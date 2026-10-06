/**
 * @file Event.cc
 * Capture event identity, exposure, triggers, vertex, centrality and
 * calorimeter context; record the established producer's event classification
 * (evt_terminal_status, recorded and never applied) and the producer-owned
 * event weights. Upstream rejection accounting remains explicit, and event
 * measurements keep their declared units and missing-value states.
 */

#include "../PhotonJetTree.h"
#include "Production.h"

#include <cmath>  // EventHeader's inline missing-value accessor uses NAN.

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/TriggerRunInfo.h>

#include <centrality/CentralityInfo.h>

#include <jetbackground/TowerBackground.h>

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>

#include <algorithm>
#include <cctype>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr int kTriggerBits = 64;

bool finiteValue(const double value) { return std::isfinite(value); }


// ============================================================================
// GL1 packet helpers
//
// The concrete packet version determines which trigger words and scaler
// accessors have distinct meaning. Gl1Packetv1 exposes only the legacy
// interface; v2/v3 expose live and scaled decision words and uint64 scalers.
//
// In v2/v3 the legacy trigger-vector accessor aliases the live word. It must
// not be interpreted as an additional independently measured raw decision.
// ============================================================================

enum TriggerAvailability : std::uint32_t
{
  kLiveWordAvailable = 1U << 0U,
  kScaledWordAvailable = 1U << 1U,
  kLegacyWordAvailable = 1U << 2U
};

int packetVersion(const Gl1Packet* packet)
{
  const std::string type = packet->ClassName();
  if (type == "Gl1Packetv1") return 1;
  if (type == "Gl1Packetv2") return 2;
  if (type == "Gl1Packetv3") return 3;
  return 0;
}

using ScalerArray = std::array<std::uint64_t, kTriggerBits>;

/**
 * Read one complete 64-bit scaler family.
 *
 * type:
 *   0 -> raw
 *   1 -> live
 *   2 -> scaled
 *
 * The accessor difference is version-dependent, but the stored representation
 * is always an unsigned 64-bit counter array.
 */
ScalerArray readScalers(Gl1Packet* packet, const int version, const int type)
{
  ScalerArray result{};

  for (int bit = 0; bit < kTriggerBits; ++bit)
  {
    result[static_cast<std::size_t>(bit)] =
        version == 1 ? static_cast<std::uint64_t>(packet->lValue(bit, type)) : packet->getScaler(bit, type);
  }

  return result;
}

/**
 * Detect a discontinuity in cumulative counters.
 *
 * Scalers are cumulative. A decrease therefore marks a reset/discontinuity
 * rather than an ordinary continuation of the previous snapshot.
 */
bool anyCounterDecreased(const ScalerArray& previous, const ScalerArray& current)
{
  for (std::size_t i = 0; i < current.size(); ++i)
  {
    if (current[i] < previous[i])
    {
      return true;
    }
  }

  return false;
}


// ============================================================================
// Calorimeter event sums
//
// A useful partial sum and complete-input validity are separate facts.
// Malformed/null/nonfinite towers are counted explicitly instead of forcing the
// entire event sum to NaN. The valid flag records whether the input census was
// complete enough to regard the sum as fully qualified.
// ============================================================================

struct CaloSum
{
  double energy = 0.0;

  std::uint64_t towerCount = 0;
  std::uint64_t acceptedCount = 0;
  std::uint64_t badQualityCount = 0;
  std::uint64_t nullCount = 0;
  std::uint64_t nonFiniteCount = 0;

  bool present = false;
  bool valid = false;
};

/**
 * Sum one calibrated calorimeter tower container.
 *
 * When requireGood is true, towers marked bad by the upstream calibration are
 * excluded from the energy sum but remain visible in the quality census.
 * Energies are summed with their stored sign; no analysis-level floor is
 * introduced here.
 */
CaloSum sumCalorimeter(PHCompositeNode* topNode, const std::string& node, const bool requireGood)
{
  CaloSum result;

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, node);
  if (!towers)
  {
    return result;
  }

  result.present = true;
  // TowerInfoContainer::size() scans the whole TClonesArray; count once.
  const auto towerCount = towers->size();
  result.towerCount = towerCount;
  result.valid = towerCount > 0;

  for (unsigned int channel = 0; channel < towerCount; ++channel)
  {
    TowerInfo* tower = towers->get_tower_at_channel(channel);

    if (!tower)
    {
      ++result.nullCount;
      result.valid = false;
      continue;
    }

    if (requireGood && !tower->get_isGood())
    {
      ++result.badQualityCount;
      continue;
    }

    const float energy = tower->get_energy();
    if (!finiteValue(energy))
    {
      ++result.nonFiniteCount;
      result.valid = false;
      continue;
    }

    result.energy += static_cast<double>(energy);  // signed, no floor
    ++result.acceptedCount;
  }

  // The retained census treats a present nonempty all-masked container as a
  // valid measured zero. Only null/nonfinite inputs or binary32 overflow lose
  // validity; the native branch widens the finished binary32 result.
  result.valid = result.valid && std::isfinite(static_cast<float>(result.energy));

  return result;
}


// ============================================================================
// Centrality helpers
// ============================================================================
// Capture the calibrated centrality inputs independently of the collaborator
// charge fallback. Preserve the reconstruction's float accumulation and its
// inclusive 0.5 charge / 25 ns selection boundaries.
template <class Pmts, class Out, class MinimumBias>
void captureCentralityInputs(photonjet::EventRecord& row, Pmts* pmts, Out* out, MinimumBias* minimumBias)
{
  using photonjet::TriState;
  row.centralityReplayVersion = 1;
  row.centralityPmtAvailable = pmts ? TriState::True : TriState::False;
  row.minimumBiasDecision =
      minimumBias ? (minimumBias->isAuAuMinimumBias() ? TriState::True : TriState::False) : TriState::Unknown;
  if (out)
  {
    row.centralityMbdZ = out->get_zvtx();
    row.centralityMbdZValid = std::isfinite(row.centralityMbdZ);
    row.centralityMbdEvent = out->get_evt();
    row.centralityMbdClock = out->get_clock();
    row.centralityMbdFemclock = out->get_femclock();
  }
  if (!pmts) return;
  if (pmts->get_npmt() != 128)
  {
    row.centralityPmtAvailable = TriState::Unknown;
    return;
  }
  bool valid = true;
  float selectedCharge = 0.f;
  for (int channel = 0; channel < 128; ++channel)
  {
    const auto* hit = pmts->get_pmt(channel);
    const float charge = hit ? hit->get_q() : std::numeric_limits<float>::quiet_NaN();
    const float time = hit ? hit->get_time() : std::numeric_limits<float>::quiet_NaN();
    const bool good = hit && hit->get_pmt() == channel && std::isfinite(charge) && std::isfinite(time);
    const bool selected = good && charge >= .5f && std::fabs(time) <= 25.f;
    row.centralityPmtId.push_back(channel);
    row.centralityPmtCharge.push_back(charge);
    row.centralityPmtTimeNs.push_back(time);
    row.centralityPmtValid.push_back(good);
    row.centralityPmtSelected.push_back(selected);
    valid = valid && good;
    if (selected) selectedCharge += charge;
  }
  if (valid) row.centralitySelectedCharge = selectedCharge;
  row.centralityInputsValid = valid && row.centralityMbdZValid && minimumBias;
}
/**
 * Map a finite centrality percentile onto the configured half-open bins.
 *
 * Values outside the configured support return -1 rather than being clipped
 * into a nearby analysis interval.
 */
int coarseCentralityBin(const double percent, const std::vector<int>& edges)
{
  if (!finiteValue(percent))
  {
    return -1;
  }

  for (std::size_t i = 0; i + 1 < edges.size(); ++i)
  {
    if (percent >= edges[i] && percent < edges[i + 1])
    {
      return static_cast<int>(i);
    }
  }

  return -1;
}

}  // namespace


// ============================================================================
// Per-event reset
// ============================================================================

/**
 * Reset only state owned by the event currently being processed.
 *
 * Source-level accounting, upstream-observer state, run trigger information,
 * scaler snapshots, weight inputs, and ROOT output handles deliberately persist
 * across events.
 */
void PhotonJetTree::resetEventState()
{
  m_event = EventRecord{};

  m_photons.clear();
  m_photonCells.clear();
  m_isolationRecords.clear();
  m_isolationConstituents.clear();

  m_jets.clear();
  m_photonJetPairs.clear();

  m_truthVertices.clear();
  m_truthPhotons.clear();
  m_truthJets.clear();

  m_photonTruthLinks.clear();
  m_photonDistanceLinks.clear();
  m_jetTruthLinks.clear();

  m_weightComponents.clear();

  m_ppg12Photons.clear();
  m_ppg12Clusters.clear();
  m_ppg12ResponseClusters.clear();
  m_ppg12ResponseClustersPresent = false;
  m_ppg12TruthJetRadiiFound.clear();
  m_ppg12EvaluatorReady = false;
}


// ============================================================================
// Upstream exposure accounting
//
// The production layer places an observer before reconstruction. Every physical
// input event is therefore observed even if a later reconstruction module
// aborts it before PhotonJetTree::process_event().
//
// A pending observation is consumed when process_event() is reached. If another
// observation arrives first, the previous event is known to have been rejected
// upstream and is recorded separately.
// ============================================================================

void PhotonJetTree::observeUpstreamEvent(const std::int64_t physicalEventSequence)
{
  if (m_upstreamPending)
  {
    const std::int64_t sourceEntry =
        m_config.source.firstEntry + static_cast<std::int64_t>(m_upstreamObservedEvents) - 1;

    recordUpstreamRejectedEvent(m_upstreamPendingRun, m_upstreamPendingSequence, sourceEntry);
  }

  m_upstreamPending = true;
  m_upstreamPendingRun = m_config.source.run;
  m_upstreamPendingSequence = physicalEventSequence;
  ++m_upstreamObservedEvents;
  photonjet::production::WriteProgress(m_config, "event_loop", static_cast<std::int64_t>(m_upstreamObservedEvents),
                                       m_config.source.firstEntry);
}

/**
 * Resolve the final observer entry after the event loop ends.
 *
 * A still-pending observation represents a physical input event that was seen
 * before reconstruction but never reached this producer.
 */
void PhotonJetTree::finishUpstreamAccounting()
{
  if (m_upstreamPending)
  {
    const std::int64_t sourceEntry =
        m_config.source.firstEntry + static_cast<std::int64_t>(m_upstreamObservedEvents) - 1;

    recordUpstreamRejectedEvent(m_upstreamPendingRun, m_upstreamPendingSequence, sourceEntry);

    m_upstreamPending = false;
  }
}


// ============================================================================
// Source and event identity
// ============================================================================

/**
 * Bind the current event to its physical source coordinates.
 *
 * The configured source/run/segment are the production-level authority. The
 * EventHeader contributes the physical event sequence when available; its run
 * was already matched to the conditions run by the exposure observer (for
 * embedded simulation, the data run the simulation is embedded in).
 *
 * sourceEntry is the original source cursor used to build the stable event
 * identity. A configured resolver handles layouts where producer encounter
 * order no longer maps directly onto that physical cursor.
 */
void PhotonJetTree::captureEventIdentity(PHCompositeNode* topNode)
{
  m_event.sourceId = m_source.sourceId;
  m_event.run = m_config.source.run;
  m_event.segment = m_config.source.segment;
  m_event.sourceFileOrdinal = m_config.source.sourceFileOrdinal;

  // --------------------------------------------------------------------------
  // Physical EventHeader
  // --------------------------------------------------------------------------

  EventHeader* header = findNode::getClass<EventHeader>(topNode, m_config.nodes.eventHeader);

  if (header && header->isValid())
  {
    m_event.physicalEventSequence = header->get_EvtSequence();
    m_event.physicalEventSequenceValid = m_event.physicalEventSequence >= 0;
  }

  // --------------------------------------------------------------------------
  // Original source entry
  //
  // Straight sequential inputs use firstEntry plus encounter position. Sources
  // affected by synchronization/skipping may supply a resolver that maps the
  // physical event identity back onto the authoritative original cursor.
  // --------------------------------------------------------------------------

  if (m_config.source.sourceEntryResolver)
  {
    m_event.sourceEntry = m_config.source.sourceEntryResolver(m_event.run, m_event.physicalEventSequence,
                                                              m_event.physicalEventSequenceValid);
  }
  else
  {
    m_event.sourceEntry = m_config.source.firstEntry + static_cast<std::int64_t>(m_upstreamObservedEvents > 0
                                                                                     ? m_upstreamObservedEvents - 1
                                                                                     : m_event.producerEventOrdinal);
  }

  // Position in the maintained global source order, when the catalogue binds
  // one; the PPG12 half-sample convention keys on it.
  m_event.sourceGlobalEntry =
      m_config.source.sourceGlobalEntryBegin >= 0 ? m_config.source.sourceGlobalEntryBegin + m_event.sourceEntry : -1;

  m_event.eventId = photonjet::makeEventIdentity(m_event.sourceId, m_event.sourceEntry, m_event.physicalEventSequence,
                                                 m_event.physicalEventSequenceValid);
}


// ============================================================================
// Trigger decisions and scaler snapshots
//
// Simulation has no detector trigger state. DATA captures the packet version,
// decision-word availability, packet metadata, and cumulative scaler state.
//
// Scalers are persisted as run-length-like snapshots: a new row is created only
// when the counters/state change or a discontinuity is observed.
// ============================================================================

void PhotonJetTree::captureTriggerInformation(PHCompositeNode* topNode)
{
  m_event.triggerPacketValid = false;
  m_event.triggerDecisionsAvailable = 0;
  m_event.triggerScalerSnapshotId = -1;
  m_event.triggerScalerDiscontinuityFlags = 8U | (!m_event.physicalEventSequenceValid ? 16U : 0U);

  if (isSimulation())
  {
    m_event.triggerDecisionState = 5;
    return;
  }

  Gl1Packet* packet = findNode::getClass<Gl1Packet>(topNode, m_config.nodes.gl1Packet);

  if (!packet)
  {
    // Some DATA files expose the GL1 packet under its numeric packet-node name.
    packet = findNode::getClass<Gl1Packet>(topNode, "14001");
  }

  if (!packet)
  {
    return;
  }

  const int version = packetVersion(packet);
  m_event.triggerPacketVersion = version;

  if (version == 0)
  {
    m_event.triggerDecisionState = 2;
    return;
  }

  m_event.triggerPacketValid = true;
  m_event.triggerPacketStatus = packet->getStatus();
  m_event.triggerDecisionState = m_event.triggerPacketStatus ? 4 : version >= 2 ? 0 : 2;
  m_event.triggerPacketNumber = packet->getPacketNumber();
  // The stored gl1_bco and snapshot endpoints are the full beam clock.
  // getBunchNumber() is the cycling bunch slot, not the monotonic GL1 BCO.
  m_event.triggerBunchNumber = packet->getBCO();

  /*
   * getTriggerVector() is the legacy accessor. Its availability is recorded
   * independently of the live/scaled words exposed by packet versions >= 2.
   */
  m_event.triggerInputBits = packet->getTriggerVector();
  m_event.triggerDecisionsAvailable |= kLegacyWordAvailable;

  if (version >= 2)
  {
    m_event.triggerLiveBits = packet->getLiveVector();
    m_event.triggerScaledBits = packet->getScaledVector();
    m_event.triggerDecisionsAvailable |= kLiveWordAvailable | kScaledWordAvailable;
  }

  // --------------------------------------------------------------------------
  // Cumulative scaler state
  // --------------------------------------------------------------------------

  const ScalerArray raw = readScalers(packet, version, 0);
  const ScalerArray live = readScalers(packet, version, 1);
  const ScalerArray scaled = readScalers(packet, version, 2);

  const std::uint64_t bunch = packet->getBCO();
  const auto unknownBunch = std::numeric_limits<std::uint64_t>::max();
  m_event.triggerScalerDiscontinuityFlags =
      (bunch == unknownBunch ? 8U : 0U) | (!m_event.physicalEventSequenceValid ? 16U : 0U);
  if (!m_triggerScalerSnapshots.empty())
  {
    const auto& previous = m_triggerScalerSnapshots.back();
    if (previous.bunchEnd != unknownBunch && bunch != unknownBunch && bunch < previous.bunchEnd)
      m_event.triggerScalerDiscontinuityFlags |= 2U;
    if (previous.physicalEventEnd >= 0 && m_event.physicalEventSequenceValid &&
        m_event.physicalEventSequence <= previous.physicalEventEnd)
      m_event.triggerScalerDiscontinuityFlags |= 4U;
    if (previous.valid && m_event.triggerPacketStatus == 0 &&
        (anyCounterDecreased(previous.raw, raw) || anyCounterDecreased(previous.live, live) ||
         anyCounterDecreased(previous.scaled, scaled)))
      m_event.triggerScalerDiscontinuityFlags |= 1U;
  }

  bool newSnapshot = m_triggerScalerSnapshots.empty();

  if (!newSnapshot)
  {
    const TriggerScalerSnapshotRecord& previous = m_triggerScalerSnapshots.back();

    newSnapshot = previous.raw != raw || previous.live != live || previous.scaled != scaled || !previous.valid ||
                  packet->getStatus() != 0 || bunch < previous.bunchEnd;
  }

  if (newSnapshot)
  {
    TriggerScalerSnapshotRecord snapshot;

    snapshot.sourceId = m_source.sourceId;
    snapshot.snapshotId = static_cast<std::int64_t>(m_triggerScalerSnapshots.size());

    snapshot.bunchStart = bunch;
    snapshot.bunchEnd = bunch;

    snapshot.sourceEntryStart = m_event.sourceEntry;
    snapshot.sourceEntryEnd = m_event.sourceEntry;
    snapshot.observedEvents = 1;
    snapshot.physicalEventEnd = m_event.physicalEventSequenceValid ? m_event.physicalEventSequence : -1;

    snapshot.raw = raw;
    snapshot.live = live;
    snapshot.scaled = scaled;

    snapshot.valid = m_event.triggerPacketStatus == 0;

    if (!m_triggerScalerSnapshots.empty())
    {
      const TriggerScalerSnapshotRecord& previous = m_triggerScalerSnapshots.back();

      snapshot.discontinuity = anyCounterDecreased(previous.raw, raw) || anyCounterDecreased(previous.live, live) ||
                               anyCounterDecreased(previous.scaled, scaled) || bunch < previous.bunchEnd;
    }

    m_triggerScalerSnapshots.push_back(snapshot);
  }
  else
  {
    TriggerScalerSnapshotRecord& snapshot = m_triggerScalerSnapshots.back();

    snapshot.bunchEnd = bunch;
    snapshot.sourceEntryEnd = m_event.sourceEntry;
    snapshot.physicalEventEnd = m_event.physicalEventSequenceValid ? m_event.physicalEventSequence : -1;
    ++snapshot.observedEvents;
  }

  m_event.triggerScalerSnapshotId = m_triggerScalerSnapshots.back().snapshotId;
}


// ============================================================================
// Run-level trigger configuration
// ============================================================================

/**
 * Capture one row per GL1 trigger bit for the current DATA run.
 *
 * Only the explicitly supported TriggerRunInfov1 interface is interpreted.
 * Missing/unsupported run information therefore yields invalid rows rather
 * than silently disappearing, preserving the distinction between unavailable
 * configuration and a configuration containing no assigned trigger.
 */
void PhotonJetTree::captureTriggerRunInformation(PHCompositeNode* topNode)
{
  auto* info = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");

  const bool supported = info && std::string(info->ClassName()) == "TriggerRunInfov1";

  for (int bit = 0; bit < kTriggerBits; ++bit)
  {
    TriggerRunInfoRecord row;

    row.sourceId = m_source.sourceId;
    row.run = m_config.source.run;
    row.bit = bit;
    row.state = info ? 2 : 1;

    if (supported)
    {
      row.name = info->getTriggerName(bit);
      row.hardwarePrescale = info->getInitialPrescaleByBit(bit);
      row.runAveragePrescale = info->getPrescaleByBit(bit);
      row.rawCount = info->getRawScalersByBit(bit);
      row.liveCount = info->getLiveScalersByBit(bit);
      row.scaledCount = info->getScalersByBit(bit);

      row.valid = !row.name.empty() && row.name != "unknown";
      row.state = row.valid ? 0 : 3;
    }

    m_triggerRunInfo.push_back(std::move(row));
  }
}


// ============================================================================
// Reconstructed vertex
// ============================================================================

/**
 * Capture the reconstructed event vertex z of the established productions.
 *
 *   pp DATA  the first MbdVertexMap vertex (the PPG12 data convention)
 *   pp SIM   the last finite MBD vertex of the first GlobalVertex (PPG12 SIM)
 *   Au+Au    the first MbdVertexMap vertex, else the first GlobalVertex
 *
 * A missing source leaves the vertex invalid (NaN); truth is never
 * substituted. The value is the stored binary32 coordinate.
 */
void PhotonJetTree::captureVertexInformation(PHCompositeNode* topNode)
{
  m_event.recoVertexValid = false;
  m_event.recoVertexSource = 0;

  auto* mbdVertices = findNode::getClass<MbdVertexMap>(topNode, m_config.nodes.mbdVertexMap);
  auto* globalVertices = findNode::getClass<GlobalVertexMap>(topNode, m_config.nodes.globalVertexMap);
  const GlobalVertex* global = globalVertices && !globalVertices->empty() ? globalVertices->begin()->second : nullptr;

  auto accept = [&](const double z, const std::int32_t source)
  {
    if (!finiteValue(z)) return false;
    m_event.recoVertexX = global ? global->get_x() : kNaN;
    m_event.recoVertexY = global ? global->get_y() : kNaN;
    m_event.recoVertexZ = static_cast<float>(z);
    m_event.recoVertexValid = true;
    m_event.recoVertexSource = source;
    return true;
  };

  if (isPP() && isSimulation())
  {
    double z = kNaN;
    if (global)
      for (auto it = global->find_vertexes(GlobalVertex::MBD); it != global->end_vertexes(); ++it)
        if (it->first == GlobalVertex::MBD)
          for (const auto* vertex : it->second)
            if (vertex && finiteValue(vertex->get_z())) z = vertex->get_z();
    accept(z, 2);
    return;
  }

  const MbdVertex* mbd = mbdVertices && !mbdVertices->empty() ? mbdVertices->begin()->second : nullptr;
  if (mbd && accept(mbd->get_z(), 1)) return;
  if (isAuAu() && global) accept(global->get_z(), 2);
}


// ============================================================================
// Minimum-bias detector
//
// Event-level MBD timing/charge and PMT-level information are stored
// independently. PMT channel order is preserved exactly; malformed or missing
// channels are represented through validity rather than removed from the row.
// ============================================================================

void PhotonJetTree::captureMbdInformation(PHCompositeNode* topNode)
{
  m_event.mbdValid = false;
  m_event.mbdPmtAvailable = TriState::Unknown;

  if (MbdOut* out = findNode::getClass<MbdOut>(topNode, m_config.nodes.mbdOut); out)
  {
    m_event.mbdAvailable = true;
    m_event.mbdNativeValid = out->isValid();
    m_event.mbdSouthNpmt = out->get_npmt(0);
    m_event.mbdNorthNpmt = out->get_npmt(1);
    m_event.mbdT0Ns = out->get_t0();
    m_event.mbdSouthTimeNs = out->get_time(0);
    m_event.mbdNorthTimeNs = out->get_time(1);

    m_event.mbdSouthCharge = out->get_q(0);
    m_event.mbdNorthCharge = out->get_q(1);

    m_event.mbdValid = finiteValue(m_event.mbdT0Ns) && finiteValue(m_event.mbdSouthTimeNs) &&
                       finiteValue(m_event.mbdNorthTimeNs) && finiteValue(m_event.mbdSouthCharge) &&
                       finiteValue(m_event.mbdNorthCharge);
  }

  MbdPmtContainer* pmts = findNode::getClass<MbdPmtContainer>(topNode, m_config.nodes.mbdPmts);

  if (!pmts && isAuAu() && !m_config.nodes.mbdPmtsFallback.empty())
    pmts = findNode::getClass<MbdPmtContainer>(topNode, m_config.nodes.mbdPmtsFallback);

  if (!pmts)
  {
    m_event.mbdPmtAvailable = TriState::False;
    return;
  }

  const int count = pmts->get_npmt();

  if (count < 0 || count > 128)
  {
    // A malformed container is neither a valid presence nor a true absence.
    m_event.mbdPmtAvailable = TriState::Unknown;
    return;
  }

  m_event.mbdPmtAvailable = TriState::True;

  MbdGeom* geometry = findNode::getClass<MbdGeom>(topNode, m_config.nodes.mbdGeometry);

  m_event.mbdPmtId.reserve(count);
  m_event.mbdPmtArm.reserve(count);
  m_event.mbdPmtCharge.reserve(count);
  m_event.mbdPmtTimeNs.reserve(count);
  m_event.mbdPmtValid.reserve(count);

  // The established event charge: the sequential binary64 sum of the finite
  // positive calibrated PMT charges, in channel order, with no time or
  // centrality threshold. A complete array has all 128 channels valid.
  double positiveCharge = 0.0;
  m_event.mbdPmtArrayComplete = count == 128;

  for (int channel = 0; channel < count; ++channel)
  {
    const MbdPmtHit* hit = pmts->get_pmt(channel);

    const double charge = hit ? hit->get_q() : kNaN;
    const double time = hit ? hit->get_time() : kNaN;
    const int arm = geometry ? geometry->get_arm(channel) : -1;

    /*
     * Preserve the original container channel as the PMT identity. The vectors
     * are never compacted, so holes and malformed channels remain visible.
     */
    m_event.mbdPmtId.push_back(channel);
    m_event.mbdPmtArm.push_back(arm == 0 || arm == 1 ? arm : -1);
    m_event.mbdPmtCharge.push_back(charge);  // zero and negative kept
    m_event.mbdPmtTimeNs.push_back(time);

    const bool valid = hit && hit->get_pmt() == channel && finiteValue(charge);
    m_event.mbdPmtValid.push_back(valid ? 1 : 0);
    m_event.mbdPmtArrayComplete = m_event.mbdPmtArrayComplete && valid;
    if (finiteValue(charge) && charge > 0.0) positiveCharge += charge;
  }
  m_event.mbdTotalCharge = positiveCharge;
}


// ============================================================================
// Au+Au centrality
//
// The reconstructed object and whether it is trustworthy are stored as separate
// facts. In DATA, the centrality value is accepted for the event only when the
// reconstructed minimum-bias decision is positive. Embedded simulation uses
// the centrality state carried by its embedded input.
//
// This prevents a stale CentralityInfo object from a prior event from being
// mistaken for a valid measurement.
// ============================================================================

void PhotonJetTree::captureCentrality(PHCompositeNode* topNode)
{
  m_event.centralityApplicable = isAuAu();
  captureCentralityInputs(
      m_event, findNode::getClass<MbdPmtContainer>(topNode, m_config.nodes.mbdPmts),
      findNode::getClass<MbdOut>(topNode, m_config.nodes.mbdOut),
      isAuAu() ? findNode::getClass<MinimumBiasInfo>(topNode, m_config.nodes.minimumBiasInfo) : nullptr);

  if (!isAuAu())
  {
    return;
  }

  auto* centrality = findNode::getClass<CentralityInfo>(topNode, m_config.nodes.centralityInfo);

  // CentralityReco fills CentralityInfo only for minimum-bias events and
  // CentralityInfov2 never resets it: in DATA the node describes this event
  // only when the minimum-bias classifier passed. Embedded simulation carries
  // its input sample's per-event centrality.
  constexpr auto property = CentralityInfo::PROP::mbd_NS;
  const bool writtenForThisEvent = isEmbeddedSimulation() || m_event.minimumBiasDecision == TriState::True;

  if (centrality && writtenForThisEvent)
  {
    const double bin = centrality->get_centrality_bin(property);
    if (finiteValue(bin) && bin >= 0.0 && bin <= 100.0) m_event.analysisCentrality = static_cast<std::int32_t>(bin);
  }

  if (!writtenForThisEvent || !centrality || !centrality->has_centile(property) ||
      !centrality->has_centrality_bin(property))
    return;
  const double centile = centrality->get_centile(property);
  const auto bin = centrality->get_centrality_bin(property);
  m_event.centralityNativeCentile = centile;
  m_event.centralityNativeBin = bin;
  m_event.centralityNativeValid = finiteValue(centile) && centile > 0. && centile <= 1. && bin >= 1 && bin <= 100;
  if (!m_event.centralityNativeValid) return;
  m_event.centralityPercent = bin;
  m_event.centralityBin = coarseCentralityBin(bin, m_config.centrality.edges);
  m_event.centralityValid = true;
}


// ============================================================================
// Calorimeter event state
//
// The producer stores both aggregate energies and the census explaining how
// those sums were formed. Availability and validity are independent masks.
//
// For Au+Au, the background object used by the reconstructed jet/isolation
// chain is also captured so the subtraction state remains inspectable later.
// ============================================================================

void PhotonJetTree::captureCalorimeterInformation(PHCompositeNode* topNode)
{
  auto* rawClusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.rawClusterNode);
  if (!rawClusters) rawClusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.photonNode);
  if (rawClusters)
  {
    m_event.maximumRawClusterEnergy = 0.;
    for (const auto& item : rawClusters->getClustersMap())
      if (item.second)
      {
        const double energy = item.second->get_energy();
        if (std::isfinite(energy) && energy >= 1.)
          m_event.maximumRawClusterEnergy = std::max(m_event.maximumRawClusterEnergy, energy);
      }
  }
  const bool requireGood = m_config.calorimeter.requireGoodTowers;

  const std::array<CaloSum, 3> sums = {sumCalorimeter(topNode, m_config.nodes.cemcTowers, requireGood),

                                       sumCalorimeter(topNode, m_config.nodes.ihcalTowers, requireGood),

                                       sumCalorimeter(topNode, m_config.nodes.ohcalTowers, requireGood)};

  m_event.caloRequiredGoodTowers = requireGood;

  m_event.caloAvailableMask = 0;
  m_event.caloValidMask = 0;

  double total = 0.0;
  bool allPresent = true;

  for (std::size_t layer = 0; layer < sums.size(); ++layer)
  {
    const CaloSum& sum = sums[layer];

    m_event.caloTowerCount[layer] = sum.towerCount;
    m_event.caloAcceptedTowerCount[layer] = sum.acceptedCount;
    m_event.caloBadQualityTowerCount[layer] = sum.badQualityCount;
    m_event.caloNullTowerCount[layer] = sum.nullCount;
    m_event.caloNonFiniteTowerCount[layer] = sum.nonFiniteCount;

    if (sum.present)
    {
      m_event.caloAvailableMask |= 1U << layer;
      total += sum.energy;
    }
    else
      allPresent = false;

    if (sum.valid)
    {
      m_event.caloValidMask |= 1U << (8U + layer);
    }
  }

  m_event.cemcEnergy = sums[0].present ? static_cast<float>(sums[0].energy) : kNaN;

  m_event.ihcalEnergy = sums[1].present ? static_cast<float>(sums[1].energy) : kNaN;

  m_event.ohcalEnergy = sums[2].present ? static_cast<float>(sums[2].energy) : kNaN;

  m_event.totalCaloEnergy = allPresent ? static_cast<float>(total) : kNaN;

  if (isAuAu())
  {
    /*
     * Retain the exact event-level background witnesses consumed by the Au+Au
     * reconstructed-object chain. The node is required in this production mode:
     * silently omitting it would make later SUB1 interpretation ambiguous.
     *
     * Historical source anchor:
     *   RecoilJets_AuAu.cc:7712-7735
     */
    auto* background = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");

    if (!background)
    {
      throw std::runtime_error("required AuAu TowerInfoBackground_Sub2 is missing");
    }

    m_event.jetBackgroundState = CaptureState::Complete;
    m_event.jetBackgroundFlowMode = m_config.calorimeter.backgroundFlowMode;

    m_event.jetBackgroundV2 = background->get_v2();
    m_event.jetBackgroundPsi2 = background->get_Psi2();

    m_event.jetBackgroundUEEmcal = background->get_UE(0);
    m_event.jetBackgroundUEIhcal = background->get_UE(1);
    m_event.jetBackgroundUEOhcal = background->get_UE(2);

    m_event.jetBackgroundNStrips = background->get_nStripsUsedForFlow();

    m_event.jetBackgroundNTowers = background->get_nTowersUsedForBkg();

    m_event.jetBackgroundFlowFailure = background->get_flow_failure_flag() ? TriState::True : TriState::False;
  }
}


// ============================================================================
// Historical event classification (pp)
//
// The established producer classified every event before its analysis fills:
// a mandatory reconstructed vertex, the active trigger streams, the PPG12
// period run filter and the |vz| cut. The tree records the outcome as
// evt_terminal_status; the PPG12 vertex/period audit occurrences it logged
// come from PPG12.cc, in their original order. It removes nothing.
// ============================================================================

/**
 * True when exactly one GL1 bit carries the named trigger in the run table
 * and that bit's scaled decision fired in a valid packet.
 */
bool PhotonJetTree::scaledTriggerFired(const std::string& name) const
{
  if (name.empty() || name == "unknown" || m_event.triggerDecisionState != 0) return false;
  int bit = -1;
  for (const auto& row : m_triggerRunInfo)
  {
    if (row.state != 0 && row.state != 3) return false;  // no supported run table
    if (row.name != name) continue;
    if (bit >= 0) return false;  // an ambiguous name is never a decision
    bit = row.bit;
  }
  return bit >= 0 && ((m_event.triggerScaledBits >> bit) & 1U);
}

void PhotonJetTree::classifyEvent()
{
  const photonjet::EventClassificationConfig& c = m_config.classification;
  m_event.terminalStatus = 0;
  m_event.activeTriggers.clear();
  m_event.nodesReady = m_event.recoVertexValid;
  if (!c.enabled || !m_event.nodesReady) return;
  if (isAuAu())
  {
    classifyAuAuEvent();
    return;
  }

  const bool outsideVertexCut = c.useVertexCut && std::fabs(m_event.recoVertexZ) >= c.vertexCutCm;

  if (isSimulation())
  {
    if (!recordPPG12SampleOwnership(true))
    {
      m_event.terminalStatus = 4;
      return;
    }
    m_event.activeTriggers.push_back("SIM");
    recordPPG12VertexAudit(true);
    if (outsideVertexCut)
    {
      m_event.terminalStatus = 2;
      return;
    }
    recordPPG12VertexAudit(false);
    recordPPG12SampleOwnership(false);
    return;
  }

  for (const auto& [name, stream] : c.triggers)
    if (scaledTriggerFired(name)) m_event.activeTriggers.push_back(stream);
  // The PPG12 scaled-bit audit stream reads the packet's scaled word as is.
  if (!c.scaledBitStream.empty() && m_event.triggerPacketVersion >= 2 &&
      ((m_event.triggerScaledBits >> c.scaledBit) & 1U) &&
      std::find(m_event.activeTriggers.begin(), m_event.activeTriggers.end(), c.scaledBitStream) ==
          m_event.activeTriggers.end())
    m_event.activeTriggers.push_back(c.scaledBitStream);
  if (m_event.activeTriggers.empty())
  {
    m_event.terminalStatus = 1;
    return;
  }

  recordPPG12VertexAudit(true);
  if (c.filterDataByPeriod &&
      !(!c.period.key.empty() && m_event.run >= c.period.runMin && m_event.run < c.period.runMaxExclusive))
  {
    m_event.terminalStatus = 3;
    return;
  }
  if (outsideVertexCut)
  {
    m_event.terminalStatus = 2;
    return;
  }
  recordPPG12VertexAudit(false);
}


// ============================================================================
// Historical event classification (Au+Au)
//
// The established Au+Au producer's gates, in its order. DATA: a configured
// scaled trigger, the minimum-bias classifier, |vz|. Embedded simulation:
// |vz|, then generator-slice stitching (logged as the STITCH occurrence).
// Both: the truncated centrality inside the analysis range. An accepted
// embedded event then carries the truth-isolation audit and the truth-photon
// misses (Truth.cc).
// ============================================================================

int photonjet::production::EmbeddedPhotonSampleCode(std::string sample)
{
  std::transform(sample.begin(), sample.end(), sample.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  for (const int code : {12, 20})
    for (const std::string& prefix : {"photonjet", "embeddedphoton", "photon"})
      if (sample.find(prefix + std::to_string(code)) != std::string::npos) return code;
  return 0;
}

// Narrower than the photon names: "PhotonJet12" is not an inclusive jet sample.
int photonjet::production::EmbeddedJetSampleCode(std::string sample)
{
  std::transform(sample.begin(), sample.end(), sample.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (sample.find("photon") != std::string::npos) return 0;
  for (const int code : {12, 20, 30, 40})
    if (sample.find("embeddedjet" + std::to_string(code)) != std::string::npos) return code;
  return 0;
}

void PhotonJetTree::classifyAuAuEvent()
{
  const photonjet::EventClassificationConfig& c = m_config.classification;
  const bool outsideVertexCut = c.useVertexCut && std::fabs(m_event.recoVertexZ) >= c.vertexCutCm;

  if (isSimulation())
  {
    m_event.activeTriggers.push_back("SIM");
    if (outsideVertexCut)
    {
      m_event.terminalStatus = 2;
      return;
    }
    if (!embeddedSampleOwnsEvent())
    {
      m_event.terminalStatus = 4;
      return;
    }
  }
  else
  {
    for (const auto& [bit, stream] : c.scaledTriggers)
      if ((m_event.triggerScaledBits >> bit) & 1U) m_event.activeTriggers.push_back(stream);
    if (m_event.activeTriggers.empty())
    {
      m_event.terminalStatus = 1;
      return;
    }
    if (c.requireMinimumBias && m_event.minimumBiasDecision != TriState::True)
    {
      m_event.terminalStatus = 3;
      return;
    }
    if (outsideVertexCut)
    {
      m_event.terminalStatus = 2;
      return;
    }
  }

  if (m_event.analysisCentrality < c.centralityMinPercent || m_event.analysisCentrality >= c.centralityMaxPercent)
  {
    m_event.terminalStatus = 5;
    return;
  }
  if (isSimulation())
  {
    recordEmbeddedTruthIsolationAudit();
    recordEmbeddedTruthMisses();
  }
}

// The embedded generator slice owns the event. Photon samples: the
// generator-filter photon in [12, 20) or above 20 GeV (logged as the STITCH
// occurrence; an unknown sample stitches nothing). Jet samples: the maximum
// R = 0.4 truth-jet pT inside the sample's configured window.
bool PhotonJetTree::embeddedSampleOwnsEvent()
{
  if (isPhotonJetSimulation())
  {
    if (!m_event.embeddedFilterPhotonCaptured)
      throw std::runtime_error("embedded generator-filter photon capture did not run");
    EmbeddedDiagnosticRecord row;
    row.kind = "STITCH";
    row.sampleCode = photonjet::production::EmbeddedPhotonSampleCode(m_config.source.sample);
    row.stitchPhotonPt = m_event.embeddedFilterPhotonPt;
    const double pt = row.stitchPhotonPt;
    bool owned = true;
    row.stitchDecision = 4;
    if (row.sampleCode == 12 || row.sampleCode == 20)
    {
      owned = std::isfinite(pt) && pt >= row.sampleCode && (row.sampleCode == 20 || pt < 20.0);
      row.stitchDecision = owned ? 1 : 2;
      if (!std::isfinite(pt) || pt <= 0.0)
      {
        owned = false;
        row.stitchDecision = 3;
      }
    }
    m_event.embeddedDiagnostics.push_back(row);
    return owned;
  }

  const auto& windows = m_config.classification.embeddedJetWindows;
  const auto window = windows.find(photonjet::production::EmbeddedJetSampleCode(m_config.source.sample));
  if (window == windows.end()) return true;
  double maxPt = -1.0;
  for (const TruthJetRecord& jet : m_truthJets)
    if (std::fabs(jet.radius - 0.4) < 1e-6 && std::isfinite(jet.pt) && jet.pt > 0.0 && jet.pt > maxPt) maxPt = jet.pt;
  const auto [low, high] = window->second;
  return maxPt > 0.0 && maxPt >= low && maxPt < high;
}


// ============================================================================
// Producer-owned event weights
//
// DATA carries weight 1. pp SIM retains the established PPG12 relative slice,
// truth-vertex, interaction and period factors in w_event and separately.
// AuAu vertex/centrality fits, source ownership and absolute normalization are
// applied downstream. Aliases w_xsec/w_exposure are not additional factors.
//
// A requested factor that cannot be determined remains invalid/NaN; it is never
// silently replaced by unity.
// ============================================================================

/**
 * Load external producer-owned weight inputs once at initialization.
 */
void PhotonJetTree::loadWeightInputs()
{
  const WeightConfig& weights = m_config.weights;

  if (!weights.applyVertexReweight)
  {
    return;
  }

  TDirectory::TContext context;

  std::unique_ptr<TFile> file(TFile::Open(weights.vertexReweightFile.c_str(), "READ"));

  if (!file || file->IsZombie())
  {
    throw std::runtime_error("cannot open the vertex-reweighting file " + weights.vertexReweightFile);
  }

  auto* histogram = dynamic_cast<TH1*>(file->Get(weights.vertexReweightHistogram.c_str()));

  if (!histogram)
  {
    throw std::runtime_error("vertex-reweighting histogram '" + weights.vertexReweightHistogram + "' is absent from " +
                             weights.vertexReweightFile);
  }

  m_vertexReweight = dynamic_cast<TH1*>(histogram->Clone("photonjet_vertex_reweight"));

  m_vertexReweight->SetDirectory(nullptr);
}


/**
 * Build the producer-owned event weight and its component records.
 *
 * The final component is always written exactly once. Individual components
 * are also retained so downstream code can audit which factor was requested,
 * whether it was valid, and how the stored final weight was constructed.
 */
void PhotonJetTree::captureEventWeights(PHCompositeNode* /*topNode*/)
{
  WeightComponentRecord final;

  final.eventId = m_event.eventId;
  final.componentType = "final";
  final.applicationCount = 1;
  // Unit factors unless the pp simulation PPG12 weight below replaces them.
  final.sliceWeight = final.crossSectionWeight = 1.0;
  final.vertexWeight = 1.0;
  final.siDiWeight = final.periodWeight = 1.0;
  final.exposureWeight = 1.0;

  if (isData())
  {
    final.finalWeight = 1.0;
    final.valid = true;

    m_weightComponents.push_back(final);

    m_event.eventWeight = 1.0;
    m_event.eventWeightValid = true;

    return;
  }

  // pp simulation carries the PPG12 weight only for an event the established
  // producer weighted: owned by its sample, accepted by the classification
  // and with its truth vertices. Every other event keeps weight 1 with the
  // PPG12 weight inactive.
  if (m_config.weights.ppg12Nominal)
  {
    const auto present = [](const double z) { return finiteValue(z) && z > -9000.0; };
    if (!m_event.ppg12Owned || m_event.terminalStatus != 0 || !present(m_event.truthHardVertexZ) ||
        (m_config.weights.doubleInteraction && !present(m_event.truthMinimumBiasVertexZ)))
    {
      final.finalWeight = 1.0;
      final.valid = true;
      m_weightComponents.push_back(final);
      m_event.eventWeight = 1.0;
      m_event.eventWeightValid = true;
      return;
    }
  }

  double product = 1.0;
  bool valid = true;

  if (m_config.weights.applyVertexReweight)
  {
    WeightComponentRecord vertex;

    vertex.eventId = m_event.eventId;
    vertex.componentType = "vertex";
    vertex.applicationCount = 1;

    if (m_vertexReweight && m_event.truthHardVertexValid)
    {
      // PPG12 interpolates the truth-vertex histogram; a bin lookup changes
      // the reference weights. pp DI has independent hard and MB vertices.
      // Au+Au is never reweighted here; its vertex weight is applied downstream.
      auto lookup = [&](double z)
      {
        if (!finiteValue(z) || z <= -9000.0) return kNaN;
        if (z < m_vertexReweight->GetXaxis()->GetXmin() || z > m_vertexReweight->GetXaxis()->GetXmax()) return 0.0;
        const double w = m_vertexReweight->Interpolate(z);
        return finiteValue(w) && w > 0.0 ? w : 0.0;
      };
      double factor = lookup(m_event.truthHardVertexZ);
      if (m_config.weights.doubleInteraction)
        factor *= m_event.truthMinimumBiasVertexValid ? lookup(m_event.truthMinimumBiasVertexZ) : kNaN;
      if (finiteValue(factor) && factor >= 0.0)
      {
        vertex.vertexWeight = factor;
        vertex.finalWeight = factor;
        vertex.valid = true;
      }
    }

    m_weightComponents.push_back(vertex);

    if (vertex.valid)
    {
      product *= vertex.vertexWeight;
      final.vertexWeight = vertex.vertexWeight;
    }
    else
    {
      valid = false;
    }
  }

  // Store the factors in the same precision and multiplication order as the
  // PPG12 producer. w_xsec aliases w_slice; w_exposure aliases w_period.
  if (m_config.weights.ppg12Nominal)
  {
    final.sliceWeight = m_config.weights.sliceFactor;
    final.crossSectionWeight = final.sliceWeight;
    final.siDiWeight = m_config.weights.interactionFactor;
    final.periodWeight = m_config.weights.periodFactor;
    final.exposureWeight = final.periodWeight;
    product = final.sliceWeight * product * final.siDiWeight * final.periodWeight;
  }
  valid = valid && finiteValue(product) && product >= 0.0;
  final.finalWeight = valid ? product : kNaN;
  final.valid = valid;

  m_weightComponents.push_back(final);

  m_event.eventWeight = final.finalWeight;
  m_event.eventWeightValid = valid;
  m_event.ppg12WeightActive = m_config.weights.ppg12Nominal && valid;
}
