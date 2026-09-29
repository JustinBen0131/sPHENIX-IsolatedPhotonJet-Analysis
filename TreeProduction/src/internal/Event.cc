/**
 * @file Event.cc
 * Capture event identity, exposure, triggers, vertex and calorimeter context.
 * Producer-owned weights and upstream rejection accounting remain explicit.
 * Event measurements retain their declared units and missing-value states.
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

bool finiteValue(const double value)
{
  return std::isfinite(value);
}


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
        version == 1 ? static_cast<std::uint64_t>(packet->lValue(bit, type))
                     : packet->getScaler(bit, type);
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
  result.valid = towerCount>0;

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
  result.valid=result.valid&&std::isfinite(static_cast<float>(result.energy));

  return result;
}


// ============================================================================
// Centrality helpers
// ============================================================================
// Capture the calibrated centrality inputs independently of the collaborator
// charge fallback. Preserve the reconstruction's float accumulation and its
// inclusive 0.5 charge / 25 ns selection boundaries.
template<class Pmts, class Out, class MinimumBias>
void captureCentralityInputs(photonjet::EventRecord& row, Pmts* pmts, Out* out, MinimumBias* minimumBias)
{
  using photonjet::TriState;
  row.centralityReplayVersion = 1;
  row.centralityPmtAvailable = pmts ? TriState::True : TriState::False;
  row.minimumBiasDecision = minimumBias
      ? (minimumBias->isAuAuMinimumBias() ? TriState::True : TriState::False) : TriState::Unknown;
  if (out) {
    row.centralityMbdZ = out->get_zvtx();
    row.centralityMbdZValid = std::isfinite(row.centralityMbdZ);
    row.centralityMbdEvent = out->get_evt();
    row.centralityMbdClock = out->get_clock();
    row.centralityMbdFemclock = out->get_femclock();
  }
  if (!pmts) return;
  if (pmts->get_npmt() != 128) {
    row.centralityPmtAvailable = TriState::Unknown;
    return;
  }
  bool valid = true;
  float selectedCharge = 0.f;
  for (int channel = 0; channel < 128; ++channel) {
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
  m_jetTruthLinks.clear();

  m_weightComponents.clear();
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

void PhotonJetTree::observeUpstreamEvent(const int run, const std::int64_t physicalEventSequence)
{
  if(run>0&&m_config.source.run>0&&run!=m_config.source.run)
    throw std::runtime_error("upstream event run disagrees with the bound source");
  if (m_upstreamPending)
  {
    const std::int64_t sourceEntry =
        m_config.source.firstEntry + static_cast<std::int64_t>(m_upstreamObservedEvents) - 1;

    recordUpstreamRejectedEvent(
        m_upstreamPendingRun,
        m_upstreamPendingSequence,
        sourceEntry);
  }

  m_upstreamPending = true;
  m_upstreamPendingRun = m_config.source.run;
  m_upstreamPendingSequence = physicalEventSequence;
  ++m_upstreamObservedEvents;
  photonjet::production::WriteProgress(m_config,"event_loop",
      static_cast<std::int64_t>(m_upstreamObservedEvents),m_config.source.firstEntry);
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

    recordUpstreamRejectedEvent(
        m_upstreamPendingRun,
        m_upstreamPendingSequence,
        sourceEntry);

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
 * EventHeader contributes the physical event sequence when available and must
 * agree with the configured run.
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

  EventHeader* header =
      findNode::getClass<EventHeader>(topNode, m_config.nodes.eventHeader);

  if (header && header->isValid())
  {
    m_event.physicalEventSequence = header->get_EvtSequence();
    m_event.physicalEventSequenceValid = m_event.physicalEventSequence>=0;

    /*
     * The configured run is the production binding. When EventHeader also
     * carries a physical run number, disagreement is a source-contract error.
     */
    const int headerRun = header->get_RunNumber();
    if (headerRun > 0 && m_event.run > 0 && headerRun != m_event.run)
    {
      throw std::runtime_error(
          "event header run " + std::to_string(headerRun) +
          " disagrees with the configured run " + std::to_string(m_event.run));
    }
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
    m_event.sourceEntry = m_config.source.sourceEntryResolver(
        m_event.run,
        m_event.physicalEventSequence,
        m_event.physicalEventSequenceValid);
  }
  else
  {
    m_event.sourceEntry =
        m_config.source.firstEntry + static_cast<std::int64_t>(
            m_upstreamObservedEvents > 0
                ? m_upstreamObservedEvents - 1
                : m_event.producerEventOrdinal);
  }

  m_event.sourceGlobalEntry = m_event.sourceEntry;

  m_event.eventId = photonjet::makeEventIdentity(
      m_event.sourceId,
      m_event.sourceEntry,
      m_event.physicalEventSequence,
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
  m_event.triggerScalerDiscontinuityFlags=8U|(!m_event.physicalEventSequenceValid?16U:0U);

  if (isSimulation())
  {
    m_event.triggerDecisionState=5;
    return;
  }

  Gl1Packet* packet =
      findNode::getClass<Gl1Packet>(topNode, m_config.nodes.gl1Packet);

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
    m_event.triggerDecisionState=2;
    return;
  }

  m_event.triggerPacketValid = true;
  m_event.triggerPacketStatus = packet->getStatus();
  m_event.triggerDecisionState=m_event.triggerPacketStatus?4:version>=2?0:2;
  m_event.triggerPacketNumber = packet->getPacketNumber();
  m_event.triggerBunchNumber = packet->getBunchNumber();

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
    m_event.triggerDecisionsAvailable |=
        kLiveWordAvailable | kScaledWordAvailable;
  }

  // --------------------------------------------------------------------------
  // Cumulative scaler state
  // --------------------------------------------------------------------------

  const ScalerArray raw = readScalers(packet, version, 0);
  const ScalerArray live = readScalers(packet, version, 1);
  const ScalerArray scaled = readScalers(packet, version, 2);

  const std::uint64_t bunch = packet->getBunchNumber();
  const auto unknownBunch=std::numeric_limits<std::uint64_t>::max();
  m_event.triggerScalerDiscontinuityFlags=(bunch==unknownBunch?8U:0U)|
      (!m_event.physicalEventSequenceValid?16U:0U);
  if(!m_triggerScalerSnapshots.empty()) {
    const auto& previous=m_triggerScalerSnapshots.back();
    if(previous.bunchEnd!=unknownBunch&&bunch!=unknownBunch&&bunch<previous.bunchEnd)
      m_event.triggerScalerDiscontinuityFlags|=2U;
    if(previous.physicalEventEnd>=0&&m_event.physicalEventSequenceValid&&
        m_event.physicalEventSequence<=previous.physicalEventEnd)
      m_event.triggerScalerDiscontinuityFlags|=4U;
    if(previous.valid&&m_event.triggerPacketStatus==0&&
        (anyCounterDecreased(previous.raw,raw)||anyCounterDecreased(previous.live,live)||
         anyCounterDecreased(previous.scaled,scaled)))m_event.triggerScalerDiscontinuityFlags|=1U;
  }

  bool newSnapshot = m_triggerScalerSnapshots.empty();

  if (!newSnapshot)
  {
    const TriggerScalerSnapshotRecord& previous =
        m_triggerScalerSnapshots.back();

    newSnapshot =
        previous.raw != raw ||
        previous.live != live ||
        previous.scaled != scaled ||
        !previous.valid ||
        packet->getStatus() != 0 ||
        bunch < previous.bunchEnd;
  }

  if (newSnapshot)
  {
    TriggerScalerSnapshotRecord snapshot;

    snapshot.sourceId = m_source.sourceId;
    snapshot.snapshotId =
        static_cast<std::int64_t>(m_triggerScalerSnapshots.size());

    snapshot.bunchStart = bunch;
    snapshot.bunchEnd = bunch;

    snapshot.sourceEntryStart = m_event.sourceEntry;
    snapshot.sourceEntryEnd = m_event.sourceEntry;
    snapshot.observedEvents = 1;
    snapshot.physicalEventEnd=m_event.physicalEventSequenceValid?m_event.physicalEventSequence:-1;

    snapshot.raw = raw;
    snapshot.live = live;
    snapshot.scaled = scaled;

    snapshot.valid = m_event.triggerPacketStatus == 0;

    if (!m_triggerScalerSnapshots.empty())
    {
      const TriggerScalerSnapshotRecord& previous =
          m_triggerScalerSnapshots.back();

      snapshot.discontinuity =
          anyCounterDecreased(previous.raw, raw) ||
          anyCounterDecreased(previous.live, live) ||
          anyCounterDecreased(previous.scaled, scaled) ||
          bunch < previous.bunchEnd;
    }

    m_triggerScalerSnapshots.push_back(snapshot);
  }
  else
  {
    TriggerScalerSnapshotRecord& snapshot =
        m_triggerScalerSnapshots.back();

    snapshot.bunchEnd = bunch;
    snapshot.sourceEntryEnd = m_event.sourceEntry;
    snapshot.physicalEventEnd=m_event.physicalEventSequenceValid?m_event.physicalEventSequence:-1;
    ++snapshot.observedEvents;
  }

  m_event.triggerScalerSnapshotId =
      m_triggerScalerSnapshots.back().snapshotId;
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
  auto* info =
      findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");

  const bool supported =
      info && std::string(info->ClassName()) == "TriggerRunInfov1";

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
 * Capture the reconstructed event vertex used by downstream reconstruction.
 *
 * The MBD vertex is preferred. GlobalVertexMap provides the fallback. Truth
 * information is deliberately never substituted for a missing reconstructed
 * vertex.
 */
void PhotonJetTree::captureVertexInformation(PHCompositeNode* topNode)
{
  m_event.recoVertexValid = false;
  m_event.recoVertexSource = 0;

  if (auto* mbdVertices =
          findNode::getClass<MbdVertexMap>(
              topNode,
              m_config.nodes.mbdVertexMap);
      mbdVertices && !mbdVertices->empty())
  {
    if (const MbdVertex* vertex = mbdVertices->begin()->second;
        vertex && finiteValue(vertex->get_z()))
    {
      m_event.recoVertexX = 0.0;
      m_event.recoVertexY = 0.0;
      m_event.recoVertexZ = vertex->get_z();
      m_event.recoVertexValid = true;
      m_event.recoVertexSource = 1;

      return;
    }
  }

  if (auto* vertices =
          findNode::getClass<GlobalVertexMap>(
              topNode,
              m_config.nodes.globalVertexMap);
      vertices && !vertices->empty())
  {
    for (auto it = vertices->begin(); it != vertices->end(); ++it)
    {
      const GlobalVertex* vertex = it->second;

      if (!vertex ||
          !finiteValue(vertex->get_x()) ||
          !finiteValue(vertex->get_y()) ||
          !finiteValue(vertex->get_z()))
      {
        continue;
      }

      m_event.recoVertexX = vertex->get_x();
      m_event.recoVertexY = vertex->get_y();
      m_event.recoVertexZ = vertex->get_z();
      m_event.recoVertexValid = true;
      m_event.recoVertexSource = 2;

      return;
    }
  }
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

  if (MbdOut* out =
          findNode::getClass<MbdOut>(
              topNode,
              m_config.nodes.mbdOut);
      out)
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

    if (finiteValue(m_event.mbdSouthCharge) &&
        finiteValue(m_event.mbdNorthCharge))
    {
      m_event.mbdTotalCharge =
          m_event.mbdSouthCharge + m_event.mbdNorthCharge;
    }

    m_event.mbdValid =
        finiteValue(m_event.mbdT0Ns) &&
        finiteValue(m_event.mbdSouthTimeNs) &&
        finiteValue(m_event.mbdNorthTimeNs) &&
        finiteValue(m_event.mbdSouthCharge) &&
        finiteValue(m_event.mbdNorthCharge);
  }

  MbdPmtContainer* pmts =
      findNode::getClass<MbdPmtContainer>(
          topNode,
          m_config.nodes.mbdPmts);

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

  MbdGeom* geometry =
      findNode::getClass<MbdGeom>(
          topNode,
          m_config.nodes.mbdGeometry);

  m_event.mbdPmtId.reserve(count);
  m_event.mbdPmtArm.reserve(count);
  m_event.mbdPmtCharge.reserve(count);
  m_event.mbdPmtTimeNs.reserve(count);
  m_event.mbdPmtValid.reserve(count);

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

    m_event.mbdPmtValid.push_back(
        hit &&
        hit->get_pmt() == channel &&
        finiteValue(charge)
            ? 1
            : 0);
  }
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
  captureCentralityInputs(m_event,
      findNode::getClass<MbdPmtContainer>(topNode, m_config.nodes.mbdPmts),
      findNode::getClass<MbdOut>(topNode, m_config.nodes.mbdOut),
      isAuAu() ? findNode::getClass<MinimumBiasInfo>(topNode, m_config.nodes.minimumBiasInfo) : nullptr);

  if (!isAuAu())
  {
    return;
  }

  if (isEmbeddedSimulation())
  {
    m_event.embeddedMinimumBias = m_event.minimumBiasDecision;
  }

  auto* centrality =
      findNode::getClass<CentralityInfo>(
          topNode,
          m_config.nodes.centralityInfo);

  if (centrality)
  {
    const double bin = centrality->get_centrality_bin(CentralityInfo::PROP::mbd_NS);
    if (finiteValue(bin) && bin >= 0.0 && bin <= std::numeric_limits<std::int32_t>::max())
      m_event.analysisCentrality = static_cast<std::int32_t>(bin);
  }

  constexpr auto property = CentralityInfo::PROP::mbd_NS;
  const bool writtenForThisEvent =
      isEmbeddedSimulation() ||
      m_event.minimumBiasDecision == TriState::True;
  // CentralityInfov2 does not reset its payload on every event. A rejected
  // DATA event must not inherit a previous event's native centrality witness.
  if (!writtenForThisEvent || !centrality || !centrality->has_centile(property) ||
      !centrality->has_centrality_bin(property)) return;
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
  auto* rawClusters=findNode::getClass<RawClusterContainer>(topNode,m_config.photon.rawClusterNode);
  if(!rawClusters)rawClusters=findNode::getClass<RawClusterContainer>(topNode,m_config.photon.photonNode);
  if(rawClusters) {
    m_event.maximumRawClusterEnergy=0.;
    for(const auto& item:rawClusters->getClustersMap())if(item.second) {
      const double energy=item.second->get_energy();
      if(std::isfinite(energy)&&energy>=1.)
        m_event.maximumRawClusterEnergy=std::max(m_event.maximumRawClusterEnergy,energy);
    }
  }
  const bool requireGood =
      m_config.calorimeter.requireGoodTowers;

  const std::array<CaloSum, 3> sums = {
      sumCalorimeter(
          topNode,
          m_config.nodes.cemcTowers,
          requireGood),

      sumCalorimeter(
          topNode,
          m_config.nodes.ihcalTowers,
          requireGood),

      sumCalorimeter(
          topNode,
          m_config.nodes.ohcalTowers,
          requireGood)};

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
    else allPresent=false;

    if (sum.valid)
    {
      m_event.caloValidMask |= 1U << (8U + layer);
    }
  }

  m_event.cemcEnergy =
      sums[0].present ? static_cast<float>(sums[0].energy) : kNaN;

  m_event.ihcalEnergy =
      sums[1].present ? static_cast<float>(sums[1].energy) : kNaN;

  m_event.ohcalEnergy =
      sums[2].present ? static_cast<float>(sums[2].energy) : kNaN;

  m_event.totalCaloEnergy =
      allPresent ? static_cast<float>(total) : kNaN;

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
    auto* background =
        findNode::getClass<TowerBackground>(
            topNode,
            "TowerInfoBackground_Sub2");

    if (!background)
    {
      throw std::runtime_error(
          "required AuAu TowerInfoBackground_Sub2 is missing");
    }

    m_event.jetBackgroundState = CaptureState::Complete;
    m_event.jetBackgroundFlowMode =
        m_config.calorimeter.backgroundFlowMode;

    m_event.jetBackgroundV2 = background->get_v2();
    m_event.jetBackgroundPsi2 = background->get_Psi2();

    m_event.jetBackgroundUEEmcal = background->get_UE(0);
    m_event.jetBackgroundUEIhcal = background->get_UE(1);
    m_event.jetBackgroundUEOhcal = background->get_UE(2);

    m_event.jetBackgroundNStrips =
        background->get_nStripsUsedForFlow();

    m_event.jetBackgroundNTowers =
        background->get_nTowersUsedForBkg();

    m_event.jetBackgroundFlowFailure =
        background->get_flow_failure_flag()
            ? TriState::True
            : TriState::False;
  }
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

  std::unique_ptr<TFile> file(
      TFile::Open(
          weights.vertexReweightFile.c_str(),
          "READ"));

  if (!file || file->IsZombie())
  {
    throw std::runtime_error(
        "cannot open the vertex-reweighting file " +
        weights.vertexReweightFile);
  }

  auto* histogram =
      dynamic_cast<TH1*>(
          file->Get(
              weights.vertexReweightHistogram.c_str()));

  if (!histogram)
  {
    throw std::runtime_error(
        "vertex-reweighting histogram '" +
        weights.vertexReweightHistogram +
        "' is absent from " +
        weights.vertexReweightFile);
  }

  m_vertexReweight =
      dynamic_cast<TH1*>(
          histogram->Clone(
              "photonjet_vertex_reweight"));

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
  final.vertexWeight = 1.0;
  final.exposureWeight = 1.0;
  if (m_config.system == CollisionSystem::PP)
  {
    final.sliceWeight = final.crossSectionWeight = 1.0;
    final.siDiWeight = final.periodWeight = 1.0;
  }

  if (isData())
  {
    final.finalWeight = 1.0;
    final.valid = true;

    m_weightComponents.push_back(final);

    m_event.eventWeight = 1.0;
    m_event.eventWeightValid = true;

    return;
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
      // Embedded AuAu needs its own qualified one-collision truth-vertex fit.
      auto lookup = [&](double z) {
        if (!finiteValue(z) || z <= -9000.0) return kNaN;
        if (z < m_vertexReweight->GetXaxis()->GetXmin() ||
            z > m_vertexReweight->GetXaxis()->GetXmax()) return 0.0;
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
}
