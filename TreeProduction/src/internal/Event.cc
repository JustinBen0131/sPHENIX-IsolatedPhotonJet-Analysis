//
// Event.cc
//
// Event-level state of the isolated-photon + jet tree producer.
//
// Owns
//   source and event identity, and the exposure accounting
//   trigger decisions, scaler snapshots and the run trigger table
//   the reconstructed vertex
//   the minimum-bias detector, event level and per tube
//   Au+Au centrality as reconstructed, with the trust rule
//   calorimeter event energy and tower census
//   producer event weights
//
// Does not own
//   photons                      -> Photons.cc
//   jets and pairs               -> Jets.cc
//   truth and associations       -> Truth.cc
//   ROOT tables                  -> Output.cc
//   reconstruction itself        -> Production.cc
//
// This file reads what the registered reconstruction put on the node tree.
// It reconstructs nothing and selects nothing.
//
#include "../PhotonJetTree.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

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
#include <cmath>
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

bool finite(const double value)
{
  return std::isfinite(value);
}

// ---- trigger packet --------------------------------------------------------
//
// The packet class version decides which accessors are meaningful. Version 1
// exposes only the legacy decision word and its scalers through the signed
// lValue interface; versions 2 and 3 expose the live and scaled words and
// unsigned 64-bit scalers. In versions 2 and 3 the legacy getter aliases the
// live word; it is not an independently measured raw word.

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

// type 0 raw, 1 live, 2 scaled.
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

// ---- calorimeter sums --------------------------------------------------------
//
// A finite partial sum and full-input validity are deliberately separate. A
// container with one malformed tower still yields an informative sum, with
// valid=false recording that the input was not complete.

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

CaloSum sumCalorimeter(PHCompositeNode* topNode, const std::string& node, const bool requireGood)
{
  CaloSum result;
  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, node);
  if (!towers)
  {
    return result;
  }
  result.present = true;
  result.valid = true;
  result.towerCount = towers->size();

  for (unsigned int channel = 0; channel < towers->size(); ++channel)
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
    if (!finite(energy))
    {
      ++result.nonFiniteCount;
      result.valid = false;
      continue;
    }
    result.energy += static_cast<double>(energy);  // signed, no floor
    ++result.acceptedCount;
  }
  if (result.acceptedCount == 0)
  {
    result.valid = false;
  }
  return result;
}

// Centrality in percent from the reconstructed object. The bin accessor is
// preferred; the centile accessor is stored as a fraction by one object
// version and as a percent by another, and is normalised accordingly.
double centralityPercent(const CentralityInfo* centrality)
{
  constexpr auto property = CentralityInfo::PROP::mbd_NS;
  if (centrality->has_centrality_bin(property))
  {
    const double bin = centrality->get_centrality_bin(property);
    if (finite(bin) && bin >= 0.0 && bin <= 100.0)
    {
      return bin;
    }
  }
  if (centrality->has_centile(property))
  {
    const double centile = centrality->get_centile(property);
    if (finite(centile) && centile >= 0.0 && centile <= 100.0)
    {
      return centile <= 1.0 ? 100.0 * centile : centile;
    }
  }
  return kNaN;
}

int coarseCentralityBin(const double percent, const std::vector<int>& edges)
{
  if (!finite(percent))
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

//
// Per-event reset
//
// Only the state of one event is cleared. Source accounting, scaler
// snapshots, the trigger table, the upstream observation and the output
// handles persist.
//
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

//
// Exposure accounting for events reconstruction stopped
//
void PhotonJetTree::observeUpstreamEvent(const int run, const std::int64_t physicalEventSequence)
{
  // A pending observation means the previous event was seen at the head of
  // the chain and never reached process_event.
  if (m_upstreamPending)
  {
    const std::int64_t sourceEntry =
        m_config.source.firstEntry + static_cast<std::int64_t>(m_upstreamObservedEvents) - 1;
    recordUpstreamRejectedEvent(m_upstreamPendingRun, m_upstreamPendingSequence, sourceEntry);
  }

  m_upstreamPending = true;
  m_upstreamPendingRun = run;
  m_upstreamPendingSequence = physicalEventSequence;
  ++m_upstreamObservedEvents;
}

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

//
// Source and event identity
//
void PhotonJetTree::captureEventIdentity(PHCompositeNode* topNode)
{
  m_event.sourceId = m_source.sourceId;
  m_event.run = m_config.source.run;
  m_event.segment = m_config.source.segment;
  m_event.sourceFileOrdinal = m_config.source.sourceFileOrdinal;

  // ---- physical event header ------------------------------------------------

  EventHeader* header = findNode::getClass<EventHeader>(topNode, m_config.nodes.eventHeader);
  if (header && header->isValid())
  {
    m_event.physicalEventSequence = header->get_EvtSequence();
    m_event.physicalEventSequenceValid = true;

    // The configured run is the production binding and must agree with the
    // physical run whenever both are meaningful.
    const int headerRun = header->get_RunNumber();
    if (headerRun > 0 && m_event.run > 0 && headerRun != m_event.run)
    {
      throw std::runtime_error("event header run " + std::to_string(headerRun) +
                               " disagrees with the configured run " + std::to_string(m_event.run));
    }
  }

  // ---- original source entry --------------------------------------------------
  //
  // Sequential sources: first entry plus this producer's encounter count.
  // Paired or skipped sources: the configured resolver maps the physical
  // sequence back to the original cursor.

  if (m_config.source.sourceEntryResolver)
  {
    m_event.sourceEntry = m_config.source.sourceEntryResolver(
        m_event.run, m_event.physicalEventSequence, m_event.physicalEventSequenceValid);
  }
  else
  {
    m_event.sourceEntry =
        m_config.source.firstEntry + static_cast<std::int64_t>(
            m_upstreamObservedEvents > 0 ? m_upstreamObservedEvents - 1 : m_event.producerEventOrdinal);
  }
  m_event.sourceGlobalEntry = m_event.sourceEntry;

  m_event.eventId = photonjet::makeEventIdentity(m_event.sourceId, m_event.sourceEntry,
                                                 m_event.physicalEventSequence,
                                                 m_event.physicalEventSequenceValid);
}

//
// Trigger decisions and scaler snapshots
//
// Simulation carries no trigger. In data the three decision words are read
// through their own accessors, the packet version and status are recorded,
// and the cumulative scalers are stored as a new snapshot only when a
// counter changed.
//
void PhotonJetTree::captureTriggerInformation(PHCompositeNode* topNode)
{
  m_event.triggerPacketValid = false;
  m_event.triggerDecisionsAvailable = 0;
  m_event.triggerScalerSnapshotId = -1;

  if (isSimulation())
  {
    return;
  }

  Gl1Packet* packet = findNode::getClass<Gl1Packet>(topNode, m_config.nodes.gl1Packet);
  if (!packet)
  {
    // Some data files carry the packet under its numeric name.
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
    return;
  }

  m_event.triggerPacketValid = true;
  m_event.triggerPacketStatus = packet->getStatus();
  m_event.triggerPacketNumber = packet->getPacketNumber();
  m_event.triggerBunchNumber = packet->getBunchNumber();

  m_event.triggerInputBits = packet->getTriggerVector();
  m_event.triggerDecisionsAvailable |= kLegacyWordAvailable;
  if (version >= 2)
  {
    m_event.triggerLiveBits = packet->getLiveVector();
    m_event.triggerScaledBits = packet->getScaledVector();
    m_event.triggerDecisionsAvailable |= kLiveWordAvailable | kScaledWordAvailable;
  }

  // ---- cumulative scalers -----------------------------------------------------

  const ScalerArray raw = readScalers(packet, version, 0);
  const ScalerArray live = readScalers(packet, version, 1);
  const ScalerArray scaled = readScalers(packet, version, 2);
  const std::uint64_t bunch = packet->getBunchNumber();

  bool newSnapshot = m_triggerScalerSnapshots.empty();
  if (!newSnapshot)
  {
    const TriggerScalerSnapshotRecord& previous = m_triggerScalerSnapshots.back();
    newSnapshot = previous.raw != raw || previous.live != live || previous.scaled != scaled ||
                  !previous.valid || packet->getStatus() != 0 || bunch < previous.bunchEnd;
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
    snapshot.raw = raw;
    snapshot.live = live;
    snapshot.scaled = scaled;
    snapshot.valid = m_event.triggerPacketStatus == 0;

    if (!m_triggerScalerSnapshots.empty())
    {
      const TriggerScalerSnapshotRecord& previous = m_triggerScalerSnapshots.back();
      snapshot.discontinuity = anyCounterDecreased(previous.raw, raw) ||
                               anyCounterDecreased(previous.live, live) ||
                               anyCounterDecreased(previous.scaled, scaled) ||
                               bunch < previous.bunchEnd;
    }
    m_triggerScalerSnapshots.push_back(snapshot);
  }
  else
  {
    TriggerScalerSnapshotRecord& snapshot = m_triggerScalerSnapshots.back();
    snapshot.bunchEnd = bunch;
    snapshot.sourceEntryEnd = m_event.sourceEntry;
    ++snapshot.observedEvents;
  }

  m_event.triggerScalerSnapshotId = m_triggerScalerSnapshots.back().snapshotId;
}

//
// Run trigger configuration
//
// One row per trigger bit, read once per run. A missing table or an
// unsupported class leaves every row invalid rather than absent, so the
// reader can distinguish "no configuration" from "no bits".
//
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
    if (supported)
    {
      row.name = info->getTriggerName(bit);
      row.hardwarePrescale = info->getInitialPrescaleByBit(bit);
      row.runAveragePrescale = info->getPrescaleByBit(bit);
      row.rawCount = info->getRawScalersByBit(bit);
      row.liveCount = info->getLiveScalersByBit(bit);
      row.scaledCount = info->getScalersByBit(bit);
      row.valid = !row.name.empty() && row.name != "unknown";
    }
    m_triggerRunInfo.push_back(std::move(row));
  }
}

//
// Reconstructed vertex
//
// The minimum-bias detector vertex is preferred; the global vertex map is
// the fallback. Both are reconstructed quantities; the truth vertex is never
// substituted for a reconstructed object.
//
void PhotonJetTree::captureVertexInformation(PHCompositeNode* topNode)
{
  m_event.recoVertexValid = false;
  m_event.recoVertexSource = 0;

  if (auto* mbdVertices = findNode::getClass<MbdVertexMap>(topNode, m_config.nodes.mbdVertexMap);
      mbdVertices && !mbdVertices->empty())
  {
    if (const MbdVertex* vertex = mbdVertices->begin()->second; vertex && finite(vertex->get_z()))
    {
      m_event.recoVertexX = 0.0;
      m_event.recoVertexY = 0.0;
      m_event.recoVertexZ = vertex->get_z();
      m_event.recoVertexValid = true;
      m_event.recoVertexSource = 1;
      return;
    }
  }

  if (auto* vertices = findNode::getClass<GlobalVertexMap>(topNode, m_config.nodes.globalVertexMap);
      vertices && !vertices->empty())
  {
    for (auto it = vertices->begin(); it != vertices->end(); ++it)
    {
      const GlobalVertex* vertex = it->second;
      if (!vertex || !finite(vertex->get_x()) || !finite(vertex->get_y()) || !finite(vertex->get_z()))
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

//
// Minimum-bias detector
//
void PhotonJetTree::captureMbdInformation(PHCompositeNode* topNode)
{
  m_event.mbdValid = false;
  m_event.mbdPmtAvailable = TriState::Unknown;

  if (MbdOut* out = findNode::getClass<MbdOut>(topNode, m_config.nodes.mbdOut); out && out->isValid())
  {
    m_event.mbdT0Ns = out->get_t0();
    m_event.mbdSouthTimeNs = out->get_time(0);
    m_event.mbdNorthTimeNs = out->get_time(1);
    m_event.mbdSouthCharge = out->get_q(0);
    m_event.mbdNorthCharge = out->get_q(1);
    if (finite(m_event.mbdSouthCharge) && finite(m_event.mbdNorthCharge))
    {
      m_event.mbdTotalCharge = m_event.mbdSouthCharge + m_event.mbdNorthCharge;
    }
    m_event.mbdValid = finite(m_event.mbdT0Ns) && finite(m_event.mbdSouthTimeNs) &&
                       finite(m_event.mbdNorthTimeNs) && finite(m_event.mbdSouthCharge) &&
                       finite(m_event.mbdNorthCharge);
  }

  MbdPmtContainer* pmts = findNode::getClass<MbdPmtContainer>(topNode, m_config.nodes.mbdPmts);
  if (!pmts)
  {
    m_event.mbdPmtAvailable = TriState::False;
    return;
  }
  const int count = pmts->get_npmt();
  if (count < 0 || count > 128)
  {
    // A malformed container is neither present nor absent.
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

  for (int channel = 0; channel < count; ++channel)
  {
    const MbdPmtHit* hit = pmts->get_pmt(channel);
    const double charge = hit ? hit->get_q() : kNaN;
    const double time = hit ? hit->get_time() : kNaN;
    const int arm = geometry ? geometry->get_arm(channel) : -1;

    // Container channel, never compacted: holes stay visible.
    m_event.mbdPmtId.push_back(channel);
    m_event.mbdPmtArm.push_back(arm == 0 || arm == 1 ? arm : -1);
    m_event.mbdPmtCharge.push_back(charge);  // zero and negative kept
    m_event.mbdPmtTimeNs.push_back(time);
    m_event.mbdPmtValid.push_back(hit && hit->get_pmt() == channel && finite(charge) && finite(time) ? 1 : 0);
  }
}

//
// Centrality
//
// The reconstructed value is trusted only when the reconstruction wrote it
// for this event: in data that means the minimum-bias decision was positive,
// in embedded simulation the value came with the input. Otherwise a stale
// object from a previous event could pass for a measurement.
//
void PhotonJetTree::captureCentrality(PHCompositeNode* topNode)
{
  m_event.centralityApplicable = isAuAu();
  m_event.minimumBiasDecision = TriState::Unknown;
  if (!isAuAu())
  {
    return;
  }

  if (auto* minimumBias = findNode::getClass<MinimumBiasInfo>(topNode, m_config.nodes.minimumBiasInfo);
      minimumBias && minimumBias->isValid())
  {
    m_event.minimumBiasDecision = minimumBias->isAuAuMinimumBias() ? TriState::True : TriState::False;
  }

  if (isEmbeddedSimulation())
  {
    m_event.embeddedMinimumBias = m_event.minimumBiasDecision;
  }

  auto* centrality = findNode::getClass<CentralityInfo>(topNode, m_config.nodes.centralityInfo);
  if (!centrality || !centrality->isValid())
  {
    return;
  }

  constexpr auto property = CentralityInfo::PROP::mbd_NS;
  if (centrality->has_quantity(property))
  {
    const double charge = centrality->get_quantity(property);
    if (finite(charge))
    {
      m_event.centralitySelectedCharge = charge;
    }
  }

  const double percent = centralityPercent(centrality);
  if (!finite(percent))
  {
    return;
  }

  const bool writtenForThisEvent =
      isEmbeddedSimulation() || m_event.minimumBiasDecision == TriState::True;

  m_event.centralityNativeBin = static_cast<std::int32_t>(std::lround(percent));
  m_event.centralityNativePercent = percent;
  m_event.centralityNativeValid = writtenForThisEvent;

  if (writtenForThisEvent)
  {
    m_event.centralityPercent = percent;
    m_event.centralityBin = coarseCentralityBin(percent, m_config.centrality.edges);
    m_event.centralityValid = true;
  }
}

//
// Calorimeter event energy and tower census
//
void PhotonJetTree::captureCalorimeterInformation(PHCompositeNode* topNode)
{
  const bool requireGood = m_config.calorimeter.requireGoodTowers;
  const std::array<CaloSum, 3> sums = {
      sumCalorimeter(topNode, m_config.nodes.cemcTowers, requireGood),
      sumCalorimeter(topNode, m_config.nodes.ihcalTowers, requireGood),
      sumCalorimeter(topNode, m_config.nodes.ohcalTowers, requireGood)};

  m_event.caloRequiredGoodTowers = requireGood;
  m_event.caloAvailableMask = 0;
  m_event.caloValidMask = 0;

  double total = 0.0;
  bool anyPresent = false;

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
      anyPresent = true;
      total += sum.energy;
    }
    if (sum.valid)
    {
      m_event.caloValidMask |= 1U << (8U + layer);
    }
  }

  m_event.cemcEnergy = sums[0].present ? sums[0].energy : kNaN;
  m_event.ihcalEnergy = sums[1].present ? sums[1].energy : kNaN;
  m_event.ohcalEnergy = sums[2].present ? sums[2].energy : kNaN;
  m_event.totalCaloEnergy = anyPresent ? total : kNaN;

  if (isAuAu())
  {
    // Frozen producer capture: RecoilJets_AuAu.cc:7712-7735.
    auto* background = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
    if (!background) throw std::runtime_error("required AuAu TowerInfoBackground_Sub2 is missing");
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

//
// Producer event weights
//
// Data carries exactly one. Simulation carries the vertex factor when the
// profile configured one, and no other producer factor: slice cross
// sections, stitching and interaction mixing are combined once downstream
// from the sample manifest. A factor that was not determined is NaN with
// valid false; it is never assumed to be one.
//
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
    throw std::runtime_error("vertex-reweighting histogram '" + weights.vertexReweightHistogram +
                             "' is absent from " + weights.vertexReweightFile);
  }
  m_vertexReweight = dynamic_cast<TH1*>(histogram->Clone("photonjet_vertex_reweight"));
  m_vertexReweight->SetDirectory(nullptr);
}

void PhotonJetTree::captureEventWeights(PHCompositeNode* /*topNode*/)
{
  WeightComponentRecord final;
  final.eventId = m_event.eventId;
  final.componentType = "final";
  final.applicationCount = 1;

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
      const int bin = m_vertexReweight->FindFixBin(m_event.truthHardVertexZ);
      const double factor = m_vertexReweight->GetBinContent(bin);
      if (finite(factor) && factor >= 0.0 && bin >= 1 && bin <= m_vertexReweight->GetNbinsX())
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

  final.finalWeight = valid ? product : kNaN;
  final.valid = valid;
  m_weightComponents.push_back(final);

  m_event.eventWeight = final.finalWeight;
  m_event.eventWeightValid = valid;
}
