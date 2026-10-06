/**
 * @file Objects.cc
 * Capture reconstructed photons, jets and their event-local relationships.
 * Photon values come from PhotonClusterBuilder; optional cells and constituents
 * preserve reconstruction witnesses without implementing a second algorithm.
 * Jet capture preserves collection order, radius and raw/calibrated identity.
 * Physics selections remain independent of the broad capture domain.
 */

#include "../PhotonJetTree.h"
#include "../PhotonClusterBuilder.h"
#include "Production.h"
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <algorithm>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TVector2.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>
#include <vector>

// Photon capture and optional reconstruction witnesses

namespace
{
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

constexpr int kCemcEtaBins = 96;
constexpr int kCemcPhiBins = 256;
constexpr int kGridHalfWidth = 3;  // 7x7 grid about the seed tower

bool finiteValue(const double value) { return std::isfinite(value); }


// ============================================================================
// CEMC coordinate helpers
// ============================================================================

/**
 * Wrap a CEMC phi-bin index onto the detector's periodic [0, 255] range.
 */
int wrapPhiIndex(int index)
{
  while (index < 0) index += kCemcPhiBins;
  while (index >= kCemcPhiBins) index -= kCemcPhiBins;

  return index;
}


// ============================================================================
// PhotonClusterBuilder-attached observables
//
// PhotonClusterBuilder stores several reconstruction products in the cluster's
// named shower-shape map. Missing entries remain NaN/absent; this layer does not
// manufacture defaults for unavailable reconstruction output.
// ============================================================================

double namedShape(const RawCluster* cluster, const std::string& name)
{
  const auto& shapes = cluster->get_all_shower_shapes();
  const auto found = shapes.find(name);

  return found == shapes.end() ? kNaN : static_cast<double>(found->second);
}

bool hasNamedShape(const RawCluster* cluster, const std::string& name)
{
  const auto& shapes = cluster->get_all_shower_shapes();

  return shapes.find(name) != shapes.end();
}


/**
 * Convert one builder model invocation into the retained witness record.
 *
 * The builder reports what it evaluated; the object identity (container key
 * and encounter ordinal) comes from this producer's own traversal. The object
 * id and the container key are distinct and both retained. Kinematic input
 * bits are available only for the frozen 25-feature NPB input order.
 */
photonjet::ModelInvocationRecord modelInvocationWitness(const photonjet::PhotonClusterBuilder::ModelInvocation& seen,
                                                        const RawCluster* cluster, const std::uint64_t clusterKey,
                                                        const std::uint32_t ordinal, const std::string& photonNode)
{
  photonjet::ModelInvocationRecord witness;
  witness.captureVersion = 1;
  witness.evaluationState = seen.outcome;
  witness.domainMode = seen.bounded ? 1 : 2;
  witness.inProducerDomain = seen.in_domain ? 1 : 0;
  witness.modelAvailable = seen.model_available ? 1 : 0;
  witness.minEt = seen.min_et;
  witness.maxEt = seen.max_et;
  witness.maxAbsEta = seen.max_abs_eta;
  witness.rawScore = seen.score;
  witness.modelFile = seen.model_file;
  witness.featureNames = seen.features;
  witness.inputOrderIsNpb25 = seen.features == photonjet::npbFeatureNames() ? 1 : 0;

  witness.inputsCaptured = !seen.features.empty() && seen.inputs.size() == seen.features.size() ? 1 : 0;
  if (witness.inputsCaptured)
  {
    witness.orderedInputs = seen.inputs;
    for (const float value : seen.inputs) witness.orderedInputFinite.push_back(std::isfinite(value) ? 1 : 0);
  }
  witness.inputsFinite =
      witness.inputsCaptured && std::all_of(witness.orderedInputFinite.begin(), witness.orderedInputFinite.end(),
                                            [](std::int32_t value) { return value == 1; })
          ? 1
          : 0;

  witness.sourcePhotonId = static_cast<std::uint32_t>(cluster->get_id());
  witness.sourceClusterEt = static_cast<float>(namedShape(cluster, "cluster_pt"));
  witness.sourceEta = static_cast<float>(namedShape(cluster, "cluster_eta"));
  witness.sourcePhi = static_cast<float>(namedShape(cluster, "cluster_phi"));
  witness.sourceIdentityAvailable = 1;
  witness.bindingPolicy = 1;  // the same object carries the score
  witness.sourceNativeClusterKey = static_cast<std::uint32_t>(clusterKey);
  witness.sourceEncounterOrdinal = ordinal;
  witness.sourceClusterNode = photonNode;

  if (witness.inputOrderIsNpb25 && witness.inputsCaptured)
  {
    witness.kinematicInputBitsAvailable = 1;
    std::memcpy(&witness.inputEtBits, &witness.orderedInputs[0], sizeof(float));
    std::memcpy(&witness.inputEtaBits, &witness.orderedInputs[1], sizeof(float));
  }
  return witness;
}


/**
 * The historical reconstructed isolation witness at each retained radius.
 *
 * The cone is the system's nominal method: pp topological clusters, Au+Au
 * SUB1 towers. pp thresholds are lines in photon ET. Au+Au thresholds are
 * lines in the truncated integer centrality percent inside the configured
 * range; outside it, or with a nonfinite centrality, the photon-ET fallback
 * line applies. Missing AuAu cones preserve the legacy 1e9 fail-safe
 * (non-isolated); pp retains its NaN invalid witness.
 */
std::vector<photonjet::IsolationWitnessRecord> isolationWitnesses(const photonjet::PhotonRecord& photon,
                                                                  const double centrality,
                                                                  const photonjet::PhotonConfig& config,
                                                                  const bool auau)
{
  using photonjet::IsolationMethod;
  const auto& witness = config.isolationWitness;
  const IsolationMethod method = auau ? IsolationMethod::CalorimeterSub1 : IsolationMethod::Topocluster;
  const int centralityBin = auau && std::isfinite(centrality) ? static_cast<int>(centrality) : -1;
  const bool routed =
      centralityBin >= 0 && centralityBin >= witness.centralityMin && centralityBin < witness.centralityMax;

  std::vector<photonjet::IsolationWitnessRecord> rows;
  for (std::size_t index = 0; index < config.isolationRadii.size(); ++index)
  {
    photonjet::IsolationWitnessRecord row;
    row.radius = config.isolationRadii[index];
    row.subtractionMethod = auau ? 2 : 1;
    row.reconstructedOrTruth = 0;
    if (auau) row.coneSum = 1e9;
    for (const photonjet::IsolationRecord& cone : photon.isolation)
      if (cone.method == method && cone.radius == row.radius && cone.valid) row.coneSum = cone.coneSum;

    const photonjet::IsolationLine& line = auau && !routed ? witness.fallback : witness.lines.at(index);
    const double x = auau && routed ? static_cast<double>(centralityBin) : photon.et;
    row.threshold = line.interceptGeV + line.slope * x;
    row.sidebandThreshold = row.threshold + line.sidebandGapGeV;
    row.passState = std::isfinite(row.coneSum) && row.coneSum < row.threshold ? 1 : 0;
    rows.push_back(row);
  }
  return rows;
}


// ============================================================================
// CEMC cell witness
//
// Cell is the temporary in-memory representation of one calorimeter location in
// the union of:
//   - towers owned by the reconstructed cluster; and
//   - towers entering any configured 7x7 shower-support grid.
//
// It keeps the calibrated tower state and the cluster-map energy separately so
// later output can distinguish detector information from cluster ownership.
// ============================================================================

struct Cell
{
  std::uint64_t towerKey = 0;

  int etaIndex = -1;
  int phiIndex = -1;

  bool present = false;
  bool good = false;
  bool saturated = false;

  double energy = kNaN;
  double timeSamples = kNaN;

  std::uint32_t status = 0;

  bool owned = false;
  double clusterMapEnergy = kNaN;

  bool inSupport = false;
};

using CellMap = std::map<std::pair<int, int>, Cell>;


/**
 * Ensure one eta/phi location exists in the temporary cell map and, when the
 * calibrated tower exists, attach its detector state.
 *
 * Phi is periodic; eta outside the physical CEMC range is ignored.
 */
void addCell(CellMap& cells, TowerInfoContainer* towers, const int etaIndex, const int phiIndex)
{
  if (etaIndex < 0 || etaIndex >= kCemcEtaBins)
  {
    return;
  }

  const int phi = wrapPhiIndex(phiIndex);

  Cell& cell = cells[{etaIndex, phi}];

  if (cell.towerKey != 0)
  {
    return;
  }

  const unsigned int key = TowerInfoDefs::encode_emcal(etaIndex, phi);

  cell.towerKey = key;
  cell.etaIndex = etaIndex;
  cell.phiIndex = phi;

  if (TowerInfo* tower = towers->get_tower_at_channel(TowerInfoDefs::decode_emcal(key)))
  {
    cell.present = true;
    cell.good = tower->get_isGood();
    cell.saturated = tower->get_isSaturated();

    cell.energy = tower->get_energy();
    cell.timeSamples = tower->get_time();
    cell.status = tower->get_status();
  }
}

}  // namespace


// ============================================================================
// PPG12 photon-yield source arbitration (pp)
//
// PPG12 counts photon yield from the leading cluster that is not on a masked
// tower. Every cluster of the photon node is considered, in encounter order:
// the mask is read at the cluster's energy-weighted tower centroid
// (cluster_ietacent/iphicent, the builder's shower-shape centroid, the
// position PPG12 used), and the leading source is the unmasked cluster with
// the highest cluster_pt of at least 5 GeV.
// ============================================================================

bool photonjet::production::PPG12TowerMasked(const TH2& mask, const double ietaCentroid, const double iphiCentroid)
{
  // Truncate toward zero exactly as PPG12's integer conversion, range-checked
  // before the conversion.
  const double ieta = std::trunc(ietaCentroid), iphi = std::trunc(iphiCentroid);
  if (!std::isfinite(ieta) || !std::isfinite(iphi) || ieta < 0 || iphi < 0 || ieta >= mask.GetNbinsX() ||
      iphi >= mask.GetNbinsY())
    return false;
  return mask.GetBinContent(static_cast<int>(ieta) + 1, static_cast<int>(iphi) + 1) > 0.0;
}

void photonjet::production::PPG12ArbitrateSource(EventRecord& event, const std::uint32_t clusterKey,
                                                 const std::int32_t encounter, const bool masked, const double et,
                                                 const double eta, const double phi)
{
  constexpr double minimumSourceEtGeV = 5.0;
  if (masked || !std::isfinite(et) || et < minimumSourceEtGeV) return;
  if (std::isfinite(event.ppg12LeadEt) && et <= event.ppg12LeadEt) return;
  event.ppg12LeadKey = clusterKey;
  event.ppg12LeadEncounter = encounter;
  event.ppg12LeadEt = et;
  event.ppg12LeadEta = eta;
  event.ppg12LeadPhi = phi;
}

/**
 * Load the pinned PPG12 tower mask once (pp). Its digest was verified with
 * the other pinned inputs before reconstruction.
 */
void PhotonJetTree::loadPPG12TowerMask()
{
  const PhotonConfig& photon = m_config.photon;
  if (photon.ppg12TowerMaskFile.empty()) return;

  TDirectory::TContext context;
  std::unique_ptr<TFile> file(TFile::Open(photon.ppg12TowerMaskFile.c_str(), "READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error("cannot open the PPG12 tower mask " + photon.ppg12TowerMaskFile);
  auto* mask = dynamic_cast<TH2*>(file->Get(photon.ppg12TowerMaskHistogram.c_str()));
  if (!mask)
    throw std::runtime_error("PPG12 tower mask '" + photon.ppg12TowerMaskHistogram + "' is absent from " +
                             photon.ppg12TowerMaskFile);
  m_ppg12TowerMask = static_cast<TH2*>(mask->Clone("photonjet_ppg12_tower_mask"));
  m_ppg12TowerMask->SetDirectory(nullptr);
}


// ============================================================================
// Reconstructed photon candidates
// ============================================================================

/**
 * Capture the complete reconstructed-photon inventory for the current event.
 *
 * Candidate admission depends only on the broad storage domain:
 *
 *   valid reconstructed event vertex
 *   configured object-vertex domain
 *   finite photon kinematics
 *   configured ET range
 *   configured eta acceptance
 *
 * No BDT score, working point, isolation requirement, or analysis-region
 * decision is applied here.
 */
void PhotonJetTree::capturePhotons(PHCompositeNode* topNode)
{
  RawClusterContainer* candidates = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.photonNode);

  if (!candidates)
  {
    throw std::runtime_error("missing photon candidate node '" + m_config.photon.photonNode + "'");
  }

  TowerInfoContainer* cemcTowers = findNode::getClass<TowerInfoContainer>(topNode, m_config.nodes.cemcTowers);

  if (!cemcTowers)
  {
    throw std::runtime_error("missing calibrated CEMC tower node '" + m_config.nodes.cemcTowers + "'");
  }

  // The workfest DATA contract limits reconstructed objects to |z| < 60 cm.
  // Its SIM population is unfiltered by that DATA storage cut. Photons use
  // the builder's vertex (notably the pp-SIM MBD component).
  m_event.recoObjectVertexInDomain =
      m_event.recoVertexValid &&
      (isSimulation() || std::abs(m_event.recoVertexZ) < m_config.photon.objectVertexAbsZMaxCm);

  if (!m_photonBuilder) throw std::runtime_error("photon capture requires the registered builder's event witness");

  const double photonVertexZ = m_photonBuilder->get_vertex_z();
  if (!finiteValue(photonVertexZ))
  {
    m_event.recoPhotonCaptureState = CaptureState::InputUnavailable;
    return;
  }
  if (isData() && !(std::abs(photonVertexZ) < m_config.photon.objectVertexAbsZMaxCm))
  {
    m_event.recoPhotonCaptureState = CaptureState::NotApplicable;
    return;
  }
  m_event.recoPhotonCaptureState = CaptureState::Complete;
  m_event.recoPhotonCaptureVersion = 1;
  // 1: ET 5-40 GeV, 2: ET 12-40 GeV (Au+Au DATA), both |eta| < 0.7; -1 otherwise.
  const bool standardEdges = m_config.photon.maxEtGeV == 40. && m_config.photon.maxAbsEta == .7;
  m_event.recoPhotonCaptureDomain = !standardEdges                    ? -1
                                    : m_config.photon.minEtGeV == 5.  ? 1
                                    : m_config.photon.minEtGeV == 12. ? 2
                                                                      : -1;

  std::uint32_t encountered = 0;

  // The PPG12 source arbitration considers every cluster, retained or not.
  const TH2* ppg12Mask = m_ppg12TowerMask;
  if (ppg12Mask) m_event.ppg12ArbitrationState = 1;

  for (const auto& [clusterKey, cluster] : candidates->getClustersMap())
  {
    const std::uint32_t ordinal = encountered++;
    ++m_event.recoPhotonInspectedCount;

    if (!cluster)
    {
      ++m_event.recoPhotonUnclassifiableCount;
      continue;
    }

    const bool ppg12Masked =
        ppg12Mask && photonjet::production::PPG12TowerMasked(*ppg12Mask, namedShape(cluster, "cluster_ietacent"),
                                                             namedShape(cluster, "cluster_iphicent"));

    if (ppg12Mask)
    {
      photonjet::production::PPG12ArbitrateSource(
          m_event, static_cast<std::uint32_t>(clusterKey), static_cast<std::int32_t>(ordinal), ppg12Masked,
          namedShape(cluster, "cluster_pt"), namedShape(cluster, "cluster_eta"), namedShape(cluster, "cluster_phi"));
    }

    PhotonRecord photon;
    photon.ppg12TowerMaskState = ppg12Mask ? (ppg12Masked ? 1 : 0) : -1;

    photon.eventId = m_event.eventId;
    photon.nativeClusterKey = clusterKey;
    photon.producerVertexZ = namedShape(cluster, "vertex_z");
    photon.clusterProbability = cluster->get_prob();
    photon.energy = cluster->get_energy();
    photon.et = namedShape(cluster, "cluster_pt");
    photon.eta = namedShape(cluster, "cluster_eta");
    photon.phi = namedShape(cluster, "cluster_phi");
    for (std::size_t index = 0; index < PhotonRecord::primaryShapeNames.size(); ++index)
      photon.primaryShapes[index] = namedShape(cluster, PhotonRecord::primaryShapeNames[index]);
    // Preselection-score invocations observed by the builder on this object.
    for (const auto& [scoreName, witness] :
         {std::pair<const char*, photonjet::ModelInvocationRecord*>{"npb_score", &photon.npbCapture},
          std::pair<const char*, photonjet::ModelInvocationRecord*>{"auau_npb_score", &photon.auauNpbCapture}})
    {
      if (const auto* seen = m_photonBuilder->get_model_invocation(cluster, scoreName))
      {
        *witness = modelInvocationWitness(*seen, cluster, clusterKey, ordinal, m_config.photon.photonNode);
      }
    }

    if (!finiteValue(photon.producerVertexZ) || photon.producerVertexZ != photonVertexZ)
      throw std::runtime_error("photon vertex does not match the current builder event");

    photon.kinematicsFinite =
        finiteValue(photon.energy) && finiteValue(photon.et) && finiteValue(photon.eta) && finiteValue(photon.phi);

    // The PPG12 simulation bookkeeping reads every candidate, stored or not.
    if (recordsPPG12Simulation())
    {
      PhotonRecord candidate = photon;
      candidate.encounterOrdinal = ordinal;
      candidate.isolation = photonIsolation(cluster, candidate);
      if (!m_config.photon.productionScore.scoreName.empty())
        candidate.productionScore = namedShape(cluster, m_config.photon.productionScore.scoreName);
      m_ppg12Photons.push_back(std::move(candidate));
    }

    if (!photonInCaptureDomain(photon))
    {
      if (!photon.kinematicsFinite)
      {
        ++m_event.recoPhotonUnclassifiableCount;
      }

      continue;
    }

    ++m_event.recoPhotonEligibleCount;
    /*
     * encounterOrdinal records where this candidate appeared in the upstream
     * candidate container. photonId combines event identity, native cluster key,
     * and that ordinal into the stable relational key used downstream.
     */
    photon.encounterOrdinal = ordinal;

    photon.photonId = photonjet::makePhotonIdentity(m_event.eventId, clusterKey, ordinal);

    capturePhotonTiming(cluster, photon);

    capturePhotonShowerShapes(cluster, cemcTowers, photon);

    // PhotonRecord carries the isolation views while the normalized Isolation
    // table receives the same records for direct relational access.
    photon.isolation = photonIsolation(cluster, photon);
    m_isolationRecords.insert(m_isolationRecords.end(), photon.isolation.begin(), photon.isolation.end());

    photon.isolationWitnesses = isolationWitnesses(photon, m_event.analysisCentrality, m_config.photon, isAuAu());

    // The configured producer score exactly as the builder wrote it.
    if (!m_config.photon.productionScore.scoreName.empty())
      photon.productionScore = namedShape(cluster, m_config.photon.productionScore.scoreName);

    capturePhotonIsolationConstituents(topNode, cluster, photon);

    /*
     * Truth.cc owns the dominant-primary calculation and final photon-truth
     * relation. Until that stage runs, simulation explicitly records the
     * evaluator as unavailable; DATA records truth as not applicable.
     */
    photon.dominantTruth.state =
        isSimulation() ? DominantTruthState::EvaluatorUnavailable : DominantTruthState::NotSimulation;

    m_photons.push_back(std::move(photon));
    ++m_event.recoPhotonWrittenCount;
  }
  if (m_event.recoPhotonInspectedCount != candidates->getClustersMap().size() ||
      m_event.recoPhotonWrittenCount != m_event.recoPhotonEligibleCount || m_event.recoPhotonUnclassifiableCount != 0)
    m_event.recoPhotonCaptureState = CaptureState::InputUnavailable;
}


// ============================================================================
// Photon storage domain
// ============================================================================

/**
 * Apply only the broad reconstruction-storage acceptance.
 *
 * This deliberately differs from a physics selection. Any photon inside this
 * domain is preserved regardless of identification or isolation so downstream
 * working points can change without regenerating the base tree.
 */
bool PhotonJetTree::photonInCaptureDomain(const PhotonRecord& photon) const
{
  if (!photon.kinematicsFinite)
  {
    return false;
  }

  if (photon.et < m_config.photon.minEtGeV || photon.et >= m_config.photon.maxEtGeV)
  {
    return false;
  }

  if (std::abs(photon.eta) >= m_config.photon.maxAbsEta)
  {
    return false;
  }

  return true;
}


// ============================================================================
// Builder-produced timing and shower views
//
// PhotonClusterBuilder owns all nominal timing and shower arithmetic. This
// capture layer only copies its named outputs and records optional calibrated
// CEMC cells as low-level witnesses. The witness census does not feed a second
// shower calculation.
// ============================================================================

void PhotonJetTree::capturePhotonTiming(const RawCluster* cluster, PhotonRecord& photon) const
{
  PhotonTimingRecord& timing = photon.timing;
  timing.energyWeightedNumerator = namedShape(cluster, "time_energy_numerator");
  timing.energyDenominator = namedShape(cluster, "time_energy_denominator");
  timing.meanTimeSamples = namedShape(cluster, "mean_time");

  const double count = namedShape(cluster, "time_contributing_towers");
  const bool countValid = finiteValue(count) && count >= 0.0 && count <= std::numeric_limits<std::uint32_t>::max() &&
                          std::floor(count) == count;
  if (countValid)
  {
    timing.contributingTowerCount = static_cast<std::uint32_t>(count);
  }

  timing.finite = finiteValue(timing.meanTimeSamples);
  timing.valid = namedShape(cluster, "time_valid") > 0.5 && countValid && finiteValue(timing.energyWeightedNumerator) &&
                 finiteValue(timing.energyDenominator) && timing.energyDenominator > 0.0 && timing.finite;
  timing.timeNs = timing.valid ? timing.meanTimeSamples * m_config.photon.timingSampleNs : kNaN;
}

void PhotonJetTree::capturePhotonShowerShapes(const RawCluster* cluster, TowerInfoContainer* cemcTowers,
                                              PhotonRecord& photon)
{
  CellMap cells;
  for (const auto& [rawKey, mapEnergy] : cluster->get_towermap())
  {
    const int etaIndex = RawTowerDefs::decode_index1(rawKey);
    const int phiIndex = wrapPhiIndex(RawTowerDefs::decode_index2(rawKey));
    addCell(cells, cemcTowers, etaIndex, phiIndex);
    if (auto found = cells.find({etaIndex, phiIndex}); found != cells.end())
    {
      found->second.owned = true;
      found->second.clusterMapEnergy = mapEnergy;
    }
  }

  const double primaryFloor = namedShape(cluster, "shower_shape_floor_gev");
  for (const std::string& definitionName : m_config.showerDefinitions)
  {
    const double floor = photonjet::showerViewDefinition(definitionName).floorGeV;
    // The builder stores binary32 floors. Match that exact representation:
    // nearby thresholds are different views, not approximate H70 aliases.
    const bool primary = primaryFloor == static_cast<double>(static_cast<float>(floor));
    const std::string prefix = primary ? "" : "shower_" + definitionName + "_";
    const auto value = [&](const std::string& key) { return namedShape(cluster, prefix + key); };

    ShowerShapeRecord shape;
    shape.definitionName = definitionName;
    shape.energyFloorGeV = floor;
    shape.et1 = value("et1");
    shape.et2 = value("et2");
    shape.et3 = value("et3");
    shape.et4 = value("et4");
    shape.e11 = value("e11");
    shape.e13 = value("e13");
    shape.e15 = value("e15");
    shape.e17 = value("e17");
    shape.e22 = value("e22");
    shape.e31 = value("e31");
    shape.e32 = value("e32");
    shape.e33 = value("e33");
    shape.e35 = value("e35");
    shape.e37 = value("e37");
    shape.e51 = value("e51");
    shape.e52 = value("e52");
    shape.e53 = value("e53");
    shape.e55 = value("e55");
    shape.e57 = value("e57");
    shape.e71 = value("e71");
    shape.e72 = value("e72");
    shape.e73 = value("e73");
    shape.e75 = value("e75");
    shape.e77 = value("e77");
    shape.weta = value("weta");
    shape.wphi = value("wphi");
    shape.wetaCog = value("weta_cog");
    shape.wphiCog = value("wphi_cog");
    shape.wetaCogX = value("weta_cogx");
    shape.wphiCogX = value("wphi_cogx");
    shape.weta33CogX = value("weta33_cogx");
    shape.wphi33CogX = value("wphi33_cogx");
    shape.w32 = value("w32");
    shape.w52 = value("w52");
    shape.w72 = value("w72");
    shape.centerOfGravityEta = value("cog_eta_local");
    shape.centerOfGravityPhi = value("cog_phi_local");
    shape.deltaEtaCog = value("detacog");
    shape.deltaPhiCog = value("dphicog");
    shape.radialSpread = value("drad");

    const double centerEta = value("center_ieta");
    const double centerPhi = value("center_iphi");
    if (finiteValue(centerEta) && centerEta >= 0.0 && centerEta < kCemcEtaBins && finiteValue(centerPhi) &&
        centerPhi >= 0.0 && centerPhi < kCemcPhiBins)
    {
      shape.centerEtaIndex = static_cast<std::int32_t>(centerEta);
      shape.centerPhiIndex = static_cast<std::int32_t>(centerPhi);
    }
    const double deltaEtaMax = value("detamax");
    const double deltaPhiMax = value("dphimax");
    const double saturated = value("nsaturated");
    if (finiteValue(deltaEtaMax) && deltaEtaMax >= 0.0) shape.deltaEtaMax = static_cast<std::int32_t>(deltaEtaMax);
    if (finiteValue(deltaPhiMax) && deltaPhiMax >= 0.0) shape.deltaPhiMax = static_cast<std::int32_t>(deltaPhiMax);
    if (finiteValue(saturated) && saturated >= 0.0) shape.saturatedTowerCount = static_cast<std::uint32_t>(saturated);

    const auto* witness = m_photonBuilder->get_shower_feature_view(cluster, definitionName);
    if (!witness) throw std::runtime_error("missing configured builder feature view: " + definitionName);
    const auto captured = [witness](const std::string& key) { return witness->at(key); };
    shape.rawCenterEta = captured("raw_center_eta");
    shape.rawCenterPhi = captured("raw_center_phi");
    shape.centerEtaIndex = static_cast<std::int32_t>(captured("center_ieta"));
    shape.centerPhiIndex = static_cast<std::int32_t>(captured("center_iphi"));
    shape.centerOfGravityEta = captured("cog_eta_local");
    shape.centerOfGravityPhi = captured("cog_phi_local");
    shape.e11 = captured("e11");
    shape.e33 = captured("e33");
    shape.e32 = captured("e32");
    shape.e35 = captured("e35");
    shape.e11OverE33 = captured("e11_over_e33");
    shape.e32OverE35 = captured("e32_over_e35");
    shape.wetaCogX = captured("weta_cogx");
    shape.wphiCogX = captured("wphi_cogx");
    shape.weta33CogX = captured("weta33_cogx");
    shape.wphi33CogX = captured("wphi33_cogx");
    shape.et1 = captured("et1");
    shape.et2 = captured("et2");
    shape.et3 = captured("et3");
    shape.et4 = captured("et4");
    shape.energySource = static_cast<std::int32_t>(captured("energy_source"));
    shape.rectangularMembership = static_cast<std::int32_t>(captured("rectangular_membership"));
    shape.momentMembership = static_cast<std::int32_t>(captured("moment_membership"));
    shape.activeSumCellCount = static_cast<std::uint32_t>(captured("active_sum_cell_count"));
    shape.activeMomentCellCount = static_cast<std::uint32_t>(captured("active_moment_cell_count"));
    const std::array<const char*, 6> momentKeys{{"moment_eta_numerator", "moment_phi_numerator", "moment_denominator",
                                                 "moment33_eta_numerator", "moment33_phi_numerator",
                                                 "moment33_denominator"}};
    for (std::size_t index = 0; index < momentKeys.size(); ++index)
      shape.momentSums[index] = captured(momentKeys[index]);

    // Cell counts are optional support-census witnesses, not reconstructed
    // shower values. Every configured view uses the builder's own center.
    if (shape.centerEtaIndex >= 0)
    {
      for (int deta = -kGridHalfWidth; deta <= kGridHalfWidth; ++deta)
      {
        for (int dphi = -kGridHalfWidth; dphi <= kGridHalfWidth; ++dphi)
        {
          const int etaIndex = shape.centerEtaIndex + deta;
          const int phiIndex = wrapPhiIndex(shape.centerPhiIndex + dphi);
          addCell(cells, cemcTowers, etaIndex, phiIndex);
          auto found = cells.find({etaIndex, phiIndex});
          if (found == cells.end()) continue;
          Cell& cell = found->second;
          cell.inSupport = true;
          if (cell.owned) ++shape.ownedCellCount;
          if (!cell.present || !finiteValue(cell.energy))
          {
            if (cell.present) ++shape.nonFiniteCellCount;
            continue;
          }
          if (cell.good) ++shape.goodCellCount;
          if (cell.energy == 0.0) ++shape.zeroCellCount;
          if (cell.energy < 0.0) ++shape.negativeCellCount;
        }
      }
    }
    // Counts come from the same arithmetic owner as the values. The optional
    // cell census above only assembles retained support, including owned tails.
    shape.ownedCellCount = static_cast<std::uint32_t>(captured("owned_cell_count"));
    shape.goodCellCount = static_cast<std::uint32_t>(captured("good_cell_count"));
    shape.zeroCellCount = static_cast<std::uint32_t>(captured("zero_cell_count"));
    shape.negativeCellCount = static_cast<std::uint32_t>(captured("negative_cell_count"));
    shape.nonFiniteCellCount = static_cast<std::uint32_t>(captured("non_finite_cell_count"));
    shape.valid = finiteValue(shape.et1) && finiteValue(shape.et2) && finiteValue(shape.et3) &&
                  finiteValue(shape.et4) && finiteValue(shape.wetaCogX) && finiteValue(shape.wphiCogX) &&
                  finiteValue(shape.weta33CogX) && finiteValue(shape.wphi33CogX) && finiteValue(shape.e11OverE33) &&
                  finiteValue(shape.e32OverE35);
    photon.showerShapes.push_back(std::move(shape));
  }

  // PhotonCells remains available for an eventual losslessness decision.
  if (!m_config.output.writePhotonCells)
  {
    return;
  }

  for (const auto& [index, cell] : cells)
  {
    (void)index;

    PhotonCellRecord row;

    row.eventId = photon.eventId;
    row.photonId = photon.photonId;

    row.towerKey = cell.towerKey;
    row.subsystem = 0;

    row.etaIndex = cell.etaIndex;
    row.phiIndex = cell.phiIndex;

    row.calibratedEnergy = cell.present ? cell.energy : kNaN;

    row.clusterMapEnergy = cell.owned ? cell.clusterMapEnergy : kNaN;

    row.towerTimeSamples = cell.present ? cell.timeSamples : kNaN;

    row.towerStatus = cell.status;
    row.towerAvailable = cell.present;

    row.ownedByCluster = cell.owned;

    row.inShowerSupport = cell.inSupport;

    row.good = cell.present && cell.good;

    row.zeroEnergy = cell.present && finiteValue(cell.energy) && cell.energy == 0.0;

    row.negativeEnergy = cell.present && finiteValue(cell.energy) && cell.energy < 0.0;

    row.finite = cell.present && finiteValue(cell.energy) && finiteValue(cell.timeSamples);

    m_photonCells.push_back(std::move(row));
  }
}


// ============================================================================
// Reconstructed isolation
//
// PhotonClusterBuilder attaches reconstructed isolation values to each photon
// candidate. TreeProduction records every configured radius/method that is
// available rather than choosing the nominal analysis isolation.
//
// The three representations are:
//
//   CalorimeterRaw
//       raw CEMC + HCALIN + HCALOUT cone components with the photon candidate
//       already removed from the electromagnetic contribution;
//
//   CalorimeterSub1
//       Au+Au underlying-event-subtracted cone components;
//
//   Topocluster
//       pp topological-cluster isolation.
//
// These are measurements. The thresholds defining isolated/non-isolated regions
// remain downstream and can therefore be recalibrated without tree production.
// ============================================================================

std::vector<photonjet::IsolationRecord> PhotonJetTree::photonIsolation(const RawCluster* cluster,
                                                                       const PhotonRecord& photon) const
{
  std::vector<IsolationRecord> cones;
  for (const double radius : m_config.photon.isolationRadii)
  {
    const int code = static_cast<int>(std::lround(radius * 100.0));

    const std::string key = code == 30 ? "03" : code == 40 ? "04" : code == 20 ? "02" : code == 10 ? "01" : "";

    if (key.empty())
    {
      continue;  // the builder attaches cones at these radii only
    }

    // ------------------------------------------------------------------------
    // Raw calorimeter cone
    // ------------------------------------------------------------------------

    {
      IsolationRecord cone;

      cone.axisEta = namedShape(cluster, "isolation_axis_eta");

      cone.axisPhi = namedShape(cluster, "isolation_axis_phi");

      cone.eventId = photon.eventId;
      cone.photonId = photon.photonId;

      cone.radius = radius;
      cone.method = IsolationMethod::CalorimeterRaw;

      cone.electromagneticComponent = namedShape(cluster, "iso_" + key + "_emcal");

      cone.innerHadronicComponent = namedShape(cluster, "iso_" + key + "_hcalin");

      cone.outerHadronicComponent = namedShape(cluster, "iso_" + key + "_hcalout");

      cone.candidateRemoved = true;

      cone.valid = namedShape(cluster, "iso_" + key + "_valid") > 0.5 && finiteValue(cone.axisEta) &&
                   finiteValue(cone.axisPhi) && finiteValue(cone.electromagneticComponent) &&
                   finiteValue(cone.innerHadronicComponent) && finiteValue(cone.outerHadronicComponent);

      if (cone.valid)
      {
        cone.coneSum = cone.electromagneticComponent + cone.innerHadronicComponent + cone.outerHadronicComponent;
      }

      cones.push_back(cone);
    }

    // ------------------------------------------------------------------------
    // Au+Au SUB1 calorimeter cone
    // ------------------------------------------------------------------------

    if (isAuAu())
    {
      IsolationRecord cone;

      cone.axisEta = namedShape(cluster, "isolation_axis_eta");

      cone.axisPhi = namedShape(cluster, "isolation_axis_phi");

      cone.eventId = photon.eventId;
      cone.photonId = photon.photonId;

      cone.radius = radius;
      cone.method = IsolationMethod::CalorimeterSub1;

      cone.electromagneticComponent = namedShape(cluster, "iso_sub_" + key + "_emcal");

      cone.innerHadronicComponent = namedShape(cluster, "iso_sub_" + key + "_hcalin");

      cone.outerHadronicComponent = namedShape(cluster, "iso_sub_" + key + "_hcalout");

      cone.candidateRemoved = true;

      cone.valid = namedShape(cluster, "iso_sub_" + key + "_valid") > 0.5 && finiteValue(cone.axisEta) &&
                   finiteValue(cone.axisPhi) && finiteValue(cone.electromagneticComponent) &&
                   finiteValue(cone.innerHadronicComponent) && finiteValue(cone.outerHadronicComponent);

      if (cone.valid)
      {
        cone.coneSum = cone.electromagneticComponent + cone.innerHadronicComponent + cone.outerHadronicComponent;
      }
      else
      {
        cone.electromagneticComponent = cone.innerHadronicComponent = cone.outerHadronicComponent = kNaN;
      }

      cones.push_back(cone);
    }

    // ------------------------------------------------------------------------
    // pp topological-cluster cone
    // ------------------------------------------------------------------------

    if (isPP() && (key == "03" || key == "04"))
    {
      IsolationRecord cone;

      cone.axisEta = namedShape(cluster, "iso_topo_axis_eta");

      cone.axisPhi = namedShape(cluster, "iso_topo_axis_phi");

      cone.eventId = photon.eventId;
      cone.photonId = photon.photonId;

      cone.radius = radius;
      cone.method = IsolationMethod::Topocluster;

      cone.candidateRemoved = true;

      const double validFlag = namedShape(cluster, "iso_topo_" + key + "_valid");

      const double value = namedShape(cluster, "iso_topo_" + key);

      cone.valid = hasNamedShape(cluster, "iso_topo_" + key) && finiteValue(validFlag) && validFlag > 0.5 &&
                   finiteValue(cone.axisEta) && finiteValue(cone.axisPhi) && finiteValue(value);

      cone.coneSum = cone.valid ? value : kNaN;

      cones.push_back(cone);
    }
  }

  return cones;
}


// ============================================================================
// Isolation constituent witnesses
//
// Constituent-level storage is optional and exists to preserve more information
// than the final cone sum.
//
// pp
// --
// Store every topological cluster within the largest configured radius,
// including its signed transverse energy and displacement from the exact
// PhotonClusterBuilder isolation axis. This permits smaller cone sums to be
// reconstructed downstream.
//
// Au+Au
// -----
// Store raw and SUB1 tower witnesses from each calorimeter layer. Raw and SUB1
// quality states remain independent, and the same reconstructed isolation axis
// used by PhotonClusterBuilder is preserved.
//
// These rows document reconstruction inputs. They do not define the eventual
// isolated/non-isolated analysis regions.
// ============================================================================

void PhotonJetTree::capturePhotonIsolationConstituents(PHCompositeNode* topNode, const RawCluster* cluster,
                                                       const PhotonRecord& photon)
{
  if (!m_config.output.writeIsolationConstituents)
  {
    return;
  }

  const double maxRadius =
      *std::max_element(m_config.photon.isolationRadii.begin(), m_config.photon.isolationRadii.end());

  // Offsets are measured from the axis each system's cone uses: the photon
  // cluster direction in pp, the builder's centre-of-gravity tower axis in
  // Au+Au (the cluster direction when that axis is undefined). Azimuthal
  // offsets use ROOT's [-pi, pi) convention.
  const double builderAxisEta = namedShape(cluster, "isolation_axis_eta");
  const double builderAxisPhi = namedShape(cluster, "isolation_axis_phi");
  const double axisEta = isAuAu() && finiteValue(builderAxisEta) ? builderAxisEta : photon.eta;
  const double axisPhi = isAuAu() && finiteValue(builderAxisPhi) ? builderAxisPhi : photon.phi;

  // ==========================================================================
  // Au+Au tower-based isolation witnesses
  // ==========================================================================

  if (isAuAu())
  {
    /*
     * Historical Au+Au isolation uses the common 24x64 tower grid for all three
     * layers. The retowered CEMC therefore uses HCALIN geometry coordinates.
     *
     * Raw tower quality and SUB1 tower quality are independent facts. SUB1
     * quality controls the subtraction witness while both states are retained.
     *
     * Historical source anchor:
     *   RecoilJets_AuAu.cc:7970-8033
     */
    struct Layer
    {
      const char* raw;
      const char* sub;
      const char* geometry;
      RawTowerDefs::CalorimeterId id;
    };

    const Layer layers[] = {
        {"TOWERINFO_CALIB_CEMC_RETOWER", "TOWERINFO_CALIB_CEMC_RETOWER_SUB1", "TOWERGEOM_HCALIN", RawTowerDefs::HCALIN},
        {"TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALIN_SUB1", "TOWERGEOM_HCALIN", RawTowerDefs::HCALIN},
        {"TOWERINFO_CALIB_HCALOUT", "TOWERINFO_CALIB_HCALOUT_SUB1", "TOWERGEOM_HCALOUT", RawTowerDefs::HCALOUT}};

    for (int layer = 0; layer < 3; ++layer)
    {
      const auto& node = layers[layer];

      auto* raw = findNode::getClass<TowerInfoContainer>(topNode, node.raw);

      auto* sub = findNode::getClass<TowerInfoContainer>(topNode, node.sub);

      auto* geometry = findNode::getClass<RawTowerGeomContainer>(topNode, node.geometry);

      if (!raw || !sub || !geometry) continue;
      const unsigned int channelCount = std::min(raw->size(), sub->size());

      for (unsigned int channel = 0; channel < channelCount; ++channel)
      {
        auto* rawTower = raw->get_tower_at_channel(channel);

        auto* subTower = sub->get_tower_at_channel(channel);

        if (!rawTower || !subTower)
        {
          continue;  // historical witness census skips unavailable pairs
        }

        const auto key = TowerInfoDefs::encode_hcal(channel);

        const int etaBin = TowerInfoDefs::getCaloTowerEtaBin(key);

        const int phiBin = TowerInfoDefs::getCaloTowerPhiBin(key);

        auto* geom = geometry->get_tower_geometry(RawTowerDefs::encode_towerid(node.id, etaBin, phiBin));

        IsolationConstituentRecord row;

        row.eventId = photon.eventId;
        row.photonId = photon.photonId;

        row.radius = maxRadius;
        row.source = "calorimeter_sub1";

        row.nativeKey = channel;
        row.subsystem = layer;

        row.rawEnergy = rawTower->get_energy();
        row.subtractedEnergy = subTower->get_energy();

        row.maskState = (rawTower->get_isGood() ? 0U : 1U) | (subTower->get_isGood() ? 0U : 2U);

        row.masked = !subTower->get_isGood();
        row.qualityState = subTower->get_isGood() ? 1 : 0;

        if (etaBin < 0 || etaBin >= 24 || phiBin < 0 || phiBin >= 64)
        {
          row.qualityState = -2;
          m_isolationConstituents.push_back(row);
          continue;
        }
        const double r = geom ? std::hypot(geom->get_center_x(), geom->get_center_y()) : 0;
        if (!geom || !(r > 0))
        {
          row.qualityState = geom ? -3 : -1;
          m_isolationConstituents.push_back(row);
          continue;
        }

        /*
         * Recompute tower eta relative to the reconstructed collision vertex.
         * The geometry object itself is detector-centered.
         */
        const double eta = std::asinh((std::sinh(geom->get_eta()) * r - photon.producerVertexZ) / r);

        row.deltaEta = eta - axisEta;
        row.deltaPhi = TVector2::Phi_mpi_pi(geom->get_phi() - axisPhi);

        row.deltaR = std::hypot(row.deltaEta, row.deltaPhi);

        if (!finiteValue(row.deltaR) || row.deltaR >= maxRadius) continue;

        row.transverseEnergy = row.rawEnergy / std::cosh(eta);

        row.subtractedTransverseEnergy = row.subtractedEnergy / std::cosh(eta);

        /*
         * This flag is a witness only. The historical nominal builder removes
         * the photon ET once from the EM cone sum; it does not form its nominal
         * isolation by masking and resumming this stored constituent core.
         */
        row.candidateRemoved = layer == 0 && row.deltaR < 0.02;

        m_isolationConstituents.push_back(std::move(row));
      }
    }

    return;
  }

  // ==========================================================================
  // pp topological-cluster isolation witnesses
  // ==========================================================================

  auto* topoclusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.topoclusterNode);

  if (!topoclusters)
  {
    return;  // builder's explicit invalid topo state remains authoritative
  }

  const CLHEP::Hep3Vector vertex(0.0, 0.0, photon.producerVertexZ);

  const auto range = topoclusters->getClusters();

  for (auto it = range.first; it != range.second; ++it)
  {
    const RawCluster* topo = it->second;

    if (!topo)
    {
      continue;
    }

    const double eta = RawClusterUtility::GetPseudorapidity(*topo, vertex);

    const double phi = RawClusterUtility::GetAzimuthAngle(*topo, vertex);

    const double energy = topo->get_energy();

    const double dEta = eta - axisEta;

    const double dPhi = TVector2::Phi_mpi_pi(phi - axisPhi);

    const double dR = std::hypot(dEta, dPhi);

    if (!finiteValue(dR) || dR >= maxRadius)
    {
      continue;
    }

    IsolationConstituentRecord row;

    row.eventId = photon.eventId;
    row.photonId = photon.photonId;

    row.radius = maxRadius;
    row.source = "topocluster";

    row.nativeKey = it->first;
    row.subsystem = 10;  // persisted topological-cluster source code

    row.deltaEta = dEta;
    row.deltaPhi = dPhi;
    row.deltaR = dR;

    row.rawEnergy = energy;

    row.transverseEnergy = energy / std::cosh(eta);  // signed, as the builder sums it
    row.subtractedTransverseEnergy = row.transverseEnergy;
    row.qualityState = finiteValue(row.transverseEnergy) ? 1 : 0;

    /*
     * Historical candidate-removal witness: a topological cluster coincident
     * with the photon axis is identified as the candidate itself.
     */
    row.candidateRemoved = dR < 0.02;

    m_isolationConstituents.push_back(std::move(row));
  }
}

// ============================================================================
// Nominal photon labels
//
// Every stored candidate receives the frozen PhotonID label and isolation
// class before writing; retention never depends on them. Inputs are the
// stored values: double kinematics, centrality and cone sums, and the view's
// binary32 model features.
// ============================================================================

double photonjet::production::ModelFeature(const std::string& name, const PhotonRecord& p, const EventRecord& e,
                                           const ShowerShapeRecord* v)
{
  constexpr double missing = std::numeric_limits<double>::quiet_NaN();
  double value = missing;
  if (name == "cluster_Et")
    value = p.et;
  else if (name == "cluster_Eta")
    value = p.eta;
  else if (name == "vertexz")
    value = e.recoVertexZ;
  else if (name == "centrality")
    value = e.analysisCentrality >= 0 ? e.analysisCentrality : missing;
  else if (v)
  {
    const std::map<std::string, double> values = {{"cluster_weta_cogx", v->wetaCogX},
                                                  {"cluster_wphi_cogx", v->wphiCogX},
                                                  {"cluster_weta33_cogx", v->weta33CogX},
                                                  {"cluster_wphi33_cogx", v->wphi33CogX},
                                                  {"e11_over_e33", v->e11OverE33},
                                                  {"e32_over_e35", v->e32OverE35},
                                                  {"cluster_et1", v->et1},
                                                  {"cluster_et2", v->et2},
                                                  {"cluster_et3", v->et3},
                                                  {"cluster_et4", v->et4}};
    value = values.at(name);
  }
  else if (name == "e11_over_e33" || name == "e32_over_e35")
  {
    const int a = name == "e11_over_e33" ? 10 : 12, b = a + 1;
    value = std::isfinite(p.primaryShapes[a]) && std::isfinite(p.primaryShapes[b]) && p.primaryShapes[b] > 0
                ? p.primaryShapes[a] / p.primaryShapes[b]
                : missing;
  }
  else
  {
    const std::string suffix = name.substr(8);
    const auto it = std::find(PhotonRecord::primaryShapeNames.begin(), PhotonRecord::primaryShapeNames.end(), suffix);
    if (it == PhotonRecord::primaryShapeNames.end()) throw std::runtime_error("unknown frozen model feature: " + name);
    value = p.primaryShapes[std::distance(PhotonRecord::primaryShapeNames.begin(), it)];
  }
  return static_cast<float>(value);  // exact evaluator input, widened for storage
}

bool photonjet::production::ReferencePreselection(const PhotonRecord& photon)
{
  const auto ratio = [](const double a, const double b)
  { return std::isfinite(a) && std::isfinite(b) && b > 0 ? a / b : std::numeric_limits<double>::quiet_NaN(); };
  const double e11OverE33 = ratio(photon.primaryShapes[10], photon.primaryShapes[11]);
  const double e32OverE35 = ratio(photon.primaryShapes[12], photon.primaryShapes[13]);
  // Every input, w_phi and the photon ET included, must be finite.
  return std::isfinite(photon.primaryShapes[1]) && std::isfinite(photon.et) && e11OverE33 < .98 &&
         photon.primaryShapes[6] > .60 && photon.primaryShapes[6] < 1. && e32OverE35 > .80 && e32OverE35 < 1. &&
         photon.primaryShapes[0] < .60;
}


int photonjet::production::ActivePreselection(const PhotonConfig& config, const PhotonRecord& photon)
{
  const auto& pre = config.preselection;
  if (pre.variant == "reference") return ReferencePreselection(photon) ? 1 : 0;
  if (pre.variant != "newPPG12") return -1;
  double cut = std::numeric_limits<double>::quiet_NaN();
  for (const auto& model : config.models)
    if (model.scoreName == pre.npbScoreName) cut = model.configuredCut;
  const double npb = (pre.npbScoreName == "auau_npb_score" ? photon.auauNpbCapture : photon.npbCapture).rawScore;
  const auto ratio = [](const double a, const double b)
  { return std::isfinite(a) && std::isfinite(b) && b > 0 ? a / b : std::numeric_limits<double>::quiet_NaN(); };
  const double et1 = photon.primaryShapes[6];
  const double e11OverE33 = ratio(photon.primaryShapes[10], photon.primaryShapes[11]);
  const double e32OverE35 = ratio(photon.primaryShapes[12], photon.primaryShapes[13]);
  return std::isfinite(npb) && std::isfinite(et1) && std::isfinite(e11OverE33) && std::isfinite(e32OverE35) &&
         npb > cut && et1 > pre.et1Min && et1 < pre.et1Max && e11OverE33 < pre.e11OverE33Max &&
         e32OverE35 > pre.e32OverE35Min && e32OverE35 < pre.e32OverE35Max;
}


bool photonjet::production::ViewShapesFinite(const ShowerShapeRecord& v)
{
  return std::isfinite(v.wetaCogX) && std::isfinite(v.wphiCogX) && std::isfinite(v.weta33CogX) &&
         std::isfinite(v.wphi33CogX) && std::isfinite(v.e11OverE33) && std::isfinite(v.e32OverE35) &&
         std::isfinite(v.et1) && std::isfinite(v.et2) && std::isfinite(v.et3) && std::isfinite(v.et4);
}

namespace
{
// PPG12 pp regions in photon ET, as frozen in PhotonID/models/nominal.json:
// tight above the tight line, non-tight strictly inside the lower band.
int ppIdClass(const double score, const double et)
{
  using photonjet::PPG12IdLines;
  const double tight = PPG12IdLines::tight(et);
  const double low = PPG12IdLines::nonTightLow(et);
  const double high = PPG12IdLines::nonTightHigh(et);
  return score > tight ? 3 : (low < score && score < high) ? 1 : score >= high ? 2 : 0;
}

// Au+Au regions: the package's WP70/WP80 lines in centrality.
int auauIdClass(const double score, const double centrality, const photonjet::NominalPhotonIdConfig& id)
{
  const double wp70 = id.wp70.interceptGeV + id.wp70.slope * centrality;
  const double wp80 = id.wp80.interceptGeV + id.wp80.slope * centrality;
  return score > wp70 ? 3 : (0.1 < score && score < wp80) ? 1 : score >= wp80 ? 2 : 0;
}

// Distance between two scores in binary32 units in the last place of the stored one.
double float32Ulp(const double recomputed, const double stored)
{
  const float magnitude = std::fabs(static_cast<float>(stored));
  const double spacing =
      static_cast<double>(std::nextafter(magnitude, std::numeric_limits<float>::infinity())) - magnitude;
  return std::fabs(recomputed - stored) / spacing;
}

int isolationClass(const photonjet::NominalIsolationConfig& c, const photonjet::PhotonRecord& p,
                   const double centrality)
{
  if (c.witnessIndex < 0 || static_cast<std::size_t>(c.witnessIndex) >= p.isolationWitnesses.size())
    throw std::runtime_error("the nominal isolation witness was not captured");
  const auto& w = p.isolationWitnesses[c.witnessIndex];
  const bool reconstructed =
      w.reconstructedOrTruth == 0 && w.subtractionMethod == c.subtractionMethod && std::isfinite(w.coneSum);
  if (c.storedThresholds)
  {
    const bool valid = reconstructed && std::isfinite(w.threshold) && std::isfinite(w.sidebandThreshold) &&
                       w.sidebandThreshold >= w.threshold;
    return !valid ? -1 : w.coneSum < w.threshold ? 2 : w.coneSum > w.sidebandThreshold ? 0 : 1;
  }
  const double threshold = c.threshold.interceptGeV + c.threshold.slope * centrality;
  const bool valid = reconstructed && std::isfinite(centrality) && centrality >= c.centralityMin &&
                     centrality < c.centralityMaxExclusive;
  return !valid ? -1 : w.coneSum < threshold ? 2 : w.coneSum > threshold + c.threshold.sidebandGapGeV ? 0 : 1;
}
}  // namespace

photonjet::NominalPhotonRecord photonjet::production::LabelNominalPhoton(const PhotonConfig& config,
                                                                         const PinnedModel& model,
                                                                         const PhotonRecord& photon,
                                                                         const EventRecord& event,
                                                                         NominalReceiptRecord& receipt)
{
  const NominalPhotonIdConfig& id = config.nominalId;
  const auto view = std::find_if(photon.showerShapes.begin(), photon.showerShapes.end(),
                                 [&](const ShowerShapeRecord& v) { return v.definitionName == id.view; });
  if (view == photon.showerShapes.end())
    throw std::runtime_error("the nominal model view was not captured: " + id.view);
  std::vector<float> features;
  bool valid = ViewShapesFinite(*view);
  for (const std::string& name : id.featureOrder)
  {
    features.push_back(static_cast<float>(ModelFeature(name, photon, event, &*view)));
    valid = valid && std::isfinite(features.back());
  }
  NominalPhotonRecord label;
  const double et = photon.et, centrality = event.analysisCentrality;
  if (id.usesProducerScore)
  {
    // pp: the PPG12 producer score where ET reaches the package minimum; a
    // recompute of the same model must agree and keep every region.
    const double producer = photon.productionScore;
    if (et >= id.etMin)
    {
      if (valid && std::isfinite(producer) && producer != -1.0)
      {
        const double recomputed = model.score(features);
        const double ulp = float32Ulp(recomputed, producer);
        if (!(ulp <= id.recomputeUlpLimit))
          throw std::runtime_error("pp producer score differs from the recompute by " + std::to_string(ulp) +
                                   " float32 ulp");
        if (ppIdClass(recomputed, et) != ppIdClass(producer, et))
          throw std::runtime_error("pp recompute changes a photon-ID region");
        receipt.maxRecomputeUlp = std::max(receipt.maxRecomputeUlp, ulp);
        label.score = producer;
        label.state = 1;
        label.idClass = ppIdClass(producer, et);
      }
      else
        label.state = 2;
    }
  }
  else
  {
    const bool inside = et >= id.etMin && et < id.etMaxExclusive && std::isfinite(centrality) &&
                        centrality >= id.centralityMin && centrality < id.centralityMaxExclusive;
    if (inside && valid)
    {
      label.score = model.score(features);
      if (!std::isfinite(label.score)) throw std::runtime_error("a scored photon has a non-finite nominal score");
      label.state = 1;
      label.idClass = auauIdClass(label.score, centrality, id);
    }
    else if (inside)
      label.state = 2;
  }
  label.isolationClass = isolationClass(config.nominalIsolation, photon, centrality);
  ++receipt.photons;
  ++receipt.states[label.state];
  ++receipt.idClasses[label.idClass];
  ++receipt.isolationClasses[label.isolationClass];
  return label;
}

// Annotate every retained candidate without filtering or reordering it.
// The score/state/ID/isolation branches let downstream code apply the declared selection.
void PhotonJetTree::labelNominalPhotons()
{
  for (PhotonRecord& photon : m_photons)
    photon.nominal =
        photonjet::production::LabelNominalPhoton(m_config.photon, *m_nominalModel, photon, m_event, m_nominalReceipt);
}

// Jet capture and photon-jet pairing

namespace
{
constexpr unsigned int kInvalidJetId = std::numeric_limits<unsigned int>::max();


// ============================================================================
// Jet capture quality flags
//
// One bit per condition used to describe the integrity of a retained JetRecord.
// The mask records what was observed during capture—finite kinematics, valid
// raw/calibrated correspondence, area availability, and storage eligibility.
// It is diagnostic metadata, not a downstream jet-selection mask.
// ============================================================================

enum JetQualityBit : std::uint64_t
{
  kFiniteCorrectedPt = 1ULL << 0U,
  kFiniteRawPt = 1ULL << 1U,
  kFiniteEta = 1ULL << 2U,
  kFinitePhi = 1ULL << 3U,
  kFiniteArea = 1ULL << 4U,
  kCalibratedFlagSet = 1ULL << 5U,
  kCorrectedIdMatchesOrdinal = 1ULL << 6U,
  kRawIdCompatible = 1ULL << 7U,
  kContainerSizesMatch = 1ULL << 8U,
  kContainerRadiusMatches = 1ULL << 9U,
  kRawPtNonNegative = 1ULL << 10U,
  kInCaptureDomain = 1ULL << 11U
};


// ============================================================================
// Small jet-geometry helpers
// ============================================================================

/**
 * Read the FastJet active area associated with one jet.
 *
 * Jetv2 may store properties at container-defined indices, so the semantic
 * area property is first translated through the owning container. The direct
 * jet property is retained as the fallback for compatible representations.
 */
double jetArea(JetContainer* container, Jet* jet)
{
  if (container->has_property(Jet::PROPERTY::prop_area))
  {
    const Jet::PROPERTY index = container->property_index(Jet::PROPERTY::prop_area);

    const float area = jet->get_property(index);

    return finiteValue(area) ? area : kNaN;
  }

  if (jet->has_property(Jet::PROPERTY::prop_area))
  {
    const float area = jet->get_property(Jet::PROPERTY::prop_area);

    return finiteValue(area) ? area : kNaN;
  }

  return kNaN;
}


/**
 * Check the configured jet radius against container metadata when available.
 *
 * Some historical containers do not persist a usable radius. That absence is
 * tolerated because the collection identity is already configuration-bound.
 * A finite positive stored radius, however, must agree with the configured
 * collection or the production contract is inconsistent.
 */
bool radiusCompatible(JetContainer* container, const double expected)
{
  const float stored = container->get_jetpar_R();

  if (!finiteValue(stored) || stored <= 0.0F)
  {
    return true;
  }

  return std::abs(static_cast<double>(stored) - expected) < 1.0e-4;
}

}  // namespace


// ============================================================================
// Configured reconstructed-jet views
// ============================================================================

/**
 * Capture every jet collection declared for this production profile.
 *
 * Event-level availability is recorded independently for each view/radius.
 * When reconstructed objects are outside the producer's vertex domain, the
 * collection is NotApplicable rather than falsely recorded as an empty valid
 * collection.
 */
void PhotonJetTree::captureJets(PHCompositeNode* topNode)
{
  // Workfest DATA keeps every event, but recoil support belongs only to
  // photon-bearing events. SIM retains its full reconstructed-jet population.
  // An inapplicable event has no collection occurrence, as in the original
  // producer; file-level collection definitions remain in metadata.
  if (!m_event.recoObjectVertexInDomain || (isData() && m_photons.empty())) return;

  for (const JetNodeConfig& collection : m_config.jets.nodes)
  {
    m_event.recoJetView.push_back(static_cast<int>(collection.view));

    m_event.recoJetRadiusCode.push_back(static_cast<int>(photonjet::radiusCode(collection.radius)));

    m_event.recoJetCaptureState.push_back(static_cast<int>(CaptureState::NotApplicable));

    m_event.recoJetCount.push_back(0);

    /*
     * Once the event is in the reconstructed-object domain, failure to obtain
     * the configured collection is an unavailable-input state until the
     * collection capture succeeds.
     */
    m_event.recoJetCaptureState.back() = static_cast<int>(CaptureState::InputUnavailable);

    const auto before = m_jets.size();

    captureJetCollection(topNode, collection);

    m_event.recoJetCaptureState.back() = static_cast<int>(CaptureState::Complete);

    m_event.recoJetCount.back() = static_cast<std::uint32_t>(m_jets.size() - before);
  }
}


// ============================================================================
// One raw/calibrated jet collection
// ============================================================================

/**
 * Capture one configured reconstructed jet view.
 *
 * The upstream calibration produces one corrected collection corresponding to
 * one raw collection. Corrected jet ordinal i is therefore expected to bind to
 * raw ordinal i, and the corrected jet id records that ordinal.
 *
 * Both collections are read so each retained row contains the raw momentum,
 * corrected momentum, native ids, area, and calibration witnesses together.
 * Multiplicity or radius disagreements indicate a broken upstream contract and
 * stop the run rather than producing ambiguously paired jets.
 */
void PhotonJetTree::captureJetCollection(PHCompositeNode* topNode, const JetNodeConfig& collection)
{
  JetContainer* corrected = findNode::getClass<JetContainer>(topNode, collection.correctedNode);

  JetContainer* raw = findNode::getClass<JetContainer>(topNode, collection.rawNode);

  if (!corrected)
  {
    throw std::runtime_error("missing corrected jet node '" + collection.correctedNode + "'");
  }

  if (!raw)
  {
    throw std::runtime_error("missing raw jet node '" + collection.rawNode + "'");
  }

  /*
   * Raw and calibrated collections are expected to have one-to-one ordinal
   * correspondence. A size mismatch means that relationship cannot be
   * established safely.
   */
  const bool sizesMatch = raw->size() == corrected->size();

  if (!sizesMatch)
  {
    throw std::runtime_error("raw and corrected jet multiplicity differ for '" + collection.rawNode + "' and '" +
                             collection.correctedNode + "'");
  }

  const bool radiusMatches = radiusCompatible(corrected, collection.radius) && radiusCompatible(raw, collection.radius);

  if (!radiusMatches)
  {
    throw std::runtime_error("configured radius disagrees with the container metadata of '" + collection.correctedNode +
                             "'");
  }

  /*
   * deterministicOrder counts retained objects only and restarts for every
   * configured collection. It intentionally differs from both native jet id
   * and input-container ordinal.
   */
  std::uint32_t retainedOrder = 0;

  // Count once: the loop never changes the container, and the repo rule
  // forbids size() calls in loop conditions (see the loop-count test).
  const auto correctedCount = corrected->size();
  for (unsigned int ordinal = 0; ordinal < correctedCount; ++ordinal)
  {
    Jet* correctedJet = corrected->get_jet(ordinal);
    Jet* rawJet = raw->get_jet(ordinal);

    /*
     * Match the reference capture semantics: invalid calibrated objects are
     * removed before a raw/calibrated row is constructed. Retained ordering
     * therefore counts only surviving corrected jets.
     */
    if (!correctedJet || !finiteValue(correctedJet->get_pt()) || correctedJet->get_pt() < 0)
    {
      continue;
    }

    if (!rawJet)
    {
      throw std::runtime_error("null jet at ordinal " + std::to_string(ordinal) + " in '" + collection.correctedNode +
                               "'");
    }

    JetRecord jet;

    jet.eventId = m_event.eventId;

    jet.view = collection.view;
    jet.radius = collection.radius;

    jet.inputIdentity = collection.inputIdentity;
    jet.subtractionIdentity = collection.subtractionIdentity;

    jet.nativeKey = correctedJet->get_id();
    jet.nativeRawKey = rawJet->get_id();

    jet.encounterOrdinal = ordinal;

    /*
     * JetCalib rebuilds each calibrated jet from the raw jet's binary32
     * (pT, eta, phi) as (c cos phi, c sin phi, c sinh eta) for the calibrated
     * pT c. The stored direction is that construction at unit scale (c equal
     * to the raw pT), so it never depends on the JES value. A negative JES
     * reverses the calibrated vector: that jet's JES is undefined, and it
     * keeps the unit-scale pT. Raw pT is retained beside both.
     */
    const float rawJetPt = rawJet->get_pt();
    const float rawJetEta = rawJet->get_eta();
    const float rawJetPhi = rawJet->get_phi();
    Jetv2 unitScale;
    unitScale.set_px(rawJetPt * std::cos(rawJetPhi));
    unitScale.set_py(rawJetPt * std::sin(rawJetPhi));
    unitScale.set_pz(rawJetPt * std::sinh(rawJetEta));
    const double alignment = static_cast<double>(correctedJet->get_px()) * unitScale.get_px() +
                             static_cast<double>(correctedJet->get_py()) * unitScale.get_py() +
                             static_cast<double>(correctedJet->get_pz()) * unitScale.get_pz();
    const bool jesDefined = alignment > 0.0;

    jet.rawPt = rawJet->get_pt();
    jet.correctedPt = jesDefined ? correctedJet->get_pt() : unitScale.get_pt();
    jet.eta = unitScale.get_eta();
    jet.phi = unitScale.get_phi();
    jet.unitScalePt = unitScale.get_pt();

    /*
     * JetCalib does not propagate the FastJet area into the corrected object,
     * so the active area is read from the corresponding raw jet.
     */
    jet.area = jetArea(raw, rawJet);
    if (collection.view != JetView::AuAuNoSub)
    {
      const auto fractionReader = m_jetFractionReaders.find(collection.rawNode);
      if (fractionReader == m_jetFractionReaders.end())
        throw std::runtime_error("no calorimeter fraction capture for raw collection " + collection.rawNode);
      jet.calorimeterFractions = fractionReader->second(rawJet);
    }

    // ------------------------------------------------------------------------
    // Calibration and raw/corrected binding
    // ------------------------------------------------------------------------

    const bool calibratedFlag = correctedJet->get_isCalib() != 0;

    const unsigned int correctedId = correctedJet->get_id();

    const unsigned int rawId = rawJet->get_id();

    const bool correctedIdMatches = correctedId == ordinal;

    /*
     * Raw native ids are retained as evidence but historically need not equal
     * their container ordinal. The corrected id is the authoritative ordinal
     * binding produced by JetCalib.
     */
    const bool rawIdCompatible = rawId == kInvalidJetId || rawId == ordinal || rawId == correctedId;

    const bool rawPtNonNegative = finiteValue(jet.rawPt) && jet.rawPt >= 0.0;

    jet.calibrationValid = jesDefined && rawPtNonNegative && finiteValue(jet.correctedPt) && finiteValue(jet.eta) &&
                           finiteValue(jet.phi) && calibratedFlag && correctedIdMatches;

    /*
     * These are production-contract requirements, not analysis-quality cuts.
     * A broken raw/calibrated binding cannot be represented as an ordinary jet
     * without losing the meaning of the stored corrected momentum.
     */
    if (!rawPtNonNegative || !calibratedFlag || !correctedIdMatches)
    {
      throw std::runtime_error("invalid calibrated/raw jet binding: " + collection.correctedNode);
    }

    if (!finiteValue(jet.area) || jet.area < 0.0)
    {
      throw std::runtime_error("missing or invalid FastJet active area: " + collection.rawNode);
    }

    // ------------------------------------------------------------------------
    // Persisted quality witnesses
    // ------------------------------------------------------------------------

    std::uint64_t quality = 0;

    if (finiteValue(jet.correctedPt)) quality |= kFiniteCorrectedPt;
    if (finiteValue(jet.rawPt)) quality |= kFiniteRawPt;
    if (finiteValue(jet.eta)) quality |= kFiniteEta;
    if (finiteValue(jet.phi)) quality |= kFinitePhi;
    if (finiteValue(jet.area)) quality |= kFiniteArea;

    if (calibratedFlag) quality |= kCalibratedFlagSet;
    if (correctedIdMatches) quality |= kCorrectedIdMatchesOrdinal;
    if (rawIdCompatible) quality |= kRawIdCompatible;

    if (sizesMatch) quality |= kContainerSizesMatch;
    if (radiusMatches) quality |= kContainerRadiusMatches;
    if (rawPtNonNegative) quality |= kRawPtNonNegative;

    // ------------------------------------------------------------------------
    // Storage domain
    //
    // TreeProduction keeps the broad reconstructed inventory. Jet pT and eta
    // selections for a measurement belong downstream in TreeToHists.
    // ------------------------------------------------------------------------

    if (!jetInCaptureDomain(jet))
    {
      continue;
    }

    quality |= kInCaptureDomain;
    jet.qualityBitmask = quality;

    /*
     * Retained order is encounter order after invalid/out-of-domain rows have
     * been skipped. It restarts independently for every configured collection.
     */
    jet.deterministicOrder = retainedOrder;

    jet.jetId =
        photonjet::makeJetIdentity(m_event.eventId, jet.view, jet.radius, jet.nativeKey, jet.deterministicOrder);

    ++retainedOrder;

    m_jets.push_back(std::move(jet));
  }
}


// ============================================================================
// Jet storage domain
// ============================================================================

/**
 * Decide whether a reconstructed jet is representable in the canonical tree.
 *
 * This is deliberately a broad storage requirement, not the measurement jet
 * selection. Any finite calibrated momentum/direction is retained; analysis
 * pT, eta, and recoil requirements remain downstream.
 */
bool PhotonJetTree::jetInCaptureDomain(const JetRecord& jet) const
{
  if (!finiteValue(jet.correctedPt) || !finiteValue(jet.eta) || !finiteValue(jet.phi))
  {
    return false;
  }

  return true;
}


// ============================================================================
// Photon x jet relationships
//
// When enabled, every retained photon is paired with every retained jet.
// TreeProduction records the geometric and momentum-balance facts but does not
// discard non-away-side pairs. Downstream analysis can therefore reproduce or
// change the recoil requirement from the stored delta-phi.
//
// Pairs are stored in the established producer's order: for each retained jet
// (views, then radii, then jets), every retained photon in candidate order.
// ============================================================================

std::vector<const photonjet::JetRecord*> photonjet::production::CanonicalJetOrder(
    const std::vector<photonjet::JetRecord>& jets)
{
  std::vector<const JetRecord*> ordered;
  ordered.reserve(jets.size());
  for (const bool unsubtracted : {false, true})
    for (const int radiusCode : {2, 3, 4})
      for (const JetRecord& jet : jets)
        if (jet.radius == radiusCode / 10. && (jet.view == JetView::AuAuNoSub) == unsubtracted) ordered.push_back(&jet);
  return ordered;
}


void PhotonJetTree::buildPhotonJetPairs()
{
  m_photonJetPairs.clear();

  if (!m_config.output.writePhotonJetPairs || m_photons.empty() || m_jets.empty())
  {
    return;
  }

  m_photonJetPairs.reserve(m_photons.size() * m_jets.size());

  // Each jet in canonical collection order, with every photon, as the
  // established producer emitted the pairs.
  for (const JetRecord* jetPointer : photonjet::production::CanonicalJetOrder(m_jets))
  {
    const JetRecord& jet = *jetPointer;
    for (const PhotonRecord& photon : m_photons)
    {
      PhotonJetPairRecord pair;

      pair.eventId = m_event.eventId;

      pair.photonId = photon.photonId;
      pair.jetId = jet.jetId;

      pair.jetView = jet.view;
      pair.jetRadius = jet.radius;

      /*
       * Historical pair semantics use photon_rank = 0 independently of the
       * candidate's position in m_photons. jetRank preserves the retained jet
       * ordinal rather than introducing a new pT ranking.
       */
      pair.photonRank = 0;
      pair.jetRank = jet.deterministicOrder;

      // Azimuthal separation in ROOT's [-pi, pi) convention, then |.|.
      pair.deltaPhi = std::fabs(TVector2::Phi_mpi_pi(jet.phi - photon.phi));

      /*
       * xJgamma uses the calibrated reconstructed jet momentum divided by the
       * reconstructed photon transverse energy. A non-positive photon
       * denominator, or a jet whose JES is undefined (negative), leaves the
       * quantity undefined.
       */
      if (jet.calibrationValid && photon.et > 0.0)
      {
        pair.xJGamma = jet.correctedPt / photon.et;
      }

      /*
       * Keep the recoil decision as a witness instead of a retention cut.
       * The boundary is inclusive: |delta phi| >= 7 pi / 8.
       */
      pair.recoilWitness = pair.deltaPhi >= photonjet::kRecoilDeltaPhiMin;

      pair.pairId = photonjet::makePairIdentity(pair.photonId, pair.jetId, pair.jetRank);

      m_photonJetPairs.push_back(std::move(pair));
    }
  }
}
