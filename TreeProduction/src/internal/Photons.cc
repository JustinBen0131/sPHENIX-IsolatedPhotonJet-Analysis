/**
 * @file Photons.cc
 * @brief Capture builder-produced photon observables and optional witnesses.
 *
 * PhotonClusterBuilder owns reconstructed timing, showers, and isolation.
 * This file captures those values in PhotonJetTree records:
 *
 *   candidate identity and kinematics
 *     -> broad tree-storage acceptance
 *     -> cluster timing
 *     -> configured shower-shape representations
 *     -> calorimeter-cell witnesses
 *     -> reconstructed isolation observables
 *     -> optional isolation-constituent witnesses
 *
 * The candidate inventory is intentionally model-independent. Photon-ID scores, tight/non-tight working points, isolated/non-isolated classifications, and
 * ABCD regions are downstream decisions and never determine whether a photon is retained here.
 *
 * Responsibility boundaries
 * -------------------------
 *   ProductionReconstruction.cc
 *       Registers PhotonClusterBuilder and the detector inputs consumed here.
 *
 *   Truth.cc
 *       Owns simulation truth and reconstruction-to-truth associations.
 *
 *   Output.cc
 *       Serializes the completed records to ROOT.
 *
 * PhotonJetTree preserves the reconstructed primitives required to evaluate
 * photon-ID models without rereading DSTs. Optional cell/constituent rows are
 * low-level audit witnesses, not another nominal reconstruction.
 */

#include "../PhotonJetTree.h"

#include <caloreco/PhotonClusterBuilder.h>

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

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPi = 3.14159265358979323846;

constexpr int kCemcEtaBins = 96;
constexpr int kCemcPhiBins = 256;
constexpr int kGridHalfWidth = 3;  // 7x7 grid about the seed tower

bool finite(const double value)
{
  return std::isfinite(value);
}


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


/**
 * Wrap a physical azimuth onto (-pi, pi].
 */
double wrapPhi(double phi)
{
  while (phi > kPi) phi -= 2.0 * kPi;
  while (phi <= -kPi) phi += 2.0 * kPi;

  return phi;
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

  return found == shapes.end()
             ? kNaN
             : static_cast<double>(found->second);
}

bool hasNamedShape(const RawCluster* cluster, const std::string& name)
{
  const auto& shapes = cluster->get_all_shower_shapes();

  return shapes.find(name) != shapes.end();
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
void addCell(CellMap& cells,
             TowerInfoContainer* towers,
             const int etaIndex,
             const int phiIndex)
{
  if (etaIndex < 0 || etaIndex >= kCemcEtaBins)
  {
    return;
  }

  const int phi = wrapPhiIndex(phiIndex);

  Cell& cell =
      cells[{etaIndex, phi}];

  if (cell.towerKey != 0)
  {
    return;
  }

  const unsigned int key =
      TowerInfoDefs::encode_emcal(
          etaIndex,
          phi);

  cell.towerKey = key;
  cell.etaIndex = etaIndex;
  cell.phiIndex = phi;

  if (TowerInfo* tower =
          towers->get_tower_at_channel(TowerInfoDefs::decode_emcal(key)))
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
  RawClusterContainer* candidates =
      findNode::getClass<RawClusterContainer>(
          topNode,
          m_config.photon.photonNode);

  if (!candidates)
  {
    throw std::runtime_error(
        "missing photon candidate node '" +
        m_config.photon.photonNode +
        "'");
  }

  TowerInfoContainer* cemcTowers =
      findNode::getClass<TowerInfoContainer>(
          topNode,
          m_config.nodes.cemcTowers);

  if (!cemcTowers)
  {
    throw std::runtime_error(
        "missing calibrated CEMC tower node '" +
        m_config.nodes.cemcTowers +
        "'");
  }

  // Keep the existing event-summary domain for jet capture. Photons have a
  // separately selected reconstructed vertex (notably the pp-SIM MBD component).
  m_event.recoObjectVertexInDomain =
      m_event.recoVertexValid &&
      std::abs(m_event.recoVertexZ) < m_config.photon.objectVertexAbsZMaxCm;

  if (!m_photonBuilder)
    throw std::runtime_error("photon capture requires the registered builder's event witness");

  const double photonVertexZ = m_photonBuilder->get_vertex_z();
  if (!finite(photonVertexZ))
  {
    m_event.recoPhotonCaptureState = CaptureState::InputUnavailable;
    return;
  }
  if (!(std::abs(photonVertexZ) < m_config.photon.objectVertexAbsZMaxCm))
  {
    m_event.recoPhotonCaptureState = CaptureState::NotApplicable;
    return;
  }
  m_event.recoPhotonCaptureState = CaptureState::Complete;

  std::uint32_t encountered = 0;

  for (const auto& [clusterKey, cluster] :
       candidates->getClustersMap())
  {
    const std::uint32_t ordinal =
        encountered++;

    if (!cluster)
    {
      ++m_event.recoPhotonUnclassifiableCount;
      continue;
    }

    PhotonRecord photon;

    photon.eventId = m_event.eventId;
    photon.nativeClusterKey = clusterKey;
    photon.producerVertexZ = namedShape(cluster, "vertex_z");
    photon.energy = cluster->get_energy();
    photon.et = namedShape(cluster, "cluster_pt");
    photon.eta = namedShape(cluster, "cluster_eta");
    photon.phi = namedShape(cluster, "cluster_phi");

    if (!finite(photon.producerVertexZ) || photon.producerVertexZ != photonVertexZ)
      throw std::runtime_error("photon vertex does not match the current builder event");

    photon.kinematicsFinite =
        finite(photon.energy) &&
        finite(photon.et) &&
        finite(photon.eta) &&
        finite(photon.phi);

    if (!photonInCaptureDomain(photon))
    {
      if (!photon.kinematicsFinite)
      {
        ++m_event.recoPhotonUnclassifiableCount;
      }

      continue;
    }

    /*
     * encounterOrdinal records where this candidate appeared in the upstream
     * candidate container. photonId combines event identity, native cluster key,
     * and that ordinal into the stable relational key used downstream.
     */
    photon.encounterOrdinal = ordinal;

    photon.photonId =
        photonjet::makePhotonIdentity(
            m_event.eventId,
            clusterKey,
            ordinal);

    capturePhotonTiming(
        cluster,
        photon);

    capturePhotonShowerShapes(
        cluster,
        cemcTowers,
        photon);

    capturePhotonIsolation(
        cluster,
        photon);

    capturePhotonIsolationConstituents(
        topNode,
        cluster,
        photon);

    /*
     * Truth.cc owns the dominant-primary calculation and final photon-truth
     * relation. Until that stage runs, simulation explicitly records the
     * evaluator as unavailable; DATA records truth as not applicable.
     */
    photon.dominantTruth.state =
        isSimulation()
            ? DominantTruthState::EvaluatorUnavailable
            : DominantTruthState::NotSimulation;

    m_photons.push_back(
        std::move(photon));
  }
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

  if (photon.et < m_config.photon.minEtGeV ||
      photon.et >= m_config.photon.maxEtGeV)
  {
    return false;
  }

  if (std::abs(photon.eta) >=
      m_config.photon.maxAbsEta)
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

void PhotonJetTree::capturePhotonTiming(
    const RawCluster* cluster,
    PhotonRecord& photon) const
{
  PhotonTimingRecord& timing = photon.timing;
  timing.energyWeightedNumerator = namedShape(cluster, "time_energy_numerator");
  timing.energyDenominator = namedShape(cluster, "time_energy_denominator");
  timing.meanTimeSamples = namedShape(cluster, "mean_time");

  const double count = namedShape(cluster, "time_contributing_towers");
  const bool countValid = finite(count) && count >= 0.0 &&
      count <= std::numeric_limits<std::uint32_t>::max() && std::floor(count) == count;
  if (countValid)
  {
    timing.contributingTowerCount = static_cast<std::uint32_t>(count);
  }

  timing.finite = finite(timing.meanTimeSamples);
  timing.valid = namedShape(cluster, "time_valid") > 0.5 && countValid &&
      finite(timing.energyWeightedNumerator) &&
      finite(timing.energyDenominator) && timing.energyDenominator > 0.0 &&
      timing.finite;
  timing.timeNs = timing.valid
      ? timing.meanTimeSamples * m_config.photon.timingSampleNs : kNaN;
}

void PhotonJetTree::capturePhotonShowerShapes(
    const RawCluster* cluster,
    TowerInfoContainer* cemcTowers,
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
    const double floor = definitionName == "H70" ? 0.070 :
                         definitionName == "H0"  ? 0.0 :
                         throw std::runtime_error("unsupported builder shower view: " + definitionName);
    // The builder stores binary32 floors. Match that exact representation:
    // nearby thresholds are different views, not approximate H70 aliases.
    const bool primary = primaryFloor == static_cast<double>(static_cast<float>(floor));
    const std::string prefix = primary ? "" : "shower_" + definitionName + "_";
    const auto value = [&](const std::string& key) { return namedShape(cluster, prefix + key); };

    ShowerShapeRecord shape;
    shape.definitionName = definitionName;
    shape.energyFloorGeV = floor;
    shape.et1 = value("et1"); shape.et2 = value("et2");
    shape.et3 = value("et3"); shape.et4 = value("et4");
    shape.e11 = value("e11"); shape.e13 = value("e13");
    shape.e15 = value("e15"); shape.e17 = value("e17");
    shape.e22 = value("e22"); shape.e31 = value("e31");
    shape.e32 = value("e32"); shape.e33 = value("e33");
    shape.e35 = value("e35"); shape.e37 = value("e37");
    shape.e51 = value("e51"); shape.e52 = value("e52");
    shape.e53 = value("e53"); shape.e55 = value("e55");
    shape.e57 = value("e57"); shape.e71 = value("e71");
    shape.e72 = value("e72"); shape.e73 = value("e73");
    shape.e75 = value("e75"); shape.e77 = value("e77");
    shape.weta = value("weta"); shape.wphi = value("wphi");
    shape.wetaCog = value("weta_cog"); shape.wphiCog = value("wphi_cog");
    shape.wetaCogX = value("weta_cogx"); shape.wphiCogX = value("wphi_cogx");
    shape.weta33CogX = value("weta33_cogx");
    shape.wphi33CogX = value("wphi33_cogx");
    shape.w32 = value("w32"); shape.w52 = value("w52");
    shape.w72 = value("w72");
    shape.centerOfGravityEta = value("cog_eta_local");
    shape.centerOfGravityPhi = value("cog_phi_local");
    shape.deltaEtaCog = value("detacog");
    shape.deltaPhiCog = value("dphicog");
    shape.radialSpread = value("drad");

    const double centerEta = value("center_ieta");
    const double centerPhi = value("center_iphi");
    if (finite(centerEta) && centerEta >= 0.0 && centerEta < kCemcEtaBins &&
        finite(centerPhi) && centerPhi >= 0.0 && centerPhi < kCemcPhiBins)
    {
      shape.centerEtaIndex = static_cast<std::int32_t>(centerEta);
      shape.centerPhiIndex = static_cast<std::int32_t>(centerPhi);
    }
    const double deltaEtaMax = value("detamax");
    const double deltaPhiMax = value("dphimax");
    const double saturated = value("nsaturated");
    if (finite(deltaEtaMax) && deltaEtaMax >= 0.0)
      shape.deltaEtaMax = static_cast<std::int32_t>(deltaEtaMax);
    if (finite(deltaPhiMax) && deltaPhiMax >= 0.0)
      shape.deltaPhiMax = static_cast<std::int32_t>(deltaPhiMax);
    if (finite(saturated) && saturated >= 0.0)
      shape.saturatedTowerCount = static_cast<std::uint32_t>(saturated);

    shape.valid = value("shower_shape_valid") > 0.5 &&
        value("moment33_valid") > 0.5 &&
        finite(shape.et1) && finite(shape.et2) &&
        finite(shape.et3) && finite(shape.et4) &&
        finite(shape.wetaCogX) && finite(shape.wphiCogX) &&
        finite(shape.weta33CogX) && finite(shape.wphi33CogX);

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
          if (!cell.present || !finite(cell.energy))
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
    photon.showerShapes.push_back(std::move(shape));
  }

  // PhotonCells remains available for an eventual losslessness decision.
  if (!m_config.output.writePhotonCells)
  {
    return;
  }

  for (const auto& [index, cell] : cells)
  {
    (void) index;

    PhotonCellRecord row;

    row.eventId = photon.eventId;
    row.photonId = photon.photonId;

    row.towerKey = cell.towerKey;
    row.subsystem = 0;

    row.etaIndex = cell.etaIndex;
    row.phiIndex = cell.phiIndex;

    row.calibratedEnergy =
        cell.present
            ? cell.energy
            : kNaN;

    row.clusterMapEnergy =
        cell.owned
            ? cell.clusterMapEnergy
            : kNaN;

    row.towerTimeSamples =
        cell.present
            ? cell.timeSamples
            : kNaN;

    row.towerStatus =
        cell.status;

    row.ownedByCluster =
        cell.owned;

    row.inShowerSupport =
        cell.inSupport;

    row.good =
        cell.present &&
        cell.good;

    row.zeroEnergy =
        cell.present &&
        finite(cell.energy) &&
        cell.energy == 0.0;

    row.negativeEnergy =
        cell.present &&
        finite(cell.energy) &&
        cell.energy < 0.0;

    row.finite =
        cell.present &&
        finite(cell.energy) &&
        finite(cell.timeSamples);

    m_photonCells.push_back(
        std::move(row));
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

void PhotonJetTree::capturePhotonIsolation(
    const RawCluster* cluster,
    PhotonRecord& photon)
{
  for (const double radius :
       m_config.photon.isolationRadii)
  {
    const int code =
        static_cast<int>(
            std::lround(
                radius *
                100.0));

    const std::string key =
        code == 30 ? "03" :
        code == 40 ? "04" :
        code == 20 ? "02" :
        code == 10 ? "01" :
                     "";

    if (key.empty())
    {
      continue;  // the builder attaches cones at these radii only
    }

    // ------------------------------------------------------------------------
    // Raw calorimeter cone
    // ------------------------------------------------------------------------

    {
      IsolationRecord cone;

      cone.axisEta =
          namedShape(
              cluster,
              "isolation_axis_eta");

      cone.axisPhi =
          namedShape(
              cluster,
              "isolation_axis_phi");

      cone.eventId = photon.eventId;
      cone.photonId = photon.photonId;

      cone.radius = radius;
      cone.method = IsolationMethod::CalorimeterRaw;

      cone.electromagneticComponent =
          namedShape(
              cluster,
              "iso_" + key + "_emcal");

      cone.innerHadronicComponent =
          namedShape(
              cluster,
              "iso_" + key + "_hcalin");

      cone.outerHadronicComponent =
          namedShape(
              cluster,
              "iso_" + key + "_hcalout");

      cone.candidateRemoved = true;

      cone.valid =
          namedShape(cluster, "iso_" + key + "_valid") > 0.5 &&
          finite(cone.axisEta) && finite(cone.axisPhi) &&
          finite(cone.electromagneticComponent) &&
          finite(cone.innerHadronicComponent) &&
          finite(cone.outerHadronicComponent);

      if (cone.valid)
      {
        cone.coneSum =
            cone.electromagneticComponent +
            cone.innerHadronicComponent +
            cone.outerHadronicComponent;
      }

      photon.isolation.push_back(cone);
    }

    // ------------------------------------------------------------------------
    // Au+Au SUB1 calorimeter cone
    // ------------------------------------------------------------------------

    if (isAuAu())
    {
      IsolationRecord cone;

      cone.axisEta =
          namedShape(
              cluster,
              "isolation_axis_eta");

      cone.axisPhi =
          namedShape(
              cluster,
              "isolation_axis_phi");

      cone.eventId = photon.eventId;
      cone.photonId = photon.photonId;

      cone.radius = radius;
      cone.method = IsolationMethod::CalorimeterSub1;

      cone.electromagneticComponent =
          namedShape(
              cluster,
              "iso_sub_" + key + "_emcal");

      cone.innerHadronicComponent =
          namedShape(
              cluster,
              "iso_sub_" + key + "_hcalin");

      cone.outerHadronicComponent =
          namedShape(
              cluster,
              "iso_sub_" + key + "_hcalout");

      cone.candidateRemoved = true;

      cone.valid =
          namedShape(cluster, "iso_sub_" + key + "_valid") > 0.5 &&
          finite(cone.axisEta) && finite(cone.axisPhi) &&
          finite(cone.electromagneticComponent) &&
          finite(cone.innerHadronicComponent) &&
          finite(cone.outerHadronicComponent);

      if (cone.valid)
      {
        cone.coneSum =
            cone.electromagneticComponent +
            cone.innerHadronicComponent +
            cone.outerHadronicComponent;
      }
      else
      {
        cone.electromagneticComponent = cone.innerHadronicComponent = cone.outerHadronicComponent = kNaN;
      }

      photon.isolation.push_back(cone);
    }

    // ------------------------------------------------------------------------
    // pp topological-cluster cone
    // ------------------------------------------------------------------------

    if (isPP() &&
        (key == "03" ||
         key == "04"))
    {
      IsolationRecord cone;

      cone.axisEta =
          namedShape(
              cluster,
              "iso_topo_axis_eta");

      cone.axisPhi =
          namedShape(
              cluster,
              "iso_topo_axis_phi");

      cone.eventId = photon.eventId;
      cone.photonId = photon.photonId;

      cone.radius = radius;
      cone.method = IsolationMethod::Topocluster;

      cone.candidateRemoved = true;

      const double validFlag =
          namedShape(
              cluster,
              "iso_topo_" + key + "_valid");

      const double value =
          namedShape(
              cluster,
              "iso_topo_" + key);

      cone.valid =
          hasNamedShape(
              cluster,
              "iso_topo_" + key) &&
          finite(validFlag) &&
          validFlag > 0.5 &&
          finite(cone.axisEta) && finite(cone.axisPhi) &&
          finite(value);

      cone.coneSum =
          cone.valid
              ? value
              : kNaN;

      photon.isolation.push_back(cone);
    }
  }

  /*
   * PhotonRecord carries the isolation views while the normalized Isolation
   * table receives the same records for direct relational access.
   */
  for (const IsolationRecord& cone :
       photon.isolation)
  {
    m_isolationRecords.push_back(cone);
  }
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

void PhotonJetTree::capturePhotonIsolationConstituents(
    PHCompositeNode* topNode,
    const RawCluster* cluster,
    const PhotonRecord& photon)
{
  if (!m_config.output.writeIsolationConstituents)
  {
    return;
  }

  const double maxRadius =
      *std::max_element(
          m_config.photon.isolationRadii.begin(),
          m_config.photon.isolationRadii.end());

  const double axisEta =
      namedShape(
          cluster,
          "isolation_axis_eta");

  const double axisPhi =
      namedShape(
          cluster,
          "isolation_axis_phi");

  if (!finite(axisEta) ||
      !finite(axisPhi))
  {
    return;  // invalid nominal axis; keep the candidate and invalid cone state
  }

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
    struct Layer { const char* raw; const char* sub; const char* geometry; RawTowerDefs::CalorimeterId id; };

    const Layer layers[] = {
      {"TOWERINFO_CALIB_CEMC_RETOWER", "TOWERINFO_CALIB_CEMC_RETOWER_SUB1", "TOWERGEOM_HCALIN", RawTowerDefs::HCALIN},
      {"TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALIN_SUB1", "TOWERGEOM_HCALIN", RawTowerDefs::HCALIN},
      {"TOWERINFO_CALIB_HCALOUT", "TOWERINFO_CALIB_HCALOUT_SUB1", "TOWERGEOM_HCALOUT", RawTowerDefs::HCALOUT}};

    for (int layer = 0; layer < 3; ++layer)
    {
      const auto& node = layers[layer];

      auto* raw =
          findNode::getClass<TowerInfoContainer>(
              topNode,
              node.raw);

      auto* sub =
          findNode::getClass<TowerInfoContainer>(
              topNode,
              node.sub);

      auto* geometry =
          findNode::getClass<RawTowerGeomContainer>(
              topNode,
              node.geometry);

      if (!raw || !sub || !geometry) continue;
      const unsigned int channelCount = std::min(raw->size(), sub->size());

      for (unsigned int channel = 0;
           channel < channelCount;
           ++channel)
      {
        auto* rawTower =
            raw->get_tower_at_channel(channel);

        auto* subTower =
            sub->get_tower_at_channel(channel);

        if (!rawTower ||
            !subTower)
        {
          continue;  // historical witness census skips unavailable pairs
        }

        const auto key =
            TowerInfoDefs::encode_hcal(channel);

        const int etaBin =
            TowerInfoDefs::getCaloTowerEtaBin(key);

        const int phiBin =
            TowerInfoDefs::getCaloTowerPhiBin(key);

        auto* geom =
            geometry->get_tower_geometry(
                RawTowerDefs::encode_towerid(
                    node.id,
                    etaBin,
                    phiBin));

        IsolationConstituentRecord row;

        row.eventId = photon.eventId; row.photonId = photon.photonId;

        row.radius = maxRadius; row.source = "calorimeter_sub1";

        row.nativeKey = channel; row.subsystem = layer;

        row.rawEnergy = rawTower->get_energy(); row.subtractedEnergy = subTower->get_energy();

        row.maskState = (rawTower->get_isGood() ? 0U : 1U) | (subTower->get_isGood() ? 0U : 2U);

        row.masked = !subTower->get_isGood(); row.qualityState = subTower->get_isGood() ? 1 : 0;

        const double r = geom ? std::hypot(geom->get_center_x(), geom->get_center_y()) : 0;

        if (!geom || !(r > 0)) { row.qualityState = -1; m_isolationConstituents.push_back(row); continue; }

        /*
         * Recompute tower eta relative to the reconstructed collision vertex.
         * The geometry object itself is detector-centered.
         */
        const double eta =
            std::asinh(
                (std::sinh(geom->get_eta()) * r -
                 photon.producerVertexZ) /
                r);

        row.deltaEta = eta - axisEta; row.deltaPhi = wrapPhi(geom->get_phi() - axisPhi);

        row.deltaR = std::hypot(row.deltaEta, row.deltaPhi);

        if (!finite(row.deltaR) || row.deltaR >= maxRadius) continue;

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

  auto* topoclusters =
      findNode::getClass<RawClusterContainer>(
          topNode,
          m_config.photon.topoclusterNode);

  if (!topoclusters)
  {
    return;  // builder's explicit invalid topo state remains authoritative
  }

  const CLHEP::Hep3Vector vertex(
      0.0,
      0.0,
      photon.producerVertexZ);

  const auto range =
      topoclusters->getClusters();

  for (auto it = range.first;
       it != range.second;
       ++it)
  {
    const RawCluster* topo =
        it->second;

    if (!topo)
    {
      continue;
    }

    const double eta =
        RawClusterUtility::GetPseudorapidity(
            *topo,
            vertex);

    const double phi =
        RawClusterUtility::GetAzimuthAngle(
            *topo,
            vertex);

    const double energy =
        topo->get_energy();

    if (!finite(eta) ||
        !finite(phi) ||
        !finite(energy))
    {
      continue;
    }

    const double dEta =
        eta -
        axisEta;

    const double dPhi =
        wrapPhi(
            phi -
            axisPhi);

    const double dR =
        std::hypot(
            dEta,
            dPhi);

    if (!finite(dR) ||
        dR >= maxRadius)
    {
      continue;
    }

    IsolationConstituentRecord row;

    row.eventId = photon.eventId;
    row.photonId = photon.photonId;

    row.radius = maxRadius;
    row.source = "topocluster";

    row.nativeKey = it->first;
    row.subsystem = -1;

    row.deltaEta = dEta;
    row.deltaPhi = dPhi;
    row.deltaR = dR;

    row.rawEnergy = energy;

    row.transverseEnergy = energy / std::cosh(eta);  // signed, as the builder sums it

    row.qualityState =
        topo->isValid()
            ? 1
            : 0;

    /*
     * Historical candidate-removal witness: a topological cluster coincident
     * with the photon axis is identified as the candidate itself.
     */
    row.candidateRemoved =
        dR < 0.02;

    m_isolationConstituents.push_back(
        std::move(row));
  }
}
