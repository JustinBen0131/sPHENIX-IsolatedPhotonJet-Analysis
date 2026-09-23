//
// Photons.cc
//
// Reconstructed photon candidates.
//
// Owns
//   candidate identity, kinematics and the storage domain
//   cluster timing from the towers the cluster owns
//   every configured shower-shape definition, evaluated from the calibrated
//   tower grid with the same arithmetic the photon builder uses
//   the calorimeter cells that entered those definitions
//   reconstructed isolation in every configured cone, with its layer
//   components, as the builder attached it to the candidate
//   constituent-level isolation witnesses
//
// Does not own
//   photon builder registration                  -> Production.cc
//   truth association of a candidate             -> Truth.cc
//   identification scores and working points     -> downstream
//   ROOT tables                                  -> Output.cc
//
// Candidates are admitted on kinematics alone. The stored primitives are the
// complete rectangular shower family, the widths, the native shower
// fractions, the timing and the isolation cones, from which the supported
// identification models' proposed inputs are reconstructible downstream.
//
#include "../PhotonJetTree.h"

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

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

int wrapPhiIndex(int index)
{
  while (index < 0) index += kCemcPhiBins;
  while (index >= kCemcPhiBins) index -= kCemcPhiBins;
  return index;
}

// Signed azimuthal index difference in [-128, 127].
int deltaPhiIndex(const int index, const int reference)
{
  int delta = index - reference;
  if (delta > kCemcPhiBins / 2) delta -= kCemcPhiBins;
  if (delta < -kCemcPhiBins / 2) delta += kCemcPhiBins;
  return delta;
}

double wrapPhi(double phi)
{
  while (phi > kPi) phi -= 2.0 * kPi;
  while (phi <= -kPi) phi += 2.0 * kPi;
  return phi;
}

// A named shower-shape value the builder attached to the cluster, or NaN.
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

// ---- shower-shape definitions -------------------------------------------------
//
// A definition is a tower energy floor and the membership rules for sums and
// moments. "H" is the hybrid the identification models used: rectangular
// sums over the full good-tower grid, moments over cluster-owned towers.

struct ShowerDefinition
{
  std::string name;
  double floorGeV = 0.0;
  bool sumsOverFullGrid = true;
  bool momentsOverOwnedOnly = true;
};

ShowerDefinition resolveDefinition(const std::string& name)
{
  if (name == "H70") return {"H70", 0.070, true, true};
  if (name == "H0") return {"H0", 0.0, true, true};
  if (name == "G70") return {"G70", 0.070, true, false};
  if (name == "G0") return {"G0", 0.0, true, false};
  if (name == "O70") return {"O70", 0.070, false, true};
  if (name == "O0") return {"O0", 0.0, false, true};
  throw std::runtime_error("unknown shower-shape definition '" + name + "'");
}

// The native shower accessor of the cluster at one floor: the four shower
// fractions, the floating seed position, and the oriented 2x2 block.
struct NativeShape
{
  bool valid = false;
  std::array<double, 4> et{{kNaN, kNaN, kNaN, kNaN}};
  double rawCenterEta = kNaN;
  double rawCenterPhi = kNaN;
  double e22 = kNaN;
};

NativeShape nativeShape(const RawCluster* cluster, const double floorGeV)
{
  NativeShape result;
  const std::vector<float> values = cluster->get_shower_shapes(static_cast<float>(floorGeV));
  if (values.size() < 12)
  {
    return result;
  }
  for (std::size_t i = 0; i < 4; ++i)
  {
    result.et[i] = values[i];
  }
  result.rawCenterEta = values[4];
  result.rawCenterPhi = values[5];
  result.e22 = static_cast<double>(values[8]) + values[9] + values[10] + values[11];
  result.valid = finite(result.rawCenterEta) && finite(result.rawCenterPhi);
  return result;
}

// One calorimeter cell in the neighbourhood of a candidate, as read from the
// calibrated tower container and from the cluster's own tower map.
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
  if (TowerInfo* tower = towers->get_tower_at_key(key))
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

//
// Candidate capture
//
void PhotonJetTree::capturePhotons(PHCompositeNode* topNode)
{
  RawClusterContainer* candidates =
      findNode::getClass<RawClusterContainer>(topNode, m_config.photon.photonNode);
  if (!candidates)
  {
    throw std::runtime_error("missing photon candidate node '" + m_config.photon.photonNode + "'");
  }

  TowerInfoContainer* cemcTowers =
      findNode::getClass<TowerInfoContainer>(topNode, m_config.nodes.cemcTowers);
  if (!cemcTowers)
  {
    throw std::runtime_error("missing calibrated CEMC tower node '" + m_config.nodes.cemcTowers + "'");
  }

  // Reconstructed objects are defined with respect to the reconstructed
  // vertex. Without one, or outside the object vertex domain, the event is
  // retained with no candidates and the capture is marked so that no truth
  // miss can be claimed from it.
  m_event.recoObjectVertexInDomain =
      m_event.recoVertexValid && std::abs(m_event.recoVertexZ) < m_config.photon.objectVertexAbsZMaxCm;
  if (!m_event.recoVertexValid)
  {
    m_event.recoPhotonCaptureState = CaptureState::InputUnavailable;
    return;
  }
  if (!m_event.recoObjectVertexInDomain)
  {
    m_event.recoPhotonCaptureState = CaptureState::NotApplicable;
    return;
  }
  m_event.recoPhotonCaptureState = CaptureState::Complete;

  const CLHEP::Hep3Vector vertex(0.0, 0.0, m_event.recoVertexZ);
  std::uint32_t encountered = 0;

  for (const auto& [clusterKey, cluster] : candidates->getClustersMap())
  {
    const std::uint32_t ordinal = encountered++;
    if (!cluster)
    {
      ++m_event.recoPhotonUnclassifiableCount;
      continue;
    }

    PhotonRecord photon;
    photon.eventId = m_event.eventId;
    photon.nativeClusterKey = clusterKey;
    photon.producerVertexZ = m_event.recoVertexZ;

    photon.energy = cluster->get_energy();
    photon.et = RawClusterUtility::GetET(*cluster, vertex);
    photon.eta = RawClusterUtility::GetPseudorapidity(*cluster, vertex);
    photon.phi = RawClusterUtility::GetAzimuthAngle(*cluster, vertex);
    photon.kinematicsFinite =
        finite(photon.energy) && finite(photon.et) && finite(photon.eta) && finite(photon.phi);

    if (!photonInCaptureDomain(photon))
    {
      if (!photon.kinematicsFinite) ++m_event.recoPhotonUnclassifiableCount;
      continue;
    }

    photon.encounterOrdinal = ordinal;
    photon.photonId = photonjet::makePhotonIdentity(m_event.eventId, clusterKey, ordinal);

    capturePhotonTiming(cluster, cemcTowers, photon);
    capturePhotonShowerShapes(cluster, cemcTowers, photon);
    capturePhotonIsolation(cluster, photon);
    capturePhotonIsolationConstituents(topNode, cluster, photon);

    // The dominant-primary witness and the truth relation are filled by
    // Truth.cc after the truth inventory exists.
    photon.dominantTruth.state =
        isSimulation() ? DominantTruthState::EvaluatorUnavailable : DominantTruthState::NotSimulation;

    m_photons.push_back(std::move(photon));
  }
}

//
// Storage domain
//
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

//
// Timing
//
// Energy-weighted mean time over every tower the cluster owns, accumulated in
// single precision as the reconstruction does, with no quality gate and no
// floor. A non-positive denominator gives no measurement.
//
void PhotonJetTree::capturePhotonTiming(const RawCluster* cluster,
                                        TowerInfoContainer* cemcTowers,
                                        PhotonRecord& photon) const
{
  float numerator = 0.0F;
  float denominator = 0.0F;
  std::uint32_t contributing = 0;

  for (const auto& [rawKey, mapEnergy] : cluster->get_towermap())
  {
    (void) mapEnergy;
    const int etaIndex = RawTowerDefs::decode_index1(rawKey);
    const int phiIndex = wrapPhiIndex(RawTowerDefs::decode_index2(rawKey));
    if (etaIndex < 0 || etaIndex >= kCemcEtaBins)
    {
      continue;
    }
    TowerInfo* tower = cemcTowers->get_tower_at_key(TowerInfoDefs::encode_emcal(etaIndex, phiIndex));
    if (!tower)
    {
      continue;
    }
    numerator += tower->get_time() * tower->get_energy();
    denominator += tower->get_energy();
    ++contributing;
  }

  PhotonTimingRecord& timing = photon.timing;
  timing.energyWeightedNumerator = numerator;
  timing.energyDenominator = denominator;
  timing.contributingTowerCount = contributing;

  if (denominator > 0.0F)
  {
    timing.meanTimeSamples = numerator / denominator;
  }
  timing.finite = finite(timing.meanTimeSamples);
  timing.valid = denominator > 0.0F && timing.finite;
  timing.timeNs = timing.valid ? timing.meanTimeSamples * m_config.photon.timingSampleNs : kNaN;
}

//
// Shower shapes and cells
//
// For each configured definition the 7x7 grid about the seed tower is read
// from the calibrated towers, with the same membership, floor and
// arithmetic as the photon builder:
//
//   rectangular sums   good towers above the floor, over the full grid
//   moments            energy-weighted second moments about the sub-tower
//                      centre of gravity, over cluster-owned towers, with the
//                      central tower excluded from the numerator
//   strip widths       over the full grid, about the seed
//   shower fractions   the cluster's own accessor at the same floor
//
// The union of the grid cells of every definition and of every owned tower
// is written to the cell table when it is enabled.
//
void PhotonJetTree::capturePhotonShowerShapes(const RawCluster* cluster,
                                              TowerInfoContainer* cemcTowers,
                                              PhotonRecord& photon)
{
  CellMap cells;

  // Cluster-owned towers.
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

  // Lead tower, for the extent of the owned towers.
  const std::pair<int, int> lead = cluster->get_lead_tower();

  for (const std::string& definitionName : m_config.showerDefinitions)
  {
    const ShowerDefinition definition = resolveDefinition(definitionName);
    const NativeShape native = nativeShape(cluster, definition.floorGeV);

    ShowerShapeRecord shape;
    shape.definitionName = definition.name;
    shape.energyFloorGeV = definition.floorGeV;
    shape.et1 = native.et[0];
    shape.et2 = native.et[1];
    shape.et3 = native.et[2];
    shape.et4 = native.et[3];
    shape.e22 = native.e22;

    if (!native.valid)
    {
      photon.showerShapes.push_back(std::move(shape));
      continue;
    }

    // Seed tower and centre of gravity, as the builder derives them.
    const float averageEta = static_cast<float>(native.rawCenterEta) + 0.5F;
    const float averagePhi = static_cast<float>(native.rawCenterPhi) + 0.5F;
    const int centerEta = static_cast<int>(std::floor(averageEta));
    const int centerPhi = wrapPhiIndex(static_cast<int>(std::floor(averagePhi)));
    const float shiftEta = averageEta - std::floor(averageEta) - 0.5F;
    const float shiftPhi = averagePhi - std::floor(averagePhi) - 0.5F;
    const float cogEta = 3.0F + shiftEta;
    const float cogPhi = 3.0F + shiftPhi;
    const int signPhi = (averagePhi - std::floor(averagePhi)) > 0.5F ? 1 : -1;

    shape.centerEtaIndex = centerEta;
    shape.centerPhiIndex = centerPhi;
    shape.centerOfGravityEta = cogEta;
    shape.centerOfGravityPhi = cogPhi;
    shape.deltaEtaCog = std::abs(centerEta - averageEta);
    shape.deltaPhiCog = std::abs(centerPhi - averagePhi);
    shape.radialSpread = std::hypot(shape.deltaEtaCog, shape.deltaPhiCog);

    if (centerEta < 0 || centerEta >= kCemcEtaBins)
    {
      photon.showerShapes.push_back(std::move(shape));
      continue;
    }

    // The 7x7 grid: energies above the floor from good towers, and ownership.
    float grid[7][7] = {{0.0F}};
    bool owned[7][7] = {{false}};
    std::uint32_t ownedCount = 0, goodCount = 0, zeroCount = 0, negativeCount = 0, nonFiniteCount = 0;
    std::uint32_t saturated = 0;
    int deltaEtaMax = 0, deltaPhiMax = 0;

    for (int deta = -kGridHalfWidth; deta <= kGridHalfWidth; ++deta)
    {
      for (int dphi = -kGridHalfWidth; dphi <= kGridHalfWidth; ++dphi)
      {
        const int etaIndex = centerEta + deta;
        const int phiIndex = wrapPhiIndex(centerPhi + dphi);
        addCell(cells, cemcTowers, etaIndex, phiIndex);
        auto found = cells.find({etaIndex, phiIndex});
        if (found == cells.end())
        {
          continue;  // outside the pseudorapidity range of the calorimeter
        }
        Cell& cell = found->second;
        cell.inSupport = true;

        const int i = deta + kGridHalfWidth;
        const int j = dphi + kGridHalfWidth;
        owned[i][j] = cell.owned;
        if (cell.owned)
        {
          ++ownedCount;
          deltaEtaMax = std::max(deltaEtaMax, std::abs(etaIndex - lead.first));
          deltaPhiMax = std::max(deltaPhiMax, std::abs(deltaPhiIndex(phiIndex, lead.second)));
          if (cell.saturated) ++saturated;
        }
        if (!cell.present)
        {
          continue;
        }
        if (!finite(cell.energy))
        {
          ++nonFiniteCount;
          continue;
        }
        if (cell.good)
        {
          ++goodCount;
          if (cell.energy == 0.0) ++zeroCount;
          if (cell.energy < 0.0) ++negativeCount;
          if (cell.energy > definition.floorGeV)
          {
            grid[i][j] = static_cast<float>(cell.energy);
          }
        }
      }
    }

    // Rectangular sums, widths and moments, in single precision as the
    // builder accumulates them.
    float e11 = grid[3][3];
    float e13 = 0, e15 = 0, e17 = 0, e31 = 0, e51 = 0, e71 = 0;
    float e33 = 0, e35 = 0, e37 = 0, e53 = 0, e55 = 0, e57 = 0, e73 = 0, e75 = 0, e77 = 0;
    float e32 = 0, e52 = 0, e72 = 0, w32 = 0, w52 = 0, w72 = 0;
    float weta = 0, wphi = 0, wetaCog = 0, wphiCog = 0, wetaCogX = 0, wphiCogX = 0;
    float weta33 = 0, wphi33 = 0;
    float momentDenominator = 0, moment33Denominator = 0;

    for (int i = 0; i < 7; ++i)
    {
      for (int j = 0; j < 7; ++j)
      {
        const int di = std::abs(i - 3);
        const int dj = std::abs(j - 3);
        const float dEta = static_cast<float>(i) - cogEta;
        const float dPhi = static_cast<float>(j) - cogPhi;
        const float energy = grid[i][j];
        const bool sumMember = definition.sumsOverFullGrid || owned[i][j];
        const bool momentMember = !definition.momentsOverOwnedOnly || owned[i][j];

        if (momentMember)
        {
          weta += energy * di * di;
          wphi += energy * dj * dj;
          wetaCog += energy * dEta * dEta;
          wphiCog += energy * dPhi * dPhi;
          momentDenominator += energy;
          if (i != 3 || j != 3)
          {
            wetaCogX += energy * dEta * dEta;
            wphiCogX += energy * dPhi * dPhi;
          }
          if (di <= 1 && dj <= 1)
          {
            moment33Denominator += energy;
            if (i != 3 || j != 3)
            {
              weta33 += energy * dEta * dEta;
              wphi33 += energy * dPhi * dPhi;
            }
          }
        }

        if (!sumMember)
        {
          continue;
        }
        e77 += energy;
        if (di <= 1 && (dj == 0 || j == 3 + signPhi)) { w32 += energy * (i - 3) * (i - 3); e32 += energy; }
        if (di <= 2 && (dj == 0 || j == 3 + signPhi)) { w52 += energy * (i - 3) * (i - 3); e52 += energy; }
        if (di <= 3 && (dj == 0 || j == 3 + signPhi)) { w72 += energy * (i - 3) * (i - 3); e72 += energy; }
        if (di <= 0 && dj <= 1) e13 += energy;
        if (di <= 0 && dj <= 2) e15 += energy;
        if (di <= 0 && dj <= 3) e17 += energy;
        if (di <= 1 && dj <= 0) e31 += energy;
        if (di <= 2 && dj <= 0) e51 += energy;
        if (di <= 3 && dj <= 0) e71 += energy;
        if (di <= 1 && dj <= 1) e33 += energy;
        if (di <= 1 && dj <= 2) e35 += energy;
        if (di <= 1 && dj <= 3) e37 += energy;
        if (di <= 2 && dj <= 1) e53 += energy;
        if (di <= 3 && dj <= 1) e73 += energy;
        if (di <= 2 && dj <= 2) e55 += energy;
        if (di <= 2 && dj <= 3) e57 += energy;
        if (di <= 3 && dj <= 2) e75 += energy;
      }
    }

    if (momentDenominator > 0.0F)
    {
      weta /= momentDenominator;
      wphi /= momentDenominator;
      wetaCog /= momentDenominator;
      wphiCog /= momentDenominator;
      wetaCogX /= momentDenominator;
      wphiCogX /= momentDenominator;
    }
    else
    {
      weta = wphi = wetaCog = wphiCog = wetaCogX = wphiCogX = kNaN;
    }
    if (moment33Denominator > 0.0F)
    {
      weta33 /= moment33Denominator;
      wphi33 /= moment33Denominator;
    }
    else
    {
      weta33 = wphi33 = kNaN;
    }
    if (e32 > 0.0F) w32 /= e32; else w32 = kNaN;
    if (e52 > 0.0F) w52 /= e52; else w52 = kNaN;
    if (e72 > 0.0F) w72 /= e72; else w72 = kNaN;

    shape.e11 = e11; shape.e13 = e13; shape.e15 = e15; shape.e17 = e17;
    shape.e31 = e31; shape.e32 = e32; shape.e33 = e33; shape.e35 = e35; shape.e37 = e37;
    shape.e51 = e51; shape.e52 = e52; shape.e53 = e53; shape.e55 = e55; shape.e57 = e57;
    shape.e71 = e71; shape.e72 = e72; shape.e73 = e73; shape.e75 = e75; shape.e77 = e77;
    shape.weta = weta; shape.wphi = wphi;
    shape.wetaCog = wetaCog; shape.wphiCog = wphiCog;
    shape.wetaCogX = wetaCogX; shape.wphiCogX = wphiCogX;
    shape.weta33CogX = weta33; shape.wphi33CogX = wphi33;
    shape.w32 = w32; shape.w52 = w52; shape.w72 = w72;
    shape.deltaEtaMax = deltaEtaMax;
    shape.deltaPhiMax = deltaPhiMax;
    shape.saturatedTowerCount = saturated;
    shape.ownedCellCount = ownedCount;
    shape.goodCellCount = goodCount;
    shape.zeroCellCount = zeroCount;
    shape.negativeCellCount = negativeCount;
    shape.nonFiniteCellCount = nonFiniteCount;
    shape.valid = finite(shape.e33) && finite(shape.wetaCogX) && finite(shape.wphiCogX) &&
                  finite(shape.et1) && finite(shape.et2) && finite(shape.et3) && finite(shape.et4);

    photon.showerShapes.push_back(std::move(shape));
  }

  // ---- cell table ----------------------------------------------------------------

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
    row.calibratedEnergy = cell.present ? cell.energy : kNaN;
    row.clusterMapEnergy = cell.owned ? cell.clusterMapEnergy : kNaN;
    row.towerTimeSamples = cell.present ? cell.timeSamples : kNaN;
    row.towerStatus = cell.status;
    row.ownedByCluster = cell.owned;
    row.inShowerSupport = cell.inSupport;
    row.good = cell.present && cell.good;
    row.zeroEnergy = cell.present && finite(cell.energy) && cell.energy == 0.0;
    row.negativeEnergy = cell.present && finite(cell.energy) && cell.energy < 0.0;
    row.finite = cell.present && finite(cell.energy) && finite(cell.timeSamples);
    m_photonCells.push_back(std::move(row));
  }
}

//
// Isolation
//
// The builder attaches, per cone radius, the raw tower cone in each
// calorimeter layer with the candidate's own transverse energy already
// removed from the electromagnetic layer, the underlying-event subtracted
// cone in Au+Au, and the topological-cluster cone in p+p. Every available
// method is stored so the nominal choice stays a downstream decision.
//
void PhotonJetTree::capturePhotonIsolation(const RawCluster* cluster, PhotonRecord& photon)
{
  for (const double radius : m_config.photon.isolationRadii)
  {
    const int code = static_cast<int>(std::lround(radius * 100.0));
    const std::string key = code == 30 ? "03" : code == 40 ? "04" : code == 20 ? "02" : code == 10 ? "01" : "";
    if (key.empty())
    {
      continue;  // the builder attaches cones at these radii only
    }

    // Raw tower cone.
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
      cone.valid = finite(cone.electromagneticComponent) && finite(cone.innerHadronicComponent) &&
                   finite(cone.outerHadronicComponent);
      if (cone.valid)
      {
        cone.coneSum = cone.electromagneticComponent + cone.innerHadronicComponent + cone.outerHadronicComponent;
      }
      photon.isolation.push_back(cone);
    }

    // Underlying-event subtracted cone, Au+Au.
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
      // The builder writes a -999 sentinel when the subtracted towers were
      // absent; that is not a measurement.
      cone.valid = finite(cone.electromagneticComponent) && cone.electromagneticComponent > -900.0 &&
                   finite(cone.innerHadronicComponent) && cone.innerHadronicComponent > -900.0 &&
                   finite(cone.outerHadronicComponent) && cone.outerHadronicComponent > -900.0;
      if (cone.valid)
      {
        cone.coneSum = cone.electromagneticComponent + cone.innerHadronicComponent + cone.outerHadronicComponent;
      }
      else
      {
        cone.electromagneticComponent = cone.innerHadronicComponent = cone.outerHadronicComponent = kNaN;
      }
      photon.isolation.push_back(cone);
    }

    // Topological-cluster cone, p+p.
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
      const double validFlag = namedShape(cluster, "iso_topo_valid");
      const double value = namedShape(cluster, "iso_topo_" + key);
      cone.valid = hasNamedShape(cluster, "iso_topo_" + key) && finite(validFlag) && validFlag > 0.5 && finite(value);
      cone.coneSum = cone.valid ? value : kNaN;
      photon.isolation.push_back(cone);
    }
  }

  for (const IsolationRecord& cone : photon.isolation)
  {
    m_isolationRecords.push_back(cone);
  }
}

//
// Isolation constituents
//
// p+p: every topological cluster inside the largest configured cone, with
// its separation and signed transverse energy, so the cone sum can be
// rebuilt with a different radius. AuAu retains raw/SUB1 tower witnesses
// with their independent quality masks and the same cone axis as the builder.
//
void PhotonJetTree::capturePhotonIsolationConstituents(PHCompositeNode* topNode,
                                                       const RawCluster* cluster,
                                                       const PhotonRecord& photon)
{
  if (!m_config.output.writeIsolationConstituents)
  {
    return;
  }

  const double maxRadius = *std::max_element(m_config.photon.isolationRadii.begin(), m_config.photon.isolationRadii.end());
  const double axisEta = namedShape(cluster, "isolation_axis_eta");
  const double axisPhi = namedShape(cluster, "isolation_axis_phi");
  if (!finite(axisEta) || !finite(axisPhi))
    throw std::runtime_error("isolation constituent capture requires the builder's exact cone axis");

  if (isAuAu())
  {
    // RecoilJets_AuAu.cc:7970-8033: all three layers use the 24x64 grid.
    // SUB1 quality controls the subtraction sum; raw quality is independent.
    struct Layer { const char* raw; const char* sub; const char* geometry; RawTowerDefs::CalorimeterId id; };
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
      if (!raw || !sub || !geometry || raw->size() != sub->size())
        throw std::runtime_error("incomplete AuAu isolation constituent input");
      for (unsigned int channel = 0; channel < sub->size(); ++channel)
      {
        auto* rawTower = raw->get_tower_at_channel(channel);
        auto* subTower = sub->get_tower_at_channel(channel);
        if (!rawTower || !subTower) throw std::runtime_error("null AuAu isolation tower");
        const auto key = TowerInfoDefs::encode_hcal(channel);
        const int etaBin = TowerInfoDefs::getCaloTowerEtaBin(key);
        const int phiBin = TowerInfoDefs::getCaloTowerPhiBin(key);
        auto* geom = geometry->get_tower_geometry(RawTowerDefs::encode_towerid(node.id, etaBin, phiBin));
        IsolationConstituentRecord row;
        row.eventId = photon.eventId; row.photonId = photon.photonId;
        row.radius = maxRadius; row.source = "calorimeter_sub1";
        row.nativeKey = channel; row.subsystem = layer;
        row.rawEnergy = rawTower->get_energy(); row.subtractedEnergy = subTower->get_energy();
        row.maskState = (rawTower->get_isGood() ? 0U : 1U) | (subTower->get_isGood() ? 0U : 2U);
        row.masked = !subTower->get_isGood(); row.qualityState = subTower->get_isGood() ? 1 : 0;
        const double r = geom ? std::hypot(geom->get_center_x(), geom->get_center_y()) : 0;
        if (!geom || !(r > 0)) { row.qualityState = -1; m_isolationConstituents.push_back(row); continue; }
        const double eta = std::asinh((std::sinh(geom->get_eta()) * r - m_event.recoVertexZ) / r);
        row.deltaEta = eta - axisEta; row.deltaPhi = wrapPhi(geom->get_phi() - axisPhi);
        row.deltaR = std::hypot(row.deltaEta, row.deltaPhi);
        if (!finite(row.deltaR) || row.deltaR >= maxRadius) continue;
        row.transverseEnergy = row.rawEnergy / std::cosh(eta);
        row.subtractedTransverseEnergy = row.subtractedEnergy / std::cosh(eta);
        // Witness only: the nominal builder subtracts the photon ET once,
        // not a sum with this core masked. Preserve the historical flag.
        row.candidateRemoved = layer == 0 && row.deltaR < 0.02;
        m_isolationConstituents.push_back(std::move(row));
      }
    }
    return;
  }
  auto* topoclusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.topoclusterNode);
  if (!topoclusters) throw std::runtime_error("required topocluster isolation input is missing");

  const CLHEP::Hep3Vector vertex(0.0, 0.0, m_event.recoVertexZ);
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
    if (!finite(eta) || !finite(phi) || !finite(energy))
    {
      continue;
    }
    const double dEta = eta - axisEta;
    const double dPhi = wrapPhi(phi - axisPhi);
    const double dR = std::hypot(dEta, dPhi);
    if (!finite(dR) || dR >= maxRadius)
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
    row.qualityState = topo->isValid() ? 1 : 0;
    // The historical candidate-removal witness: the cluster coincident with
    // the candidate itself.
    row.candidateRemoved = dR < 0.02;
    m_isolationConstituents.push_back(std::move(row));
  }
}
