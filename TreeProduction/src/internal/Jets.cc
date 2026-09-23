//
// Jets.cc
//
// Reconstructed jets and photon-jet relationships.
//
// Owns
//   reading the configured jet collections, raw and corrected
//   the raw-to-corrected binding and the calibration state of every jet
//   the measured jet area
//   deterministic retained order within each view and radius
//   the photon x jet relationships with their azimuthal separation and
//   momentum balance
//
// Does not own
//   jet clustering, underlying-event subtraction, the jet energy scale
//                                                     -> Production.cc
//   truth jets and reconstruction-to-truth matching   -> Truth.cc
//   ROOT tables                                       -> Output.cc
//
// The jet energy scale is applied exactly once, upstream, from the raw node
// into the corrected node. This file never recalibrates; it reads both nodes
// so the tree keeps the raw quantity, the corrected quantity, the collection
// identity and the calibration state side by side.
//
// Ordering. deterministic_order is the encounter order in the corrected
// container after invalid calibrated jets are skipped, restarting at zero in
// every view and radius. It is not a momentum ranking, and downstream code
// must not read it as one.
//
#include "../PhotonJetTree.h"

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPi = 3.14159265358979323846;
constexpr unsigned int kInvalidJetId = std::numeric_limits<unsigned int>::max();

// Bits of JetRecord::qualityBitmask. They describe how the stored row was
// assembled; they are not analysis cuts.
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

bool finite(const double value)
{
  return std::isfinite(value);
}

// |phi1 - phi2| wrapped into [0, pi].
double absoluteDeltaPhi(const double phi1, const double phi2)
{
  double delta = phi1 - phi2;
  while (delta > kPi) delta -= 2.0 * kPi;
  while (delta <= -kPi) delta += 2.0 * kPi;
  return std::abs(delta);
}

// The FastJet area. Jetv2 stores properties at container-defined indices, so
// the semantic property has to be translated through the container first.
double jetArea(JetContainer* container, Jet* jet)
{
  if (container->has_property(Jet::PROPERTY::prop_area))
  {
    const Jet::PROPERTY index = container->property_index(Jet::PROPERTY::prop_area);
    const float area = jet->get_property(index);
    return finite(area) ? area : kNaN;
  }
  if (jet->has_property(Jet::PROPERTY::prop_area))
  {
    const float area = jet->get_property(Jet::PROPERTY::prop_area);
    return finite(area) ? area : kNaN;
  }
  return kNaN;
}

// A missing stored radius is tolerated because some historical containers do
// not persist it; a stored radius that disagrees with the configuration is
// a production error.
bool radiusCompatible(JetContainer* container, const double expected)
{
  const float stored = container->get_jetpar_R();
  if (!finite(stored) || stored <= 0.0F)
  {
    return true;
  }
  return std::abs(static_cast<double>(stored) - expected) < 1.0e-4;
}

}  // namespace

//
// All configured collections
//
void PhotonJetTree::captureJets(PHCompositeNode* topNode)
{
  for (const JetNodeConfig& collection : m_config.jets.nodes)
  {
    m_event.recoJetView.push_back(static_cast<int>(collection.view));
    m_event.recoJetRadiusCode.push_back(static_cast<int>(photonjet::radiusCode(collection.radius)));
    m_event.recoJetCaptureState.push_back(static_cast<int>(CaptureState::NotApplicable));
    m_event.recoJetCount.push_back(0);
    if (!m_event.recoObjectVertexInDomain) continue;
    m_event.recoJetCaptureState.back() = static_cast<int>(CaptureState::InputUnavailable);
    const auto before = m_jets.size();
    captureJetCollection(topNode, collection);
    m_event.recoJetCaptureState.back() = static_cast<int>(CaptureState::Complete);
    m_event.recoJetCount.back() = static_cast<std::uint32_t>(m_jets.size() - before);
  }
}

//
// One collection
//
// The calibration creates corrected jets in the order of the raw container
// and sets each corrected jet's id to that ordinal, so corrected ordinal i
// is raw ordinal i. Both containers are read; a size mismatch means the
// upstream chain did not run as configured and stops the run.
//
void PhotonJetTree::captureJetCollection(PHCompositeNode* topNode,
                                         const JetNodeConfig& collection)
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

  const bool sizesMatch = raw->size() == corrected->size();
  if (!sizesMatch)
  {
    throw std::runtime_error("raw and corrected jet multiplicity differ for '" +
                             collection.rawNode + "' and '" + collection.correctedNode + "'");
  }

  const bool radiusMatches =
      radiusCompatible(corrected, collection.radius) && radiusCompatible(raw, collection.radius);
  if (!radiusMatches)
  {
    throw std::runtime_error("configured radius disagrees with the container metadata of '" +
                             collection.correctedNode + "'");
  }

  std::uint32_t retainedOrder = 0;

  for (unsigned int ordinal = 0; ordinal < corrected->size(); ++ordinal)
  {
    Jet* correctedJet = corrected->get_jet(ordinal);
    Jet* rawJet = raw->get_jet(ordinal);

    // Invalid calibrated objects are skipped before raw binding, as in the
    // reference capture. Retained order counts only the surviving objects.
    if (!correctedJet || !finite(correctedJet->get_pt()) || correctedJet->get_pt() < 0) continue;
    if (!rawJet)
    {
      throw std::runtime_error("null jet at ordinal " + std::to_string(ordinal) + " in '" +
                               collection.correctedNode + "'");
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

    // The corrected jet defines the stored momentum; the raw transverse
    // momentum is kept independently. Jet mass is not part of the record.
    jet.rawPt = rawJet->get_pt();
    jet.correctedPt = correctedJet->get_pt();
    jet.eta = correctedJet->get_eta();
    jet.phi = correctedJet->get_phi();

    // The calibration does not copy the area to the corrected jet.
    jet.area = jetArea(raw, rawJet);

    // ---- calibration state -------------------------------------------------

    const bool calibratedFlag = correctedJet->get_isCalib() != 0;
    const unsigned int correctedId = correctedJet->get_id();
    const unsigned int rawId = rawJet->get_id();
    const bool correctedIdMatches = correctedId == ordinal;
    const bool rawIdCompatible = rawId == kInvalidJetId || rawId == ordinal || rawId == correctedId;
    const bool rawPtNonNegative = finite(jet.rawPt) && jet.rawPt >= 0.0;

    jet.calibrationValid = rawPtNonNegative && finite(jet.correctedPt) && finite(jet.eta) &&
                           finite(jet.phi) && calibratedFlag && correctedIdMatches;
    // Raw object ids need not equal their container ordinals. JetCalib binds
    // its output id to the raw ordinal; retain the native raw id as evidence.
    if (!rawPtNonNegative || !calibratedFlag || !correctedIdMatches)
      throw std::runtime_error("invalid calibrated/raw jet binding: " + collection.correctedNode);
    if (!finite(jet.area) || jet.area < 0.0)
      throw std::runtime_error("missing or invalid FastJet active area: " + collection.rawNode);

    std::uint64_t quality = 0;
    if (finite(jet.correctedPt)) quality |= kFiniteCorrectedPt;
    if (finite(jet.rawPt)) quality |= kFiniteRawPt;
    if (finite(jet.eta)) quality |= kFiniteEta;
    if (finite(jet.phi)) quality |= kFinitePhi;
    if (finite(jet.area)) quality |= kFiniteArea;
    if (calibratedFlag) quality |= kCalibratedFlagSet;
    if (correctedIdMatches) quality |= kCorrectedIdMatchesOrdinal;
    if (rawIdCompatible) quality |= kRawIdCompatible;
    if (sizesMatch) quality |= kContainerSizesMatch;
    if (radiusMatches) quality |= kContainerRadiusMatches;
    if (rawPtNonNegative) quality |= kRawPtNonNegative;

    // ---- storage domain ----------------------------------------------------
    //
    // A jet with non-finite corrected kinematics is not a retained physics
    // object. Retain all finite nonnegative calibrated jets; pT and eta
    // measurement cuts belong to TreeToHists, not this inventory.

    if (!jetInCaptureDomain(jet))
    {
      continue;
    }
    quality |= kInCaptureDomain;
    jet.qualityBitmask = quality;

    // The retained order is the encounter order of the jets that survived,
    // restarting in every collection.
    jet.deterministicOrder = retainedOrder;
    jet.jetId = photonjet::makeJetIdentity(m_event.eventId, jet.view, jet.radius,
                                           jet.nativeKey, jet.deterministicOrder);
    ++retainedOrder;

    m_jets.push_back(std::move(jet));
  }
}

//
// Storage domain
//
bool PhotonJetTree::jetInCaptureDomain(const JetRecord& jet) const
{
  if (!finite(jet.correctedPt) || !finite(jet.eta) || !finite(jet.phi))
  {
    return false;
  }
  return true;
}

//
// Photon x jet relationships
//
// Every retained photon is paired with every retained jet. Nothing is
// dropped for failing the away-side requirement; the recoil witness records
// it and the downstream analysis reproduces any requirement from the stored
// separation.
//
void PhotonJetTree::buildPhotonJetPairs()
{
  m_photonJetPairs.clear();

  if (!m_config.output.writePhotonJetPairs || m_photons.empty() || m_jets.empty())
  {
    return;
  }

  m_photonJetPairs.reserve(m_photons.size() * m_jets.size());

  for (std::size_t photonIndex = 0; photonIndex < m_photons.size(); ++photonIndex)
  {
    const PhotonRecord& photon = m_photons[photonIndex];

    for (const JetRecord& jet : m_jets)
    {
      PhotonJetPairRecord pair;
      pair.eventId = m_event.eventId;
      pair.photonId = photon.photonId;
      pair.jetId = jet.jetId;
      pair.jetView = jet.view;
      pair.jetRadius = jet.radius;

      pair.photonRank = 0;  // reference pair contract, independent of candidate index
      pair.jetRank = jet.deterministicOrder;

      if (finite(photon.phi) && finite(jet.phi))
      {
        pair.deltaPhi = absoluteDeltaPhi(jet.phi, photon.phi);
      }

      // Momentum balance: corrected jet momentum over photon transverse
      // energy, NaN unless the denominator is positive.
      if (finite(jet.correctedPt) && finite(photon.et) && photon.et > 0.0)
      {
        pair.xJGamma = jet.correctedPt / photon.et;
      }

      // Inclusive comparison: at least 7 pi / 8.
      pair.recoilWitness = finite(pair.deltaPhi) && pair.deltaPhi >= photonjet::kRecoilDeltaPhiMin;

      pair.pairId = photonjet::makePairIdentity(pair.photonId, pair.jetId, pair.jetRank);

      m_photonJetPairs.push_back(std::move(pair));
    }
  }
}
