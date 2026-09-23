/**
 * @file Jets.cc
 * @brief Reconstructed jet capture and photon-jet relationship construction.
 *
 * This file owns the reconstructed-jet information written by PhotonJetTree:
 *
 *   configured raw and calibrated jet collections
 *     -> raw-to-corrected jet correspondence
 *     -> calibration and quality witnesses
 *     -> FastJet active area
 *     -> deterministic retained jet order
 *     -> photon x jet relationships
 *     -> delta-phi, xJgamma, and recoil witnesses
 *
 * The reconstruction itself is deliberately upstream of this file.
 *
 *   Production.cc
 *       jet clustering, Au+Au underlying-event subtraction, and JES application
 *
 *   Truth.cc
 *       truth jets and reconstructed-to-truth jet associations
 *
 *   Output.cc
 *       ROOT booking and serialization
 *
 * Jet energy calibration is applied exactly once before PhotonJetTree reads the collections. Jets.cc never recalibrates a jet. It captures both the raw and
 * corrected quantities so the stored record makes the calibration relationship explicit and auditable.
 *
 * Ordering
 * --------
 * JetRecord::deterministicOrder is the encounter order of retained jets in the
 * corrected container after invalid calibrated objects have been skipped. It
 * restarts at zero for every configured view/radius. It is an identity/order
 * witness, not a transverse-momentum ranking.
 *
 * Analysis selections such as jet-pT thresholds, eta windows, and away-side
 * requirements are not applied to the jet inventory here. The retained records
 * preserve the inputs needed to make those selections downstream.
 */

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

bool finite(const double value)
{
  return std::isfinite(value);
}

// ============================================================================
// Small jet-geometry helpers
// ============================================================================

/**
 * Return |phi1 - phi2| wrapped onto [0, pi].
 */
double absoluteDeltaPhi(const double phi1, const double phi2)
{
  double delta = phi1 - phi2;

  while (delta > kPi) delta -= 2.0 * kPi;
  while (delta <= -kPi) delta += 2.0 * kPi;

  return std::abs(delta);
}


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
    const Jet::PROPERTY index =
        container->property_index(Jet::PROPERTY::prop_area);

    const float area = jet->get_property(index);

    return finite(area) ? area : kNaN;
  }

  if (jet->has_property(Jet::PROPERTY::prop_area))
  {
    const float area =
        jet->get_property(Jet::PROPERTY::prop_area);

    return finite(area) ? area : kNaN;
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

  if (!finite(stored) || stored <= 0.0F)
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
  for (const JetNodeConfig& collection : m_config.jets.nodes)
  {
    m_event.recoJetView.push_back(
        static_cast<int>(collection.view));

    m_event.recoJetRadiusCode.push_back(
        static_cast<int>(photonjet::radiusCode(collection.radius)));

    m_event.recoJetCaptureState.push_back(
        static_cast<int>(CaptureState::NotApplicable));

    m_event.recoJetCount.push_back(0);

    if (!m_event.recoObjectVertexInDomain) continue;

    /*
     * Once the event is in the reconstructed-object domain, failure to obtain
     * the configured collection is an unavailable-input state until the
     * collection capture succeeds.
     */
    m_event.recoJetCaptureState.back() =
        static_cast<int>(CaptureState::InputUnavailable);

    const auto before = m_jets.size();

    captureJetCollection(topNode, collection);

    m_event.recoJetCaptureState.back() =
        static_cast<int>(CaptureState::Complete);

    m_event.recoJetCount.back() =
        static_cast<std::uint32_t>(m_jets.size() - before);
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
void PhotonJetTree::captureJetCollection(PHCompositeNode* topNode,
                                         const JetNodeConfig& collection)
{
  JetContainer* corrected =
      findNode::getClass<JetContainer>(
          topNode,
          collection.correctedNode);

  JetContainer* raw =
      findNode::getClass<JetContainer>(
          topNode,
          collection.rawNode);

  if (!corrected)
  {
    throw std::runtime_error(
        "missing corrected jet node '" +
        collection.correctedNode +
        "'");
  }

  if (!raw)
  {
    throw std::runtime_error(
        "missing raw jet node '" +
        collection.rawNode +
        "'");
  }

  /*
   * Raw and calibrated collections are expected to have one-to-one ordinal
   * correspondence. A size mismatch means that relationship cannot be
   * established safely.
   */
  const bool sizesMatch =
      raw->size() == corrected->size();

  if (!sizesMatch)
  {
    throw std::runtime_error(
        "raw and corrected jet multiplicity differ for '" +
        collection.rawNode +
        "' and '" +
        collection.correctedNode +
        "'");
  }

  const bool radiusMatches =
      radiusCompatible(corrected, collection.radius) &&
      radiusCompatible(raw, collection.radius);

  if (!radiusMatches)
  {
    throw std::runtime_error(
        "configured radius disagrees with the container metadata of '" +
        collection.correctedNode +
        "'");
  }

  /*
   * deterministicOrder counts retained objects only and restarts for every
   * configured collection. It intentionally differs from both native jet id
   * and input-container ordinal.
   */
  std::uint32_t retainedOrder = 0;

  for (unsigned int ordinal = 0; ordinal < corrected->size(); ++ordinal)
  {
    Jet* correctedJet = corrected->get_jet(ordinal);
    Jet* rawJet = raw->get_jet(ordinal);

    /*
     * Match the reference capture semantics: invalid calibrated objects are
     * removed before a raw/calibrated row is constructed. Retained ordering
     * therefore counts only surviving corrected jets.
     */
    if (!correctedJet ||
        !finite(correctedJet->get_pt()) ||
        correctedJet->get_pt() < 0)
    {
      continue;
    }

    if (!rawJet)
    {
      throw std::runtime_error(
          "null jet at ordinal " +
          std::to_string(ordinal) +
          " in '" +
          collection.correctedNode +
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
     * The corrected collection defines the stored reconstructed four-vector
     * quantities. Raw pT is retained independently so the JES transformation
     * remains inspectable. Jet mass is intentionally not part of this record.
     */
    jet.rawPt = rawJet->get_pt();
    jet.correctedPt = correctedJet->get_pt();

    jet.eta = correctedJet->get_eta();
    jet.phi = correctedJet->get_phi();

    /*
     * JetCalib does not propagate the FastJet area into the corrected object,
     * so the active area is read from the corresponding raw jet.
     */
    jet.area = jetArea(raw, rawJet);

    // ------------------------------------------------------------------------
    // Calibration and raw/corrected binding
    // ------------------------------------------------------------------------

    const bool calibratedFlag =
        correctedJet->get_isCalib() != 0;

    const unsigned int correctedId =
        correctedJet->get_id();

    const unsigned int rawId =
        rawJet->get_id();

    const bool correctedIdMatches =
        correctedId == ordinal;

    /*
     * Raw native ids are retained as evidence but historically need not equal
     * their container ordinal. The corrected id is the authoritative ordinal
     * binding produced by JetCalib.
     */
    const bool rawIdCompatible =
        rawId == kInvalidJetId ||
        rawId == ordinal ||
        rawId == correctedId;

    const bool rawPtNonNegative =
        finite(jet.rawPt) &&
        jet.rawPt >= 0.0;

    jet.calibrationValid =
        rawPtNonNegative &&
        finite(jet.correctedPt) &&
        finite(jet.eta) &&
        finite(jet.phi) &&
        calibratedFlag &&
        correctedIdMatches;

    /*
     * These are production-contract requirements, not analysis-quality cuts.
     * A broken raw/calibrated binding cannot be represented as an ordinary jet
     * without losing the meaning of the stored corrected momentum.
     */
    if (!rawPtNonNegative ||
        !calibratedFlag ||
        !correctedIdMatches)
    {
      throw std::runtime_error(
          "invalid calibrated/raw jet binding: " +
          collection.correctedNode);
    }

    if (!finite(jet.area) || jet.area < 0.0)
    {
      throw std::runtime_error(
          "missing or invalid FastJet active area: " +
          collection.rawNode);
    }

    // ------------------------------------------------------------------------
    // Persisted quality witnesses
    // ------------------------------------------------------------------------

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

    jet.jetId = photonjet::makeJetIdentity(
        m_event.eventId,
        jet.view,
        jet.radius,
        jet.nativeKey,
        jet.deterministicOrder);

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
  if (!finite(jet.correctedPt) ||
      !finite(jet.eta) ||
      !finite(jet.phi))
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
// ============================================================================

void PhotonJetTree::buildPhotonJetPairs()
{
  m_photonJetPairs.clear();

  if (!m_config.output.writePhotonJetPairs ||
      m_photons.empty() ||
      m_jets.empty())
  {
    return;
  }

  m_photonJetPairs.reserve(
      m_photons.size() *
      m_jets.size());

  for (std::size_t photonIndex = 0;
       photonIndex < m_photons.size();
       ++photonIndex)
  {
    const PhotonRecord& photon =
        m_photons[photonIndex];

    for (const JetRecord& jet : m_jets)
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

      if (finite(photon.phi) && finite(jet.phi))
      {
        pair.deltaPhi =
            absoluteDeltaPhi(
                jet.phi,
                photon.phi);
      }

      /*
       * xJgamma uses the calibrated reconstructed jet momentum divided by the
       * reconstructed photon transverse energy. A non-positive or invalid
       * photon denominator leaves the quantity undefined.
       */
      if (finite(jet.correctedPt) &&
          finite(photon.et) &&
          photon.et > 0.0)
      {
        pair.xJGamma =
            jet.correctedPt /
            photon.et;
      }

      /*
       * Keep the recoil decision as a witness instead of a retention cut.
       * The boundary is inclusive: |delta phi| >= 7 pi / 8.
       */
      pair.recoilWitness =
          finite(pair.deltaPhi) &&
          pair.deltaPhi >= photonjet::kRecoilDeltaPhiMin;

      pair.pairId =
          photonjet::makePairIdentity(
              pair.photonId,
              pair.jetId,
              pair.jetRank);

      m_photonJetPairs.push_back(
          std::move(pair));
    }
  }
}
