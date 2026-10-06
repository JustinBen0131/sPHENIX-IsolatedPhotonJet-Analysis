/**
 * @file Truth.cc
 * Capture simulation truth and reconstruction-to-truth associations, and the
 * established producer's Au+Au embedded truth diagnostics (the truth-isolation
 * audit and the truth-photon misses). Retain source identity, ancestry, truth
 * isolation and matching validity. Data has no invented truth; absent and
 * ambiguous matches remain explicit.
 */

#include "../PhotonJetTree.h"
#include "Production.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4eval/CaloRawClusterEval.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <TLorentzVector.h>
#include <TVector2.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPi = 3.14159265358979323846;
constexpr double kRadiusTolerance = 1.0e-6;

constexpr double kTruthIsolationR03 = 0.30;
constexpr double kTruthIsolationR04 = 0.40;

/**
 * Embedding identifier of the hard-scatter subevent in the archived
 * double-interaction simulation.
 */
constexpr int kDoubleInteractionHardEmbedding = 2;

// The unsplit cluster node of the established Fig. 8 response targets; the
// builder's input clusters serve when a DST has none.
constexpr const char* kPPG12ResponseClusterNode = "CLUSTERINFO_CEMC_NO_SPLIT";

bool finiteValue(const double value) { return std::isfinite(value); }


// ============================================================================
// Particle-kinematic helpers
// ============================================================================

/**
 * Read transverse momentum, pseudorapidity, and azimuth from a Geant4 particle,
 * in the established producer's arithmetic (TLorentzVector).
 *
 * Invalid values are returned as NaN and the function returns false whenever
 * the momentum is non-finite or has non-positive transverse momentum.
 */
bool kinematics(const PHG4Particle* particle, double& pt, double& eta, double& phi)
{
  pt = eta = phi = kNaN;

  const TLorentzVector momentum(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e());

  if (!finiteValue(momentum.Pt()) || !finiteValue(momentum.Eta()) || !finiteValue(momentum.Phi()) ||
      !(momentum.Pt() > 0.0))
  {
    return false;
  }

  pt = momentum.Pt();
  eta = momentum.Eta();
  phi = momentum.Phi();

  return true;
}


// ============================================================================
// Geant4 -> HepMC generator association
//
// Geant4 track identity and HepMC generator identity are distinct. Association
// proceeds through the embedding-id-specific generator event and then through
// the Geant4 particle barcode. Failure at either step remains an explicit
// missing generator association.
// ============================================================================

/**
 * The archived AuAu photon stitching witness. Preserve its event-map priority
 * (0, then 1, then the first entry), inclusive |eta| <= 1.5 and parent-PDG
 * requirement. This is independent of the Geant4 analysis photon inventory.
 */
double embeddedFilterPhotonPt(PHCompositeNode* topNode, const std::string& node)
{
  auto* events = findNode::getClass<PHHepMCGenEventMap>(topNode, node);
  PHHepMCGenEvent* wrapper = nullptr;
  if (events)
  {
    wrapper = events->get(0);
    if (!wrapper) wrapper = events->get(1);
    if (!wrapper && !events->empty()) wrapper = events->begin()->second;
  }
  const auto* event = wrapper ? wrapper->getEvent() : nullptr;
  if (!event) return -1.0;

  double result = -1.0;
  for (auto it = event->particles_begin(); it != event->particles_end(); ++it)
  {
    const auto* particle = *it;
    if (!particle || particle->pdg_id() != 22) continue;
    const double pt = std::hypot(particle->momentum().px(), particle->momentum().py());
    const double eta = particle->momentum().pseudoRapidity();
    if (!std::isfinite(pt) || !std::isfinite(eta) || pt <= 0.0 || eta < -1.5 || eta > 1.5) continue;
    const auto* vertex = particle->production_vertex();
    if (!vertex) continue;
    bool parentOk = false;
    for (auto parent = vertex->particles_in_const_begin(); parent != vertex->particles_in_const_end(); ++parent)
    {
      if (!*parent) continue;
      const int pdg = std::abs((*parent)->pdg_id());
      if (pdg >= 1 && pdg <= 22)
      {
        parentOk = true;
        break;
      }
    }
    if (parentOk && pt > result) result = pt;
  }
  return result;
}

const HepMC::GenEvent* generatorEvent(const PHHepMCGenEventMap* events, const int embeddingId)
{
  if (!events)
  {
    return nullptr;
  }

  const PHHepMCGenEvent* wrapper = events->get(embeddingId);

  if (!wrapper)
  {
    return nullptr;
  }

  if (wrapper->get_embedding_id() != embeddingId)
  {
    throw std::runtime_error("generator event map key disagrees with the event's embedding id");
  }

  return wrapper->getEvent();
}


const HepMC::GenParticle* generatorParticle(const HepMC::GenEvent* event, const int barcode)
{
  if (!event)
  {
    return nullptr;
  }

  for (auto it = event->particles_begin(); it != event->particles_end(); ++it)
  {
    if (*it && (*it)->barcode() == barcode)
    {
      return *it;
    }
  }

  return nullptr;
}


// ============================================================================
// Prompt-photon classification
//
// The classifier follows the photon backward through single-photon
// continuations until it reaches the first meaningful production vertex.
//
// Classification
// --------------
//   Direct
//     two incoming and two outgoing legs, all with |PDG| <= 22
//
//   Fragmentation
//     one incoming quark/lepton with |PDG| <= 11 reappears beside the photon
//
//   Hadronic
//     one incoming particle with |PDG| > 37
//
//   Unclassified
//     a usable production vertex exists but none of the above patterns match
//
//   Unknown
//     no usable photon/history is available
//
// Generator status codes do not enter the classification. This reproduces the
// predecessor classifier, with only a guard against malformed cyclic histories.
// ============================================================================

int promptClassCode(const HepMC::GenParticle* pho)
{
  if (!pho || pho->pdg_id() != 22) return -1;

  const HepMC::GenVertex* vertex = pho->production_vertex();
  if (!vertex) return -1;

  std::vector<const HepMC::GenParticle*> incomingParticles;

  for (auto inItr = vertex->particles_in_const_begin(); inItr != vertex->particles_in_const_end(); ++inItr)
  {
    if (*inItr) incomingParticles.push_back(*inItr);
  }

  std::set<const HepMC::GenVertex*> seen{vertex};

  while (incomingParticles.size() == 1 && incomingParticles[0] && incomingParticles[0]->pdg_id() == 22)
  {
    vertex = incomingParticles[0]->production_vertex();

    if (!vertex || !seen.insert(vertex).second) return -1;

    incomingParticles.clear();

    for (auto inItr = vertex->particles_in_const_begin(); inItr != vertex->particles_in_const_end(); ++inItr)
    {
      if (*inItr) incomingParticles.push_back(*inItr);
    }
  }

  std::vector<const HepMC::GenParticle*> outgoingParticles;

  for (auto outItr = vertex->particles_out_const_begin(); outItr != vertex->particles_out_const_end(); ++outItr)
  {
    if (*outItr) outgoingParticles.push_back(*outItr);
  }

  bool hasOutgoingPhoton = false;

  for (const auto* outgoing : outgoingParticles)
  {
    if (outgoing && outgoing->pdg_id() == 22)
    {
      hasOutgoingPhoton = true;
      break;
    }
  }

  if (!hasOutgoingPhoton) return -1;

  if (incomingParticles.size() == 2 && outgoingParticles.size() == 2)
  {
    const int in0 = incomingParticles[0] ? incomingParticles[0]->pdg_id() : 0;

    const int in1 = incomingParticles[1] ? incomingParticles[1]->pdg_id() : 0;

    const int out0 = outgoingParticles[0] ? outgoingParticles[0]->pdg_id() : 0;

    const int out1 = outgoingParticles[1] ? outgoingParticles[1]->pdg_id() : 0;

    if (std::abs(in0) <= 22 && std::abs(in1) <= 22 && std::abs(out0) <= 22 && std::abs(out1) <= 22)
    {
      return 1;
    }
  }
  else if (incomingParticles.size() == 1 && incomingParticles[0])
  {
    const int inPid = incomingParticles[0]->pdg_id();

    if (std::abs(inPid) <= 11 && outgoingParticles.size() == 2)
    {
      for (const auto* outgoing : outgoingParticles)
      {
        if (outgoing && outgoing->pdg_id() == inPid) return 2;
      }
    }

    if (std::abs(inPid) > 37) return 3;
  }

  return 0;
}


photonjet::TruthPhotonClass promptClass(const HepMC::GenParticle* photon)
{
  return static_cast<photonjet::TruthPhotonClass>(promptClassCode(photon));
}


// ============================================================================
// Reconstructed-to-truth jet matching helpers
//
// Candidate edges are built geometrically. Greedy one-to-one assignment then
// follows a deterministic ordering so the same event always produces the same
// selected matches independent of container iteration accidents.
// ============================================================================

// The separation of an associated photon pair, hypot(delta eta, delta phi)
// with delta phi reduced to [-pi, pi]. It is evaluated with the algorithm of
// CPython 3.10 math.hypot, the arithmetic of the established association
// repair, so the stored diagnostic separation reproduces bit for bit; a
// platform std::hypot can differ in the last place. Inputs are differences.
double associationSeparation(const double deltaEta, const double deltaPhi)
{
  const double values[2] = {std::fabs(deltaEta), std::fabs(deltaPhi)};
  double max = 0.0;
  bool foundNaN = false;
  for (const double value : values)
  {
    foundNaN |= std::isnan(value);
    if (value > max) max = value;
  }
  if (std::isinf(max)) return max;
  if (foundNaN) return std::numeric_limits<double>::quiet_NaN();
  if (max == 0.0) return max;

  // Sum the squares exactly (Dekker splitting), then correct the root once.
  constexpr double kSplit = 134217729.0;  // 2^27 + 1
  double sum = 1.0, fraction1 = 0.0, fraction2 = 0.0, fraction3 = 0.0;
  const auto accumulate = [&sum](const double term, double& fraction)
  {
    const double previous = sum;
    sum += term;
    fraction += (previous - sum) + term;
  };
  int maxExponent = 0;
  std::frexp(max, &maxExponent);
  if (maxExponent < -1023)
  {
    for (double value : values)
    {
      value /= max;
      accumulate(value * value, fraction1);
    }
    return max * std::sqrt(sum - 1.0 + fraction1);
  }
  const double scale = std::ldexp(1.0, -maxExponent);
  for (double value : values)
  {
    value *= scale;
    const double t = value * kSplit;
    const double high = t - (t - value);
    const double low = value - high;
    accumulate(high * high, fraction1);
    accumulate(2.0 * high * low, fraction2);
    fraction3 += low * low;
  }
  const double root = std::sqrt(sum - 1.0 + (fraction1 + fraction2 + fraction3));
  const double t = root * kSplit;
  const double high = t - (t - root);
  const double low = root - high;
  accumulate(-high * high, fraction1);
  accumulate(-2.0 * high * low, fraction2);
  accumulate(-low * low, fraction3);
  const double residual = sum - 1.0 + (fraction1 + fraction2 + fraction3);
  return (root + residual / (2.0 * root)) / scale;
}

// The established producer's photon links (plink_*). Each reconstructed photon
// links to the nearest selected truth photon closer than Delta R = 0.1, and is
// a fake recording that search radius otherwise; each selected truth photon
// that no reconstructed photon reached is a miss. Distance alone decides, and
// several reconstructed photons may link one truth photon. The identity
// association (matchPhotonRecords) is the analysis relation.
std::vector<photonjet::PhotonTruthLinkRecord> linkPhotonsByDistance(
    const photonjet::EventRecord& event, const std::vector<photonjet::PhotonRecord>& reco,
    const std::vector<photonjet::TruthPhotonRecord>& truth, const std::vector<bool>& selectedTruth)
{
  using namespace photonjet;
  constexpr double kSearchRadius = 0.1;
  std::vector<PhotonTruthLinkRecord> result;
  std::vector<bool> reached(truth.size(), false);
  for (const auto& photon : reco)
  {
    PhotonTruthLinkRecord row;
    row.eventId = event.eventId;
    row.photonId = photon.photonId;
    row.deltaR = kSearchRadius;
    std::size_t nearest = truth.size();
    for (std::size_t t = 0; t < truth.size(); ++t)
    {
      if (!selectedTruth[t]) continue;
      const double distance = std::hypot(photon.eta - truth[t].eta, TVector2::Phi_mpi_pi(photon.phi - truth[t].phi));
      if (distance < row.deltaR)
      {
        row.deltaR = distance;
        nearest = t;
      }
    }
    row.state = nearest < truth.size() ? AssociationState::Matched : AssociationState::Fake;
    if (nearest < truth.size())
    {
      row.truthPhotonId = truth[nearest].truthPhotonId;
      reached[nearest] = true;
    }
    row.linkId = makeLinkIdentity(event.eventId, row.photonId, row.truthPhotonId, LinkKind::PhotonToTruthPhoton);
    result.push_back(std::move(row));
  }
  for (std::size_t t = 0; t < truth.size(); ++t)
  {
    if (!selectedTruth[t] || reached[t]) continue;
    PhotonTruthLinkRecord row;
    row.eventId = event.eventId;
    row.truthPhotonId = truth[t].truthPhotonId;
    row.state = AssociationState::Miss;
    row.linkId = makeLinkIdentity(event.eventId, row.photonId, row.truthPhotonId, LinkKind::PhotonToTruthPhoton);
    result.push_back(std::move(row));
  }
  return result;
}

// One identity join supplies both the raw association and the producer's
// selected-signal links. Several clusters can match one truth photon. A raw
// match outside the signal subset remains a raw match and a selected fake.
std::vector<photonjet::PhotonTruthLinkRecord> matchPhotonRecords(const photonjet::EventRecord& event,
                                                                 const std::vector<photonjet::PhotonRecord>& reco,
                                                                 const std::vector<photonjet::TruthPhotonRecord>& truth,
                                                                 const std::vector<bool>& selectedTruth)
{
  using namespace photonjet;
  if (!event.eventId.valid() || selectedTruth.size() != truth.size())
    throw std::invalid_argument("photon association needs an event and explicit truth subset");
  std::map<std::pair<int, int>, std::size_t> byNativeKey;
  std::set<std::pair<std::uint64_t, std::uint64_t>> truthIds, recoIds;
  for (std::size_t t = 0; t < truth.size(); ++t)
  {
    const auto& row = truth[t];
    if (row.eventId != event.eventId || !row.truthPhotonId.valid() || !row.geantValid || row.trackId <= 0 ||
        !truthIds.emplace(row.truthPhotonId.hi, row.truthPhotonId.lo).second ||
        !byNativeKey.emplace(std::make_pair(row.embeddingId, row.trackId), t).second)
      throw std::invalid_argument("invalid or ambiguous event-scoped truth photon inventory");
  }
  std::vector<PhotonTruthLinkRecord> result;
  std::set<std::size_t> matched;
  bool unresolved = false;
  const auto append = [&](PhotonTruthLinkRecord row)
  {
    row.linkId = makeLinkIdentity(event.eventId, row.photonId, row.truthPhotonId, LinkKind::PhotonToTruthPhoton);
    result.push_back(std::move(row));
  };
  for (const auto& photon : reco)
  {
    if (photon.eventId != event.eventId || !photon.photonId.valid() ||
        !recoIds.emplace(photon.photonId.hi, photon.photonId.lo).second)
      throw std::invalid_argument("invalid or duplicate event-scoped reco photon identity");
    PhotonTruthLinkRecord row;
    row.eventId = event.eventId;
    row.photonId = photon.photonId;
    const auto& witness = photon.dominantTruth;
    const bool known = witness.evaluator == DominantTruthEvaluator::TowerInfo ||
                       witness.evaluator == DominantTruthEvaluator::LegacyRawTower;
    if (known && witness.state == DominantTruthState::NoPrimary)
    {
      row.associationState = 1;
      row.state = AssociationState::Fake;
    }
    else if (known && witness.state == DominantTruthState::ValidPrimary && witness.trackId > 0 && witness.pid != 0 &&
             std::isfinite(witness.energyContribution) && witness.energyContribution >= 0.)
    {
      row.energyContribution = witness.energyContribution;
      if (witness.pid != 22)
      {
        row.associationState = 3;
        row.state = AssociationState::Fake;
      }
      else
      {
        const auto found = byNativeKey.find({witness.embeddingId, witness.trackId});
        if (found != byNativeKey.end())
        {
          const auto t = found->second;
          row.associationState = 2;
          row.associationTruthPhotonId = truth[t].truthPhotonId;
          if (std::isfinite(photon.eta) && std::isfinite(photon.phi) && std::isfinite(truth[t].eta) &&
              std::isfinite(truth[t].phi))
            row.associationDeltaR =
                associationSeparation(photon.eta - truth[t].eta, std::remainder(photon.phi - truth[t].phi, 2. * kPi));
          if (selectedTruth[t])
          {
            row.truthPhotonId = truth[t].truthPhotonId;
            row.state = AssociationState::Matched;
            row.deltaR = row.associationDeltaR;
            matched.insert(t);
          }
          else
            row.state = AssociationState::Fake;
        }
      }
    }
    if (row.associationState == 0) unresolved = true;
    append(std::move(row));
  }
  const bool complete =
      event.recoPhotonCaptureState == CaptureState::Complete && event.recoPhotonUnclassifiableCount == 0;
  for (std::size_t t = 0; t < truth.size(); ++t)
  {
    if (!selectedTruth[t] || matched.count(t)) continue;
    PhotonTruthLinkRecord row;
    row.eventId = event.eventId;
    row.truthPhotonId = truth[t].truthPhotonId;
    row.state = complete && !unresolved ? AssociationState::Miss : AssociationState::Unknown;
    append(std::move(row));
  }
  return result;
}

// Pure matching on captured records. The persisted contract stores every
// candidate edge first, then selected edges, fakes and misses. Truth jets are
// shared across reconstructed views; matching them once per view changes that
// contract. Identity tie breaks follow canonical collection/encounter order.
std::vector<photonjet::JetTruthLinkRecord> matchJetRecords(const photonjet::EventRecord& event,
                                                           const std::vector<photonjet::JetRecord>& reco,
                                                           const std::vector<photonjet::TruthJetRecord>& truth,
                                                           const double maximumDeltaR)
{
  using namespace photonjet;
  if (!std::isfinite(maximumDeltaR) || maximumDeltaR <= 0.)
    throw std::invalid_argument("jet matching radius must be finite and positive");
  struct Edge
  {
    std::size_t reco, truth;
    double distance;
  };
  std::vector<Edge> edges;
  std::vector<double> nearestReco(reco.size(), kNaN), nearestTruth(truth.size(), kNaN);
  const auto recoKey = [&](std::size_t i)
  { return std::make_tuple(reco[i].view == JetView::AuAuNoSub, reco[i].radius, reco[i].deterministicOrder); };
  const auto truthKey = [&](std::size_t i) { return std::make_pair(truth[i].radius, truth[i].nativeKey); };
  for (std::size_t r = 0; r < reco.size(); ++r)
    for (std::size_t t = 0; t < truth.size(); ++t)
    {
      if (!reco[r].jetId.valid() || !truth[t].truthJetId.valid() ||
          std::abs(reco[r].radius - truth[t].radius) > kRadiusTolerance)
        continue;
      const double phi = reco[r].phi - truth[t].phi;
      const double distance = std::hypot(reco[r].eta - truth[t].eta, std::atan2(std::sin(phi), std::cos(phi)));
      if (!std::isfinite(distance)) continue;
      if (!std::isfinite(nearestReco[r]) || distance < nearestReco[r]) nearestReco[r] = distance;
      if (!std::isfinite(nearestTruth[t]) || distance < nearestTruth[t]) nearestTruth[t] = distance;
      if (distance < maximumDeltaR) edges.push_back({r, t, distance});
    }
  // Equal distances: the larger raw (pre-JES) jet pT first. The established
  // output orders every tie this way (with the calibrated pT, sub-GeV Au+Au
  // jets, whose JES factors differ between views, would swap).
  std::sort(edges.begin(), edges.end(),
            [&](const Edge& a, const Edge& b)
            {
              if (a.distance != b.distance) return a.distance < b.distance;
              if (reco[a.reco].rawPt != reco[b.reco].rawPt) return reco[a.reco].rawPt > reco[b.reco].rawPt;
              if (recoKey(a.reco) != recoKey(b.reco)) return recoKey(a.reco) < recoKey(b.reco);
              return truthKey(a.truth) < truthKey(b.truth);
            });
  std::vector<JetTruthLinkRecord> result;
  const auto append = [&](std::size_t r, std::size_t t, double distance, AssociationState state)
  {
    JetTruthLinkRecord row;
    row.eventId = event.eventId;
    if (r < reco.size())
    {
      row.jetId = reco[r].jetId;
      row.jetView = reco[r].view;
      row.jetRadius = reco[r].radius;
    }
    else
    {
      row.jetRadius = truth[t].radius;
      row.jetView = event.recoJetView.empty() ? JetView::PP : static_cast<JetView>(event.recoJetView.front());
    }
    if (t < truth.size()) row.truthJetId = truth[t].truthJetId;
    row.deltaR = distance;
    row.state = state;
    row.selectedMatch = state == AssociationState::Matched;
    row.linkId = makeLinkIdentity(event.eventId, row.jetId, row.truthJetId, LinkKind::JetToTruthJet);
    row.linkId.lo = mixIdentityWord(row.linkId.lo ^ static_cast<std::uint64_t>(state));
    result.push_back(row);
  };
  for (const auto& edge : edges) append(edge.reco, edge.truth, edge.distance, AssociationState::Unmatched);
  std::set<std::size_t> assignedReco, assignedTruth;
  for (const auto& edge : edges)
    if (!assignedReco.count(edge.reco) && !assignedTruth.count(edge.truth))
    {
      assignedReco.insert(edge.reco);
      assignedTruth.insert(edge.truth);
      append(edge.reco, edge.truth, edge.distance, AssociationState::Matched);
    }
  std::vector<std::size_t> recoOrder(reco.size()), truthOrder(truth.size());
  std::iota(recoOrder.begin(), recoOrder.end(), 0);
  std::iota(truthOrder.begin(), truthOrder.end(), 0);
  std::sort(recoOrder.begin(), recoOrder.end(), [&](auto a, auto b) { return recoKey(a) < recoKey(b); });
  std::sort(truthOrder.begin(), truthOrder.end(), [&](auto a, auto b) { return truthKey(a) < truthKey(b); });
  // As in the established producer, an unmatched reconstructed jet is a fake
  // and an unmatched truth jet a miss, including when the other side has no
  // jets in the event (e.g. no reconstructed vertex).
  for (auto r : recoOrder)
    if (!assignedReco.count(r)) append(r, truth.size(), nearestReco[r], AssociationState::Fake);
  for (auto t : truthOrder)
    if (!assignedTruth.count(t)) append(reco.size(), t, nearestTruth[t], AssociationState::Miss);
  return result;
}


/**
 * Native truth-photon identity inside one event:
 *
 *   embedding id + Geant4 track id
 */
using TruthPhotonKey = std::pair<int, int>;

// Pseudorapidity as the established stitch computed it from the momentum.
double safeEta(const double px, const double py, const double pz)
{
  const double p2 = px * px + py * py + pz * pz;
  if (!std::isfinite(p2) || p2 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double p = std::sqrt(p2);
  if (p - pz <= 0.0 || p + pz <= 0.0) return (pz >= 0.0 ? 1.0 : -1.0) * std::numeric_limits<double>::infinity();
  return 0.5 * std::log((p + pz) / (p - pz));
}

// Maximum pT of the stored Geant4 primary photons (optionally embedded) with
// pT > 1 GeV, E > 1 GeV and |eta| <= etaAbsMax; -1 without one.
double maxStoredTruthPhotonPt(PHG4TruthInfoContainer* truth, const bool requireEmbedded, const double etaAbsMax)
{
  if (!truth) return -1.0;
  double maxPt = -1.0;
  const auto primaries = truth->GetPrimaryParticleRange();
  for (auto it = primaries.first; it != primaries.second; ++it)
  {
    const PHG4Particle* p = it->second;
    if (!p || (requireEmbedded && truth->isEmbeded(p->get_track_id()) < 1) || p->get_pid() != 22) continue;
    const double pt = std::hypot(p->get_px(), p->get_py());
    const double eta = safeEta(p->get_px(), p->get_py(), p->get_pz());
    if (!std::isfinite(pt) || pt <= 1.0 || !std::isfinite(p->get_e()) || p->get_e() <= 1.0 || !std::isfinite(eta) ||
        std::fabs(eta) > etaAbsMax)
      continue;
    if (pt > maxPt) maxPt = pt;
  }
  return maxPt;
}

// Maximum pT of the generator photons with |eta| <= 1.5, optionally final-state
// and optionally produced by a quark, gluon or photon (|pdg| 1-22); -1 without.
double maxGeneratorPhotonPt(PHHepMCGenEventMap* events, const bool requireParent, const bool requireStable)
{
  PHHepMCGenEvent* generated = nullptr;
  if (events)
  {
    generated = events->get(0);
    if (!generated) generated = events->get(1);
    if (!generated && !events->empty()) generated = events->begin()->second;
  }
  const HepMC::GenEvent* event = generated ? generated->getEvent() : nullptr;
  if (!event) return -1.0;
  double maxPt = -1.0;
  for (auto it = event->particles_begin(); it != event->particles_end(); ++it)
  {
    const HepMC::GenParticle* p = *it;
    if (!p || p->pdg_id() != 22 || (requireStable && p->status() != 1)) continue;
    const double pt = std::hypot(p->momentum().px(), p->momentum().py());
    const double eta = p->momentum().pseudoRapidity();
    if (!std::isfinite(pt) || !std::isfinite(eta) || pt <= 0.0 || eta < -1.5 || eta > 1.5) continue;
    if (requireParent)
    {
      const HepMC::GenVertex* vertex = p->production_vertex();
      bool parent = false;
      if (vertex)
        for (auto in = vertex->particles_in_const_begin(); in != vertex->particles_in_const_end() && !parent; ++in)
          parent = *in && std::abs((*in)->pdg_id()) >= 1 && std::abs((*in)->pdg_id()) <= 22;
      if (!parent) continue;
    }
    if (pt > maxPt) maxPt = pt;
  }
  return maxPt;
}

}  // namespace


// ============================================================================
// Simulation truth stage
//
// Build the native truth inventories before constructing any relation to
// reconstructed objects. This keeps the truth denominator independent of the
// success or failure of reconstruction matching.
// ============================================================================

void PhotonJetTree::captureSimulationTruth(PHCompositeNode* topNode)
{
  if (!isSimulation())
  {
    return;
  }

  if (m_config.system == photonjet::CollisionSystem::AuAu &&
      m_config.simulationRole == photonjet::SimulationRole::PhotonJet)
  {
    m_event.embeddedFilterPhotonPt = embeddedFilterPhotonPt(topNode, m_config.nodes.hepmcEventMap);
    m_event.embeddedFilterPhotonCaptured = true;
  }

  // The complete truth denominator is captured before any matching.
  captureTruthVertices(topNode);
  captureTruthPhotons(topNode);
  captureTruthJets(topNode);
  if (m_config.weights.ppg12Nominal) capturePPG12OwnershipInputs(topNode);
}


// ============================================================================
// PPG12 sample-ownership inputs (pp simulation)
//
// The maximum stored truth photon pT of each established stitch definition and
// the leading truth R=0.4 jet pT, as the producer computed them. The ownership
// decision and its audit rows belong to the event classification (Event.cc).
// ============================================================================

/**
 * The truth inputs of the ownership decision.
 */
void PhotonJetTree::capturePPG12OwnershipInputs(PHCompositeNode* topNode)
{
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);
  auto* events = findNode::getClass<PHHepMCGenEventMap>(topNode, m_config.nodes.hepmcEventMap);
  m_event.ppg12StitchPhotonPt = {{maxStoredTruthPhotonPt(truth, true, 1.5), maxGeneratorPhotonPt(events, true, false),
                                  maxGeneratorPhotonPt(events, true, true), maxGeneratorPhotonPt(events, false, false),
                                  maxStoredTruthPhotonPt(truth, false, 1.5)}};
  m_event.ppg12SourcePhotonPt = maxStoredTruthPhotonPt(truth, true, 0.7);
  for (const auto& collection : m_config.jets.nodes)
  {
    if (collection.truthNode.empty() || std::fabs(collection.radius - 0.4) > 1e-6) continue;
    if (auto* jets = findNode::getClass<JetContainer>(topNode, collection.truthNode))
      for (const Jet* jet : *jets)
        if (jet && std::isfinite(jet->get_pt()) && jet->get_pt() > m_event.ppg12MaxTruthJetPtR04)
          m_event.ppg12MaxTruthJetPtR04 = jet->get_pt();
    break;
  }
}


// ============================================================================
// Truth vertices
//
// Store the primary Geant4 vertices together with explicit event-level
// witnesses for the hard-scatter and minimum-bias vertices.
//
// The PHG4TruthInfoContainer primary vertex defines the hard-scatter witness.
// In a double-interaction event, the first other finite primary vertex is
// retained as the minimum-bias companion.
// ============================================================================

void PhotonJetTree::captureTruthVertices(PHCompositeNode* topNode)
{
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);

  if (!truth)
  {
    m_event.truthVertexCaptureState = CaptureState::InputUnavailable;

    return;
  }

  m_event.truthVertexCaptureState = CaptureState::Complete;

  // --------------------------------------------------------------------------
  // Map native Geant4 vertex ids onto embedding ids.
  // --------------------------------------------------------------------------

  std::map<int, int> embeddingByVertex;

  {
    const auto embedded = truth->GetEmbeddedVtxIds();

    for (auto it = embedded.first; it != embedded.second; ++it)
    {
      embeddingByVertex.emplace(it->first, it->second);
    }
  }

  // --------------------------------------------------------------------------
  // Hard-scatter vertex witness.
  // --------------------------------------------------------------------------

  const int primaryVertexId = truth->GetPrimaryVertexIndex();

  m_event.truthPrimaryVertexId = primaryVertexId;

  if (const PHG4VtxPoint* primary = truth->GetVtx(primaryVertexId); primary && finiteValue(primary->get_z()))
  {
    m_event.truthHardVertexZ = primary->get_z();

    m_event.truthHardVertexValid = true;
  }

  // --------------------------------------------------------------------------
  // Primary truth vertices, the established producer's vertex census: the
  // hard scatter and, in double-interaction samples, its companion. Secondary
  // Geant4 vertices are not recorded.
  // --------------------------------------------------------------------------

  const auto vertices = truth->GetPrimaryVtxRange();

  for (auto it = vertices.first; it != vertices.second; ++it)
  {
    TruthVertexRecord row;

    row.eventId = m_event.eventId;
    row.vertexId = it->first;

    const auto embedding = embeddingByVertex.find(it->first);

    row.embeddingValid = embedding != embeddingByVertex.end();

    row.embeddingId = row.embeddingValid ? embedding->second : 0;

    row.z = it->second ? it->second->get_z() : kNaN;

    row.valid = it->second && finiteValue(row.z);

    m_truthVertices.push_back(row);
  }

  // --------------------------------------------------------------------------
  // Minimum-bias companion vertex witness.
  // --------------------------------------------------------------------------

  const auto primaries = truth->GetPrimaryVtxRange();

  for (auto it = primaries.first; it != primaries.second; ++it)
  {
    if (it->first == primaryVertexId || !it->second || !finiteValue(it->second->get_z()))
    {
      continue;
    }

    m_event.truthMinimumBiasVertexZ = it->second->get_z();

    m_event.truthMinimumBiasVertexValid = true;

    break;
  }
}


// ============================================================================
// Truth photons
//
// The truth-photon inventory begins from embedded Geant4 primary photons.
// Geant4 identity, generator association, prompt class, and truth isolation are
// recorded independently so a missing generator history does not erase an
// otherwise valid simulated photon.
//
// Any failure that compromises denominator completeness is counted explicitly.
// ============================================================================

void PhotonJetTree::captureTruthPhotons(PHCompositeNode* topNode)
{
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);

  if (!truth)
  {
    m_event.truthPhotonCaptureState = CaptureState::InputUnavailable;

    m_event.truthDenominatorComplete = false;

    return;
  }

  m_event.truthPhotonCaptureState = CaptureState::Complete;

  auto* generatorEvents = findNode::getClass<PHHepMCGenEventMap>(topNode, m_config.nodes.hepmcEventMap);

  const double signalRadius = m_config.photon.truthIsolationSignalRadius;

  std::set<TruthPhotonKey> seen;

  const auto particles = truth->GetPrimaryParticleRange();

  for (auto it = particles.first; it != particles.second; ++it)
  {
    PHG4Particle* particle = it->second;

    if (!particle)
    {
      ++m_event.truthPhotonIsolationIncompleteCount;
      continue;
    }

    const int trackId = particle->get_track_id();

    const int embeddingId = truth->isEmbeded(trackId);

    /*
     * Truth-photon population:
     *
     *   - embedded primary subevent
     *   - photon PID
     *
     * No reconstructed-object requirement enters this census.
     */
    if (embeddingId < 1 || particle->get_pid() != 22)
    {
      continue;
    }

    ++m_event.truthPhotonEmbeddedPrimaryCount;

    /*
     * (embedding id, Geant4 track id) must uniquely identify one truth photon
     * inside the event. Duplicate native identities invalidate denominator
     * completeness rather than being silently merged.
     */
    if (!seen.insert({embeddingId, trackId}).second)
    {
      ++m_event.truthPhotonDuplicateTrackCount;
      continue;
    }

    double pt, eta, phi;

    if (!kinematics(particle, pt, eta, phi))
    {
      ++m_event.truthPhotonRejectedCount;
      ++m_event.truthPhotonRejectedKinematicsCount;

      continue;
    }

    // ------------------------------------------------------------------------
    // Truth isolation
    // ------------------------------------------------------------------------

    double isolationR03 = kNaN;
    double isolationR04 = kNaN;

    bool isolationComplete = false;

    const bool isolationValid =
        calculateTruthPhotonIsolation(topNode, particle, isolationR03, isolationR04, isolationComplete);

    if (!isolationComplete)
    {
      ++m_event.truthPhotonIsolationIncompleteCount;
      ++m_event.truthPhotonRejectedIsolationCount;
    }

    // ------------------------------------------------------------------------
    // Generator association and prompt classification
    // ------------------------------------------------------------------------

    const int barcode = particle->get_barcode();

    const HepMC::GenParticle* generatorPhoton =
        generatorParticle(generatorEvent(generatorEvents, embeddingId), barcode);

    const TruthPhotonClass photonClass = promptClass(generatorPhoton);

    /*
     * In archived double-interaction production, the hard subevent alone owns
     * the signal role. Source ownership is stored separately from the native
     * photon truth-signal predicate below.
     */
    const bool ownership = !isArchivedDoubleInteraction() || embeddingId == kDoubleInteractionHardEmbedding;

    const double signalIsolation =
        std::abs(signalRadius - kTruthIsolationR03) < kRadiusTolerance ? isolationR03 : isolationR04;

    /*
     * Native truth-photon signal predicate.
     *
     * This records the established producer reference definition. Source-role
     * ownership remains separate in signalSourceRole.
     */
    const bool analysisSignal =
        std::abs(eta) < m_config.photon.maxAbsEta && photonjet::isAnalysisSignalClass(photonClass) && isolationValid &&
        finiteValue(signalIsolation) && signalIsolation < m_config.photon.truthIsolationMaxEtGeV;

    // ------------------------------------------------------------------------
    // Persisted truth-photon record
    // ------------------------------------------------------------------------

    TruthPhotonRecord row;

    row.eventId = m_event.eventId;

    row.truthPhotonId = photonjet::makeTruthPhotonIdentity(m_event.eventId, embeddingId, trackId);

    row.trackId = trackId;

    row.vertexId = particle->get_vtx_id();

    row.barcode = barcode;

    row.embeddingId = embeddingId;

    row.pid = 22;

    row.pt = pt;

    row.eta = eta;

    row.phi = phi;

    row.energy = particle->get_e();

    row.geantValid = true;

    row.generatorAssociationValid = generatorPhoton && generatorPhoton->pdg_id() == 22;

    row.promptClass = photonClass;

    row.isolationR03 = isolationR03;

    row.isolationR04 = isolationR04;

    row.isolationValid = isolationValid;

    row.analysisSignal = analysisSignal;

    row.sampleRole = m_config.simulationRole;

    row.signalSourceRole = ownership;  // separate from the native photon signal predicate

    m_truthPhotons.push_back(std::move(row));

    ++m_event.truthPhotonWrittenCount;

    if (analysisSignal)
    {
      ++m_event.truthPhotonAnalysisSignalCount;
    }
  }

  /*
   * A denominator is complete only when every candidate identity was unique,
   * truth-isolation inputs were complete, and no primary photon was lost
   * because its kinematics could not be represented.
   */
  m_event.truthDenominatorComplete = m_event.truthPhotonDuplicateTrackCount == 0 &&
                                     m_event.truthPhotonIsolationIncompleteCount == 0 &&
                                     m_event.truthPhotonRejectedKinematicsCount == 0;
}


// ============================================================================
// Truth-photon isolation
//
// Reproduce the predecessor isolation convention using transverse energy from
// TLorentzVector::Et(). Every embedded primary particle contributes according
// to its geometric distance from the photon axis.
//
// The central merged core is accumulated separately and removed from both cone
// sums. Input completeness and finite resulting sums are tracked independently.
// ============================================================================

bool PhotonJetTree::calculateTruthPhotonIsolation(PHCompositeNode* topNode, const PHG4Particle* photon,
                                                  double& isolationR03, double& isolationR04, bool& complete) const
{
  isolationR03 = isolationR04 = kNaN;
  complete = false;

  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);

  if (!truth || !photon)
  {
    return false;
  }

  double photonPt, photonEta, photonPhi;

  if (!kinematics(photon, photonPt, photonEta, photonPhi))
  {
    return false;
  }

  /*
   * RecoilObjects.cc:17735-17779 uses TLorentzVector::Et, not pT. In
   * particular massive primary particles have a different transverse energy.
   */
  const TLorentzVector axis(photon->get_px(), photon->get_py(), photon->get_pz(), photon->get_e());

  const double core = m_config.photon.truthIsolationCoreRadius;

  double sumR03 = 0.0;
  double sumR04 = 0.0;
  double coreSum = 0.0;

  complete = true;

  const auto particles = truth->GetPrimaryParticleRange();

  for (auto it = particles.first; it != particles.second; ++it)
  {
    const PHG4Particle* particle = it->second;

    if (!particle)
    {
      complete = false;
      continue;
    }

    if (truth->isEmbeded(particle->get_track_id()) < 1)
    {
      continue;
    }

    const TLorentzVector momentum(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e());

    const double et = momentum.Et();

    if (!finiteValue(et) || et <= 0.0)
    {
      if (!finiteValue(et) || et < 0.0)
      {
        complete = false;
      }

      continue;
    }

    const double dr = axis.DeltaR(momentum);

    if (!finiteValue(dr))
    {
      complete = false;
      continue;
    }

    if (dr < kTruthIsolationR04) sumR04 += et;
    if (dr < kTruthIsolationR03) sumR03 += et;
    if (dr < core) coreSum += et;
  }

  isolationR03 = sumR03 - coreSum;

  isolationR04 = sumR04 - coreSum;

  return finiteValue(isolationR03) && finiteValue(isolationR04);
}


// ============================================================================
// Truth jets
//
// Capture every unique configured truth-jet collection once. Several
// reconstructed views may point to the same truth node, so (radius, node name)
// is used to suppress duplicate reads.
//
// Truth jets are stored exactly as provided by the truth collection. No JES or
// reconstructed-jet correction is applied.
// ============================================================================

void PhotonJetTree::captureTruthJets(PHCompositeNode* topNode)
{
  m_event.truthJetCaptureState = CaptureState::Complete;

  std::set<std::pair<int, std::string>> visited;

  for (const JetNodeConfig& collection : m_config.jets.nodes)
  {
    if (collection.truthNode.empty())
    {
      continue;
    }

    const int radiusCode = static_cast<int>(std::lround(collection.radius * 100.0));

    if (!visited.insert({radiusCode, collection.truthNode}).second)
    {
      continue;  // several reconstructed views share one truth collection
    }

    JetContainer* container = findNode::getClass<JetContainer>(topNode, collection.truthNode);

    m_event.truthJetRadiusCode.push_back(radiusCode);

    m_event.truthJetContainerValid.push_back(container ? 1 : 0);

    if (!container)
    {
      m_event.truthJetCaptureState = CaptureState::InputUnavailable;
    }

    if (!container)
    {
      continue;
    }
    m_ppg12TruthJetRadiiFound.push_back(radiusCode);

    for (Jet* jet : *container)
    {
      if (!jet)
      {
        m_event.truthJetContainerValid.back() = 0;
        continue;
      }

      const double pt = jet->get_pt();

      const double eta = jet->get_eta();

      const double phi = jet->get_phi();

      if (!finiteValue(pt) || !finiteValue(eta) || !finiteValue(phi) || pt < 0.0)
      {
        m_event.truthJetContainerValid.back() = 0;
        continue;
      }

      TruthJetRecord row;

      row.eventId = m_event.eventId;

      row.radius = collection.radius;

      row.nativeKey = jet->get_id();

      row.truthJetId = photonjet::makeTruthJetIdentity(m_event.eventId, collection.radius, row.nativeKey);

      row.pt = pt;

      row.eta = eta;

      row.phi = phi;

      row.containerValid = true;

      m_truthJets.push_back(std::move(row));
    }
  }
}


// ============================================================================
// Reconstructed-photon -> truth-photon relations
//
// For every reconstructed photon, CaloRawClusterEval identifies the Geant4
// primary contributing the most energy to the calorimeter cluster.
//
// Association is then performed by native truth identity:
//
//   dominant primary
//       -> embedding id + Geant4 track id
//       -> retained TruthPhotonRecord
//
// Delta-R is computed only after an identity match and is retained as a
// diagnostic. It is never used as the photon-matching criterion.
//
// Relation semantics
// ------------------
//   Matched
//     dominant primary is a retained truth photon in the native signal subset
//
//   Fake
//     evaluator proved no primary, the primary is not a photon, or the raw
//     association is outside the native signal subset
//
//   Unknown
//     evaluator/input is incomplete, the primary is
//     malformed, or a photon primary exists without a retained truth row
//
//   Miss
//     selected truth photon matched by no reconstructed photon, asserted only
//     when reconstructed-photon capture and all candidate associations are
//     complete
// ============================================================================

void PhotonJetTree::buildPhotonTruthAssociations(PHCompositeNode* topNode)
{
  m_photonTruthLinks.clear();
  m_photonDistanceLinks.clear();
  if (!isSimulation()) return;
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);
  auto* clusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.photonNode);

  // The established evaluator: tower-info truth nodes when the event has
  // them, otherwise raw towers, chosen once per event and loaded for it.
  CaloRawClusterEval evaluator(topNode, "CEMC");
  evaluator.set_usetowerinfo(true);
  evaluator.next_event(topNode);
  DominantTruthEvaluator mode = DominantTruthEvaluator::TowerInfo;
  if (!evaluator.has_reduced_node_pointers())
  {
    evaluator.set_usetowerinfo(false);
    evaluator.next_event(topNode);
    mode = DominantTruthEvaluator::LegacyRawTower;
  }
  m_ppg12EvaluatorReady = evaluator.has_reduced_node_pointers();

  // The maximum-energy primary of one photon cluster, captured once; the
  // record helper below owns every association and selection decision. A
  // photon-only search cannot reproduce this detector witness.
  const auto dominantPrimary = [&](const std::uint64_t key)
  {
    DominantTruthRecord witness;
    witness.state = DominantTruthState::EvaluatorUnavailable;
    witness.evaluator = mode;
    RawCluster* cluster = clusters ? clusters->getCluster(static_cast<unsigned int>(key)) : nullptr;
    if (!truth || !cluster || !evaluator.has_reduced_node_pointers()) return witness;
    PHG4Particle* primary = evaluator.max_truth_primary_particle_by_energy(cluster);
    witness.state = DominantTruthState::NoPrimary;
    if (!primary) return witness;
    witness.trackId = primary->get_track_id();
    witness.pid = primary->get_pid();
    witness.barcode = primary->get_barcode();
    witness.vertexId = primary->get_vtx_id();
    witness.energyContribution = evaluator.get_energy_contribution(cluster, primary);
    witness.state = DominantTruthState::InvalidPrimary;
    if (witness.trackId <= 0 || witness.pid == 0 || !std::isfinite(witness.energyContribution) ||
        witness.energyContribution < 0.)
      return witness;
    witness.embeddingId = truth->isEmbeded(witness.trackId);
    witness.state = DominantTruthState::ValidPrimary;
    return witness;
  };
  for (auto& photon : m_photons) photon.dominantTruth = dominantPrimary(photon.nativeClusterKey);
  for (auto& photon : m_ppg12Photons) photon.dominantTruth = dominantPrimary(photon.nativeClusterKey);
  if (recordsPPG12Simulation()) capturePPG12Clusters(topNode, evaluator);

  std::vector<bool> selected;
  for (const auto& row : m_truthPhotons) selected.push_back(row.analysisSignal && row.signalSourceRole);
  m_photonTruthLinks = matchPhotonRecords(m_event, m_photons, m_truthPhotons, selected);
  m_photonDistanceLinks = linkPhotonsByDistance(m_event, m_photons, m_truthPhotons, selected);
}

// ============================================================================
// PPG12 simulation bookkeeping clusters
//
// The established producer read the builder's input clusters, and the unsplit
// cluster node when the DST had one, directly: each cluster's energy, its
// energy vector and pseudorapidity from the reconstructed vertex, and the
// Geant4 primary depositing most of its energy. Its truth tags only consult
// that primary when the event has a signal photon.
// ============================================================================

void PhotonJetTree::capturePPG12Clusters(PHCompositeNode* topNode, CaloRawClusterEval& evaluator)
{
  bool signal = false;
  for (const auto& row : m_truthPhotons) signal |= row.analysisSignal;
  const CLHEP::Hep3Vector vertex(0.0, 0.0, m_event.recoVertexZ);
  const auto capture = [&](RawClusterContainer* container, std::vector<PPG12ClusterRecord>& records)
  {
    for (const auto& [key, cluster] : container->getClustersMap())
    {
      if (!cluster) continue;
      PPG12ClusterRecord record;
      record.energy = cluster->get_energy();
      const CLHEP::Hep3Vector energy = RawClusterUtility::GetEVec(*cluster, vertex);
      record.directionEnergy = energy.mag();
      if (std::isfinite(record.directionEnergy) && record.directionEnergy > 0.0)
      {
        TLorentzVector momentum;
        momentum.SetPxPyPzE(energy.x(), energy.y(), energy.z(), record.directionEnergy);
        record.pt = momentum.Pt();
        record.eta = momentum.Eta();
        record.phi = TVector2::Phi_mpi_pi(momentum.Phi());
      }
      record.vertexEta = RawClusterUtility::GetPseudorapidity(*cluster, vertex);
      if (signal && evaluator.has_reduced_node_pointers())
        if (PHG4Particle* primary = evaluator.max_truth_primary_particle_by_energy(cluster))
        {
          record.primaryFound = true;
          record.truthTrackId = primary->get_track_id();
          record.energyContribution = evaluator.get_energy_contribution(cluster, primary);
        }
      records.push_back(record);
    }
  };
  if (auto* clusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.rawClusterNode))
    capture(clusters, m_ppg12Clusters);
  if (auto* unsplit = findNode::getClass<RawClusterContainer>(topNode, kPPG12ResponseClusterNode))
  {
    m_ppg12ResponseClustersPresent = true;
    capture(unsplit, m_ppg12ResponseClusters);
  }
}

// ============================================================================
// Au+Au embedded truth-isolation audit (embdiag AUDIT and TRUTH_PHOTON)
//
// Once per accepted embedded event the established producer audited the
// truth photons behind its fixed R = 0.3 truth-isolation diagnostic: |vz| below
// 30 cm, a centrality class (0-20, 20-50, 50-80%) and its truth-photon map;
// then, per photon in track order, its class, |eta| < 0.7 and isolation,
// logging each accepted photon. Its map held the embedded primary photons
// with finite kinematics and finite R = 0.3 and 0.4 truth isolation.
// ============================================================================

void PhotonJetTree::recordEmbeddedTruthIsolationAudit()
{
  const int sample = isPhotonJetSimulation() ? photonjet::production::EmbeddedPhotonSampleCode(m_config.source.sample)
                                             : photonjet::production::EmbeddedJetSampleCode(m_config.source.sample);
  if (sample == 0) return;
  const auto audit = [&](const int bin)
  {
    EmbeddedDiagnosticRecord row;
    row.kind = "AUDIT";
    row.auditBin = bin;
    row.sampleCode = sample;
    m_event.embeddedDiagnostics.push_back(row);
  };

  if (!std::isfinite(m_event.recoVertexZ) || std::fabs(m_event.recoVertexZ) >= 30.0)
  {
    audit(8);
    return;
  }
  const int percent = m_event.analysisCentrality;
  const int centrality = percent >= 0 && percent < 20    ? 0
                         : percent >= 20 && percent < 50 ? 1
                         : percent >= 50 && percent < 80 ? 2
                                                         : -1;
  if (centrality < 0)
  {
    audit(9);
    return;
  }
  if (m_event.truthPhotonCaptureState == CaptureState::InputUnavailable)
  {
    audit(10);
    return;
  }
  audit(1);
  audit(sample == 12 ? 12 : 13);
  audit(14 + centrality);
  for (const TruthPhotonRecord& photon : m_truthPhotons)
  {
    if (!std::isfinite(photon.isolationR03) || !std::isfinite(photon.isolationR04)) continue;
    audit(2);
    const int photonClass = static_cast<int>(photon.promptClass);
    if (photonClass != 1 && photonClass != 2)
    {
      audit(5);
      continue;
    }
    audit(photonClass == 1 ? 3 : 4);
    if (std::fabs(photon.eta) >= 0.7)
    {
      audit(7);
      continue;
    }
    EmbeddedDiagnosticRecord row;
    row.kind = "TRUTH_PHOTON";
    row.sampleCode = sample;
    row.centralityBinIndex = centrality;
    row.photonClass = photonClass;
    row.truthPt = photon.pt;
    row.truthIsolationWitness = photon.isolationR03;
    m_event.embeddedDiagnostics.push_back(row);
    if (photon.pt >= 10.0 && photon.pt < 12.0) audit(11);
  }
}


// ============================================================================
// Au+Au embedded truth-photon misses (truthmiss_*)
//
// The established producer's photon-count unfolding recorded a miss for an
// accepted embedded event with a truth signal photon when it selected no
// reconstructed photon, or selected one that is not that truth photon. The
// selected photon is the leading isolated tight photon: |eta| < 0.7,
// 15 <= ET < 35 GeV, the active preselection, a passing R = 0.4 isolation
// witness, and the production's cut-based tight identification on the
// photon's own shower shapes. It is the truth photon when its dominant truth
// primary is the leading truth signal photon.
// ============================================================================

void PhotonJetTree::recordEmbeddedTruthMisses()
{
  const photonjet::TruthMissConfig& c = m_config.truthMisses;
  if (!c.enabled) return;

  const TruthPhotonRecord* lead = nullptr;
  for (const TruthPhotonRecord& photon : m_truthPhotons)
    if (photon.analysisSignal && (!lead || photon.pt > lead->pt)) lead = &photon;
  if (!lead) return;

  const auto inside = [](const double value, const double low, const double high)
  { return low < value && value < high; };
  const auto ratio = [](const double a, const double b)
  { return std::isfinite(a) && std::isfinite(b) && b > 0 ? a / b : std::numeric_limits<double>::quiet_NaN(); };
  const PhotonRecord* selected = nullptr;
  for (const PhotonRecord& photon : m_photons)
  {
    if (!(std::fabs(photon.eta) < c.absEtaMax && photon.et >= c.etMinGeV && photon.et < c.etMaxGeV)) continue;
    if (photonjet::production::ActivePreselection(m_config.photon, photon) != 1) continue;
    const auto witness =
        std::find_if(photon.isolationWitnesses.begin(), photon.isolationWitnesses.end(),
                     [&](const photonjet::IsolationWitnessRecord& w) { return w.radius == c.isolationRadius; });
    if (witness == photon.isolationWitnesses.end() || witness->passState != 1) continue;
    const auto& shape = photon.primaryShapes;  // weta, wphi, ..., et1 [6], e11 [10], e33 [11], e32 [12], e35 [13]
    const double widthMax = c.widthMaxIntercept + c.widthMaxSlopePerGeV * photon.et;
    const bool tight = inside(shape[0], c.widthMin, widthMax) && inside(shape[1], c.widthMin, widthMax) &&
                       inside(ratio(shape[10], shape[11]), c.e11OverE33[0], c.e11OverE33[1]) &&
                       inside(shape[6], c.et1[0], c.et1[1]) &&
                       inside(ratio(shape[12], shape[13]), c.e32OverE35[0], c.e32OverE35[1]);
    if (tight && (!selected || photon.et > selected->et)) selected = &photon;
  }
  if (selected)
  {
    const DominantTruthRecord& truth = selected->dominantTruth;
    const bool primaryFound =
        truth.state == DominantTruthState::ValidPrimary || truth.state == DominantTruthState::InvalidPrimary;
    if (primaryFound && truth.trackId == lead->trackId) return;
  }

  const double centrality = m_event.analysisCentrality;
  photonjet::TruthMissRecord miss;
  miss.centralityBin = centrality >= 0 && centrality < 20    ? 0
                       : centrality >= 20 && centrality < 50 ? 1
                       : centrality >= 50 && centrality < 80 ? 2
                                                             : -1;
  miss.truthPt = lead->pt;
  miss.decisionPath = selected ? "SELECTED_RECO_NOT_TRUTH_SIGNAL" : "NO_SELECTED_RECO";
  m_event.truthMisses.push_back(miss);
}


// ============================================================================
// Reconstructed-jet -> truth-jet relations
//
// Same-radius edges share one truth assignment across the retained views.
// Candidate edges, selected edges and unmatched objects remain separate rows:
// an unmatched reconstructed jet is a fake, an unmatched truth jet a miss.
// ============================================================================

void PhotonJetTree::buildJetTruthAssociations()
{
  m_jetTruthLinks.clear();
  if (!isSimulation()) return;
  m_jetTruthLinks = matchJetRecords(m_event, m_jets, m_truthJets, m_config.jets.truthMatchMaxDeltaR);
}
