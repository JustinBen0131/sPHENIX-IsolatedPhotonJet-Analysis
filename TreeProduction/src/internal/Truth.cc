/**
 * @file Truth.cc
 * Capture simulation truth and reconstruction-to-truth associations.
 * Retain source identity, ancestry, truth isolation and matching validity.
 * Data has no invented truth; absent and ambiguous matches remain explicit.
 */

#include "../PhotonJetTree.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4eval/CaloRawClusterEval.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TLorentzVector.h>

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

bool finiteValue(const double value)
{
  return std::isfinite(value);
}


// ============================================================================
// Angular and particle-kinematic helpers
// ============================================================================

double wrapPhi(double phi)
{
  while (phi > kPi) phi -= 2.0 * kPi;
  while (phi <= -kPi) phi += 2.0 * kPi;

  return phi;
}


double deltaR(const double eta1,
              const double phi1,
              const double eta2,
              const double phi2)
{
  return std::hypot(
      eta1 - eta2,
      wrapPhi(phi1 - phi2));
}


/**
 * Read transverse momentum, pseudorapidity, and azimuth from a Geant4 particle.
 *
 * Invalid values are returned as NaN and the function returns false whenever
 * the momentum is non-finite or has non-positive transverse momentum.
 */
bool kinematics(const PHG4Particle* particle,
                double& pt,
                double& eta,
                double& phi)
{
  pt = eta = phi = kNaN;

  const double px = particle->get_px();
  const double py = particle->get_py();
  const double pz = particle->get_pz();

  if (!finiteValue(px) || !finiteValue(py) || !finiteValue(pz))
  {
    return false;
  }

  pt = std::hypot(px, py);

  if (!finiteValue(pt) || !(pt > 0.0))
  {
    return false;
  }

  eta = std::asinh(pz / pt);
  phi = std::atan2(py, px);

  return finiteValue(eta) && finiteValue(phi);
}


// ============================================================================
// Geant4 -> HepMC generator association
//
// Geant4 track identity and HepMC generator identity are distinct. Association
// proceeds through the embedding-id-specific generator event and then through
// the Geant4 particle barcode. Failure at either step remains an explicit
// missing generator association.
// ============================================================================

const HepMC::GenEvent* generatorEvent(const PHHepMCGenEventMap* events,
                                      const int embeddingId)
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
    throw std::runtime_error(
        "generator event map key disagrees with the event's embedding id");
  }

  return wrapper->getEvent();
}


const HepMC::GenParticle* generatorParticle(const HepMC::GenEvent* event,
                                            const int barcode)
{
  if (!event)
  {
    return nullptr;
  }

  for (auto it = event->particles_begin();
       it != event->particles_end();
       ++it)
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

  for (auto inItr = vertex->particles_in_const_begin();
       inItr != vertex->particles_in_const_end();
       ++inItr)
  {
    if (*inItr) incomingParticles.push_back(*inItr);
  }

  std::set<const HepMC::GenVertex*> seen{vertex};

  while (incomingParticles.size() == 1 &&
         incomingParticles[0] &&
         incomingParticles[0]->pdg_id() == 22)
  {
    vertex = incomingParticles[0]->production_vertex();

    if (!vertex || !seen.insert(vertex).second) return -1;

    incomingParticles.clear();

    for (auto inItr = vertex->particles_in_const_begin();
         inItr != vertex->particles_in_const_end();
         ++inItr)
    {
      if (*inItr) incomingParticles.push_back(*inItr);
    }
  }

  std::vector<const HepMC::GenParticle*> outgoingParticles;

  for (auto outItr = vertex->particles_out_const_begin();
       outItr != vertex->particles_out_const_end();
       ++outItr)
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

  if (incomingParticles.size() == 2 &&
      outgoingParticles.size() == 2)
  {
    const int in0 =
        incomingParticles[0]
            ? incomingParticles[0]->pdg_id()
            : 0;

    const int in1 =
        incomingParticles[1]
            ? incomingParticles[1]->pdg_id()
            : 0;

    const int out0 =
        outgoingParticles[0]
            ? outgoingParticles[0]->pdg_id()
            : 0;

    const int out1 =
        outgoingParticles[1]
            ? outgoingParticles[1]->pdg_id()
            : 0;

    if (std::abs(in0) <= 22 &&
        std::abs(in1) <= 22 &&
        std::abs(out0) <= 22 &&
        std::abs(out1) <= 22)
    {
      return 1;
    }
  }
  else if (incomingParticles.size() == 1 &&
           incomingParticles[0])
  {
    const int inPid =
        incomingParticles[0]->pdg_id();

    if (std::abs(inPid) <= 11 &&
        outgoingParticles.size() == 2)
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
  return static_cast<photonjet::TruthPhotonClass>(
      promptClassCode(photon));
}


// ============================================================================
// Reconstructed-to-truth jet matching helpers
//
// Candidate edges are built geometrically. Greedy one-to-one assignment then
// follows a deterministic ordering so the same event always produces the same
// selected matches independent of container iteration accidents.
// ============================================================================

// One identity join supplies both the raw association and the producer's
// selected-signal links. Several clusters can match one truth photon. A raw
// match outside the signal subset remains a raw match and a selected fake.
std::vector<photonjet::PhotonTruthLinkRecord> matchPhotonRecords(
    const photonjet::EventRecord& event,
    const std::vector<photonjet::PhotonRecord>& reco,
    const std::vector<photonjet::TruthPhotonRecord>& truth,
    const std::vector<bool>& selectedTruth)
{
  using namespace photonjet;
  if (!event.eventId.valid() || selectedTruth.size() != truth.size())
    throw std::invalid_argument("photon association needs an event and explicit truth subset");
  std::map<std::pair<int,int>, std::size_t> byNativeKey;
  std::set<std::pair<std::uint64_t,std::uint64_t>> truthIds, recoIds;
  for (std::size_t t = 0; t < truth.size(); ++t) {
    const auto& row = truth[t];
    if (row.eventId != event.eventId || !row.truthPhotonId.valid() ||
        !row.geantValid || row.trackId <= 0 || !truthIds.emplace(row.truthPhotonId.hi,row.truthPhotonId.lo).second ||
        !byNativeKey.emplace(std::make_pair(row.embeddingId,row.trackId),t).second)
      throw std::invalid_argument("invalid or ambiguous event-scoped truth photon inventory");
  }
  std::vector<PhotonTruthLinkRecord> result;
  std::set<std::size_t> matched;
  bool unresolved = false;
  const auto append = [&](PhotonTruthLinkRecord row) {
    row.linkId = makeLinkIdentity(event.eventId,row.photonId,row.truthPhotonId,LinkKind::PhotonToTruthPhoton);
    result.push_back(std::move(row));
  };
  for (const auto& photon : reco) {
    if (photon.eventId != event.eventId || !photon.photonId.valid() ||
        !recoIds.emplace(photon.photonId.hi,photon.photonId.lo).second)
      throw std::invalid_argument("invalid or duplicate event-scoped reco photon identity");
    PhotonTruthLinkRecord row; row.eventId = event.eventId; row.photonId = photon.photonId;
    const auto& witness = photon.dominantTruth;
    const bool known = witness.evaluator == DominantTruthEvaluator::TowerInfo ||
                       witness.evaluator == DominantTruthEvaluator::LegacyRawTower;
    if (known && witness.state == DominantTruthState::NoPrimary) {
      row.associationState = 1; row.state = AssociationState::Fake;
    } else if (known && witness.state == DominantTruthState::ValidPrimary &&
               witness.trackId > 0 && witness.pid != 0 &&
               std::isfinite(witness.energyContribution) && witness.energyContribution >= 0.) {
      row.energyContribution = witness.energyContribution;
      if (witness.pid != 22) {
        row.associationState = 3; row.state = AssociationState::Fake;
      } else {
        const auto found = byNativeKey.find({witness.embeddingId,witness.trackId});
        if (found != byNativeKey.end()) {
          const auto t = found->second;
          row.associationState = 2; row.associationTruthPhotonId = truth[t].truthPhotonId;
          if (std::isfinite(photon.eta) && std::isfinite(photon.phi) &&
              std::isfinite(truth[t].eta) && std::isfinite(truth[t].phi))
            row.associationDeltaR = std::hypot(photon.eta-truth[t].eta,
                std::remainder(photon.phi-truth[t].phi,2.*kPi));
          if (selectedTruth[t]) {
            row.truthPhotonId = truth[t].truthPhotonId; row.state = AssociationState::Matched;
            row.deltaR = row.associationDeltaR; matched.insert(t);
          } else row.state = AssociationState::Fake;
        }
      }
    }
    if (row.associationState == 0) unresolved = true;
    append(std::move(row));
  }
  const bool complete = event.recoPhotonCaptureState == CaptureState::Complete &&
                        event.recoPhotonUnclassifiableCount == 0;
  for (std::size_t t = 0; t < truth.size(); ++t) {
    if (!selectedTruth[t] || matched.count(t)) continue;
    PhotonTruthLinkRecord row; row.eventId = event.eventId; row.truthPhotonId = truth[t].truthPhotonId;
    row.state = complete && !unresolved ? AssociationState::Miss : AssociationState::Unknown;
    append(std::move(row));
  }
  return result;
}

// Pure matching on captured records. The persisted contract stores every
// candidate edge first, then selected edges, fakes and misses. Truth jets are
// shared across reconstructed views; matching them once per view changes that
// contract. Identity tie breaks follow canonical collection/encounter order.
std::vector<photonjet::JetTruthLinkRecord> matchJetRecords(
    const photonjet::EventRecord& event,
    const std::vector<photonjet::JetRecord>& reco,
    const std::vector<photonjet::TruthJetRecord>& truth,
    const double maximumDeltaR)
{
  using namespace photonjet;
  if (!std::isfinite(maximumDeltaR) || maximumDeltaR <= 0.)
    throw std::invalid_argument("jet matching radius must be finite and positive");
  struct Edge { std::size_t reco, truth; double distance; };
  std::vector<Edge> edges;
  std::vector<double> nearestReco(reco.size(), kNaN), nearestTruth(truth.size(), kNaN);
  const auto recoKey = [&](std::size_t i) {
    return std::make_tuple(reco[i].view == JetView::AuAuNoSub, reco[i].radius, reco[i].deterministicOrder);
  };
  const auto truthKey = [&](std::size_t i) { return std::make_pair(truth[i].radius, truth[i].nativeKey); };
  for (std::size_t r = 0; r < reco.size(); ++r)
    for (std::size_t t = 0; t < truth.size(); ++t)
    {
      if (!reco[r].jetId.valid() || !truth[t].truthJetId.valid() ||
          std::abs(reco[r].radius - truth[t].radius) > kRadiusTolerance) continue;
      const double phi = reco[r].phi - truth[t].phi;
      const double distance = std::hypot(reco[r].eta - truth[t].eta, std::atan2(std::sin(phi), std::cos(phi)));
      if (!std::isfinite(distance)) continue;
      if (!std::isfinite(nearestReco[r]) || distance < nearestReco[r]) nearestReco[r] = distance;
      if (!std::isfinite(nearestTruth[t]) || distance < nearestTruth[t]) nearestTruth[t] = distance;
      if (distance < maximumDeltaR) edges.push_back({r,t,distance});
    }
  std::sort(edges.begin(), edges.end(), [&](const Edge& a, const Edge& b) {
    if (a.distance != b.distance) return a.distance < b.distance;
    if (reco[a.reco].correctedPt != reco[b.reco].correctedPt)
      return reco[a.reco].correctedPt > reco[b.reco].correctedPt;
    if (recoKey(a.reco) != recoKey(b.reco)) return recoKey(a.reco) < recoKey(b.reco);
    return truthKey(a.truth) < truthKey(b.truth);
  });
  std::vector<JetTruthLinkRecord> result;
  const auto append = [&](std::size_t r, std::size_t t, double distance, AssociationState state) {
    JetTruthLinkRecord row;
    row.eventId = event.eventId;
    if (r < reco.size()) {
      row.jetId = reco[r].jetId; row.jetView = reco[r].view; row.jetRadius = reco[r].radius;
    } else {
      row.jetRadius = truth[t].radius;
      row.jetView = event.recoJetView.empty() ? JetView::PP : static_cast<JetView>(event.recoJetView.front());
    }
    if (t < truth.size()) row.truthJetId = truth[t].truthJetId;
    row.deltaR = distance; row.state = state; row.selectedMatch = state == AssociationState::Matched;
    row.linkId = makeLinkIdentity(event.eventId, row.jetId, row.truthJetId, LinkKind::JetToTruthJet);
    row.linkId.lo = mixIdentityWord(row.linkId.lo ^ static_cast<std::uint64_t>(state));
    result.push_back(row);
  };
  for (const auto& edge : edges) append(edge.reco, edge.truth, edge.distance, AssociationState::Unmatched);
  std::set<std::size_t> assignedReco, assignedTruth;
  for (const auto& edge : edges)
    if (!assignedReco.count(edge.reco) && !assignedTruth.count(edge.truth)) {
      assignedReco.insert(edge.reco); assignedTruth.insert(edge.truth);
      append(edge.reco, edge.truth, edge.distance, AssociationState::Matched);
    }
  std::vector<std::size_t> recoOrder(reco.size()), truthOrder(truth.size());
  std::iota(recoOrder.begin(), recoOrder.end(), 0);
  std::iota(truthOrder.begin(), truthOrder.end(), 0);
  std::sort(recoOrder.begin(), recoOrder.end(), [&](auto a, auto b) { return recoKey(a) < recoKey(b); });
  std::sort(truthOrder.begin(), truthOrder.end(), [&](auto a, auto b) { return truthKey(a) < truthKey(b); });
  for (auto r : recoOrder) if (!assignedReco.count(r)) {
    bool complete = false;
    for (std::size_t i = 0; i < event.truthJetRadiusCode.size(); ++i)
      if (event.truthJetRadiusCode[i] == static_cast<int>(radiusCode(reco[r].radius)))
        complete = event.truthJetContainerValid.at(i) == 1;
    append(r, truth.size(), nearestReco[r], complete ? AssociationState::Fake : AssociationState::Unknown);
  }
  for (auto t : truthOrder) if (!assignedTruth.count(t)) {
    bool observed = false, complete = true;
    for (std::size_t i = 0; i < event.recoJetRadiusCode.size(); ++i)
      if (event.recoJetRadiusCode[i] == static_cast<int>(radiusCode(truth[t].radius))) {
        observed = true; complete &= event.recoJetCaptureState.at(i) == static_cast<int>(CaptureState::Complete);
      }
    append(reco.size(), t, nearestTruth[t], observed && complete ? AssociationState::Miss : AssociationState::Unknown);
  }
  return result;
}


/**
 * Native truth-photon identity inside one event:
 *
 *   embedding id + Geant4 track id
 */
using TruthPhotonKey = std::pair<int, int>;

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

  // The complete truth denominator is captured before any matching.
  captureTruthVertices(topNode);
  captureTruthPhotons(topNode);
  captureTruthJets(topNode);
}


// ============================================================================
// Truth vertices
//
// Store the complete Geant4 vertex census together with explicit event-level
// witnesses for the hard-scatter and minimum-bias vertices.
//
// The PHG4TruthInfoContainer primary vertex defines the hard-scatter witness.
// In a double-interaction event, the first other finite primary vertex is
// retained as the minimum-bias companion.
// ============================================================================

void PhotonJetTree::captureTruthVertices(PHCompositeNode* topNode)
{
  auto* truth =
      findNode::getClass<PHG4TruthInfoContainer>(
          topNode,
          m_config.nodes.truthInfo);

  if (!truth)
  {
    m_event.truthVertexCaptureState =
        CaptureState::InputUnavailable;

    return;
  }

  m_event.truthVertexCaptureState =
      CaptureState::Complete;

  // --------------------------------------------------------------------------
  // Map native Geant4 vertex ids onto embedding ids.
  // --------------------------------------------------------------------------

  std::map<int, int> embeddingByVertex;

  {
    const auto embedded =
        truth->GetEmbeddedVtxIds();

    for (auto it = embedded.first;
         it != embedded.second;
         ++it)
    {
      embeddingByVertex.emplace(
          it->first,
          it->second);
    }
  }

  // --------------------------------------------------------------------------
  // Hard-scatter vertex witness.
  // --------------------------------------------------------------------------

  const int primaryVertexId =
      truth->GetPrimaryVertexIndex();

  m_event.truthPrimaryVertexId =
      primaryVertexId;

  if (const PHG4VtxPoint* primary =
          truth->GetVtx(primaryVertexId);
      primary &&
      finiteValue(primary->get_z()))
  {
    m_event.truthHardVertexZ =
        primary->get_z();

    m_event.truthHardVertexValid =
        true;
  }

  // --------------------------------------------------------------------------
  // Complete truth-vertex census.
  // --------------------------------------------------------------------------

  const auto vertices =
      truth->GetVtxRange();

  for (auto it = vertices.first;
       it != vertices.second;
       ++it)
  {
    TruthVertexRecord row;

    row.eventId = m_event.eventId;
    row.vertexId = it->first;

    const auto embedding =
        embeddingByVertex.find(it->first);

    row.embeddingValid =
        embedding != embeddingByVertex.end();

    row.embeddingId =
        row.embeddingValid
            ? embedding->second
            : 0;

    row.z =
        it->second
            ? it->second->get_z()
            : kNaN;

    row.valid =
        it->second &&
        finiteValue(row.z);

    m_truthVertices.push_back(row);
  }

  // --------------------------------------------------------------------------
  // Minimum-bias companion vertex witness.
  // --------------------------------------------------------------------------

  const auto primaries =
      truth->GetPrimaryVtxRange();

  for (auto it = primaries.first;
       it != primaries.second;
       ++it)
  {
    if (it->first == primaryVertexId ||
        !it->second ||
        !finiteValue(it->second->get_z()))
    {
      continue;
    }

    m_event.truthMinimumBiasVertexZ =
        it->second->get_z();

    m_event.truthMinimumBiasVertexValid =
        true;

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
  auto* truth =
      findNode::getClass<PHG4TruthInfoContainer>(
          topNode,
          m_config.nodes.truthInfo);

  if (!truth)
  {
    m_event.truthPhotonCaptureState =
        CaptureState::InputUnavailable;

    m_event.truthDenominatorComplete =
        false;

    return;
  }

  m_event.truthPhotonCaptureState =
      CaptureState::Complete;

  auto* generatorEvents =
      findNode::getClass<PHHepMCGenEventMap>(
          topNode,
          m_config.nodes.hepmcEventMap);

  const double signalRadius =
      m_config.photon.truthIsolationSignalRadius;

  std::set<TruthPhotonKey> seen;

  const auto particles =
      truth->GetPrimaryParticleRange();

  for (auto it = particles.first;
       it != particles.second;
       ++it)
  {
    PHG4Particle* particle =
        it->second;

    if (!particle)
    {
      ++m_event.truthPhotonIsolationIncompleteCount;
      continue;
    }

    const int trackId =
        particle->get_track_id();

    const int embeddingId =
        truth->isEmbeded(trackId);

    /*
     * Truth-photon population:
     *
     *   - embedded primary subevent
     *   - photon PID
     *
     * No reconstructed-object requirement enters this census.
     */
    if (embeddingId < 1 ||
        particle->get_pid() != 22)
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

    if (!kinematics(
            particle,
            pt,
            eta,
            phi))
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
        calculateTruthPhotonIsolation(
            topNode,
            particle,
            isolationR03,
            isolationR04,
            isolationComplete);

    if (!isolationComplete)
    {
      ++m_event.truthPhotonIsolationIncompleteCount;
      ++m_event.truthPhotonRejectedIsolationCount;
    }

    // ------------------------------------------------------------------------
    // Generator association and prompt classification
    // ------------------------------------------------------------------------

    const int barcode =
        particle->get_barcode();

    const HepMC::GenParticle* generatorPhoton =
        generatorParticle(
            generatorEvent(
                generatorEvents,
                embeddingId),
            barcode);

    const TruthPhotonClass photonClass =
        promptClass(generatorPhoton);

    /*
     * In archived double-interaction production, the hard subevent alone owns
     * the signal role. Source ownership is stored separately from the native
     * photon truth-signal predicate below.
     */
    const bool ownership =
        !isArchivedDoubleInteraction() ||
        embeddingId == kDoubleInteractionHardEmbedding;

    const double signalIsolation =
        std::abs(
            signalRadius -
            kTruthIsolationR03) <
                kRadiusTolerance
            ? isolationR03
            : isolationR04;

    /*
     * Native truth-photon signal predicate.
     *
     * This records the established producer reference definition. Source-role
     * ownership remains separate in signalSourceRole.
     */
    const bool analysisSignal =
        std::abs(eta) < m_config.photon.maxAbsEta &&
        photonjet::isAnalysisSignalClass(photonClass) &&
        isolationValid &&
        finiteValue(signalIsolation) &&
        signalIsolation < m_config.photon.truthIsolationMaxEtGeV;

    // ------------------------------------------------------------------------
    // Persisted truth-photon record
    // ------------------------------------------------------------------------

    TruthPhotonRecord row;

    row.eventId =
        m_event.eventId;

    row.truthPhotonId =
        photonjet::makeTruthPhotonIdentity(
            m_event.eventId,
            embeddingId,
            trackId);

    row.trackId =
        trackId;

    row.vertexId =
        particle->get_vtx_id();

    row.barcode =
        barcode;

    row.embeddingId =
        embeddingId;

    row.pid =
        22;

    row.pt =
        pt;

    row.eta =
        eta;

    row.phi =
        phi;

    row.geantValid =
        true;

    row.generatorAssociationValid =
        generatorPhoton &&
        generatorPhoton->pdg_id() == 22;

    row.promptClass =
        photonClass;

    row.isolationR03 =
        isolationR03;

    row.isolationR04 =
        isolationR04;

    row.isolationValid =
        isolationValid;

    row.analysisSignal =
        analysisSignal;

    row.sampleRole =
        m_config.simulationRole;

    row.signalSourceRole =
        ownership;  // separate from the native photon signal predicate

    m_truthPhotons.push_back(
        std::move(row));

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
  m_event.truthDenominatorComplete =
      m_event.truthPhotonDuplicateTrackCount == 0 &&
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

bool PhotonJetTree::calculateTruthPhotonIsolation(
    PHCompositeNode* topNode,
    const PHG4Particle* photon,
    double& isolationR03,
    double& isolationR04,
    bool& complete) const
{
  isolationR03 = isolationR04 = kNaN;
  complete = false;

  auto* truth =
      findNode::getClass<PHG4TruthInfoContainer>(
          topNode,
          m_config.nodes.truthInfo);

  if (!truth || !photon)
  {
    return false;
  }

  double photonPt, photonEta, photonPhi;

  if (!kinematics(
          photon,
          photonPt,
          photonEta,
          photonPhi))
  {
    return false;
  }

  /*
   * RecoilObjects.cc:17735-17779 uses TLorentzVector::Et, not pT. In
   * particular massive primary particles have a different transverse energy.
   */
  const TLorentzVector axis(
      photon->get_px(),
      photon->get_py(),
      photon->get_pz(),
      photon->get_e());

  const double core =
      m_config.photon.truthIsolationCoreRadius;

  double sumR03 = 0.0;
  double sumR04 = 0.0;
  double coreSum = 0.0;

  complete = true;

  const auto particles =
      truth->GetPrimaryParticleRange();

  for (auto it = particles.first;
       it != particles.second;
       ++it)
  {
    const PHG4Particle* particle =
        it->second;

    if (!particle)
    {
      complete = false;
      continue;
    }

    if (truth->isEmbeded(
            particle->get_track_id()) < 1)
    {
      continue;
    }

    const TLorentzVector momentum(
        particle->get_px(),
        particle->get_py(),
        particle->get_pz(),
        particle->get_e());

    const double et =
        momentum.Et();

    if (!finiteValue(et) || et <= 0.0)
    {
      if (!finiteValue(et) || et < 0.0)
      {
        complete = false;
      }

      continue;
    }

    const double dr =
        axis.DeltaR(momentum);

    if (!finiteValue(dr))
    {
      complete = false;
      continue;
    }

    if (dr < kTruthIsolationR04) sumR04 += et;
    if (dr < kTruthIsolationR03) sumR03 += et;
    if (dr < core) coreSum += et;
  }

  isolationR03 =
      sumR03 - coreSum;

  isolationR04 =
      sumR04 - coreSum;

  return finiteValue(isolationR03) &&
         finiteValue(isolationR04);
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
  m_event.truthJetCaptureState =
      CaptureState::Complete;

  std::set<std::pair<int, std::string>> visited;

  for (const JetNodeConfig& collection :
       m_config.jets.nodes)
  {
    if (collection.truthNode.empty())
    {
      continue;
    }

    const int radiusCode =
        static_cast<int>(
            std::lround(
                collection.radius * 100.0));

    if (!visited.insert(
            {radiusCode, collection.truthNode})
             .second)
    {
      continue;  // several reconstructed views share one truth collection
    }

    JetContainer* container =
        findNode::getClass<JetContainer>(
            topNode,
            collection.truthNode);

    m_event.truthJetRadiusCode.push_back(
        radiusCode);

    m_event.truthJetContainerValid.push_back(
        container ? 1 : 0);

    if (!container)
    {
      m_event.truthJetCaptureState =
          CaptureState::InputUnavailable;
    }

    if (!container)
    {
      continue;
    }

    for (Jet* jet : *container)
    {
      if (!jet)
      {
        m_event.truthJetContainerValid.back() = 0;
        continue;
      }

      const double pt =
          jet->get_pt();

      const double eta =
          jet->get_eta();

      const double phi =
          jet->get_phi();

      if (!finiteValue(pt) ||
          !finiteValue(eta) ||
          !finiteValue(phi) ||
          pt < 0.0)
      {
        m_event.truthJetContainerValid.back() = 0;
        continue;
      }

      TruthJetRecord row;

      row.eventId =
          m_event.eventId;

      row.radius =
          collection.radius;

      row.nativeKey =
          jet->get_id();

      row.truthJetId =
          photonjet::makeTruthJetIdentity(
              m_event.eventId,
              collection.radius,
              row.nativeKey);

      row.pt =
          pt;

      row.eta =
          eta;

      row.phi =
          phi;

      row.containerValid =
          true;

      m_truthJets.push_back(
          std::move(row));
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
  if (!isSimulation()) return;
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode,m_config.nodes.truthInfo);
  auto* clusters = findNode::getClass<RawClusterContainer>(topNode,m_config.photon.rawClusterNode);

  CaloRawClusterEval towerInfoEvaluator(topNode,"CEMC");
  towerInfoEvaluator.set_strict(false);
  towerInfoEvaluator.do_caching(true);
  towerInfoEvaluator.set_usetowerinfo(true);
  CaloRawClusterEval legacyEvaluator(topNode,"CEMC");
  legacyEvaluator.set_strict(false);
  legacyEvaluator.do_caching(true);
  legacyEvaluator.set_usetowerinfo(false);

  // Capture the maximum-energy primary once; the record helper below owns
  // every association and selection decision. A photon-only search cannot
  // reproduce this detector witness.
  for (auto& photon : m_photons) {
    auto& witness = photon.dominantTruth;
    witness = DominantTruthRecord{};
    witness.state = DominantTruthState::EvaluatorUnavailable;
    witness.evaluator = DominantTruthEvaluator::Unavailable;
    RawCluster* cluster = clusters ? clusters->getCluster(static_cast<unsigned int>(photon.nativeClusterKey)) : nullptr;
    if (!truth || !cluster || (!towerInfoEvaluator.has_reduced_node_pointers() &&
                              !legacyEvaluator.has_reduced_node_pointers())) continue;
    auto* evaluator = &towerInfoEvaluator;
    witness.evaluator = DominantTruthEvaluator::TowerInfo;
    PHG4Particle* primary = towerInfoEvaluator.has_reduced_node_pointers()
        ? towerInfoEvaluator.max_truth_primary_particle_by_energy(cluster) : nullptr;
    if (!primary && legacyEvaluator.has_reduced_node_pointers()) {
      evaluator = &legacyEvaluator;
      witness.evaluator = DominantTruthEvaluator::LegacyRawTower;
      primary = legacyEvaluator.max_truth_primary_particle_by_energy(cluster);
    }
    witness.state = DominantTruthState::NoPrimary;
    if (!primary) continue;
    witness.trackId = primary->get_track_id();
    witness.pid = primary->get_pid();
    witness.barcode = primary->get_barcode();
    witness.vertexId = primary->get_vtx_id();
    witness.energyContribution = evaluator->get_energy_contribution(cluster,primary);
    witness.state = DominantTruthState::InvalidPrimary;
    if (witness.trackId <= 0 || witness.pid == 0 ||
        !std::isfinite(witness.energyContribution) || witness.energyContribution < 0.) continue;
    witness.embeddingId = truth->isEmbeded(witness.trackId);
    witness.state = DominantTruthState::ValidPrimary;
  }
  std::vector<bool> selected;
  for (const auto& row : m_truthPhotons) selected.push_back(row.analysisSignal && row.signalSourceRole);
  m_photonTruthLinks = matchPhotonRecords(m_event,m_photons,m_truthPhotons,selected);
}

// ============================================================================
// Reconstructed-jet -> truth-jet relations
//
// Same-radius edges share one truth assignment across the retained views.
// Candidate edges, selected edges and unmatched objects remain separate rows.
// Missing opposite-side capture produces Unknown rather than a fake or miss.
// ============================================================================

void PhotonJetTree::buildJetTruthAssociations()
{
  m_jetTruthLinks.clear();
  if (!isSimulation()) return;
  m_jetTruthLinks = matchJetRecords(m_event, m_jets, m_truthJets, m_config.jets.truthMatchMaxDeltaR);
}
