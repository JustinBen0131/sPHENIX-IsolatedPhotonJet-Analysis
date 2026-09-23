//
// Truth.cc
//
// Simulation truth and its relation to reconstruction.
//
// Owns
//   the Geant4 vertex census and the hard and minimum-bias vertex witnesses
//   the embedded primary photon inventory with generator association,
//   prompt classification and truth isolation
//   truth jets
//   the dominant-primary witness of every reconstructed photon and the
//   photon-to-truth relation built from it
//   the deterministic jet-to-truth-jet relation
//
// Does not own
//   truth-jet clustering                 -> Production.cc
//   reconstructed photons and jets       -> Photons.cc, Jets.cc
//   ROOT tables                          -> Output.cc
//
// Rules that must not drift
//   Truth capture is selection neutral: every embedded primary photon is
//   written, and the reasons for any exclusion are counted.
//   Geant4 identity and generator identity are kept separately; a missing
//   generator association is recorded, never invented.
//   The prompt class is the production-vertex walk; the signal predicate
//   admits every class below Hadronic, including Unknown and Unclassified.
//   Truth isolation sums all embedded primaries in the cone, not photons
//   only, with the merged core removed.
//   The photon relation is by energy-dominant primary identity, never by
//   nearest separation; the stored separation is diagnostic.
//   Truth jets carry no energy scale. Jet matching is one-to-one, within the
//   configured radius, with a deterministic edge order.
//   Miss is asserted only when reconstructed capture is known complete.
//
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

// The embedding identifier of the hard-scatter subevent in the archived
// double-interaction samples.
constexpr int kDoubleInteractionHardEmbedding = 2;

bool finite(const double value)
{
  return std::isfinite(value);
}

double wrapPhi(double phi)
{
  while (phi > kPi) phi -= 2.0 * kPi;
  while (phi <= -kPi) phi += 2.0 * kPi;
  return phi;
}

double deltaR(const double eta1, const double phi1, const double eta2, const double phi2)
{
  return std::hypot(eta1 - eta2, wrapPhi(phi1 - phi2));
}

// Transverse momentum, pseudorapidity and azimuth of a Geant4 particle.
// Returns false for a non-finite or non-positive transverse momentum.
bool kinematics(const PHG4Particle* particle, double& pt, double& eta, double& phi)
{
  pt = eta = phi = kNaN;
  const double px = particle->get_px();
  const double py = particle->get_py();
  const double pz = particle->get_pz();
  if (!finite(px) || !finite(py) || !finite(pz))
  {
    return false;
  }
  pt = std::hypot(px, py);
  if (!finite(pt) || !(pt > 0.0))
  {
    return false;
  }
  eta = std::asinh(pz / pt);
  phi = std::atan2(py, px);
  return finite(eta) && finite(phi);
}

// ---- generator association ------------------------------------------------------

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

// The prompt class from the production-vertex walk.
//
// Walk backwards through single-photon continuations until the first
// vertex with more than a photon-to-photon relation. Two incoming and two
// outgoing legs, all with |PDG| at most 22, is a direct photon. One incoming
// lepton or quark with |PDG| at most 11 that emerges beside a photon is a
// fragmentation photon. One incoming particle with |PDG| above 37 is a
// hadronic decay. Anything else is unclassified; a missing history is
// unknown. No generator status code enters this decision.
// Exact predecessor classifier; only the malformed-cycle guard is added.
int promptClassCode(const HepMC::GenParticle* pho)
{
  if (!pho || pho->pdg_id() != 22) return -1;

  const HepMC::GenVertex* vertex = pho->production_vertex();
  if (!vertex) return -1;

  std::vector<const HepMC::GenParticle*> incomingParticles;
  for (auto inItr = vertex->particles_in_const_begin();
       inItr != vertex->particles_in_const_end(); ++inItr)
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
         inItr != vertex->particles_in_const_end(); ++inItr)
    {
      if (*inItr) incomingParticles.push_back(*inItr);
    }
  }

  std::vector<const HepMC::GenParticle*> outgoingParticles;
  for (auto outItr = vertex->particles_out_const_begin();
       outItr != vertex->particles_out_const_end(); ++outItr)
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
    if (std::abs(in0) <= 22 && std::abs(in1) <= 22 &&
        std::abs(out0) <= 22 && std::abs(out1) <= 22)
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

// ---- jet matching ---------------------------------------------------------------

struct MatchEdge
{
  std::size_t recoIndex = 0;
  std::size_t truthIndex = 0;
  double deltaR = kNaN;
  double recoPt = kNaN;
};

// Smallest separation first, then highest reconstructed momentum, then the
// identities, so the greedy assignment is reproducible.
bool edgeBefore(const MatchEdge& a, const MatchEdge& b,
                const std::vector<photonjet::JetRecord>& reco,
                const std::vector<photonjet::TruthJetRecord>& truth)
{
  if (a.deltaR != b.deltaR) return a.deltaR < b.deltaR;
  if (a.recoPt != b.recoPt) return a.recoPt > b.recoPt;
  const auto& ra = reco[a.recoIndex].jetId;
  const auto& rb = reco[b.recoIndex].jetId;
  if (ra.hi != rb.hi) return ra.hi < rb.hi;
  if (ra.lo != rb.lo) return ra.lo < rb.lo;
  const auto& ta = truth[a.truthIndex].truthJetId;
  const auto& tb = truth[b.truthIndex].truthJetId;
  if (ta.hi != tb.hi) return ta.hi < tb.hi;
  return ta.lo < tb.lo;
}

using TruthPhotonKey = std::pair<int, int>;  // embedding id, Geant4 track id

}  // namespace

//
// Truth stage
//
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

//
// Truth vertices
//
// The hard-scatter vertex is the container's primary vertex. In a
// double-interaction sample the first other primary vertex is the
// minimum-bias companion. The full census is kept beside these witnesses.
//
void PhotonJetTree::captureTruthVertices(PHCompositeNode* topNode)
{
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);
  if (!truth)
  {
    m_event.truthVertexCaptureState = CaptureState::InputUnavailable;
    return;
  }
  m_event.truthVertexCaptureState = CaptureState::Complete;

  std::map<int, int> embeddingByVertex;
  {
    const auto embedded = truth->GetEmbeddedVtxIds();
    for (auto it = embedded.first; it != embedded.second; ++it)
    {
      embeddingByVertex.emplace(it->first, it->second);
    }
  }

  const int primaryVertexId = truth->GetPrimaryVertexIndex();
  m_event.truthPrimaryVertexId = primaryVertexId;
  if (const PHG4VtxPoint* primary = truth->GetVtx(primaryVertexId); primary && finite(primary->get_z()))
  {
    m_event.truthHardVertexZ = primary->get_z();
    m_event.truthHardVertexValid = true;
  }

  const auto vertices = truth->GetVtxRange();
  for (auto it = vertices.first; it != vertices.second; ++it)
  {
    TruthVertexRecord row;
    row.eventId = m_event.eventId;
    row.vertexId = it->first;
    const auto embedding = embeddingByVertex.find(it->first);
    row.embeddingValid = embedding != embeddingByVertex.end();
    row.embeddingId = row.embeddingValid ? embedding->second : 0;
    row.z = it->second ? it->second->get_z() : kNaN;
    row.valid = it->second && finite(row.z);
    m_truthVertices.push_back(row);
  }

  const auto primaries = truth->GetPrimaryVtxRange();
  for (auto it = primaries.first; it != primaries.second; ++it)
  {
    if (it->first == primaryVertexId || !it->second || !finite(it->second->get_z()))
    {
      continue;
    }
    m_event.truthMinimumBiasVertexZ = it->second->get_z();
    m_event.truthMinimumBiasVertexValid = true;
    break;
  }
}

//
// Truth photons
//
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

    // The truth population is the embedded primary subevents.
    if (embeddingId < 1 || particle->get_pid() != 22)
    {
      continue;
    }
    ++m_event.truthPhotonEmbeddedPrimaryCount;

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

    double isolationR03 = kNaN, isolationR04 = kNaN;
    bool isolationComplete = false;
    const bool isolationValid = calculateTruthPhotonIsolation(topNode, particle, isolationR03, isolationR04, isolationComplete);
    if (!isolationComplete)
    {
      ++m_event.truthPhotonIsolationIncompleteCount;
      ++m_event.truthPhotonRejectedIsolationCount;
    }

    const int barcode = particle->get_barcode();
    const HepMC::GenParticle* generatorPhoton =
        generatorParticle(generatorEvent(generatorEvents, embeddingId), barcode);
    const TruthPhotonClass photonClass = promptClass(generatorPhoton);

    // In an archived double-interaction sample only the hard subevent owns
    // the signal.
    const bool ownership = !isArchivedDoubleInteraction() || embeddingId == kDoubleInteractionHardEmbedding;
    const double signalIsolation = std::abs(signalRadius - kTruthIsolationR03) < kRadiusTolerance
                                       ? isolationR03 : isolationR04;

    const bool analysisSignal =
                                std::abs(eta) < m_config.photon.maxAbsEta &&
                                photonjet::isAnalysisSignalClass(photonClass) &&
                                isolationValid && finite(signalIsolation) &&
                                signalIsolation < m_config.photon.truthIsolationMaxEtGeV;

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

  m_event.truthDenominatorComplete =
      m_event.truthPhotonDuplicateTrackCount == 0 && m_event.truthPhotonIsolationIncompleteCount == 0 &&
      m_event.truthPhotonRejectedKinematicsCount == 0;
}

// Truth isolation uses embedded primary transverse energies (TLorentzVector::Et).
// Finite sums and complete input coverage are recorded independently.
bool PhotonJetTree::calculateTruthPhotonIsolation(PHCompositeNode* topNode,
                                                  const PHG4Particle* photon,
                                                  double& isolationR03,
                                                  double& isolationR04,
                                                  bool& complete) const
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

  // RecoilJets.cc:17735-17779 uses TLorentzVector::Et, not pT. In
  // particular massive primary particles have a different transverse energy.
  const TLorentzVector axis(photon->get_px(), photon->get_py(), photon->get_pz(), photon->get_e());
  const double core = m_config.photon.truthIsolationCoreRadius;
  double sumR03 = 0.0, sumR04 = 0.0, coreSum = 0.0;
  complete = true;
  const auto particles = truth->GetPrimaryParticleRange();
  for (auto it = particles.first; it != particles.second; ++it)
  {
    const PHG4Particle* particle = it->second;
    if (!particle) { complete = false; continue; }
    if (truth->isEmbeded(particle->get_track_id()) < 1) continue;
    const TLorentzVector momentum(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e());
    const double et = momentum.Et();
    if (!finite(et) || et <= 0.0)
    {
      if (!finite(et) || et < 0.0) complete = false;
      continue;
    }
    const double dr = axis.DeltaR(momentum);
    if (!finite(dr)) { complete = false; continue; }
    if (dr < kTruthIsolationR04) sumR04 += et;
    if (dr < kTruthIsolationR03) sumR03 += et;
    if (dr < core) coreSum += et;
  }
  isolationR03 = sumR03 - coreSum;
  isolationR04 = sumR04 - coreSum;
  return finite(isolationR03) && finite(isolationR04);
}

//
// Truth jets
//
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
    if (!container) m_event.truthJetCaptureState = CaptureState::InputUnavailable;
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
      const double pt = jet->get_pt();
      const double eta = jet->get_eta();
      const double phi = jet->get_phi();
      if (!finite(pt) || !finite(eta) || !finite(phi) || pt < 0.0)
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

//
// Photon-to-truth relation
//
// For every reconstructed photon the calorimeter truth evaluator names the
// Geant4 primary that deposited the most energy in its cluster. That primary
// is looked up in the truth-photon inventory by embedding and track
// identity. A primary that is not a retained truth photon makes the
// candidate a fake, unless it is a photon whose row was not written, which
// stays unknown. Truth photons matched by nobody are misses only when
// reconstructed capture is complete and every candidate was resolved.
//
void PhotonJetTree::buildPhotonTruthAssociations(PHCompositeNode* topNode)
{
  m_photonTruthLinks.clear();
  if (!isSimulation())
  {
    return;
  }

  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_config.nodes.truthInfo);
  auto* clusters = findNode::getClass<RawClusterContainer>(topNode, m_config.photon.rawClusterNode);

  std::map<TruthPhotonKey, std::size_t> truthByKey;
  for (std::size_t index = 0; index < m_truthPhotons.size(); ++index)
  {
    const TruthPhotonRecord& row = m_truthPhotons[index];
    if (!truthByKey.emplace(TruthPhotonKey{row.embeddingId, row.trackId}, index).second)
    {
      throw std::runtime_error("duplicate truth-photon identity in the retained inventory");
    }
  }

  auto pushLink = [&](PhotonTruthLinkRecord link) {
    link.linkId = photonjet::makeLinkIdentity(m_event.eventId, link.photonId, link.truthPhotonId,
                                              LinkKind::PhotonToTruthPhoton);
    m_photonTruthLinks.push_back(std::move(link));
  };

  // The TowerInfo evaluator is preferred; the legacy tower evaluator is the
  // fallback when the first finds no primary.
  CaloRawClusterEval towerInfoEvaluator(topNode, "CEMC");
  towerInfoEvaluator.set_strict(false);
  towerInfoEvaluator.do_caching(true);
  towerInfoEvaluator.set_usetowerinfo(true);

  CaloRawClusterEval legacyEvaluator(topNode, "CEMC");
  legacyEvaluator.set_strict(false);
  legacyEvaluator.do_caching(true);
  legacyEvaluator.set_usetowerinfo(false);

  std::set<std::size_t> matchedTruth;
  bool unresolved = false;

  for (PhotonRecord& photon : m_photons)
  {
    photon.dominantTruth = DominantTruthRecord{};
    photon.dominantTruth.state = DominantTruthState::EvaluatorUnavailable;
    photon.dominantTruth.evaluator = DominantTruthEvaluator::Unavailable;

    PhotonTruthLinkRecord link;
    link.eventId = m_event.eventId;
    link.photonId = photon.photonId;
    link.state = AssociationState::Unknown;

    RawCluster* cluster = clusters ? clusters->getCluster(static_cast<unsigned int>(photon.nativeClusterKey)) : nullptr;
    if (!truth || !cluster || (!towerInfoEvaluator.has_reduced_node_pointers() && !legacyEvaluator.has_reduced_node_pointers()))
    {
      unresolved = true;
      pushLink(link);
      continue;
    }

    CaloRawClusterEval* evaluator = &towerInfoEvaluator;
    photon.dominantTruth.evaluator = DominantTruthEvaluator::TowerInfo;
    PHG4Particle* primary = towerInfoEvaluator.has_reduced_node_pointers()
        ? towerInfoEvaluator.max_truth_primary_particle_by_energy(cluster) : nullptr;
    if (!primary && legacyEvaluator.has_reduced_node_pointers())
    {
      evaluator = &legacyEvaluator;
      photon.dominantTruth.evaluator = DominantTruthEvaluator::LegacyRawTower;
      primary = legacyEvaluator.max_truth_primary_particle_by_energy(cluster);
    }

    if (!primary)
    {
      photon.dominantTruth.state = DominantTruthState::NoPrimary;
      unresolved = true;
      pushLink(link);
      continue;
    }

    photon.dominantTruth.trackId = primary->get_track_id();
    photon.dominantTruth.pid = primary->get_pid();
    photon.dominantTruth.barcode = primary->get_barcode();
    photon.dominantTruth.vertexId = primary->get_vtx_id();
    photon.dominantTruth.energyContribution = evaluator->get_energy_contribution(cluster, primary);

    if (photon.dominantTruth.trackId <= 0 || photon.dominantTruth.pid == 0 ||
        !finite(photon.dominantTruth.energyContribution) || photon.dominantTruth.energyContribution < 0.0)
    {
      photon.dominantTruth.state = DominantTruthState::InvalidPrimary;
      unresolved = true;
      pushLink(link);
      continue;
    }

    photon.dominantTruth.embeddingId = truth->isEmbeded(photon.dominantTruth.trackId);
    photon.dominantTruth.state = DominantTruthState::ValidPrimary;
    link.energyContribution = photon.dominantTruth.energyContribution;

    const auto found = truthByKey.find({photon.dominantTruth.embeddingId, photon.dominantTruth.trackId});
    if (found == truthByKey.end())
    {
      if (photon.dominantTruth.pid == 22)
      {
        // A photon primary without a retained truth row: the inventory was
        // incomplete, so the relation stays unknown rather than fake.
        link.state = AssociationState::Unknown;
        unresolved = true;
      }
      else
      {
        link.state = AssociationState::Fake;
      }
      pushLink(link);
      continue;
    }

    const TruthPhotonRecord& truthPhoton = m_truthPhotons[found->second];
    link.truthPhotonId = truthPhoton.truthPhotonId;
    link.state = AssociationState::Matched;
    if (finite(photon.eta) && finite(photon.phi) && finite(truthPhoton.eta) && finite(truthPhoton.phi))
    {
      link.deltaR = deltaR(photon.eta, photon.phi, truthPhoton.eta, truthPhoton.phi);  // diagnostic only
    }
    pushLink(link);
    matchedTruth.insert(found->second);
  }

  const bool recoCaptureComplete =
      m_event.recoPhotonCaptureState == CaptureState::Complete && m_event.recoPhotonUnclassifiableCount == 0;

  for (std::size_t index = 0; index < m_truthPhotons.size(); ++index)
  {
    if (matchedTruth.count(index) != 0)
    {
      continue;
    }
    PhotonTruthLinkRecord link;
    link.eventId = m_event.eventId;
    link.truthPhotonId = m_truthPhotons[index].truthPhotonId;
    link.state = recoCaptureComplete && !unresolved ? AssociationState::Miss : AssociationState::Unknown;
    pushLink(link);
  }
}

//
// Jet-to-truth-jet relation
//
// Within each reconstructed view and radius, every candidate edge inside
// the matching radius is retained; the edges are sorted deterministically
// and assigned greedily one to one. Unassigned reconstructed jets are fakes
// and unassigned truth jets are misses, relative to that view.
//
void PhotonJetTree::buildJetTruthAssociations()
{
  m_jetTruthLinks.clear();
  if (!isSimulation() || (m_jets.empty() && m_truthJets.empty()))
  {
    return;
  }

  const double maxDeltaR = m_config.jets.truthMatchMaxDeltaR;

  auto pushLink = [&](JetTruthLinkRecord link) {
    link.linkId = photonjet::makeLinkIdentity(m_event.eventId, link.jetId, link.truthJetId, LinkKind::JetToTruthJet);
    link.linkId.hi = photonjet::mixIdentityWord(link.linkId.hi ^ static_cast<std::uint64_t>(link.jetView));
    link.linkId.lo = photonjet::mixIdentityWord(link.linkId.lo ^ photonjet::radiusCode(link.jetRadius));
    m_jetTruthLinks.push_back(std::move(link));
  };

  std::set<std::pair<JetView, int>> views;
  for (const JetRecord& jet : m_jets)
  {
    views.insert({jet.view, static_cast<int>(std::lround(jet.radius * 100.0))});
  }
  for (const JetNodeConfig& collection : m_config.jets.nodes)
  {
    views.insert({collection.view, static_cast<int>(std::lround(collection.radius * 100.0))});
  }

  for (const auto& [view, radiusCode] : views)
  {
    const double radius = radiusCode / 100.0;
    bool recoComplete = false, truthComplete = false;
    for (std::size_t i = 0; i < m_event.recoJetView.size(); ++i)
      if (m_event.recoJetView[i] == static_cast<int>(view) && m_event.recoJetRadiusCode[i] == radiusCode)
        recoComplete = m_event.recoJetCaptureState[i] == static_cast<int>(CaptureState::Complete);
    for (std::size_t i = 0; i < m_event.truthJetRadiusCode.size(); ++i)
      if (m_event.truthJetRadiusCode[i] == radiusCode)
        truthComplete = m_event.truthJetContainerValid[i] == 1;

    std::vector<std::size_t> recoIndices, truthIndices;
    for (std::size_t i = 0; i < m_jets.size(); ++i)
    {
      if (m_jets[i].view == view && std::abs(m_jets[i].radius - radius) < kRadiusTolerance)
      {
        recoIndices.push_back(i);
      }
    }
    for (std::size_t i = 0; i < m_truthJets.size(); ++i)
    {
      if (std::abs(m_truthJets[i].radius - radius) < kRadiusTolerance)
      {
        truthIndices.push_back(i);
      }
    }
    if (recoIndices.empty() && truthIndices.empty())
    {
      continue;
    }

    std::vector<MatchEdge> edges;
    for (const std::size_t r : recoIndices)
    {
      const JetRecord& reco = m_jets[r];
      if (!finite(reco.eta) || !finite(reco.phi) || !finite(reco.correctedPt))
      {
        continue;
      }
      for (const std::size_t t : truthIndices)
      {
        const TruthJetRecord& truthJet = m_truthJets[t];
        const double dr = deltaR(reco.eta, reco.phi, truthJet.eta, truthJet.phi);
        if (finite(dr) && dr < maxDeltaR)
        {
          edges.push_back({r, t, dr, reco.correctedPt});
        }
      }
    }
    std::stable_sort(edges.begin(), edges.end(),
                     [&](const MatchEdge& a, const MatchEdge& b) { return edgeBefore(a, b, m_jets, m_truthJets); });

    std::set<std::size_t> assignedReco, assignedTruth;
    for (const MatchEdge& edge : edges)
    {
      const bool selected = assignedReco.count(edge.recoIndex) == 0 && assignedTruth.count(edge.truthIndex) == 0;
      if (selected)
      {
        assignedReco.insert(edge.recoIndex);
        assignedTruth.insert(edge.truthIndex);
      }
      JetTruthLinkRecord link;
      link.eventId = m_event.eventId;
      link.jetId = m_jets[edge.recoIndex].jetId;
      link.truthJetId = m_truthJets[edge.truthIndex].truthJetId;
      link.jetView = view;
      link.jetRadius = radius;
      link.state = selected ? AssociationState::Matched : AssociationState::Unmatched;
      link.deltaR = edge.deltaR;
      link.selectedMatch = selected;
      pushLink(link);
    }

    for (const std::size_t r : recoIndices)
    {
      if (assignedReco.count(r) != 0) continue;
      JetTruthLinkRecord link;
      link.eventId = m_event.eventId;
      link.jetId = m_jets[r].jetId;
      link.jetView = view;
      link.jetRadius = radius;
      link.state = truthComplete ? AssociationState::Fake : AssociationState::Unknown;
      pushLink(link);
    }
    for (const std::size_t t : truthIndices)
    {
      if (assignedTruth.count(t) != 0) continue;
      JetTruthLinkRecord link;
      link.eventId = m_event.eventId;
      link.truthJetId = m_truthJets[t].truthJetId;
      link.jetView = view;
      link.jetRadius = radius;
      link.state = recoComplete ? AssociationState::Miss : AssociationState::Unknown;
      pushLink(link);
    }
  }
}
