// Truth side of the photon + recoil-jet response on eventTree v1 simulation: the truth photon, its truth recoil
// jets, the reco-truth jet match, the nearby-truth-jet test of combinatoric reco jets, and the response categories.
//
// Companion of PhotonJetSelection.h, which holds every reco decision. The truth recoil cuts are the reco recoil
// cuts of the same Profile. Functions take plain values, or indexable arrays (TTreeReaderArray, std::vector) with
// their length, so the ROOT reader (FillPhotonJetResponse.C) and any other caller make identical decisions.
//
// Order: eventTree collections keep the row order of the production tables, and a row identity is (domain,
// row + 1), so within an event "identity order" is "index order". The reference breaks ties by identity; the
// functions below break them by index, which is the same decision.
#pragma once

#include "PhotonJetSelection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace photonjet {

// Truth photon acceptance of the reference response: |eta| below this. The reference keeps it apart from the reco
// photon window (Profile::photonAbsEtaMax); both are 0.7.
constexpr double kTruthPhotonAbsEtaMax = 0.7;
// Analysis jet radius (collections jet04 and truthjet04). The nearby-truth-jet test uses it as its Delta R.
constexpr double kAnalysisJetRadius = 0.4;
// truthpho_source_role of an analysis-signal truth photon.
constexpr int kSignalSourceRole = 1;
// The analysis truth signal, recorded by the producer as source role 1: a valid Geant4 photon of a prompt class below
// hadronic, |eta| < kTruthPhotonAbsEtaMax, and truth isolation in R = 0.3 below 4 GeV.
constexpr int kHadronicPromptClass = 3;
constexpr double kTruthIsolationMaxGeV = 4.0;
inline bool analysisTruthSignal(int g4Valid, int promptClass, double eta, double isolationR03) {
  return g4Valid == 1 && promptClass < kHadronicPromptClass && std::fabs(eta) < kTruthPhotonAbsEtaMax &&
         isolationR03 < kTruthIsolationMaxGeV;
}
// plink_class / jlink_class of a final reco-truth match.
constexpr int kLinkMatch = 0;
// jlink_class of a candidate edge: a reco and a truth jet of the same radius within the producer's match distance.
constexpr int kLinkMatchCandidate = 5;

// An analysis event enters the response when its vertex is finite and inside the profile window. This applies to
// every event, with or without a reco photon, before the truth side is filled.
inline bool responseEvent(const Profile& p, double vz) { return std::isfinite(vz) && std::fabs(vz) < p.vertexAbsZMax; }

// The truth photon of the event: the highest-pT truth photon with source role 1, finite kinematics, pT > 0 and
// |eta| < kTruthPhotonAbsEtaMax. Returns its index, or -1 when there is none.
template <class Values, class Ints>
long truthPhoton(const Values& pt, const Values& eta, const Values& phi, const Ints& sourceRole, std::size_t n) {
  long best = -1;
  for (std::size_t t = 0; t < n; ++t) {
    if (sourceRole[t] != kSignalSourceRole) continue;
    if (!(std::isfinite(pt[t]) && std::isfinite(eta[t]) && std::isfinite(phi[t]))) continue;
    if (!(pt[t] > 0.0 && std::fabs(eta[t]) < kTruthPhotonAbsEtaMax)) continue;
    // an exact pT tie goes to the later photon: the reference takes the maximum of (pT, identity)
    if (best < 0 || pt[t] >= pt[best]) best = static_cast<long>(t);
  }
  return best;
}

// A truth recoil jet of the truth photon: the reco recoil cuts of the profile (recoilJet) applied to a truth jet.
// Both angles take the jet first, as the reference does; the wrapped Delta phi is not bitwise symmetric in its
// arguments, so the order is kept.
inline bool truthRecoilJet(const Profile& p, double photonEta, double photonPhi, double pt, double eta, double phi) {
  if (!(std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi))) return false;
  if (pt < p.jetPtMin) return false;
  if (deltaR(eta, phi, photonEta, photonPhi) < p.overlapDeltaR) return false;
  if (std::fabs(eta) >= p.jetAbsEtaMax) return false;
  return wrappedAbsDeltaPhi(phi, photonPhi) >= p.deltaPhiMin;
}

// The indices i < n with accept(i), by descending pt[i]; equal pT keep index order. This is the reference's order of
// recoil jets: (-pT, identity) for truth jets, and (-pT, deterministic_order, identity) for reco jets, whose
// deterministic_order is the index within the collection.
template <class Accept, class Values>
std::vector<long> orderedByPt(std::size_t n, Accept accept, const Values& pt) {
  std::vector<long> out;
  for (std::size_t i = 0; i < n; ++i)
    if (accept(i)) out.push_back(static_cast<long>(i));
  std::stable_sort(out.begin(), out.end(), [&pt](long a, long b) { return pt[a] > pt[b]; });
  return out;
}

// A jet candidate edge: event-local reco and truth jet indices, its Delta R and the reco jet's corrected pT.
struct JetEdge {
  std::size_t reco, truth;
  double deltaR, recoPt;
};

// The producer's one-to-one jet match over a set of candidate edges: edges in order of Delta R, then reco pT
// descending, reco index, truth index (index order is identity order); an edge is taken when neither of its jets is
// taken yet. Returns the (reco, truth) pairs taken, sorted.
inline std::vector<std::pair<std::size_t, std::size_t>> greedyJetMatch(std::vector<JetEdge> edges) {
  std::sort(edges.begin(), edges.end(), [](const JetEdge& a, const JetEdge& b) {
    if (a.deltaR != b.deltaR) return a.deltaR < b.deltaR;
    if (a.recoPt != b.recoPt) return a.recoPt > b.recoPt;
    if (a.reco != b.reco) return a.reco < b.reco;
    return a.truth < b.truth;
  });
  std::vector<std::pair<std::size_t, std::size_t>> taken;
  auto free = [&taken](std::size_t reco, std::size_t truth) {
    for (const auto& pair : taken)
      if (pair.first == reco || pair.second == truth) return false;
    return true;
  };
  for (const JetEdge& e : edges)
    if (free(e.reco, e.truth)) taken.emplace_back(e.reco, e.truth);
  std::sort(taken.begin(), taken.end());
  return taken;
}

// True when a truth jet (truthjet04) above the threshold lies within matchDeltaR of the reco jet. A selected reco
// recoil of a matched photon without such a jet is combinatoric. Fiducial and recoil cuts are not applied to the
// truth jet: they belong to the observable, not to the label.
template <class Values>
bool nearbyTruthJet(double threshold, double recoEta, double recoPhi, const Values& pt, const Values& eta,
                    const Values& phi, std::size_t n, double matchDeltaR = kAnalysisJetRadius) {
  for (std::size_t t = 0; t < n; ++t) {
    if (!(std::isfinite(pt[t]) && std::isfinite(eta[t]) && std::isfinite(phi[t]))) continue;
    if (pt[t] >= threshold && deltaR(eta[t], phi[t], recoEta, recoPhi) < matchDeltaR) return true;
  }
  return false;
}

// Photon-level response category. A matched pair (the region-A photon's final match is the truth photon) fills the
// response when both photons are inside their pT axes. Otherwise the reco photon is a fake and the truth photon a
// miss, each only when inside its axis; so a matched photon whose truth pT is outside the truth axis is a reco fake.
struct PhotonResponse {
  bool matched, fake, miss;
};

inline PhotonResponse photonResponse(bool linked, bool recoInAxis, bool truthInAxis) {
  if (linked && recoInAxis && truthInAxis) return {true, false, false};
  return {false, recoInAxis, truthInAxis};
}

// Category of one selected reco recoil jet (inside the xJ axis) of the region-A photon:
//   kMatched       photon matched, and the jet's final match is a selected truth recoil not yet taken
//   kFakePhoton    photon not matched
//   kCombinatoric  photon matched, no truth jet nearby (nearbyTruthJet)
//   kDetectorFake  photon matched, a truth jet nearby, but no free linked truth recoil
// Reco recoils are classified in orderedByPt order and each truth recoil is taken by the first jet linked to it; a
// later jet linked to the same truth recoil is a fake. The result therefore depends on that order. truthJetNearby
// is only called for a matched photon without a match, as in the reference.
enum class RecoilClass { kMatched, kFakePhoton, kCombinatoric, kDetectorFake };

template <class NearbyTest>
RecoilClass recoilClass(bool photonMatched, bool linkedToFreeTruthRecoil, NearbyTest truthJetNearby) {
  if (photonMatched && linkedToFreeTruthRecoil) return RecoilClass::kMatched;
  if (!photonMatched) return RecoilClass::kFakePhoton;
  return truthJetNearby() ? RecoilClass::kDetectorFake : RecoilClass::kCombinatoric;
}

// Signal leakage: one region per signal truth photon. The leader of region r (A, B, C, D) is a witness of the truth
// photon of its final match when that photon has source role 1 (any truth photon of the event; no acceptance cut).
// A truth photon keeps its smallest witness by (match metric, -ET, region, photon index). Returns which regions'
// leaders are kept. photonTruth / photonMetric: final match of each reco photon (-1 / NaN when none).
// Throws when a witness has no finite match metric, as the reference does.
template <class Values, class Ints>
std::array<bool, 4> leakageWitnesses(const std::array<long, 4>& leader, const std::vector<long>& photonTruth,
                                     const std::vector<double>& photonMetric, const Values& et,
                                     const Ints& truthSourceRole) {
  using Order = std::tuple<double, double, int, long>;
  std::array<long, 4> truth{-1, -1, -1, -1};
  std::array<Order, 4> order{};
  for (int r = 0; r < 4; ++r) {
    const long i = leader[r];
    if (i < 0 || photonTruth[i] < 0 || truthSourceRole[photonTruth[i]] != kSignalSourceRole) continue;
    if (!std::isfinite(photonMetric[i])) throw std::runtime_error("final photon match without a finite metric");
    truth[r] = photonTruth[i];
    order[r] = Order(photonMetric[i], -et[i], r, i);
  }
  std::array<bool, 4> keep{false, false, false, false};
  for (int r = 0; r < 4; ++r) {
    if (truth[r] < 0) continue;
    keep[r] = true;
    for (int s = 0; s < 4; ++s)
      if (s != r && truth[s] == truth[r] && order[s] < order[r]) keep[r] = false;
  }
  return keep;
}

}  // namespace photonjet
