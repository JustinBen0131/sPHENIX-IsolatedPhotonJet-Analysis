#pragma once
/* Shared, unchanged native match validation and analysis-view matching. */
#include "PhotonJet.h"
namespace photonjet { namespace response {
/* Reco-truth matches of one event. Links index the event's concatenated collections in meta/collections order: reco
   jet02, jet03, jet04 (then jetunsub02..04 in Au+Au) and truth truthjet02, truthjet03, truthjet04. */

/* Photons: the final matches (link class kLinkMatch). A reco photon with a second final match throws, as the
   reference does. */

/* Jets: the producer makes its final jet match in one greedy pass over the candidate edges (kLinkMatchCandidate) of
   every reconstructed jet view at once, so in Au+Au an unsubtracted jet can take the truth jet of its analysis twin;
   the analysis jet would then count as a fake and its truth jet as a miss. The analysis match is the same greedy
   (greedyJetMatch) rerun over the analysis view alone, for R = 0.4 jet04 x truthjet04. It is used only when the greedy
   over every view of the analysis radius (jet04 and jetunsub04) reproduces the producer's final matches of that
   radius exactly, which proves the
   ordering, the stored precision, and that the file keeps every candidate edge. In pp (one view) the analysis match
   is the producer's. Throws on a jet link joining different radii (the producer never makes one), on a candidate
   edge without a finite Delta R or reco pT, and on a reco jet with a second final match, as the reference does. */
struct EventMatches {
  std::vector<long> photonTruth;     /* truth photon of each reco photon, -1 when none */
  std::vector<double> photonMetric;  /* its match metric */
  std::vector<long> jetTruth;        /* truth jet of each analysis jet (same radius), -1 when none */
};

/* Where the collections of the analysis radius sit among the event's concatenated reco and truth jet collections:
   the subtracted reco jets (jet), the unsubtracted ones (unsub) and the truth jets. */
struct JetRanges {
  size_t nReco, jetStart, nJet, unsubStart, nUnsub, nTruth, truthStart, nTruthJet;
  bool jet(size_t i) const { return i >= jetStart && i < jetStart + nJet; }
  bool unsub(size_t i) const { return i >= unsubStart && i < unsubStart + nUnsub; }
  bool truth(size_t i) const { return i >= truthStart && i < truthStart + nTruthJet; }
  /* From the sizes of the three radii (0.2, 0.3, 0.4) of each kind and the index of the analysis radius. */
  static JetRanges of(const std::array<size_t, 3>& jets, const std::array<size_t, 3>& unsubs,
                      const std::array<size_t, 3>& truths, int radius) {
    auto before = [radius](const std::array<size_t, 3>& n) { size_t sum = 0; for (int i = 0; i < radius; ++i) sum += n[i]; return sum; };
    const size_t allJets = jets[0] + jets[1] + jets[2];
    return {allJets + unsubs[0] + unsubs[1] + unsubs[2], before(jets), jets[radius], allJets + before(unsubs), unsubs[radius],
            truths[0] + truths[1] + truths[2], before(truths), truths[radius]};
  }
};

template <class Ints, class Values>
EventMatches eventMatches(size_t nPhoton, size_t nTruthPhoton, const Ints& pReco, const Ints& pTruth,
                          const Ints& pClass, const Values& pMetric, size_t nPhotonLinks, const JetRanges& r,
                          const Values& jetRawPt, const Values& unsubRawPt, const Ints& jReco, const Ints& jTruth,
                          const Ints& jClass, const Values& jMetric, size_t nJetLinks) {
  EventMatches m{std::vector<long>(nPhoton, -1),
                 std::vector<double>(nPhoton, std::numeric_limits<double>::quiet_NaN()),
                 std::vector<long>(r.nJet, -1)};
  for (size_t l = 0; l < nPhotonLinks; ++l) {
    if (pClass[l] != kLinkMatch || pReco[l] < 0 || pTruth[l] < 0) continue;
    const size_t reco = pReco[l], truth = pTruth[l];
    if (reco >= nPhoton || truth >= nTruthPhoton) throw std::runtime_error("photon link outside its collections");
    if (m.photonTruth[reco] >= 0) throw std::runtime_error("reco photon with two final matches");
    m.photonTruth[reco] = truth;
    m.photonMetric[reco] = pMetric[l];
  }
  std::vector<char> matched(r.nReco, 0);
  std::vector<std::pair<size_t, size_t>> produced;
  std::vector<JetEdge> everyView, analysisView;
  for (size_t l = 0; l < nJetLinks; ++l) {
    const int linkClass = jClass[l];
    if (linkClass != kLinkMatch && linkClass != kLinkMatchCandidate) continue;
    if (jReco[l] < 0 || jTruth[l] < 0) {
      if (linkClass == kLinkMatchCandidate) throw std::runtime_error("jet candidate edge without both ends");
      continue;
    }
    const size_t reco = jReco[l], truth = jTruth[l];
    if (reco >= r.nReco || truth >= r.nTruth) throw std::runtime_error("jet link outside its collections");
    const bool analysisRadius = r.jet(reco) || r.unsub(reco);
    if (analysisRadius != r.truth(truth)) throw std::runtime_error("jet link joining different radii");
    if (linkClass == kLinkMatch) {
      if (matched[reco]) throw std::runtime_error("reco jet with two final matches");
      matched[reco] = 1;
      if (analysisRadius) produced.emplace_back(reco, truth);
      continue;
    }
    if (!analysisRadius) continue;
    const double pt = r.jet(reco) ? jetRawPt[reco - r.jetStart] : unsubRawPt[reco - r.unsubStart];
    if (!std::isfinite(jMetric[l]) || !std::isfinite(pt))
      throw std::runtime_error("jet candidate edge without a finite Delta R or reco pT");
    everyView.push_back({reco, truth, jMetric[l], pt});
    if (r.jet(reco)) analysisView.push_back(everyView.back());
  }
  std::sort(produced.begin(), produced.end());
  if (greedyJetMatch(everyView) != produced)
    throw std::runtime_error("the jet candidate edges do not reproduce the producer's final matches");
  for (const auto& pair : greedyJetMatch(analysisView)) m.jetTruth[pair.first - r.jetStart] = pair.second - r.truthStart;
  return m;
}

} }
