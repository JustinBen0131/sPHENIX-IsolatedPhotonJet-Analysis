// Fill the photon + recoil-jet response inputs from eventTree v1 simulation files.
//
//   root -l -q 'TreeToHists/src/FillPhotonJetResponse.C("a.eventtree.root b.eventtree.root", "out.root", "nominal")'
//
// Analysis events (evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && evt_terminal_status == 0) whose vertex
// is finite and inside the profile window. The reco photon is the region-A event-leading photon (PhotonJetSelection.h),
// the truth photon and the truth recoil jets come from PhotonJetTruth.h, and the two are paired by the reco-truth
// matches (plink_, jlink_; see eventMatches: jets use the match inside the analysis jet view, rebuilt from the
// producer's candidate edges). Per centrality class <c> (Au+Au: cent0_20, cent20_50, cent50_80; pp: all) and jet-pT
// threshold <t> (jet5, jet7, jet10, jet12 GeV), every name below carries the suffix _<c>_<t>:
//   photon_reco                             TH1D  reco pT                  region-A photon
//   photon_truth                            TH1D  truth pT                 truth photon
//   photon_response_truth_x_reco            TH2D  truth pT x reco pT       matched photons
//   photon_fakes, photon_misses             TH1D  reco pT, truth pT        unmatched reco, truth photons
//   xj_reco, xj_truth                       TH2D  pT x xJ                  reco, truth recoil jets
//   xj_response_truth_global_x_reco_global  TH2D  truth cell x reco cell   matched recoil jets
//   xj_fakes                                TH2D  reco pT x xJ             reco recoil jets without a match, split into
//   xj_fake_photon_reco, xj_combinatoric_reco, xj_detector_fakes_reco      (RecoilClass in PhotonJetTruth.h)
//   xj_misses                               TH2D  truth pT x xJ            truth recoil jets without a match
//   signal_leakage_truth_matched_abcd       TH2D  region x reco pT         one region leader per signal truth photon
//   dphi_truthjet_<side>                    TH3D  reco pT x xJ x |dphi|    every jet of a matched region-A photon passing
//   dphi_combinatoric_<side>                                               the recoil cuts except dphi, with and without
//                                                                          a truth jet nearby (nearbyTruthJet); side near
//                                                                          (|deta| <= 0.8) or far
// The dphi fills are the simulation closure of the data-driven combinatoric estimate (FinalAnalysis combinatoric.py):
// their sum is what FillPhotonJet.C fills on data, and the recoil-window integral of dphi_combinatoric is
// xj_combinatoric_reco.
// Axes (arguments; the defaults are the reference analysis's): reco pT 10 15 20 25 35 40 GeV; truth pT 5 10 15 20 25
// 35 40 GeV; xJ in 0.05 steps to 2, then 3; dphi deltaPhiEdges() (EventTreeInput.h). A response cell is pT bin *
// (xJ bins) + xJ bin (0-based). A value outside its axis is never filled (no under- or overflow) and decides the
// category as described in PhotonJetTruth.h, so only the axis ranges change the categories. The axes used (dphi
// excepted) are written as the JSON object "grids". Weights are w_event
// times, in Au+Au, the centrality weight of the event's 5% slice (an input, CentralityWeights in EventTreeInput.h; 1
// when the product is not finite), written as "centrality_weights". All histograms keep Sumw2. "sources" lists every
// input with its production input identity (meta/source input_file_sha256), which the half closure splits on.
// "variation" records the reco energy scale and resolution variation of a systematic run (Variation; "" is nominal).
#include "SelectionInput.h"
#include "PhotonJetTruth.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libTreePlayer)
#endif

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TObjString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace photonjet {
namespace response {

// The response axes. The response categories depend only on the outer range of each axis (a value outside it is
// never filled and decides the category as described in PhotonJetTruth.h), not on the inner edges. So a response
// filled on fine axes sums exactly onto any coarser axes with the same ranges whose edges lie on the fine ones.
// The defaults are the reference analysis's axes.
struct ResponseGrids {
  std::vector<double> reco{10.0, 15.0, 20.0, 25.0, 35.0, 40.0};
  std::vector<double> truth{5.0, 10.0, 15.0, 20.0, 25.0, 35.0, 40.0};
  std::vector<double> xj = xjEdges();
};

// Whitespace-separated, strictly increasing edges; an empty string keeps the default.
inline std::vector<double> parseEdges(const char* text, const std::vector<double>& fallback) {
  std::vector<double> edges;
  std::istringstream stream(text ? text : "");
  for (double value; stream >> value;) edges.push_back(value);
  if (edges.empty()) return fallback;
  if (edges.size() < 2) throw std::runtime_error("an axis needs at least two edges");
  for (size_t i = 1; i < edges.size(); ++i)
    if (!(edges[i] > edges[i - 1])) throw std::runtime_error("axis edges must be strictly increasing");
  return edges;
}

inline std::string jsonArray(const std::vector<double>& edges) {
  std::ostringstream s;
  s.precision(17);
  s << '[';
  for (size_t i = 0; i < edges.size(); ++i) s << (i ? "," : "") << edges[i];
  return s.str() + ']';
}

// A variation of the reco energy scale and resolution, for the systematic uncertainties: every reco jet (jet04) and
// photon of the event is taken at value x scale x (1 + smear g) before any selection, g a standard Gaussian drawn per
// object from a counter-based generator keyed by (seed, the file's production input, entry, object), so the draws of
// one file do not depend on the other inputs or their order. The reco-truth jet match stays the producer's, and the
// photon classes and isolation thresholds are the stored ones. The text is "" (nominal) or whitespace-separated
// key=value pairs of jet_scale, jet_smear, photon_scale, photon_smear, seed.
struct Variation {
  double jetScale = 1.0, jetSmear = 0.0, photonScale = 1.0, photonSmear = 0.0;
  std::uint64_t seed = 0;
};

inline Variation parseVariation(const char* text) {
  Variation v;
  std::vector<std::string> seen;
  std::istringstream stream(text ? text : "");
  for (std::string item; stream >> item;) {
    const size_t at = item.find('=');
    if (at == std::string::npos) throw std::runtime_error("a variation item is key=value: " + item);
    const std::string key = item.substr(0, at), value = item.substr(at + 1);
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) throw std::runtime_error("a variation key appears twice: " + key);
    seen.push_back(key);
    size_t used = 0;
    if (key == "seed") {
      if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("a seed is a non-negative integer: " + item);
      v.seed = std::stoull(value);
      continue;
    }
    double x = 0.0;
    try {
      x = std::stod(value, &used);
    } catch (const std::exception&) {
      used = 0;
    }
    if (used == 0 || used != value.size()) throw std::runtime_error("a variation value is one number: " + item);
    if (!std::isfinite(x)) throw std::runtime_error("a variation value must be finite: " + item);
    if (key == "jet_scale" || key == "photon_scale") {
      if (!(x > 0.0)) throw std::runtime_error("a scale must be positive: " + item);
      (key == "jet_scale" ? v.jetScale : v.photonScale) = x;
    } else if (key == "jet_smear" || key == "photon_smear") {
      if (x < 0.0) throw std::runtime_error("a smear must not be negative: " + item);
      (key == "jet_smear" ? v.jetSmear : v.photonSmear) = x;
    } else {
      throw std::runtime_error("unknown variation key " + key);
    }
  }
  return v;
}

inline std::string variationJson(const Variation& v) {
  std::ostringstream s;
  s.precision(17);
  s << "{\"jet_scale\":" << v.jetScale << ",\"jet_smear\":" << v.jetSmear << ",\"photon_scale\":" << v.photonScale
    << ",\"photon_smear\":" << v.photonSmear << ",\"seed\":" << v.seed << '}';
  return s.str();
}

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

inline std::uint64_t textKey(const std::string& text) {  // FNV-1a
  std::uint64_t h = 14695981039346656037ULL;
  for (unsigned char c : text) h = (h ^ c) * 1099511628211ULL;
  return h;
}

// Standard Gaussian of one key: Box-Muller on two 53-bit uniforms.
inline double gaussian(std::uint64_t key) {
  const std::uint64_t a = splitmix64(key), b = splitmix64(a);
  const double unit = 1.0 / 9007199254740992.0;  // 2^-53
  const double u1 = ((a >> 11) + 1.0) * unit, u2 = (b >> 11) * unit;
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

// The reco values of one collection of one entry under a variation; scale 1 without smearing keeps them exactly.
class Varied {
 public:
  template <class Values>
  Varied(const Values& values, double scale, double smear, std::uint64_t seed, std::uint64_t source, std::uint64_t entry,
         std::uint64_t collection)
      : v_(values.GetSize()) {
    for (size_t i = 0; i < v_.size(); ++i) {
      double factor = scale;
      if (smear != 0.0) factor *= 1.0 + smear * gaussian(splitmix64(splitmix64(splitmix64(seed ^ source) ^ entry) ^ collection) ^ i);
      v_[i] = factor == 1.0 ? values[i] : values[i] * factor;
    }
  }
  size_t GetSize() const { return v_.size(); }
  double operator[](size_t i) const { return v_[i]; }

 private:
  std::vector<double> v_;
};

// Bin of value on edges, -1 when it is not finite or outside [front, back).
inline int binIndex(double value, const std::vector<double>& edges) {
  if (!std::isfinite(value) || value < edges.front() || value >= edges.back()) return -1;
  for (size_t i = 0; i + 1 < edges.size(); ++i)
    if (edges[i] <= value && value < edges[i + 1]) return static_cast<int>(i);
  return -1;
}

// The JSON object of meta/collections whose name is `name` (compact JSON, Format tables.dumps).
inline std::string collectionEntry(const std::string& collections, const std::string& name) {
  const size_t at = collections.find("\"name\":\"" + name + "\"");
  if (at == std::string::npos) return "";
  const size_t open = collections.rfind('{', at), close = collections.find('}', at);
  return open == std::string::npos || close == std::string::npos ? "" : collections.substr(open, close - open + 1);
}

// System of a simulation file. Fails unless it is simulation with anti-kT jet04 and truthjet04 collections, the jet04
// view being the analysis view of the system.
inline System simulationSystem(TFile& f) {
  const System system = fileSystem(f);
  if (!isSimulation(f)) throw std::runtime_error(std::string("not a simulation file: ") + f.GetName());
  const std::string jet04 = collectionEntry(metaText(f, "meta/collections"), "jet04");
  const std::string view = system == System::pp ? "towerinfo_calib" : "towerinfo_sub1";
  if (jet04.find("\"algorithm\":\"antikt\"") == std::string::npos ||
      jet04.find("\"input_identity\":\"" + view + "\"") == std::string::npos)
    throw std::runtime_error(std::string("jet04 is not the anti-kT analysis view in ") + f.GetName());
  if (metaText(f, "meta/constants").find("\"truthJets.algorithm\":\"antikt\"") == std::string::npos)
    throw std::runtime_error(std::string("truth jets are not anti-kT in ") + f.GetName());
  return system;
}

// The accumulators of one centrality class and jet threshold.
struct Accumulators {
  TH1D *photonReco, *photonTruth, *photonFakes, *photonMisses;
  TH2D *photonResponse, *leakage;
  TH2D *xjReco, *xjTruth, *xjResponse, *xjFakes, *xjFakePhoton, *xjDetectorFakes, *xjCombinatoric, *xjMisses;
  std::array<TH3D*, 2> dphiTruthJet{}, dphiCombinatoric{};  // side 0 near, 1 far
  std::vector<std::unique_ptr<TH1>> owned;  // booking order, which is the write order

  Accumulators(const std::string& suffix, const ResponseGrids& g) {
    const std::vector<double>& xj = g.xj;
    const int nXj = xj.size() - 1, nReco = g.reco.size() - 1, nTruth = g.truth.size() - 1;
    const double* reco = g.reco.data();
    const double* truth = g.truth.data();
    const std::string recoAxes = ";reco p_{T}^{#gamma} [GeV];reco x_{J#gamma}";
    const std::string truthAxes = ";truth p_{T}^{#gamma} [GeV];truth x_{J#gamma}";
    auto h1 = [&](const std::string& key, const std::string& title, const std::vector<double>& edges) {
      return keep(new TH1D((key + suffix).c_str(), title.c_str(), edges.size() - 1, edges.data()));
    };
    auto h2 = [&](const std::string& key, const std::string& title, int nx, const double* x, int ny, const double* y) {
      return keep(new TH2D((key + suffix).c_str(), title.c_str(), nx, x, ny, y));
    };
    photonReco = h1("photon_reco", "region-A photons;reco p_{T}^{#gamma} [GeV]", g.reco);
    photonTruth = h1("photon_truth", "truth photons;truth p_{T}^{#gamma} [GeV]", g.truth);
    photonResponse = h2("photon_response_truth_x_reco",
                        "matched photons;truth p_{T}^{#gamma} [GeV];reco p_{T}^{#gamma} [GeV]", nTruth, truth, nReco, reco);
    photonFakes = h1("photon_fakes", "unmatched region-A photons;reco p_{T}^{#gamma} [GeV]", g.reco);
    photonMisses = h1("photon_misses", "unmatched truth photons;truth p_{T}^{#gamma} [GeV]", g.truth);
    xjReco = h2("xj_reco", "reco recoil jets" + recoAxes, nReco, reco, nXj, xj.data());
    xjTruth = h2("xj_truth", "truth recoil jets" + truthAxes, nTruth, truth, nXj, xj.data());
    std::vector<double> truthCells(nTruth * nXj + 1), recoCells(nReco * nXj + 1);
    for (size_t i = 0; i < truthCells.size(); ++i) truthCells[i] = i;
    for (size_t i = 0; i < recoCells.size(); ++i) recoCells[i] = i;
    xjResponse = h2("xj_response_truth_global_x_reco_global",
                    "matched recoil jets;truth cell (p_{T} bin #times 41 + x_{J} bin);"
                    "reco cell (p_{T} bin #times 41 + x_{J} bin)",
                    nTruth * nXj, truthCells.data(), nReco * nXj, recoCells.data());
    xjFakes = h2("xj_fakes", "unmatched reco recoil jets" + recoAxes, nReco, reco, nXj, xj.data());
    xjFakePhoton = h2("xj_fake_photon_reco", "reco recoil jets of an unmatched photon" + recoAxes, nReco, reco, nXj,
                      xj.data());
    xjDetectorFakes = h2("xj_detector_fakes_reco", "unmatched reco recoil jets with a truth jet nearby" + recoAxes,
                         nReco, reco, nXj, xj.data());
    xjMisses = h2("xj_misses", "unmatched truth recoil jets" + truthAxes, nTruth, truth, nXj, xj.data());
    xjCombinatoric = h2("xj_combinatoric_reco", "reco recoil jets without a truth jet nearby" + recoAxes, nReco, reco,
                        nXj, xj.data());
    const double regions[5] = {0, 1, 2, 3, 4};
    leakage = h2("signal_leakage_truth_matched_abcd",
                 "one region leader per signal truth photon;region (A,B,C,D);reco p_{T}^{#gamma} [GeV]", 4, regions,
                 nReco, reco);
    const std::vector<double> dphi = deltaPhiEdges();
    for (int side = 0; side < 2; ++side) {
      auto h3 = [&](const std::string& key, const std::string& title) {
        return keep(new TH3D((key + (side ? "_far" : "_near") + suffix).c_str(),
                             (title + ";reco p_{T}^{#gamma} [GeV];reco x_{J#gamma};|#Delta#phi|").c_str(), nReco, reco,
                             nXj, xj.data(), dphi.size() - 1, dphi.data()));
      };
      dphiTruthJet[side] = h3("dphi_truthjet", "jets of a matched photon with a truth jet nearby");
      dphiCombinatoric[side] = h3("dphi_combinatoric", "jets of a matched photon without a truth jet nearby");
    }
  }
  void write() {
    for (auto& h : owned) h->Write();
  }

 private:
  template <class H>
  H* keep(H* h) {
    h->Sumw2();
    h->SetDirectory(nullptr);
    owned.emplace_back(h);
    return h;
  }
};

// Reco-truth matches of one event. Links index the event's concatenated collections in meta/collections order: reco
// jet02, jet03, jet04 (then jetunsub02..04 in Au+Au) and truth truthjet02, truthjet03, truthjet04.
//
// Photons: the final matches (link class kLinkMatch). A reco photon with a second final match throws, as the
// reference does.
//
// Jets: the producer makes its final jet match in one greedy pass over the candidate edges (kLinkMatchCandidate) of
// every reconstructed jet view at once, so in Au+Au an unsubtracted jet can take the truth jet of its analysis twin;
// the analysis jet would then count as a fake and its truth jet as a miss. The analysis match is the same greedy
// (greedyJetMatch) rerun over the analysis view alone, jet04 x truthjet04. It is used only when the greedy over every
// view at R = 0.4 (jet04 and jetunsub04) reproduces the producer's final R = 0.4 matches exactly, which proves the
// ordering, the stored precision, and that the file keeps every candidate edge. In pp (one view) the analysis match
// is the producer's. Throws on a jet link joining different radii (the producer never makes one), on a candidate
// edge without a finite Delta R or reco pT, and on a reco jet with a second final match, as the reference does.
struct EventMatches {
  std::vector<long> photonTruth;     // truth photon of each reco photon, -1 when none
  std::vector<double> photonMetric;  // its match metric
  std::vector<long> jetTruth;        // truthjet04 jet of each jet04 jet, -1 when none
};

// Where the R = 0.4 collections sit among the event's concatenated reco and truth jet collections.
struct JetRanges {
  size_t nReco, jet04Start, nJet04, unsub04Start, nUnsub04, nTruth, truth04Start, nTruth04;
  bool jet04(size_t i) const { return i >= jet04Start && i < jet04Start + nJet04; }
  bool unsub04(size_t i) const { return i >= unsub04Start && i < unsub04Start + nUnsub04; }
  bool truth04(size_t i) const { return i >= truth04Start && i < truth04Start + nTruth04; }
};

template <class Ints, class Values>
EventMatches eventMatches(size_t nPhoton, size_t nTruthPhoton, const Ints& pReco, const Ints& pTruth,
                          const Ints& pClass, const Values& pMetric, size_t nPhotonLinks, const JetRanges& r,
                          const Values& jet04Pt, const Values& unsub04Pt, const Ints& jReco, const Ints& jTruth,
                          const Ints& jClass, const Values& jMetric, size_t nJetLinks) {
  EventMatches m{std::vector<long>(nPhoton, -1),
                 std::vector<double>(nPhoton, std::numeric_limits<double>::quiet_NaN()),
                 std::vector<long>(r.nJet04, -1)};
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
    const bool radius04 = r.jet04(reco) || r.unsub04(reco);
    if (radius04 != r.truth04(truth)) throw std::runtime_error("jet link joining different radii");
    if (linkClass == kLinkMatch) {
      if (matched[reco]) throw std::runtime_error("reco jet with two final matches");
      matched[reco] = 1;
      if (radius04) produced.emplace_back(reco, truth);
      continue;
    }
    if (!radius04) continue;
    const double pt = r.jet04(reco) ? jet04Pt[reco - r.jet04Start] : unsub04Pt[reco - r.unsub04Start];
    if (!std::isfinite(jMetric[l]) || !std::isfinite(pt))
      throw std::runtime_error("jet candidate edge without a finite Delta R or reco pT");
    everyView.push_back({reco, truth, jMetric[l], pt});
    if (r.jet04(reco)) analysisView.push_back(everyView.back());
  }
  std::sort(produced.begin(), produced.end());
  if (greedyJetMatch(everyView) != produced)
    throw std::runtime_error("the jet candidate edges do not reproduce the producer's final matches");
  for (const auto& pair : greedyJetMatch(analysisView)) m.jetTruth[pair.first - r.jet04Start] = pair.second - r.truth04Start;
  return m;
}

inline void fillFile(const std::string& path, const Profile& requested, const std::vector<CentralityClass>& classes,
                     const ResponseGrids& grids, const CentralityWeights& centralityWeight, const Variation& variation,
                     std::vector<std::vector<Accumulators>>& acc) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  const Profile profile = validateSelectionInput(*file, requested);
  if (simulationSystem(*file) != profile.system) throw std::runtime_error("profile system differs from " + path);
  const std::uint64_t source = textKey(sourceSha(*file));
  TTreeReader reader("eventTree", file.get());
  const VertexEventWeight vertexWeight(*file,reader);
  TTreeReaderValue<int> record(reader, "evt_record"), ready(reader, "evt_nodes_ready"), vtxValid(reader, "vtx_valid"),
      terminal(reader, "evt_terminal_status");
  TTreeReaderValue<double> vz(reader, "vtx_z"), centrality(reader, "cent_value"), weightValue(reader, "w_event");
  TTreeReaderArray<double> storedEt(reader, "pho_et"), eta(reader, "pho_eta"), phi(reader, "pho_phi"),
      score(reader, profile.scoreBranch.c_str());
  TTreeReaderArray<int> presel(reader, "pho_presel_active"), state(reader, profile.stateBranch.c_str());
  TTreeReaderArray<double> cone(reader, "isow_cone_sum"), threshold(reader, "isow_threshold"),
      sideband(reader, "isow_sideband_threshold");
  TTreeReaderArray<int> method(reader, "isow_subtraction_method"), role(reader, "isow_reconstructed_or_truth");
  TTreeReaderArray<int> idClass(reader, "pho_id_class"), isoClass(reader, "pho_iso_class");
  TTreeReaderArray<double> storedJetPt(reader, "jet04_pt"), jetEta(reader, "jet04_eta"), jetPhi(reader, "jet04_phi");
  TTreeReaderValue<int> nJet02(reader, "n_jet02"), nJet03(reader, "n_jet03"), nJet04(reader, "n_jet04"),
      nUnsub02(reader, "n_jetunsub02"), nUnsub03(reader, "n_jetunsub03"), nUnsub04(reader, "n_jetunsub04");
  TTreeReaderArray<double> truthPt(reader, "truthpho_pt"), truthEta(reader, "truthpho_eta"),
      truthPhi(reader, "truthpho_phi");
  TTreeReaderArray<int> truthRole(reader, "truthpho_source_role"), truthG4Valid(reader, "truthpho_g4_photon_valid"),
      truthPromptClass(reader, "truthpho_prompt_class");
  TTreeReaderArray<double> truthIsolation(reader, "truthpho_truth_isolation_r03");
  TTreeReaderArray<double> truthJetPt(reader, "truthjet04_pt"), truthJetEta(reader, "truthjet04_eta"),
      truthJetPhi(reader, "truthjet04_phi");
  TTreeReaderValue<int> nTruthJet02(reader, "n_truthjet02"), nTruthJet03(reader, "n_truthjet03"),
      nTruthJet04(reader, "n_truthjet04");
  TTreeReaderArray<int> plinkReco(reader, "plink_reco"), plinkTruth(reader, "plink_truth"),
      plinkClass(reader, "plink_class");
  TTreeReaderArray<double> plinkMetric(reader, "plink_metric");
  TTreeReaderArray<int> jlinkReco(reader, "jlink_reco"), jlinkTruth(reader, "jlink_truth"),
      jlinkClass(reader, "jlink_class");
  TTreeReaderArray<double> jlinkMetric(reader, "jlink_metric"), unsubPt(reader, "jetunsub04_pt");

  const std::vector<double>& xj = grids.xj;
  const long nXj = xj.size() - 1;
  // Sums: entries are in source order, the reference's event order. Every fill of one event carries the same weight,
  // so the order of the fills inside an event (region, truth photon, recoil order) cannot change any sum.
  while (reader.Next()) {
    reportProgress(reader.GetCurrentEntry()+1);
    validateSelectionEvent(profile, *centrality, storedEt, score, state, cone, threshold, sideband, method, role, idClass, isoClass);
    if (*record != 0) continue;  // upstream-rejected entries carry only evt_ fields
    // the reco energies every selection below uses (Variation); the reco-truth match uses the stored ones
    const std::uint64_t entry = static_cast<std::uint64_t>(reader.GetCurrentEntry());
    const Varied et(storedEt, variation.photonScale, variation.photonSmear, variation.seed, source, entry, 0);
    const Varied jetPt(storedJetPt, variation.jetScale, variation.jetSmear, variation.seed, source, entry, 1);
    const size_t nPhoton = et.GetSize(), nTruthJets04 = truthJetPt.GetSize();
    if (score.GetSize() != nPhoton || state.GetSize() != nPhoton)
      throw std::runtime_error("photon-ID branches pho_bdt_score, pho_bdt_state are not filled in " + path);
    if (profile.augmentedClasses && (idClass.GetSize() != nPhoton || isoClass.GetSize() != nPhoton))
      throw std::runtime_error("the augmented classes pho_id_class, pho_iso_class are not filled in " + path);
    // Matches are read for every event: the reference refuses a duplicate final match anywhere in its input.
    const JetRanges ranges{static_cast<size_t>(*nJet02 + *nJet03 + *nJet04 + *nUnsub02 + *nUnsub03 + *nUnsub04),
                           static_cast<size_t>(*nJet02 + *nJet03), storedJetPt.GetSize(),
                           static_cast<size_t>(*nJet02 + *nJet03 + *nJet04 + *nUnsub02 + *nUnsub03), unsubPt.GetSize(),
                           static_cast<size_t>(*nTruthJet02 + *nTruthJet03 + *nTruthJet04),
                           static_cast<size_t>(*nTruthJet02 + *nTruthJet03), nTruthJets04};
    if (ranges.nJet04 != static_cast<size_t>(*nJet04) || ranges.nUnsub04 != static_cast<size_t>(*nUnsub04) ||
        ranges.nTruth04 != static_cast<size_t>(*nTruthJet04))
      throw std::runtime_error("jet collection counters differ from their arrays in " + path);
    const EventMatches match =
        eventMatches(nPhoton, truthPt.GetSize(), plinkReco, plinkTruth, plinkClass, plinkMetric, plinkReco.GetSize(),
                     ranges, storedJetPt, unsubPt, jlinkReco, jlinkTruth, jlinkClass, jlinkMetric, jlinkReco.GetSize());
    // The producer's signal flag must be the analysis truth signal of the stored primitives.
    for (size_t t = 0; t < truthPt.GetSize(); ++t)
      if (analysisTruthSignal(truthG4Valid[t], truthPromptClass[t], truthEta[t], truthIsolation[t]) !=
          (truthRole[t] == kSignalSourceRole))
        throw std::runtime_error("truth signal flag differs from its definition in " + path);
    if (*ready != 1 || *vtxValid != 1 || *terminal != 0) continue;
    if (!responseEvent(profile, *vz)) continue;
    // w_event x the centrality weight; a non-finite product counts 1, as the reference does
    const double product = *weightValue * vertexWeight() * centralityWeight(*centrality);
    if (strictReduction && !std::isfinite(product)) throw std::runtime_error("nonfinite response weight");
    const double w = std::isfinite(product) ? product : 1.0;
    if (cone.GetSize() != 2 * nPhoton) throw std::runtime_error("isolation witnesses are not 2 per photon in " + path);

    const long tp = truthPhoton(truthPt, truthEta, truthPhi, truthRole, truthPt.GetSize());
    const int truthBin = tp >= 0 ? binIndex(truthPt[tp], grids.truth) : -1;
    for (size_t k = 0; k < classes.size(); ++k) {
      if (profile.system == System::auau && !(classes[k].low <= *centrality && *centrality < classes[k].high)) continue;
      const std::array<long, 4> leader = regionLeaders(profile, *vz, *centrality, nPhoton, et, eta, phi, score, presel,
                                                       state, cone, threshold, sideband, method, role, idClass, isoClass);
      const long a = leader[kA];
      const int recoBin = a >= 0 ? binIndex(et[a], grids.reco) : -1;
      const bool photonMatched = a >= 0 && tp >= 0 && match.photonTruth[a] == tp;
      const PhotonResponse photon = photonResponse(photonMatched, recoBin >= 0, truthBin >= 0);
      const std::array<bool, 4> witness =
          leakageWitnesses(leader, match.photonTruth, match.photonMetric, et, truthRole);

      for (size_t t = 0; t < kJetThresholds.size(); ++t) {
        Accumulators& h = acc[k][t];
        Profile p = profile;
        p.jetPtMin = kJetThresholds[t];
        if (truthBin >= 0) h.photonTruth->Fill(truthPt[tp], w);
        if (recoBin >= 0) h.photonReco->Fill(et[a], w);
        if (photon.matched) h.photonResponse->Fill(truthPt[tp], et[a], w);
        if (photon.fake) h.photonFakes->Fill(et[a], w);
        if (photon.miss) h.photonMisses->Fill(truthPt[tp], w);
        for (int r = 0; r < 4; ++r)
          if (witness[r] && binIndex(et[leader[r]], grids.reco) >= 0) h.leakage->Fill(r + 0.5, et[leader[r]], w);

        // Truth recoils of the truth photon inside the xJ axis, in orderedByPt order.
        std::vector<int> truthXjBin(nTruthJets04, -1);  // xJ bin of each selected truth recoil
        std::vector<long> selected;
        if (truthBin >= 0) {
          auto truthRecoil = [&](size_t j) {
            return truthRecoilJet(p, truthEta[tp], truthPhi[tp], truthJetPt[j], truthJetEta[j], truthJetPhi[j]);
          };
          for (long j : orderedByPt(nTruthJets04, truthRecoil, truthJetPt)) {
            const double x = truthJetPt[j] / truthPt[tp];
            const int b = binIndex(x, xj);
            if (b < 0) continue;
            truthXjBin[j] = b;
            selected.push_back(j);
            h.xjTruth->Fill(truthPt[tp], x, w);
          }
        }

        // Reco recoils of the region-A photon, in orderedByPt order. The order decides which of several reco jets
        // linked to one truth recoil is matched: the first; the others are fakes (RecoilClass).
        std::vector<char> taken(nTruthJets04, 0);
        if (a >= 0) {
          auto recoRecoil = [&](size_t j) { return recoilJet(p, eta[a], phi[a], jetPt[j], jetEta[j], jetPhi[j]); };
          for (long j : orderedByPt(jetPt.GetSize(), recoRecoil, jetPt)) {
            const double x = jetPt[j] / et[a];
            const int b = binIndex(x, xj);
            if (recoBin < 0 || b < 0) continue;
            h.xjReco->Fill(et[a], x, w);
            const long tj = match.jetTruth[j];
            const bool available = tj >= 0 && truthXjBin[tj] >= 0 && !taken[tj];
            const RecoilClass category = recoilClass(photonMatched, available, [&] {
              return nearbyTruthJet(p.jetPtMin, jetEta[j], jetPhi[j], truthJetPt, truthJetEta, truthJetPhi,
                                    nTruthJets04);
            });
            if (category == RecoilClass::kMatched) {
              h.xjResponse->Fill(truthBin * nXj + truthXjBin[tj] + 0.5, recoBin * nXj + b + 0.5, w);
              taken[tj] = 1;
              continue;
            }
            h.xjFakes->Fill(et[a], x, w);
            if (category == RecoilClass::kFakePhoton) h.xjFakePhoton->Fill(et[a], x, w);
            else if (category == RecoilClass::kCombinatoric) h.xjCombinatoric->Fill(et[a], x, w);
            else h.xjDetectorFakes->Fill(et[a], x, w);
          }
        }
        for (long j : selected)
          if (!taken[j]) h.xjMisses->Fill(truthPt[tp], truthJetPt[j] / truthPt[tp], w);

        // Every jet of a matched photon at any Delta phi, by whether a truth jet is nearby.
        if (photonMatched && recoBin >= 0) {
          Profile anyAngle = p;
          anyAngle.deltaPhiMin = 0.0;
          for (size_t j = 0; j < jetPt.GetSize(); ++j) {
            if (!recoilJet(anyAngle, eta[a], phi[a], jetPt[j], jetEta[j], jetPhi[j])) continue;
            const double x = jetPt[j] / et[a];
            if (binIndex(x, xj) < 0) continue;
            const int side = std::fabs(jetEta[j] - eta[a]) > kFarDeltaEta ? 1 : 0;
            const bool nearby =
                nearbyTruthJet(p.jetPtMin, jetEta[j], jetPhi[j], truthJetPt, truthJetEta, truthJetPhi, nTruthJets04);
            (nearby ? h.dphiTruthJet : h.dphiCombinatoric)[side]->Fill(et[a], x, wrappedAbsDeltaPhi(jetPhi[j], phi[a]), w);
          }
        }
      }
    }
  }
  // a missing branch or a read error ends the loop like the end of the tree; only the end is accepted
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
}

inline int run(const char* inputs, const char* output, const char* profileName, const char* recoEdges,
               const char* truthEdges, const char* xjAxisEdges, const char* centralityWeightText, const char* variationText,
               const char* selection) {
  const Variation variation = parseVariation(variationText);
  const ResponseGrids defaults;
  ResponseGrids grids;
  grids.reco = parseEdges(recoEdges, defaults.reco);
  grids.truth = parseEdges(truthEdges, defaults.truth);
  grids.xj = parseEdges(xjAxisEdges, defaults.xj);
  const std::vector<std::string> paths = inputPaths(inputs);
  System system;
  {
    std::unique_ptr<TFile> first(TFile::Open(paths.front().c_str(), "READ"));
    if (!first || first->IsZombie()) throw std::runtime_error("cannot open " + paths.front());
    system = simulationSystem(*first);
  }
  const CentralityWeights centralityWeight = centralityWeights(centralityWeightText, system);
  const Profile profile = configuredProfile(profileName, system, selection);
  const std::vector<CentralityClass> classes = centralityClasses(system);
  std::vector<std::vector<Accumulators>> acc(classes.size());
  for (size_t k = 0; k < classes.size(); ++k)
    for (double t : kJetThresholds) acc[k].emplace_back("_" + classes[k].name + "_" + thresholdName(t), grids);
  for (const auto& path : paths) fillFile(path, profile, classes, grids, centralityWeight, variation, acc);
  const std::string sources = sourceRecords(paths);  // before the output is opened (sourceRecord)
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  for (auto& perClass : acc)
    for (auto& h : perClass) h.write();
  TObjString(profile.name.c_str()).Write("profile");
  writeSelection(profile);
  TObjString(("{\"reco_pt\":" + jsonArray(grids.reco) + ",\"truth_pt\":" + jsonArray(grids.truth) + ",\"xj\":" +
              jsonArray(grids.xj) + "}").c_str()).Write("grids");
  TObjString(("{\"per_slice\":" + jsonArray(centralityWeight.perSlice) + "}").c_str()).Write("centrality_weights");
  TObjString(variationJson(variation).c_str()).Write("variation");
  TObjString(sources.c_str()).Write("sources");
  out->Close();
  return 0;
}

}  // namespace response
}  // namespace photonjet

// inputs: whitespace-separated eventTree simulation paths. profile: "nominal" (the analysis), a variant of it
// (EventTreeInput.h namedProfile), or "historical" (the reference reproduction only).
// Optional axes (whitespace-separated edges; empty keeps the reference axes): reco photon pT, truth photon pT, xJ.
// centralityWeights: Au+Au only and required there, 16 weights for the 5% slices of 0-80% (CentralityWeights).
// variation: "" (nominal) or the reco energy scale and resolution variation of the systematics (Variation).
int FillPhotonJetResponse(const char* inputs, const char* output, const char* profileName = "nominal",
                          const char* recoPtEdges = "", const char* truthPtEdges = "", const char* xjEdges = "",
                          const char* centralityWeights = "", const char* variation = "", const char* selection = "") {
  return photonjet::response::run(inputs, output, profileName, recoPtEdges, truthPtEdges, xjEdges, centralityWeights,
                                  variation, selection);
}
