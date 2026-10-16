#include "NativeModelDefinition.h"
/* TreeToHists: eventTree v1 files -> one histogram package per DATA collection or simulated sample. */

/* tree_to_hists plan       --settings CFG --inputs INPUTS.yaml --work DIR [--sample SYSTEM/NAME]...
   tree_to_hists fill       --work DIR --job N
   tree_to_hists merge      --work DIR --sample SYSTEM/NAME --out DIR
   tree_to_hists centrality --settings CFG --inputs INPUTS.yaml --data AuAu/data.root --out weights.json
   tree_to_hists isolation  population|derive|freeze|check ... (Isolation.cc) */

/* Sections: the settings and the resolved selection; the package (inventory and merge); the catalogue, the plan and
   one fill job; the Au+Au centrality weights; main. The physics lives in the fills (FillPhotonJet.C,
   FillPhotonJetResponse.C, FillPhotonQA.C) and PhotonJet.h. */
#include "TreeToHists.h"
#include "CandidateCachePackage.h"
#include "BuildInfo.h"
#include "ScientificProducerEvidence.h"
#include "Configuration.h"
#include "FillPhotonQA.C"
#include "FillPhotonJet.C"
#include "FillPhotonJetResponse.C"
#include "EmsComposition.h"
#include "NativeResponseGrid.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unistd.h>

#include "TEnv.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "THnSparse.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

namespace treetohists {

/* The settings and the resolved selection */
struct Settings {
  std::string path;                 /* the settings file */
  std::string home;                 /* the TreeToHists directory (the canonical packages are found from it) */
  Json selection;                   /* resolved selection (PhotonJetSelectionV1), written to every product and package */
  Json emsSelectionContract;
  Json responseGrid;                /* opt-in response components, outside common histogram settings */
  Json axes;                        /* native_xj_edges, native_reco_response_edges_gev, native_truth_response_edges_gev */
  bool extensions = true, requireStage3Fractions = false, candidateCache = false;
  long candidateCacheMaxRows = 2000000;
  double candidateCacheMinEt = 0.;
  std::string responseVariation;    /* "" nominal, or the reco energy scale and resolution variation of a systematic run */
  double jetRadius = 0.4;           /* Fill.JetRadius: the jet collection of the fill, 0.2, 0.3 or 0.4 (photonjet::configureJetRadius) */
  double dataJetScalePp = 1.0, dataJetScaleAuAu = 1.0;  /* Fill.DataJetScale.pp, .auau: photonjet::configuredDataJetScale */
  static Settings read(const std::string &path, const std::string &home);
};
/* Plans and timing receipts bind the same resolved histogram settings. */
static Json histogramSettings(const Settings &s) {
  Json histograms=s.axes;
  auto edges=photonjet::photonEtEdges();
  for (const auto &[system,selection] : s.selection.at("systems").items()) if (selection.contains("optional_model")) {
    const auto &model=selection.at("optional_model"); const auto &domain=model.at("domain");
    const int low=std::floor(domain.at("et_min").get<double>()),high=std::ceil(domain.at("et_max_exclusive").get<double>());
    for (int e=low;e<=high;++e) edges.push_back(e);
    if (model.contains("et_ranges_gev")) histograms["model_et_ranges_gev"][system]=model.at("et_ranges_gev");
    for (const char *axis : {"native_reco_response_edges_gev","native_truth_response_edges_gev"}) {
      auto response=histograms.at(axis).get<std::vector<double>>();
      if (model.contains("et_ranges_gev")) for (const char *range : {"application","analysis"})
        for (const auto &edge : model.at("et_ranges_gev").at(range)) response.push_back(edge.get<double>());
      std::sort(response.begin(),response.end());response.erase(std::unique(response.begin(),response.end()),response.end());histograms[axis]=response;
    }
  }
  std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
  histograms.update(Json{{"native_measurement_et_edges_gev",edges}, {"jet_thresholds_gev",photonjet::kJetThresholds},
    {"reported_photon_et_edges_gev",{17,23,33}}, {"require_stage3_fractions",s.requireStage3Fractions},
    {"candidate_records",s.candidateCache}, {"candidate_records_max_rows",s.candidateCacheMaxRows}, {"candidate_records_min_et",s.candidateCacheMinEt}, {"response_variation",s.responseVariation}, {"jet_radius",s.jetRadius},
    {"data_jet_scale",{{"pp",s.dataJetScalePp},{"auau",s.dataJetScaleAuAu}}}});
  if(!s.emsSelectionContract.is_null()){
    histograms["ems_selection_contract"]=s.emsSelectionContract;
    histograms["ems_qa_extensions"]="nominal_archival_only";
    histograms["ems_candidate_archive"]="original_nominal_inputs_replay_unqualified";
  }
  return histograms;
}
std::string sha256Text(const std::string &text);
std::string canonical(const Json &value);  /* compact, sorted keys, ASCII: the digest form */
std::string digest(const Json &value);     /* sha256 of canonical(value) */
Json readJson(const std::string &path);
/* A completed fill is reusable only with the code that produced its plan.
   Old work directories remain evidence; a changed recipe needs a new plan. */
Json reductionCode() {
  Json code = Json::object();
#ifdef TREETOHISTS_SOURCE_DIR
  const std::filesystem::path src = TREETOHISTS_SOURCE_DIR;
#else
  const std::filesystem::path src = std::filesystem::path(__FILE__).parent_path();
#endif
  for (const char *name : {"TreeToHists.cc", "TreeToHists.h", "CandidateCache.h", "CandidateCachePackage.h", "MeasurementAccumulator.h", "CandidateMeasurementReplay.h", "CandidateResponseReplay.h", "ResponseMatching.h", "Campaign.cc", "Package.cc", "Isolation.cc", "Corrections.cc", "Configuration.h", "PhotonJet.h", "Preselection.h", "MbdTimeReference.h", "LeadingAssociation.h", "FillPhotonJet.C", "FillPhotonQA.C", "FillPhotonJetResponse.C", "EnergyVariation.h", "CentralityWeightSource.h", "CentralityReplay.h", "PhotonEnergySelection.h", "PhotonEnergyInputs.h", "NativeEnergyVariation.h", "PurityInputs.h", "SelectionReplay.h", "ResponseOnlyResolution.h", "CollectionContract.h", "EmsComposition.h", "EmsNativeSelection.h", "NativeResponseGrid.h", "SelectionViewContract.h", "SelectionEfficiency.h", "SelectionFanout.cc", "ScientificTrace.h", "HistogramScientificDefinitions.h", "QAScientificContext.h", "MeasurementScientificDefinitions.h", "CandidateScientificContext.h", "ReplayScientificContext.h", "NativeModelDefinition.h", "ScoreReviewScientific.h", "ScoreReview.C", "ScientificProducerEvidence.h"})
    code[name] = sha256File((src / name).string());
  code["TreeProduction/Configuration.h"] = sha256File((src / "../../TreeProduction/src/internal/Configuration.h").lexically_normal().string());
  return code;
}
void checkReductionCode(const Json &plan) {
  if (plan.value("producer_code", Json()) != reductionCode() || plan.at("producer_code")!=Json::parse(TREETOHISTS_COMPILED_CODE))
    throw std::runtime_error("reduction code differs from the plan; prepare a new work directory");
}

/* Planned model bytes may not drift between sealing the snapshot and a
   native fill. Saved packages remain self-contained after transfer. */
void checkNativeModelFiles(const Json &plan,const std::string &system) {
  const auto &selection=plan.at("selection").at("systems").at(system=="PP" ? "pp" : "auau");
  if (!selection.contains("optional_model") || !selection.at("optional_model").contains("native_evaluation")) return;
  for (const auto &cell : selection.at("optional_model").at("native_evaluation").at("cells"))
    if (sha256File(cell.at("file"))!=cell.at("sha256").get<std::string>()) throw std::runtime_error("native model file changed after input sealing");
}
Json treeCopySummary(const Json &sample) {
  if (sample.at("system")!="AuAu" || sample.at("family")=="data") return nullptr;
  Json counts={{"matching_receipt",0},{"other_receipt_ignored",0},{"unweighted",0}};
  for (const auto &source : sample.at("sources")) {
    const Json w=source.at("metadata").value("weighting",Json::object());
    const std::string key=!photonjet::declaredAuAuCopies(w) ? "unweighted" :
      !sample.at("vertex_reweighting").is_null() && w.at("receipt_sha256")==sample.at("vertex_reweighting").at("evidence_sha256") &&
      w.at("stitching_sha256")==sample.at("normalization").at("stitching_table_sha256") ? "matching_receipt" : "other_receipt_ignored";
    counts[key]=counts.at(key).get<long long>()+1;
  }
  return Json{{"sources",counts},{"applied_in_histograms",false}};
}
/* The vertex histogram is external to the eventTrees. Recheck its planned
   hash at execution, on completed-job reuse and before package publication. */
void checkVertexFit(const Json &sample) {
  const Json fit = sample.value("vertex_reweighting", Json());
  if (fit.is_null()) return;
  const std::string path = fit.at("path");
  if (!std::filesystem::is_regular_file(path) || sha256File(path) != fit.at("sha256").get<std::string>())
    throw std::runtime_error(sample.at("key").get<std::string>() + ": vertex fit is missing or changed since planning");
  const std::string evidence=fit.at("evidence");
  if (sha256File(evidence)!=fit.at("evidence_sha256").get<std::string>())
    throw std::runtime_error("AuAu weight qualification receipt changed since planning");
  const Json receipt=readJson(evidence);
  if (receipt.value("schema", "")!="AuAuSimWeightFreezeV1" || receipt.value("status", "")!="QUALIFIED" ||
      receipt.at("selection_sha256")!=fit.at("selection_sha256") || receipt.at("data_trigger")!=fit.at("data_trigger") ||
      receipt.at("centrality").at("weights")!=sample.at("weights") ||
      receipt.at("centrality").at("weights_sha256")!=digest(sample.at("weights")) ||
      receipt.at("vertex").at("sha256")!=fit.at("sha256") || receipt.at("vertex").at("histogram")!=fit.at("histogram") ||
      receipt.at("vertex").at("coordinate")!="truth_vtx_z")
    throw std::runtime_error("AuAu simulation needs a matching QUALIFIED weight receipt; TEST_ONLY is not production input");
}
/* The p+p simulation periods (photonjet::configurePpPeriods, the PPG12 two-period sum): TreeProduction's ppg12_weights,
   the constants and truth-vertex fits the producer used for each source's own period. Mode all (the default) sums
   every period per event; mode source keeps w_event as written. */
Json ppSimulationPeriods(const std::filesystem::path &home, const std::string &mode,const std::string &fitDirectory="") {
  if (mode != "all" && mode != "source") throw std::runtime_error("pp.SimulationPeriods is all or source");
  if (mode == "source") {
    if(!fitDirectory.empty())throw std::runtime_error("pp.VertexFitDirectory requires SimulationPeriods: all");
    return Json{{"mode", mode}};
  }
  const std::filesystem::path config = (home / "../TreeProduction/config/tree_production.yaml").lexically_normal();
  const Json w = photonjet_config::readYaml(config.string()).at("ppg12_weights");
  const double total = w.at("total_luminosity_pb_inv").get<double>();
  Json periods = Json::object();
  double shares = 0;
  for (const auto &[name, p] : w.at("periods").items()) {
    const Json &fit = p.at("vertex_reweighting");
    if (!fit.value("enabled", false)) throw std::runtime_error("p+p period " + name + " has no enabled truth-vertex fit");
    const double share = p.at("luminosity_pb_inv").get<double>() / total;
    shares += share;
    periods[name] = {{"share", share}, {"double_fraction", p.at("double_fraction")}, {"luminosity_pb_inv", p.at("luminosity_pb_inv")},
                     {"runs", p.at("runs")}, {"fit", {{"path", fit.at("file")}, {"sha256", fit.at("sha256")}, {"histogram", fit.at("histogram")}}}};
    // A portable location override changes neither the fit bytes nor its expected hash.
    if(!fitDirectory.empty())periods[name]["fit"]["path"]=(std::filesystem::path(fitDirectory)/name/"reweight.root").string();
  }
  if (periods.empty() || std::abs(shares - 1.0) > 1e-9) throw std::runtime_error("the p+p period luminosities do not add up to the total");
  return Json{{"mode", mode}, {"rule", "w = w_slice x sum_p share_p x mix_p x v_p(z_truth)"}, {"total_luminosity_pb_inv", total},
              {"periods", periods}, {"source", "TreeProduction/config/tree_production.yaml ppg12_weights"}};
}
/* The period fits are external to the eventTrees: recheck their hashes at planning, at execution, on completed-job
   reuse and before package publication. */
void checkPpPeriodFits(const Json &periods) {
  if (!periods.is_object() || periods.value("mode", "") != "all") return;
  for (const auto &[name, p] : periods.at("periods").items()) {
    const std::string path = p.at("fit").at("path");
    if (!std::filesystem::is_regular_file(path) || sha256File(path) != p.at("fit").at("sha256").get<std::string>())
      throw std::runtime_error("the " + name + " p+p truth-vertex fit is missing or changed: " + path);
  }
}
/* The p+p period definition a sample's jobs use: the selection's for p+p simulation, none otherwise. */
Json ppPeriodsOf(const Json &plan, const Json &sample) {
  if (sample.at("system") != "PP" || sample.at("family") == "data") return Json();
  return plan.at("selection").at("systems").at("pp").value("simulation_periods", Json());
}
/* The package: fill-job directories, the inventory and the merge */
std::string productDirectory(const std::string &work, const std::string &sample, long source);  /* sharded below 1000 */
Json inventory(const std::string &system, const std::string &family);
void merge(const Json &plan, const Json &sample, const std::string &work, const std::string &output);
/* The canonical Au+Au embedded stitching (config/auau_embedded_stitching.json): a sample's normalisation record, and
   the events of one eventTree file in its denominator ({"kind", "events", "missing"}) */
Json stitchingNormalization(const std::string &home, const std::string &system, const std::string &name, const std::string &family);
Json stitchingDenominator(const std::string &path, const Json &normalization);
/* The Au+Au simulation centrality weights from the DATA package and the simulated distributions */
Json centralityWeights(const Settings &s, const std::string &catalogue, const std::string &dataPackage);

}  /* namespace treetohists */

/* === Settings ========================================================================================================
   Settings: the settings file and the one resolved selection of every DATA and simulation file. */

/* The selection starts from two frozen packages: the photon-ID definition (PhotonID/models/nominal.json: models,
   working points, regions) and the Au+Au isolation line (config/auau_isolation_wp80.json). A settings file may override
   a score branch, a working-point line or an isolation parameter, or select a sealed PhotonID release instead of the
   nominal score. The resolved selection (PhotonJetSelectionV1) is what the fillers read; its sha256 is the profile
   every product and package carries, so DATA and simulation of one result always share one selection. */
namespace treetohists {
namespace {

/* SHA-256 (FIPS 180-4). */
class Sha256 {
 public:
  void update(const unsigned char *p, size_t n) {
    length_ += n;
    while (n) {
      const size_t take = std::min(n, 64 - used_);
      std::memcpy(block_ + used_, p, take);
      used_ += take, p += take, n -= take;
      if (used_ == 64) compress(), used_ = 0;
    }
  }
  std::string hex() {
    const uint64_t bits = length_ * 8;
    const unsigned char one = 0x80, zero = 0;
    update(&one, 1);
    while (used_ != 56) update(&zero, 1);
    unsigned char size[8];
    for (int i = 0; i < 8; ++i) size[i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
    update(size, 8);
    std::ostringstream out;
    for (uint32_t v : h_) for (int s = 28; s >= 0; s -= 4) out << "0123456789abcdef"[(v >> s) & 15];
    return out.str();
  }

 private:
  static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
  void compress() {
    static const uint32_t k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98,
        0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
        0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8,
        0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819,
        0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
        0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2};
    uint32_t w[64];
    for (int i = 0; i < 16; ++i)
      w[i] = uint32_t(block_[4 * i]) << 24 | uint32_t(block_[4 * i + 1]) << 16 | uint32_t(block_[4 * i + 2]) << 8 | block_[4 * i + 3];
    for (int i = 16; i < 64; ++i)
      w[i] = w[i - 16] + (rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3)) + w[i - 7] +
             (rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10));
    uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3], e = h_[4], f = h_[5], g = h_[6], h = h_[7];
    for (int i = 0; i < 64; ++i) {
      const uint32_t t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      h = g, g = f, f = e, e = d + t1, d = c, c = b, b = a, a = t1 + t2;
    }
    const uint32_t add[8] = {a, b, c, d, e, f, g, h};
    for (int i = 0; i < 8; ++i) h_[i] += add[i];
  }
  uint32_t h_[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  unsigned char block_[64] = {};
  size_t used_ = 0;
  uint64_t length_ = 0;
};

/* A setting's value: the text after "Key:", up to a trailing "# comment". */
std::string value(TEnv &env, const std::string &key) {
  std::string v = env.GetValue(key.c_str(), "");
  v = v.substr(0, v.find('#'));
  const auto first = v.find_first_not_of(" \t"), last = v.find_last_not_of(" \t");
  return first == std::string::npos ? "" : v.substr(first, last - first + 1);
}

std::string required(TEnv &env, const std::string &key) {
  const std::string v = value(env, key);
  if (v.empty()) throw std::runtime_error("setting " + key + " is missing");
  return v;
}

std::vector<double> numbers(const std::string &key, const std::string &text) {
  std::istringstream in(text);
  std::vector<double> out;
  for (std::string item; in >> item;) {
    size_t used = 0;
    double x = 0;
    try { x = std::stod(item, &used); } catch (const std::exception &) { used = 0; }
    if (used != item.size() || !std::isfinite(x)) throw std::runtime_error("setting " + key + ": '" + item + "' is not a number");
    out.push_back(x);
  }
  return out;
}

Json line(const std::string &key, const std::string &text) {  /* "intercept slope" */
  const auto v = numbers(key, text);
  if (v.size() != 2) throw std::runtime_error("setting " + key + " is a line: intercept slope");
  return Json{{"intercept", v[0]}, {"slope", v[1]}};
}

Json line(const Json &v, const std::string &label) {
  if (!v.is_object() || v.size() != 2 || !v.contains("intercept") || !v.contains("slope"))
    throw std::runtime_error(label + " needs intercept and slope");
  for (const char *k : {"intercept", "slope"})
    if (!v.at(k).is_number() || !std::isfinite(v.at(k).get<double>())) throw std::runtime_error(label + " must be finite");
  return Json{{"intercept", v.at("intercept").get<double>()}, {"slope", v.at("slope").get<double>()}};
}

double finite(const Json &v, const std::string &label) {
  if (!v.is_number() || v.is_boolean() || !std::isfinite(v.get<double>())) throw std::runtime_error(label + " must be a finite number");
  return v.get<double>();
}

Json package(const std::string &path, const std::string &schema) {
  Json p = readJson(path);
  if (p.value("schema", "") != schema) throw std::runtime_error(path + ": expected " + schema);
  return p;
}

const std::regex kBranch("[A-Za-z][A-Za-z0-9_]*"), kSha("[0-9a-f]{64}");

/* A sealed PhotonID release (PhotonIDSelectionV1 + PhotonIDReleaseV1 + PhotonIDWorkingPointsV1) as the selection's
   score: its single score branch, applicability domain and routed WP70/WP80 lines. */
Json optionalSelection(const std::filesystem::path &path, const std::string &system, const Json &features,
                       const std::string &tightWP="WP70", const std::string &nonTightWP="WP80",
                       double nonTightMinimum=0.1) {
  for (const auto &name : {tightWP, nonTightWP})
    if (name!="WP70" && name!="WP80" && name!="WP90")
      throw std::runtime_error("choose a released working point: WP70, WP80 or WP90");
  if (!std::isfinite(nonTightMinimum) || nonTightMinimum<0 || nonTightMinimum>=1)
    throw std::runtime_error("NonTightMinimum must be finite in [0,1)");
  const auto folder = path.parent_path();
  const Json spec = package(path.string(), "PhotonIDSelectionV1");
  const Json release = package((folder / "manifest.json").string(), "PhotonIDReleaseV1");
  Json body = release;
  body.erase("release_sha256");
  if (release.value("status", "") != "PASS" || digest(body) != release.value("release_sha256", ""))
    throw std::runtime_error("optional model release is not validated or its manifest changed");
  for (auto it = release.at("artifacts").begin(); it != release.at("artifacts").end(); ++it)
    if (sha256File((folder / it.key()).string()) != it.value().get<std::string>())
      throw std::runtime_error("optional model artifact changed: " + it.key());
  if (!release.at("artifacts").contains(path.filename().string()) || !release.at("artifacts").contains(spec.at("working_points_file").get<std::string>()))
    throw std::runtime_error("selection and working points must be sealed release artifacts");
  const std::string systemName = system == "pp" ? "PP" : "AuAu";
  if (spec.at("system").get<std::string>() != systemName || spec.at("system") != release.at("system")) throw std::runtime_error("optional selection collision system differs");
  if (!spec.at("state_branch").is_null() || spec.at("score_branch") != release.at("score_branch"))
    throw std::runtime_error("optional selection must use the released single score");
  const std::string score = spec.at("score_branch");
  if (!std::regex_match(score, kBranch) || score == "pho_bdt_score") throw std::runtime_error("optional score collides with nominal or has an invalid name");
  if (spec.at("cells") != release.at("cells") || spec.at("domain") != release.at("recipe").at("domain"))
    throw std::runtime_error("selection routing/domain differs from its released model");
  if (spec.at("non_tight_lower") != 0.1 || spec.at("tight") != "score > WP70" || spec.at("non_tight") != "0.1 < score < WP80")
    throw std::runtime_error("unsupported optional photon-ID inequalities");
  const auto wpPath = folder / spec.at("working_points_file").get<std::string>();
  const Json wp = package(wpPath.string(), "PhotonIDWorkingPointsV1");
  const Json &d = spec.at("domain");
  const double etMin = finite(d.at("et_min"), "ET min"), etMax = finite(d.at("et_max_exclusive"), "ET max");
  const double centMin = d.contains("centrality_min") ? finite(d.at("centrality_min"), "centrality min") : 0.0;
  const double centMax = d.contains("centrality_max_exclusive") ? finite(d.at("centrality_max_exclusive"), "centrality max") : 80.0;
  if (etMin >= etMax || centMin >= centMax) throw std::runtime_error("invalid optional model domain");
  Json cells = Json::array();
  std::set<std::string> ids;
  std::vector<double> ets{etMin, etMax}, cents{centMin, centMax};
  for (const auto &cell : spec.at("cells")) {
    const Json &bounds = cell.at("bounds");
    for (auto it = bounds.begin(); it != bounds.end(); ++it)
      if (it.key() != "et" && it.key() != "centrality_coarse" && it.key() != "centrality_fine") throw std::runtime_error("model cell: unknown bound " + it.key());
    if (bounds.contains("centrality_coarse") && bounds.contains("centrality_fine")) throw std::runtime_error("cell has ambiguous centrality routing");
    const Json et = bounds.contains("et") ? bounds.at("et") : Json{etMin, etMax};
    const Json cent = bounds.contains("centrality_coarse") ? bounds.at("centrality_coarse")
                      : bounds.contains("centrality_fine") ? bounds.at("centrality_fine") : Json{centMin, centMax};
    for (const auto &[axis, lo, hi] : {std::tuple<Json, double, double>{et, etMin, etMax}, {cent, centMin, centMax}})
      if (axis.size() != 2 || !(lo <= axis[0].get<double>() && axis[0].get<double>() < axis[1].get<double>() && axis[1].get<double>() <= hi))
        throw std::runtime_error("model cell lies outside its domain");
    const std::string id = cell.at("id");
    const Json &w = wp.at("cells").at(id);
    if (w.value("status", "") != "PASS" || w.at("axis") != (system == "pp" ? "pho_et" : "cent_value"))
      throw std::runtime_error("model cell has no applicable working-point validation");
    // Released cuts also retain calibration bins and efficiency diagnostics.
    // The selection consumes only the sealed affine coefficients.
    const auto coefficients=[&](const std::string &name) {
      const auto &cut=w.at("cuts").at(name);
      return line(Json{{"intercept",cut.at("intercept")},{"slope",cut.at("slope")}},name);
    };
    const Json tight=coefficients(tightWP), high=coefficients(nonTightWP);
    for (const auto &x : system == "pp" ? et : cent) {
      const double t = tight["intercept"].get<double>() + tight["slope"].get<double>() * x.get<double>();
      const double hh = high["intercept"].get<double>() + high["slope"].get<double>() * x.get<double>();
      if (!(nonTightMinimum < hh && hh <= t && t <= 1.0)) throw std::runtime_error("selected working-point bands overlap or are invalid");
    }
    if (!ids.insert(id).second) throw std::runtime_error("duplicate working-point cells");
    ets.insert(ets.end(), {et[0].get<double>(), et[1].get<double>()});
    cents.insert(cents.end(), {cent[0].get<double>(), cent[1].get<double>()});
    cells.push_back(Json{{"id", id}, {"et", et}, {"centrality", cent}, {"tight", tight}, {"non_tight_high", high}});
  }
  if (ids.size() != wp.at("cells").size()) throw std::runtime_error("duplicate or missing working-point cells");
  /* every elementary rectangle is covered exactly once (half-open bounds, as in scoring) */
  std::sort(ets.begin(), ets.end()), ets.erase(std::unique(ets.begin(), ets.end()), ets.end());
  std::sort(cents.begin(), cents.end()), cents.erase(std::unique(cents.begin(), cents.end()), cents.end());
  for (size_t i = 0; i + 1 < ets.size(); ++i)
    for (size_t j = 0; j + 1 < cents.size(); ++j) {
      const double e = (ets[i] + ets[i + 1]) / 2, c = (cents[j] + cents[j + 1]) / 2;
      int n = 0;
      for (const auto &cell : cells)
        n += cell["et"][0].get<double>() <= e && e < cell["et"][1].get<double>() && cell["centrality"][0].get<double>() <= c &&
             c < cell["centrality"][1].get<double>();
      if (n != 1) throw std::runtime_error("optional model cells overlap or leave a gap");
    }
  Json result={{"score_branch", score}, {"state_branch", nullptr}, {"score_metadata", "meta/photon_models"},
              {"model_name", release.at("name")}, {"release_sha256", release.at("release_sha256")}, {"domain", d},
              {"cells", cells}, {"selection_sha256", sha256File(path.string())}, {"working_points_sha256", sha256File(wpPath.string())},
              {"working_point_choice", {{"tight",tightWP},{"non_tight_upper",nonTightWP},{"non_tight_lower",nonTightMinimum}}}};
  if (spec.contains("et_ranges_gev")) result["et_ranges_gev"]=spec.at("et_ranges_gev");
  if (wp.contains("calibration_domain") && wp.at("calibration_domain")!=d) throw std::runtime_error("working points were calibrated on another model domain");
  if (release.at("recipe").at("algorithm")=="bdt") {
    if (release.at("recipe").at("features")!=features.at("order") || features.at("view")!="H70")
      throw std::runtime_error("native model feature order differs from the retained H70 contract");
    Json nativeCells=Json::array(); bool native=true;
    for (size_t i=0;i<release.at("cells").size();++i) {
      const auto &cell=release.at("cells").at(i);
      const std::filesystem::path artifact=folder/cell.at("artifact").get<std::string>();
      const Json model=readJson(artifact.string());
      if (!model.contains("native_evaluator")) {native=false;break;}
      const auto &evaluation=model.at("native_evaluator");
      const std::filesystem::path filename=std::filesystem::canonical(artifact.parent_path()/evaluation.at("file").get<std::string>());
      const std::string relativeModel=filename.lexically_relative(std::filesystem::canonical(folder)).string();
      if (evaluation.at("type")!="ROOT_RBDT" || evaluation.at("inputs")!=features.at("count") ||
          model.at("features")!=features.at("order") || model.at("native_export_validation").at("status")!="PASS" ||
          relativeModel.empty() || relativeModel.rfind("..",0)==0 || !release.at("artifacts").contains(relativeModel))
        throw std::runtime_error("native model artifact is not bound to the released H70 inputs");
      Json row=cells.at(i);row["file"]=std::filesystem::canonical(filename).string();row["key"]=evaluation.at("key");
      row["sha256"]=release.at("artifacts").at(relativeModel);
      if(sha256File(filename.string())!=row.at("sha256").get<std::string>())throw std::runtime_error("native model bytes changed before scientific capture");
      row["scientific_function"]=photonjet_science::readRBDTDefinition(filename.string(),evaluation.at("key"),features.at("order"));
      nativeCells.push_back(row);
    }
    if (native) result["native_evaluation"]={{"type","ROOT_RBDT"},{"feature_contract",features},{"domain",d},{"cells",nativeCells}};
  }
  return result;
}

}  /* namespace */

std::string sha256Text(const std::string &text) {
  Sha256 h;
  h.update(reinterpret_cast<const unsigned char *>(text.data()), text.size());
  return h.hex();
}

std::string sha256File(const std::string &path) {
  /* Retry only read completion, before trusting any digest. Scientific
     failures and mismatched hashes never enter this recovery path. */
  for (int attempt=0;attempt<3;++attempt) {
    std::ifstream in(path,std::ios::binary); Sha256 h; std::vector<char> buffer(1<<22);
    if (in) {
      while (in) {
        ReductionProgress::pulse("hashing_payload");
        in.read(buffer.data(),buffer.size()); h.update(reinterpret_cast<const unsigned char *>(buffer.data()),in.gcount());
        if (!photonjet::progressPath.empty() && std::chrono::duration<double>(std::chrono::steady_clock::now()-photonjet::progressLast).count()>=2)
          photonjet::reportProgress(photonjet::entryBegin+photonjet::progressRowProcessed);
      }
      if (in.eof()) return h.hex();
    }
    if (attempt<2) std::this_thread::sleep_for(std::chrono::milliseconds(250));
  }
  throw std::runtime_error("cannot complete input hash after three read attempts: " + path);
}

std::string canonical(const Json &value) { return value.dump(-1, ' ', true); }
static void writeCanonical(std::ostream &out,const Json &value) {
  if (value.is_object()) {
    out<<'{'; bool first=true;
    for (const auto &[key,item] : value.items()) {
      if (!first) out<<','; first=false;
      out<<Json(key).dump(-1,' ',true)<<':'; writeCanonical(out,item);
    }
    out<<'}';
  } else if (value.is_array()) {
    out<<'['; bool first=true;
    for (const auto &item : value) {if (!first) out<<','; first=false; writeCanonical(out,item);}
    out<<']';
  } else out<<value.dump(-1,' ',true);
}
std::string digest(const Json &value) {
  return digestStream([&](std::ostream &out) {writeCanonical(out,value);});
}

/* Hash canonical metadata as it is written. Large source inventories need
   neither a second JSON tree nor one giant serialized string in memory. */
std::string digestStream(const std::function<void(std::ostream &)> &write) {
  class Buffer : public std::streambuf {
   public:
    Sha256 hash;
    std::streamsize xsputn(const char *data,std::streamsize size) override {
      hash.update(reinterpret_cast<const unsigned char *>(data),size); return size;
    }
    int_type overflow(int_type ch) override {
      if (!traits_type::eq_int_type(ch,traits_type::eof())) {
        const unsigned char byte=traits_type::to_char_type(ch); hash.update(&byte,1);
      }
      return traits_type::not_eof(ch);
    }
  } buffer;
  std::ostream out(&buffer); write(out);
  if (!out) throw std::runtime_error("cannot hash canonical metadata");
  return buffer.hash.hex();
}

Json readJson(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  return Json::parse(in);
}

Settings Settings::read(const std::string &path, const std::string &home) {
  TEnv env;
  photonjet_config::readSettings(path, env);
  const std::filesystem::path here = std::filesystem::absolute(path).parent_path(), root = std::filesystem::absolute(home);
  auto relative = [&](const std::string &p) { return std::filesystem::path(p).is_absolute() ? std::filesystem::path(p) : here / p; };
  Settings s;
  s.path = std::filesystem::absolute(path).string();
  s.home = root.string();
  std::set<std::string> known = {"Selection.Name", "Selection.PhotonIDPackage", "Selection.AuAuIsolationPackage", "Selection.PhotonIDModels", "Axes.XjEdges",
                                 "Axes.RecoResponseEdges", "Axes.TruthResponseEdges", "Fill.Preselection", "Fill.RequireStage3Fractions", "Fill.ResponseVariation", "Fill.EMSContract", "Fill.PhotonJetRecords", "Fill.RecordMaxRows", "Fill.RecordMinEt", "Fill.PhotonJetRecordsMatchMode", "Fill.JetRadius",
                                 "Fill.DataJetScale.pp", "Fill.DataJetScale.auau", "Fill.OutputMode", "Fill.ResponseGrid"};

  /* the frozen packages */
  const auto canonicalId = (root / "../PhotonID/models/nominal.json").lexically_normal(), canonicalIso = root / "config/auau_isolation_wp80.json";
  const auto idPath = value(env, "Selection.PhotonIDPackage").empty() ? canonicalId : relative(value(env, "Selection.PhotonIDPackage"));
  const auto isoPath = value(env, "Selection.AuAuIsolationPackage").empty() ? canonicalIso : relative(value(env, "Selection.AuAuIsolationPackage"));
  const Json pid = package(idPath.string(), "CanonicalPhotonIDFreezeV1"), iso = package(isoPath.string(), "CanonicalIsolationFreezeV1");
  const std::string pidSha = sha256File(idPath.string()), isoSha = sha256File(isoPath.string());
  if (iso.at("system") != "auau" || iso.at("witness").at("reconstructed") != true)
    throw std::runtime_error("the Au+Au isolation package must use a reconstructed witness");
  const Json canonicalPid = package(canonicalId.string(), "CanonicalPhotonIDFreezeV1");

  Json resolved = {{"schema", "PhotonJetSelectionV1"}, {"version", 1}, {"systems", Json::object()}};
  for (const std::string system : {"pp", "auau"}) {
    const std::string prefix = system + ".";
    for (const char *k : {"ScoreBranch", "StateBranch", "ScoreMetadata", "Tight", "NonTightLow", "NonTightHigh", "IsolationRadius",
                          "IsolationThreshold", "IsolationSidebandGap", "PhotonIDSelection", "Model", "ScoreMode",
                          "TightWP", "NonTightWP", "NonTightMinimum", "DataTrigger", "SimulationPeriods", "VertexFitDirectory"})
      known.insert(prefix + k);
    auto has = [&](const char *k) { return !value(env, prefix + k).empty(); };
    const Json &spec = pid.at(system);
    if (!std::regex_match(spec.at("model").at("sha256").get<std::string>(), kSha)) throw std::runtime_error("invalid model sha256 for " + system);
    Json cuts, isolation;
    if (system == "auau") {
      cuts = {{"tight", line(spec.at("working_points").at("WP70"), "WP70")}, {"non_tight_low", {{"intercept", 0.1}, {"slope", 0.0}}},
              {"non_tight_high", line(spec.at("working_points").at("WP80"), "WP80")}};
      isolation = {{"radius", finite(iso.at("witness").at("radius"), "radius")}, {"subtraction_method", 2},
                   {"threshold", {{"intercept", finite(iso.at("threshold").at("intercept_gev"), "isolation intercept")},
                                  {"slope", finite(iso.at("threshold").at("slope_gev_per_percent"), "isolation slope")}}},
                   {"sideband_gap_gev", finite(iso.at("sideband_gap_gev"), "isolation gap")},
                   {"centrality_min", finite(iso.at("domain").at("centrality_min"), "centrality min")},
                   {"centrality_max", finite(iso.at("domain").at("centrality_max_exclusive"), "centrality max")},
                   {"package_sha256", isoSha}};
      if (iso.at("witness").at("subtraction_method") != 2 || iso.at("threshold").at("variable") != "cent_value")
        throw std::runtime_error("the Au+Au package requires SUB1 and a centrality threshold");
    } else {
      /* the frozen pp contract states the PPG12 lines as formulae; a replacement package may give explicit working_points */
      if (spec.contains("working_points") && !spec.at("working_points").empty()) {
        for (const char *k : {"tight", "non_tight_low", "non_tight_high"}) cuts[k] = line(spec.at("working_points").at(k), k);
      } else {
        if (spec.at("regions") != canonicalPid.at("pp").at("regions")) throw std::runtime_error("a different pp ID definition needs explicit working_points");
        cuts = {{"tight", {{"intercept", 0.815625}, {"slope", -0.0015625}}},
                {"non_tight_low", {{"intercept", 0.7333333333333333}, {"slope", -0.01333333333333333}}},
                {"non_tight_high", {{"intercept", 0.684375}, {"slope", 0.0015625}}}};
      }
      isolation = {{"radius", 0.4}, {"subtraction_method", 1}, {"threshold", "stored"}, {"sideband_gap_gev", nullptr},
                   {"centrality_min", 0.0}, {"centrality_max", 80.0}, {"package_sha256", nullptr},
                   {"simulation_scale",1.2},{"simulation_shift_gev",0.1},{"non_isolated_max_gev",20.0}};
    }
    const Json originalCuts = cuts, originalIsolation = isolation;
    const std::string trigger = has("DataTrigger") ? value(env, prefix + "DataTrigger") : system == "auau" ? "photon12" : "photon4";
    if ((system == "auau" && trigger != "photon12") || (system == "pp" && trigger != "photon4"))
      throw std::runtime_error("nominal DATA triggers are photon12 for Au+Au and photon4 for p+p");
    for (const auto &[key, name] : {std::pair<const char *, const char *>{"Tight", "tight"}, {"NonTightLow", "non_tight_low"}, {"NonTightHigh", "non_tight_high"}})
      if (has(key)) cuts[name] = line(prefix + key, value(env, prefix + key));
    if (has("IsolationRadius")) isolation["radius"] = numbers(prefix + "IsolationRadius", value(env, prefix + "IsolationRadius")).at(0);
    if (has("IsolationThreshold")) isolation["threshold"] = line(prefix + "IsolationThreshold", value(env, prefix + "IsolationThreshold"));
    if (has("IsolationSidebandGap")) isolation["sideband_gap_gev"] = numbers(prefix + "IsolationSidebandGap", value(env, prefix + "IsolationSidebandGap")).at(0);
    const double radius = isolation["radius"].get<double>();
    if (radius != 0.3 && radius != 0.4) throw std::runtime_error("isolation radius must be a stored witness: 0.3 or 0.4");
    if (system == "auau" && isolation["radius"] != originalIsolation["radius"]) throw std::runtime_error("a new Au+Au radius requires a matching isolation package");
    if (system == "pp" && radius != 0.4 && !has("IsolationThreshold")) throw std::runtime_error("a changed pp radius needs its explicit threshold");
    if (!isolation["sideband_gap_gev"].is_null() && isolation["sideband_gap_gev"].get<double>() < 0) throw std::runtime_error("isolation sideband gap cannot be negative");
    if (!(0 <= isolation["centrality_min"].get<double>() && isolation["centrality_min"].get<double>() < isolation["centrality_max"].get<double>() &&
          isolation["centrality_max"].get<double>() <= 80))
      throw std::runtime_error("isolation centrality domain is outside the retained nominal domain");
    for (double x : system == "pp" ? std::vector<double>{15.0, 35.0} : std::vector<double>{0.0, 80.0}) {
      auto at = [&](const char *k) { return cuts[k]["intercept"].get<double>() + cuts[k]["slope"].get<double>() * x; };
      if (!(at("non_tight_low") < at("non_tight_high") && at("non_tight_high") <= at("tight")))
        throw std::runtime_error(system + ": tight and non-tight cuts overlap or have an empty band");
    }
    Json branches = {{"score_branch", has("ScoreBranch") ? value(env, prefix + "ScoreBranch") : "pho_bdt_score"},
                     {"state_branch", has("StateBranch") ? value(env, prefix + "StateBranch") : "pho_bdt_state"},
                     {"score_metadata", has("ScoreMetadata") ? value(env, prefix + "ScoreMetadata") : "meta/photon_id"}};
    for (const char *k : {"score_branch", "state_branch"})
      if (!std::regex_match(branches[k].get<std::string>(), kBranch)) throw std::runtime_error("score/state must name eventTree branches");
    if (!std::regex_match(branches["score_metadata"].get<std::string>(), std::regex("meta/[A-Za-z0-9_/]+")))
      throw std::runtime_error("score_metadata must name a meta/ receipt");
    const Json defaults = {{"score_branch", "pho_bdt_score"}, {"state_branch", "pho_bdt_state"}, {"score_metadata", "meta/photon_id"}};
    if (branches != defaults && branches["score_metadata"] == "meta/photon_id" &&
        (branches["score_branch"] != "pho_bdt_score" || branches["state_branch"] != "pho_bdt_state"))
      throw std::runtime_error("an alternative score needs its own augmentation metadata");
    const bool nominalScore = branches == defaults && spec.at("model").at("sha256") == canonicalPid.at(system).at("model").at("sha256");
    Json entry = branches;
    entry.update(cuts);
    entry.update(Json{{"model_sha256", spec.at("model").at("sha256")}, {"photon_id_package_sha256", pidSha}, {"isolation", isolation},
                      {"data_trigger", trigger},
                      {"event_selection",
                       {{"vertex_abs_z_max_cm", system == "auau" ? 10.0 : 60.0}, {"require_minimum_bias", system == "auau"},
                        {"trigger_word", "scaled"},
                        {"trigger_name", system == "auau" ? "Photon 12 GeV + MBD NS >= 2, vtx < 150 cm" : "Photon 4 GeV"},
                        {"data_run_periods", system == "auau" ? Json::object() : Json{{"0mrad", {47289, 51274}}, {"1p5mrad", {51274, 54000}}}}}},
                      {"verify_nominal_id", nominalScore && cuts == originalCuts && pidSha == sha256File(canonicalId.string())},
                      // pp stored classes precede the SIM transform and20GeV end;
                      // derive the requested classes from raw witnesses at fill time.
                      {"verify_nominal_isolation", system == "auau" && isolation == originalIsolation && isoSha == sha256File(canonicalIso.string())}});
    if (has("Model") && has("PhotonIDSelection"))
      throw std::runtime_error("choose Model or the explicit PhotonIDSelection path, not both");
    std::filesystem::path modelSelection;
    if (has("PhotonIDSelection")) modelSelection=relative(value(env,prefix+"PhotonIDSelection"));
    const std::string modelName=has("Model") ? value(env,prefix+"Model") : std::string();
    const bool originalModel=(system=="pp" && modelName=="ppg12_original") ||
                             (system=="auau" && modelName=="schema10_original");
    if (has("Model") && !originalModel) {
      if (!std::regex_match(modelName,kBranch)) throw std::runtime_error("invalid PhotonID model name");
      const auto models=value(env,"Selection.PhotonIDModels").empty() ? (root/"../PhotonID/models").lexically_normal()
                       : relative(value(env,"Selection.PhotonIDModels"));
      modelSelection=models/(system=="pp" ? "PP" : "AuAu")/modelName/"selection.json";
      if (!std::filesystem::is_regular_file(modelSelection))
        throw std::runtime_error("selected model release is unavailable: "+modelSelection.string());
    }
    if (!modelSelection.empty()) {
      for (const char *k : {"ScoreBranch", "StateBranch", "ScoreMetadata", "NonTightLow"})
        if (has(k)) throw std::runtime_error("an optional release owns its score and working points; do not override " + prefix + k);
      const std::string tight=has("TightWP") ? value(env,prefix+"TightWP") : "WP70";
      const std::string high=has("NonTightWP") ? value(env,prefix+"NonTightWP") : "WP80";
      double low=0.1;
      if (has("NonTightMinimum")) {
        const auto values=numbers(prefix+"NonTightMinimum",value(env,prefix+"NonTightMinimum"));
        if (values.size()!=1) throw std::runtime_error("NonTightMinimum needs one numeric score");
        low=values[0];
      }
      Json optional = optionalSelection(modelSelection, system,spec.at("features"),tight,high,low);
      const bool explicitLines=has("Tight") || has("NonTightHigh");
      if (explicitLines) {
        if (system!="auau" || !has("Tight") || !has("NonTightHigh") || has("TightWP") || has("NonTightWP"))
          throw std::runtime_error("explicit release working points need AuAu Tight and NonTightHigh together and no WP aliases");
        // The non-tight lower edge is a selection nuisance. optionalSelection
        // validates its finite score range; each explicit band must still be disjoint.
        const Json tightLine=cuts.at("tight"), highLine=cuts.at("non_tight_high");
        for (auto &cell : optional["cells"]) {
          for (const auto &x : cell.at("centrality")) {
            const double t=tightLine.at("intercept").get<double>()+tightLine.at("slope").get<double>()*x.get<double>();
            const double h=highLine.at("intercept").get<double>()+highLine.at("slope").get<double>()*x.get<double>();
            if (!(low<h && h<=t && t<=1.0)) throw std::runtime_error("explicit AuAu working-point bands overlap or are invalid");
          }
          cell["tight"]=tightLine;cell["non_tight_high"]=highLine;
        }
        if(optional.contains("native_evaluation")) for(auto &cell : optional["native_evaluation"]["cells"]) {
          cell["tight"]=tightLine;cell["non_tight_high"]=highLine;
        }
        optional["release_working_points_sha256"]=optional.at("working_points_sha256");
        optional.erase("working_points_sha256");
        optional["working_point_choice"]={{"mode","explicit"},{"variable","cent_value"},
          {"tight",tightLine},{"non_tight_upper",highLine},{"non_tight_lower",low}};
      }
      // Named models use stored scores normally. The explicit selection-path
      // interface retains native evaluation for pinned production recipes.
      const std::string mode=has("ScoreMode") ? value(env,prefix+"ScoreMode") : (has("Model") ? "stored" : "native");
      if (mode!="stored" && mode!="native") throw std::runtime_error("ScoreMode must be stored or native");
      if (mode=="native" && !optional.contains("native_evaluation"))
        throw std::runtime_error("selected model has no qualified native evaluator");
      if (mode=="stored") optional.erase("native_evaluation");
      optional["score_mode"]=mode;
      for (const char *k : {"score_branch", "state_branch", "score_metadata"}) entry[k] = optional.at(k);
      entry["non_tight_low"]={{"intercept",low},{"slope",0.0}};
      entry["optional_model"] = optional;
      entry["verify_nominal_id"] = false;
    } else if (originalModel) {
      if (has("ScoreMode") && value(env,prefix+"ScoreMode")!="stored")
        throw std::runtime_error("the original frozen model uses its verified stored score");
      if (has("TightWP") || has("NonTightWP") || has("NonTightMinimum"))
        throw std::runtime_error("the original frozen model uses its established cuts; use explicit line overrides for a variation");
      entry["model_choice"]=modelName;
      entry["verify_nominal_id"]=false;
    } else if (has("ScoreMode") || has("TightWP") || has("NonTightWP") || has("NonTightMinimum")) {
      throw std::runtime_error("model/working-point options require a selected Model or PhotonIDSelection");
    }
    /* p+p simulation: the PPG12 sum over crossing-angle periods (photonjet::configurePpPeriods) */
    if (system == "pp") {
      const auto directory=value(env,prefix+"VertexFitDirectory");
      const std::string fitDirectory=directory.empty()?"":std::filesystem::weakly_canonical(std::filesystem::path(path).parent_path()/directory).string();
      entry["simulation_periods"] = ppSimulationPeriods(root, has("SimulationPeriods") ? value(env, prefix + "SimulationPeriods") : "all",fitDirectory);
    } else if (has("SimulationPeriods") || has("VertexFitDirectory")) throw std::runtime_error("simulation-period and vertex-fit locations apply to p+p only");
    resolved["systems"][system] = entry;
  }
  const auto preselection=value(env,"Fill.Preselection");
  if(!preselection.empty() && preselection!="newPPG12")throw std::runtime_error("Fill.Preselection must be canonical newPPG12");
  resolved["preselection"]=photonjet::canonicalPreselection();
  resolved["timing_reference"]=photonjet::timingReferenceContract();
  resolved["fill_diagnostics"]={
    {"data_timecut",{{"formula","isfinite(pho_time_ns) && abs(pho_time_ns+2)<6"},{"units","ns"},{"scope","DATA_only_both_systems"},{"nominal_changed",false},{"leader_recomputed",true},{"event_records",true}}},
    {"leading_global",{{"chronology","global_stored_leading_after_finite_pt_overlap_then_eta_recoil_veto_no_fallback"},{"match","matched_signal_photon_and_direct_leading_jet_deltaR_le_0.3"},{"response_grid_restriction",false},{"et_bins_gev",1},{"xj_bins",0.05}}},
    {"reco_pair_truth_reason",{{"schema","RecoPairTruthReasonV1"},{"bit_semantics","independent_positive_facts"},{"truth_fa_et_window_gev",{17,33}},
      {"bits",{"photon_matched_to_response_truth","final_analysis_jet_link","truth_photon_in_native_response_grid","valid_response_truth_photon","linked_truth_jet_pt_pass","linked_truth_jet_eta_pass","linked_truth_jet_overlap_pass","linked_truth_jet_recoil_pass","linked_truth_xj_in_native_grid","linked_partner_available","linked_partner_already_taken","any_stored_truth_jet_nearby","above_threshold_stored_truth_jet_nearby","stored_truth_jet_collection_empty","linked_truth_jet_nonfinite","truth_photon_in_17_33"}}}}
  };
  resolved["selection_sha256"] = digest(resolved);
  resolved["profile"] = "selection:" + resolved["selection_sha256"].get<std::string>();
  resolved["name"] = required(env, "Selection.Name");
  s.selection = resolved;

  /* storage axes: xJ on [0, 3]; the response grids keep their support and every photon interval edge */
  for (const auto &[key, name] : {std::pair<const char *, const char *>{"Axes.XjEdges", "native_xj_edges"},
                                  {"Axes.RecoResponseEdges", "native_reco_response_edges_gev"},
                                  {"Axes.TruthResponseEdges", "native_truth_response_edges_gev"}}) {
    const auto e = numbers(key, required(env, key));
    if (e.size() < 2) throw std::runtime_error(std::string(key) + " needs at least two edges");
    for (size_t i = 1; i < e.size(); ++i) if (!(e[i] > e[i - 1])) throw std::runtime_error(std::string(key) + " must increase");
    s.axes[name] = e;
  }
  auto has = [](const Json &edges, double x) { for (const auto &e : edges) if (e.get<double>() == x) return true; return false; };
  if (s.axes["native_xj_edges"].front() != 0.0 || s.axes["native_xj_edges"].back() != 3.0) throw std::runtime_error("xJ storage support must be [0, 3]");
  for (const auto &[name, low, high] : {std::tuple<const char *, double, double>{"native_reco_response_edges_gev", 10, 40},
                                        {"native_truth_response_edges_gev", 5, 40}}) {
    const Json &e = s.axes[name];
    if (e.front() != low || e.back() != high) throw std::runtime_error(std::string(name) + " must span its response support");
    for (double edge : {15.0, 17.0, 23.0, 33.0, 35.0})
      if (!has(e, edge)) throw std::runtime_error(std::string(name) + " must keep every photon interval edge");
  }
  const std::string stage3 = required(env, "Fill.RequireStage3Fractions");
  if (stage3 != "0" && stage3 != "1") throw std::runtime_error("Fill.RequireStage3Fractions is 0 or 1");
  s.requireStage3Fractions = stage3 == "1";
  s.responseVariation = value(env, "Fill.ResponseVariation");
  const auto selectedVariation=photonjet::response::parseVariation(s.responseVariation.c_str());
  const auto emsPath=value(env,"Fill.EMSContract");
  if(photonjet::response::emsSelectionResponse(selectedVariation)){
    if(emsPath.empty())throw std::runtime_error("EMS requires explicit Fill.EMSContract");
    s.emsSelectionContract=readJson(relative(emsPath).string());
    photonjet::ems_native::validateContract(s.emsSelectionContract);
    if(s.emsSelectionContract.at("selection")!=s.selection)throw std::runtime_error("EMS contract selection differs");
  }else if(!emsPath.empty())throw std::runtime_error("Fill.EMSContract requires explicit EMS composition mode");
  photonjet::response::requireRunnableVariation(selectedVariation,!s.emsSelectionContract.is_null());
  const auto cache=value(env,"Fill.PhotonJetRecords");
  const auto gridPath=value(env,"Fill.ResponseGrid");
  if(!gridPath.empty()) {
    if(!s.responseVariation.empty() || !s.emsSelectionContract.is_null())throw std::runtime_error("Fill.ResponseGrid owns its response recipes; leave Fill.ResponseVariation/EMSContract empty");
    s.responseGrid=readJson(relative(gridPath).string());
    photonjet::response::native_grid::validateSettings(s.responseGrid);
  }
  if(!cache.empty() && cache!="0" && cache!="1")throw std::runtime_error("Fill.PhotonJetRecords is 0 or 1");
  s.candidateCache=cache=="1";
  const auto cacheRows=value(env,"Fill.RecordMaxRows");
  if(!cacheRows.empty()){if(cacheRows.find_first_not_of("0123456789")!=std::string::npos)throw std::runtime_error("invalid candidate cache row limit");s.candidateCacheMaxRows=std::stol(cacheRows);}
  if(s.candidateCacheMaxRows<1 || s.candidateCacheMaxRows>10000000)throw std::runtime_error("candidate cache row limit outside 1..10000000");
  const auto cacheEt=value(env,"Fill.RecordMinEt");
  if(!cacheEt.empty()){size_t used=0;s.candidateCacheMinEt=std::stod(cacheEt,&used);if(used!=cacheEt.size())throw std::runtime_error("invalid candidate cache ET floor");}
  if(!std::isfinite(s.candidateCacheMinEt)||s.candidateCacheMinEt<0||s.candidateCacheMinEt>12)throw std::runtime_error("candidate cache ET floor outside 0..12 GeV");
  /* The jet radius is a choice of the fill, not of the photon selection: the selection identity, the weights bound
     to it and the header inspection stay valid for every radius. A package records its radius in its histogram
     settings, so a collection cannot mix radii. */
  if (!value(env, "Fill.JetRadius").empty()) {
    const auto radius = numbers("Fill.JetRadius", value(env, "Fill.JetRadius"));
    if (radius.size() != 1 || (radius[0] != 0.2 && radius[0] != 0.3 && radius[0] != 0.4))
      throw std::runtime_error("Fill.JetRadius is 0.2, 0.3 or 0.4");
    s.jetRadius = radius[0];
  }
  /* The in-situ jet energy scale of DATA, per system: a calibration applied to DATA jets when they are read. */
  for (const auto &[key, target] : {std::pair<const char *, double *>{"Fill.DataJetScale.pp", &s.dataJetScalePp},
                                    {"Fill.DataJetScale.auau", &s.dataJetScaleAuAu}}) {
    if (value(env, key).empty()) continue;
    const auto scale = numbers(key, value(env, key));
    if (scale.size() != 1 || !(scale[0] > 0.8 && scale[0] < 1.25)) throw std::runtime_error(std::string(key) + " is one factor between 0.8 and 1.25");
    *target = scale[0];
  }
  TIter next(env.GetTable());
  while (auto *record = static_cast<TEnvRec *>(next()))
    if (!known.count(record->GetName())) throw std::runtime_error(path + ": unknown setting " + record->GetName());
  return s;
}

}  /* namespace treetohists */

/* === Package =========================================================================================================
   Package: the object inventory of a histogram package and the merge of one sample's fill jobs into it. */

/* A package holds, for one DATA collection or simulated sample, every histogram of the inventory summed over its
   source files (explicit Sumw2, flow bins included), the event records of the whole-event bootstrap (records/events,
   photons, jets, with 64-bit event references), the sparse per-source contributions of the measurement and response
   cells (records/source_contributions) with each source's identity and metadata (records/sources) for the source
   bootstrap, and meta/package, meta/inventory. The DATA trigger, runs and sources are chosen once, by the selection:
   the package holds only the analysis population. */
namespace treetohists {
namespace fs = std::filesystem;

/* One fill job's products: WORK/<system>/<sample>/<source / 1000>/<source>/, so no directory holds 1000 jobs. */
std::string productDirectory(const std::string &work, const std::string &sample, long source) {
  char shard[32], leaf[32];
  std::snprintf(shard, sizeof shard, "%03ld", source / 1000);
  std::snprintf(leaf, sizeof leaf, "%06ld", source);
  return (fs::path(work) / sample / shard / leaf).string();
}

/* Every object of a package: its path, axes, fill job product ("kind") and the key the filler writes it under. */
Json inventory(const std::string &system, const std::string &family) {
  const bool auau = system == "AuAu", sim = family != "data";
  const std::vector<std::string> classes = auau ? std::vector<std::string>{"cent0_20", "cent20_50", "cent50_80"} : std::vector<std::string>{"all"};
  const std::vector<std::string> labels = !sim ? std::vector<std::string>{"data"}
                                          : family == "photonjet" ? std::vector<std::string>{"sig", "phoother"} : std::vector<std::string>{"bkg", "jetprompt"};
  const std::string F = !sim ? "data" : family == "photonjet" ? "phosim" : "jetsim";
  const std::vector<std::string> shapes = {"e11_over_e33", "e32_over_e35", "bdt", "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx", "weta35_cogx", "wphi53_cogx", "et1", "et2", "et3", "et4"};
  const std::vector<std::string> widths = {"weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx", "weta35_cogx", "wphi53_cogx"};
  const std::vector<std::pair<std::string, std::string>> response = {
      {"photon_reco", "reco_et"}, {"photon_truth", "truth_et"}, {"photon_response_truth_x_reco", "truth_et,reco_et"},
      {"photon_fakes", "reco_et"}, {"photon_misses", "truth_et"}, {"xj_reco", "reco_et,xj"}, {"xj_truth", "truth_et,xj"},
      {"xj_response_truth_global_x_reco_global", "truth_cell,reco_cell"}, {"xj_fakes", "reco_et,xj"}, {"xj_fake_photon_reco", "reco_et,xj"},
      {"xj_detector_fakes_reco", "reco_et,xj"}, {"xj_misses", "truth_et,xj"}, {"xj_combinatoric_reco", "reco_et,xj"},
      {"xj_leading_reco", "reco_et,xj"}, {"xj_leading_matched_reco", "reco_et,xj"},
      {"xj_leading_global_reco","reco_et_1gev,xj_005"},{"xj_leading_global_matched_reco","reco_et_1gev,xj_005"},
      {"reco_pair_truth_reason","reco_et,xj,reason_bitmask"},
      {"signal_leakage_truth_matched_abcd", "region_abcd,reco_et"}, {"dphi_truthjet_near", "reco_et,xj,dphi"},
      {"dphi_truthjet_far", "reco_et,xj,dphi"}, {"dphi_combinatoric_near", "reco_et,xj,dphi"}, {"dphi_combinatoric_far", "reco_et,xj,dphi"}};
  Json rows = Json::array();
  auto add = [&](const std::string &path, const std::string &axes, const std::string &kind, const std::string &key = "") {
    std::string k = key;
    if (k.empty()) for (char c : path) k += c == '/' ? std::string("__") : std::string(1, c);
    rows.push_back(Json{{"path", path}, {"axes", axes}, {"kind", kind}, {"key", k}});
  };
  for (const auto &c : classes) {
    add("measurement/" + c + "/events", "count", "spectra", "events_" + c);
    add("measurement/" + c + "/photons", "et", "spectra", "photons_" + c);
    add("measurement/" + c + "/abcd", "region_abcd,et", "spectra", "abcd_" + c);
    for (int t : {5, 7, 10, 12}) {
      const std::string jet = "jet" + std::to_string(t);
      add("measurement/"+c+"/"+jet+"/leading_global","reco_et_1gev,xj_005","spectra","leading_global_"+c+"_"+jet);
      add("measurement/" + c + "/" + jet + "/recoil", "region_ac,et,xj", "spectra", "recoil_" + c + "_" + jet);
      for (const char *region : {"A", "C"})
        for (const char *side : {"near", "far"})
          add("measurement/" + c + "/" + jet + "/dphi" + region + "_" + side, "et,xj,dphi", "spectra",
              std::string("dphi") + region + "_" + c + "_" + jet + "_" + side);
      if (sim) for (const auto &[n, axes] : response) add("response/" + c + "/" + jet + "/" + n, axes, "response", n + "_" + c + "_" + jet);
      if(sim && !auau)for(const std::string n:{"xj_leading_reco_u","xj_leading_matched_reco_u","xj_leading_matched_any_reco_u"})
        add("response/"+c+"/"+jet+"/"+n,"reco_et_1gev,xj_005","response",n+"_"+c+"_"+jet);
    }
    for (const auto &label : labels) {
      const std::string dest = "qa/photons/" + c + "/" + label, tag = label + "_" + c;
      for (const auto &v : shapes) add(dest + "/shape_" + v, v + ",et,stage", "qa", "ss_" + v + "_" + tag);
      for (const auto &[n, axes] : std::vector<std::pair<std::string, std::string>>{
               {"preselection_compare", "stored_flag,canonical_npb,et"}, {"ncb", "npb,et"}, {"iso", "iso,et,id"}, {"abcd", "region_abcd,et"}, {"recoil", "region_ac,et,xj"},
               {"scoreiso", "score,iso,et"}, {"scoreiso_fine", "score,iso,et"}, {"clu_et", "wide_et,stage"}, {"clu_etaphi", "eta,phi,stage"}})
        add(dest + "/" + n, axes, "qa", n + "_" + tag);
      if (auau) add(dest + "/isonoue", "iso,et,id", "qa", "isonoue_" + tag);
      if (sim) {
        add(dest + "/abcdown", "region_abcd,et,generator_window", "qa", "abcdown_" + tag);
        add(dest + "/recoilown", "region_ownership,et,xj", "qa", "recoilown_" + tag);
      }
    }
    for (const auto &[n, axes] : std::vector<std::pair<std::string, std::string>>{
             {"evt_vz", "vz"}, {"jet_pt", "jet_pt,event_selection"}, {"jet_etaphi", "jet_eta,jet_phi,event_selection"}})
      add("qa/events/" + c + "/" + n, axes, "qa", n + "_" + F + "_" + c);
    if (auau && sim) add("qa/events/"+c+"/evt_vz_novertexweight","vz","qa","evt_vz_novertexweight_"+F+"_"+c);
    if (!sim) {
      add("qa/photons/" + c + "/data/e11e33_dt", "e11_over_e33,et,time_minus_mbd", "qa", "e11e33_dt_data_" + c);
      add("qa/photons/" + c + "/data/etatime", "eta,time_minus_mbd", "qa", "etatime_data_" + c);
      add("qa/photons/"+c+"/data/time_stage","corrected_time,et,stage","qa","clu_time_stage_data_"+c);
      add("qa/photons/"+c+"/data/etatime_stage","eta,corrected_time,stage","qa","clu_etatime_stage_data_"+c);
      add("qa/photons/"+c+"/data/etatime_npb","eta,corrected_time,npb_score","qa","clu_etatime_npb_data_"+c);
      add("qa/photons/"+c+"/data/rawtime_stage","raw_time,et,stage","qa","clu_rawtime_stage_data_"+c);
    }
    if (family == "photonjet") {
      add("qa/truth/" + c + "/truth_iso", "truth_iso,truth_et,prompt_class", "qa", "truth_iso_" + c);
      add("qa/truth/" + c + "/eff", "truth_et,efficiency_stage,reco_et", "qa", "eff_" + c);
    }
    for (const char *radius : {"02", "03", "04"}) add("qa/jets/" + c + "/R" + radius + "/calo_fractions", "fraction,jet_pt,calo_component", "extensions");
    if (sim) for (const char *obj : {"photon", "jet"}) add("qa/matching/" + c + "/" + obj, "delta_r,reco_pt,match_class", "extensions");
    else add("qa/photons/" + c + "/data/e11e33_npb_timing", "e11_over_e33,et,time_minus_mbd", "extensions");
  }
  for (const auto &label : labels) {
    const std::string dest = "qa/slices/" + label;
    add(dest + "/bdt", "score,et,centrality_slice", "qa", "bdt_" + label);
    add(dest + "/isoslice", "iso,et,centrality_slice", "qa", "isoslice_" + label);
    if (auau) add(dest + "/isonoueslice", "iso,et,centrality_slice", "qa", "isonoueslice_" + label);
    for (const auto &v : widths) {
      add(dest + "/width_" + v, v + ",et,centrality_slice", "qa", "width_" + v + "_" + label);
      add(dest + "/width_tight_iso_" + v, v + ",et,centrality_slice", "qa", "width_tight_iso_" + v + "_" + label);
    }
  }
  add("qa/events/calo_total", "calo_total,centrality,mb_state", "qa", "evt_calo_" + F);
  if (auau) {
    add("qa/events/centrality", "centrality,weight_mode", "qa", "evt_cent_" + F);
    if(sim)add("qa/events/vertex","vz,weight_mode","qa","evt_vzw_"+F);
    if (!sim) {
      add("qa/events/centrality_scaled_bits", "centrality,scaled_bit", "qa", "cent_bits_data");
      add("qa/events/vertex_selected", "vz", "qa", "evt_vz_selected_data");
      add("qa/events/vertex_selected_wide","vz","qa","evt_vz_selected_wide_data");
    }
  }
  if (sim) {
    add("qa/generator/leading_pt", "generator_pt,generator_kind", "qa", "gen_" + F);
    add("qa/generator/raw_leading_pt", "generator_pt,generator_kind", "qa", "gen_raw_" + F);
    if (!auau) {
      add("qa/stitching/reco_tag_vz30", "photon_et", "qa", "reco_tag_vz30_" + F);
    }
    add("qa/photon_id/score_roles", "score,et,centrality_slice,truth_role", "extensions");
  } else {
    add("qa/events/trigger_bits_raw", "max_cluster_energy,scaled_bit,live_bit", "qa", "trig_data");
    add("qa/events/trigger_study_bits_raw", "max_cluster_energy,scaled_bit,live_bit", "qa", "trig_study_data");
    add("qa/events/trigger_study_bits_corrected", "max_cluster_energy,scaled_bit,live_bit", "extensions", "derived_trigger_study");
    add("qa/events/trigger_max_cluster_energy", "max_cluster_energy,semantic_trigger,trigger_state", "extensions");
  }
  for (const auto &[n, axes] : std::vector<std::pair<std::string, std::string>>{
           {"cutflow", "gate,centrality,weight_mode"}, {"weights", "event_weight,weight_state"}, {"mbd_charge", "mbd_charge,centrality"}})
    add("qa/events/" + n, axes, "extensions");
  for (const std::string d : {"emcal", "ihcal", "ohcal", "total"}) {
    add("qa/events/" + d + "_vs_mbd_charge", "calo_" + d + ",mbd_charge,mb_state", "extensions");
    if (d != "total") add("qa/events/" + d + "_vs_centrality", "calo_" + d + ",centrality,mb_state", "extensions");
  }
  if(!sim){const Json original=rows;for(const auto&row:original)if(row.at("kind")=="spectra")
    add("variations/data_timecut/"+row.at("path").get<std::string>(),row.at("axes"),"spectra","timecut/"+row.at("key").get<std::string>());}
  return rows;
}

Json packageInventory(const Json &plan,const Json &sample) {
  Json rows=inventory(sample.at("system"),sample.at("family"));
  if(sample.at("family")=="data" || !plan.contains("response_component_grid"))return rows;
  const auto &config=plan.at("response_component_grid");photonjet::response::native_grid::validateSettings(config);
  const Json nominal=rows;
  for(const auto &variant:config.at("variants"))if(variant.at("name")!="nominal")for(auto row:nominal)if(row.at("kind")=="response") {
    const std::string name=variant.at("name");row["path"]="response_components/"+name+"/"+row.at("path").get<std::string>();
    row["kind"]="response__"+name;rows.push_back(row);
  }
  return rows;
}

namespace {

std::string text(TFile &f, const std::string &key) {
  auto *o = dynamic_cast<TObjString *>(f.Get(key.c_str()));
  if (!o) throw std::runtime_error(std::string(f.GetName()) + ": no " + key);
  return o->GetString().Data();
}

TDirectory *directory(TFile &file, const std::string &path) {
  if (path.empty()) return &file;
  if (!file.GetDirectory(path.c_str())) file.mkdir(path.c_str(), "", true);
  auto *d = file.GetDirectory(path.c_str());
  if (!d) throw std::runtime_error("cannot create package directory " + path);
  return d;
}

void sameAxes(const TH1 &a, const TH1 &b) {
  if (a.GetDimension() != b.GetDimension() || a.GetNcells() != b.GetNcells()) throw std::runtime_error("histogram shape mismatch");
  const TAxis *aa[] = {a.GetXaxis(), a.GetYaxis(), a.GetZaxis()}, *bb[] = {b.GetXaxis(), b.GetYaxis(), b.GetZaxis()};
  for (int d = 0; d < a.GetDimension(); ++d) {
    if (aa[d]->GetNbins() != bb[d]->GetNbins()) throw std::runtime_error("histogram axis mismatch");
    for (int i = 1; i <= aa[d]->GetNbins() + 1; ++i)
      if (aa[d]->GetBinLowEdge(i) != bb[d]->GetBinLowEdge(i)) throw std::runtime_error("histogram edges mismatch");
  }
}

/* The saved package: every inventory object present (optional extensions excepted) with finite values and variances,
   event references inside the events, one source record per source. */
void check(const std::string &path, const Json &rows, size_t sources) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot reopen " + path);
  for (const auto &row : rows) {
    TObject *o = f->Get(row.at("path").get<std::string>().c_str());
    if (!o) {
      if (row.at("kind") == "extensions") continue;
      throw std::runtime_error("missing package object " + row.at("path").get<std::string>());
    }
    if (auto *h = dynamic_cast<TH1 *>(o)) {
      if (h->GetSumw2N() != h->GetNcells()) throw std::runtime_error("no explicit Sumw2: " + row.at("path").get<std::string>());
      for (int c = 0; c < h->GetNcells(); ++c)
        if (!std::isfinite(h->GetBinContent(c)) || !std::isfinite(h->GetSumw2()->At(c)))
          throw std::runtime_error("nonfinite cell in " + row.at("path").get<std::string>());
    }
  }
  auto *events = dynamic_cast<TTree *>(f->Get("records/events"));
  if (!events) throw std::runtime_error("no event records");
  for (const char *name : {"records/photons", "records/jets"}) {
    auto *t = dynamic_cast<TTree *>(f->Get(name));
    if (!t) throw std::runtime_error(std::string("no ") + name);
    t->SetBranchStatus("*", false);
    t->SetBranchStatus("event", true);
    TLeaf *leaf = t->GetLeaf("event");
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      const Long64_t e = leaf->GetValueLong64();
      if (e < 0 || e >= events->GetEntries()) throw std::runtime_error(std::string("broken event reference in ") + name);
    }
  }
  auto *s = dynamic_cast<TTree *>(f->Get("records/sources"));
  if (!s || s->GetEntries() != static_cast<Long64_t>(sources)) throw std::runtime_error("source records lost");
}

}  /* namespace */

/* Sample publication combines histograms AND event/source records from the
   pinned fill jobs. Source contributions survive this step because
   FinalAnalysis resamples whole simulated sources, not individual bins. */
void merge(const Json &plan, const Json &sample, const std::string &work, const std::string &outputDir) {
  checkReductionCode(plan);
  checkVertexFit(sample);
  checkPpPeriodFits(ppPeriodsOf(plan, sample));
  const std::string key = sample.at("key"), system = sample.at("system"), family = sample.at("family");
  const fs::path target = fs::path(outputDir) / (key + ".root"), partial = target.string() + ".partial";
  if (sample.contains("run_packages")) {
    std::vector<std::string> inputs;
    for (const auto &run : sample.at("run_packages")) inputs.push_back((fs::path(outputDir)/(run.get<std::string>()+".root")).string());
    combinePackages(plan,sample,inputs,target.string());
    return;
  }
  ReductionProgress progress("merge",(sample.contains("_source_count") ? sample.at("_source_count").get<size_t>() : sample.at("sources").size())+1);
  if (fs::exists(target)) {
    validatePackage(target.string(), plan, sample);
    progress.finish(); return;  /* exact validated publication is reusable */
  }
  fs::create_directories(target.parent_path());
  if (sample.contains("reused_package")) {
    const auto &saved=sample.at("reused_package");
    const std::string original=saved.at("path");
    if (sha256File(original)!=saved.at("sha256").get<std::string>()) throw std::runtime_error("reused package changed: " + original);
    validatePackage(original,plan,sample);
    std::error_code error;
    fs::create_hard_link(original,target,error);
    if (error) fs::copy_file(original,target);
    validatePackage(target.string(),plan,sample);
    std::cout<<"[TreeToHists] reused "<<key<<" from "<<original<<'\n';
    progress.finish(); return;
  }
  fs::remove(partial);
  const Json rows = packageInventory(plan,sample);
  const std::string selectionText = plan.at("selection").dump(), profile = plan.at("selection").at("profile");
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(), "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot create " + partial.string());
  auto *records = directory(*out, "records");
  records->cd();
  Long64_t source = 0, bin = 0, contributionBegin = 0, contributionEnd = 0;
  double value = 0, variance = 0;
  std::string histogram, sourceJson;
  TTree contributions("source_contributions", "Sparse response, measurement and run-dependent QA cells, including ROOT flow bins");
  contributions.Branch("source", &source);
  contributions.Branch("histogram", &histogram);
  contributions.Branch("bin", &bin);
  contributions.Branch("value", &value);
  contributions.Branch("variance", &variance);
  contributions.SetAutoFlush(-8 * 1024 * 1024);
  TTree sources("sources", "Source identities, input hashes, producer metadata and normalization");
  sources.Branch("source", &source);
  sources.Branch("json", &sourceJson);
  sources.Branch("contribution_begin", &contributionBegin);
  sources.Branch("contribution_end", &contributionEnd);
  std::map<std::string, std::unique_ptr<TH1>> sums;
  std::map<std::string, std::unique_ptr<THnSparse>> sparse;
  RecordTable events, photons, jets, timeEvents, timePhotons, timeJets;
  CandidateCachePackage candidateCache;
  Long64_t eventOffset = 0, timeEventOffset = 0;
  Json available = Json::object();
  Json nativeAncestry=Json::array();
  std::string grids;
  Json dataJesFusion;
  std::vector<std::string> kinds = {"spectra", "qa", "extensions"};
  if (family != "data") kinds.push_back("response");
  if(family!="data" && plan.contains("response_component_grid"))for(const auto &v:plan.at("response_component_grid").at("variants"))
    if(v.at("name")!="nominal")kinds.push_back("response__"+v.at("name").get<std::string>());
  long long denominator = 0, missing = 0;  /* the canonical Au+Au embedded stitching: the events the fill jobs counted */
  for (const auto &src : sample.at("sources")) {
    const long sourceIndex = src.value("index", source);
    std::vector<Json> pieces;
    for (const auto &job : plan.at("jobs"))
      if (job.at("sample") == key && job.at("source").get<long>() == sourceIndex) pieces.push_back(job);
    if (pieces.empty()) throw std::runtime_error("no fill rows for original source");
    std::sort(pieces.begin(), pieces.end(), [](const Json &a, const Json &b) {
      return a.value("range_begin", 0LL) < b.value("range_begin", 0LL);
    });
    long long next = 0;
    std::map<std::pair<std::string, long>, std::pair<double, double>> cells;
    contributionBegin = contributions.GetEntries();
    Json ranges = Json::array();
    for (const auto &job : pieces) {
      const long long begin = job.value("range_begin", 0LL), end = job.value("range_end", src.at("metadata").at("entries").get<long long>());
      if (begin != next || end < begin || (end==begin && src.at("metadata").at("entries")!=0)) throw std::runtime_error("gapped or overlapping source ranges");
      next = end;
      if (end>begin) ranges.push_back({begin, end});
      const fs::path dir = jobDirectory(work, job);
      if (!fs::exists(dir / "done.json")) throw std::runtime_error("fill job incomplete: " + dir.string());
      const Json done = readJson((dir / "done.json").string());
      if (done.at("plan_sha256") != plan.at("plan_sha256") || done.at("input") != src.at("path") ||
          done.at("sha256") != src.at("expected_sha256") ||
          done.value("producer_code_sha256", "") != digest(plan.at("producer_code")) ||
          done.value("range_begin", 0LL) != begin || done.value("range_end", end) != end)
        throw std::runtime_error("fill belongs to another source, range, plan or code: " + dir.string());
      if (done.contains("stitching_denominator")) {
        denominator += done.at("stitching_denominator").at("events").get<long long>();
        missing += done.at("stitching_denominator").at("missing").get<long long>();
      }
      for (const auto &kind : kinds) {
        const fs::path product = dir / (kind + ".root");
        if (!done.contains("products") || done.at("products").at(kind).get<std::string>() != sha256File(product.string()))
          throw std::runtime_error("filled product hash differs: " + product.string());
        std::unique_ptr<TFile> in(TFile::Open(product.c_str(), "READ"));
        if (!in || in->IsZombie()) throw std::runtime_error("cannot open product " + product.string());
        nativeAncestry.push_back(photonjet_science::producerInput(*in));
        if (text(*in, "profile") != profile || text(*in, "selection") != selectionText)
          throw std::runtime_error("product selection differs: " + dir.string());
        if(plan.contains("response_component_grid") && (kind=="response" || kind.rfind("response__",0)==0)) {
          const std::string label=kind=="response"?"nominal":kind.substr(10);
          const auto provenance=Json::parse(text(*in,"native_response_grid"));
          const auto &spec=provenance.at("grid");
          if(spec.value("ems_selection_contract",Json())!=plan.at("response_component_grid").value("ems_selection_contract",Json()))
            throw std::runtime_error("response component EMS contract differs from planned configuration");
          const Json expectedSource={{"path",src.at("path")},{"sha256",src.at("expected_sha256")},{"input_file_sha256",src.at("identity")},{"entries",src.at("metadata").at("entries")},{"entry_begin",begin},{"entry_end",end}};
          if(spec.at("sources")!=Json::array({expectedSource}))throw std::runtime_error("response component original source/range differs");
          Json resolved=Json::array();for(const auto &v:plan.at("response_component_grid").at("variants")) {
            Json row=v;row["plan"]=photonjet::response::native_grid::resolvedPlan(spec,row);row["plan_sha256"]=digest(row["plan"]);resolved.push_back(row);
          }
          if(spec.at("variants")!=resolved)throw std::runtime_error("response component resolved plan roster differs");
          if(provenance.at("scope")!="RESPONSE_COMPONENT_ONLY" || provenance.at("grid_sha256")!=digest(spec) ||
             spec.at("producer_code")!=plan.at("producer_code") || spec.at("limits")!=plan.at("response_component_grid").at("limits"))
            throw std::runtime_error("response component code/configuration binding differs");
          bool matched=false;
          for(const auto &v:plan.at("response_component_grid").at("variants"))if(v.at("name")==label) {
            matched=provenance.at("variant").at("variation")==v.at("variation") &&
              Json::parse(text(*in,"variation"))==Json::parse(photonjet::response::variationJson(photonjet::response::parseVariation(v.at("variation").get<std::string>().c_str())));
          }
          if(!matched)throw std::runtime_error("response component variant binding differs");
        }
        if(kind=="spectra" && in->Get("data_jes_fusion")) {
          const auto provenance=Json::parse(text(*in,"data_jes_fusion"));
          if(!done.contains("data_jes_fusion") || done.at("data_jes_fusion").at("science_traversals")!=1 || provenance.at("qa_extensions")!="nominal_reused" || provenance.at("scale")!=plan.at("histograms").at("data_jet_scale").at("pp"))
            throw std::runtime_error("DATA JES spectrum lacks shared-traversal lineage");
          auto binding=provenance;binding.erase("science_readers");binding.erase("original_entries");
          if(!dataJesFusion.is_null() && dataJesFusion!=binding)throw std::runtime_error("DATA JES view lineage differs");
          dataJesFusion=binding;
        }
        if (kind == "response") {
          const std::string g = text(*in, "grids");
          if (!grids.empty() && grids != g) throw std::runtime_error("response grids differ between sources");
          grids = g;
        }
        if (kind == "extensions") {
          const Json a = Json::parse(text(*in, "availability"));
          for (const auto &[name, state] : a.items()) {
            if (!available.contains(name)) available[name] = state;
            else if (available.at(name) != state) available[name] = "PARTIAL";
          }
        }
        for (const auto &row : rows) {
          if (row.at("kind").get<std::string>() != kind) continue;
          const std::string path = row.at("path");
          TObject *obj = in->Get(row.at("key").get<std::string>().c_str());
          if (!obj) {
            if (kind != "extensions") throw std::runtime_error("missing filled histogram " + row.at("key").get<std::string>());
            continue;
          }
          if (auto *h = dynamic_cast<TH1 *>(obj)) {
            if (h->GetSumw2N() != h->GetNcells()) throw std::runtime_error("histogram without explicit Sumw2: " + path);
            const bool runQA = kind == "qa" && (path == "qa/events/centrality_scaled_bits" || path == "qa/events/trigger_bits_raw" || path == "qa/events/trigger_study_bits_raw" || path == "qa/events/centrality" || path == "qa/events/vertex_selected" || path == "qa/events/vertex_selected_wide" || path == "qa/events/vertex");
            for (int cell = 0; cell < h->GetNcells(); ++cell) {
              if (cell%10000==0) ReductionProgress::pulse("merging_histograms");
              const double v = h->GetBinContent(cell), e = h->GetSumw2()->At(cell);
              if (!std::isfinite(v) || !std::isfinite(e) || e < 0) throw std::runtime_error("invalid histogram cell in " + path);
              if ((kind == "response" || kind.rfind("response__",0)==0 || kind == "spectra" || runQA) && (v != 0 || e != 0)) {
                auto &sum = cells[{path, cell}];
                sum.first += v, sum.second += e;
              }
            }
            if (!sums.count(path)) {
              sums[path].reset(static_cast<TH1 *>(h->Clone())); sums[path]->SetDirectory(nullptr);
            } else {
              sameAxes(*sums[path], *h);
              if (!sums[path]->Add(h)) throw std::runtime_error("histogram addition failed: " + path);
              photonjet_science::mergeDefinitions(*sums[path],*h);
            }
          } else if (auto *h = dynamic_cast<THnSparse *>(obj)) {
            if (!sparse.count(path)) sparse[path].reset(static_cast<THnSparse *>(h->Clone()));
            else sparse[path]->Add(h);
          } else throw std::runtime_error("unsupported object " + path);
        }
        if (kind == "spectra") {
          auto *ev = dynamic_cast<TTree *>(in->Get("table_events"));
          auto *ph = dynamic_cast<TTree *>(in->Get("table_photons"));
          auto *jt = dynamic_cast<TTree *>(in->Get("table_jets"));
          if (!ev || !ph || !jt) throw std::runtime_error("event records missing from " + dir.string());
          events.append(*ev, "events", source, eventOffset, records);
          photons.append(*ph, "photons", source, eventOffset, records);
          jets.append(*jt, "jets", source, eventOffset, records);
          if(plan.at("histograms").value("candidate_records",false))candidateCache.append(*in,"table_candidates_",source,records);
          eventOffset += ev->GetEntries();
          if (family == "data") {
            auto *te = dynamic_cast<TTree *>(in->Get("timecut/table_events"));
            auto *tp = dynamic_cast<TTree *>(in->Get("timecut/table_photons"));
            auto *tj = dynamic_cast<TTree *>(in->Get("timecut/table_jets"));
            if (!te || !tp || !tj) throw std::runtime_error("timing event records missing from " + dir.string());
            auto *timing = directory(*out, "variations/data_timecut/records");
            timeEvents.append(*te, "events", source, timeEventOffset, timing);
            timePhotons.append(*tp, "photons", source, timeEventOffset, timing);
            timeJets.append(*tj, "jets", source, timeEventOffset, timing);
            timeEventOffset += te->GetEntries();
          }
        }
      }
    }
    if (next != src.at("metadata").at("entries").get<long long>()) throw std::runtime_error("source range does not cover the complete tree");
    if (sample.value("data_run",Json::object()).value("trigger_study",Json::object()).value("status","")=="AVAILABLE") {
      const auto &factors=sample.at("data_run").at("trigger_study").at("factors");
      const std::string raw="qa/events/trigger_study_bits_raw", corrected="qa/events/trigger_study_bits_corrected";
      const auto *shape=dynamic_cast<TH3 *>(sums.at(raw).get());
      if (!shape) throw std::runtime_error("trigger study histogram has the wrong shape");
      if (!sums.count(corrected)) {
        sums[corrected].reset(static_cast<TH1 *>(shape->Clone())); sums[corrected]->SetDirectory(nullptr); sums[corrected]->Reset();
      }
      std::map<std::pair<std::string,long>,std::pair<double,double>> derived;
      for (const auto &[cell,sum] : cells) if (cell.first==raw) {
        int x,y,z;shape->GetBinXYZ(cell.second,x,y,z);
        const std::string bit=std::to_string(y-1);
        if (!factors.contains(bit)) continue;
        const double factor=factors.at(bit), v=sum.first*factor, e=sum.second*factor*factor;
        derived[{corrected,cell.second}]={v,e};
        auto &h=*sums.at(corrected);h.AddBinContent(cell.second,v);h.GetSumw2()->AddAt(h.GetSumw2()->At(cell.second)+e,cell.second);
      }
      cells.insert(derived.begin(),derived.end());
    }
    // DATA readers use only these QA cells. Filter after the full raw study has
    // built its corrected aggregate above; changing that input would change physics.
    // Keep all study cells for unavailable run metadata (conservative fallback).
    const auto study=sample.value("data_run",Json::object()).value("trigger_study",Json::object());
    const bool studyAvailable=study.value("status","")=="AVAILABLE";
    const auto *studyShape=sums.count("qa/events/trigger_study_bits_raw")
      ? dynamic_cast<TH3 *>(sums.at("qa/events/trigger_study_bits_raw").get()) : nullptr;
    std::set<int> readerBits;
    if (family=="data" && studyAvailable && studyShape) {
      // Match FinalAnalysis sourceStudies: exact name, else FIRST preferred prefix.
      // dataRunRecord's factor lookup has different duplicate-prefix precedence.
      const auto &names=sample.at("data_run").at("trigger_record").at("name");
      const std::vector<std::string> targets=system=="AuAu"
        ? std::vector<std::string>{"MBD N&S >= 2, vtx < 150 cm","Photon 10 GeV + MBD NS >= 2","Photon 12 GeV + MBD NS >= 2"}
        : std::vector<std::string>{"MBD N&S >= 1","Photon 3 GeV + MBD NS >= 1","Photon 4 GeV + MBD NS >= 1","Photon 5 GeV + MBD NS >= 1"};
      for (const auto &target : targets) {
        int bit=-1,first=-1,preferred=-1;
        for (size_t i=0;i<names.size();++i) if (names[i]==target) {bit=i;break;}
        if (bit<0) {
          for (size_t i=0;i<names.size();++i) {
            const std::string name=names[i];
            if (name.rfind(target,0)!=0) continue;
            if (first<0) first=i;
            if (preferred<0 && name.find("vtx < 150 cm")!=std::string::npos) preferred=i;
          }
          bit=preferred>=0 ? preferred : first;
        }
        if (bit<0) throw std::runtime_error("AVAILABLE trigger study lacks its reader bit");
        readerBits.insert(bit);
      }
    }
    for (const auto &[cell, sum] : cells) {
      if (family=="data") {
        const auto &path=cell.first;
        bool keep=path=="qa/events/centrality_scaled_bits" || path=="qa/events/vertex_selected" || path=="qa/events/centrality";
        if (path=="qa/events/trigger_bits_raw" && !studyShape) keep=true;
        if (path=="qa/events/trigger_study_bits_raw") {
          keep=true;
          if (studyAvailable && studyShape) {
            int x,y,z;studyShape->GetBinXYZ(cell.second,x,y,z);
            keep=readerBits.count(y-1) && (z==65 || readerBits.count(z-1));
          }
        }
        if (!keep) continue;
      }
      histogram = cell.first, bin = cell.second, value = sum.first, variance = sum.second;
      contributions.Fill();  /* one cell per original source, independent of job partitioning */
    }
    sourceJson = Json{{"path", src.at("path")}, {"sha256", src.at("expected_sha256")}, {"identity", src.at("identity")},
                      {"metadata", src.at("metadata")}, {"entry_ranges", ranges}, {"reduced_entries", next}}.dump();
    contributionEnd = contributions.GetEntries(); sources.Fill(); ++source; progress.advance(source);
  }
  for (auto &[path, h] : sums) {
    const auto slash = path.rfind('/');
    directory(*out, path.substr(0, slash))->cd();
    h->Write(path.substr(slash + 1).c_str());
  }
  for (auto &[path, h] : sparse) {
    const auto slash = path.rfind('/');
    directory(*out, path.substr(0, slash))->cd();
    h->Write(path.substr(slash + 1).c_str());
  }
  records->cd();
  contributions.Write("", TObject::kOverwrite);
  sources.Write("", TObject::kOverwrite);
  events.write(), photons.write(), jets.write();
  timeEvents.write(), timePhotons.write(), timeJets.write();
  candidateCache.write();
  out->cd();
  TObjString(profile.c_str()).Write("profile");
  TObjString(selectionText.c_str()).Write("selection");
  if (!grids.empty()) TObjString(grids.c_str()).Write("grids");
  const std::string sys = system == "PP" ? "pp" : "auau";
  Json normalization = sample.at("normalization");
  if (normalization.contains("ownership")) {  /* cross section over the sample's events in its denominator */
    if (missing > 0) throw std::runtime_error(key + ": " + std::to_string(missing) + " events without the ownership variable; the stitching denominator is not exact");
    if (denominator <= 0 && !sample.value("intermediate", false)) throw std::runtime_error(key + ": no events in the stitching denominator");
    normalization["denominator_events"] = denominator;
    normalization["pb_per_generated_event"] = denominator > 0 ? Json(normalization.at("cross_section_pb").get<double>() / denominator) : Json();
  }
  Json meta = {{"schema", "PhotonJetSampleV1"}, {"system", system}, {"name", sample.at("name")}, {"family", family},
               {"profile", profile}, {"selection", plan.at("selection")}, {"histograms", plan.at("histograms")},
               {"normalization", normalization}, {"centrality_weights", sample.at("weights")},
               {"centrality_evidence", sample.at("centrality_evidence")}, {"centrality_data_trigger", sample.at("centrality_data_trigger")},
               {"vertex_reweighting", sample.at("vertex_reweighting")}, {"plan_sha256", plan.at("plan_sha256")},
               {"runtime", {{"root_version", gROOT->GetVersion()}}}, {"histograms_written", sums.size() + sparse.size()},
               {"availability", available}, {"coverage", sample.value("coverage", Json::object())}, {"source_count", source}, {"event_records", eventOffset},
               {"acceptance", "STRUCTURE_VALIDATED_PHYSICS_NOT_EVALUATED"}};
  if(plan.contains("selection_fanout"))meta["selection_fanout"]=plan.at("selection_fanout");
  if(family!="data" && plan.contains("response_component_grid"))meta["response_components"]={{"scope","RESPONSE_COMPONENT_ONLY"},{"configuration",plan.at("response_component_grid")},{"ems_measurement_qa_fanout","UNQUALIFIED"}};
  if(!dataJesFusion.is_null())meta["data_jes_fusion"]=dataJesFusion;
  if(plan.at("histograms").value("candidate_records",false))meta["candidate_records"]=candidateRecordMetadata(candidateCache,plan,system);
  meta["package_binding_sha256"]=packageBinding(plan,sample);
  /* These values come from the same resolved Profile used by every fill. They are descriptive metadata,
     not a second selection: changing a cut changes the fill and this record together. */
  photonjet::configureJetRadius(plan.at("histograms").at("jet_radius").get<double>());
  const auto cuts = photonjet::configuredProfile("nominal", system == "PP" ? photonjet::System::pp : photonjet::System::auau,
                                                (fs::path(work) / "selection.json").c_str());
  Json qaClasses = Json::array();
  for (const auto &cls : photonjet::centralityClasses(cuts.system))
    qaClasses.push_back({{"id",cls.name},{"centrality_percent",cuts.system==photonjet::System::pp ? Json() : Json{cls.low,cls.high}}});
  Json angle = Json();
  for (int denominator=1;denominator<=64;++denominator) {
    const int numerator=std::lround(cuts.deltaPhiMin/M_PI*denominator);
    if (std::fabs(cuts.deltaPhiMin-M_PI*numerator/denominator)<1e-14) {angle={numerator,denominator};break;}
  }
  meta["render_selection"] = {
      {"schema", "PhotonJetRenderSelectionV1"}, {"selection_sha256", plan.at("selection").at("selection_sha256")},
      {"producer_code", plan.at("producer_code")},
      {"collision", {{"system", system}, {"sqrt_s_nn_gev", 200.0}, {"contract", "sPHENIX_eventTree_v1_200GeV"}}},
      {"photon", {{"et_min_gev", cuts.photonEtMin}, {"et_max_gev", cuts.photonEtMax}, {"abs_eta_max", cuts.photonAbsEtaMax}}},
      {"event", plan.at("selection").at("systems").at(sys).at("event_selection")},
      {"jet", {{"collection", cuts.jetCollection("jet")}, {"radius", cuts.jetRadius}, {"abs_eta_max", cuts.jetAbsEtaMax},
               {"data_pt_scale", family == "data" ? plan.at("histograms").at("data_jet_scale").at(sys).get<double>() : 1.0},
               {"pt_min_gev", cuts.jetPtMin}, {"overlap_delta_r", cuts.overlapDeltaR}, {"delta_phi_min_rad", cuts.deltaPhiMin},
               {"delta_phi_pi_fraction",angle}, {"pt_operator",">="}, {"delta_phi_operator",">="}, {"overlap_operator",">="}}},
      {"truth_photon", {{"isolation_radius", 0.3}, {"abs_eta_max", photonjet::kTruthPhotonAbsEtaMax},
                        {"isolation_max_gev", photonjet::kTruthIsolationMaxGeV}, {"signal_source_role", photonjet::kSignalSourceRole}}},
      {"qa", {{"cluster_map_et_min_gev",12.0}, {"et_min_gev", photonjet::qa::kQaEtMin}, {"et_max_gev", photonjet::qa::kQaEtMax},
              {"classes",qaClasses}, {"stored_abs_eta_max", photonjet::qa::kStoredAbsEtaMax}, {"jet_map_pt_min_gev", photonjet::qa::kQaJetMapPtMin}}},
      {"photon_id", plan.at("selection").at("systems").at(sys)},
      {"histogram_settings", plan.at("histograms")}};
  if (family == "data") {
    meta["data_trigger"] = plan.at("selection").at("systems").at(sys).at("data_trigger");
    meta["trigger_selection"] = "applied_in_package";  /* the measurement and records hold only this trigger's population */
    if (sample.contains("run")) meta["run"]=sample.at("run");
    if (sample.contains("data_run")) meta["data_run"]=sample.at("data_run");
  }
  directory(*out, "meta")->cd();
  meta["package_kind"] = sample.value("intermediate", false) ? "INTERMEDIATE" : "FINAL";
  if (system=="PP" && family!="data") meta["si_di_role"]=sample.at("sources").at(0).at("metadata").at("source").at("si_di_role");
  meta["tree_weight_copies"]=treeCopySummary(sample);
  TObjString(meta.dump().c_str()).Write("package");
  if (family != "data" && !normalization.at("pb_per_generated_event").is_null())  /* FinalAnalysis applies it once */
    TParameter<double>("pb_per_generated_event", normalization.at("pb_per_generated_event").get<double>()).Write();
  TObjString(rows.dump().c_str()).Write("inventory");
  photonjet_science::writeProducerEvidence(*out,nativeAncestry);
  contributions.SetDirectory(nullptr);
  sources.SetDirectory(nullptr);
  out->Close();
  sums.clear(); sparse.clear(); events.output.reset(); photons.output.reset(); jets.output.reset();
  publishValidatedPackage(partial.string(),target.string(),plan,sample);
  progress.finish();
  std::cout << "[TreeToHists] " << target.string() << ": " << source << " sources, " << meta.at("histograms_written") << " histograms, "
            << eventOffset << " event records" << std::endl;
}

}  /* namespace treetohists */

/* === Catalogue, plan and fill ========================================================================================
   plan resolves the selection, expands every catalogued sample into its eventTree files (loose files, or the files of a
   published stage-2 group), checks each file's system, family, source identity and score binding, and writes
   DIR/plan.json, DIR/selection.json and DIR/jobs.txt (one fill job per file). fill runs the fillers on one file into
   DIR/<system>/<sample>/<job>/ and records the input's sha256; a job whose products are complete is not redone. merge
   sums the jobs of one sample into OUT/<system>/<sample>.root (above). */
namespace treetohists {
namespace fs = std::filesystem;
/* Minimal, idempotent sourceInfo contract. Fill cuts, model evaluator bytes,
   timing and histogram choices cannot invalidate already-read headers. */
Json inspectionSelection(const Json &selection) {
  Json systems=Json::object();
  for(const auto &[system,spec] : selection.at("systems").items()) {
    const bool native=spec.contains("optional_model") && spec.at("optional_model").contains("native_evaluation");
    Json out={{"state_branch",spec.at("state_branch")},{"verify_nominal_isolation",false}};
    if(native) out["optional_model"]={{"native_evaluation",Json::object()}};
    else {
      out["score_branch"]=spec.at("score_branch");out["score_metadata"]=spec.at("score_metadata");
      if(spec.contains("optional_model")) {
        Json model=Json::object();
        for(const char *key : {"model_name","release_sha256","score_branch","domain"}) model[key]=spec.at("optional_model").at(key);
        out["optional_model"]=model;
      } else {
        out["model_sha256"]=spec.at("model_sha256");out["photon_id_package_sha256"]=spec.at("photon_id_package_sha256");
      }
    }
    if(system=="pp")out["simulation_periods"]={{"mode",spec.value("simulation_periods",Json::object()).value("mode","")=="all" ? "all" : "single"}};
    if(system=="auau" && spec.at("verify_nominal_isolation")==true) {
      out["verify_nominal_isolation"]=true;out["isolation"]={{"package_sha256",spec.at("isolation").at("package_sha256")}};
    }
    systems[system]=out;
  }
  return {{"systems",systems}};
}
bool inspectionCompatible(const Json &oldSnapshot,const Json &newSnapshot) {
  if(oldSnapshot.value("root_version",Json())!=newSnapshot.at("root_version"))return false;
  Json old=inspectionSelection(oldSnapshot.at("selection"));const Json now=inspectionSelection(newSnapshot.at("selection"));
  // A previous successful nominal-isolation check also satisfies a new cut
  // override that no longer requires that check. The reverse is forbidden.
  for(const auto &[system,spec] : now.at("systems").items())
    if(spec.at("verify_nominal_isolation")==false && old.at("systems").contains(system)) {
      old["systems"][system]["verify_nominal_isolation"]=false;old["systems"][system].erase("isolation");
    }
  return old==now;
}

namespace {

void writeText(const fs::path &path, const std::string &text) {  /* atomic: a reader never sees a partial file */
  char host[256] = "host";
  ::gethostname(host, sizeof(host) - 1);
  /* this writer's own temporary name: two jobs that write one file (a repeat of a slow fill group) never share it */
  const fs::path partial = path.string() + ".partial." + host + "." + std::to_string(::getpid());
  {
    std::ofstream out(partial);
    out << text;
    if (!out) throw std::runtime_error("cannot write " + partial.string());
  }
  fs::rename(partial, path);
}

/* YAML names the produced ROOT files or published groups and their supplied
   normalization/weighting bindings. The existing planner checks them below.
   Previously prepared JSON catalogues retain their original reader. */
Json readInputCollection(const std::string &path) {
  if (!photonjet_config::isYaml(path)) return readJson(path);
  Json collection = photonjet_config::readYaml(path);
  if (!collection.is_object()) throw std::runtime_error(path + ": inputs must be a YAML mapping");
  if (!collection.contains("schema")) collection["schema"] = "PhotonJetInputCollectionV1";
  return collection;
}

/* The files of a published stage-2 group: GROUP/CURRENT -> versions/<version>/ with a hash-bound MANIFEST.json. */
std::vector<std::pair<fs::path, std::string>> groupFiles(const fs::path &group) {
  const fs::path current = group / "CURRENT";
  if (!fs::is_symlink(current)) throw std::runtime_error("group is not published: " + group.string());
  const fs::path directory = fs::canonical(current);
  if (directory.parent_path() != fs::canonical(group / "versions")) throw std::runtime_error("group pointer escapes its version directory");
  const Json manifest = readJson((directory / "MANIFEST.json").string());
  Json body = manifest;
  body.erase("version");
  if (digest(body) != manifest.at("version").get<std::string>() || directory.filename() != manifest.at("version").get<std::string>())
    throw std::runtime_error("group publication manifest changed: " + group.string());
  std::vector<std::pair<fs::path, std::string>> out;
  std::set<std::string> names, identities;
  for (const auto &f : manifest.at("files")) {
    const fs::path name = f.at("relative_path").get<std::string>();
    if (name.is_absolute() || name.lexically_normal().string().rfind("..", 0) == 0 || name.extension() != ".root")
      throw std::runtime_error("unsafe publication file name: " + name.string());
    if (!names.insert(name.string()).second || !identities.insert(f.at("logical_id").get<std::string>()).second)
      throw std::runtime_error("duplicate file or source in group manifest: " + group.string());
    out.emplace_back(directory / name, f.at("sha256").get<std::string>());
  }
  return out;
}

/* One eventTree file: its source metadata (the record FinalAnalysis and the source bootstrap read) and its family. */
Json sourceInfo(const fs::path &path, const Json &selection, std::string &family, std::string &system) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path.string());
  const Json schema = Json::parse(photonjet::metaText(*f, "meta/schema")), source = Json::parse(photonjet::metaText(*f, "meta/source"));
  system = schema.at("system");
  const std::string dataset = source.at("dataset");
  family = dataset.rfind("isSim", 0) != 0 ? "data" : dataset.find("Inclusive") != std::string::npos ? "jetsim" : "phosim";
  auto *tree = dynamic_cast<TTree *>(f->Get("eventTree"));
  if (!tree) throw std::runtime_error(path.string() + ": no eventTree");
  const Json &spec = selection.at("systems").at(system);
  const bool native=spec.contains("optional_model") && spec.at("optional_model").contains("native_evaluation");
  std::vector<std::string> needed = {"pho_et", "isow_cone_sum", "isow_threshold", "isow_sideband_threshold"};
  if (native) {needed.push_back("shv_features");needed.push_back("shv_finite_feature_state");}
  else needed.push_back(spec.at("score_branch"));
  if (system == "pp" && family != "data") {  /* the canonical weight state the fills check (photonjet::PpSimulationGate) */
    photonjet::ppWeightComponent(source);
    for (const char *name : {"ppg12_weight_active_state", "ppg12_weight_lane_code", "w_application_count", "w_event", "w_slice", "w_vertex", "w_sidi", "w_period"})
      needed.push_back(name);
    if (spec.value("simulation_periods", Json::object()).value("mode", "") == "all")  /* the truth vertices of the period sum */
      for (const char *name : {"truth_vtx_z", "truth_mbvtx_z"}) needed.push_back(name);
  }
  if (!spec.at("state_branch").is_null()) needed.push_back(spec.at("state_branch"));
  for (const auto &name : needed)
    if (!tree->GetBranch(name.c_str())) throw std::runtime_error(path.string() + ": missing selection input " + name);
  const Json receipt = native ? Json::object() : Json::parse(photonjet::metaText(*f, spec.at("score_metadata")));
  if (spec.contains("optional_model") && !native) {
    const Json &o = spec.at("optional_model");
    const Json score = receipt.value(o.at("model_name").get<std::string>(), Json::object());
    for (const char *k : {"release_sha256", "score_branch", "domain"})
      if (score.value(k, Json()) != o.at(k)) throw std::runtime_error(path.string() + ": optional score is not bound to the selected release");
    if (score.value("status", "") != "PASS") throw std::runtime_error(path.string() + ": optional score is not validated");
  } else if (!native) {
    if (receipt.at("model_sha256") != spec.at("model_sha256") || receipt.value("freeze_sha256", "") != spec.at("photon_id_package_sha256").get<std::string>())
      throw std::runtime_error(path.string() + ": saved score differs from the selected photon-ID package");
  }
  if (spec.at("verify_nominal_isolation") == true && system == "auau" &&
      Json::parse(photonjet::metaText(*f, "meta/isolation")).at("freeze_sha256") != spec.at("isolation").at("package_sha256"))
    throw std::runtime_error(path.string() + ": saved nominal Au+Au isolation package differs");
  Json meta = {{"source", source}, {"entries", tree->GetEntries()}};
  for (const char *key : {"trigger_selection", "trigger_run_info", "normalization", "weighting", "event_accounting", "accounting", "consolidation"}) {
    const std::string name = std::string("meta/") + key;
    if (!f->Get(name.c_str())) continue;
    Json v = Json::parse(photonjet::metaText(*f, name));
    if (std::string(key) == "consolidation") {  /* keep the source and chunk identity, not the repeated chunk dictionaries */
      const std::string whole = digest(v);
      v.erase("chunks");
      v["input_metadata_sha256"] = whole;
    }
    meta[key] = v;
  }
  const std::string identity = source.value("logical_source_id", source.value("input_file_sha256", ""));
  if (identity.empty()) throw std::runtime_error("missing original source identity: " + path.string());
  return Json{{"path", path.string()}, {"identity", identity}, {"metadata", meta}};
}

/* Run-level GL1 records are shared by all sources of a run, never summed per
   source. Only trigger-study spectra use the legacy six-decimal live/scaled factor. */
Json dataRunRecord(const Json &sources, const Json &luminosity, const std::string &system) {
  Json info=Json(), study={{"status","UNAVAILABLE"},{"reason","Run GL1 live/scaled counters are not retained"}};
  for (const auto &src : sources) {
    const auto &m=src.at("metadata");
    const Json table=m.contains("trigger_run_info") ? m.at("trigger_run_info") : m.value("trigger_selection",Json());
    if (info.is_null()) info=table;
    else if (table!=info) throw std::runtime_error("DATA sources disagree on their run trigger record");
  }
  if (info.is_object() && info.contains("run_live") && info.contains("run_scaled")) {
    const auto &names=info.at("name");
    const std::vector<std::string> targets=system=="AuAu" ? std::vector<std::string>{"MBD N&S >= 2, vtx < 150 cm","Photon 10 GeV + MBD NS >= 2","Photon 12 GeV + MBD NS >= 2"} :
      std::vector<std::string>{"MBD N&S >= 1","Photon 3 GeV + MBD NS >= 1","Photon 4 GeV + MBD NS >= 1","Photon 5 GeV + MBD NS >= 1"};
    Json factors=Json::object(); bool valid=system!="PP" || sources.at(0).at("metadata").at("source").value("run",-1)>=47289;
    for (const auto &target : targets) {
      int bit=-1;
      for (size_t i=0;i<names.size();++i) if (names[i]==target) {bit=i;break;}
      if (bit<0) for (size_t i=0;i<names.size();++i) if (names[i].get<std::string>().rfind(target,0)==0) {
        if (bit<0 || names[i].get<std::string>().find("vtx < 150 cm")!=std::string::npos) bit=i;
      }
      if (bit<0 || info.at("initial_prescale").at(bit).get<int>()<0 || info.at("run_live").at(bit).is_null() || info.at("run_scaled").at(bit).is_null()) {valid=false;continue;}
      const double live=info.at("run_live").at(bit), scaled=info.at("run_scaled").at(bit);
      if (!(live>0 && scaled>0 && std::isfinite(live/scaled))) {valid=false;continue;}
      char text[64];std::snprintf(text,sizeof text,"%.6f",live/scaled);
      factors[std::to_string(bit)]=std::stod(text);
    }
    if (valid) study={{"status","AVAILABLE"},{"rule","run_live/run_scaled, rounded to six decimals"},{"factors",factors},{"vertex_abs_z_cm",30},{"energy_min_gev",1}};
    else study["reason"]="The study's common trigger group is disabled or has missing/zero GL1 counters";
  }
  return {{"trigger_record",info},{"luminosity",luminosity},{"trigger_study",study}};
}

/* Campaign planning keeps source metadata in pinned, ordered inventories.
   The scientific binding is the same as an interactive plan; only its
   storage changes. A DATA run is the largest metadata unit held in memory. */
Json sourceInventory(const fs::path &path,const Json &sources) {
  fs::create_directories(path.parent_path());
  const fs::path next=path.string()+".partial";
  std::ofstream out(next);
  for (const auto &source : sources) out<<source.dump(-1,' ',true)<<'\n';
  out.close();
  if (!out) throw std::runtime_error("cannot write source inventory " + path.string());
  const std::string sha=sha256File(next.string());
  if (fs::exists(path)) {
    if (sha256File(path.string())!=sha) throw std::runtime_error("existing source inventory differs: " + path.string());
    fs::remove(next);
  } else fs::rename(next,path);
  return {{"path",path.string()},{"sha256",sha},{"count",sources.size()}};
}
void externalizeSample(Json &sample,const Json &common,const fs::path &directory) {
  const std::string binding=packageBinding(common,sample);
  if (sample.at("system")=="PP" && sample.at("family")!="data") {
    std::set<std::string> roles;
    for (const auto &src : sample.at("sources")) roles.insert(src.at("metadata").at("source").at("si_di_role"));
    if (roles.size()!=1) throw std::runtime_error("mixed p+p SI/DI source roles in one sample");
    sample["si_di_role"]=*roles.begin();
  }
  sample["_source_file"]=sourceInventory(directory/(sample.at("key").get<std::string>()+".jsonl"),sample.at("sources"));
  sample["_source_count"]=sample.at("sources").size();
  sample["_package_binding_sha256"]=binding;
  sample.erase("sources");
}
Json readSourceInventory(const fs::path &path) {
  Json rows=Json::array(); std::ifstream in(path); std::string line;
  if (!in) throw std::runtime_error("cannot read source inventory " + path.string());
  while (std::getline(in,line)) rows.push_back(Json::parse(line));
  if (!in.eof()) throw std::runtime_error("source inventory read failed: " + path.string());
  return rows;
}

/* Published parts share a small number of directories. Resolve each parent
   once instead of walking the same remote directory hierarchy for every part.
   Leaf symlinks and filesystem errors retain weakly_canonical's behavior. */
class SourcePaths {
 public:
  fs::path resolve(const fs::path &path) {
    const auto parent=path.parent_path();
    auto found=parents_.find(parent);
    if (found==parents_.end()) found=parents_.emplace(parent,fs::weakly_canonical(parent)).first;
    const fs::path resolved=found->second/path.filename();
    std::error_code error;
    const auto status=fs::symlink_status(resolved,error);
    if (fs::is_symlink(status) || (error && error!=std::errc::no_such_file_or_directory))
      return fs::weakly_canonical(path);
    return resolved.lexically_normal();
  }
 private:
  std::map<fs::path,fs::path> parents_;
};

/* Read one bounded header-inspection result at a time, in frozen input order.
   Fills still verify the actual ROOT bytes and headers before reading events. */
class InspectedSources {
 public:
  InspectedSources(const Json &parts,const Json &snapshot): parts_(parts),snapshot_(snapshot) {}
  bool enabled() const {return !parts_.empty();}
  Json next(const fs::path &path,const std::string &sha,std::string &family,std::string &system) {
    if (row_==rows_.size()) {
      if (part_==parts_.size()) throw std::runtime_error("input inspection is truncated");
      const auto &pin=parts_.at(part_++);
      if (sha256File(pin.at("path"))!=pin.at("sha256").get<std::string>()) throw std::runtime_error("inspection descriptor changed");
      const Json descriptor=readJson(pin.at("path"));
      Json result=readJson(pin.at("result"));const std::string claimed=result.at("sha256");result.erase("sha256");
      if (digest(result)!=claimed || result.at("schema")!="TreeToHistsSourceInspectionV1" ||
          result.at("descriptor_sha256")!=pin.at("sha256") || result.at("snapshot")!=snapshot_ ||
          result.at("rows").size()!=descriptor.at("files").size()) throw std::runtime_error("inspection result has a different binding");
      rows_=std::move(result.at("rows"));expected_=descriptor.at("files");row_=0;
    }
    const Json &row=rows_.at(row_),&expected=expected_.at(row_++);
    if (row.at("input")!=expected || expected.at("path")!=path.string() || expected.at("sha256")!=sha ||
        row.at("info").at("path")!=path.string()) throw std::runtime_error("inspection source order or hash differs");
    family=row.at("family");system=row.at("system");return row.at("info");
  }
  void finish() const {
    if (enabled() && (part_!=parts_.size() || row_!=rows_.size())) throw std::runtime_error("unused source inspection records");
  }
 private:
  const Json &parts_,&snapshot_;
  Json rows_=Json::array(),expected_=Json::array();size_t part_=0,row_=0;
};

/* Expand and check every catalogued sample. Ordinary plans remain complete
   in-memory documents. Campaigns spool DATA metadata per run and retain
   compact sample contracts instead of a second full combined inventory. */
Json plan(const Settings &s, const std::string &catalogue, const std::vector<std::string> &only,
          long long maxEvents = 100000,const fs::path &sourceDirectory = {},const Json &inspections=Json::array()) {
  if (maxEvents <= 0 || maxEvents > 100000) throw std::runtime_error("max-events must be between 1 and 100000");
  const Json c = readInputCollection(catalogue);
  if (c.value("schema", "") != "PhotonJetInputCollectionV1" || c.value("samples", Json::array()).empty())
    throw std::runtime_error("expected PhotonJetInputCollectionV1 with samples");
  const fs::path base = fs::absolute(catalogue).parent_path();
  const bool bounded=!sourceDirectory.empty();
  Json common={{"selection",s.selection},{"histograms",histogramSettings(s)},
    {"producer_code",reductionCode()},{"root_version",gROOT->GetVersion()}};
  if(!s.responseGrid.is_null())common["response_component_grid"]=s.responseGrid;
  /* a header inspection is bound to the selection and ROOT only (Campaign.cc inspectionSnapshot) */
  const Json inspectionBinding={{"selection",inspectionSelection(s.selection)},{"root_version",gROOT->GetVersion()}};
  InspectedSources inspected(inspections,inspectionBinding);
  if (inspected.enabled() && !only.empty()) throw std::runtime_error("campaign inspection requires its complete frozen input order");
  Json samples = Json::array(), jobs = Json::array();
  std::set<std::string> keys, identities;
  std::map<std::string,std::vector<std::pair<fs::path,std::string>>> inputFiles;
  SourcePaths sourcePaths;
  long long totalFiles=0, processedFiles=0;
  for (const auto &sample : c.at("samples")) {
    const std::string key=sample.at("system").get<std::string>()+"/"+sample.at("name").get<std::string>();
    if (!only.empty() && std::find(only.begin(),only.end(),key)==only.end()) continue;
    auto &files=inputFiles[key];
    for (const auto &f : sample.value("files",Json::array())) {
      const fs::path supplied=base/(f.is_string() ? f.get<std::string>() : f.at("path").get<std::string>());
      // Inspected campaigns retain the exact path bound by their frozen list.
      files.emplace_back(inspected.enabled() ? supplied : sourcePaths.resolve(supplied),
                         f.is_string() ? "" : f.at("sha256").get<std::string>());
    }
    for (const auto &g : sample.value("groups",Json::array()))
      for (auto &f : groupFiles(fs::weakly_canonical(base/g.get<std::string>()))) files.push_back(f);
    totalFiles+=files.size();
  }
  if (!photonjet::progressPath.empty()) {
    photonjet::progressTotal=photonjet::progressRowTotal=totalFiles; photonjet::progressOffset=0;
    photonjet::entryBegin=0; photonjet::entryEnd=-1; photonjet::progressRowProcessed=0;
    photonjet::progressIdentity={{"campaign",sha256File(catalogue)},{"row","PLAN"},{"unit","sources"}};
    photonjet::progressPhase="inspecting_input_snapshot"; photonjet::progressStart=std::chrono::steady_clock::now();
    photonjet::reportProgress(0);
  }
  for (const auto &sample : c.at("samples")) {
    const std::string system = sample.at("system"), name = sample.at("name"), family = sample.at("family");
    if ((system != "PP" && system != "AuAu") || (family != "data" && family != "photonjet" && family != "inclusivejet") ||
        !std::regex_match(name, std::regex("[A-Za-z0-9_-]+")))
      throw std::runtime_error("invalid system, sample name or family");
    if ((name == "data") != (family == "data")) throw std::runtime_error("each system has one data sample named data");
    const std::string key = system + "/" + name;
    if (!only.empty() && std::find(only.begin(), only.end(), key) == only.end()) continue;
    if (!keys.insert(key).second) throw std::runtime_error("duplicate sample " + key);
    const auto &files=inputFiles.at(key);
    if (files.empty()) throw std::runtime_error(key + " has no input files");
    const std::string expected = family == "data" ? "data" : family == "photonjet" ? "phosim" : "jetsim";
    Json sources = Json::array();
    size_t selectedSources=0;
    std::map<int,fs::path> runInventories;
    std::ofstream runStream; int openRun=-1;
    std::set<std::string> populations;
    for (const auto &[path, sha] : files) {
      std::string fileFamily, fileSystem;
      Json info = inspected.enabled() ? inspected.next(path,sha,fileFamily,fileSystem) : sourceInfo(path,s.selection,fileFamily,fileSystem);
      if (fileSystem != (system == "PP" ? "pp" : "auau") || fileFamily != expected)
        throw std::runtime_error("sample identity disagrees with eventTree metadata: " + key + " " + path.string());
      if (!identities.insert(info.at("identity")).second) throw std::runtime_error("duplicate original source identity: " + info.at("identity").get<std::string>());
      populations.insert(info.at("metadata").at("source").value("sample", ""));
      /* Published bytes are checked by the fill that actually reads them.
         Its metadata must also match this inspected header before any event
         loop. Interactive/unpublished inputs are still hashed here. */
      const std::string observedSha = bounded && !sha.empty() ? sha : sha256File(path.string());
      if (!sha.empty() && observedSha != sha) throw std::runtime_error("published input differs from its group manifest: " + path.string());
      info["expected_sha256"] = observedSha;
      info["index"] = selectedSources;
      const long long entries = info.at("metadata").at("entries");
      if (entries < 0) throw std::runtime_error("invalid eventTree entry count: " + path.string());
      if (entries==0) jobs.push_back(Json{{"sample",key},{"source",selectedSources},{"input",path.string()},
        {"range_begin",0},{"range_end",0},{"source_entries",0}});
      for (long long begin = 0; begin < entries; begin += maxEvents)
        jobs.push_back(Json{{"sample", key}, {"source", selectedSources}, {"input", path.string()},
          {"range_begin", begin}, {"range_end", std::min(entries, begin + maxEvents)}, {"source_entries", entries}});
      if (bounded && family=="data") {
        const int run=info.at("metadata").at("source").value("run",-1);
        if (run<=0) throw std::runtime_error(key+": DATA source has no valid run number");
        const fs::path spool=sourceDirectory/"inspection"/system/("run_"+std::to_string(run)+".partial");
        const bool first=runInventories.emplace(run,spool).second;
        if (openRun!=run) {
          if (runStream.is_open()) {runStream.close(); if (!runStream) throw std::runtime_error("cannot spool DATA metadata");}
          fs::create_directories(spool.parent_path());
          runStream.clear(); runStream.open(spool,first ? std::ios::out : std::ios::out|std::ios::app); openRun=run;
        }
        runStream<<info.dump(-1,' ',true)<<'\n';
        if (!runStream) throw std::runtime_error("cannot spool DATA metadata");
      } else sources.push_back(std::move(info));
      ++selectedSources;
      photonjet::reportProgress(++processedFiles);
    }
    if (runStream.is_open()) {runStream.close(); if (!runStream) throw std::runtime_error("cannot spool DATA metadata");}
    if (family != "data" && (populations.size() != 1 || populations.count("")))
      throw std::runtime_error(key + ": a simulation package must hold exactly one named source population");
    if (sample.contains("source_sample") && populations != std::set<std::string>{sample.at("source_sample").get<std::string>()})
      throw std::runtime_error(key + ": source_sample binding differs");
    if (system == "PP" && family != "data") {  /* each source's own period is configured, and the producer weighted it with that period's fit */
      const Json periods = s.selection.at("systems").at("pp").value("simulation_periods", Json());
      checkPpPeriodFits(periods);
      if (periods.is_object() && periods.value("mode", "") == "all")
        for (const auto &src : sources) {
          const Json &m = src.at("metadata");
          const std::string own = m.at("source").value("period", "");
          const Json weighting = m.value("weighting", Json::object());
          if (!periods.at("periods").contains(own)) throw std::runtime_error(key + ": source period '" + own + "' is not a configured p+p period");
          if (weighting.value("composition", "") != "slice*vertex*sidi*period" ||
              weighting.value("vertex_file_sha256", "") != periods.at("periods").at(own).at("fit").at("sha256").get<std::string>())
            throw std::runtime_error(key + ": the producer did not weight this source with its period's truth-vertex fit");
        }
    }
    const Json weights = sample.value("centrality_weights", Json::array());
    Json vertex = nullptr;
    if (system == "AuAu" && family != "data") {
      if (weights.size() != 16) throw std::runtime_error(key + " needs its 16 production-derived centrality_weights");
      for (const auto &w : weights) if (!w.is_number() || !(w.get<double>() > 0) || !std::isfinite(w.get<double>())) throw std::runtime_error(key + ": invalid centrality weight");
      if (sample.value("centrality_evidence", "").empty()) throw std::runtime_error(key + " needs centrality_evidence identifying its derivation");
      if (sample.value("centrality_data_trigger", "") != s.selection.at("systems").at("auau").at("event_selection").at("trigger_name").get<std::string>())
        throw std::runtime_error(key + ": centrality weights must bind the selected DATA trigger");
      vertex = sample.value("vertex_reweighting", Json());
      if (!vertex.is_object() || vertex.value("evidence", "").empty()) throw std::runtime_error(key + " needs a qualified Au+Au vertex_reweighting fit");
      if (vertex.value("coordinate", "") != "truth_vtx_z" || vertex.value("selection_sha256", "") != s.selection.at("selection_sha256").get<std::string>() ||
          vertex.value("centrality_weights_sha256", "") != digest(weights))
        throw std::runtime_error(key + ": the vertex fit must bind truth_vtx_z, the selection and the centrality weights");
      const fs::path fit = fs::weakly_canonical(base / vertex.value("path", ""));
      if (!fs::is_regular_file(fit) || sha256File(fit.string()) != vertex.value("sha256", "")) throw std::runtime_error(key + ": vertex fit is missing or changed");
      vertex["path"] = fit.string();
      const fs::path evidence=fs::weakly_canonical(base/vertex.at("evidence").get<std::string>());
      vertex["evidence"]=evidence.string(); vertex["evidence_sha256"]=sha256File(evidence.string());
      vertex["data_trigger"]=sample.at("centrality_data_trigger");
      checkVertexFit(Json{{"key",key},{"weights",weights},{"vertex_reweighting",vertex}});
      for (const auto &src : sources)
        if (src.at("metadata").value("weighting", Json::object()).value("vertex_enabled", false) && !photonjet::declaredAuAuCopies(src.at("metadata").at("weighting"))) throw std::runtime_error(key + ": input already includes vertex weights");
    } else if (!weights.empty() || sample.contains("vertex_reweighting")) {
      throw std::runtime_error(key + ": centrality weights and vertex fits apply only to Au+Au simulation");
    }
    Json normalization = sample.value("normalization", Json::object());
    if (system=="PP" && family!="data") {
      if (normalization.value("source_factor","")!="family_reference")
        throw std::runtime_error(key+": p+p relative slice weights require family_reference normalization");
      const std::string reference=family=="photonjet" ? "photon_jet" : "inclusive_jet";
      if (normalization.value("reference",reference)!=reference) throw std::runtime_error(key+": normalization reference differs from its family");
      const Json producer=photonjet_config::readYaml(s.home+"/../TreeProduction/config/tree_production.yaml").at("ppg12_weights").at(reference);
      const double numerator=producer.at("reference_cross_section_pb");
      if (c.value("reference_cross_section_pb",Json::object()).contains(reference) && c.at("reference_cross_section_pb").at(reference)!=numerator)
        throw std::runtime_error(key+": reference cross section differs from the producer");
      long long generated=0;
      for (const auto &src : sources) generated+=src.at("metadata").at("entries").get<long long>();
      if (generated<=0) throw std::runtime_error(key+": no generated exposure in the selected trees");
      if (normalization.contains("generated_events") && normalization.at("generated_events").get<double>()!=generated)
        throw std::runtime_error(key+": handwritten generated_events differs from the exact selected tree entries");
      normalization["reference"]=reference; normalization["numerator_pb"]=numerator;
      normalization["generated_events"]=generated; normalization["denominator_events"]=generated;
      normalization["denominator"]= "all_selected_eventTree_entries";
      normalization["pb_per_generated_event"]=numerator/generated; normalization["applied_in_histograms"]=false;
      const std::string sliceName=producer.at("samples").at(name);
      const Json &slice=producer.at("slices").at(sliceName);
      normalization["generator_diagnostic"]={{"source_moment","raw_counts"},
        {"cross_section_pb",slice.at("cross_section_pb")},{"generated_events",generated},
        {"window_gev",slice.at("window_gev")},{"sample_slice",sliceName},
        {"photon_definition","maximum finite positive truthpho_pt with truthpho_source_role == 1"},
        {"jet_definition","maximum finite positive truthjet04_pt"},
        {"reco_comparison_vertex_abs_z_max_cm",30.0}};
    } else if (family != "data" && normalization.value("source_factor", "") == "auau_embedded_stitching") {
      /* the fill jobs count the events of the denominator and the merge sets the normalisation */
      normalization = stitchingNormalization(s.home, system, name, family);
    } else if (family != "data") {
      const double generated = normalization.value("generated_events", -1.0);
      double numerator = -1;
      if (normalization.value("source_factor", "") == "family_reference")
        numerator = c.value("reference_cross_section_pb", Json::object()).value(normalization.value("reference", ""), -1.0);
      else if (normalization.value("source_factor", "") == "cross_section_over_generated")
        numerator = normalization.value("cross_section_pb", -1.0);
      else throw std::runtime_error(key + " has an unsupported normalization rule");
      if (!(generated > 0) || !(numerator > 0)) throw std::runtime_error(key + " needs a positive cross section and generated-event count");
      normalization["numerator_pb"] = numerator;
      normalization["pb_per_generated_event"] = numerator / generated;
      normalization["applied_in_histograms"] = false;
    }
    const long population = sample.value("population_sources", -1L);
    if (population != -1 && population < static_cast<long>(selectedSources)) throw std::runtime_error(key + ": population_sources is smaller than the selected source count");
    const Json coverage = {{"selected_sources", selectedSources}, {"population_sources", population == -1 ? Json(nullptr) : Json(population)},
                           {"status", population == -1 ? "UNKNOWN" : population == static_cast<long>(selectedSources) ? "COMPLETE" : "PARTIAL"}};
    Json planned={{"key", key}, {"system", system}, {"name", name}, {"family", family}, {"sources", std::move(sources)}, {"coverage", coverage},
                           {"weights", weights}, {"vertex_reweighting", vertex}, {"normalization", normalization},
                           {"centrality_evidence", sample.value("centrality_evidence", Json())},
                           {"centrality_data_trigger", sample.value("centrality_data_trigger", Json())}};
    if (family=="data") {
      /* Fill each original once. Run files precede the combined analysis file;
         local indices change while source identity and metadata remain exact. */
      std::map<int,Json> runs;
      for (const auto &src : planned.at("sources")) {
        const int run=src.at("metadata").at("source").value("run",-1);
        if (run<=0) throw std::runtime_error(key+": DATA source has no valid run number");
        if (!runs.count(run)) runs[run]=Json::array();
        runs[run].push_back(src);
      }
      if (bounded) for (const auto &[run,path] : runInventories) runs[run]=Json();
      Json runTemplate=sampleHeader(planned); runTemplate["sources"]=Json::array();
      std::map<size_t,std::pair<std::string,size_t>> assignments;
      planned["sources"]=Json::array(); planned["run_packages"]=Json::array(); planned["run_coverages"]=Json::object(); planned["run_records"]=Json::object();
      const auto excluded=sample.value("exclude_runs",std::vector<int>{});
      planned["excluded_runs"]=excluded;
      size_t combinedSources=0;
      for (const auto &[run,inMemoryRows] : runs) {
        Json loaded; if (bounded) loaded=readSourceInventory(runInventories.at(run));
        const Json &rows=bounded ? loaded : inMemoryRows;
        Json perRun=runTemplate;
        const std::string runName="run_"+std::to_string(run), runKey=system+"/"+runName;
        if (!keys.insert(runKey).second) throw std::runtime_error("duplicate DATA run package " + runKey);
        perRun["name"]=runName; perRun["key"]=runKey; perRun["run"]=run; perRun["sources"]=rows;
        const auto runPop=sample.value("run_population_sources",Json::object());
        const long expectedRun=runPop.value(std::to_string(run),-1L);
        if (expectedRun!=-1 && expectedRun<static_cast<long>(rows.size())) throw std::runtime_error(runKey+": invalid run source population");
        perRun["coverage"]={{"selected_sources",rows.size()},{"population_sources",expectedRun==-1 ? Json() : Json(expectedRun)},
          {"status",expectedRun==-1 ? "UNKNOWN" : expectedRun==static_cast<long>(rows.size()) ? "COMPLETE" : "PARTIAL"}};
        perRun["data_run"]=dataRunRecord(rows,sample.value("run_luminosity",Json::object()).value(std::to_string(run),Json()),system);
        for (size_t i=0;i<rows.size();++i) {
          const size_t originalIndex=rows.at(i).at("index");
          perRun["sources"][i]["index"]=i;
          assignments[originalIndex]={runKey,i};
        }
        if (bounded) externalizeSample(perRun,common,sourceDirectory);
        samples.push_back(perRun);
        if (std::find(excluded.begin(),excluded.end(),run)!=excluded.end()) continue;
        planned["run_packages"].push_back(runKey);
        planned["run_coverages"][runKey]=perRun.at("coverage");
        planned["run_records"][runKey]=perRun.at("data_run");
        if (bounded) planned["_source_files"].push_back(perRun.at("_source_file"));
        else for (const auto &src : perRun.at("sources")) planned["sources"].push_back(src);
        combinedSources+=rows.size();
      }
      for (auto &job : jobs) if (job.at("sample")==key) {
        const auto &assigned=assignments.at(job.at("source").get<size_t>());
        job["sample"]=assigned.first; job["source"]=assigned.second;
      }
      if (planned.at("run_packages").empty()) throw std::runtime_error(key+": no analysis runs selected");
      planned["coverage"]["selected_sources"]=combinedSources;
      if (!excluded.empty()) planned["coverage"]["status"]="PARTIAL";
      if (bounded) {
        planned["_source_count"]=combinedSources;
        planned["_package_binding_sha256"]=packageBindingFromSources(common,planned,[&](std::ostream &out) {
          out<<'['; size_t count=0;
          for (const auto &pin : planned.at("_source_files")) {
            std::ifstream in(pin.at("path").get<std::string>()); std::string line;
            while (std::getline(in,line)) {
              const Json src=Json::parse(line); if (count++) out<<',';
              out<<Json{{"identity",src.at("identity")},{"sha256",src.at("expected_sha256")},{"metadata",src.at("metadata")}}.dump(-1,' ',true);
            }
            if (!in.eof()) throw std::runtime_error("cannot read planned DATA inventory");
          }
          if (count!=combinedSources) throw std::runtime_error("combined DATA inventory count differs");
          out<<']';
        });
        planned.erase("sources");
        for (const auto &[run,path] : runInventories) fs::remove(path);
      }
    }
    if (bounded && family!="data") externalizeSample(planned,common,sourceDirectory);
    samples.push_back(planned);
  }
  inspected.finish();
  if (samples.empty()) throw std::runtime_error("no catalogued sample matches");
  const Json histograms = histogramSettings(s);
  const auto responseVariation = photonjet::response::parseVariation(s.responseVariation.c_str());
  photonjet::response::requireRunnableVariation(responseVariation,!s.emsSelectionContract.is_null());
  for (const auto& sample : samples) if (sample.at("family") != "data")
    photonjet::response::validateVariationSystem(responseVariation,
        sample.at("system") == "PP" ? photonjet::System::pp : photonjet::System::auau);
  Json out = {{"schema", "TreeToHistsPlanV1"}, {"settings", s.path}, {"catalogue", fs::absolute(catalogue).string()},
              {"catalogue_sha256", sha256File(catalogue)}, {"selection", s.selection}, {"histograms", histograms},
              {"root_version", gROOT->GetVersion()}, {"producer_code", reductionCode()}, {"max_events_per_job", maxEvents},
              {"samples", std::move(samples)}, {"jobs", std::move(jobs)}};
  if (bounded) out["source_inventory_layout"]="TreeToHistsSourceInventoryV1";
  if(!s.responseGrid.is_null()) {
    photonjet::response::native_grid::validateSettings(s.responseGrid);
    for(const auto &sample:out.at("samples"))if(sample.at("family")!="data")for(const auto &variant:s.responseGrid.at("variants"))
      photonjet::response::validateVariationSystem(photonjet::response::parseVariation(variant.at("variation").get<std::string>().c_str()),sample.at("system")=="PP"?photonjet::System::pp:photonjet::System::auau);
    out["response_component_grid"]=s.responseGrid;
  }
  out["plan_sha256"] = digest(out);
  return out;
}

const Json &sampleOf(const Json &p, const std::string &key) {
  for (const auto &s : p.at("samples")) if (s.at("key").get<std::string>() == key) return s;
  throw std::runtime_error("no sample " + key + " in the plan");
}

std::string listed(const Json &values) {
  std::ostringstream out;
  out.precision(17);
  for (const auto &v : values) out << v.get<double>() << ' ';
  return out.str();
}

/* One job: every filler on one eventTree file. Products are written under .partial names and moved into place, and
   done.json (with the input's sha256) is written last, so a resumed plan redoes exactly the incomplete jobs. */
}  /* namespace */

/* The canonical Au+Au embedded stitching of a sample (config/auau_embedded_stitching.json): its ownership window and
   denominator kind, its cross section and the table's sha256; pb_per_generated_event is set from the counted
   denominator. */
Json stitchingNormalization(const std::string &home, const std::string &system, const std::string &name, const std::string &family) {
  const std::string table = home + "/config/auau_embedded_stitching.json";
  const Json entry = readJson(table).at("samples").value(name, Json());
  if (system != "AuAu" || entry.is_null() || entry.at("family").get<std::string>() != family)
    throw std::runtime_error(system + "/" + name + ": no Au+Au embedded " + family + " sample of that name in " + table);
  return Json{{"source_factor", "auau_embedded_stitching"},
              {"ownership", {{"owner", entry.at("owner")}, {"window_gev", entry.at("window_gev")}, {"denominator", entry.at("denominator")}}},
              {"cross_section_pb", entry.at("cross_section_pb")}, {"stitching_table_sha256", sha256File(table)},
              {"pb_per_generated_event", nullptr}, {"applied_in_histograms", false}};
}

/* The ownership of a sample's normalisation record (none without one). */
photonjet::Ownership ownershipOf(const Json &normalization) {
  photonjet::Ownership o;
  if (!normalization.contains("ownership")) return o;
  const Json &own = normalization.at("ownership");
  const std::string owner = own.at("owner");
  if (owner == "leading_truth_jet_r04_pt") o.owner = photonjet::Ownership::Owner::leadingTruthJet;
  else if (owner == "generator_filter_photon_pt") o.owner = photonjet::Ownership::Owner::filterPhoton;
  else throw std::runtime_error("unknown ownership variable " + owner);
  o.low = own.at("window_gev")[0].get<double>();
  if (!own.at("window_gev")[1].is_null()) o.high = own.at("window_gev")[1].get<double>();
  return o;
}

/* The events of one eventTree file in a sample's stitching denominator: every entry (all_events), or the entries whose
   owner variable lies in the window (events_in_window); entries without the owner variable are counted as missing. */
Json stitchingDenominator(const std::string &path, const Json &normalization) {
  const photonjet::Ownership o = ownershipOf(normalization);
  const std::string kind = normalization.at("ownership").at("denominator");
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  TTreeReader reader("eventTree", f.get());
  photonjet::applyEntryRange(reader);
  photonjet::OwnershipGate ownership(reader, o, true);
  long long events = 0, missing = 0;
  while (reader.Next()) {
    photonjet::reportProgress(reader.GetCurrentEntry()+1);
    const double x = ownership.value();
    if (!std::isfinite(x) && !(kind == "all_events" && ownership.rejectedBeforeOwnership())) ++missing;
    if (kind == "all_events" || o.contains(x)) ++events;
  }
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
  return Json{{"kind", kind}, {"events", events}, {"missing", missing}};
}

namespace {

/* Execute one source job from the sealed plan. Each product kind sees the
   same resolved selection, source ownership and event-weight context.
   done.json is reusable only while those bindings and input bytes match. */
void fill(const fs::path &work, int index, const char* jesUp=nullptr,const char* jesDown=nullptr) {
  // Preflight hashes precede this row's counters. They must not pulse a prior row.
  photonjet::progressPath.clear();
  const fs::path descriptor = work / "rows" / std::to_string(index / 1000) / (std::to_string(index) + ".json");
  Json p = readJson((fs::exists(descriptor) ? descriptor : work / "plan.json").string());
  if (fs::exists(descriptor)) {
    const std::string claimed = p.at("descriptor_sha256");
    p.erase("descriptor_sha256");
    if (digest(p) != claimed || readJson((work / "binding.json").string()).at("plan_sha256") != p.at("plan_sha256"))
      throw std::runtime_error("fill descriptor differs from its immutable plan binding");
  } else {
    Json body = p; body.erase("plan_sha256");
    if (digest(body) != p.at("plan_sha256").get<std::string>()) throw std::runtime_error("plan digest differs");
  }
  checkReductionCode(p);
  if (p.at("root_version").get<std::string>() != gROOT->GetVersion()) throw std::runtime_error("ROOT version differs from the plan");
  const Json &job = p.at("jobs").at(fs::exists(descriptor) ? 0 : index);
  const Json &sample = sampleOf(p, job.at("sample"));
  checkVertexFit(sample);
  const Json ppPeriods = ppPeriodsOf(p, sample);
  checkPpPeriodFits(ppPeriods);
  const Json &source = sample.at("sources").at(job.at("source").get<size_t>());
  const fs::path dir = jobDirectory(work.string(), job);
  if(p.contains("selection_fanout") && !fs::exists(dir/"done.json"))throw std::runtime_error("selection view requires fill-selection-group to preserve parent/source composition and gap witnesses");
  fs::create_directories(dir);
  photonjet::entryBegin = job.value("range_begin", 0LL);
  photonjet::entryEnd = job.value("range_end", source.at("metadata").at("entries").get<long long>());
  photonjet::progressRowTotal = photonjet::entryEnd - photonjet::entryBegin;
  std::vector<std::string> kinds = {"spectra", "qa", "extensions"};
  if (sample.at("family") != "data") kinds.push_back("response");
  const bool countDenominator = sample.at("normalization").contains("ownership");
  photonjet::progressTotal = photonjet::progressRowTotal * (kinds.size() + (countDenominator ? 1 : 0));
  photonjet::progressOffset = 0;
  photonjet::progressPath = (dir / "progress.json").string();
  photonjet::progressIdentity = {{"campaign", p.at("plan_sha256")}, {"row", index}, {"sample", job.at("sample")},
    {"source", source.at("identity")}, {"entry_begin", photonjet::entryBegin}, {"entry_end", photonjet::entryEnd}};
  photonjet::progressStart = std::chrono::steady_clock::now();
  photonjet::progressPhase = "checking_input";
  photonjet::reportProgress(photonjet::entryBegin);
  /* The tree file is not hashed: its embedded identity and its metadata are compared with the plan below, and its
     size and modification time must be the same after the reduction as before. `sha` is the planned label. */
  const std::string input = job.at("input"), sha = source.at("expected_sha256");
  const auto stamp = [&]() { return std::make_pair(fs::file_size(input), fs::last_write_time(input)); };
  const auto before = stamp();
  checkNativeModelFiles(p,sample.at("system"));
  std::string fileFamily,fileSystem;
  const Json inspected=sourceInfo(input,p.at("selection"),fileFamily,fileSystem);
  if (inspected.at("identity")!=source.at("identity") || inspected.at("metadata")!=source.at("metadata"))
    throw std::runtime_error("input metadata differs from the planned source: " + input);
  const std::string codeSha = digest(p.at("producer_code"));
  if (fs::exists(dir / "done.json")) {
    const Json saved = readJson((dir / "done.json").string());
    if (saved.at("plan_sha256") != p.at("plan_sha256") || saved.at("sha256") != sha ||
        saved.value("producer_code_sha256", "") != codeSha || saved.at("range_begin") != photonjet::entryBegin ||
        saved.at("range_end") != photonjet::entryEnd || saved.at("root_version").get<std::string>() != gROOT->GetVersion())
      throw std::runtime_error("completed fill differs from its input, range, runtime or plan");
    for(const auto &kind:kinds)if(!saved.at("products").contains(kind))throw std::runtime_error("completed fill lacks a required product");
    if(p.contains("response_component_grid") && sample.at("family")!="data")for(const auto &v:p.at("response_component_grid").at("variants"))
      if(v.at("name")!="nominal" && !saved.at("products").contains("response__"+v.at("name").get<std::string>()))throw std::runtime_error("completed fill lacks a response component");
    for (const auto &[kind,hash] : saved.at("products").items())
      if (hash.get<std::string>() != sha256File((dir / (kind + ".root")).string()))
        throw std::runtime_error("completed product hash differs");
    photonjet::progressOffset = photonjet::progressTotal - photonjet::progressRowTotal;
    photonjet::progressPhase = "reused"; photonjet::reportProgress(photonjet::entryEnd);
    return;
  }
  const auto start = std::chrono::steady_clock::now();
  const std::string selection = (work / "selection.json").string(), weights = listed(sample.at("weights"));
  photonjet::configuredXjEdges = p.at("histograms").at("native_xj_edges").get<std::vector<double>>();
  photonjet::configureJetRadius(p.at("histograms").at("jet_radius").get<double>());
  photonjet::configuredDataJetScale = p.at("histograms").at("data_jet_scale").at(sample.at("system") == "PP" ? "pp" : "auau").get<double>();
  photonjet::configuredMeasurementEtEdges=p.at("histograms").at("native_measurement_et_edges_gev").get<std::vector<double>>();
  photonjet::candidate_cache::enabled=p.at("histograms").value("candidate_records",false);
  photonjet::candidate_cache::maxRows=p.at("histograms").value("candidate_records_max_rows",2000000L);
  photonjet::candidate_cache::minRecoEt=p.at("histograms").value("candidate_records_min_et",0.);
  photonjet::strictReduction = true;
  photonjet::configureVertexFit(sample.at("vertex_reweighting").is_null() ? Json::object() : sample.at("vertex_reweighting"));
  photonjet::auauReceiptSha=sample.at("vertex_reweighting").is_null() ? "" : sample.at("vertex_reweighting").at("evidence_sha256").get<std::string>();
  photonjet::auauStitchingSha=sample.at("normalization").value("stitching_table_sha256", "");
  photonjet::configurePpPeriods(ppPeriods.is_null() ? Json::object() : ppPeriods);
  photonjet::configuredOwnership = ownershipOf(sample.at("normalization"));
  if (readJson(selection) != p.at("selection")) throw std::runtime_error("resolved selection differs from plan");
  const Json &h = p.at("histograms");
  const std::string emsContractText=h.contains("ems_selection_contract")?h.at("ems_selection_contract").dump():"";
  const std::string variationText=h.at("response_variation").get<std::string>();
  Json products = Json::object();
  Json gridReceipt;
  for (const auto &kind : kinds) {
    const fs::path out = dir / (kind + ".root"), partial = dir / (kind + ".root.partial");
    const fs::path checkpoint = dir / (kind + ".json");
    Json binding = {{"plan_sha256", p.at("plan_sha256")}, {"input_sha256", sha}, {"producer_code_sha256", codeSha},
      {"range_begin", photonjet::entryBegin}, {"range_end", photonjet::entryEnd}, {"root_version", gROOT->GetVersion()}};
    photonjet::progressPhase = kind; photonjet::reportProgress(photonjet::entryBegin);
    if (fs::exists(checkpoint)) {
      Json saved = readJson(checkpoint.string()); const std::string expected = saved.at("sha256"); saved.erase("sha256");
      if (saved != binding || expected != sha256File(out.string())) throw std::runtime_error("fill checkpoint differs: " + checkpoint.string());
      products[kind] = expected;
    } else {
      fs::remove(partial);
      int code = 0;
      if (kind == "spectra") code = FillPhotonJet(input.c_str(), partial.c_str(), "nominal", selection.c_str(), weights.c_str(),variationText.c_str(),emsContractText.c_str(),sha256File,jesUp,jesDown);
      else if (kind == "qa") code = FillPhotonQA(input.c_str(), partial.c_str(), "nominal", weights.c_str(), selection.c_str(),variationText.c_str(),emsContractText.c_str(),sha256File);
      else if (kind == "extensions") code = FillExtensions(input.c_str(), partial.c_str(), selection.c_str(), weights.c_str(), h.at("require_stage3_fractions").get<bool>(), photonjet::response::parseVariation(variationText.c_str()).centralityOff);
      else if(p.contains("response_component_grid")) {
        namespace ng=photonjet::response::native_grid;
        const auto &config=p.at("response_component_grid");ng::validateSettings(config);
        const auto system=sample.at("system")=="PP"?photonjet::System::pp:photonjet::System::auau;
        const auto profile=photonjet::configuredProfile("nominal",system,selection.c_str());
        photonjet::response::ResponseGrids grids;grids.reco=h.at("native_reco_response_edges_gev").get<std::vector<double>>();
        grids.truth=h.at("native_truth_response_edges_gev").get<std::vector<double>>();grids.xj=h.at("native_xj_edges").get<std::vector<double>>();
        const auto cw=photonjet::centralityWeights(weights,system);
        Json spec={{"schema","PhotonJetNativeResponseGridV1"},{"identity",{{"task","maintained-campaign"},{"workstream","response-components"},{"campaign",p.at("plan_sha256")},{"row",std::to_string(index)}}},
          {"sources",Json::array({{{"path",input},{"sha256",sha},{"input_file_sha256",source.at("identity")},
            {"entries",source.at("metadata").at("entries")},{"entry_begin",photonjet::entryBegin},{"entry_end",photonjet::entryEnd}}})},
          {"variants",config.at("variants")},{"producer_code",p.at("producer_code")},{"runtime_context",ng::context(profile,grids,cw)},
          {"limits",config.at("limits")},{"output_directory",fs::absolute(dir/"response_grid").lexically_normal().string()}};
        ng::copyEmsContract(config,spec);
        for(auto &v:spec["variants"]){v["plan"]=ng::resolvedPlan(spec,v);v["plan_sha256"]=digest(v["plan"]);}
        const auto previousObserver=photonjet::progressObserver;
        photonjet::progressObserver=[&](const Json &value) {
          if(!previousObserver)return;
          Json mapped=value;
          if(value.value("schema","")=="PhotonJetNativeResponseGridProgressV1") {
            mapped["processed"]=photonjet::progressOffset+value.at("processed").get<Long64_t>();
            mapped["total"]=photonjet::progressTotal;mapped["remaining"]=photonjet::progressTotal-mapped.at("processed").get<Long64_t>();
          }
          previousObserver(mapped);
        };
        try {gridReceipt=ng::run(spec,profile,grids,cw);}catch(...){photonjet::progressObserver=previousObserver;throw;}
        photonjet::progressObserver=previousObserver;
        if(gridReceipt.at("original_entries").get<Long64_t>()!=photonjet::progressRowTotal)throw std::runtime_error("native grid did not finish the public row");
        photonjet::reportProgress(photonjet::entryEnd);
        for(const auto &product:gridReceipt.at("products")) {
          const std::string name=product.at("name"),key=name=="nominal"?"response":"response__"+name;
          const auto from=dir/"response_grid"/product.at("path").get<std::string>();
          const auto to=name=="nominal"?partial:dir/(key+".root");fs::create_hard_link(from,to);
          if(name!="nominal")products[key]=product.at("sha256");
        }
      }
      else code = FillPhotonJetResponse(input.c_str(), partial.c_str(), "nominal", listed(h.at("native_reco_response_edges_gev")).c_str(),
        listed(h.at("native_truth_response_edges_gev")).c_str(), listed(h.at("native_xj_edges")).c_str(),
        weights.c_str(), variationText.c_str(), selection.c_str(),emsContractText.c_str(),sha256File);
      if (code) throw std::runtime_error("the " + kind + " filler failed on " + input);
      if (photonjet::progressRowProcessed != photonjet::progressRowTotal) throw std::runtime_error("filler stopped before the end of its range");
      fs::rename(partial, out);
      products[kind] = sha256File(out.string()); binding["sha256"] = products.at(kind);
      writeText(checkpoint, binding.dump() + "\n");
    }
    photonjet::reportProgress(photonjet::entryEnd);
    photonjet::progressOffset += photonjet::progressRowTotal;
  }
  if(p.contains("response_component_grid") && sample.at("family")!="data") {
    if(gridReceipt.is_null())gridReceipt=readJson((dir/"response_grid/receipt.json").string());
    if(gridReceipt.at("status")!="PASS_RESPONSE_COMPONENT_ONLY")throw std::runtime_error("response grid did not publish a PASS component receipt");
    for(const auto &product:gridReceipt.at("products")) {
      const std::string name=product.at("name"),key=name=="nominal"?"response":"response__"+name;
      if(sha256File((dir/(key+".root")).string())!=product.at("sha256").get<std::string>())throw std::runtime_error("response component differs from its receipt");
      products[key]=product.at("sha256");
    }
  }
  Json done = {{"input", input}, {"sha256", sha}, {"root_version", gROOT->GetVersion()}, {"plan_sha256", p.at("plan_sha256")},
    {"producer_code_sha256", codeSha}, {"range_begin", photonjet::entryBegin}, {"range_end", photonjet::entryEnd}, {"products", products}};
  if(!gridReceipt.is_null())done["response_grid"]=gridReceipt;
  if (countDenominator) {
    photonjet::progressPhase = "stitching_denominator"; photonjet::reportProgress(photonjet::entryBegin);
    done["stitching_denominator"] = stitchingDenominator(input, sample.at("normalization"));
    photonjet::reportProgress(photonjet::entryEnd);
  }
  // Hash heartbeats include progressRowProcessed: retain the last pass's start,
  // rather than the next-pass offset left by the product loop.
  photonjet::progressOffset = photonjet::progressTotal - photonjet::progressRowTotal;
  photonjet::progressPhase = "validating_input_and_models";
  photonjet::reportProgress(photonjet::entryEnd);
  if (stamp() != before) throw std::runtime_error("input changed while it was being reduced: " + input);
  checkNativeModelFiles(p,sample.at("system"));
  checkVertexFit(sample); checkPpPeriodFits(ppPeriods);
  done["seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  done["processed_entries"] = photonjet::progressRowTotal;
  done["sample"] = sample.at("key");
  done["sample_class"] = sample.at("system").get<std::string>() + "/" + sample.at("family").get<std::string>();
  done["reduction_settings_sha256"] = digest(Json{{"selection",p.at("selection")},{"histograms",p.at("histograms")}});
  writeText(dir / "done.json", done.dump(1) + "\n");
  photonjet::progressOffset = photonjet::progressTotal - photonjet::progressRowTotal;
  photonjet::progressPhase = "complete"; photonjet::reportProgress(photonjet::entryEnd);
}

/* Three independently sealed ordinary plans keep each DATA product's scale,
   settings and package identity truthful. Only the science traversal is fused;
   QA/extensions are explicitly nominal and ordinary merges need no new path. */
void fillDataJes(const fs::path& nominal,const fs::path& up,const fs::path& down,int index,bool muteReusedProgress=false) {
  const std::array<fs::path,3> works={fs::canonical(nominal),fs::canonical(up),fs::canonical(down)};
  if(works[0]==works[1] || works[0]==works[2] || works[1]==works[2])throw std::runtime_error("DATA JES work identities must be distinct");
  std::array<Json,3> plans;std::array<fs::path,3> dirs;
  const std::array<double,3> scales={1./.916,1./.907,1./.923};
  Json baseline;
  for(size_t i=0;i<3;++i){
    const auto descriptor=works[i]/"rows"/std::to_string(index/1000)/(std::to_string(index)+".json");
    const bool scoped=muteReusedProgress && fs::exists(descriptor);
    plans[i]=readJson((scoped?descriptor:works[i]/"plan.json").string());auto compare=plans[i];
    if(scoped){
      const auto claimed=compare.at("descriptor_sha256");compare.erase("descriptor_sha256");
      if(digest(compare)!=claimed.get<std::string>() || compare.at("plan_sha256")!=readJson((works[i]/"binding.json").string()).at("plan_sha256"))
        throw std::runtime_error("DATA JES row descriptor binding differs");
      compare.erase("plan_sha256");
    }else{
      const auto claimed=compare.at("plan_sha256");compare.erase("plan_sha256");
      if(digest(compare)!=claimed.get<std::string>())throw std::runtime_error("DATA JES plan digest differs");
    }
    checkReductionCode(plans[i]);
    const auto& h=plans[i].at("histograms");
    if(h.at("data_jet_scale").at("pp")!=scales[i] || h.at("data_jet_scale").at("auau")!=1. || h.at("response_variation")!="")
      throw std::runtime_error("DATA JES requires exact absolute scales and nominal other-system/response settings");
    for(const auto& sample:plans[i].at("samples")){
      if(sample.at("system")!="PP" || sample.at("family")!="data" || sample.contains("reused_package") || sample.contains("_source_files"))
        throw std::runtime_error("DATA JES requires ordinary original-only pp DATA plans");
      std::set<std::string> sources;
      for(const auto& source:sample.at("sources"))if(!sources.insert(source.at("identity").get<std::string>()).second)
        throw std::runtime_error("duplicate original source identity in DATA JES plan");
    }
    compare.erase("settings");compare["histograms"]["data_jet_scale"]["pp"]=scales[0];
    if(i==0)baseline=compare;else if(compare!=baseline)throw std::runtime_error("DATA JES plans differ beyond their exact pp DATA scale");
    const auto& job=plans[i].at("jobs").at(scoped?0:index);dirs[i]=jobDirectory(works[i].string(),job);
    for(const char* name:{"done.json","spectra.json","spectra.root","spectra.root.partial","qa.root","extensions.root"})
      if(fs::exists(dirs[i]/name))throw std::runtime_error("DATA JES refuses existing row products; reuse or reconcile the exact row first");
  }
  for(const auto& dir:dirs)fs::create_directories(dir);
  const auto upPartial=(dirs[1]/"spectra.root.partial").string(),downPartial=(dirs[2]/"spectra.root.partial").string();
  fill(works[0],index,upPartial.c_str(),downPartial.c_str());
  const auto nominalDone=readJson((dirs[0]/"done.json").string());
  for(size_t i=1;i<3;++i){
    fs::rename(dirs[i]/"spectra.root.partial",dirs[i]/"spectra.root");
    for(const char* kind:{"spectra","qa","extensions"}){
      const std::string name=kind;
      if(name!="spectra")fs::create_hard_link(dirs[0]/(name+".root"),dirs[i]/(name+".root"));
      Json binding=readJson((dirs[0]/(name+".json")).string());
      binding["plan_sha256"]=plans[i].at("plan_sha256");binding["sha256"]=sha256File((dirs[i]/(name+".root")).string());
      writeText(dirs[i]/(name+".json"),binding.dump()+"\n");
    }
    const auto observer=photonjet::progressObserver;
    if(muteReusedProgress)photonjet::progressObserver={};
    try{fill(works[i],index);}catch(...){photonjet::progressObserver=observer;throw;}
    photonjet::progressObserver=observer;
    auto done=readJson((dirs[i]/"done.json").string());
    done["data_jes_fusion"]={{"nominal_plan_sha256",plans[0].at("plan_sha256")},{"science_traversals",1},{"qa_extensions","nominal_reused"}};
    writeText(dirs[i]/"done.json",done.dump(1)+"\n");
  }
}

std::string option(int argc, char **argv, const std::string &name, bool requiredOption = true) {
  for (int i = 2; i + 1 < argc; ++i) if (name == argv[i]) return argv[i + 1];
  if (requiredOption) throw std::runtime_error("missing " + name);
  return "";
}

}  /* namespace */

std::string jobDirectory(const std::string &work, const Json &job) {
  const long original = job.value("original_source", job.at("source").get<long>());
  std::string path = productDirectory(work, job.at("sample"), original);
  if (job.contains("range_begin")) path += "/" + std::to_string(job.at("range_begin").get<long long>()) + "_" + std::to_string(job.at("range_end").get<long long>());
  return path;
}
Json preparePlan(const std::string &settings, const std::string &inputs, const std::string &home,
                 const std::vector<std::string> &samples, long long maxEvents, const std::string &work,bool bounded,const Json &inspections) {
  if (bounded && work.empty()) throw std::runtime_error("bounded planning requires a work directory");
  if (!work.empty()) {fs::create_directories(work); photonjet::progressPath=(fs::path(work)/"plan_progress.json").string();}
  Json p=plan(Settings::read(settings, home), inputs, samples, maxEvents,bounded ? fs::path(work)/"source_inventory" : fs::path{},inspections);
  checkReductionCode(p);
  const Json collection=readInputCollection(inputs);
  if (collection.contains("reuse_collection")) {
    reusePackages(p,fs::canonical(fs::absolute(inputs).parent_path()/collection.at("reuse_collection").get<std::string>()).string());
    p.erase("plan_sha256"); p["plan_sha256"]=digest(p);
  }
  if (!work.empty()) {
    photonjet::progressPhase="planned"; photonjet::progressIdentity["campaign"]=p.at("plan_sha256");
    photonjet::reportProgress(photonjet::progressRowTotal); photonjet::progressPath.clear();
  }
  return p;
}
Json settingsSnapshot(const std::string &settings, const std::string &home) {
  const auto s=Settings::read(settings,home);
  Json result={{"selection",s.selection},{"histograms",histogramSettings(s)}};
  if(!s.responseGrid.is_null())result["response_component_grid"]=s.responseGrid;
  return result;
}
Json inspectSources(const Json &descriptor,const Json &snapshot) {
  const auto &files=descriptor.at("files");
  if (files.empty() || files.size()>10000) throw std::runtime_error("invalid bounded inspection size");
  Json rows=Json::array();ReductionProgress progress("inspect",files.size());
  for (const auto &file : files) {
    std::string family,system;
    Json info=sourceInfo(file.at("path").get<std::string>(),snapshot.at("selection"),family,system);
    rows.push_back({{"input",file},{"info",std::move(info)},{"family",family},{"system",system}});
    if (rows.size()<files.size()) progress.advance(rows.size());
  }
  progress.finish();return rows;
}
Json freezeInputs(const std::string &inputs) {
  Json c=readInputCollection(inputs);
  const fs::path base=fs::absolute(inputs).parent_path();
  c["configuration_sha256"]=sha256File(inputs);
  c["publication_bindings"]=Json::array();
  if (c.contains("reuse_collection")) {
    c["reuse_collection"]=fs::canonical(base/c.at("reuse_collection").get<std::string>()).string();
    c["publication_bindings"].push_back({{"path",c.at("reuse_collection")},{"sha256",sha256File(c.at("reuse_collection"))}});
  }
  for (auto &sample : c.at("samples")) {
    Json files=Json::array();
    for (const auto &f : sample.value("files",Json::array())) {
      if (!f.is_object() || !f.contains("sha256"))
        throw std::runtime_error("campaign inputs need published groups or files with an explicit sha256; use plan on a compute node to inspect other files");
      // Freeze names and hashes without serial filesystem lookups. Parallel
      // header inspection checks readability; each fill checks the ROOT bytes.
      files.push_back({{"path",(base/f.at("path").get<std::string>()).string()},{"sha256",f.at("sha256")}});
    }
    for (const auto &g : sample.value("groups",Json::array())) {
      const fs::path group=fs::weakly_canonical(base/g.get<std::string>()), current=fs::canonical(group/"CURRENT");
      for (const auto &[path,sha] : groupFiles(group)) files.push_back({{"path",path.string()},{"sha256",sha}});
      c["publication_bindings"].push_back({{"path",(current/"MANIFEST.json").string()}, {"sha256",sha256File((current/"MANIFEST.json").string())}});
    }
    if (files.empty()) throw std::runtime_error("empty campaign input sample");
    sample["files"]=files; sample.erase("groups");
  }
  return c;
}
void savePlanRow(const Json &p, long index, const std::string &directory,const Json &providedSource) {
  const fs::path work=directory;
  Json job=p.at("jobs").at(index), sample;
  const Json &whole=sampleOf(p,job.at("sample"));
  const Json source=providedSource.is_null() ? whole.at("sources").at(job.at("source").get<size_t>()) : providedSource;
  for (const auto &[key,value] : whole.items()) if (key!="sources") sample[key]=value;
  for (const char *key : {"_source_count","_package_binding_sha256","_source_file","_source_files"}) sample.erase(key);
  sample["sources"]=Json::array({source}); job["original_source"]=job.at("source"); job["source"]=0;
  Json row={{"schema",p.at("schema")},{"plan_sha256",p.at("plan_sha256")},{"root_version",p.at("root_version")},
    {"selection",p.at("selection")},{"histograms",p.at("histograms")},{"producer_code",p.at("producer_code")},
    {"samples",Json::array({sample})},{"jobs",Json::array({job})}};
  if(p.contains("response_component_grid"))row["response_component_grid"]=p.at("response_component_grid");
  if(p.contains("selection_fanout"))row["selection_fanout"]=p.at("selection_fanout");
  row["descriptor_sha256"]=digest(row);
  const fs::path path=work/"rows"/std::to_string(index/1000)/(std::to_string(index)+".json");
  fs::create_directories(path.parent_path()); writeText(path,row.dump()+"\n");
}
void savePlan(const Json &p, const std::string &directory, bool writeRows) {
  const fs::path work = directory; fs::create_directories(work);
  if (fs::exists(work / "plan.json") && readJson((work / "plan.json").string()).at("plan_sha256") != p.at("plan_sha256"))
    throw std::runtime_error(work.string() + " holds another plan; choose a new work directory");
  writeText(work / "selection.json", p.at("selection").dump() + "\n");
  const fs::path next=work/"plan.json.partial";
  std::ofstream stream(next); stream<<std::setw(1)<<p<<'\n'; stream.close();
  if (!stream) throw std::runtime_error("cannot write plan"); fs::rename(next,work/"plan.json");
  std::ostringstream sampleList;
  for (const auto &sample : p.at("samples")) sampleList<<sample.at("key").get<std::string>()<<'\n';
  writeText(work/"samples.txt",sampleList.str());
  writeText(work / "binding.json", Json{{"plan_sha256",p.at("plan_sha256")}}.dump() + "\n");
  std::ostringstream list;
  for (size_t i = 0; i < p.at("jobs").size(); ++i) {
    const Json &job=p.at("jobs")[i];
    if (writeRows) savePlanRow(p,i,directory);
    list << i << " " << job.at("sample").get<std::string>() << " " << job.at("input").get<std::string>() << "\n";
  }
  writeText(work / "jobs.txt", list.str());
}
void fillRow(const std::string &work, long index) { fill(work, index); }


/* A farm row reduces a bounded, ordered batch in private worker scratch, then
   publishes one intermediate package. Existing fills and source-aware merges
   retain each original's contribution and event identity. Scratch is removed
   only after the chunk and completion receipt have both been validated.
   A group may be run by more than one job (a repeat of a slow one): the package
   takes its name by a hard link, which only one job can make, so the first
   to finish publishes and every other one keeps that result and stops. */
struct DataJesGroup {
  std::array<fs::path,3> homes,outputs;
  std::array<Json,3> plans;
  std::array<std::string,3> descriptorHashes;
};
void fillGroupImpl(const std::string &work, const std::string &group, const std::string &scratch,
                   const DataJesGroup* jes=nullptr) {
  if (!std::regex_match(group,std::regex("G[0-9]+"))) throw std::runtime_error("invalid fill group");
  const std::string shard=std::to_string(std::stoul(group.substr(1))/1000);
  const fs::path home=work, descriptor=home/"fill_groups"/shard/(group+".json"), output=home/"chunks"/shard/group;
  Json d=readJson(descriptor.string()); const std::string descriptorSha=d.at("descriptor_sha256");
  d.erase("descriptor_sha256");
  if (digest(d)!=descriptorSha || d.at("plan").at("plan_sha256")!=readJson((home/"binding.json").string()).at("plan_sha256"))
    throw std::runtime_error("fill group differs from its immutable plan");
  const Json &p=d.at("plan"), &sample=p.at("samples").at(0);
  checkReductionCode(p); checkVertexFit(sample); checkPpPeriodFits(ppPeriodsOf(p,sample));
  if (p.at("root_version").get<std::string>()!=gROOT->GetVersion()) throw std::runtime_error("ROOT version differs from the plan");
  long long entries=0;
  for (const auto &src : sample.at("sources")) entries+=src.at("metadata").at("entries").get<long long>();
  const long passes=sample.at("family")=="data" ? 3 : 4+(sample.at("normalization").contains("ownership") ? 1 : 0);
  const long long total=entries*passes;
  fs::create_directories(output);
  const std::string package=(output/(sample.at("key").get<std::string>()+".root")).string();
  const auto started=std::chrono::steady_clock::now();
  long long completed=0,highestProcessed=0;
  /* While the group runs its progress is written to the job's own scratch (live), every few seconds. The shared
     output folder gets one record, at the end: thousands of jobs rewriting a file there every two seconds is the
     small-file load a shared file system must not carry. */
  fs::path live;
  auto progress=[&](const std::string &phase,const Json &child=Json::object(),const Json &error=nullptr) {
    long long processed=completed+child.value("processed",0LL);
    if(jes){highestProcessed=std::max(highestProcessed,processed);processed=highestProcessed;}
    if (processed<0 || processed>total) throw std::runtime_error("group progress exceeds its sealed event total");
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    const double rate=seconds>0 ? processed/seconds : 0;
    Json record={{"campaign",p.at("plan_sha256")},{"row",group},{"sample",sample.at("key")},{"source",group},
      {"phase",phase},{"unit","event_passes"},{"processed",processed},{"total",total},{"remaining",total-processed},
      {"row_processed",processed},{"row_total",total},{"row_remaining",total-processed},
      {"deterministic_skip_processed",0},{"deterministic_skip_total",0},{"deterministic_skip_remaining",0},
      {"source_progress",child},{"input_entries",entries},
      {"elapsed_seconds",seconds},{"rate",rate},{"eta_seconds",rate>0 ? (total-processed)/rate : -1},
      {"updated_at_unix",std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},
      {"latest_error",error}};
    if(jes)record["data_jes_fusion"]={{"views",3},{"physical_passes_per_source",3},{"physical_work_owner",p.at("plan_sha256")}};
    if (phase=="complete" || phase=="failed") writeText(output/"progress.json",record.dump()+"\n");
    else if (!live.empty()) writeText(live,record.dump()+"\n");
  };
  auto finish=[&](const std::string &expected,const Json &seconds,const std::string &state) {
    writeText(output/"done.json",Json{{"plan_sha256",p.at("plan_sha256")},{"descriptor_sha256",descriptorSha},
      {"sha256",expected},{"status",state},{"processed_entries",entries},{"source_count",sample.at("sources").size()},
      {"sample_class",sample.at("system").get<std::string>()+"/"+sample.at("family").get<std::string>()},
      {"seconds",seconds},{"root_version",gROOT->GetVersion()},{"producer_code_sha256",digest(p.at("producer_code"))},
      {"reduction_settings_sha256",digest(Json{{"selection",p.at("selection")},{"histograms",p.at("histograms")}})}}.dump(1)+"\n");
    completed=total; progress("complete");
  };
  if (fs::exists(output/"done.json")) {
    const Json saved=readJson((output/"done.json").string());
    if (saved.at("descriptor_sha256")!=descriptorSha || saved.at("plan_sha256")!=p.at("plan_sha256") ||
        saved.at("sha256").get<std::string>()!=sha256File(package) || saved.at("processed_entries")!=entries)
      throw std::runtime_error("completed group differs from its input, descriptor or package");
    validatePackage(package,p,sample); completed=total; progress("complete"); return;
  }
  if (fs::exists(package)) {
    const Json intent=readJson((output/"publication.json").string());
    if (intent.at("descriptor_sha256")!=descriptorSha || intent.at("plan_sha256")!=p.at("plan_sha256") ||
        intent.at("sha256").get<std::string>()!=sha256File(package)) throw std::runtime_error("interrupted chunk publication differs");
    validatePackage(package,p,sample);
    finish(intent.at("sha256"),nullptr,"RECOVERED_VALIDATED"); return; // No event loop or invented timing.
  }
  std::string scratchHome=scratch;
  if (scratchHome.empty()) {const char *s=std::getenv("_CONDOR_SCRATCH_DIR"); if (s) scratchHome=s;}
  if (scratchHome.empty() || !fs::is_directory(scratchHome)) throw std::runtime_error("fill-group requires worker scratch or --scratch DIR");
  const std::string nonce=digest(Json{{"time",std::chrono::steady_clock::now().time_since_epoch().count()},
    {"thread",std::hash<std::thread::id>{}(std::this_thread::get_id())}}).substr(0,24);
  const fs::path temporary=fs::canonical(scratchHome)/("treetohists-"+group+"-"+nonce);
  if (!fs::create_directory(temporary)) throw std::runtime_error("cannot create unique row scratch");
  live=temporary/"progress.json";
  progress("starting");
  try {
    savePlan(p,temporary.string());
    if(jes)for(size_t v=1;v<3;++v)savePlan(jes->plans[v],(temporary/(v==1?"energy_up":"energy_down")).string());
    photonjet::progressObserver=[&](const Json &value) {
      const std::string phase=value.at("phase");
      // Only the group's own finish publishes a shared terminal progress record.
      progress(phase=="complete" || phase=="reused" ? "source_complete" : phase,value);
    };
    for (size_t i=0;i<p.at("jobs").size();++i) {
      if(jes)fillDataJes(temporary,temporary/"energy_up",temporary/"energy_down",i,true);
      else fillRow(temporary.string(),i);
      completed+=photonjet::progressTotal;
      photonjet::progressPath.clear();
    }
    photonjet::progressObserver={}; progress("merging_chunk");
    if(jes){
      std::array<std::string,3> hashes;
      long long actualScienceReaders=0,actualScienceEntries=0;
      for(const auto& job:p.at("jobs")){
        const auto raw=fs::path(jobDirectory((temporary/"energy_up").string(),job))/"spectra.root";
        TFile source(raw.c_str(),"READ");const auto observed=Json::parse(text(source,"data_jes_fusion"));
        const auto& src=sample.at("sources").at(job.at("source").get<size_t>());
        const auto expectedEntries=job.value("range_end",src.at("metadata").at("entries").get<long long>())-job.value("range_begin",0LL);
        if(observed.at("science_readers")!=1 || observed.at("original_entries")!=expectedEntries)
          throw std::runtime_error("DATA JES original traversal counters differ from sealed row range");
        actualScienceReaders+=observed.at("science_readers").get<long long>();actualScienceEntries+=observed.at("original_entries").get<long long>();
      }
      if(actualScienceEntries!=entries)throw std::runtime_error("DATA JES original group entry counters differ");
      for(size_t v=0;v<3;++v){
        const auto& vp=jes->plans[v];const auto& vs=vp.at("samples").at(0);
        const fs::path local=v==0?temporary:temporary/(v==1?"energy_up":"energy_down");
        merge(vp,vs,local.string(),(local/"packages").string());
        const auto original=local/"packages"/(vs.at("key").get<std::string>()+".root");
        const auto target=jes->outputs[v]/(vs.at("key").get<std::string>()+".root");
        const auto partial=target.string()+".partial."+nonce;hashes[v]=sha256File(original.string());
        fs::create_directories(target.parent_path());fs::copy_file(original,partial);
        if(sha256File(partial)!=hashes[v])throw std::runtime_error("DATA JES group publication copy differs");
        if(::link(partial.c_str(),target.c_str())!=0)throw std::runtime_error("DATA JES group publication collided; exact reconciliation required");
        fs::remove(partial);
        writeText(jes->outputs[v]/"publication.json",Json{{"plan_sha256",vp.at("plan_sha256")},
          {"descriptor_sha256",jes->descriptorHashes[v]},{"sha256",hashes[v]}}.dump()+"\n");
        validatePackage(target.string(),vp,vs);
      }
      for(size_t v=1;v<3;++v){
        const auto& vp=jes->plans[v];const auto& vs=vp.at("samples").at(0);
        writeText(jes->outputs[v]/"done.json",Json{{"plan_sha256",vp.at("plan_sha256")},{"descriptor_sha256",jes->descriptorHashes[v]},
          {"sha256",hashes[v]},{"status","COMPLETE"},{"processed_entries",entries},{"source_count",vs.at("sources").size()},
          {"sample_class","PP/data"},{"seconds",nullptr},{"root_version",gROOT->GetVersion()},
          {"producer_code_sha256",digest(vp.at("producer_code"))},
          {"reduction_settings_sha256",digest(Json{{"selection",vp.at("selection")},{"histograms",vp.at("histograms")}})},
          {"data_jes_fusion",{{"physical_work_owner",p.at("plan_sha256")},{"physical_event_passes",0},{"shared_science_readers",actualScienceReaders},{"shared_science_entries",actualScienceEntries},{"qa_extensions","nominal_reused"}}}}.dump(1)+"\n");
        writeText(jes->outputs[v]/"progress.json",Json{{"campaign",vp.at("plan_sha256")},{"row",group},{"sample",vs.at("key")},
          {"phase","complete"},{"unit","event_passes"},{"processed",0},{"total",0},{"remaining",0},
          {"row_processed",0},{"row_total",0},{"row_remaining",0},
          {"deterministic_skip_processed",0},{"deterministic_skip_total",0},{"deterministic_skip_remaining",0},
          {"rate",0},{"eta_seconds",0},{"updated_at_unix",std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},
          {"latest_error",nullptr},{"input_entries",entries},{"source",group},
          {"data_jes_fusion",{{"physical_work_owner",p.at("plan_sha256")},{"physical_progress",(output/"progress.json").string()},{"physical_event_passes",0}}}}.dump()+"\n");
      }
      finish(hashes[0],std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count(),"COMPLETE");
      auto ownerDone=readJson((output/"done.json").string());
      ownerDone["data_jes_fusion"]={{"views",3},{"physical_event_passes",total},{"science_readers",actualScienceReaders},{"science_entries",actualScienceEntries},{"qa_passes",1},{"extensions_passes",1}};
      writeText(output/"done.json",ownerDone.dump(1)+"\n");
      fs::remove_all(temporary);return;
    }
    merge(p,sample,temporary.string(),(temporary/"packages").string());
    const std::string original=(temporary/"packages"/(sample.at("key").get<std::string>()+".root")).string();
    const std::string expected=sha256File(original), partial=package+".partial."+nonce;
    fs::create_directories(fs::path(package).parent_path());
    fs::copy_file(original,partial,fs::copy_options::overwrite_existing);
    if (sha256File(partial)!=expected) throw std::runtime_error("chunk publication copy differs");
    const bool first=::link(partial.c_str(),package.c_str())==0;  /* fails when another job has published the group */
    fs::remove(partial);
    if (!first) {
      for (int waited=0; waited<300 && !fs::exists(output/"done.json"); ++waited) std::this_thread::sleep_for(std::chrono::seconds(1));
      if (!fs::exists(output/"done.json")) throw std::runtime_error("another job published this group and left no receipt");
      validatePackage(package,p,sample); completed=total; fs::remove_all(temporary); return;
    }
    writeText(output/"publication.json",Json{{"plan_sha256",p.at("plan_sha256")},
      {"descriptor_sha256",descriptorSha},{"sha256",expected}}.dump()+"\n");
    validatePackage(package,p,sample);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    finish(expected,seconds,"COMPLETE");
    fs::remove_all(temporary);
  } catch (const std::exception &error) {
    photonjet::progressObserver={};
    const Json failedProgress=jes && !photonjet::progressPath.empty() ?
      Json{{"processed",photonjet::progressOffset+photonjet::progressRowProcessed}} : Json::object();
    progress("failed",failedProgress,error.what());
    throw; // Failed row scratch is retained for diagnosis; no other row is touched.
  }
}

/* Default group execution retains its original recipe and lifecycle. */
void fillGroup(const std::string& work,const std::string& group,const std::string& scratch){fillGroupImpl(work,group,scratch);}

void fillDataJesGroup(const std::string& nominal,const std::string& up,const std::string& down,
                     const std::string& group,const std::string& scratch){
  if(!std::regex_match(group,std::regex("G[0-9]+")))throw std::runtime_error("invalid fill group");
  const std::string shard=std::to_string(std::stoul(group.substr(1))/1000);
  DataJesGroup views;views.homes={fs::canonical(nominal),fs::canonical(up),fs::canonical(down)};
  if(views.homes[0]==views.homes[1] || views.homes[0]==views.homes[2] || views.homes[1]==views.homes[2])
    throw std::runtime_error("DATA JES group work identities must be distinct");
  const std::array<double,3> scales={1./.916,1./.907,1./.923};Json baseline;
  for(size_t v=0;v<3;++v){
    auto d=readJson((views.homes[v]/"fill_groups"/shard/(group+".json")).string());
    views.descriptorHashes[v]=d.at("descriptor_sha256");d.erase("descriptor_sha256");
    if(digest(d)!=views.descriptorHashes[v])throw std::runtime_error("DATA JES group descriptor digest differs");
    views.plans[v]=d.at("plan");const auto& p=views.plans[v];
    if(p.at("plan_sha256")!=readJson((views.homes[v]/"binding.json").string()).at("plan_sha256"))
      throw std::runtime_error("DATA JES group immutable plan binding differs");
    checkReductionCode(p);
    if(p.at("root_version").get<std::string>()!=gROOT->GetVersion())throw std::runtime_error("DATA JES group ROOT version differs");
    if(p.at("samples").size()!=1 || p.at("jobs").empty())throw std::runtime_error("DATA JES group must contain one bounded sample and jobs");
    const auto& sample=p.at("samples").at(0);const auto& h=p.at("histograms");
    if(sample.at("system")!="PP" || sample.at("family")!="data" || sample.contains("reused_package") ||
       sample.contains("_source_files") || h.at("response_variation")!="" || h.at("data_jet_scale").at("pp")!=scales[v] || h.at("data_jet_scale").at("auau")!=1.)
      throw std::runtime_error("DATA JES group requires exact pp DATA scales and original-only nominal response inputs");
    std::set<std::string> sources;
    for(const auto& source:sample.at("sources"))if(!sources.insert(source.at("identity").get<std::string>()).second)
      throw std::runtime_error("duplicate original source in DATA JES group");
    auto compare=p;compare.erase("plan_sha256");compare["histograms"]["data_jet_scale"]["pp"]=scales[0];
    if(v==0)baseline=compare;else if(compare!=baseline)throw std::runtime_error("DATA JES group source/range/selection/settings bindings differ");
    views.outputs[v]=fs::weakly_canonical(views.homes[v]/"chunks"/shard/group);
    if(fs::exists(views.outputs[v]) && !fs::is_empty(views.outputs[v]))
      throw std::runtime_error("DATA JES group refuses existing results; exact reconciliation required");
  }
  if(views.outputs[0]==views.outputs[1] || views.outputs[0]==views.outputs[2] || views.outputs[1]==views.outputs[2])
    throw std::runtime_error("DATA JES group output identities alias");
  fillGroupImpl(views.homes[0].string(),group,scratch,&views);
}
}  /* namespace treetohists */


/* === Centrality weights ==============================================================================================
   Centrality: the Au+Au simulation centrality weights, w(c_i) = p_DATA(c_i) / p_SIM(c_i) for the 5% slices c_i of 0-80%. */

/* tree_to_hists centrality --settings CFG --inputs INPUTS.yaml --data AuAu/data.root --out weights.json */

/* p_X(c_i) is the fraction of the centrality distribution of X inside 0-80% in slice i. DATA: the events that fired the
   selection's Au+Au trigger, minimum-bias classified and inside the vertex window (qa/events/centrality_scaled_bits of
   the DATA package), per source, found by name in the source's own trigger table: a source where the trigger is absent
   or disabled (initial prescale < 0) contributes nothing and is listed. Simulation: hist/SIM/h_centrality of every file
   of the catalogue's Au+Au photon samples (the response's), summed with each sample's cross section per generated event. Slices are
   [5i, 5i + 5) on the histograms' own edges (floor(c / 5), as the fills apply them). The weights are the catalogue's
   centrality_weights; weights.json also records every count and input. Run it after the DATA package, before the
   simulation is planned. */
namespace treetohists {
namespace fs = std::filesystem;
namespace {

/* Intermediate packages retain the full plan's coverage. Corrections report
   only the source records they actually used, including an unknown total. */
Json dataCorrectionCoverage(const Json &metadata, std::size_t sources) {
  Json coverage=metadata.value("coverage",Json::object());
  if(!coverage.is_object()) coverage=Json::object();
  coverage["selected_sources"]=sources;
  if(coverage.contains("population_sources") && !coverage.at("population_sources").is_null()) {
    const auto population=coverage.at("population_sources").get<long long>();
    if(population<static_cast<long long>(sources)) throw std::runtime_error("centrality DATA coverage exceeds its population");
    coverage["status"]=population==static_cast<long long>(sources) ? "COMPLETE" : "PARTIAL";
  } else coverage["status"]="UNKNOWN";
  return coverage;
}


/* Counts in the slices [5i, 5i + 5) of 0-80% of a histogram axis; every slice edge must be a bin edge. */
std::vector<double> sliceCounts(const std::vector<double> &values, const std::vector<double> &edges) {
  std::vector<double> out(16);
  std::vector<int> at;
  for (int i = 0; i <= 16; ++i) {
    const auto it = std::find_if(edges.begin(), edges.end(), [&](double e) { return std::fabs(e - 5.0 * i) < 1e-9; });
    if (it == edges.end()) throw std::runtime_error("the 5% slice edges of 0-80% must be bin edges of the centrality histogram");
    at.push_back(it - edges.begin());
  }
  for (int i = 0; i < 16; ++i) for (int b = at[i]; b < at[i + 1]; ++b) out[i] += values[b];
  return out;
}

std::vector<double> axis(const TAxis *a) {
  std::vector<double> e;
  for (int i = 1; i <= a->GetNbins() + 1; ++i) e.push_back(a->GetBinLowEdge(i));
  return e;
}

}  /* namespace */

bool finalDataAggregatesValidated(const Json &receipt, const Json &m, const std::string &fileHash) {
  if (!receipt.is_object() || !receipt.contains("validation") ||
      m.value("schema", "") != "PhotonJetSampleV1" || m.value("package_kind", "") != "FINAL" ||
      m.value("system", "") != "AuAu" || m.value("family", "") != "data") return false;
  const auto &v = receipt.at("validation");
  const Json required = {"identity_and_settings", "histograms_and_sumw2",
                         "source_ranges_and_contributions", "event_references"};
  return v.value("schema", "") == "PhotonJetPackageValidationV1" && v.value("status", "") == "PASS" &&
      v.value("checks", Json()) == required && receipt.value("sha256", "") == fileHash &&
      receipt.at("sources") == m.at("source_count") && receipt.at("plan_sha256") == m.at("plan_sha256") &&
      receipt.at("package_binding_sha256") == m.at("package_binding_sha256") &&
      receipt.at("selection_sha256") == m.at("selection").at("selection_sha256") &&
      v.at("producer_code_sha256") == digest(m.at("render_selection").at("producer_code"));
}

Json centralityWeights(const Settings &s, const std::string &catalogue, const std::string &dataPackage) {
  const std::string trigger = s.selection.at("systems").at("auau").at("event_selection").at("trigger_name");
  /* DATA: the named trigger's column of every source's centrality x scaled-bit cells */
  std::unique_ptr<TFile> d(TFile::Open(dataPackage.c_str(), "READ"));
  if (!d || d->IsZombie()) throw std::runtime_error("cannot open " + dataPackage);
  const Json meta = Json::parse(static_cast<TObjString *>(d->Get("meta/package"))->GetString().Data());
  if (meta.at("system") != "AuAu" || meta.at("family") != "data") throw std::runtime_error("the DATA target needs the Au+Au DATA package");
  const std::string key = "qa/events/centrality_scaled_bits";
  auto *bits = dynamic_cast<TH2 *>(d->Get(key.c_str()));
  if (!bits) throw std::runtime_error(dataPackage + ": no " + key);
  const std::vector<double> centralityEdges = axis(bits->GetXaxis());
  std::vector<double> data(centralityEdges.size() - 1);
  Json runs = Json::array();
  auto *sources=d->Get<TTree>("records/sources"),*cells=d->Get<TTree>("records/source_contributions");
  if(!sources || !cells) throw std::runtime_error("centrality DATA source records are missing");
  // The existing hash loop pulses this counter while verifying large inputs.
  ReductionProgress dataProgress("centrality-data",2*sources->GetEntries()+1);
  const std::string dataHash=sha256File(dataPackage),configurationHash=sha256File(catalogue);
  fs::path sidecar=dataPackage;sidecar.replace_extension(".json");
  const std::string receiptHash=fs::exists(sidecar)?sha256File(sidecar.string()):"";
  bool aggregate= !receiptHash.empty() && finalDataAggregatesValidated(readJson(sidecar.string()),meta,dataHash);
  int commonBit=-1;
  if(aggregate && sources->GetEntries()!=meta.at("source_count").get<long long>())
    throw std::runtime_error("centrality DATA source records are missing or incomplete");
  TTreeReader rs(sources), rc(cells);
  TTreeReaderValue<std::string> json(rs, "json"), histogram(rc, "histogram");
  TTreeReaderValue<Long64_t> begin(rs, "contribution_begin"), end(rs, "contribution_end"), bin(rc, "bin");
  TTreeReaderValue<double> value(rc, "value");
  while (rs.Next()) {
    const Json source = Json::parse(*json);
    const Json &m = source.at("metadata");
    const Json table = m.contains("trigger_selection") ? m.at("trigger_selection") : m.value("trigger_run_info", Json());
    if (table.is_null()) throw std::runtime_error("missing per-source trigger names in " + dataPackage);
    int bit = -1;
    for (size_t b = 0; b < table.at("name").size(); ++b)
      if (table.at("name")[b].get<std::string>() == trigger) {
        if (bit >= 0) throw std::runtime_error("ambiguous trigger name in a source trigger table");
        bit = b;
      }
    const bool enabled = bit >= 0 && table.at("initial_prescale")[bit].get<int>() >= 0;
    Json record = {{"source", runs.size()}, {"input_sha256", source.at("sha256")}, {"bit", bit >= 0 ? Json(bit) : Json()}, {"enabled", enabled}};
    if (enabled && commonBit<0) commonBit=bit;
    if (enabled && bit!=commonBit) aggregate=false;
    runs.push_back(record);
    dataProgress.advance(runs.size());
  }
  aggregate=aggregate && commonBit>=0;
  /* A shared enabled trigger column commutes with source summation. Subtract
     that column for disabled/absent sources, where the global bit may mean
     another trigger. Mixed enabled mappings keep per-source extraction. */
  if (aggregate) {
    if(bits->GetSumw2N()!=bits->GetNcells()) throw std::runtime_error("centrality DATA has no retained Sumw2");
    // DATA entries have unit weight. Sumw2 stores their exact double-precision
    // counts even when the TH2F contents rounded during a large final merge.
    for (size_t b=0;b<data.size();++b) {
      const int cell=bits->GetBin(b+1,commonBit+1);
      const double count=bits->GetSumw2()->At(cell),content=bits->GetBinContent(cell);
      if(!std::isfinite(count) || count<0 || count!=std::floor(count) ||
         std::abs(content-count)>1e-6*std::max(1.,count))
        throw std::runtime_error("validated DATA centrality needs integer unit counts and matching Sumw2");
      data[b]=count;
    }
  }
  rs.Restart();size_t ordinal=0;
  while (rs.Next()) {
      auto &record=runs.at(ordinal++);
      const bool enabled=record.at("enabled").get<bool>();
      dataProgress.advance(sources->GetEntries()+ordinal-1);
      if (aggregate?enabled:!enabled) continue;
      const int bit=aggregate?commonBit:record.at("bit").get<int>();
      std::vector<double> column(centralityEdges.size() - 1);
      for (Long64_t i = *begin; i < *end; ++i) {
        if (rc.SetEntry(i) != TTreeReader::kEntryValid) throw std::runtime_error("cannot read a source contribution");
        if (*histogram != key) continue;
        int bx, by, bz;
        bits->GetBinXYZ(static_cast<int>(*bin), bx, by, bz);
        if (by - 1 == bit && bx >= 1 && bx <= static_cast<int>(column.size())) column[bx - 1] += *value;
      }
      for (size_t i = 0; i < column.size(); ++i) data[i] += (aggregate?-1:1)*column[i];
      if(enabled) record["counts"] = sliceCounts(column, centralityEdges);
  }
  dataProgress.advance(2*sources->GetEntries());
  const std::vector<double> dataSlices = sliceCounts(data, centralityEdges);
  /* simulation: every Au+Au simulated sample of the catalogue, weighted by its cross section per generated event */
  const Json c = readInputCollection(catalogue);
  const fs::path base = fs::absolute(catalogue).parent_path();
  std::vector<double> simulation(16);
  Json samples = Json::object();
  Json snapshot=Json::array(),publications=Json::array();
  std::set<std::string> allInputPaths;
  for (const auto &sample : c.at("samples")) {
    if (sample.at("system") != "AuAu" || sample.at("family") != "photonjet") continue;  /* the response's samples */
    std::map<std::string,std::string> inputFiles;
    for(const auto &file:sample.value("files",Json::array())) {
      const auto path=fs::canonical(base/(file.is_string() ? file.get<std::string>() : file.at("path").get<std::string>())).string();
      const auto hash=file.is_object() ? file.at("sha256").get<std::string>() : sha256File(path);
      if(!inputFiles.emplace(path,hash).second)throw std::runtime_error("duplicate centrality input file: "+path);
    }
    for(const auto &name:sample.value("groups",Json::array())) {
      const auto group=fs::canonical(base/name.get<std::string>());
      const auto manifest=fs::canonical(group/"CURRENT")/"MANIFEST.json";
      publications.push_back({{"path",manifest.string()},{"sha256",sha256File(manifest.string())}});
      for(const auto &[path,hash]:groupFiles(group)) if(!inputFiles.emplace(path.string(),hash).second)
        throw std::runtime_error("duplicate centrality group/file input: "+path.string());
    }
    if(inputFiles.empty())throw std::runtime_error("empty centrality simulation sample");
    for(const auto &[path,hash]:inputFiles) {
      if(!allInputPaths.insert(path).second)throw std::runtime_error("centrality source is repeated across samples: "+path);
      if(sha256File(path)!=hash)throw std::runtime_error("centrality input changed: "+path);
    }
    const Json &n = sample.at("normalization");
    double weight = -1;
    if (n.value("source_factor", "") == "auau_embedded_stitching") {
      /* the canonical stitching: cross section over the sample's counted denominator (its files' events). The
         producer's h_centrality holds every event of a file; the embedded background, and so the centrality, does not
         depend on the generator-filter photon that decides the ownership. */
      const Json record = stitchingNormalization(s.home, "AuAu", sample.at("name"), "photonjet");
      long long events = 0;
      for (const auto &[path,hash] : inputFiles)
        events += stitchingDenominator(path, record).at("events").get<long long>();
      weight = events > 0 ? record.at("cross_section_pb").get<double>() / events : -1;
    } else {
      const double numerator = n.value("source_factor", "") == "family_reference"
                                   ? c.value("reference_cross_section_pb", Json::object()).value(n.value("reference", ""), -1.0)
                                   : n.value("cross_section_pb", -1.0);
      weight = numerator / n.value("generated_events", -1.0);
    }
    if (!(weight > 0)) throw std::runtime_error(sample.at("name").get<std::string>() + ": no positive cross section per generated event");
    std::vector<double> summed;
    std::vector<double> edges;
    Json files = Json::array();
    for (const auto &[path,hash] : inputFiles) {
      std::unique_ptr<TFile> in(TFile::Open(path.c_str(), "READ"));
      auto *h = in ? dynamic_cast<TH1 *>(in->Get("hist/SIM/h_centrality")) : nullptr;
      if (!h) throw std::runtime_error(path + " has no hist/SIM/h_centrality");
      const std::vector<double> these = axis(h->GetXaxis());
      if (edges.empty()) edges = these, summed.assign(these.size() - 1, 0.0);
      else if (these != edges) throw std::runtime_error(path + ": hist/SIM/h_centrality has other bin edges");
      for (int b = 1; b <= h->GetNbinsX(); ++b) summed[b - 1] += h->GetBinContent(b);
      if(sha256File(path)!=hash)throw std::runtime_error("centrality input changed during reading: "+path);
      files.push_back(Json{{"path", path}, {"sha256", hash}});
    }
    if (summed.empty()) continue;
    const std::vector<double> counts = sliceCounts(summed, edges);
    for (int i = 0; i < 16; ++i) simulation[i] += weight * counts[i];
    samples[sample.at("name").get<std::string>()] = Json{{"weight_pb_per_event", weight}, {"counts", counts}, {"files", files}};
    snapshot.push_back({{"sample",sample.at("name")},{"files",files},{"normalization",n},{"weight_pb_per_event",weight}});
  }
  double dataSum = 0, simulationSum = 0;
  for (int i = 0; i < 16; ++i) dataSum += dataSlices[i], simulationSum += simulation[i];
  const std::string counts = " (DATA per slice " + Json(dataSlices).dump() + ", simulation " + Json(simulation).dump() + ")";
  if (!(dataSum > 0) || !(simulationSum > 0)) throw std::runtime_error("an empty centrality distribution" + counts);
  std::vector<double> weights(16);
  for (int i = 0; i < 16; ++i) {
    if (!(simulation[i] > 0)) throw std::runtime_error("a simulated slice without events: " + std::to_string(i) + counts);
    if (!(dataSlices[i] > 0)) throw std::runtime_error("a data slice without events: " + std::to_string(i) + counts);
    weights[i] = (dataSlices[i] / dataSum) / (simulation[i] / simulationSum);
  }
  int with = 0;
  for (const auto &r : runs) with += r.at("enabled").get<bool>();
  for(const auto &publication:publications)if(sha256File(publication.at("path"))!=publication.at("sha256").get<std::string>())throw std::runtime_error("centrality publication manifest changed");
  if(sha256File(dataPackage)!=dataHash || sha256File(catalogue)!=configurationHash ||
     (!receiptHash.empty() && sha256File(sidecar.string())!=receiptHash))
    throw std::runtime_error("centrality DATA/configuration changed during derivation");
  dataProgress.finish();
  return Json{{"schema", "AuAuCentralityWeightsV1"}, {"status","DERIVED_UNQUALIFIED"},
              {"input_configuration_sha256",configurationHash},{"input_snapshot_sha256",digest(snapshot)},
              {"publication_bindings",publications},{"data_coverage",dataCorrectionCoverage(meta,runs.size())},
              {"formula", "w(c_i) = p_DATA(c_i) / p_SIM(c_i), 5% slices of 0-80%"},
              {"per_slice", weights}, {"centrality_weights_sha256", digest(weights)},
              {"data", {{"package", fs::absolute(dataPackage).string()}, {"package_sha256", dataHash}, {"histogram", key},
                        {"trigger", trigger}, {"counts", dataSlices}, {"sources_with_trigger", with},
                        {"aggregation",aggregate?"validated_final_common_trigger":"per_source_trigger_columns"},
                        {"validation_receipt_sha256",receiptHash},
                        {"sources_without_trigger", static_cast<int>(runs.size()) - with}, {"sources", runs}}},
              {"simulation", {{"histogram", "hist/SIM/h_centrality"}, {"counts", simulation}, {"samples", samples}}},
              {"selection_sha256", s.selection.at("selection_sha256")}};
}

}  /* namespace treetohists */

#include "SelectionFanout.cc"

/* Command dispatch. Settings/plan validation precedes filling; isolation's
   derivation commands live in Isolation.cc and use the shared definitions. */
int main(int argc, char **argv) {
  using namespace treetohists;
  gROOT->SetBatch(true);
  const std::string command = argc > 1 ? argv[1] : "";
  try {
    if(command=="selection-groups"){selection_fanout::localGroups(option(argc,argv,"--work"));return 0;}
    if(command=="selection-plan"){
      selection_fanout::planView(option(argc,argv,"--work"),option(argc,argv,"--out"),option(argc,argv,"--variant"),option(argc,argv,"--prescriptions"));return 0;
    }
    if(command=="fill-selection-group" || command=="validate-selection-group"){
      std::vector<fs::path> works={option(argc,argv,"--work")};
      for(int i=2;i+1<argc;++i)if(std::string(argv[i])=="--view-work")works.emplace_back(argv[i+1]);
      if(command=="validate-selection-group")std::cout<<selection_fanout::validateSaved(works,option(argc,argv,"--group")).dump(2)<<std::endl;
      else selection_fanout::group(works,option(argc,argv,"--group"),option(argc,argv,"--scratch"));return 0;
    }
    if (command == "plan") {
      std::vector<std::string> only;
      for (int i = 2; i + 1 < argc; ++i) if (std::string(argv[i]) == "--sample") only.push_back(argv[i + 1]);
      const fs::path work = option(argc, argv, "--work");
      const std::string max = option(argc, argv, "--max-events", false);
      const Json p = preparePlan(option(argc, argv, "--settings"), option(argc, argv, "--inputs"), option(argc, argv, "--home"), only,
        max.empty() ? 100000 : std::stoll(max), work.string());
      savePlan(p, work.string());
      std::cout << "[TreeToHists] " << p.at("samples").size() << " samples, " << p.at("jobs").size() << " fill jobs, selection "
                << p.at("selection").at("profile").get<std::string>() << ", plan " << (work / "plan.json").string() << std::endl;
      return 0;
    }
    if (command == "fill-data-jes") {
      fillDataJes(option(argc,argv,"--work"),option(argc,argv,"--up-work"),option(argc,argv,"--down-work"),std::stoi(option(argc,argv,"--job")));
      return 0;
    }
    if (command == "fill") {
      fill(option(argc, argv, "--work"), std::stoi(option(argc, argv, "--job")));
      return 0;
    }
    if (command == "fill-data-jes-group") {
      fillDataJesGroup(option(argc,argv,"--work"),option(argc,argv,"--up-work"),option(argc,argv,"--down-work"),option(argc,argv,"--group"),option(argc,argv,"--scratch",false));
      return 0;
    }
    if (command == "fill-group") {
      fillGroup(option(argc,argv,"--work"),option(argc,argv,"--group"),option(argc,argv,"--scratch",false));
      return 0;
    }
    if (command == "centrality") {
      const Settings s = Settings::read(option(argc, argv, "--settings"), option(argc, argv, "--home"));
      const fs::path out = option(argc, argv, "--out");
      if (fs::exists(out)) throw std::runtime_error(out.string() + " exists");
      const Json w = centralityWeights(s, option(argc, argv, "--inputs"), option(argc, argv, "--data"));
      writeText(out, w.dump(1) + "\n");
      std::cout << "[TreeToHists] centrality weights (" << w.at("data").at("sources_with_trigger") << " DATA sources with the trigger): ";
      for (const auto &x : w.at("per_slice")) std::cout << x.get<double>() << " ";
      std::cout << std::endl;
      return 0;
    }
    if (command == "merge") {
      const fs::path work = option(argc, argv, "--work");
      const Json p = readJson((work / "plan.json").string());
      merge(p, sampleOf(p, option(argc, argv, "--sample")), work.string(), option(argc, argv, "--out"));
      return 0;
    }
    if (command == "campaign" || command == "campaign-inspect" || command == "campaign-plan" || command == "campaign-submit-check" || command == "merge-group" || command == "collection")
      return campaign(command, std::vector<std::string>(argv + 2, argv + argc));
    if (command == "isolation") return isolation(std::vector<std::string>(argv + 2, argv + argc));
    if (command == "vertex") {
      if(reductionCode()!=Json::parse(TREETOHISTS_COMPILED_CODE)) throw std::runtime_error("vertex executable differs from current sources; rebuild first");
      return vertex(std::vector<std::string>(argv + 2, argv + argc));
    }
    std::cerr << "usage: tree_to_hists plan|fill|merge|centrality|isolation|vertex ... (see tree_to_hists.sh)" << std::endl;
    return 2;
  } catch (const std::exception &e) {
    if (!photonjet::progressPath.empty()) {
      try {
        Json progress = readJson(photonjet::progressPath); progress["latest_error"] = e.what(); progress["phase"] = "failed";
        writeText(photonjet::progressPath, progress.dump() + "\n");
      } catch (...) { /* retain the original failure */ }
    }
    std::cerr << "TreeToHists: " << e.what() << std::endl;
    return 1;
  }
}
