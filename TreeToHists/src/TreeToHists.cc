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
#include "Configuration.h"
#include "FillPhotonQA.C"
#include "FillPhotonJet.C"
#include "FillPhotonJetResponse.C"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

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
  Json axes;                        /* native_xj_edges, native_reco_response_edges_gev, native_truth_response_edges_gev */
  bool extensions = true, requireStage3Fractions = false;
  std::string responseVariation;    /* "" nominal, or the reco energy scale and resolution variation of a systematic run */
  static Settings read(const std::string &path, const std::string &home);
};
std::string sha256Text(const std::string &text);
std::string canonical(const Json &value);  /* compact, sorted keys, ASCII: the digest form */
std::string digest(const Json &value);     /* sha256 of canonical(value) */
Json readJson(const std::string &path);
/* A completed fill is reusable only with the code that produced its plan.
   Old work directories remain evidence; a changed recipe needs a new plan. */
Json reductionCode() {
  Json code = Json::object();
  const std::filesystem::path src = std::filesystem::path(__FILE__).parent_path();
  for (const char *name : {"TreeToHists.cc", "Isolation.cc", "Configuration.h", "PhotonJet.h", "FillPhotonJet.C", "FillPhotonQA.C", "FillPhotonJetResponse.C"})
    code[name] = sha256File((src / name).string());
  return code;
}
void checkReductionCode(const Json &plan) {
  if (plan.value("producer_code", Json()) != reductionCode())
    throw std::runtime_error("reduction code differs from the plan; prepare a new work directory");
}
/* The vertex histogram is external to the eventTrees. Recheck its planned
   hash at execution, on completed-job reuse and before package publication. */
void checkVertexFit(const Json &sample) {
  const Json fit = sample.value("vertex_reweighting", Json());
  if (fit.is_null()) return;
  const std::string path = fit.at("path");
  if (!std::filesystem::is_regular_file(path) || sha256File(path) != fit.at("sha256").get<std::string>())
    throw std::runtime_error(sample.at("key").get<std::string>() + ": vertex fit is missing or changed since planning");
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
Json optionalSelection(const std::filesystem::path &path, const std::string &system) {
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
    const Json tight = line(w.at("cuts").at("WP70"), "WP70"), high = line(w.at("cuts").at("WP80"), "WP80");
    for (const auto &x : system == "pp" ? et : cent) {
      const double t = tight["intercept"].get<double>() + tight["slope"].get<double>() * x.get<double>();
      const double hh = high["intercept"].get<double>() + high["slope"].get<double>() * x.get<double>();
      if (!(0.1 < hh && hh <= t && t <= 1.0)) throw std::runtime_error("optional working-point bands overlap or are invalid");
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
  return Json{{"score_branch", score}, {"state_branch", nullptr}, {"score_metadata", "meta/photon_models"},
              {"model_name", release.at("name")}, {"release_sha256", release.at("release_sha256")}, {"domain", d},
              {"cells", cells}, {"selection_sha256", sha256File(path.string())}, {"working_points_sha256", sha256File(wpPath.string())}};
}

}  /* namespace */

std::string sha256Text(const std::string &text) {
  Sha256 h;
  h.update(reinterpret_cast<const unsigned char *>(text.data()), text.size());
  return h.hex();
}

std::string sha256File(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  Sha256 h;
  std::vector<char> buffer(1 << 22);
  while (in) {
    in.read(buffer.data(), buffer.size());
    h.update(reinterpret_cast<const unsigned char *>(buffer.data()), in.gcount());
  }
  return h.hex();
}

std::string canonical(const Json &value) { return value.dump(-1, ' ', true); }
std::string digest(const Json &value) { return sha256Text(canonical(value)); }

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
  std::set<std::string> known = {"Selection.Name", "Selection.PhotonIDPackage", "Selection.AuAuIsolationPackage", "Axes.XjEdges",
                                 "Axes.RecoResponseEdges", "Axes.TruthResponseEdges", "Fill.RequireStage3Fractions", "Fill.ResponseVariation"};

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
                          "IsolationThreshold", "IsolationSidebandGap", "PhotonIDSelection", "DataTrigger"})
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
                   {"centrality_min", 0.0}, {"centrality_max", 80.0}, {"package_sha256", nullptr}};
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
                        {"trigger_name", system == "auau" ? "Photon 12 GeV + MBD NS >= 2, vtx < 150 cm" : "Photon 4 GeV + MBD NS >= 1"},
                        {"data_run_periods", system == "auau" ? Json::object() : Json{{"0mrad", {47289, 51274}}, {"1p5mrad", {51274, 54000}}}}}},
                      {"verify_nominal_id", nominalScore && cuts == originalCuts && pidSha == sha256File(canonicalId.string())},
                      {"verify_nominal_isolation", isolation == originalIsolation && (system == "pp" || isoSha == sha256File(canonicalIso.string()))}});
    if (has("PhotonIDSelection")) {
      for (const char *k : {"ScoreBranch", "StateBranch", "ScoreMetadata", "Tight", "NonTightLow", "NonTightHigh"})
        if (has(k)) throw std::runtime_error("an optional release owns its score and working points; do not override " + prefix + k);
      const Json optional = optionalSelection(relative(value(env, prefix + "PhotonIDSelection")), system);
      for (const char *k : {"score_branch", "state_branch", "score_metadata"}) entry[k] = optional.at(k);
      entry["optional_model"] = optional;
      entry["verify_nominal_id"] = false;
    }
    resolved["systems"][system] = entry;
  }
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
  const std::vector<std::string> shapes = {"e11_over_e33", "e32_over_e35", "bdt", "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx", "et1", "et2", "et3", "et4"};
  const std::vector<std::string> widths = {"weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx"};
  const std::vector<std::pair<std::string, std::string>> response = {
      {"photon_reco", "reco_et"}, {"photon_truth", "truth_et"}, {"photon_response_truth_x_reco", "truth_et,reco_et"},
      {"photon_fakes", "reco_et"}, {"photon_misses", "truth_et"}, {"xj_reco", "reco_et,xj"}, {"xj_truth", "truth_et,xj"},
      {"xj_response_truth_global_x_reco_global", "truth_cell,reco_cell"}, {"xj_fakes", "reco_et,xj"}, {"xj_fake_photon_reco", "reco_et,xj"},
      {"xj_detector_fakes_reco", "reco_et,xj"}, {"xj_misses", "truth_et,xj"}, {"xj_combinatoric_reco", "reco_et,xj"},
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
      add("measurement/" + c + "/" + jet + "/recoil", "region_ac,et,xj", "spectra", "recoil_" + c + "_" + jet);
      for (const char *region : {"A", "C"})
        for (const char *side : {"near", "far"})
          add("measurement/" + c + "/" + jet + "/dphi" + region + "_" + side, "et,xj,dphi", "spectra",
              std::string("dphi") + region + "_" + c + "_" + jet + "_" + side);
      if (sim) for (const auto &[n, axes] : response) add("response/" + c + "/" + jet + "/" + n, axes, "response", n + "_" + c + "_" + jet);
    }
    for (const auto &label : labels) {
      const std::string dest = "qa/photons/" + c + "/" + label, tag = label + "_" + c;
      for (const auto &v : shapes) add(dest + "/shape_" + v, v + ",et,stage", "qa", "ss_" + v + "_" + tag);
      for (const auto &[n, axes] : std::vector<std::pair<std::string, std::string>>{
               {"ncb", "npb,et"}, {"iso", "iso,et,id"}, {"abcd", "region_abcd,et"}, {"recoil", "region_ac,et,xj"},
               {"scoreiso", "score,iso,et"}, {"clu_et", "wide_et,stage"}, {"clu_etaphi", "eta,phi,stage"}})
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
    if (!sim) {
      add("qa/photons/" + c + "/data/e11e33_dt", "e11_over_e33,et,time_minus_mbd", "qa", "e11e33_dt_data_" + c);
      add("qa/photons/" + c + "/data/etatime", "eta,time_minus_mbd", "qa", "etatime_data_" + c);
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
    for (const auto &v : widths) add(dest + "/width_" + v, v + ",et,centrality_slice", "qa", "width_" + v + "_" + label);
  }
  add("qa/events/calo_total", "calo_total,centrality,mb_state", "qa", "evt_calo_" + F);
  if (auau) {
    add("qa/events/centrality", "centrality,weight_mode", "qa", "evt_cent_" + F);
    if (!sim) add("qa/events/centrality_scaled_bits", "centrality,scaled_bit", "qa", "cent_bits_data");
  }
  if (sim) {
    add("qa/generator/leading_pt", "generator_pt,generator_kind", "qa", "gen_" + F);
    add("qa/photon_id/score_roles", "score,et,centrality_slice,truth_role", "extensions");
  } else {
    add("qa/events/trigger_bits_raw", "max_cluster_energy,scaled_bit,live_bit", "qa", "trig_data");
    add("qa/events/trigger_max_cluster_energy", "max_cluster_energy,semantic_trigger,trigger_state", "extensions");
  }
  for (const auto &[n, axes] : std::vector<std::pair<std::string, std::string>>{
           {"cutflow", "gate,centrality,weight_mode"}, {"weights", "event_weight,weight_state"}, {"mbd_charge", "mbd_charge,centrality"}})
    add("qa/events/" + n, axes, "extensions");
  for (const std::string d : {"emcal", "ihcal", "ohcal", "total"}) {
    add("qa/events/" + d + "_vs_mbd_charge", "calo_" + d + ",mbd_charge,mb_state", "extensions");
    if (d != "total") add("qa/events/" + d + "_vs_centrality", "calo_" + d + ",centrality,mb_state", "extensions");
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

/* An event-record table appended source by source: integers become Long64_t, event references move by the events
   already written, "file" and "source" carry the source index. */
struct Table {
  std::unique_ptr<TTree> output;
  std::map<std::string, Long64_t> ints;
  std::map<std::string, Double_t> doubles;
  void append(TTree &input, const std::string &name, Long64_t source, Long64_t offset, TDirectory *dest) {
    if (!output) {
      dest->cd();
      output.reset(new TTree(name.c_str(), input.GetTitle()));
      output->SetAutoFlush(-8 * 1024 * 1024);
      for (auto *obj : *input.GetListOfLeaves()) {
        auto *leaf = static_cast<TLeaf *>(obj);
        const std::string n = leaf->GetName(), t = leaf->GetTypeName();
        if (t == "Double_t" || t == "Float_t") output->Branch(n.c_str(), &(doubles[n] = 0), (n + "/D").c_str());
        else output->Branch(n.c_str(), &(ints[n] = 0), (n + "/L").c_str());
      }
      output->Branch("source", &(ints["source"] = 0), "source/L");
    }
    for (Long64_t row = 0; row < input.GetEntries(); ++row) {
      if (input.GetEntry(row) <= 0) throw std::runtime_error("cannot read an event record");
      for (auto &[n, v] : doubles) {
        auto *l = input.GetLeaf(n.c_str());
        if (!l) throw std::runtime_error("event-record schema changed");
        v = l->GetValue();
      }
      for (auto &[n, v] : ints) {
        if (n == "source") { v = source; continue; }
        auto *l = input.GetLeaf(n.c_str());
        if (!l) throw std::runtime_error("event-record schema changed");
        v = l->GetValueLong64();
        if (n == "event") v += offset;
        if (n == "file") v = source;
      }
      output->Fill();
    }
  }
  void write() {
    if (!output) return;
    output->GetDirectory()->cd();
    output->Write("", TObject::kOverwrite);
    output->SetDirectory(nullptr);
  }
};

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
  const std::string key = sample.at("key"), system = sample.at("system"), family = sample.at("family");
  const fs::path target = fs::path(outputDir) / (key + ".root"), partial = target.string() + ".partial";
  if (fs::exists(target)) throw std::runtime_error(target.string() + " exists; a package is never replaced");
  fs::create_directories(target.parent_path());
  fs::remove(partial);
  const Json rows = inventory(system, family);
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
  Table events, photons, jets;
  Long64_t eventOffset = 0;
  Json available = Json::object();
  std::string grids;
  std::vector<std::string> kinds = {"spectra", "qa", "extensions"};
  if (family != "data") kinds.push_back("response");
  long long denominator = 0, missing = 0;  /* the canonical Au+Au embedded stitching: the events the fill jobs counted */
  for (const auto &src : sample.at("sources")) {
    const fs::path dir = productDirectory(work, key, source);
    if (!fs::exists(dir / "done.json")) throw std::runtime_error("fill job incomplete: " + dir.string());
    const Json done = readJson((dir / "done.json").string());
    if (done.at("plan_sha256") != plan.at("plan_sha256") || done.at("input") != src.at("path") ||
        done.value("producer_code_sha256", "") != digest(plan.at("producer_code")))
      throw std::runtime_error("fill job belongs to another plan or reduction code: " + dir.string());
    if (done.contains("stitching_denominator"))
      denominator += done.at("stitching_denominator").at("events").get<long long>(), missing += done.at("stitching_denominator").at("missing").get<long long>();
    sourceJson = Json{{"path", src.at("path")}, {"sha256", done.at("sha256")}, {"identity", src.at("identity")}, {"metadata", src.at("metadata")}}.dump();
    contributionBegin = contributions.GetEntries();
    for (const auto &kind : kinds) {
      std::unique_ptr<TFile> in(TFile::Open((dir / (kind + ".root")).c_str(), "READ"));
      if (!in || in->IsZombie()) throw std::runtime_error("cannot open product " + (dir / kind).string());
      if (text(*in, "profile") != profile || text(*in, "selection") != selectionText) throw std::runtime_error("product selection differs: " + dir.string());
      if (kind == "response") {
        const std::string g = text(*in, "grids");
        if (!grids.empty() && grids != g) throw std::runtime_error("response grids differ between sources");
        grids = g;
      }
      if (kind == "extensions") {
        const Json a = Json::parse(text(*in, "availability"));
        if (!available.empty() && available != a) throw std::runtime_error("mixed extension availability; use a uniform input release");
        available = a;
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
          const bool runQA = kind == "qa" && (path == "qa/events/centrality_scaled_bits" || path == "qa/events/trigger_bits_raw" || path == "qa/events/centrality");
          for (int cell = 0; cell < h->GetNcells(); ++cell) {
            const double v = h->GetBinContent(cell), e = h->GetBinError(cell);
            if (!std::isfinite(v) || !std::isfinite(e)) throw std::runtime_error("nonfinite histogram cell in " + path);
            if ((kind == "response" || kind == "spectra" || runQA) && (v != 0 || e != 0)) {
              bin = cell, histogram = path, value = v, variance = h->GetSumw2()->At(cell);
              contributions.Fill();
            }
          }
          if (!sums.count(path)) {
            sums[path].reset(static_cast<TH1 *>(h->Clone()));
            sums[path]->SetDirectory(nullptr);
          } else {
            sameAxes(*sums[path], *h);
            if (!sums[path]->Add(h)) throw std::runtime_error("histogram addition failed: " + path);
          }
        } else if (auto *h = dynamic_cast<THnSparse *>(obj)) {
          if (!sparse.count(path)) sparse[path].reset(static_cast<THnSparse *>(h->Clone()));
          else sparse[path]->Add(h);
        } else {
          throw std::runtime_error("unsupported object " + path);
        }
      }
      if (kind == "spectra") {
        auto *ev = dynamic_cast<TTree *>(in->Get("table_events"));
        auto *ph = dynamic_cast<TTree *>(in->Get("table_photons"));
        auto *jt = dynamic_cast<TTree *>(in->Get("table_jets"));
        if (!ev || !ph || !jt) throw std::runtime_error("event records missing from " + dir.string());
        events.append(*ev, "events", source, eventOffset, records);
        photons.append(*ph, "photons", source, eventOffset, records);
        jets.append(*jt, "jets", source, eventOffset, records);
        eventOffset += ev->GetEntries();
      }
    }
    contributionEnd = contributions.GetEntries();
    sources.Fill();
    ++source;
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
  out->cd();
  TObjString(profile.c_str()).Write("profile");
  TObjString(selectionText.c_str()).Write("selection");
  if (!grids.empty()) TObjString(grids.c_str()).Write("grids");
  const std::string sys = system == "PP" ? "pp" : "auau";
  Json normalization = sample.at("normalization");
  if (normalization.contains("ownership")) {  /* cross section over the sample's events in its denominator */
    if (missing > 0) throw std::runtime_error(key + ": " + std::to_string(missing) + " events without the ownership variable; the stitching denominator is not exact");
    if (denominator <= 0) throw std::runtime_error(key + ": no events in the stitching denominator");
    normalization["denominator_events"] = denominator;
    normalization["pb_per_generated_event"] = normalization.at("cross_section_pb").get<double>() / denominator;
  }
  Json meta = {{"schema", "PhotonJetSampleV1"}, {"system", system}, {"name", sample.at("name")}, {"family", family},
               {"profile", profile}, {"selection", plan.at("selection")}, {"histograms", plan.at("histograms")},
               {"normalization", normalization}, {"centrality_weights", sample.at("weights")},
               {"centrality_evidence", sample.at("centrality_evidence")}, {"centrality_data_trigger", sample.at("centrality_data_trigger")},
               {"vertex_reweighting", sample.at("vertex_reweighting")}, {"plan_sha256", plan.at("plan_sha256")},
               {"runtime", {{"root_version", gROOT->GetVersion()}}}, {"histograms_written", sums.size() + sparse.size()},
               {"availability", available}, {"coverage", sample.value("coverage", Json::object())}, {"source_count", source}, {"event_records", eventOffset},
               {"acceptance", "LOCAL_STRUCTURE_ONLY_REQUIRES_PHYSICS_QUALIFICATION"}};
  /* These values come from the same resolved Profile used by every fill. They are descriptive metadata,
     not a second selection: changing a cut changes the fill and this record together. */
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
      {"jet", {{"collection", "jet04"}, {"radius", 0.4}, {"abs_eta_max", cuts.jetAbsEtaMax},
               {"pt_min_gev", cuts.jetPtMin}, {"overlap_delta_r", cuts.overlapDeltaR}, {"delta_phi_min_rad", cuts.deltaPhiMin},
               {"delta_phi_pi_fraction",angle}, {"pt_operator",">="}, {"delta_phi_operator",">="}, {"overlap_operator",">="}}},
      {"truth_photon", {{"isolation_radius", 0.3}, {"abs_eta_max", photonjet::kTruthPhotonAbsEtaMax},
                        {"isolation_max_gev", photonjet::kTruthIsolationMaxGeV}, {"signal_source_role", photonjet::kSignalSourceRole}}},
      {"qa", {{"et_min_gev", photonjet::qa::kQaEtMin}, {"et_max_gev", photonjet::qa::kQaEtMax},
              {"classes",qaClasses}, {"stored_abs_eta_max", photonjet::qa::kStoredAbsEtaMax}, {"jet_map_pt_min_gev", photonjet::qa::kQaJetMapPtMin}}},
      {"photon_id", plan.at("selection").at("systems").at(sys)},
      {"histogram_settings", plan.at("histograms")}};
  if (family == "data") {
    meta["data_trigger"] = plan.at("selection").at("systems").at(sys).at("data_trigger");
    meta["trigger_selection"] = "applied_in_package";  /* the measurement and records hold only this trigger's population */
  }
  directory(*out, "meta")->cd();
  TObjString(meta.dump().c_str()).Write("package");
  if (family != "data")  /* the normalisation FinalAnalysis applies once when it sums the samples of a result */
    TParameter<double>("pb_per_generated_event", normalization.at("pb_per_generated_event").get<double>()).Write();
  TObjString(rows.dump().c_str()).Write("inventory");
  contributions.SetDirectory(nullptr);
  sources.SetDirectory(nullptr);
  out->Close();
  check(partial.string(), rows, source);
  checkVertexFit(sample);
  fs::rename(partial, target);
  const Json sidecar = {{"package", target.filename().string()}, {"sha256", sha256File(target.string())}, {"sources", source},
                        {"histograms", meta.at("histograms_written")}, {"event_records", eventOffset}, {"availability", available},
                        {"plan_sha256", plan.at("plan_sha256")}};
  std::ofstream((target.parent_path() / (sample.at("name").get<std::string>() + ".json")).string()) << sidecar.dump(1) << "\n";
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
namespace {

void writeText(const fs::path &path, const std::string &text) {  /* atomic: a reader never sees a partial file */
  const fs::path partial = path.string() + ".partial";
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
  std::vector<std::string> needed = {spec.at("score_branch"), "pho_et", "isow_cone_sum", "isow_threshold", "isow_sideband_threshold"};
  if (system == "pp" && family != "data") {  /* the canonical weight state the fills check (photonjet::PpSimulationGate) */
    photonjet::ppWeightComponent(source);
    for (const char *name : {"ppg12_weight_active_state", "ppg12_weight_lane_code", "w_application_count", "w_event", "w_slice", "w_vertex", "w_sidi", "w_period"})
      needed.push_back(name);
  }
  if (!spec.at("state_branch").is_null()) needed.push_back(spec.at("state_branch"));
  for (const auto &name : needed)
    if (!tree->GetBranch(name.c_str())) throw std::runtime_error(path.string() + ": missing selection input " + name);
  const Json receipt = Json::parse(photonjet::metaText(*f, spec.at("score_metadata")));
  if (spec.contains("optional_model")) {
    const Json &o = spec.at("optional_model");
    const Json score = receipt.value(o.at("model_name").get<std::string>(), Json::object());
    for (const char *k : {"release_sha256", "score_branch", "domain"})
      if (score.value(k, Json()) != o.at(k)) throw std::runtime_error(path.string() + ": optional score is not bound to the selected release");
    if (score.value("status", "") != "PASS") throw std::runtime_error(path.string() + ": optional score is not validated");
  } else {
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

/* Expand and check every catalogued sample. */
Json plan(const Settings &s, const std::string &catalogue, const std::vector<std::string> &only) {
  const Json c = readInputCollection(catalogue);
  if (c.value("schema", "") != "PhotonJetInputCollectionV1" || c.value("samples", Json::array()).empty())
    throw std::runtime_error("expected PhotonJetInputCollectionV1 with samples");
  const fs::path base = fs::absolute(catalogue).parent_path();
  Json samples = Json::array(), jobs = Json::array();
  std::set<std::string> keys, identities;
  for (const auto &sample : c.at("samples")) {
    const std::string system = sample.at("system"), name = sample.at("name"), family = sample.at("family");
    if ((system != "PP" && system != "AuAu") || (family != "data" && family != "photonjet" && family != "inclusivejet") ||
        !std::regex_match(name, std::regex("[A-Za-z0-9_-]+")))
      throw std::runtime_error("invalid system, sample name or family");
    if ((name == "data") != (family == "data")) throw std::runtime_error("each system has one data sample named data");
    const std::string key = system + "/" + name;
    if (!only.empty() && std::find(only.begin(), only.end(), key) == only.end()) continue;
    if (!keys.insert(key).second) throw std::runtime_error("duplicate sample " + key);
    std::vector<std::pair<fs::path, std::string>> files;
    for (const auto &f : sample.value("files", Json::array())) files.emplace_back(fs::weakly_canonical(base / f.get<std::string>()), "");
    for (const auto &g : sample.value("groups", Json::array())) for (auto &f : groupFiles(fs::weakly_canonical(base / g.get<std::string>()))) files.push_back(f);
    if (files.empty()) throw std::runtime_error(key + " has no input files");
    const std::string expected = family == "data" ? "data" : family == "photonjet" ? "phosim" : "jetsim";
    Json sources = Json::array();
    std::set<std::string> populations;
    for (const auto &[path, sha] : files) {
      std::string fileFamily, fileSystem;
      Json info = sourceInfo(path, s.selection, fileFamily, fileSystem);
      if (fileSystem != (system == "PP" ? "pp" : "auau") || fileFamily != expected)
        throw std::runtime_error("sample identity disagrees with eventTree metadata: " + key + " " + path.string());
      if (!identities.insert(info.at("identity")).second) throw std::runtime_error("duplicate original source identity: " + info.at("identity").get<std::string>());
      populations.insert(info.at("metadata").at("source").value("sample", ""));
      const std::string observedSha = sha256File(path.string());
      if (!sha.empty() && observedSha != sha) throw std::runtime_error("published input differs from its group manifest: " + path.string());
      info["expected_sha256"] = observedSha;
      jobs.push_back(Json{{"sample", key}, {"source", sources.size()}, {"input", path.string()}});
      sources.push_back(info);
    }
    if (family != "data" && (populations.size() != 1 || populations.count("")))
      throw std::runtime_error(key + ": a simulation package must hold exactly one named source population");
    if (sample.contains("source_sample") && populations != std::set<std::string>{sample.at("source_sample").get<std::string>()})
      throw std::runtime_error(key + ": source_sample binding differs");
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
      for (const auto &src : sources)
        if (src.at("metadata").value("weighting", Json::object()).value("vertex_enabled", false)) throw std::runtime_error(key + ": input already includes vertex weights");
    } else if (!weights.empty() || sample.contains("vertex_reweighting")) {
      throw std::runtime_error(key + ": centrality weights and vertex fits apply only to Au+Au simulation");
    }
    Json normalization = sample.value("normalization", Json::object());
    if (family != "data" && normalization.value("source_factor", "") == "auau_embedded_stitching") {
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
    if (population != -1 && population < static_cast<long>(sources.size())) throw std::runtime_error(key + ": population_sources is smaller than the selected source count");
    const Json coverage = {{"selected_sources", sources.size()}, {"population_sources", population == -1 ? Json(nullptr) : Json(population)},
                           {"status", population == -1 ? "UNKNOWN" : population == static_cast<long>(sources.size()) ? "COMPLETE" : "PARTIAL"}};
    samples.push_back(Json{{"key", key}, {"system", system}, {"name", name}, {"family", family}, {"sources", sources}, {"coverage", coverage},
                           {"weights", weights}, {"vertex_reweighting", vertex}, {"normalization", normalization},
                           {"centrality_evidence", sample.value("centrality_evidence", Json())},
                           {"centrality_data_trigger", sample.value("centrality_data_trigger", Json())}});
  }
  if (samples.empty()) throw std::runtime_error("no catalogued sample matches");
  Json histograms = s.axes;
  histograms.update(Json{{"native_measurement_et_edges_gev", photonjet::photonEtEdges()}, {"jet_thresholds_gev", photonjet::kJetThresholds},
                         {"reported_photon_et_edges_gev", {17, 23, 33}}, {"require_stage3_fractions", s.requireStage3Fractions},
                         {"response_variation", s.responseVariation}});
  photonjet::response::parseVariation(s.responseVariation.c_str());  /* a malformed variation fails here, before any job */
  Json out = {{"schema", "TreeToHistsPlanV1"}, {"settings", s.path}, {"catalogue", fs::absolute(catalogue).string()},
              {"catalogue_sha256", sha256File(catalogue)}, {"selection", s.selection}, {"histograms", histograms},
              {"root_version", gROOT->GetVersion()}, {"producer_code", reductionCode()}, {"samples", samples}, {"jobs", jobs}};
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
  photonjet::OwnershipGate ownership(reader, o, true);
  long long events = 0, missing = 0;
  while (reader.Next()) {
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
void fill(const fs::path &work, int index) {
  const Json p = readJson((work / "plan.json").string());
  checkReductionCode(p);
  const Json &job = p.at("jobs").at(index);
  const Json &sample = sampleOf(p, job.at("sample"));
  checkVertexFit(sample);
  const Json &source = sample.at("sources").at(job.at("source").get<size_t>());
  const fs::path dir = productDirectory(work.string(), job.at("sample"), job.at("source").get<long>());
  if (fs::exists(dir / "done.json")) {
    const Json saved = readJson((dir / "done.json").string());
    if (saved.at("plan_sha256") != p.at("plan_sha256") || saved.at("sha256").get<std::string>() != sha256File(job.at("input")) ||
        saved.value("producer_code_sha256", "") != digest(p.at("producer_code")))
      throw std::runtime_error("completed fill differs from its input or plan");
    return;
  }
  fs::create_directories(dir);
  const auto start = std::chrono::steady_clock::now();
  const std::string input = job.at("input"), sha = sha256File(input);
  if (!source.at("expected_sha256").get<std::string>().empty() && sha != source.at("expected_sha256").get<std::string>())
    throw std::runtime_error("published input differs from its group manifest: " + input);
  const std::string selection = (work / "selection.json").string(), weights = listed(sample.at("weights"));
  photonjet::configuredXjEdges = p.at("histograms").at("native_xj_edges").get<std::vector<double>>();
  photonjet::strictReduction = true;
  photonjet::configureVertexFit(sample.at("vertex_reweighting").is_null() ? Json::object() : sample.at("vertex_reweighting"));
  photonjet::configuredOwnership = ownershipOf(sample.at("normalization"));
  const Json &h = p.at("histograms");
  std::vector<std::string> kinds = {"spectra", "qa", "extensions"};
  if (sample.at("family") != "data") kinds.push_back("response");
  for (const auto &kind : kinds) {
    const fs::path out = dir / (kind + ".root"), partial = dir / (kind + ".root.partial");
    fs::remove(partial);
    photonjet::progressPath = (dir / ("progress_" + kind + ".json")).string();
    photonjet::progressTotal = source.at("metadata").at("entries").get<long long>();
    int code = 0;
    if (kind == "spectra") code = FillPhotonJet(input.c_str(), partial.c_str(), "nominal", selection.c_str(), weights.c_str());
    else if (kind == "qa") code = FillPhotonQA(input.c_str(), partial.c_str(), "nominal", weights.c_str(), selection.c_str());
    else if (kind == "extensions") code = FillExtensions(input.c_str(), partial.c_str(), selection.c_str(), weights.c_str(), h.at("require_stage3_fractions").get<bool>());
    else code = FillPhotonJetResponse(input.c_str(), partial.c_str(), "nominal", listed(h.at("native_reco_response_edges_gev")).c_str(),
                                      listed(h.at("native_truth_response_edges_gev")).c_str(), listed(h.at("native_xj_edges")).c_str(),
                                      weights.c_str(), h.at("response_variation").get<std::string>().c_str(), selection.c_str());
    if (code) throw std::runtime_error("the " + kind + " filler failed on " + input);
    fs::rename(partial, out);
  }
  Json done = {{"input", input}, {"sha256", sha}, {"root_version", gROOT->GetVersion()}, {"plan_sha256", p.at("plan_sha256")},
               {"producer_code_sha256", digest(p.at("producer_code"))}};
  if (sample.at("normalization").contains("ownership")) done["stitching_denominator"] = stitchingDenominator(input, sample.at("normalization"));
  if (sha256File(input) != sha) throw std::runtime_error("input changed while it was being reduced: " + input);
  checkVertexFit(sample);
  done["seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  writeText(dir / "done.json", done.dump(1) + "\n");
}

std::string option(int argc, char **argv, const std::string &name, bool requiredOption = true) {
  for (int i = 2; i + 1 < argc; ++i) if (name == argv[i]) return argv[i + 1];
  if (requiredOption) throw std::runtime_error("missing " + name);
  return "";
}

}  /* namespace */
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
  TTreeReader rs(dynamic_cast<TTree *>(d->Get("records/sources"))), rc(dynamic_cast<TTree *>(d->Get("records/source_contributions")));
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
    if (enabled) {
      std::vector<double> column(centralityEdges.size() - 1);
      for (Long64_t i = *begin; i < *end; ++i) {
        if (rc.SetEntry(i) != TTreeReader::kEntryValid) throw std::runtime_error("cannot read a source contribution");
        if (*histogram != key) continue;
        int bx, by, bz;
        bits->GetBinXYZ(static_cast<int>(*bin), bx, by, bz);
        if (by - 1 == bit && bx >= 1 && bx <= static_cast<int>(column.size())) column[bx - 1] += *value;
      }
      for (size_t i = 0; i < column.size(); ++i) data[i] += column[i];
      record["counts"] = sliceCounts(column, centralityEdges);
    }
    runs.push_back(record);
  }
  const std::vector<double> dataSlices = sliceCounts(data, centralityEdges);
  /* simulation: every Au+Au simulated sample of the catalogue, weighted by its cross section per generated event */
  const Json c = readInputCollection(catalogue);
  const fs::path base = fs::absolute(catalogue).parent_path();
  std::vector<double> simulation(16);
  Json samples = Json::object();
  for (const auto &sample : c.at("samples")) {
    if (sample.at("system") != "AuAu" || sample.at("family") != "photonjet") continue;  /* the response's samples */
    const Json &n = sample.at("normalization");
    double weight = -1;
    if (n.value("source_factor", "") == "auau_embedded_stitching") {
      /* the canonical stitching: cross section over the sample's counted denominator (its files' events). The
         producer's h_centrality holds every event of a file; the embedded background, and so the centrality, does not
         depend on the generator-filter photon that decides the ownership. */
      const Json record = stitchingNormalization(s.home, "AuAu", sample.at("name"), "photonjet");
      long long events = 0;
      for (const auto &f : sample.value("files", Json::array()))
        events += stitchingDenominator(fs::weakly_canonical(base / f.get<std::string>()).string(), record).at("events").get<long long>();
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
    for (const auto &f : sample.value("files", Json::array())) {
      const std::string path = fs::weakly_canonical(base / f.get<std::string>()).string();
      std::unique_ptr<TFile> in(TFile::Open(path.c_str(), "READ"));
      auto *h = in ? dynamic_cast<TH1 *>(in->Get("hist/SIM/h_centrality")) : nullptr;
      if (!h) throw std::runtime_error(path + " has no hist/SIM/h_centrality");
      const std::vector<double> these = axis(h->GetXaxis());
      if (edges.empty()) edges = these, summed.assign(these.size() - 1, 0.0);
      else if (these != edges) throw std::runtime_error(path + ": hist/SIM/h_centrality has other bin edges");
      for (int b = 1; b <= h->GetNbinsX(); ++b) summed[b - 1] += h->GetBinContent(b);
      files.push_back(Json{{"path", path}, {"sha256", sha256File(path)}});
    }
    if (!sample.value("groups", Json::array()).empty()) throw std::runtime_error("centrality: list the simulated files (groups are not read here)");
    if (summed.empty()) continue;
    const std::vector<double> counts = sliceCounts(summed, edges);
    for (int i = 0; i < 16; ++i) simulation[i] += weight * counts[i];
    samples[sample.at("name").get<std::string>()] = Json{{"weight_pb_per_event", weight}, {"counts", counts}, {"files", files}};
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
  return Json{{"schema", "AuAuCentralityWeightsV1"}, {"formula", "w(c_i) = p_DATA(c_i) / p_SIM(c_i), 5% slices of 0-80%"},
              {"per_slice", weights}, {"centrality_weights_sha256", digest(weights)},
              {"data", {{"package", fs::absolute(dataPackage).string()}, {"package_sha256", sha256File(dataPackage)}, {"histogram", key},
                        {"trigger", trigger}, {"counts", dataSlices}, {"sources_with_trigger", with},
                        {"sources_without_trigger", static_cast<int>(runs.size()) - with}, {"sources", runs}}},
              {"simulation", {{"histogram", "hist/SIM/h_centrality"}, {"counts", simulation}, {"samples", samples}}},
              {"selection_sha256", s.selection.at("selection_sha256")}};
}

}  /* namespace treetohists */

/* Command dispatch. Settings/plan validation precedes filling; isolation's
   derivation commands live in Isolation.cc and use the shared definitions. */
int main(int argc, char **argv) {
  using namespace treetohists;
  gROOT->SetBatch(true);
  const std::string command = argc > 1 ? argv[1] : "";
  try {
    if (command == "plan") {
      std::vector<std::string> only;
      for (int i = 2; i + 1 < argc; ++i) if (std::string(argv[i]) == "--sample") only.push_back(argv[i + 1]);
      const fs::path work = option(argc, argv, "--work");
      const Settings s = Settings::read(option(argc, argv, "--settings"), option(argc, argv, "--home"));
      const Json p = plan(s, option(argc, argv, "--inputs"), only);
      fs::create_directories(work);
      if (fs::exists(work / "plan.json") && readJson((work / "plan.json").string()).at("plan_sha256") != p.at("plan_sha256"))
        throw std::runtime_error(work.string() + " holds another plan; choose a new work directory");
      writeText(work / "selection.json", s.selection.dump() + "\n");
      writeText(work / "plan.json", p.dump(1) + "\n");
      std::ostringstream list;
      for (size_t i = 0; i < p.at("jobs").size(); ++i) list << i << " " << p.at("jobs")[i].at("sample").get<std::string>() << " " << p.at("jobs")[i].at("input").get<std::string>() << "\n";
      writeText(work / "jobs.txt", list.str());
      std::cout << "[TreeToHists] " << p.at("samples").size() << " samples, " << p.at("jobs").size() << " fill jobs, selection "
                << s.selection.at("profile").get<std::string>() << ", plan " << (work / "plan.json").string() << std::endl;
      return 0;
    }
    if (command == "fill") {
      fill(option(argc, argv, "--work"), std::stoi(option(argc, argv, "--job")));
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
    if (command == "isolation") return isolation(std::vector<std::string>(argv + 2, argv + argc));
    std::cerr << "usage: tree_to_hists plan|fill|merge|centrality|isolation ... (see tree_to_hists.sh)" << std::endl;
    return 2;
  } catch (const std::exception &e) {
    std::cerr << "TreeToHists: " << e.what() << std::endl;
    return 1;
  }
}
