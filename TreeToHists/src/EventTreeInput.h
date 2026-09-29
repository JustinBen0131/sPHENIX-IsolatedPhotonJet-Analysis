// What every TreeToHists reader shares: the input file list, the system and metadata of an eventTree file, the
// centrality classes, slices and simulation weights, the recoil-jet thresholds and the native photon E_T and xJ grids.
#pragma once

#include "PhotonJetSelection.h"

#include <TFile.h>
#include <TObjString.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>
#include <map>
#include <fstream>
#include <chrono>

namespace photonjet {

// Set only by the versioned public lifecycle, before booking any histograms.
// Empty keeps the historical macro interface and its reference axes unchanged.
inline std::vector<double> configuredXjEdges;
inline bool strictReduction = false;
inline std::string progressPath;
inline long long progressTotal = 0;
inline void reportProgress(long long processed) {
  if (progressPath.empty()) return;
  static auto start=std::chrono::steady_clock::now(), last=start-std::chrono::seconds(2);
  const auto now=std::chrono::steady_clock::now();
  if (processed!=progressTotal && now-last<std::chrono::seconds(1)) return;
  last=now;
  const double seconds=std::chrono::duration<double>(now-start).count();
  const double rate=seconds>0 ? processed/seconds : 0;
  std::ofstream out(progressPath+".next");
  out.precision(17);
  out << "{\"processed\":" << processed << ",\"total\":" << progressTotal
      << ",\"remaining\":" << progressTotal-processed << ",\"elapsed_seconds\":" << seconds
      << ",\"rate\":" << rate << ",\"eta_seconds\":" << (rate>0 ? (progressTotal-processed)/rate : -1)
      << ",\"updated_at_unix\":" << std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()
      << ",\"first_error\":null}";
  out.close(); std::filesystem::rename(progressPath+".next",progressPath);
}

struct CentralityClass {
  std::string name;
  double low, high;  // low <= centrality < high; ignored in pp
};

inline std::vector<CentralityClass> centralityClasses(System system) {
  if (system == System::auau) return {{"cent0_20", 0, 20}, {"cent20_50", 20, 50}, {"cent50_80", 50, 80}};
  return {{"all", 0, 0}};
}

// The 5% centrality slice of an Au+Au event, floor(c / 5): 0..15 inside the classes (0-80%).
inline int centralitySlice(double centrality) { return static_cast<int>(std::floor(centrality / 5.0)); }

// Au+Au simulation centrality weights, one per 5% slice of 0-80% and 1 outside: w(c) = p_DATA(c) / p_SIM(c), which
// weights the simulated centrality distribution to the data's. They are derived for each production from its data
// and simulation centrality distributions, so they are an input. Every simulation fill is weighted by w_event x w(c).
struct CentralityWeights {
  std::vector<double> perSlice;  // empty: no centrality weight (pp)
  double operator()(double centrality) const {
    if (perSlice.empty() || !(centrality >= 0.0 && centrality < 80.0)) return 1.0;
    return perSlice[centralitySlice(centrality)];
  }
};

// Parse the weights of a system: 16 finite positive numbers for Au+Au, none for pp.
inline CentralityWeights centralityWeights(const std::string& text, System system) {
  CentralityWeights w;
  std::istringstream stream(text);
  for (double v; stream >> v;) {
    if (!(std::isfinite(v) && v > 0.0)) throw std::runtime_error("centrality weights must be finite and positive");
    w.perSlice.push_back(v);
  }
  if (!stream.eof()) throw std::runtime_error("centrality weights are not numbers: " + text);
  if (system == System::pp && !w.perSlice.empty()) throw std::runtime_error("pp takes no centrality weights");
  if (system == System::auau && w.perSlice.size() != 16)
    throw std::runtime_error("Au+Au simulation needs 16 centrality weights, one per 5% slice of 0-80%");
  return w;
}

const std::vector<double> kJetThresholds = {5.0, 7.0, 10.0, 12.0};

inline std::string thresholdName(double t) { return "jet" + std::to_string(static_cast<int>(t)); }

inline std::vector<double> photonEtEdges() {
  std::vector<double> e;
  for (int v = 15; v <= 35; ++v) e.push_back(v);  // native 1 GeV grid
  return e;
}

inline std::vector<double> xjEdges() {
  if (!configuredXjEdges.empty()) return configuredXjEdges;
  std::vector<double> e;
  for (int i = 0; i <= 40; ++i) e.push_back(std::round(i * 0.05 * 1e6) / 1e6);  // 0.05 steps to 2.0
  e.push_back(3.0);
  return e;
}

// |Delta eta| above which a jet is on the far side of the photon: the shape of the data-driven combinatoric estimate.
const double kFarDeltaEta = 0.8;

// pi/32 steps with 0.8 and 2.5 rad added, so the recoil window (dphi >= 7pi/8) and the sideband and fit
// ranges of the combinatorial subtraction are bin edges.
inline std::vector<double> deltaPhiEdges() {
  std::vector<double> e;
  for (int k = 0; k <= 32; ++k) e.push_back(k * M_PI / 32.0);
  e.push_back(0.8);
  e.push_back(2.5);
  std::sort(e.begin(), e.end());
  return e;
}

inline std::vector<std::string> inputPaths(const char* inputs) {
  std::vector<std::string> paths;
  std::map<std::string, std::string> pinnedGroups;
  std::istringstream stream(inputs);
  for (std::string path; stream >> path;) {
    const auto at = path.find("/CURRENT/");
    if (at != std::string::npos) {
      const std::string group = path.substr(0, at);
      auto found = pinnedGroups.find(group);
      if (found == pinnedGroups.end())
        found = pinnedGroups.emplace(group, std::filesystem::canonical(group + "/CURRENT").string()).first;
      path = found->second + path.substr(at + 8);
    }
    paths.push_back(path);
  }
  if (paths.empty()) throw std::runtime_error("no input files");
  return paths;
}

inline std::string metaText(TFile& f, const std::string& key) {
  auto* text = dynamic_cast<TObjString*>(f.Get(key.c_str()));
  if (!text) throw std::runtime_error("no " + key + " in " + f.GetName());
  return text->GetString().Data();
}

// The system of an eventTree file (meta/schema is compact JSON, Format tables.dumps).
inline System fileSystem(TFile& f) {
  const std::string text = metaText(f, "meta/schema");
  if (text.find("\"system\":\"auau\"") != std::string::npos) return System::auau;
  if (text.find("\"system\":\"pp\"") != std::string::npos) return System::pp;
  throw std::runtime_error(std::string("unknown system in ") + f.GetName());
}

inline System fileSystem(const std::string& path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  return fileSystem(*f);
}

// True when the file is simulation: meta/source dataset isSim, isSimInclusive, isSimEmbedded or isSimEmbeddedInclusive.
inline bool isSimulation(TFile& f) {
  const std::string source = metaText(f, "meta/source");
  for (const char* dataset : {"isSim", "isSimInclusive", "isSimEmbedded", "isSimEmbeddedInclusive"})
    if (source.find(std::string("\"dataset\":\"") + dataset + "\"") != std::string::npos) return true;
  return false;
}

// The identity of a file's production input: meta/source input_file_sha256.
inline std::string sourceSha(TFile& f) {
  const std::string source = metaText(f, "meta/source");
  // Consolidated/augmented versions share one original logical source identity.
  // Their container hash changes with augmentation and must not become a new
  // normalization or resampling unit. Original chunk hashes remain in provenance.
  const std::string logical = "\"logical_source_id\":\"";
  const std::string key = source.find(logical) != std::string::npos ? logical : "\"input_file_sha256\":\"";
  const size_t at = source.find(key);
  if (at == std::string::npos) throw std::runtime_error(std::string("no input_file_sha256 in meta/source of ") + f.GetName());
  const size_t begin = at + key.size(), end = source.find('"', begin);
  return source.substr(begin, end - begin);
}

// {"path":...,"input_file_sha256":...} of one input. Opening the file moves the current ROOT directory, so outputs
// call this before they open their own file.
inline std::string sourceRecord(const std::string& path) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  return "{\"path\":\"" + path + "\",\"input_file_sha256\":\"" + sourceSha(*file) + "\"}";
}

// The JSON list of sourceRecord for every input.
inline std::string sourceRecords(const std::vector<std::string>& paths) {
  std::string out = "[";
  for (size_t i = 0; i < paths.size(); ++i) out += (i ? "," : "") + sourceRecord(paths[i]);
  return out + "]";
}

// A profile by name: "nominal", "historical", or a variant of one for the ABCD sideband systematic,
// "<base>:key=value[,key=value]", e.g. "nominal:nontight_min=0.2" or "nominal:sideband_gap=0". Keys:
//   nontight_min        Au+Au: the non-tight lower bound on the score (nominal 0.1); with augmented classes it is read
//                       from pho_id_class and the stored score, so the non-tight upper bound stays the package's
//   nontight_low_shift  pp: a shift of the non-tight lower score line
//   sideband_gap        GeV from the isolation threshold to the non-isolated sideband (pp nominal: the stored sideband);
//                       with augmented classes only 0, which joins the gap class to the non-isolated side
// The name, variant included, is written to every output, and FinalAnalysis requires every file of a result to carry
// the configured one.
inline Profile namedProfile(const std::string& name, System system) {
  const size_t colon = name.find(':');
  const std::string base = name.substr(0, colon);
  Profile p;
  if (base == "nominal") {
    p = nominal(system);
  } else if (base == "historical") {
    p = historical(system);
  } else {
    throw std::runtime_error("unknown profile " + name);
  }
  if (colon == std::string::npos) return p;
  const std::string variant = name.substr(colon + 1);
  std::vector<std::string> seen;
  size_t start = 0;
  while (true) {
    const size_t comma = variant.find(',', start);
    const std::string item = variant.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
    const size_t at = item.find('=');
    if (at == std::string::npos) throw std::runtime_error("a profile variant item is key=value: '" + item + "' in " + name);
    const std::string key = item.substr(0, at), value = item.substr(at + 1);
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) throw std::runtime_error("a profile variant key appears twice: " + name);
    seen.push_back(key);
    size_t used = 0;
    double x = 0.0;
    try {
      x = std::stod(value, &used);
    } catch (const std::exception&) {
      used = 0;
    }
    if (used == 0 || used != value.size() || !std::isfinite(x))
      throw std::runtime_error("a profile variant value is one finite number: " + item);
    if (key == "nontight_min") {
      if (system != System::auau) throw std::runtime_error("nontight_min is the Au+Au bound; pp shifts its line (nontight_low_shift)");
      if (!(0.0 <= x && x < 1.0)) throw std::runtime_error("nontight_min is a score in [0, 1): " + item);
      p.nonTightMin = x;
    } else if (key == "nontight_low_shift") {
      if (system != System::pp) throw std::runtime_error("nontight_low_shift is the pp line shift; Au+Au sets nontight_min");
      if (!(std::fabs(x) < 1.0)) throw std::runtime_error("nontight_low_shift is a score shift below 1 in size: " + item);
      p.nonTightLowShift = x;
    } else if (key == "sideband_gap") {
      if (x < 0.0) throw std::runtime_error("sideband_gap is not negative: " + item);
      if (p.augmentedClasses && x != 0.0)
        throw std::runtime_error("with augmented isolation classes the sideband gap can only be 0; another gap needs its own isolation package");
      p.sidebandGap = x;
    } else {
      throw std::runtime_error("unknown profile variant key " + key + " in " + name);
    }
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  p.name = name;
  return p;
}

}  // namespace photonjet
