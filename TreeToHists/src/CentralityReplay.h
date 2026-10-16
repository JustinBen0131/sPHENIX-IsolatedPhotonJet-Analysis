#pragma once
// Shared admission for the one centrality on/off comparison. The historic
// frozen target is a pinned mechanical fixture, not current-data acceptance.
#include "CentralityWeightSource.h"
#include "Configuration.h"
#include <TFile.h>
#include <TH3.h>
#include <TObjString.h>
#include <fstream>
#include <map>
#include <memory>
#include <set>

namespace photonjet_centrality_replay {
using Json = nlohmann::json;
inline void require(bool ok, const std::string& why) {
  if (!ok) throw std::runtime_error("centrality replay: " + why);
}
inline Json read(const std::string& path) {
  std::ifstream stream(path); require(bool(stream), "unreadable input " + path);
  Json value; stream >> value; return value;
}
inline Json object(TFile& file, const char* path) {
  const auto* value = file.Get<TObjString>(path);
  require(value, std::string("missing object ") + path);
  return Json::parse(value->GetString().Data());
}
inline bool validate(const Json& recipe, const std::string& sample, const std::string& family) {
  const auto& effects = recipe.at("sample_effects");
  if (!recipe.contains("centrality_weight")) {
    for (const auto& effect : effects)
      require(effect.get<std::string>().find("centrality_weight=") == std::string::npos,
              "centrality effect lacks its scope declaration");
    return false;
  }
  const auto& c = recipe.at("centrality_weight");
  require(!recipe.contains("selection_change") && !recipe.contains("response_only_resolution"), "mixed source scope");
  require(recipe.at("admission") == "MECHANICAL_FIXTURE_ONLY" &&
          recipe.at("cone_prescription") == "unchanged-nominal-measurement-and-isolation", "unqualified admission or cone policy");
  require(c.at("schema") == "PhotonJetCentralityWeightReplayV1" &&
          c.at("collection_subset") == "all_auau_parent_samples" &&
          c.at("measurement") == "nominal" && c.at("classification") == "original_unsmeared" &&
          c.at("response") == "all_native_accumulators" && c.at("vertex") == "nominal" &&
          c.at("source_normalization") == "nominal_once", "altered centrality source policy");
  const std::string mode = c.at("mode");
  require(mode == "on" || mode == "off", "unknown centrality mode");
  const auto& roster = c.at("sample_families");
  require(roster.is_object() && !roster.empty() && roster.size() == effects.size(), "incomplete sample effects");
  require(roster.contains(sample) && roster.at(sample) == family, "sample is outside the declared family roster");
  std::set<std::string> affected, expected;
  for (auto it = roster.begin(); it != roster.end(); ++it) {
    require(it.key().rfind("AuAu/", 0) == 0, "centrality comparison is AuAu-only");
    const std::string f = it.value();
    require(f == "data" || f == "photonjet" || f == "inclusivejet", "unknown AuAu family");
    if (f != "data") expected.insert(it.key());
    require(effects.contains(it.key()) && effects.at(it.key()) == (f == "data" ? "" : "centrality_weight=" + mode),
            "undeclared, missing or mixed energy effect");
  }
  for (const auto& value : c.at("affected_samples"))
    require(affected.insert(value.get<std::string>()).second, "duplicate affected sample");
  require(!expected.empty() && affected == expected, "every SIM response component must vary together");
  return mode == "off" && family != "data";
}

// Source identities and counts are checked before the producer creates output,
// and again by the actual consumer. Paths may move; the frozen bytes may not.
template<class Hash>
inline Json admit(const Json& recipe, const Json& selection, Hash hash) {
  if (!recipe.contains("centrality_weight")) return Json();
  const auto& c = recipe.at("centrality_weight");
  const std::map<std::string, std::string> expected = {
    {"photonjet", "48f10242d219be747a2808f78fb577a0da87ac4581687f3ce26a87d9111e4de2"},
    {"inclusivejet", "bba92948b34a09d45b2109f76650e67ffe0b6c5b926b550c53233d351b600856"},
    {"qualification", "3d7f8ca875f098b5a41d52529aca914c90b8b149ce6a396779622a295d511633"},
    {"resolution", "11d82d834a25615350b2bc504c57bf447f98679bf17548d7d8ae38fbbd7d503e"},
    {"data", "9d1b62fe92c90313fcf5dfe8252dea4e9c1d4f14c747ae09254282ad29f2cf5c"},
    {"vertex", "ffaeaf54ab4cccf6e044e6cd2b77d347431dd9c5652038b578b679be44083be5"}
  };
  const auto& pins = c.at("inputs");
  require(pins.size() == expected.size(), "missing or extra frozen inputs");
  for (const auto& [name, digest] : expected) {
    const auto& pin = pins.at(name);
    require(pin.at("sha256") == digest && hash(pin.at("path").template get<std::string>()) == digest,
            "frozen input bytes differ: " + name);
  }
  const auto load = [&](const char* name) { return read(pins.at(name).at("path")); };
  const auto photon = load("photonjet"), inclusive = load("inclusivejet"), q = load("qualification");
  const auto selected = load("resolution").at("site_selection");
  for (const char* key : {"schema", "version", "preselection", "timing_reference"})
    require(selected.at(key) == selection.at(key), std::string("source selection differs: ") + key);
  auto sourceAuAu=selection.at("systems").at("auau");
  const auto& frozenAuAu=selected.at("systems").at("auau");
  // The fixture has a local copy of the same model. Permit only that exact
  // hash-verified file relocation; every classifier/domain/WP field remains.
  auto& sourceCells=sourceAuAu.at("optional_model").at("native_evaluation").at("cells");
  const auto& frozenCells=frozenAuAu.at("optional_model").at("native_evaluation").at("cells");
  require(sourceCells.size()==frozenCells.size(), "source model cell roster differs");
  for(size_t i=0;i<sourceCells.size();++i) {
    const std::string path=sourceCells.at(i).at("file");
    require(sourceCells.at(i).at("sha256")==frozenCells.at(i).at("sha256") &&
            hash(path)==frozenCells.at(i).at("sha256").template get<std::string>(), "source model bytes differ");
    sourceCells.at(i)["file"]=frozenCells.at(i).at("file");
  }
  require(frozenAuAu == sourceAuAu, "source AuAu selection differs");
  std::unique_ptr<TFile> data(TFile::Open(pins.at("data").at("path").template get<std::string>().c_str(), "READ"));
  require(data && !data->IsZombie(), "unreadable frozen DATA");
  const auto* h = data->Get<TH3>("qa/slices/data/width_tight_iso_weta_cogx");
  require(h && h->GetNbinsZ() == 16, "missing frozen target histogram");
  const int first = h->GetYaxis()->FindBin(15. + 1e-9), last = h->GetYaxis()->FindBin(35. - 1e-9);
  require(h->GetYaxis()->GetBinLowEdge(first) == 15. && h->GetYaxis()->GetBinUpEdge(last) == 35., "target ET edges differ");
  Json counts = Json::array(); double outside = 0.;
  for (int z = 1; z <= 16; ++z) {
    double n = 0.;
    for (int y = 0; y <= h->GetNbinsY() + 1; ++y) for (int x = 0; x <= h->GetNbinsX() + 1; ++x)
      if (y >= first && y <= last) n += h->GetBinContent(x, y, z); else outside += h->GetBinContent(x, y, z);
    counts.push_back(n);
  }
  require(outside == 0. && c.at("target_counts") == counts, "selected target counts differ");
  photonjet::centrality_source::validateTarget(photon, inclusive, selected,
      object(*data, "meta/package").at("selection"), counts, q);
  const auto& parentPin = c.at("parent_collection");
  require(hash(parentPin.at("path").template get<std::string>()) == parentPin.at("sha256").template get<std::string>(), "parent collection bytes differ");
  const auto parent = photonjet_config::readYaml(parentPin.at("path"));
  require(parent.at("schema") == "PhotonJetHistogramCollectionV1" && parent.at("status") == "VALIDATED" &&
          !parent.contains("measurement_replay") && !parent.contains("compact_of"), "invalid original parent collection");
  Json families = Json::object(), packages = Json::object();
  for (const auto& row : parent.at("samples")) if (row.at("system") == "AuAu") {
    const std::string key = row.at("key");
    require(!families.contains(key), "duplicate parent sample");
    families[key] = row.at("family"); packages[key] = row.at("sha256");
  }
  require(families == c.at("sample_families") && packages == c.at("parent_packages"), "incomplete AuAu parent roster");
  require(parent.at("selection_sha256") == selection.at("selection_sha256"), "source selection identity differs from parent");
  return {{"photonjet", photon}, {"inclusivejet", inclusive}};
}
inline void validatePackage(const Json& recipe, const Json& meta, const Json& freezes) {
  if (!recipe.contains("centrality_weight")) return;
  const std::string family = meta.at("family");
  require(meta.at("system") == "AuAu", "centrality package must be AuAu");
  if (family == "data") return;
  const auto& f = freezes.at(family);
  require(meta.at("centrality_weights") == f.at("centrality").at("weights"), "package frozen weight map differs");
  require(meta.at("vertex_reweighting").at("sha256") == f.at("vertex").at("sha256") &&
          meta.at("vertex_reweighting").at("coordinate") == "truth_vtx_z" &&
          meta.at("centrality_data_trigger") == f.at("data_trigger"), "package vertex or target trigger differs");
  require(meta.at("normalization").at("applied_in_histograms") == false, "source normalization already applied");
}
inline void recordedFactor(double percentile, double recorded, const Json& weights) {
  require(weights.is_array() && weights.size() == 16 && std::isfinite(percentile), "invalid recorded centrality domain");
  const double expected = percentile >= 0. && percentile < 80. ? weights.at(static_cast<size_t>(percentile / 5.)).get<double>() : 1.;
  require(std::isfinite(recorded) && recorded > 0. && recorded == expected, "recorded centrality factor differs from admitted map");
}
}
