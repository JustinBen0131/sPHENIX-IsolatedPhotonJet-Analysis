/* Save the selection, units, bins and calculation settings used by the measurement.
   Figures and tables read these values from the result, so their labels describe the saved calculation.
   The accompanying receipt binds this metadata to its ROOT product. */
#include "FinalAnalysis.h"
#include "RenderMetadata.h"
#include "MeasurementReplay.h"
#include "NamedResponse.h"
#include <set>
#include <memory>
#include <filesystem>
#include "TFile.h"
#include "TObjString.h"

namespace fa {
std::string renderMetadata(const Settings &s, bool allowMissingData) {
  using namespace photonjet_render;
  Json body = {{"schema", "PhotonJetRenderMetadataV1"}, {"packages", Json::object()}, {"results", Json::object()},
               {"analysis_program", Json::object()}};
  for (const auto &file : std::filesystem::directory_iterator(std::filesystem::path(__FILE__).parent_path()))
    if (file.is_regular_file() && (file.path().extension()==".cc" || file.path().extension()==".h"))
      body["analysis_program"][file.path().filename().string()] = sha256(file.path().string(), true);
  for (const auto &file : measurementReplayCodeFiles())
    body["analysis_program"][std::filesystem::relative(file, std::filesystem::path(__FILE__).parent_path()).generic_string()] =
        sha256(file.string(), true);
  body["analysis_program"]["../../TreeToHists/src/CollectionContract.h"]=sha256(named_response::codeFile().string(),true);
  std::set<std::string> seen;
  for (const auto &sys : s.systems) {
    Json classes = Json::array();
    for (const auto &other : s.systems) if (other.package == sys.package) {
      Json range = other.slices.empty() || other.slices.front()<0 ? Json() : Json{5*other.slices.front(), 5*(other.slices.back()+1)};
      classes.push_back({{"id", other.centrality}, {"centrality_percent", range}});
    }
    const Json centrality = sys.slices.empty() || sys.slices.front()<0 ? Json() : Json{5*sys.slices.front(), 5*(sys.slices.back()+1)};
    body["results"][sys.name] = {{"package", sys.package}, {"centrality", sys.centrality}, {"centrality_percent", centrality}, {"photon_edges_gev", sys.photonEdges},
      {"recoil_min_gev", s.jetThreshold}, {"reported_xj", {sys.reportedLow, s.reportedHigh}}, {"data_toy_model", s.dataToyModel},
      {"data_seed", s.seed}, {"response_seed", s.responseSeed}, {"closure_seed", s.closureSeed}, {"prior", s.priorForm},
      {"prior_fraction", s.priorFraction}, {"combinatoric_method", s.combinatoricMethod},
      {"closure_scope", "prior_fit_and_unfolding_after_corrections"},
      {"full_correction_chain_closure", "REQUIRES_INDEPENDENT_PHYSICS_QUALIFICATION"}};
    if(named_response::active(sys))body["results"][sys.name]["named_response_component"]=named_response::selection(sys);
    if (!seen.insert(sys.package).second) continue;
    Json samples = Json::object(), weightPackages = Json::array();
    std::set<std::string> names{"data"};
    for (const auto &other : s.systems) if (other.package == sys.package) names.insert(other.simulation.begin(), other.simulation.end());
    const auto inclusive = s.inclusive.find(sys.package);
    if (inclusive != s.inclusive.end()) names.insert(inclusive->second.begin(), inclusive->second.end());
    for (const auto &name : names) {
      const auto path = s.inputDirectory+"/"+sys.package+"/"+name+".root";
      if (allowMissingData && name == "data" && !std::filesystem::exists(path)) continue;
      std::unique_ptr<TFile> file(TFile::Open(path.c_str()));
      auto *saved = file && !file->IsZombie() ? dynamic_cast<TObjString *>(file->Get("meta/package")) : nullptr;
      if (!saved) throw std::runtime_error("render metadata needs sample metadata: " + path);
      const Json meta = Json::parse(saved->GetString().Data());
      if (meta.value("schema", "") != "PhotonJetSampleV1" || !meta.contains("render_selection") ||
          !meta.at("render_selection").contains("producer_code"))
        throw std::runtime_error("physics rendering requires a native histogram package: " + path);
      const std::string digest = sha256(path, true);  /* once: a combined DATA package is tens of GB */
      weightPackages.push_back({{"schema", meta.at("schema")}, {"name", sys.package+"/"+name},
        {"sha256", digest}, {"family", meta.at("family")}, {"source_count", meta.at("source_count")},
        {"selection_sha256", meta.at("selection").at("selection_sha256")},
        {"producer_sha256", sha256(meta.at("render_selection").at("producer_code").dump())},
        {"normalization", name == "data" ? Json() : meta.at("normalization")}});
      samples[name] = {{"path", path}, {"sha256", digest}, {"family", meta.at("family")},
                       {"normalization", meta.value("normalization", Json())}, {"availability", meta.value("availability", Json())}};
      if(!measurementReplayRecipe(meta).is_null())samples[name]["replay_provenance"]=meta.at("replay_provenance");
      if (name == "data" || !body["packages"].contains(sys.package)) {
        body["packages"][sys.package] = {{"classes", classes}, {"selection", meta.value("render_selection", Json())},
          {"annotation_status", meta.contains("render_selection") ? "COMPILED_SELECTION" : "LEGACY_UNVERIFIED_ANNOTATION"}};
      }
      if (meta.contains("render_selection") && meta.at("render_selection").at("selection_sha256") != meta.at("selection").at("selection_sha256"))
        throw std::runtime_error("render selection does not bind the filled selection: " + path);
    }
    Json weights = {{"schema", "PhotonJetHistogramPlotWeightsV1"}, {"role", "physics_event_source"},
      {"training_weights_used", false}, {"packages", weightPackages}};
    requirePhysicsPlotWeights(weights, "native histogram collection");
    body["packages"][sys.package]["weight_provenance"] = weights;
    body["packages"][sys.package]["samples"] = samples;
    auto &qa=body["packages"][sys.package]["selection"]["qa"];
    qa["abcd_edges_gev"]=s.qaAbcdEdges;
    qa["windows"]={{"shower_shapes_gev",s.qaShowerShapes},{"isolation_overlay_gev",s.qaIsolationOverlay},
      {"truth_isolation_gev",s.qaTruthIsolation},{"ncb_score_gev",s.qaNcbScore}};
  }
  return seal(body).dump();
}

void renderProductReceipt(const std::string &path, const std::string &metadata) {
  using namespace photonjet_render;
  const auto body = verify(Json::parse(metadata));
  for (const auto &package : body.at("packages")) for (const auto &sample : package.at("samples"))
    if (sha256(sample.at("path").get<std::string>(), true) != sample.at("sha256").get<std::string>())
      throw std::runtime_error("sample changed while producing analysis output");
  write(path+".render.json", {{"schema", "PhotonJetRenderProductV1"}, {"product_sha256", sha256(path, true)},
                            {"metadata_sha256", Json::parse(metadata).at("sha256")}});
}
}  /* namespace fa */
