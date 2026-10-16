#ifndef PHOTONJET_MEASUREMENT_REPLAY_H
#define PHOTONJET_MEASUREMENT_REPLAY_H
#include "RenderMetadata.h"
#include "../../TreeToHists/src/SelectionReplay.h"
#include "../../TreeToHists/src/ResponseOnlyResolution.h"
#include "../../TreeToHists/src/CentralityReplay.h"
#include <vector>

namespace fa {
/* External code used by the replay consumer, including the centrality
   header's own local dependencies. Saved analysis and pipeline receipts pin
   these alongside src/. Keep this list aligned with the compiler closure. */
inline std::vector<std::filesystem::path> measurementReplayCodeFiles() {
  const auto source = std::filesystem::path(__FILE__).parent_path();
  std::vector<std::filesystem::path> files;
  for (const char *name : {"SelectionReplay.h", "ResponseOnlyResolution.h", "CentralityReplay.h",
                          "CentralityWeightSource.h", "Configuration.h"})
    files.push_back((source / "../../TreeToHists/src" / name).lexically_normal());
  return files;
}
/* THE482/owner368: this response-only prescription varies the simulated
   sideband transfer with the response. The legacy template instead uses the
   nominal ABCD leakage as its normalization and has not been qualified. */
inline void validateMeasurementReplayConsumerScope(const nlohmann::json& recipe,
    const std::string& combinatoricMethod,const std::string& combinatoricShape) {
  if(recipe.is_null() || (!recipe.contains("response_only_resolution") && !recipe.contains("centrality_weight")))return;
  if(combinatoricMethod!="data_dphi" || combinatoricShape!="simulation")
    throw std::runtime_error("response-only resolution consumer requires Combinatoric.Method=data_dphi and Combinatoric.Shape=simulation; other correction modes are unqualified");
}
/* A saved-record consumer preserves the original producer. Its separate
   recipe and per-sample effect describe every varied measurement input. */
inline nlohmann::json measurementReplayRecipe(const nlohmann::json& meta) {
  using Json=nlohmann::json;
  const bool declared=meta.at("histograms").contains("measurement_replay");
  if(declared!=meta.contains("replay_provenance"))
    throw std::runtime_error("measurement replay declaration/provenance differs");
  if(!declared)return Json();
  const auto& p=meta.at("replay_provenance");
  const auto& r=meta.at("histograms").at("measurement_replay");
  const std::string key=meta.at("system").get<std::string>()+"/"+meta.at("name").get<std::string>();
  const std::string parent=p.at("parent_package_sha256");
  if(p.at("schema")!="PhotonJetMeasurementReplayV1" || p.at("scope")!="measurement_response_purity_inputs" ||
     p.at("nominal_diagnostics")!="parent_only" || p.at("normalization_applied")!=false ||
     p.at("candidate_archive")!="parent_only" || r.at("candidate_archive")!="parent_only" ||
     meta.contains("candidate_records") || !meta.contains("parent_candidate_records") ||
     r.at("schema")!="MeasurementReplayRecipeV1" || r.at("admission")!="MECHANICAL_FIXTURE_ONLY" ||
     p.at("recipe")!=r || p.at("recipe_sha256")!=photonjet_render::sha256(r.dump()) ||
     p.at("consumer_code")!=r.at("consumer_code") || r.at("consumer_code").empty() ||
     p.at("applied_effect")!=r.at("sample_effects").at(key) ||
     p.at("original_producer_code")!=meta.at("render_selection").at("producer_code") ||
     parent.size()!=64 || parent.find_first_not_of("0123456789abcdef")!=std::string::npos)
    throw std::runtime_error("measurement replay recipe, lineage or producer binding differs");
  photonjet_response_only::validate(r,key,meta.at("family"));
  photonjet_centrality_replay::validate(r,key,meta.at("family"));
  const auto centralityFreezes=photonjet_centrality_replay::admit(r,meta.at("selection"),
      [](const std::string& path){return photonjet_render::sha256(path,true);});
  photonjet_centrality_replay::validatePackage(r,meta,centralityFreezes);
  if(r.contains("centrality_weight") && r.at("centrality_weight").at("parent_packages").at(key)!=parent)
    throw std::runtime_error("centrality replay parent package identity differs");
  if(r.contains("selection_change")) {
    photonjet_selection_replay::validate(r,meta.at("parent_selection"),meta.at("selection"),
      [](const Json& j){return photonjet_render::sha256(j.dump());});
    if(meta.at("render_selection").at("selection_sha256")!=meta.at("selection").at("selection_sha256") ||
       meta.at("render_selection").at("photon_id")!=meta.at("selection").at("systems").at(meta.at("system")=="PP"?"pp":"auau"))
      throw std::runtime_error("selection replay render binding differs");
  } else if(meta.contains("parent_selection"))throw std::runtime_error("undeclared original selection");
  return r;
}
}
#endif
