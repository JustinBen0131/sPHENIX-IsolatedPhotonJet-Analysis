/* A collection explicitly binds one simulation response and, optionally, one system's DATA scale recipe.
   This changes no package metadata and grants no scientific acceptance to that recipe.
   TreeToHists assembly and FinalAnalysis ingestion use this same contract. */
#ifndef PHOTONJET_COLLECTION_CONTRACT_H
#define PHOTONJET_COLLECTION_CONTRACT_H
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

namespace photonjet_collection {
using Json = nlohmann::json;
inline Json simulationResponse(const std::string &recipe) {
  if (recipe.empty() || recipe.find_first_of("\r\n\t") != std::string::npos)
    throw std::runtime_error("simulation response recipe must be explicit and nonempty");
  return {{"schema", "PhotonJetSimulationResponseCollectionV1"},
          {"data_response_variation", ""}, {"simulation_response_variation", recipe}};
}
inline bool keys(const Json &value, std::initializer_list<const char*> names) {
  if (!value.is_object() || value.size()!=names.size()) return false;
  for (const auto *name:names) if (!value.contains(name)) return false;
  return true;
}
inline bool positive(const Json &value) {
  return value.is_number() && std::isfinite(value.get<double>()) && value.get<double>()>0.;
}
inline bool finite(const Json &value) {
  return value.is_number() && std::isfinite(value.get<double>());
}
inline std::map<std::string,std::string> variationItems(const std::string &recipe) {
  std::map<std::string,std::string> out; std::istringstream input(recipe);
  for (std::string item;input>>item;) {
    const auto at=item.find('=');
    if (at==std::string::npos || at==0 || at+1==item.size() ||
        !out.emplace(item.substr(0,at),item.substr(at+1)).second)
      throw std::runtime_error("invalid or duplicate collection recipe item");
  }
  return out;
}
inline bool emsRecipe(const std::string &recipe) {
  const auto items=variationItems(recipe); const auto mode=items.find("photon_selection");
  return mode!=items.end() && mode->second=="ems_selection_response";
}
inline void requireEmsRecipe(const std::string &recipe) {
  const auto items=variationItems(recipe);
  if (!emsRecipe(recipe) || !items.count("photon_resolution") || !items.count("jet_resolution"))
    throw std::runtime_error("EMS requires its exact selection token and named resolutions");
  for (const auto &[key,value]:items) {
    bool valid=false;
    if (key=="photon_selection") valid=value=="ems_selection_response";
    else if (key=="photon_resolution") valid=value=="ppg18_nominal" || value=="ppg18_high" || value=="ppg18_low";
    else if (key=="jet_resolution") valid=value=="sphenix_pp2026_nominal" || value=="sphenix_pp2026_high" || value=="sphenix_pp2026_low";
    else if (key=="centrality_weight") valid=value=="on" || value=="off";
    else if (key=="seed") {
      valid=!value.empty() && value.find_first_not_of("0123456789")==std::string::npos;
      try {std::stoull(value);} catch (...) {valid=false;}
    }
    else if (key=="em_scale" || key=="photon_scale" || key=="jet_scale" || key=="photon_smear" || key=="jet_smear") {
      size_t used=0; double number=0.; try {number=std::stod(value,&used);} catch (...) {}
      valid=used==value.size() && std::isfinite(number) &&
        (key=="em_scale" ? number>0. : (key=="photon_scale" || key=="jet_scale" ? number==1. : number==0.));
    }
    if (!valid) throw std::runtime_error("invalid EMS collection recipe item: "+key);
  }
}
inline void requireEmsContract(const Json &ems) {
  const Json cone={{"rule","retain_captured_analysis_witness"},{"status","MECHANICAL_UNADOPTED"},
                   {"raw_auau_isolation","ARCHIVAL_ONLY"}};
  if (!keys(ems,{"schema","source","cone","selection","models","pp_witness_line"}) ||
      ems.at("schema")!="EmsNativeSelectionContractV1" ||
      ems.at("source")!=Json{{"sha256","cabc3a5eb081ce77e0e1d83936e04bcc427bed6f5be5066183983edd4dc32105"}} ||
      ems.at("cone")!=cone || !ems.at("selection").is_object() || !keys(ems.at("models"),{"pp","npb"}) ||
      !keys(ems.at("pp_witness_line"),{"intercept","slope","gap"}))
    throw std::runtime_error("invalid source-bound EMS collection contract");
  for (const auto *name:{"pp","npb"}) {
    const auto &model=ems.at("models").at(name);
    if (!keys(model,{"file","key","sha256"}) || !model.at("file").is_string() || model.at("file")=="" ||
        !model.at("key").is_string() || model.at("key")=="" || !model.at("sha256").is_string())
      throw std::runtime_error("invalid EMS model binding");
    const auto pin=model.at("sha256").get<std::string>();
    if (pin.size()!=64 || pin.find_first_not_of("0123456789abcdef")!=std::string::npos)
      throw std::runtime_error("invalid EMS model hash");
  }
  const auto &line=ems.at("pp_witness_line");
  if (!finite(line.at("intercept")) || !finite(line.at("slope")) || !positive(line.at("gap")))
    throw std::runtime_error("invalid EMS witness line");
}
inline void requireDataScale(const Json &scale) {
  if (!keys(scale,{"system","nominal","varied"}) ||
      (scale.at("system")!="pp" && scale.at("system")!="auau") ||
      !positive(scale.at("nominal")) || !positive(scale.at("varied")) || scale.at("nominal")==scale.at("varied"))
    throw std::runtime_error("DATA scale recipe requires one system and distinct positive finite nominal/varied values");
}
inline Json collectionRecipe(const std::string &recipe,const Json &ems=Json(),const Json &dataScale=Json()) {
  if (ems.is_null() && dataScale.is_null()) return simulationResponse(recipe);
  if (!recipe.empty()) simulationResponse(recipe);
  if (!ems.is_null()) {requireEmsRecipe(recipe);requireEmsContract(ems);}
  else if (emsRecipe(recipe)) throw std::runtime_error("EMS recipe requires an explicit EMS selection contract");
  if (!dataScale.is_null()) requireDataScale(dataScale);
  return {{"schema","PhotonJetRecipeCollectionV2"},{"data_response_variation",""},
          {"simulation_response_variation",recipe},{"ems_selection_contract",ems},{"data_jet_scale",dataScale}};
}
inline bool recipeV2(const Json &contract) {
  return contract.is_object() && contract.value("schema","")=="PhotonJetRecipeCollectionV2";
}
inline std::string dataScaleSystem(const Json &contract) {
  return recipeV2(contract) && !contract.at("data_jet_scale").is_null()
    ? contract.at("data_jet_scale").at("system").get<std::string>() : "";
}
inline void requireContract(const Json &contract, const Json &reference) {
  if (contract.is_null()) return;
  if (recipeV2(contract)) {
    if (!keys(contract,{"schema","data_response_variation","simulation_response_variation","ems_selection_contract","data_jet_scale"}) ||
        !contract.at("simulation_response_variation").is_string() ||
        contract!=collectionRecipe(contract.at("simulation_response_variation"),contract.at("ems_selection_contract"),contract.at("data_jet_scale")) ||
        reference.value("response_variation",Json())!="" || reference.contains("measurement_replay") ||
        reference.contains("ems_selection_contract") || reference.contains("ems_qa_extensions") || reference.contains("ems_candidate_archive"))
      throw std::runtime_error("invalid composed collection contract");
    const auto system=dataScaleSystem(contract);
    if (!system.empty()) {
      const auto &scales=reference.at("data_jet_scale");
      if (!keys(scales,{"pp","auau"}) || !positive(scales.at("pp")) || !positive(scales.at("auau")) ||
          scales.at(system)!=contract.at("data_jet_scale").at("varied"))
        throw std::runtime_error("collection reference must declare the effective varied DATA scale");
    }
    return;
  }
  if (!contract.is_object() || !contract.contains("simulation_response_variation") ||
      !contract.at("simulation_response_variation").is_string() ||
      contract != simulationResponse(contract.at("simulation_response_variation").get<std::string>()) ||
      !reference.contains("response_variation") || reference.at("response_variation") != "" ||
      reference.contains("measurement_replay"))
    throw std::runtime_error("invalid simulation response collection contract");
}
inline void requireSelection(const Json &contract,const Json &selection) {
  if (recipeV2(contract) && !contract.at("ems_selection_contract").is_null() &&
      contract.at("ems_selection_contract").at("selection")!=selection)
    throw std::runtime_error("EMS collection selection differs");
}
inline bool histogramsAgree(const Json &filled, const Json &reference, const Json &contract,
                            const std::string &family,const std::string &system="") {
  requireContract(contract, reference);
  if (contract.is_null()) return filled == reference;
  Json expected = reference;
  if (recipeV2(contract)) {
    const bool simulation=family=="photonjet" || family=="inclusivejet";
    if ((!simulation && family!="data") || (system!="PP" && system!="AuAu")) return false;
    if (simulation && !contract.at("simulation_response_variation").get<std::string>().empty()) {
      const auto items=variationItems(contract.at("simulation_response_variation"));
      if (items.count("centrality_weight") && items.at("centrality_weight")=="off" && system!="AuAu") return false;
    }
    const auto variedSystem=dataScaleSystem(contract);
    if (!variedSystem.empty()) {
      if (system!="PP" && system!="AuAu") return false;
      const std::string actual=system=="PP" ? "pp" : "auau";
      if (simulation || actual!=variedSystem)
        expected["data_jet_scale"][variedSystem]=contract.at("data_jet_scale").at("nominal");
    }
    if (simulation) {
      expected["response_variation"]=contract.at("simulation_response_variation");
      const auto &ems=contract.at("ems_selection_contract");
      if (!ems.is_null()) {
        expected["ems_selection_contract"]=ems;
        expected["ems_qa_extensions"]="nominal_archival_only";
        expected["ems_candidate_archive"]="original_nominal_inputs_replay_unqualified";
      }
    }
    return filled==expected;
  }
  if (family == "photonjet" || family == "inclusivejet")
    expected["response_variation"] = contract.at("simulation_response_variation");
  else if (family != "data") return false;
  return filled == expected;
}
}  // namespace photonjet_collection
#endif
