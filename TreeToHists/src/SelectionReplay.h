#ifndef PHOTONJET_SELECTION_REPLAY_H
#define PHOTONJET_SELECTION_REPLAY_H

#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <string>

/* A changed cut uses the same retained score, witness and event population.
   Changes to models, event cuts, weights or NPB require their own qualification.
   This contract is shared by publication and the analysis consumer. */
namespace photonjet_selection_replay {
using Json=nlohmann::json;
inline void require(bool test,const char* message) {
  if(!test)throw std::runtime_error(std::string("selection replay: ")+message);
}
inline bool digestString(const Json& value) {
  if(!value.is_string())return false;
  const auto s=value.get<std::string>();
  return s.size()==64 && s.find_first_not_of("0123456789abcdef")==std::string::npos;
}
inline Json identityBody(Json selection) {
  selection.erase("selection_sha256");selection.erase("profile");selection.erase("name");
  return selection;
}
template<class Digest> inline void identity(const Json& selection,Digest digest) {
  require(selection.at("schema")=="PhotonJetSelectionV1" && selection.at("version")==1,"unknown selection schema");
  require(digestString(selection.at("selection_sha256")) && digest(identityBody(selection))==selection.at("selection_sha256").template get<std::string>(),"selection digest differs");
  require(selection.at("profile")=="selection:"+selection.at("selection_sha256").get<std::string>(),"selection profile differs");
}
inline double line(const Json& value,double x) {
  require(value.is_object() && value.size()==2 && value.contains("intercept") && value.contains("slope"),"invalid cut line");
  const double a=value.at("intercept"),b=value.at("slope");
  require(std::isfinite(a) && std::isfinite(b),"nonfinite cut line");
  return a+b*x;
}
inline Json invariantBody(Json selection) {
  selection=identityBody(selection);
  for(const char* name:{"pp","auau"}) {
    auto& s=selection.at("systems").at(name);
    for(const char* key:{"tight","non_tight_low","non_tight_high","verify_nominal_id","verify_nominal_isolation"})s.erase(key);
    auto& iso=s.at("isolation");iso.erase("threshold");iso.erase("sideband_gap_gev");
    if(s.contains("optional_model")) {
      auto& m=s.at("optional_model");m.erase("working_point_choice");
      for(auto& cell:m.at("cells")){cell.erase("tight");cell.erase("non_tight_high");}
      if(m.contains("native_evaluation"))for(auto& cell:m.at("native_evaluation").at("cells")){
        cell.erase("tight");cell.erase("non_tight_high");
      }
    }
  }
  return selection;
}
inline void bands(const Json& s,const Json& tight,const Json& high,double lo,double hi) {
  require(std::isfinite(lo) && std::isfinite(hi) && lo<hi,"invalid cut domain");
  for(double x:{lo,hi}) {
    const double l=line(s.at("non_tight_low"),x),h=line(high,x),t=line(tight,x);
    require(0.<=l && l<h && h<=t && t<=1.,"ID bands overlap or leave score domain");
  }
}
inline void cuts(const Json& selection) {
  for(const char* name:{"pp","auau"}) {
    const auto& s=selection.at("systems").at(name);
    require(s.at("verify_nominal_id")==false && s.at("verify_nominal_isolation")==false,"varied cuts cannot assert nominal producer-label equality");
    const bool pp=std::string(name)=="pp";
    bands(s,s.at("tight"),s.at("non_tight_high"),pp?15.:0.,pp?35.:80.);
    if(s.contains("optional_model")) {
      const auto& m=s.at("optional_model");
      const auto& choice=m.at("working_point_choice");
      require(choice.at("mode")=="explicit" && choice.at("tight")==s.at("tight") &&
        choice.at("non_tight_upper")==s.at("non_tight_high") &&
        choice.at("non_tight_lower")==s.at("non_tight_low").at("intercept") &&
        s.at("non_tight_low").at("slope")==0.,"explicit working-point description differs");
      for(size_t i=0;i<m.at("cells").size();++i) {
        const auto& cell=m.at("cells").at(i);
        const auto& bounds=cell.at(pp?"et":"centrality");
        bands(s,cell.at("tight"),cell.at("non_tight_high"),bounds.at(0),bounds.at(1));
        require(cell.at("tight")==s.at("tight") && cell.at("non_tight_high")==s.at("non_tight_high"),"cell and explicit cut lines differ");
        if(m.contains("native_evaluation")) {
          const auto& native=m.at("native_evaluation").at("cells").at(i);
          require(native.at("id")==cell.at("id") && native.at("tight")==cell.at("tight") &&
            native.at("non_tight_high")==cell.at("non_tight_high"),"native model routing cut lines differ");
        }
      }
    }
    const auto& iso=s.at("isolation");
    require(iso.at("radius")== (pp?.4:.3) && iso.at("subtraction_method")== (pp?1:2),"unsupported or unsubtracted isolation");
    if(!iso.at("threshold").is_string()) {
      line(iso.at("threshold"),pp?15.:0.);line(iso.at("threshold"),pp?35.:80.);
    } else require(iso.at("threshold")=="stored" && pp,"unqualified stored isolation cut");
    if(!iso.at("sideband_gap_gev").is_null()) {
      const double gap=iso.at("sideband_gap_gev");require(std::isfinite(gap) && gap>=0,"invalid isolation gap");
    }
  }
}
template<class Digest> inline void validate(const Json& recipe,const Json& parent,const Json& effective,Digest digest) {
  if(!recipe.contains("selection_change")){require(parent==effective,"undeclared selection change");return;}
  identity(parent,digest);identity(effective,digest);
  const auto& change=recipe.at("selection_change");
  require(change.at("schema")=="PhotonJetSelectionReplayV1" &&
    change.at("parent_selection_sha256")==parent.at("selection_sha256") &&
    change.at("effective_selection_sha256")==effective.at("selection_sha256"),"selection lineage differs");
  require(change.at("correction_policy")=="frozen_parent_weights_mechanical_only" &&
    change.at("score_policy")=="retain_exact_selected_score" &&
    recipe.at("admission")=="MECHANICAL_FIXTURE_ONLY", "unqualified corrections or score policy");
  require(recipe.at("cone_prescription")=="retained-subtracted-witness-selection-only","unknown isolation witness prescription");
  for(const auto& effect:recipe.at("sample_effects"))require(effect=="","combined energy and cut variations are not qualified");
  require(invariantBody(parent)==invariantBody(effective),"changed model, event selection, NPB, weight or other invariant");
  require(change.at("definition_sources").is_array() && !change.at("definition_sources").empty(),"missing cut-definition source");
  for(const auto& source:change.at("definition_sources")) {
    require(digestString(source.at("sha256")) && source.at("path").is_string() &&
      !source.at("path").get<std::string>().empty() && source.at("locator").is_string() &&
      !source.at("locator").get<std::string>().empty(),"invalid cut-definition source");
  }
  cuts(effective);
}
}
#endif
