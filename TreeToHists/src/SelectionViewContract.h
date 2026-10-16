/* Cut-only original-tree views. Existing arrays/models, event gates and weights
   are invariant. Explicit efficiency payloads are accepted only with a pinned
   derivation that binds this exact parent selection and supplies both lines. */
#pragma once
#include "PhotonJet.h"
#include "SelectionReplay.h"
#include "EnergyVariation.h"
#include "Configuration.h"
namespace photonjet { namespace selection_views {
using Json=nlohmann::json;
inline void require(bool ok,const char* why){if(!ok)throw std::runtime_error(std::string("native selection views: ")+why);}
/* Selection alternatives share either the historical nominal or exact
   corrected nominal response. Detector endpoints stay separate sinks. */
inline bool cutComposition(const response::Variation& v){
 if(response::variationJson(v)==response::variationJson(response::Variation{}))return true;
 return response::variationJson(v)==response::variationJson(response::parseVariation(
  "photon_resolution=ppg18_nominal photon_selection=response_only jet_resolution=sphenix_pp2026_nominal seed=480"));
}
inline void compatible(const Profile& parent,const Profile& effective,double expectedStoredGap=std::numeric_limits<double>::quiet_NaN()){
  require(parent.configured && effective.configured,"resolved selections required");
  require(parent.system==effective.system,"collision system differs");
  if(std::isnan(parent.sidebandGap) && std::isfinite(effective.sidebandGap))
    require(std::isfinite(expectedStoredGap) && expectedStoredGap>=0.,"stored nominal gap requires an explicit checked witness value");
  const auto invariants=[](const Profile& p){return std::tie(p.photonEtMin,p.photonEtMax,p.photonAbsEtaMax,p.vertexAbsZMax,
    p.jetRadius,p.jetPtMin,p.jetAbsEtaMax,p.overlapDeltaR,p.deltaPhiMin,p.augmentedClasses,p.configured,p.ppIsolationParity,
    p.optionalID,p.nativeID,p.scoreBranch,p.stateBranch,p.isolationRadiusIndex,p.isolationSubtraction,p.isolationCentMin,p.isolationCentMax);};
  require(invariants(parent)==invariants(effective),"runtime profile invariants differ");
  const auto a=Json::parse(parent.selectionJson),b=Json::parse(effective.selectionJson);
  require(photonjet_selection_replay::invariantBody(a)==photonjet_selection_replay::invariantBody(b),"model, score, preselection, population or weight invariant changed");
  require(!effective.verifyNominalID && !effective.verifyNominalIsolation,"varied profile cannot verify nominal labels");
  photonjet_selection_replay::cuts(b);
  const auto& spec=b.at("systems").at(effective.system==System::pp?"pp":"auau");
  const auto same=[](const Line& l,const Json& j){return l.intercept==j.at("intercept").get<double>() && l.slope==j.at("slope").get<double>();};
  require(effective.selectionId==b.at("selection_sha256").get<std::string>() && effective.name==b.at("profile").get<std::string>(),"profile identity differs from metadata");
  require(same(effective.tight,spec.at("tight")) && same(effective.nonTightLow,spec.at("non_tight_low")) && same(effective.nonTightHigh,spec.at("non_tight_high")),"runtime ID cuts differ from metadata");
  const auto& iso=spec.at("isolation");
  require(effective.isolationFromWitness==iso.at("threshold").is_string(),"runtime isolation mode differs from metadata");
  if(!effective.isolationFromWitness)require(same(effective.isolationThreshold,iso.at("threshold")),"runtime isolation line differs from metadata");
  require(iso.at("sideband_gap_gev").is_null()?std::isnan(effective.sidebandGap):effective.sidebandGap==iso.at("sideband_gap_gev").get<double>(),"runtime gap differs from metadata");
  if(effective.optionalID){
    const auto& cells=spec.at("optional_model").at("cells");require(cells.size()==effective.idCells.size(),"runtime cell count differs");
    for(size_t i=0;i<cells.size();++i){const auto& c=cells.at(i);const auto& p=effective.idCells[i];
      require(same(p.tight,c.at("tight"))&&same(p.nonTightHigh,c.at("non_tight_high")) &&
        p.etMin==c.at("et").at(0).get<double>() && p.etMax==c.at("et").at(1).get<double>() &&
        p.centMin==c.at("centrality").at(0).get<double>() && p.centMax==c.at("centrality").at(1).get<double>(),"runtime cell cuts/ranges differ from metadata");
    }
  }
}
/* A null parent gap preserves per-photon sideband-threshold, not an
   assumed0.8GeV. A changed absolute gap must prove that declared nominal on
   every valid retained witness before filling the changed view. */
template<class Values,class Ints> void validateStoredGap(const Profile& parent,double expected,size_t n,
 const Values& cone,const Values& threshold,const Values& sideband,const Ints& method,const Ints& role){
 if(std::isnan(expected))return;
 require(std::isfinite(expected)&&expected>=0.,"invalid stored-gap declaration");
 require(cone.GetSize()==2*n && threshold.GetSize()==2*n && sideband.GetSize()==2*n &&
  method.GetSize()==2*n && role.GetSize()==2*n,"stored-gap witness arrays differ");
 for(size_t i=0;i<n;++i){const auto k=2*i+parent.isolationRadiusIndex;
  if(role[k]!=0 || method[k]!=parent.isolationSubtraction || !std::isfinite(cone[k]) || !std::isfinite(threshold[k]) || !std::isfinite(sideband[k]))continue;
  require(std::abs(sideband[k]-threshold[k]-expected)<=1.e-12*std::max(1.,std::abs(expected)),"source's nominal isolation gap differs from declaration");
 }
}
template<class Digest,class FileDigest> Json resolve(const Json& parent,const Json& payload,Digest digest,FileDigest fileDigest){
  namespace replay=photonjet_selection_replay;
  replay::identity(parent,digest);
  require(payload.at("schema")=="PhotonJetNativeSelectionViewV1","unknown payload schema");
  require(payload.at("parent_selection_sha256")==parent.at("selection_sha256"),"parent selection differs");
  const std::string system=payload.at("system"),kind=payload.at("kind");
  require(system=="pp" || system=="auau","unknown system");
  const auto& source=payload.at("definition_source");
  require(replay::digestString(source.at("sha256")) && source.at("locator").is_string() && !source.at("locator").get<std::string>().empty(),"source pin/locator absent");
  require(fileDigest(source.at("path").template get<std::string>())==source.at("sha256").template get<std::string>(),"definition source bytes differ");
  Json result=parent;auto& s=result["systems"][system];
  if(kind=="isolation_efficiency" || kind=="tight_efficiency"){
    const auto& binding=payload.at("efficiency_binding");
    std::ifstream derivation(source.at("path").template get<std::string>());
    Json derivationBody;derivation>>derivationBody;
    require(derivationBody.at(Json::json_pointer(source.at("locator").template get<std::string>()))==binding,"efficiency binding is absent from pinned derivation");
    require(binding.at("parent_selection_sha256")==parent.at("selection_sha256"),"efficiency derivation is not this pass");
    require(binding.at("system")==system && binding.at("coordinate")== (system=="pp"?"photon_et_gev":"centrality_percent"),"efficiency coordinate/system differs");
    require(binding.at("score_branch")==s.at("score_branch") && binding.at("nominal_line")== (kind=="tight_efficiency"?s.at("tight"):s.at("isolation").at("threshold")),"nominal efficiency line/score differs");
    const double nominal=binding.at("nominal_efficiency"),target=binding.at("target_efficiency");
    require(std::isfinite(nominal)&&std::isfinite(target)&&nominal>0&&nominal<1&&target>0&&target<1&&target!=nominal,"invalid efficiencies");
    require(binding.at("line")==payload.at("line"),"target efficiency line differs");
    replay::line(payload.at("line"),0.);
    if(kind=="isolation_efficiency")s["isolation"]["threshold"]=payload.at("line");
    else{
      require(std::abs(std::abs(target-nominal)-.05)<1.e-12,"tight variation must be five efficiency points");
      s["tight"]=payload.at("line");
      if(s.contains("optional_model")){
        auto& m=s["optional_model"];m["working_point_choice"]["tight"]=s["tight"];
        for(auto& cell:m["cells"])cell["tight"]=s["tight"];
        if(m.contains("native_evaluation"))for(auto& cell:m["native_evaluation"]["cells"])cell["tight"]=s["tight"];
      }
    }
  }else if(kind=="sideband_gap"){
    if(s.at("isolation").at("sideband_gap_gev").is_null()){
      require(payload.at("nominal_gap_policy")=="check_each_retained_witness","stored gap needs a per-source witness check");
      const double nominal=payload.at("nominal_gap_gev");require(std::isfinite(nominal)&&nominal>=0.,"stored nominal gap is absent");
    }else require(payload.at("nominal_gap_gev")==s.at("isolation").at("sideband_gap_gev"),"gap proposal does not match actual nominal");
    const double gap=payload.at("gap_gev");require(std::isfinite(gap)&&gap>=0,"invalid gap");
    s["isolation"]["sideband_gap_gev"]=gap;
  }else if(kind=="non_tight_lower"){
    require(payload.at("nominal_line")==s.at("non_tight_low"),"lower-edge proposal does not match actual nominal");
    const double low=payload.at("lower_score");require(std::isfinite(low)&&low>=0&&low<=1,"invalid score edge");
    s["non_tight_low"]={{"intercept",low},{"slope",0.}};
    if(s.contains("optional_model"))s["optional_model"]["working_point_choice"]["non_tight_lower"]=low;
  }else throw std::runtime_error("native selection views: unsupported variation kind");
  // The frozen other-system cuts remain byte-identical except nominal-label checks.
  for(auto& systemSpec:result["systems"]){systemSpec["verify_nominal_id"]=false;systemSpec["verify_nominal_isolation"]=false;}
  replay::cuts(result);
  result["selection_sha256"]=digest(replay::identityBody(result));
  result["profile"]="selection:"+result.at("selection_sha256").template get<std::string>();
  replay::identity(result,digest);return result;
}

template<class Digest,class FileDigest> Json adopted(const Json& parent,const std::string& name,
 const std::string& prescriptionPath,Digest digest,FileDigest fileDigest){
 namespace replay=photonjet_selection_replay;replay::identity(parent,digest);
 // The YAML and its relative retained-line table travel with the package.
 // File hashes bind derivation bytes; nominal equality binds applicability.
 const Json source=photonjet_config::readYaml(prescriptionPath);
 require(source.at("schema")=="PhotonJetSelectionVariationsV1","unknown selection-variation schema");
 const auto keys=[](const Json& object,std::initializer_list<const char*> expected){
  require(object.is_object() && object.size()==expected.size(),"unknown or missing selection-variation setting");
  for(const auto* key:expected)require(object.contains(key),"missing selection-variation setting");
 };
 keys(source,{"schema","sources","prescriptions"});
 keys(source.at("sources"),{"ppg12","auau_pass_lines"});
 keys(source.at("sources").at("ppg12"),{"version","text_sha256","pdf_sha256"});
 for(const char* key:{"text_sha256","pdf_sha256"})require(replay::digestString(source.at("sources").at("ppg12").at(key)),"invalid reference identity");
 const auto& prescriptions=source.at("prescriptions");
 keys(prescriptions,{"pp_tight_intercept_up","pp_non_tight_lower_down","pp_nonisolated_gap","auau_non_tight_lower","auau_non_tight_upper"});
 keys(prescriptions.at("pp_tight_intercept_up"),{"value"});
 keys(prescriptions.at("pp_non_tight_lower_down"),{"value"});
 keys(prescriptions.at("pp_nonisolated_gap"),{"nominal","values"});
 keys(prescriptions.at("auau_non_tight_lower"),{"nominal","value"});
 keys(prescriptions.at("auau_non_tight_upper"),{"nominal_efficiency","value_efficiency"});
 const auto finite=[](const Json& value){return value.is_number() && std::isfinite(value.template get<double>());};
 require(finite(prescriptions.at("pp_tight_intercept_up").at("value")) &&
  finite(prescriptions.at("pp_non_tight_lower_down").at("value")),"invalid pp score shift");
 const auto& gaps=prescriptions.at("pp_nonisolated_gap");
 require(gaps.at("values").is_array() && gaps.at("values").size()==2 && finite(gaps.at("nominal")) && gaps.at("nominal").template get<double>()>=0.,"invalid isolation gap definition");
 for(const auto& gap:gaps.at("values"))require(finite(gap) && gap.template get<double>()>=0.,"invalid isolation gap");
 const auto& lower=prescriptions.at("auau_non_tight_lower");
 require(finite(lower.at("nominal")) && finite(lower.at("value")) &&
  lower.at("nominal").template get<double>()>=0. && lower.at("value").template get<double>()>=0. &&
  lower.at("nominal").template get<double>()<=1. && lower.at("value").template get<double>()<=1.,"invalid lower score definition");
 const auto& upper=prescriptions.at("auau_non_tight_upper");
 require(upper.at("nominal_efficiency")==.8 && upper.at("value_efficiency")==.9,"upper efficiency does not identify the retained80/90 lines");
 Json result=parent;const bool pp=name.rfind("pp_",0)==0;auto& s=result["systems"][pp?"pp":"auau"];
 Json evidence=Json::array({{{"path",prescriptionPath},{"sha256",fileDigest(prescriptionPath)},{"locator","/prescriptions"}}});
 Json pass;
 if(!pp){
  const auto& pin=source.at("sources").at("auau_pass_lines");keys(pin,{"path","sha256"});
  require(replay::digestString(pin.at("sha256")),"invalid retained-line source hash");
  const auto relative=std::filesystem::path(pin.at("path").template get<std::string>());
  require(!relative.empty() && relative.is_relative() && std::none_of(relative.begin(),relative.end(),[](const auto& p){return p=="..";}),"retained-line source must be package-relative");
  const std::string path=(std::filesystem::path(prescriptionPath).parent_path()/relative).lexically_normal().string();
  require(fileDigest(path)==pin.at("sha256").template get<std::string>(),"AuAu pass line source bytes differ");
  std::ifstream file(path);file>>pass;require(pass.at("schema")=="PhotonJetWorkingPointLinesV1" && pass.at("coordinate")=="centrality_percent","unknown retained-line schema or coordinate");
  for(const char* key:{"qa_source_sha256","qa_cpp_sha256","native_helpers_sha256","existing_derivation_sha256","retained_source_sha256"})
   require(replay::digestString(pass.at(key)),"invalid retained derivation identity");
  for(const char* key:{"bdt60","bdt70","bdt80","bdt90","isolation70","isolation80","isolation90"}){
   const auto& l=pass.at("lines").at(key);
   require(l.size()==2 && finite(l.at("intercept")) && finite(l.at("slope_per_percent")),"invalid retained working-point line");
  }
  evidence.push_back({{"path",path},{"sha256",pin.at("sha256")},{"locator","/lines"},{"qa_source_sha256",pass.at("qa_source_sha256")}});
 }
 const auto line=[&](const char* key){const auto& l=pass.at("lines").at(key);return Json{{"intercept",l.at("intercept")},{"slope",l.at("slope_per_percent")}};};
 if(!pp)require(s.at("tight")==line("bdt70") && s.at("non_tight_high")==line("bdt80") &&
  s.at("isolation").at("threshold")==line("isolation80"),"AuAu nominal differs from retained70/80 working points");
 if(name=="pp_tight_plus005")s["tight"]["intercept"]=s.at("tight").at("intercept").template get<double>()+source.at("prescriptions").at("pp_tight_intercept_up").at("value").template get<double>();
 else if(name=="pp_lower_minus010")s["non_tight_low"]["intercept"]=s.at("non_tight_low").at("intercept").template get<double>()+source.at("prescriptions").at("pp_non_tight_lower_down").at("value").template get<double>();
 else if(name=="pp_gap01" || name=="pp_gap10")s["isolation"]["sideband_gap_gev"]=source.at("prescriptions").at("pp_nonisolated_gap").at("values").at(name=="pp_gap01"?0:1);
 else if(name=="auau_iso70" || name=="auau_iso90"){
  require(s.at("isolation").at("threshold")==line("isolation80"),"AuAu nominal isolation is not the source-bound80 line");
  s["isolation"]["threshold"]=line(name=="auau_iso70"?"isolation70":"isolation90");
 }else if(name=="auau_tight60" || name=="auau_tight80"){
  require(s.at("tight")==line("bdt70") && s.at("non_tight_high")==line("bdt80"),"AuAu nominal ID differs from source-bound70/80 lines");
  s["tight"]=line(name=="auau_tight60"?"bdt60":"bdt80");
 }else if(name=="auau_lower02"){
  require(s.at("non_tight_low")==Json({{"intercept",source.at("prescriptions").at("auau_non_tight_lower").at("nominal")},{"slope",0.}}),"AuAu nominal lower edge differs");
  s["non_tight_low"]={{"intercept",source.at("prescriptions").at("auau_non_tight_lower").at("value")},{"slope",0.}};
 }else if(name=="auau_upper90"){
  require(s.at("non_tight_high")==line("bdt80"),"AuAu nominal upper edge differs");s["non_tight_high"]=line("bdt90");
 }else throw std::runtime_error("native selection views: unsupported adopted prescription "+name);
 if(s.contains("optional_model")){
  auto& m=s["optional_model"];m["working_point_choice"]={{"mode","explicit"},{"variable",pp?"pho_et":"cent_value"},{"tight",s.at("tight")},{"non_tight_upper",s.at("non_tight_high")},{"non_tight_lower",s.at("non_tight_low").at("intercept")}};
  for(auto& c:m["cells"]){c["tight"]=s.at("tight");c["non_tight_high"]=s.at("non_tight_high");}
  if(m.contains("native_evaluation"))for(auto& c:m["native_evaluation"]["cells"]){c["tight"]=s.at("tight");c["non_tight_high"]=s.at("non_tight_high");}
 }
 for(auto& sys:result["systems"]){sys["verify_nominal_id"]=false;sys["verify_nominal_isolation"]=false;}
 replay::cuts(result);result["selection_sha256"]=digest(replay::identityBody(result));result["profile"]="selection:"+result.at("selection_sha256").template get<std::string>();replay::identity(result,digest);
 Json contract={{"schema","PhotonJetNativeSelectionPlanV2"},{"prescription",name},{"parent_selection_sha256",parent.at("selection_sha256")},{"effective_selection_sha256",result.at("selection_sha256")},{"selection",result},{"definition_sources",evidence},{"admission","SOURCE_BOUND_CANDIDATE_NOT_FINAL_SCIENCE_ADOPTION"},{"correction_policy","frozen_parent_weights_first_pass_candidate"}};
 if(name=="pp_gap01"||name=="pp_gap10")contract["expected_parent_stored_gap_gev"]=source.at("prescriptions").at("pp_nonisolated_gap").at("nominal");
 return contract;
}
} }
