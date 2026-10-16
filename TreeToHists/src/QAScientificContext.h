/* Effective values observed by FillPhotonQA, with unavailable upstream branch
   ancestry kept explicit. This performs no event loop and changes no cuts. */
#ifndef PHOTONJET_QA_SCIENTIFIC_CONTEXT_H
#define PHOTONJET_QA_SCIENTIFIC_CONTEXT_H
#include "HistogramScientificDefinitions.h"
#include "PhotonJet.h"
namespace photonjet_science {
// The eta-time maps have no ET axis: apply their common comparison window
// during filling and record it separately from the ET-resolved QA window.
inline constexpr double kNcbEtaTimeEtMinGeV = 12.0, kNcbEtaTimeEtMaxGeV = 40.0;
inline bool ncbEtaTimeEnergyAccepted(double et) {
  return kNcbEtaTimeEtMinGeV <= et && et < kNcbEtaTimeEtMaxGeV;
}
inline Json metadata(TFile& file,const std::string& key) {
  if(!file.Get(key.c_str()))return Json::object();
  return Json::parse(photonjet::metaText(file,key.c_str()));
}
inline Json commonSourceDefinitions(const std::vector<Json>& inputs) {
  if(inputs.empty())return unknown("no contributing native source definition");
  if(std::all_of(inputs.begin(),inputs.end(),[&](const Json& x){return x==inputs.front();}))return inputs.front();
  if(std::all_of(inputs.begin(),inputs.end(),[](const Json& x){return x.is_object() && !x.contains("state");})) {
    std::set<std::string> keys;for(const auto& input:inputs)for(const auto& item:input.items())keys.insert(item.key());
    Json result=Json::object();
    for(const auto& key:keys) {
      std::vector<Json> values;for(const auto& input:inputs)values.push_back(input.value(key,unknown("definition absent in a contributing source")));
      result[key]=commonSourceDefinitions(values);
    }
    return result;
  }
  Json value=unknown("contributing sources have different or incomplete effective definitions");
  std::set<std::string> unique;for(const auto& input:inputs)unique.insert(input.dump());
  value["observed_definitions"]=Json::array();for(const auto& input:unique)value["observed_definitions"].push_back(Json::parse(input));
  return value;
}
inline Json upstreamDefinitions(TFile& file) {
  // Consolidation preserves original metadata per native entry interval. A
  // copied first-file header cannot certify the remaining intervals. Resolve
  // each field independently so an unrelated branch does not taint it.
  if(file.Get("meta/consolidation")) {
    const auto consolidation=metadata(file,"meta/consolidation");std::vector<Json> sources;
    for(const auto& chunk:consolidation.at("chunks"))if(chunk.at("entry_end").get<long long>()>chunk.at("entry_begin").get<long long>())
      sources.push_back(chunk.at("original_meta").value("scientific_definitions",unknown("contributing chunk has no scientific definitions")));
    return commonSourceDefinitions(sources);
  }
  return metadata(file,"meta/scientific_definitions");
}
inline Json lookup(const TH1* h) {
  if(!h)return absent("vertex weighting disabled in this native invocation");
  Json centers=Json::array(),values=Json::array();
  for(int i=1;i<=h->GetNbinsX();++i) {centers.push_back(h->GetBinCenter(i));values.push_back(h->GetBinContent(i));}
  return portable({{"lookup","TH1::Interpolate"},{"centers",centers},{"values",values},
      {"support",{h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax()}},{"outside",0.0}});
}
inline Json modelFeatureDefinitions(const Json& upstream,const Json& function,const std::string& view="") {
  const auto features=function.value("feature_order",Json::array());
  if(!features.is_array() || features.empty())return unknown("model feature order not recorded consistently");
  Json values=Json::object();
  for(const auto& feature:features) {
    const std::string name=feature.get<std::string>();
    if(name=="cluster_Et" || name=="cluster_Eta")values[name]=child(upstream,"/candidate_geometry/kinematics");
    else if(name=="vertexz")values[name]=child(upstream,"/vertex");
    else if(name=="centrality")values[name]=child(upstream,"/event_acceptance/centrality_measurement");
    else if(!view.empty())values[name]=child(upstream,"/shower_feature_views/"+view);
    else values[name]=child(upstream,"/branch_definitions/"+(starts(name,"cluster_")?name.substr(8):name));
  }
  return values;
}
inline Json qaContext(TFile& file,const photonjet::Profile& p,const photonjet::CentralityWeights& weights,
                      const photonjet::DataEventSelection& selected,double mbdReference,bool centralityOff,
                      const Json& loadedFunctions=Json()) {
  using namespace photonjet;
  const bool sim=isSimulation(file),auau=p.system==System::auau;
  const Json upstream=upstreamDefinitions(file),source=metadata(file,"meta/source");
  Json out={{"system",auau?"AuAu":"PP"},{"simulation",sim},{"collision",{{"system",auau?"AuAu":"PP"}}},
    {"source_channel",{{"sample",source.value("sample","")},{"period",source.value("period","")},{"si_di_role",source.value("si_di_role","")}}},
    {"event_record_selection",{{"evt_record",0},{"evt_nodes_ready",1},{"terminal_rule","acceptedEventClass"},
       {"upstream_acceptance",child(upstream,"/event_acceptance")}}},
    {"raw_generator_selection",{{"evt_record",0},{"auau_terminal",{0,4}},{"auau_centrality",{0.,80.}},
       {"upstream_acceptance",child(upstream,"/event_acceptance")}}},
    {"trigger_study_selection",{{"finite_vz",true},{"abs_vz_lt_cm",30.},{"max_cluster_e_ge_gev",1.},
       {"upstream_acceptance",child(upstream,"/event_acceptance")}}},
    {"valid_vertex",{{"vtx_valid",1}}},{"vertex_window",{{"abs_vz_lt_cm",p.vertexAbsZMax},{"finite",true}}},
    {"data_event_selection",sim?absent("simulation does not apply the DATA trigger selector"):
       Json{{"active",selected.active},{"scaled_mask",selected.mask},{"minimum_bias",selected.minimumBias},{"run_windows",selected.runWindows}}},
    {"candidate_geometry",{{"finite_et_eta_phi",true},{"abs_eta_lt",p.photonAbsEtaMax},
       {"upstream_candidates",child(upstream,"/candidate_geometry")}}},
    {"qa_energy_window",{{"et_min_gev",10.},{"et_max_exclusive_gev",40.}}},
    {"ncb_eta_time_energy_window",{{"et_min_gev",kNcbEtaTimeEtMinGeV},{"et_max_exclusive_gev",kNcbEtaTimeEtMaxGeV}}},
    {"jet_selection",{{"radius",p.jetRadius},{"finite_pt_eta_phi",true},{"abs_eta_lt",p.jetAbsEtaMax},
       {"pt_min_gev",p.jetPtMin},{"overlap_delta_r",p.overlapDeltaR},{"delta_phi_min",p.deltaPhiMin},
       {"upstream_jets",child(upstream,"/jets")}}},
    {"truth_classification",child(upstream,"/truth_classification")},
    {"matched_truth",{{"plink_class",kLinkMatch},{"truth_signal_role",kSignalSourceRole}}},
    {"region_leader",{{"photon_et",{p.photonEtMin,p.photonEtMax}},{"selection","highest_ET_in_each_ABCD_region"}}},
    {"ownership",{{"owner",static_cast<int>(configuredOwnership.owner)},{"low",configuredOwnership.low},{"high",std::isfinite(configuredOwnership.high)?Json(configuredOwnership.high):Json("unbounded")},
                  {"upstream",child(upstream,"/ownership")}}},
    {"branch_definitions",child(upstream,"/branch_definitions")},{"ncb_model",child(upstream,"/models/npb")},
    {"cluster_time_reference",{{"mbd_reference_ns",mbdReference},{"upstream",child(upstream,"/cluster_timing")}}}};
  out["centrality_channels"]=Json::object();int centralitySlot=0;
  for(const auto& c:centralityClasses(p.system))out["centrality_channels"][c.name]={{"slot",centralitySlot++},{"low",c.low},{"high_exclusive",c.high}};
  Json pre=canonicalPreselection();pre.erase("npb_model_sha256");
  pre["primary_observables"]={{"et1",child(upstream,"/branch_definitions/et1")},
    {"e11_over_e33",child(upstream,"/branch_definitions/e11_over_e33")},
    {"e32_over_e35",child(upstream,"/branch_definitions/e32_over_e35")}};
  out["preselection"]=pre;
  out["ncb_model"]={{"function",child(upstream,"/models/npb")},
    {"features",modelFeatureDefinitions(upstream,child(upstream,"/models/npb"))}};
  Json cells=Json::array();for(const auto& c:p.idCells)cells.push_back({{"et",{c.etMin,c.etMax}},
      {"centrality",{c.centMin,c.centMax}},{"tight",{c.tight.intercept,c.tight.slope}},
      {"non_tight_high",{c.nonTightHigh.intercept,c.nonTightHigh.slope}}});
  out["tight_id"]={{"score_branch",p.scoreBranch},{"augmented_classes",p.augmentedClasses},{"optional",p.optionalID},
      {"state_branch",p.optionalID?"finite_score":p.stateBranch},{"line",{p.tight.intercept,p.tight.slope}},{"cells",cells},
      {"augmented_definition",p.augmentedClasses?child(upstream,"/nominal_id"):absent("configured score selector")}};
  out["non_tight_id"]={{"low",{p.nonTightLow.intercept,p.nonTightLow.slope}},{"high",{p.nonTightHigh.intercept,p.nonTightHigh.slope}},
      {"low_constant",p.nonTightLowIsConstant},{"low_value",p.nonTightLowConstant}};
  out["isolation"]={{"radius_index",p.isolationRadiusIndex},{"subtraction_method",p.isolationSubtraction},
      {"from_witness",p.isolationFromWitness},{"line",{p.isolationThreshold.intercept,p.isolationThreshold.slope}},
      {"centrality_domain",{p.isolationCentMin,p.isolationCentMax}},
      {"witness_definition",child(upstream,"/isolation")}};
  out["isolation_sideband"]={{"gap",std::isfinite(p.sidebandGap)?Json(p.sidebandGap):child(upstream,"/isolation/sideband_gap_gev")}};
  out["isolation_observable"]={{"radius_index",p.isolationRadiusIndex},{"subtraction_method",p.isolationSubtraction},
      {"definition",child(upstream,"/isolation")}};
  const auto selectedView=p.scoreBranch=="pho_bdt_score"?child(upstream,"/nominal_id/view"):Json("");
  out["score_model"]={{"function",child(upstream,"/models/"+p.scoreBranch)},
    {"features",selectedView.is_string()?modelFeatureDefinitions(upstream,child(upstream,"/models/"+p.scoreBranch),selectedView.get<std::string>()):
      unknown("stored nominal score feature view absent")}};
  if(p.nativeID) {
    const auto resolved=Json::parse(p.selectionJson);
    const auto& native=resolved.at("systems").at(auau?"auau":"pp").at("optional_model").at("native_evaluation");
    Json functions=Json::array();for(const auto& c:native.at("cells"))functions.push_back({{"et",c.at("et")},{"centrality",c.at("centrality")},
      {"function",c.value("scientific_function",unknown("native model numeric function not exported"))}});
    if(!loadedFunctions.is_null())functions=loadedFunctions.empty()?unknown("no loaded native model witnessed"):loadedFunctions;
    const auto& feature=native.at("feature_contract");
    out["score_model"]={{"routing","half_open_ET_centrality"},{"functions",functions},
      {"feature_order",feature.at("order")},{"input_precision","binary32"},
      {"feature_calculation",modelFeatureDefinitions(upstream,Json{{"feature_order",feature.at("order")}},feature.at("view").get<std::string>())}};
  }
  out["base_event_weight"]=sim?child(upstream,"/event_weight"):Json{{"rule","unit"}};
  if(sim && auau && (strictReduction||declaredAuAuCopies(metadata(file,"meta/weighting"))))out["base_event_weight"]={{"rule","unit_base_for_recomputed_AuAu_weights"}};
  out["vertex_weight"]=lookup(auauVertexFit.get());
  out["centrality_weights"]=Json::array();for(double w:weights.perSlice)out["centrality_weights"].push_back(centralityOff?1.:w);
  out["definition_references"]=upstream.value("definition_references",Json::object());
  out["definition_origins"]=upstream.value("definition_origins",Json::object());
  return out;
}
}  // namespace photonjet_science
#endif
