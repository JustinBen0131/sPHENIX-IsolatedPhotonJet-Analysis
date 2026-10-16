/* Typed definitions for measurement and response accumulators. These terms
   are copied from the validated effective Profile and from observed upstream
   metadata, not from a run/source digest. Old upstream metadata stays unknown. */
#pragma once
#include "QAScientificContext.h"

namespace photonjet_science {
inline Json measurementProfileWitness(const photonjet::Profile& p) {
  Json cells=Json::array();for(const auto& c:p.idCells)cells.push_back({{"et",{c.etMin,c.etMax}},{"centrality",{c.centMin,c.centMax}},
    {"tight",{c.tight.intercept,c.tight.slope}},{"non_tight_high",{c.nonTightHigh.intercept,c.nonTightHigh.slope}}});
  Json native=absent("native model evaluation disabled");
  if(p.nativeID) {
    const auto selection=Json::parse(p.selectionJson);
    const auto& evaluation=selection.at("systems").at(p.system==photonjet::System::pp?"pp":"auau").at("optional_model").at("native_evaluation");
    native={{"feature_contract",evaluation.at("feature_contract")},{"domain",evaluation.at("domain")},{"cells",Json::array()}};
    for(const auto& c:evaluation.at("cells"))native["cells"].push_back({{"et",c.at("et")},{"centrality",c.at("centrality")},
      {"function",c.value("scientific_function",unknown("native model function absent"))}});
  }
  Json preselection=photonjet::canonicalPreselection();preselection.erase("npb_model_sha256");
  return portable({{"system",p.system==photonjet::System::pp?"PP":"AuAu"},{"simulation",p.simulation},
    {"native_model",native},{"preselection",preselection},
    {"vertex",p.vertexAbsZMax},{"photon",{p.photonEtMin,p.photonEtMax,p.photonAbsEtaMax}},
    {"jets",{p.jetRadius,p.jetAbsEtaMax,p.jetPtMin,p.overlapDeltaR,p.deltaPhiMin}},
    {"id",{{"score_branch",p.scoreBranch},{"state_branch",p.stateBranch},{"native",p.nativeID},{"optional",p.optionalID},
      {"augmented",p.augmentedClasses},{"tight",{p.tight.intercept,p.tight.slope}},
      {"non_tight_low",{p.nonTightLow.intercept,p.nonTightLow.slope}},
      {"non_tight_high",{p.nonTightHigh.intercept,p.nonTightHigh.slope}},
      {"constant_low",p.nonTightLowIsConstant},{"low_constant",p.nonTightLowConstant},{"cells",cells}}},
    {"isolation",{{"radius_index",p.isolationRadiusIndex},{"subtraction",p.isolationSubtraction},
      {"from_witness",p.isolationFromWitness},{"threshold",{p.isolationThreshold.intercept,p.isolationThreshold.slope}},
      {"centrality",{p.isolationCentMin,p.isolationCentMax}},
      {"sideband_gap",std::isfinite(p.sidebandGap)?Json(p.sidebandGap):absent("stored witness sideband")}}}});
}
inline Json measurementContext(TFile& file,const photonjet::Profile& profile,
    const photonjet::CentralityWeights& weights,const Json& variation=Json::object(),
    const Json& axes=Json::object(),double jetScale=1.,bool timing=false,const std::string& trigger="",
    const Json& loadedFunctions=Json()) {
  const photonjet::DataEventSelection selected(file,profile);
  const bool centralityOff=variation.value("centrality_weight",std::string())=="off";
  Json c=qaContext(file,profile,weights,selected,0.,centralityOff,loadedFunctions);
  c["profile_witness"]=measurementProfileWitness(profile);
  const auto upstream=upstreamDefinitions(file);
  c["measurement_photon_selection"]={{"et_min",profile.photonEtMin},{"et_max_exclusive",profile.photonEtMax},
    {"abs_eta_max",profile.photonAbsEtaMax},{"geometry",child(c,"/candidate_geometry")}};
  c["truth_photon_selection"]={{"pt_min_exclusive",0.},{"abs_eta_max",photonjet::kTruthPhotonAbsEtaMax},
    {"isolation_max",photonjet::kTruthIsolationMaxGeV},{"source_role",photonjet::kSignalSourceRole},
    {"leader","highest_pt; exact ties select later index"},{"upstream",child(upstream,"/truth_classification")}};
  c["response_matching"]=child(upstream,"/response_matching");
  c["energy_variation"]=variation;
  c["measurement_axes"]=axes;
  c["data_jet_scale"]=jetScale;
  c["timing_selection"]=timing?Json{{"rule","finite(time_ns) and abs(time_ns-center_ns)<half_width_ns"},{"center_ns",-2.0},{"half_width_ns",6.0},{"time_definition",child(upstream,"/cluster_timing")}}:
      absent("no DATA timing cut in this accumulation");
  c["explicit_trigger"]=trigger.empty()?absent("no extra trigger selector"):Json{{"name",trigger},{"scaled_mask",photonjet::dataTriggerMask(file,trigger)}};
  return c;
}

inline Json measurementDefinition(const std::string& name,const Json& context,bool response=false) {
  std::smatch match;std::string stem=name,centrality="all";double threshold=-1;
  const std::regex suffix("(.*)_(cent[0-9]+_[0-9]+|all)(?:_jet([0-9]+))?(?:_(near|far))?$");
  std::string side;
  if(std::regex_match(name,match,suffix)) {
    stem=match[1];centrality=match[2];if(match[3].matched)threshold=std::stod(match[3]);if(match[4].matched)side=match[4];
  }
  const bool event=stem=="events",kinematic=stem=="photons";
  const bool pureTruth=response && (stem=="photon_truth" || stem=="xj_truth");
  const bool jets=starts(stem,"xj_")||starts(stem,"recoil")||starts(stem,"dphi")||starts(stem,"leading_")||starts(stem,"reco_pair_");
  const bool dphi=starts(stem,"dphi");
  const std::set<std::string> knownResponse={"photon_reco","photon_truth","photon_fakes","photon_misses",
    "photon_response_truth_x_reco","signal_leakage_truth_matched_abcd","xj_reco","xj_truth",
    "xj_response_truth_global_x_reco_global","xj_fakes","xj_fake_photon_reco","xj_combinatoric_reco",
    "xj_detector_fakes_reco","xj_misses","xj_leading_reco","xj_leading_matched_reco","xj_leading_reco_u",
    "xj_leading_matched_reco_u","xj_leading_matched_any_reco_u","xj_leading_global_reco",
    "xj_leading_global_matched_reco","reco_pair_truth_reason","dphi_truthjet_near","dphi_truthjet_far",
    "dphi_combinatoric_near","dphi_combinatoric_far"};
  const bool recognized=response?knownResponse.count(stem):event||kinematic||stem=="abcd"||stem=="recoil"||stem=="leading_global"||stem=="dphiA"||stem=="dphiC";
  Json terms=Json::array(),declarations=Json::object();
  for(const auto& category:categories())declarations[category]=recognized?known("native measurement term roster"):
      unknown("measurement observable has no typed term roster: "+stem);
  const auto add=[&](const std::string& category,const std::string& id,const Json& value) {
    const auto refs=context.value("definition_references",Json::object());const std::string reference=refs.value(id,std::string());
    const auto origins=context.value("definition_origins",Json::object());
    terms.push_back(term(category,id,value,Json::array(),reference,reference.empty()?Json():origins.value(reference,Json())));
  };
  const auto read=[&](const std::string& category,const std::string& key){add(category,key,child(context,"/"+key));};
  for(const char* key:{"collision","event_record_selection","ownership","valid_vertex","vertex_window","data_event_selection","explicit_trigger"})read("selection",key);
  if(centrality!="all") {
    std::smatch range;
    if(std::regex_match(centrality,range,std::regex("cent([0-9]+)_([0-9]+)")))
      add("selection","centrality_class",{{"low",std::stoi(range[1])},{"high_exclusive",std::stoi(range[2])}});
  }
  if(!event && !pureTruth) {
    read("selection","measurement_photon_selection");
    if(!kinematic) {
      for(const char* key:{"preselection","tight_id","non_tight_id","isolation","isolation_sideband","region_leader","timing_selection"})read("selection",key);
      read("models","ncb_model");read("models","score_model");
    }
  }
  if(jets) {
    Json selection=child(context,"/jet_selection");
    // Delta-phi controls a recoil cut only for the recoil window, not for a
    // histogram whose observable is that very angle before the cut.
    if(dphi && selection.is_object())selection.erase("delta_phi_min");
    add("selection","jet_selection",selection);add("selection","jet_threshold_gev",threshold);
    if(!side.empty())add("selection","eta_side",{{"side",side},{"boundary_abs_deta",.8}});
    read("projection","data_jet_scale");
  }
  if(response) {
    if(pureTruth || stem!="photon_reco")read("selection","truth_photon_selection");
    if(!pureTruth && stem!="photon_reco" && stem!="xj_reco")read("selection","response_matching");
    Json axes=child(context,"/measurement_axes");if(!jets && axes.is_object())axes.erase("xj");
    add("projection","measurement_axes",axes);
  }
  if(!event && !pureTruth) {
    Json variation=child(context,"/energy_variation");
    if(!jets && variation.is_object())for(const char* key:{"jet_scale","jet_smear","jet_resolution"})variation.erase(key);
    // The random seed is irrelevant when every smearing operation is off.
    if(variation.is_object() && variation.value("jet_smear",0.)==0. && variation.value("photon_smear",0.)==0. &&
       variation.value("jet_resolution",std::string("off"))=="off" && variation.value("photon_resolution",std::string("off"))=="off")variation.erase("seed");
    add("projection","energy_variation",variation);
  }
  read("weights","base_event_weight");read("weights","vertex_weight");
  if(context.value("simulation",false) && context.value("system","")=="AuAu") {
    const auto weights=child(context,"/centrality_weights");std::smatch range;
    const bool bounded=std::regex_match(centrality,range,std::regex("cent([0-9]+)_([0-9]+)"));
    if(weights.is_array())for(std::size_t i=0;i<weights.size();++i) {
      if(!bounded || (5*int(i)>=std::stoi(range[1]) && 5*int(i)<std::stoi(range[2])))
        add("weights","centrality_slice_"+std::to_string(i),weights[i]);
    } else add("weights","centrality_weights",weights);
  }
  if(context.contains("replay_normalization"))add("normalization","replay_factor",context.at("replay_normalization"));
  add("normalization","sample_normalization",{{"applied",false},{"rule","FinalAnalysis applies sample cross-section factor once"}});
  add("projection","accumulation",{{"algorithm",response?"accumulateResponse/v1":"accumulateMeasurement/v1"},
    {"variance","sum_w_squared"},{"histogram",stem}});
  int thresholdSlot=-1;for(std::size_t i=0;i<photonjet::kJetThresholds.size();++i)if(photonjet::kJetThresholds[i]==threshold)thresholdSlot=i;
  Json centralitySlot=-1;
  if(centrality!="all") {
    const auto declared=child(context,"/centrality_channels/"+centrality);
    if(declared.contains("slot"))centralitySlot=declared.at("slot");
    else {centralitySlot="unmapped";add("selection","centrality_slot_mapping",unknown("native centrality slot was not declared"));}
  }
  return definition(stem,{{"system",context.value("system","")},{"centrality_slot",centralitySlot},{"threshold_slot",thresholdSlot},
    {"source_channel",context.value("source_channel",Json::object())}},terms,declarations);
}

inline void observeMeasurementDefinition(TH1& h,const Json& context,bool response=false) {
  auto* functions=h.GetListOfFunctions();
  Json definitions=Json::array();
  if(auto* old=dynamic_cast<TNamed*>(functions->FindObject("physics_effective_definitions"))) {
    definitions=Json::parse(old->GetTitle());functions->Remove(old);
  }
  const auto definition=measurementDefinition(h.GetName(),context,response);
  if(std::find(definitions.begin(),definitions.end(),definition)==definitions.end())definitions.push_back(definition);
  functions->Add(new TNamed("physics_effective_definitions",definitions.dump().c_str()));
  if(auto* old=functions->FindObject("physics_covariance"))functions->Remove(old);
  functions->Add(new TNamed("physics_covariance",covariance("not_computed",{},0,
      "native accumulation stores Sumw2 diagonal variances; event correlations retained in record tables").dump().c_str()));
}
}
