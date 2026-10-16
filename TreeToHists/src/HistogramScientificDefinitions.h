/* Effective QA histogram definitions, emitted by the fill that used them.
   These are not run/configuration fingerprints. Stage-dependent selections
   remain conditional terms so a stage-0 projection does not depend on an ID
   cut applied only at later stages. Unknown branch ancestry remains explicit. */
#ifndef PHOTONJET_HISTOGRAM_SCIENTIFIC_DEFINITIONS_H
#define PHOTONJET_HISTOGRAM_SCIENTIFIC_DEFINITIONS_H
#include "ScientificTrace.h"
#include <regex>
#include <TH1.h>
#include <TNamed.h>
#include <TList.h>

namespace photonjet_science {
inline Json child(const Json& root,const std::string& pointer) {
  const Json::json_pointer p(pointer);
  return root.contains(p)?root.at(p):unknown("producer did not record "+pointer);
}
inline bool starts(const std::string& value,const std::string& prefix) { return value.rfind(prefix,0)==0; }

/* context is supplied from the validated Profile and native source metadata
   at fill time. It contains no output filename, run id, source digest or
   general producer-code digest. Unrecognized names cannot become complete. */
inline Json qaDefinition(const std::string& histogramName,const Json& context) {
  std::string stem=histogramName, population="all", centrality="all";
  std::smatch match;
  if(std::regex_match(stem,match,std::regex("(.*)_(cent[0-9]+_[0-9]+|all)$"))) {
    stem=match[1];centrality=match[2];
  }
  for(const char* label:{"jetprompt","phoother","phosim","jetsim","data","sig","bkg"}) {
    const std::string suffix=std::string("_")+label;
    if(stem.size()>=suffix.size() && stem.compare(stem.size()-suffix.size(),suffix.size(),suffix)==0) {
      population=label;stem.resize(stem.size()-suffix.size());break;
    }
  }
  static const std::set<std::string> exact={"ncb","iso","isonoue","abcd","recoil","scoreiso","scoreiso_fine",
    "clu_et","clu_etaphi","preselection_compare","abcdown","recoilown","evt_vz","evt_vz_novertexweight",
    "jet_pt","jet_etaphi","e11e33_dt","etatime","clu_time_stage","clu_etatime_stage","clu_etatime_npb",
    "clu_rawtime_stage","truth_iso","eff","bdt","isoslice","isonoueslice","evt_calo","evt_cent","evt_vzw",
    "cent_bits","evt_vz_selected","evt_vz_selected_wide","gen","gen_raw","reco_tag_vz30","trig","trig_study"};
  const bool shape=starts(stem,"ss_"),width=starts(stem,"width_");
  const bool recognized=exact.count(stem)||shape||width;
  Json declarations=Json::object();
  for(const auto& category:categories())declarations[category]=recognized?known("native per-observable term roster"):
    unknown("histogram family has no effective-definition rule: "+stem);
  Json terms=Json::array();
  const auto add=[&](const std::string& category,const std::string& id,const Json& value,const Json& ranges=Json::array()) {
    const auto refs=context.value("definition_references",Json::object());
    const std::string reference=refs.value(id,std::string());
    const auto origins=context.value("definition_origins",Json::object());
    terms.push_back(term(category,id,value,ranges,reference,reference.empty()?Json():origins.value(reference,Json())));
  };
  const auto read=[&](const std::string& category,const std::string& id,const std::string& pointer,const Json& ranges=Json::array()) {
    add(category,id,child(context,pointer),ranges);
  };
  const bool generator=starts(stem,"gen"),raw=stem=="gen_raw";
  const bool event=starts(stem,"evt_")||starts(stem,"cent_bits")||starts(stem,"trig");
  const bool jet=starts(stem,"jet_")||stem=="recoil"||stem=="recoilown";
  const bool truth=stem=="truth_iso"||stem=="eff"||generator;
  const bool slice=starts(stem,"width_")||stem=="bdt"||stem=="isoslice"||stem=="isonoueslice";
  const bool staged=shape||stem=="clu_et"||stem=="clu_etaphi"||stem=="clu_time_stage"||stem=="clu_etatime_stage"||stem=="clu_rawtime_stage";
  const int stageAxis=stem=="clu_et"?1:2;
  const bool region=stem=="abcd"||stem=="abcdown"||stem=="recoil"||stem=="recoilown"||stem=="reco_tag_vz30";
  const bool tight=starts(stem,"width_tight_iso_")||region;
  add("projection","fill_algorithm",{{"id","photon_qa_accumulation_v1"},{"accumulator","content=sum(w); variance=sum(w*w)"},
      {"cumulative_stage_axis",staged},{"finite_observable_required",true}});
  read("selection","collision","/collision");
  if(raw) {
    read("selection","generator_record_acceptance","/raw_generator_selection");
    read("selection","ownership","/ownership",window(1,3,4));
  } else if(stem=="trig_study")read("selection","trigger_study_selection","/trigger_study_selection");
  else {
    read("selection","event_record_acceptance","/event_record_selection");
    if(stem!="trig")read("selection","ownership","/ownership");
    if(!generator && stem!="trig")read("selection","valid_vertex","/valid_vertex");
    if(!generator && stem!="trig" && stem!="evt_vz" && stem!="evt_vz_novertexweight" && stem!="evt_calo" && stem!="evt_vz_selected_wide")
      read("selection","vertex_window","/vertex_window");
    if(!event && !generator || stem=="evt_vz_selected" || stem=="evt_vz_selected_wide")
      read("selection","data_event_selection","/data_event_selection");
  }
  if(!event && !generator && !truth && !jet)read("selection","candidate_geometry","/candidate_geometry");
  if(!event && !generator && !truth && !jet && stem!="clu_et" && stem!="clu_etaphi" && stem!="abcd" && stem!="abcdown") {
    if(stem=="etatime" || stem=="clu_etatime_stage" || stem=="clu_etatime_npb")
      read("selection","ncb_eta_time_energy_window","/ncb_eta_time_energy_window");
    else read("selection","qa_energy_window","/qa_energy_window");
  }
  if(stem=="clu_etaphi")add("selection","cluster_map_energy_min",{{"et_min_gev",12.0}});
  if(jet)read("selection","jet_selection","/jet_selection");
  if(truth||population!="data")read("selection","truth_classification","/truth_classification");
  Json centralityChannel=-1;
  if(centrality!="all") {
    const auto declared=child(context,"/centrality_channels/"+centrality);
    if(declared.contains("slot") && declared.contains("low") && declared.contains("high_exclusive")) {
      centralityChannel=declared.at("slot");
      add("selection","centrality_class",{{"low",declared.at("low")},{"high_exclusive",declared.at("high_exclusive")}});
    } else {centralityChannel="unmapped";add("selection","centrality_class",unknown("native centrality slot was not declared"));}
  }
  if(staged) {
    read("selection","preselection","/preselection",window(stageAxis,1,4));
    read("models","ncb_model","/ncb_model",window(stageAxis,1,4));
    read("selection","tight_id","/tight_id",window(stageAxis,2,4));
    read("models","selected_score_model","/score_model",window(stageAxis,2,4));
    read("selection","isolation","/isolation",window(stageAxis,3,4));
  } else if(stem=="eff") {
    read("selection","matched_truth","/matched_truth",window(1,1,6));
    read("selection","preselection","/preselection",window(1,2,6));
    read("models","ncb_model","/ncb_model",window(1,2,6));
    read("selection","tight_id","/tight_id",window(1,4,6));
    read("models","selected_score_model","/score_model",window(1,4,6));
    read("selection","isolation","/isolation",window(1,3,4));
    read("selection","isolation_stage5","/isolation",window(1,5,6));
  } else if(starts(stem,"jet_")) {
    const auto selected=window(stem=="jet_pt"?1:2,1,2);
    read("selection","region_leader","/region_leader",selected);
    read("selection","preselection","/preselection",selected);
    read("models","ncb_model","/ncb_model",selected);
    read("selection","tight_id","/tight_id",selected);
    read("models","selected_score_model","/score_model",selected);
    read("selection","isolation","/isolation",selected);
  } else if(!event && !generator && stem!="ncb" && stem!="etatime" && stem!="e11e33_dt" && stem!="truth_iso") {
    read("selection","preselection","/preselection");
    read("models","ncb_model","/ncb_model");
    if(tight || stem=="iso" || stem=="isonoue" || stem=="eff") {
      read("selection","tight_id","/tight_id");read("models","selected_score_model","/score_model");
    }
    if(tight)read("selection","isolation","/isolation");
  }
  if(region) {
    read("selection","region_leader","/region_leader");
    read("selection","non_tight_id","/non_tight_id");
    read("selection","isolation_sideband","/isolation_sideband");
  }
  if(stem=="ncb")read("models","ncb_model","/ncb_model");
  if(stem=="bdt" || stem=="ss_bdt" || starts(stem,"scoreiso"))read("models","selected_score_model_observable","/score_model");
  if(shape && stem!="ss_bdt")read("projection","shower_observable","/branch_definitions/"+stem.substr(3));
  if(width)read("projection","shower_observable","/branch_definitions/"+stem.substr(starts(stem,"width_tight_iso_")?16:6));
  if(stem.find("time")!=std::string::npos || stem=="e11e33_dt")read("projection","cluster_time_reference","/cluster_time_reference");
  if(starts(stem,"iso")||stem=="truth_iso")read("projection","isolation_observable","/isolation_observable");
  const bool unitWeight=raw||starts(stem,"trig")||stem=="cent_bits"||stem=="evt_vz_selected"||stem=="evt_vz_selected_wide";
  if(unitWeight)add("weights","event_weight",{{"rule","unit"}});
  else {
    read("weights","base_event_weight","/base_event_weight");
    if(stem!="evt_vz_novertexweight" && stem!="evt_cent")
      read("weights","vertex_weight","/vertex_weight",stem=="evt_vzw"?window(1,1,2):Json::array());
    if(stem=="evt_cent")read("weights","vertex_weight","/vertex_weight");
    if(!generator && context.value("simulation",false) && context.value("system","")=="AuAu") {
      const auto weights=child(context,"/centrality_weights");
      if(weights.is_array())for(std::size_t i=0;i<weights.size();++i) {
        // Slice coordinates in a slices histogram are 0..16, not percent.
        if(slice)add("weights","centrality_slice_"+std::to_string(i),weights[i],window(2,i,i+1));
        else if(centrality!="all" && std::regex_match(centrality,match,std::regex("cent([0-9]+)_([0-9]+)"))) {
          const int lo=std::stoi(match[1]),hi=std::stoi(match[2]);
          if(5*int(i)>=lo && 5*int(i)<hi)add("weights","centrality_slice_"+std::to_string(i),weights[i]);
        } else add("weights","centrality_slice_"+std::to_string(i),weights[i],stem=="evt_cent"?window(1,1,2):Json::array());
      } else add("weights","centrality_weights",weights);
    }
  }
  add("normalization","sample_normalization",{{"applied",false},{"consumer","FinalAnalysis applies the exact saved sample factor once"}});
  return definition(stem,{{"system",context.value("system","")},{"population",population},{"centrality_slot",centralityChannel},{"source_channel",context.value("source_channel",Json::object())}},terms,declarations);
}
inline void putDefinition(TH1& histogram,const Json& definition) {
  auto* functions=histogram.GetListOfFunctions();
  if(auto* old=functions->FindObject("physics_effective_definition"))functions->Remove(old);
  functions->Add(new TNamed("physics_effective_definition",definition.dump().c_str()));
  if(auto* old=functions->FindObject("physics_covariance"))functions->Remove(old);
  functions->Add(new TNamed("physics_covariance",covariance("not_computed",{},0,
      "native photon QA accumulation stores Sumw2 diagonal variances; no cross-bin covariance calculated").dump().c_str()));
}
inline Json storedDefinitions(const TH1& histogram) {
  auto* functions=histogram.GetListOfFunctions();
  if(auto* many=dynamic_cast<TNamed*>(functions->FindObject("physics_effective_definitions")))return Json::parse(many->GetTitle());
  if(auto* one=dynamic_cast<TNamed*>(functions->FindObject("physics_effective_definition")))return Json::array({Json::parse(one->GetTitle())});
  return Json::array({Json()}); // one contributor has no recorded ancestry
}
inline void mergeDefinitions(TH1& destination,const TH1& input) {
  Json rows=storedDefinitions(destination);
  for(const auto& row:storedDefinitions(input))if(std::find(rows.begin(),rows.end(),row)==rows.end())rows.push_back(row);
  std::sort(rows.begin(),rows.end(),[](const Json& a,const Json& b){return a.dump()<b.dump();});
  auto* functions=destination.GetListOfFunctions();
  for(const char* name:{"physics_effective_definition","physics_effective_definitions"})
    if(auto* old=functions->FindObject(name))functions->Remove(old);
  functions->Add(new TNamed("physics_effective_definitions",rows.dump().c_str()));
}
}  // namespace photonjet_science
#endif
