#pragma once
#include "TreeToHists.h"
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
namespace treetohists {
inline std::vector<std::string> candidateCacheNames(){return {"events","photons","isolation","jets","truth_photons","truth_jets","photon_links","jet_links","unsub_jets"};}
struct CandidateCachePackage {
 std::map<std::string,RecordTable> rows;Long64_t offset=0;
 void append(TFile&in,const std::string&prefix,Long64_t source,TDirectory*dest,const std::vector<Long64_t>*map=nullptr){
  if(!dest)throw std::runtime_error("candidate record destination missing");
  auto* target=dest->GetDirectory("candidates");if(!target)target=dest->mkdir("candidates");
  auto*events=dynamic_cast<TTree*>(in.Get((prefix+"events").c_str()));if(!events)throw std::runtime_error("candidate cache event records missing");
  for(const auto&name:candidateCacheNames()){
   auto*t=dynamic_cast<TTree*>(in.Get((prefix+name).c_str()));if(!t)throw std::runtime_error("candidate cache table missing: "+name);
   rows[name].append(*t,name,source,offset,target,map);
  }offset+=events->GetEntries();
 }
 Json summary(double etFloor,double radius) const {
  Json result={{"schema","PhotonJetRecordsV1"},{"jet_radius",radius},{"reco_et_min_gev",etFloor},
   {"photon_index","dense_retained_per_event"},{"input_index","original_per_event"},
   {"photon_link_reco","remapped_dense_or_minus_one"},{"isolation_order","two_reco_witnesses_R03_R04_per_retained_photon"},
   {"truth_and_jets","unfiltered_selected_radius;unsubtracted_match_inputs_verified_then_omitted_for_SIM"},
   {"events",offset},{"rows",Json::object()},
   {"response_match_recipe","native_all_view_validation_then_selected_radius_analysis_view_greedy_v1"},
   {"response_match_witness_scope","SIM_verified_matches_1;DATA_witness_0"},
   {"omitted_after_SIM_match_validation",{"photon_links","jet_links","unsub_jets"}},
   {"raw_isolation_use","archival_only_forbidden_downstream"}};
  for(const auto&[n,t]:rows)result["rows"][n]=t.output?t.output->GetEntries():0;return result;
 }
 void write(){for(auto&[n,t]:rows)t.write();}
};
inline Json candidateRecordMetadata(const CandidateCachePackage& records,const Json& plan,const std::string& system){
 const auto& h=plan.at("histograms");
 Json result=records.summary(h.value("candidate_records_min_et",0.0),h.at("jet_radius").get<double>());
 result["jet_radius"]=h.at("jet_radius");
 const auto& selection=plan.at("selection").at("systems").at(system=="PP"?"pp":"auau");
 if(selection.contains("optional_model") && selection.at("optional_model").contains("native_evaluation")){
  const auto& contract=selection.at("optional_model").at("native_evaluation").at("feature_contract");
  result["native_model_input_order"]=contract.at("order");
  result["native_model_input_view"]=contract.at("view");
  result["native_model_input_semantic_sha256"]=contract.at("view_semantic_sha256");
 }
 return result;
}
inline void validateCandidateCache(TFile&f,Long64_t sources){
 auto*events=dynamic_cast<TTree*>(f.Get("records/candidates/events"));if(!events)throw std::runtime_error("missing candidate cache events");
 std::vector<std::pair<Long64_t,Long64_t>> ranges(sources,{0,0});Long64_t previousSource=-1;
 for(const auto&name:candidateCacheNames()){
  auto*t=dynamic_cast<TTree*>(f.Get(("records/candidates/"+name).c_str()));if(!t||!t->GetLeaf("source"))throw std::runtime_error("candidate cache source schema differs");
  t->SetBranchStatus("*",false);t->SetBranchStatus("source",true);
  if(name!="events"){if(!t->GetLeaf("event")||!t->GetLeaf("index"))throw std::runtime_error("candidate cache event schema differs");t->SetBranchStatus("event",true);t->SetBranchStatus("index",true);}
  Long64_t previousEvent=-1,index=-1;
  for(Long64_t i=0;i<t->GetEntries();++i){
   if(i%10000==0)ReductionProgress::pulse("validating_candidate_cache");
   if(t->GetEntry(i)<=0)throw std::runtime_error("unreadable candidate cache row");
   const auto source=t->GetLeaf("source")->GetValueLong64();if(source<0||source>=sources)throw std::runtime_error("candidate cache source outside inventory");
   if(name=="events"){
    if(source<previousSource)throw std::runtime_error("candidate cache source order differs");
    if(source!=previousSource)ranges[source].first=i;ranges[source].second=i+1;previousSource=source;
   }else{
    const auto event=t->GetLeaf("event")->GetValueLong64();const auto j=t->GetLeaf("index")->GetValueLong64();
    if(event<ranges[source].first||event>=ranges[source].second||event<previousEvent)throw std::runtime_error("candidate cache foreign event reference");
    if(event!=previousEvent)index=-1;
    if(j!=index+1)throw std::runtime_error("candidate cache candidate identity gap");index=j;previousEvent=event;
   }
  }t->SetBranchStatus("*",true);
 }
}
}
