/* Optional candidate cache: finite reco photons above the recorded ET floor, plus complete jets/truth of accepted events.
   Missing optional fields have has_<name>=0 and NaN; local replay must gate their use. */
#pragma once
#include "ResponseMatching.h"
#include "CandidateScientificContext.h"
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TDirectory.h>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include <functional>
#include <cmath>
#include <numeric>
namespace photonjet { namespace candidate_cache {
inline bool enabled=false;
inline Long64_t maxRows=2000000;
inline double minRecoEt=0.;
inline int radiusIndex(double radius){
 if(radius==.2)return 0;if(radius==.3)return 1;if(radius==.4)return 2;
 throw std::runtime_error("photon-jet records require jet radius 0.2, 0.3 or 0.4");
}
inline std::string radiusSuffix(double radius){return std::array<std::string,3>{"02","03","04"}.at(radiusIndex(radius));}
inline std::string sourceField(const std::string& field,double radius){
 for(const std::string kind:{"jetunsub","truthjet","jet"}){
  const std::string prefix=kind+"_";
  if(field.rfind(prefix,0)==0)return kind+radiusSuffix(radius)+field.substr(kind.size());
 }return field;
}
struct WeightFactors {double stored=1,base=1,vertex=1,centrality=1;};
using PhotonValues = std::map<std::string,double>;
struct TableSpec {std::string name,count;std::vector<std::string> fields;};
inline const std::vector<TableSpec>& specs(){static const std::vector<TableSpec> s={
{"photons","pho_et",{"selected_score","pho_et","pho_eta","pho_phi","pho_bdt_score","pho_bdt_state","pho_presel_active","pho_presel_bits","pho_id_class","pho_iso_class","pho_ppg12_tower_mask_state","pho_ppg12_source_eligible","pho_npb_score","pho_npb_score_state","pho_npb_pass","pho_npb_pass_state","pho_npb_score_native","pho_auau_npb_score_native","pho_npbcap_raw_score","pho_npbcap_evaluation_state","pho_aunpbcap_raw_score","pho_aunpbcap_evaluation_state","pho_time_ns","pho_time_raw","pho_time_raw_state","pho_time_ntowers","pho_time_den_e","pho_ntime_num","pho_ntime_den","pho_ntime_ntowers","pho_ntime_finite","pho_bdt_score_cent_coarse_v1","pho_bdt_score_cent_fine_v1","pho_bdt_score_depth_10_v1","pho_bdt_score_depth_2_v1","pho_bdt_score_depth_3_v1","pho_bdt_score_depth_5_v1","pho_bdt_score_depth_6_v1","pho_bdt_score_depth_7_v1","pho_bdt_score_depth_8_v1","pho_bdt_score_depth_9_v1","pho_bdt_score_et12to40_v1","pho_bdt_score_et_cent_coarse_v1","pho_bdt_score_et_cent_fine_v1","pho_bdt_score_nominal_retrained_v1","pho_bdt_score_ppg12_retrained_v1","pho_presel_reference","pho_below15_state","pho_finite_feature_state","pho_et1","pho_e11_over_e33","pho_e32_over_e35","pho_weta_cogx","pho_wphi_cogx","pho_weta33_cogx","pho_wphi33_cogx","pho_weta35_cogx","pho_wphi53_cogx","pho_npb_cluster_et2","pho_npb_cluster_et3","pho_npb_cluster_et4"}},
{"isolation","isow_cone_sum",{"isow_cone_sum","isow_threshold","isow_sideband_threshold","isow_subtraction_method","isow_reconstructed_or_truth"}},
{"jets","jet_pt",{"jet_pt","jet_eta","jet_phi","jet_jes_valid","jet_emcal_frac","jet_ihcal_frac","jet_ohcal_frac","jet_rawpt","jet_time_ns","jet_time_ntowers"}},
{"truth_photons","truthpho_pt",{"truthpho_pt","truthpho_eta","truthpho_phi","truthpho_source_role","truthpho_prompt_class","truthpho_g4_photon_valid","truthpho_truth_isolation_valid","truthpho_truth_isolation_r03","truthpho_truth_isolation_r04","truthpho_analysis_signal_r03","truthpho_sample_source_role"}},
{"truth_jets","truthjet_pt",{"truthjet_pt","truthjet_eta","truthjet_phi"}},
{"photon_links","plink_reco",{"plink_reco","plink_reco_original","plink_truth","plink_class","plink_metric","plink_policy"}},
{"jet_links","jlink_reco",{"jlink_reco","jlink_truth","jlink_class","jlink_metric"}},
{"unsub_jets","jetunsub_pt",{"jetunsub_pt","jetunsub_rawpt","jetunsub_eta","jetunsub_phi","jetunsub_time_ns","jetunsub_time_ntowers"}},
};static const auto extended=[] {auto result=s;auto &fields=result.front().fields;
 for(const char *key:{"applied_id_class","applied_iso_class","applied_region","applied_preselection","applied_kinematics","analysis_isolation"})fields.push_back(key);
 fields.push_back("response_truth_index");fields.push_back("response_match_metric");
 fields.push_back("native_model_feature_state");fields.push_back("native_model_encoded_feature_state");for(int i=0;i<14;++i)fields.push_back("native_model_input_"+std::to_string(i));
 // Exact original evaluator vectors, not reconstructed raw-shape substitutes.
 // The producer vector is 11 inputs in pp; its per-row count is explicit.
 fields.push_back("producer_model_input_count");
 for(int i=0;i<14;++i)fields.push_back("producer_model_input_"+std::to_string(i));
 for(const char* key:{"pho_prodbdt_score","pho_prodbdt_finite","pho_prodbdt_applicability",
     "pho_npbcap_inputs_captured","pho_npbcap_inputs_finite","pho_npbcap_input_order_is_npb25",
     "pho_npbcap_model_available","pho_npbcap_in_producer_domain","pho_npbcap_domain_mode",
     "pho_npbcap_min_et","pho_npbcap_max_et","pho_npbcap_max_abs_eta",
     "pho_npbcap_kinematic_input_bits_available","pho_npbcap_input_et_bits","pho_npbcap_input_eta_bits"})fields.push_back(key);
 for(int i=0;i<25;++i){fields.push_back("npb_model_input_"+std::to_string(i));fields.push_back("npb_model_input_finite_"+std::to_string(i));}
 result[2].fields.push_back("response_truth_index");
 result[1].fields.push_back("diagnostic_raw_cone_sum");
 return result;}();return extended;}
struct Array {virtual ~Array()=default;virtual size_t size() const=0;virtual double at(size_t)const=0;};
template<class T> struct ReadArray final:Array {
 TTreeReaderArray<T> value;
 ReadArray(TTreeReader&r,const std::string&n):value(r,n.c_str()){}
 size_t size()const override{return value.GetSize();}double at(size_t i)const override{const auto v=value[i];if constexpr(std::is_integral_v<T> && sizeof(T)>=8){if(v>static_cast<T>(9007199254740992ULL) || (std::is_signed_v<T> && v<static_cast<T>(-9007199254740992LL)))throw std::runtime_error("candidate cache: integer cannot be represented exactly");}return v;}
};
inline std::unique_ptr<Array> reader(TTree &tree,TTreeReader&r,const std::string&n){
 auto *leaf=tree.GetLeaf(n.c_str());if(!leaf)return {};const std::string t=leaf->GetTypeName();
 if(t=="Double_t")return std::make_unique<ReadArray<Double_t>>(r,n);
 if(t=="Float_t")return std::make_unique<ReadArray<Float_t>>(r,n);
 if(t=="Int_t")return std::make_unique<ReadArray<Int_t>>(r,n);
 if(t=="UInt_t")return std::make_unique<ReadArray<UInt_t>>(r,n);
 if(t=="Long64_t")return std::make_unique<ReadArray<Long64_t>>(r,n);
 if(t=="ULong64_t")return std::make_unique<ReadArray<ULong64_t>>(r,n);
 throw std::runtime_error("candidate cache: unsupported field type "+n+":"+t);
}
// A photon-aligned view into one already captured flattened vector. Retain
// binary32 values exactly by widening to double; never evaluate a new model.
struct CapturedColumn final:Array {
 std::shared_ptr<Array> input;size_t stride,column;bool count;
 CapturedColumn(std::shared_ptr<Array> a,size_t n,size_t i,bool c=false):input(std::move(a)),stride(n),column(i),count(c){}
 size_t size()const override{if(input->size()%stride)throw std::runtime_error("ragged captured model vector");return input->size()/stride;}
 double at(size_t i)const override{if(i>=size())throw std::runtime_error("captured model index outside photons");return count?double(stride):input->at(i*stride+column);}
};
struct Table {
 std::unique_ptr<TTree> tree;Long64_t event=0,index=0,inputIndex=0;
 std::map<std::string,Double_t> values;std::map<std::string,Long64_t> present;
 Table(TDirectory&out,const TableSpec&s){out.cd();tree=std::make_unique<TTree>(("table_candidates_"+s.name).c_str(),"PhotonJetRecordsV1");
  tree->SetAutoFlush(-8*1024*1024);tree->Branch("event",&event,"event/L");tree->Branch("index",&index,"index/L");tree->Branch("input_index",&inputIndex,"input_index/L");
  for(const auto&n:s.fields){tree->Branch(n.c_str(),&(values[n]=0),(n+"/D").c_str());tree->Branch(("has_"+n).c_str(),&(present[n]=0),("has_"+n+"/L").c_str());}}
 void write(){tree->GetDirectory()->cd();tree->Write("",TObject::kOverwrite);tree->SetDirectory(nullptr);}
};
// A complete original-file timing reference, never an accepted-event median.
// State: 0 not required (SIM), 1 complete and finite, 2 complete with no finite
// samples, 3 partial source (unavailable for exact timing QA replay).
struct QaSourceWitness { double mbdReference=std::numeric_limits<double>::quiet_NaN(); Long64_t state=0,observed=0,expected=0,finite=0; };
struct Writer {
 std::unique_ptr<TTree> events;std::map<std::string,std::unique_ptr<Table>> tables;
 Long64_t rows=0,run=0,seq=0,entry=0,file=0;Double_t centrality=0,vz=0,weight=0;
 Double_t mbdT0=std::numeric_limits<double>::quiet_NaN();Long64_t hasMbdT0=0;
 Double_t truthVz=std::numeric_limits<double>::quiet_NaN();Long64_t hasTruthVz=0;
 Long64_t inputPhotons=0,retainedPhotons=0,droppedBelowFloor=0,droppedInvalidEt=0;
 Long64_t responseMatchWitness=0,originalJetLinks=0;
 Double_t etFloor=0;
 Double_t qaFilterPhoton=std::numeric_limits<double>::quiet_NaN(),qaTruthJet04=-1;
 Long64_t hasQaFilterPhoton=0,hasQaTruthJet04=0;
 std::map<Long64_t,QaSourceWitness> qaSources;
 WeightFactors factors;
 std::map<std::string,Long64_t> counts;
 Double_t jetRadius;
 explicit Writer(TDirectory&out,double radius):jetRadius(radius){radiusIndex(radius);out.cd();events=std::make_unique<TTree>("table_candidates_events","Accepted events including those with no nominal leader");events->SetAutoFlush(-8*1024*1024);
  events->Branch("jet_radius",&jetRadius,"jet_radius/D");
  events->Branch("run",&run,"run/L");events->Branch("seq",&seq,"seq/L");events->Branch("entry",&entry,"entry/L");events->Branch("file",&file,"file/L");
  events->Branch("centrality",&centrality,"centrality/D");events->Branch("vtx_z",&vz,"vtx_z/D");events->Branch("weight",&weight,"weight/D");
  events->Branch("mbd_t0",&mbdT0,"mbd_t0/D");events->Branch("has_mbd_t0",&hasMbdT0,"has_mbd_t0/L");
  events->Branch("truth_vtx_z",&truthVz,"truth_vtx_z/D");events->Branch("has_truth_vtx_z",&hasTruthVz,"has_truth_vtx_z/L");
  events->Branch("qa_filter_photon_pt",&qaFilterPhoton,"qa_filter_photon_pt/D");
  events->Branch("has_qa_filter_photon_pt",&hasQaFilterPhoton,"has_qa_filter_photon_pt/L");
  events->Branch("qa_leading_truthjet04_pt",&qaTruthJet04,"qa_leading_truthjet04_pt/D");
  events->Branch("has_qa_leading_truthjet04_pt",&hasQaTruthJet04,"has_qa_leading_truthjet04_pt/L");
  events->Branch("weight_stored",&factors.stored,"weight_stored/D");events->Branch("weight_base",&factors.base,"weight_base/D");
  events->Branch("weight_vertex",&factors.vertex,"weight_vertex/D");events->Branch("weight_centrality",&factors.centrality,"weight_centrality/D");
  events->Branch("photon_input_count",&inputPhotons,"photon_input_count/L");events->Branch("photon_retained_count",&retainedPhotons,"photon_retained_count/L");
  events->Branch("photon_below_floor_count",&droppedBelowFloor,"photon_below_floor_count/L");events->Branch("photon_invalid_et_count",&droppedInvalidEt,"photon_invalid_et_count/L");
  events->Branch("reco_et_min_gev",&etFloor,"reco_et_min_gev/D");
  events->Branch("response_match_witness",&responseMatchWitness,"response_match_witness/L");
  events->Branch("jet_link_input_count",&originalJetLinks,"jet_link_input_count/L");
  for(const auto&k:{"jet","jetunsub","truthjet"})for(const auto&r:{"02","03","04"}){const std::string n=std::string(k)+r+"_count";events->Branch(n.c_str(),&(counts[n]=0),(n+"/L").c_str());}
  for(const auto&s:specs())tables[s.name]=std::make_unique<Table>(out,s);
 }
 void scientificContext(Long64_t source,const photonjet_science::Json& context) {
  auto rows=photonjet_science::candidateContexts(*events);const auto key=std::to_string(source);
  photonjet_science::require(!rows.contains(key)||rows.at(key)==context,"candidate source context changed");rows[key]=context;
  photonjet_science::candidateContexts(*events,rows);
 }
 void write(){
  // The existing event pass measures the source reference before event cuts.
  // Append these columns only after every source has ended; no extra source
  // scan, event reweighting, planner change or merge specialization is needed.
  if(events->GetBranch("qa_file_mbd_reference"))throw std::runtime_error("candidate cache already finalized");
  QaSourceWitness witness;
  std::vector<TBranch*> columns;
  columns.push_back(events->Branch("qa_file_mbd_reference",&witness.mbdReference,"qa_file_mbd_reference/D"));
  columns.push_back(events->Branch("qa_file_mbd_reference_state",&witness.state,"qa_file_mbd_reference_state/L"));
  columns.push_back(events->Branch("qa_source_observed_entries",&witness.observed,"qa_source_observed_entries/L"));
  columns.push_back(events->Branch("qa_source_expected_entries",&witness.expected,"qa_source_expected_entries/L"));
  columns.push_back(events->Branch("qa_source_finite_mbd_entries",&witness.finite,"qa_source_finite_mbd_entries/L"));
  auto* files=events->GetBranch("file");
  for(Long64_t i=0;i<events->GetEntries();++i){
   if(files->GetEntry(i)<=0||!qaSources.count(file))throw std::runtime_error("candidate cache lacks finalized source QA witness");
   witness=qaSources.at(file);for(auto* column:columns)if(column->Fill()<0)throw std::runtime_error("candidate QA witness column write failed");
  }
  events->GetDirectory()->cd();events->Write("",TObject::kOverwrite);events->SetDirectory(nullptr);
  // Columns refer to the local witness only while finalizing, never afterward.
  for(auto* column:columns)column->ResetAddress();
  for(auto&[n,t]:tables)t->write();
 }
};
struct Source {
 Writer &out;std::map<std::string,std::unique_ptr<Array>> arrays;
 bool simulation=false;
 Long64_t qaExpected=0,qaObserved=0,qaFirst=-1,qaLast=-1;
 bool qaContiguous=true,qaFinished=false;
 std::vector<double> qaMbdSamples;
 std::unique_ptr<TTreeReaderValue<Double_t>> mbdT0;
 std::unique_ptr<TTreeReaderValue<Double_t>> truthVz;
 Source(Writer&w,TTree&tree,TTreeReader&r,bool sim=false):out(w),simulation(sim),qaExpected(tree.GetEntries()){
  if(auto*leaf=tree.GetLeaf("mbd_t0")){
   if(std::string(leaf->GetTypeName())!="Double_t")throw std::runtime_error("candidate cache: mbd_t0 type differs");
   mbdT0=std::make_unique<TTreeReaderValue<Double_t>>(r,"mbd_t0");
  }
  if(auto*leaf=tree.GetLeaf("truth_vtx_z")){
   if(std::string(leaf->GetTypeName())!="Double_t")throw std::runtime_error("candidate cache: truth_vtx_z type differs");
   truthVz=std::make_unique<TTreeReaderValue<Double_t>>(r,"truth_vtx_z");
  }
  for(const auto&s:specs())for(const auto&n:s.fields)if(!arrays.count(n))arrays[n]=reader(tree,r,sourceField(n,out.jetRadius));
  const auto capture=[&](const std::string& branch,const std::string& prefix,size_t stride,const std::string& count=""){
   std::shared_ptr<Array> raw=reader(tree,r,branch);if(!raw)return;
   if(stride==0||stride>25)throw std::runtime_error("unsupported captured model stride");
   for(size_t i=0;i<stride;++i)arrays[prefix+std::to_string(i)]=std::make_unique<CapturedColumn>(raw,stride,i);
   if(!count.empty())arrays[count]=std::make_unique<CapturedColumn>(raw,stride,0,true);
  };
  if(tree.GetLeaf("pho_prodbdt_inputs")){
   auto* schema=tree.GetCurrentFile()->Get<TObjString>("meta/schema");
   if(!schema)throw std::runtime_error("producer model vector needs source schema");
   const auto metadata=nlohmann::json::parse(schema->GetString().Data());
   const size_t stride=metadata.at("strides").at("pho_prodbdt_inputs");
   if(stride==0||stride>14)throw std::runtime_error("producer model vector exceeds retained capacity");
   capture("pho_prodbdt_inputs","producer_model_input_",stride,"producer_model_input_count");
  }
  capture("pho_npbcap_inputs","npb_model_input_",25);
  capture("pho_npbcap_input_finite","npb_model_input_finite_",25);
  for(const auto&k:{"jet","jetunsub","truthjet"})for(const auto&rad:{"02","03","04"}){const std::string n=std::string(k)+rad+"_pt";if(!arrays.count(n))arrays[n]=reader(tree,r,n);}
  for(const char *n:{"isoc_pho","isoc_delta_r","isoc_delta_eta","isoc_calibrated_energy","isoc_sub1_energy"})arrays[n]=reader(tree,r,n);
  arrays["embdiag_stitch_photon_pt"]=reader(tree,r,"embdiag_stitch_photon_pt");
 }
 // Call once for every source event, before ownership, vertex or trigger cuts.
 void observeSource(Long64_t entry){
  if(qaFinished)throw std::runtime_error("source QA observation after finalization");
  if(!qaObserved)qaFirst=entry;else if(entry!=qaLast+1)qaContiguous=false;
  qaLast=entry;++qaObserved;
  if(!simulation&&mbdT0&&std::isfinite(**mbdT0))qaMbdSamples.push_back(**mbdT0);
 }
 void finishSource(Long64_t file){
  if(qaFinished||out.qaSources.count(file))throw std::runtime_error("duplicate source QA finalization");
  QaSourceWitness witness;witness.observed=qaObserved;witness.expected=qaExpected;witness.finite=qaMbdSamples.size();
  if(!simulation){
   const bool complete=qaContiguous&&qaObserved==qaExpected&&(qaExpected==0||(qaFirst==0&&qaLast==qaExpected-1));
   witness.state=complete?(qaMbdSamples.empty()?2:1):3;
   if(witness.state==1){auto middle=qaMbdSamples.begin()+qaMbdSamples.size()/2;
    std::nth_element(qaMbdSamples.begin(),middle,qaMbdSamples.end());witness.mbdReference=*middle;}
  }
  out.qaSources.emplace(file,witness);qaFinished=true;
 }
 void append(Long64_t run,Long64_t seq,Long64_t entry,Long64_t file,double cent,double vz,double weight,const std::function<double(size_t)>& selectedScore={},
             const WeightFactors& factors={},const std::function<PhotonValues(size_t)>& applied={}){
  out.run=run;out.seq=seq;out.entry=entry;out.file=file;out.centrality=cent;out.vz=vz;out.weight=weight;
  out.hasMbdT0=bool(mbdT0);out.mbdT0=mbdT0?**mbdT0:std::numeric_limits<double>::quiet_NaN();
  out.hasTruthVz=bool(truthVz);out.truthVz=truthVz?**truthVz:std::numeric_limits<double>::quiet_NaN();
  out.qaFilterPhoton=std::numeric_limits<double>::quiet_NaN();out.hasQaFilterPhoton=0;
  if(const auto& filter=arrays.at("embdiag_stitch_photon_pt"))for(size_t i=0;i<filter->size();++i){
   const double value=filter->at(i);if(!std::isfinite(value))continue;
   if(out.hasQaFilterPhoton&&value!=out.qaFilterPhoton)throw std::runtime_error("candidate QA has inconsistent filter photons");
   out.qaFilterPhoton=value;out.hasQaFilterPhoton=1;
  }
  const auto& jets04=arrays.at("truthjet04_pt");out.hasQaTruthJet04=simulation&&bool(jets04);out.qaTruthJet04=-1;
  if(out.hasQaTruthJet04)for(size_t i=0;i<jets04->size();++i)
   if(std::isfinite(jets04->at(i)))out.qaTruthJet04=std::max(out.qaTruthJet04,jets04->at(i));
  out.factors=factors;
  if(applied && (!std::isfinite(factors.base)||!std::isfinite(factors.vertex)||!std::isfinite(factors.centrality)||
      weight!=factors.base*factors.vertex*factors.centrality))throw std::runtime_error("candidate cache weight factors do not reproduce event weight");
  for(auto&[n,count]:out.counts){const auto key=n.substr(0,n.size()-6)+"_pt";count=arrays.at(key)?arrays.at(key)->size():0;}
  const auto event=out.events->GetEntries();
  if(!std::isfinite(minRecoEt)||minRecoEt<0||minRecoEt>12)throw std::runtime_error("candidate cache ET floor outside 0..12 GeV");
  const auto&et=arrays.at("pho_et");const size_t nPhotons=et?et->size():0;
  std::vector<Long64_t> remap(nPhotons,-1);std::vector<size_t> retained;
  out.inputPhotons=nPhotons;out.droppedBelowFloor=0;out.droppedInvalidEt=0;out.etFloor=minRecoEt;
  for(size_t i=0;i<nPhotons;++i){const auto value=et->at(i);
   if(!std::isfinite(value)){++out.droppedInvalidEt;continue;}
   if(value<minRecoEt){++out.droppedBelowFloor;continue;}
   remap[i]=retained.size();retained.push_back(i);
  }out.retainedPhotons=retained.size();
  out.originalJetLinks=arrays.at("jlink_reco")?arrays.at("jlink_reco")->size():0;
  out.responseMatchWitness=simulation ? 1 : 0;
  response::EventMatches matches;
  if(out.responseMatchWitness){
   const auto required=[&](const char* name)->const Array& {
    const auto& a=arrays.at(name);if(!a)throw std::runtime_error(std::string("compact response cache lacks ")+name);return *a;};
   struct Floats {const Array& a;size_t GetSize()const{return a.size();}double operator[](size_t i)const{return a.at(i);}};
   struct Ints {const Array& a;size_t GetSize()const{return a.size();}int operator[](size_t i)const{
    const double x=a.at(i);if(!std::isfinite(x)||x!=std::floor(x)||x<std::numeric_limits<int>::min()||x>std::numeric_limits<int>::max())
      throw std::runtime_error("compact response cache has an invalid match index/class");return static_cast<int>(x);}};
   const auto counts=[&](const char* kind){std::array<size_t,3> result{};int k=0;for(const char* radius:{"02","03","04"})
    result[k++]=out.counts.at(std::string(kind)+radius+"_count");return result;};
   const auto ranges=response::JetRanges::of(counts("jet"),counts("jetunsub"),counts("truthjet"),radiusIndex(out.jetRadius));
   if(ranges.nJet!=required("jet_rawpt").size()||ranges.nUnsub!=required("jetunsub_rawpt").size())
     throw std::runtime_error("compact response cache raw jet arrays differ");
   for(const char* name:{"plink_truth","plink_class","plink_metric"})if(required(name).size()!=required("plink_reco").size())
     throw std::runtime_error("compact response cache photon links differ");
   for(const char* name:{"jlink_truth","jlink_class","jlink_metric"})if(required(name).size()!=required("jlink_reco").size())
     throw std::runtime_error("compact response cache jet links differ");
   matches=response::eventMatches(nPhotons,required("truthpho_pt").size(),
    Ints{required("plink_reco")},Ints{required("plink_truth")},Ints{required("plink_class")},Floats{required("plink_metric")},required("plink_reco").size(),
    ranges,Floats{required("jet_rawpt")},Floats{required("jetunsub_rawpt")},
    Ints{required("jlink_reco")},Ints{required("jlink_truth")},Ints{required("jlink_class")},Floats{required("jlink_metric")},required("jlink_reco").size());
  }
  const auto&iso=arrays.at("isow_cone_sum");
  if(iso&&iso->size()!=0&&iso->size()!=2*nPhotons)throw std::runtime_error("candidate cache requires two ordered reco isolation witnesses per photon");
  // Archival diagnostic ONLY. Justin explicitly excludes raw isolation from
  // downstream selections, tuning, response and uncertainty calculations.
  std::vector<double> ue(2*nPhotons,0.);
  bool hasUE=true;size_t nConstituents=arrays.at("isoc_pho")?arrays.at("isoc_pho")->size():0;
  for(const char *key:{"isoc_pho","isoc_delta_r","isoc_delta_eta","isoc_calibrated_energy","isoc_sub1_energy"})
   hasUE=hasUE && arrays.at(key) && arrays.at(key)->size()==nConstituents;
  hasUE=hasUE && arrays.at("pho_eta") && arrays.at("pho_eta")->size()==nPhotons;
  if(hasUE)for(size_t q=0;q<nConstituents;++q){
   const double p=arrays.at("isoc_pho")->at(q);if(!std::isfinite(p)||p!=std::floor(p)||p<0||p>=nPhotons)continue;
   const auto i=static_cast<size_t>(p);const double dr=arrays.at("isoc_delta_r")->at(q);
   const double ch=std::cosh(arrays.at("pho_eta")->at(i)+arrays.at("isoc_delta_eta")->at(q));
   const double removed=(arrays.at("isoc_calibrated_energy")->at(q)-arrays.at("isoc_sub1_energy")->at(q))/ch;
   for(size_t r=0;r<2;++r)if(dr<(r==0?.3:.4))ue[2*i+r]+=std::isfinite(removed)?removed:std::numeric_limits<double>::quiet_NaN();
  }
  auto recoIndex=[&](size_t i){const double value=arrays.at("plink_reco")->at(i);
   if(!std::isfinite(value)||value!=std::floor(value)||value< -1||value>=static_cast<double>(nPhotons))throw std::runtime_error("candidate cache invalid reconstructed photon link");
   return static_cast<Long64_t>(value);};
  for(const auto&s:specs()){auto&t=*out.tables.at(s.name);const auto&count=arrays.at(s.count);const size_t n=count?count->size():0;
   for(const auto&key:s.fields)if(arrays.at(key)&&arrays.at(key)->size()!=0&&arrays.at(key)->size()!=n)throw std::runtime_error("candidate cache: ragged field "+key);
   std::vector<size_t> selected;
   if(out.responseMatchWitness && (s.name=="jet_links"||s.name=="photon_links"||s.name=="unsub_jets")){}
   else if(s.name=="photons")selected=retained;
   else if(s.name=="isolation"){if(n)for(auto i:retained){selected.push_back(2*i);selected.push_back(2*i+1);}}
   else if(s.name=="photon_links"){for(size_t i=0;i<n;++i){const auto old=recoIndex(i);if(old<0||remap[old]>=0)selected.push_back(i);}}
   else{selected.resize(n);std::iota(selected.begin(),selected.end(),0);}
   Long64_t dense=0;
   for(size_t i:selected){const auto derived=s.name=="photons" && applied ? applied(i) : PhotonValues{};
    t.event=event;t.index=dense++;t.inputIndex=i;for(const auto&key:s.fields){const auto&a=arrays.at(key);t.present[key]=a&&a->size()==n;t.values[key]=t.present[key]?a->at(i):std::numeric_limits<double>::quiet_NaN();
     if(key=="selected_score"&&selectedScore){t.present[key]=1;t.values[key]=selectedScore(i);}
     if(derived.count(key)){t.present[key]=1;t.values[key]=derived.at(key);}
     if(out.responseMatchWitness && key=="response_truth_index"){
      t.present[key]=1;t.values[key]=s.name=="photons"?matches.photonTruth.at(i):matches.jetTruth.at(i);
     }
     if(out.responseMatchWitness && s.name=="photons" && key=="response_match_metric"){
      t.present[key]=1;t.values[key]=matches.photonMetric.at(i);
     }
     if(key=="diagnostic_raw_cone_sum"){
      const auto &method=arrays.at("isow_subtraction_method");
      if(method && method->size()==n && (method->at(i)==1 || (method->at(i)==2 && hasUE))){
       const double raw=iso->at(i)+(method->at(i)==2?ue[i]:0.);
       t.present[key]=std::isfinite(raw);t.values[key]=raw;
      }
     }
     if(s.name=="photon_links"&&(key=="plink_reco"||key=="plink_reco_original")){const auto old=recoIndex(i);t.present[key]=1;t.values[key]=key=="plink_reco_original"?old:(old<0?-1:remap[old]);}
    }if(++out.rows>maxRows)throw std::runtime_error("candidate cache row ceiling exceeded; no truncation");t.tree->Fill();}}
  if(++out.rows>maxRows)throw std::runtime_error("candidate cache row ceiling exceeded; no truncation");
  out.events->Fill();
 }
};
} }
