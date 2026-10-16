#pragma once
/* Streaming measurement replay through the same accumulator as native fills.
   Caller binds the package, selection and complete normalization denominator.
   This kernel does not grant package, response or scientific acceptance. */
#include "MeasurementAccumulator.h"
#include "CandidateCache.h"
#include "PhotonEnergyInputs.h"
#include "NativeEnergyVariation.h"
#include "PurityInputs.h"
#include "ReplayScientificContext.h"
#include <TLeaf.h>
#include <limits>
#include <set>

namespace photonjet { namespace candidate_replay {
inline void require(bool good,const std::string& why) {
  if(!good)throw std::runtime_error("candidate replay: "+why);
}
template<class T> struct Values:std::vector<T> {
  using std::vector<T>::vector;
  size_t GetSize()const{return this->size();}
};
inline TLeaf* leaf(TTree& tree,const std::string& name,const std::string& type) {
  auto* result=tree.GetLeaf(name.c_str());
  require(result && result->GetTypeName()==type,"missing or mistyped "+name);
  tree.SetBranchStatus(name.c_str(),true);
  return result;
}
inline Long64_t integer(double value,const std::string& name) {
  require(std::isfinite(value) && value==std::floor(value) &&
          value>=std::numeric_limits<int>::min() && value<=std::numeric_limits<int>::max(),
          "inexact integer "+name);
  return static_cast<Long64_t>(value);
}
struct Slice {
  Long64_t source=-1;
  std::map<std::string,std::pair<Long64_t,Long64_t>> rows;
};
inline std::pair<Long64_t,Long64_t> rangeOf(TTree& table,const std::string& name,const Slice* slice) {
  const auto range=slice?slice->rows.at(name):std::make_pair(Long64_t(0),table.GetEntries());
  require(range.first>=0 && range.second>=range.first && range.second<=table.GetEntries(),"invalid source table range "+name);
  if(slice)require(slice->source>=0,"invalid sliced source identity");
  return range;
}
inline void sourceIdentities(const std::vector<std::string>& identities,const Slice* slice) {
  const auto check=[](const std::string& id){require(id.size()==64 && id.find_first_not_of("0123456789abcdef")==std::string::npos,"invalid variation source identity");};
  if(slice){require(slice->source>=0 && slice->source<static_cast<Long64_t>(identities.size()),"slice source outside identity inventory");check(identities.at(slice->source));}
  else for(const auto& id:identities)check(id);
}
class Rows {
  TTree* tree_;
  Long64_t row_=0,end_=0;
  TLeaf *event_,*index_,*input_,*source_=nullptr;
  std::map<std::string,std::pair<TLeaf*,TLeaf*>> fields_;
  std::set<std::string> optional_;
 public:
  Rows(TFile& file,const std::string& prefix,const std::string& name,
       const std::vector<std::string>& fields,bool packaged,const std::vector<std::string>& optionalFields={},const Slice* slice=nullptr) {
    tree_=file.Get<TTree>((prefix+name).c_str());require(tree_,"missing "+name);
    const auto range=rangeOf(*tree_,name,slice);row_=range.first;end_=range.second;
    tree_->SetBranchStatus("*",false);
    event_=leaf(*tree_,"event","Long64_t");index_=leaf(*tree_,"index","Long64_t");
    input_=leaf(*tree_,"input_index","Long64_t");
    if(packaged)source_=leaf(*tree_,"source","Long64_t");
    for(const auto& key:fields)fields_[key]={leaf(*tree_,key,"Double_t"),leaf(*tree_,"has_"+key,"Long64_t")};
    for(const auto& key:optionalFields){fields_[key]={leaf(*tree_,key,"Double_t"),leaf(*tree_,"has_"+key,"Long64_t")};optional_.insert(key);}
    if(row_<end_)require(tree_->GetEntry(row_)>0,"unreadable first "+name+" row");
  }
  std::vector<std::map<std::string,double>> take(Long64_t event,Long64_t source,Long64_t maxCount) {
    std::vector<std::map<std::string,double>> result;
    Long64_t previousInput=-1;
    while(row_<end_) {
      const auto key=event_->GetValueLong64();require(key>=event,"regressing event reference");
      if(key!=event)break;
      require(index_->GetValueLong64()==static_cast<Long64_t>(result.size()),"candidate index gap");
      require(input_->GetValueLong64()>previousInput,"original candidate order changed");previousInput=input_->GetValueLong64();
      require(!source_ || source_->GetValueLong64()==source,"foreign source reference");
      require(static_cast<Long64_t>(result.size())<maxCount,"per-event candidate ceiling");
      std::map<std::string,double> values{{"__input_index",static_cast<double>(previousInput)}};
      for(const auto& [name,leaves]:fields_){
        if(optional_.count(name)&&leaves.second->GetValueLong64()==0)continue;
        require(leaves.second->GetValueLong64()==1,"unavailable required field "+name);
        values[name]=leaves.first->GetValue();
      }
      result.push_back(std::move(values));++row_;
      if(row_<end_)require(tree_->GetEntry(row_)>0,"unreadable candidate row");
    }
    return result;
  }
  void finish()const{require(row_==end_,"orphan rows after final event");}
};
struct Summary {
  Long64_t cacheEvents=0,measurementEvents=0;
  std::vector<Long64_t> measurementSources;
};
/* Full, finite event weights are applied once here. The complete per-sample
   pb factor is deliberately an argument; fragment-level normalization is not
   inferred from the selected cached event count. */
inline Summary measurement(TFile& input,const std::string& prefix,const Profile& profile,
    const std::vector<CentralityClass>& classes,Histograms& histograms,EventTable& output,
    double normalization,Long64_t sourceCount,bool timing=false,Long64_t maxCandidates=2000000,
    const response::Variation* energyVariation=nullptr,const energy_selection::Context* energySelection=nullptr,
    const std::vector<std::string>& sourceIdentities={},purity_inputs::Book* purity=nullptr,const Slice* slice=nullptr) {
  const response::Variation chosenVariation=energyVariation?*energyVariation:response::Variation{};
  const bool varied=energyVariation||energySelection,truthNeeded=energy_selection::needsTruth(chosenVariation);
  const bool purityTruth=purity && profile.simulation,photonTruthNeeded=truthNeeded||purityTruth;
  require(!energy_selection::changesPhotons(chosenVariation)||energySelection,"photon measurement replay requires whole-selection evaluator");
  require(!truthNeeded||profile.simulation,"truth-based energy resolution requires SIM");
  if(varied){require(sourceIdentities.size()==static_cast<size_t>(sourceCount),"measurement variation needs original source identities");
    candidate_replay::sourceIdentities(sourceIdentities,slice);}
  require(prefix=="records/candidates/" || prefix=="table_candidates_","unsupported table layout");
  require(profile.system!=System::auau || (profile.isolationRadiusIndex==0 && profile.isolationSubtraction==2),
          "AuAu analysis requires R03 SUB1 isolation");
  require(profile.system!=System::pp || (profile.isolationRadiusIndex==1 && profile.isolationSubtraction==1),
          "pp analysis requires R04 topocluster isolation");
  require(std::isfinite(normalization) && normalization>0 && sourceCount>0 && maxCandidates>0,"invalid replay envelope");
  require(profile.simulation || normalization==1.,"DATA normalization must be one");
  const bool packaged=prefix=="records/candidates/";
  auto* events=input.Get<TTree>((prefix+"events").c_str());require(events,"missing event table");
  for(const auto& context:photonjet_science::replayScientificContexts(*events,profile,slice?slice->source:-1,
      nlohmann::json::parse(response::variationJson(chosenVariation)),nlohmann::json::object(),normalization,timing)) {
    histograms.scientificContext(context);output.scientificContext(context);
  }
  events->SetBranchStatus("*",false);
  std::map<std::string,TLeaf*> ei,ed;
  for(const char* name:{"run","seq","entry","file","photon_input_count","photon_retained_count",
      "photon_below_floor_count","photon_invalid_et_count","jet02_count","jet03_count","jet04_count"})ei[name]=leaf(*events,name,"Long64_t");
  if(packaged)ei["source"]=leaf(*events,"source","Long64_t");
  if(photonTruthNeeded){ei["response_match_witness"]=leaf(*events,"response_match_witness","Long64_t");
    const std::string name="truthjet"+candidate_cache::radiusSuffix(profile.jetRadius)+"_count";ei[name]=leaf(*events,name,"Long64_t");}
  for(const char* name:{"centrality","vtx_z","weight","weight_base","weight_vertex","weight_centrality","reco_et_min_gev","jet_radius"})
    ed[name]=leaf(*events,name,"Double_t");
  std::vector<std::string> photonFields={"pho_et","pho_eta","pho_phi","selected_score","applied_preselection"};
  if(!profile.optionalID)photonFields.push_back(profile.stateBranch);
  if(profile.augmentedClasses){photonFields.push_back("pho_id_class");photonFields.push_back("pho_iso_class");}
  if(timing)photonFields.push_back("pho_time_ns");
  if(photonTruthNeeded)photonFields.push_back("response_truth_index");
  Rows photons(input,prefix,"photons",photonFields,packaged,energySelection?energy_selection::capturedFields(profile):std::vector<std::string>{},slice);
  // Raw diagnostic isolation is intentionally never opened by this reader.
  Rows isolation(input,prefix,"isolation",{"isow_cone_sum","isow_threshold","isow_sideband_threshold",
      "isow_subtraction_method","isow_reconstructed_or_truth"},packaged,{},slice);
  std::vector<std::string> jetFields{"jet_pt","jet_eta","jet_phi"};
  if(chosenVariation.emScale!=1.)jetFields.push_back("jet_emcal_frac");
  if(truthNeeded)jetFields.push_back("response_truth_index");
  Rows jets(input,prefix,"jets",jetFields,packaged,{},slice);
  std::unique_ptr<Rows> truths,truthJets;
  if(photonTruthNeeded){std::vector<std::string> fields{"truthpho_pt"};if(purityTruth)fields.push_back("truthpho_source_role");
    truths=std::make_unique<Rows>(input,prefix,"truth_photons",fields,packaged,std::vector<std::string>{},slice);}
  if(truthNeeded)truthJets=std::make_unique<Rows>(input,prefix,"truth_jets",std::vector<std::string>{"truthjet_pt"},packaged,std::vector<std::string>{},slice);
  Summary summary;Long64_t previousSource=-1,previousEntry=-1,previousFile=-1;
  const auto eventRange=rangeOf(*events,"events",slice);
  for(Long64_t e=eventRange.first;e<eventRange.second;++e){
    require(events->GetEntry(e)>0,"unreadable event row");
    const auto i=[&](const std::string& key){return ei.at(key)->GetValueLong64();};
    const auto d=[&](const std::string& key){return ed.at(key)->GetValue();};
    const Long64_t source=packaged?i("source"):i("file");
    require(source>=0 && source<sourceCount && source>=previousSource,"invalid source inventory/order");
    require(!slice || source==slice->source,"event outside requested source slice");
    require(i("file")>=0 && i("file")<=std::numeric_limits<int>::max() && i("entry")>=0,"invalid original event identity");
    if(source==previousSource && i("file")==previousFile)require(i("entry")>previousEntry,"duplicate or regressing original event");
    previousSource=source;previousEntry=i("entry");previousFile=i("file");
    require(i("run")>=0 && i("run")<=std::numeric_limits<int>::max(),"run outside native integer range");
    require(i("photon_input_count")>=0 && i("photon_retained_count")>=0 && i("photon_below_floor_count")>=0 && i("photon_invalid_et_count")>=0 &&
        i("photon_input_count")==i("photon_retained_count")+i("photon_below_floor_count")+i("photon_invalid_et_count"),"photon retention accounting differs");
    require(std::isfinite(d("reco_et_min_gev")) && d("reco_et_min_gev")>=0 && profile.photonEtMin>=d("reco_et_min_gev"),"replay requests discarded photon energies");
    const double weight=d("weight"),product=d("weight_base")*d("weight_vertex")*d("weight_centrality");
    require(std::isfinite(weight) && std::isfinite(product) && weight==product && std::isfinite(weight*normalization),"weight factors differ or overflow");
    auto p=photons.take(e,source,maxCandidates),q=isolation.take(e,source,2*maxCandidates),j=jets.take(e,source,maxCandidates);
    std::vector<std::map<std::string,double>> ts,tjs;
    if(photonTruthNeeded){require(i("response_match_witness")==1,"measurement physics inputs lack verified matches");ts=truths->take(e,source,maxCandidates);}
    if(truthNeeded){tjs=truthJets->take(e,source,maxCandidates);
      require(tjs.size()==static_cast<size_t>(i("truthjet"+candidate_cache::radiusSuffix(profile.jetRadius)+"_count")),"truth jet cardinality differs");}
    require(p.size()==static_cast<size_t>(i("photon_retained_count")) && q.size()==2*p.size(),"photon or isolation cardinality differs");
    require(d("jet_radius")==profile.jetRadius,"record jet radius differs from replay profile");
    const auto jetCount=i("jet"+candidate_cache::radiusSuffix(profile.jetRadius)+"_count");
    require(jetCount>=0 && j.size()==static_cast<size_t>(jetCount),"selected-radius jet cardinality differs");
    for(size_t k=0;k<p.size();++k){
      const double original=p[k].at("__input_index");
      require(original<i("photon_input_count") && q[2*k].at("__input_index")==2*original &&
          q[2*k+1].at("__input_index")==2*original+1,"isolation witness no longer belongs to its photon");
    }
    for(size_t k=0;k<j.size();++k)require(j[k].at("__input_index")==static_cast<double>(k),"unfiltered jet identity changed");
    Values<double> et,eta,phi,score,cone,threshold,sideband,pt,je,jp;
    Values<int> state,method,role,id,iso;std::vector<int> accepted,visible;
    for(const auto& row:p){
      et.push_back(row.at("pho_et"));eta.push_back(row.at("pho_eta"));phi.push_back(row.at("pho_phi"));score.push_back(row.at("selected_score"));
      require(std::isfinite(et.back()) && et.back()>=d("reco_et_min_gev"),"cached photon below its recorded floor");
      const auto pre=integer(row.at("applied_preselection"),"preselection");require(pre==0 || pre==1,"invalid preselection value");
      const bool keep=!timing || dataTimeCut(row.at("pho_time_ns"));accepted.push_back(pre && keep);visible.push_back(keep);
      state.push_back(profile.optionalID?0:integer(row.at(profile.stateBranch),"score state"));
      id.push_back(profile.augmentedClasses?integer(row.at("pho_id_class"),"ID class"):0);
      iso.push_back(profile.augmentedClasses?integer(row.at("pho_iso_class"),"isolation class"):0);
    }
    for(const auto& row:q){cone.push_back(row.at("isow_cone_sum"));threshold.push_back(row.at("isow_threshold"));sideband.push_back(row.at("isow_sideband_threshold"));
      method.push_back(integer(row.at("isow_subtraction_method"),"isolation method"));role.push_back(integer(row.at("isow_reconstructed_or_truth"),"isolation role"));}
    for(const auto& row:j){pt.push_back(row.at("jet_pt"));je.push_back(row.at("jet_eta"));jp.push_back(row.at("jet_phi"));}
    ++summary.cacheEvents;reportProgress(summary.cacheEvents);
    const double cent=d("centrality"),vz=d("vtx_z");
    require(std::isfinite(vz) && std::abs(vz)<profile.vertexAbsZMax,"cached event outside original vertex selection");
    size_t k=0;if(profile.system==System::auau)while(k<classes.size() && !(classes[k].low<=cent && cent<classes[k].high))++k;
    if(k==classes.size())continue;
    const double w=weight*normalization;histograms.h1("events_"+classes[k].name).Fill(.5,w);
    if(varied && !p.empty()) {
      const double scale=chosenVariation.photonScale*chosenVariation.emScale;
      require(!i("photon_below_floor_count") || (chosenVariation.photonSmear==0 && chosenVariation.photonResolution=="off" && d("reco_et_min_gev")*scale<=profile.photonEtMin),
              "discarded photons could enter the varied measurement");
      const auto total=i("photon_input_count");require(total>=0&&total<=maxCandidates,"measurement photon ceiling");
      Values<double> original(total,std::numeric_limits<double>::quiet_NaN());
      for(const auto& row:p)original.at(integer(row.at("__input_index"),"original photon index"))=row.at("pho_et");
      const auto key=response::textKey(sourceIdentities.at(source));
      energy_selection::ChangedMomenta momenta(original,pt,chosenVariation,key,i("entry"));
      if(truthNeeded) {
        Values<double> tp,tj;std::vector<long> pm(total,-1),jm(j.size(),-1);
        for(size_t t=0;t<ts.size();++t){require(ts[t].at("__input_index")==static_cast<double>(t),"truth photon order changed");tp.push_back(ts[t].at("truthpho_pt"));}
        for(size_t t=0;t<tjs.size();++t){require(tjs[t].at("__input_index")==static_cast<double>(t),"truth jet order changed");tj.push_back(tjs[t].at("truthjet_pt"));}
        for(const auto& row:p){const auto match=integer(row.at("response_truth_index"),"photon truth match");
          require(match>=-1&&match<static_cast<Long64_t>(tp.size()),"photon match outside truth");pm.at(integer(row.at("__input_index"),"photon index"))=match;}
        std::set<long> used;
        for(size_t t=0;t<j.size();++t){jm[t]=integer(j[t].at("response_truth_index"),"jet truth match");
          require(jm[t]>=-1&&jm[t]<static_cast<Long64_t>(tj.size()),"jet match outside truth");require(jm[t]<0||used.insert(jm[t]).second,"duplicate jet truth match");}
        if(chosenVariation.photonResolution!="off"){momenta.photons.smearPhotons(tp,pm,chosenVariation,key,i("entry"));momenta.photons.scaleAll(scale);}
        if(chosenVariation.jetResolution!="off"){momenta.jets.smearJets(tj,jm,chosenVariation,profile.jetRadius,key,i("entry"));momenta.jets.scaleAll(chosenVariation.jetScale);}
      }
      if(chosenVariation.emScale!=1.){Values<double> em;for(const auto& row:j)em.push_back(row.at("jet_emcal_frac"));momenta.jets.scaleEmFraction(em,chosenVariation.emScale);}
      for(size_t t=0;t<j.size();++t)pt[t]=momenta.jets[t];
      for(size_t t=0;t<p.size();++t){et[t]=momenta.photons[integer(p[t].at("__input_index"),"original photon index")];
        if(energySelection){auto input=energy_selection::recordInput(p[t],{q[2*t],q[2*t+1]},vz,cent);const auto result=energySelection->apply(input,et[t]);
          score[t]=result.score;state[t]=result.scoreState;accepted[t]=result.preselection&&visible[t];
          for(size_t r=0;r<2;++r){cone[2*t+r]=input.prescribedCone[r];threshold[2*t+r]=result.threshold[r];sideband[2*t+r]=result.sideband[r];}}
      }
    }
    if(purity) {
      std::vector<int> signal(p.size(),0);
      if(purityTruth) {
        for(size_t t=0;t<ts.size();++t)require(ts[t].at("__input_index")==static_cast<double>(t),"purity truth photon order changed");
        for(size_t t=0;t<p.size();++t) {
          const auto match=integer(p[t].at("response_truth_index"),"purity photon truth match");
          require(match>=-1 && match<static_cast<Long64_t>(ts.size()),"purity photon match outside truth collection");
          if(match>=0)signal[t]=integer(ts[match].at("truthpho_source_role"),"purity truth role")==kSignalSourceRole;
        }
      }
      purity->accumulate(profile,classes[k],vz,cent,w,et,eta,phi,score,accepted,state,cone,threshold,sideband,method,role,id,iso,pt,je,jp,signal,
          profile.simulation?1.:configuredDataJetScale);
    }
    const auto before=output.run.size();
    accumulateMeasurement(profile,classes[k],static_cast<int>(i("run")),i("seq"),i("entry"),static_cast<int>(i("file")),cent,vz,w,
        histograms,&output,et,eta,phi,score,accepted,visible,state,cone,threshold,sideband,method,role,id,iso,pt,je,jp,
        profile.simulation?1.:configuredDataJetScale);
    if(output.run.size()!=before){require(output.run.size()==before+1,"unexpected event multiplicity");summary.measurementSources.push_back(source);++summary.measurementEvents;}
  }
  photons.finish();isolation.finish();jets.finish();if(truths)truths->finish();if(truthJets)truthJets->finish();
  return summary;
}
} } // namespace photonjet::candidate_replay
