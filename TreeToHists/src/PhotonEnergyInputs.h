/* Native and retained-record inputs for PhotonEnergySelection.
   Adapters preserve original photon indices and exact captured float inputs.
   Cone modification is an explicit caller-supplied prescription; this module
   never opens the archival raw AuAu diagnostic. */
#pragma once
#include "PhotonEnergySelection.h"
#include "CandidateCache.h"

namespace photonjet { namespace energy_selection {
inline std::vector<std::string> capturedFields(const Profile& p) {
  std::vector<std::string> names{"pho_et1","pho_e11_over_e33","pho_e32_over_e35",
    "pho_npbcap_input_order_is_npb25","pho_npbcap_model_available","pho_npbcap_domain_mode",
    "pho_npbcap_min_et","pho_npbcap_max_et","pho_npbcap_max_abs_eta",
    "pho_npbcap_inputs_captured","pho_npbcap_inputs_finite"};
  for(int j=0;j<25;++j){names.push_back("npb_model_input_"+std::to_string(j));names.push_back("npb_model_input_finite_"+std::to_string(j));}
  const bool pp=p.system==System::pp;
  names.push_back(pp?"producer_model_input_count":"native_model_feature_state");
  for(int j=0;j<(pp?11:14);++j)names.push_back(std::string(pp?"producer_model_input_":"native_model_input_")+std::to_string(j));
  return names;
}

template<class T>struct Values:std::vector<T> {
  using std::vector<T>::vector;
  size_t GetSize()const{return this->size();}
};

struct Counters {
  std::uint64_t candidates=0,incomingId=0,outgoingId=0,incomingNpb=0,outgoingNpb=0,preselectionChanged=0;
};
struct Context {
  const Evaluator& evaluator;
  std::function<std::array<double,2>(const Input&,double)> conePrescription;
  std::string conePrescriptionIdentity;
  mutable Counters counters;
  Result apply(Input& input,double et)const {
    if(!conePrescription||conePrescriptionIdentity.empty())throw std::runtime_error("energy selection requires an explicit cone prescription");
    const auto nominal=evaluator.evaluate(input,input.et);
    input.prescribedCone=conePrescription(input,et);
    auto varied=evaluator.evaluate(input,et);
    ++counters.candidates;counters.incomingId+=varied.idApplicable&&!nominal.idApplicable;
    counters.outgoingId+=nominal.idApplicable&&!varied.idApplicable;
    counters.incomingNpb+=varied.npbApplicable&&!nominal.npbApplicable;
    counters.outgoingNpb+=nominal.npbApplicable&&!varied.npbApplicable;
    counters.preselectionChanged+=varied.preselection!=nominal.preselection;
    return varied;
  }
};

class NativeInputs {
  const Profile& profile_;
  const PhotonScores& scores_;
  std::map<std::string,std::unique_ptr<candidate_cache::Array>> arrays_;
  size_t producerStride_=0;
  double get(const std::string& name,size_t i)const {
    const auto found=arrays_.find(name);
    if(found==arrays_.end()||!found->second||i>=found->second->size())throw std::runtime_error("native energy input missing/misaligned: "+name);
    return found->second->at(i);
  }
public:
  NativeInputs(TTree& tree,TTreeReader& reader,const Profile& profile,const PhotonScores& scores)
      :profile_(profile),scores_(scores) {
    auto fields=capturedFields(profile);
    for(const char* name:{"pho_et","pho_eta","pho_phi","isow_cone_sum","isow_threshold","isow_sideband_threshold",
                          "isow_subtraction_method","isow_reconstructed_or_truth","pho_prodbdt_inputs","pho_npbcap_inputs","pho_npbcap_input_finite"})fields.push_back(name);
    for(const auto& name:fields)arrays_[name]=candidate_cache::reader(tree,reader,name);
    auto* schema=tree.GetCurrentFile()->Get<TObjString>("meta/schema");
    if(!schema)throw std::runtime_error("native energy inputs need source schema");
    producerStride_=SelectionJson::parse(schema->GetString().Data()).at("strides").at("pho_prodbdt_inputs");
    if(profile.system==System::pp&&producerStride_!=11)throw std::runtime_error("native pp input stride differs");
  }
  Input photon(size_t i,double vertex,double centrality)const {
    Input input;input.et=get("pho_et",i);input.eta=get("pho_eta",i);input.phi=get("pho_phi",i);
    input.vertex=vertex;input.centrality=centrality;input.acceptedEvent=true;
    for(const auto& name:capturedFields(profile_)) {
      const auto& array=arrays_.at(name);if(array&&i<array->size())input.fields[name]=array->at(i);
    }
    const auto vector=[&](const char* branch,const char* prefix,size_t stride){
      const auto& a=arrays_.at(branch);
      if(!a)return;
      if(a->size()!=arrays_.at("pho_et")->size()*stride)throw std::runtime_error("native classifier vector alignment differs");
      for(size_t j=0;j<stride;++j)input.fields[std::string(prefix)+std::to_string(j)]=a->at(i*stride+j);
    };
    vector("pho_npbcap_inputs","npb_model_input_",25);vector("pho_npbcap_input_finite","npb_model_input_finite_",25);
    if(profile_.system==System::pp){vector("pho_prodbdt_inputs","producer_model_input_",producerStride_);input.fields["producer_model_input_count"]=producerStride_;}
    else {const auto features=scores_.cacheInputs(i);input.fields.insert(features.begin(),features.end());}
    for(size_t r=0;r<2;++r){const auto at=2*i+r;
      input.nominalThreshold[r]=get("isow_threshold",at);input.nominalSideband[r]=get("isow_sideband_threshold",at);
      input.prescribedCone[r]=get("isow_cone_sum",at);input.method[r]=get("isow_subtraction_method",at);input.role[r]=get("isow_reconstructed_or_truth",at);}
    return input;
  }
};

inline Input recordInput(const std::map<std::string,double>& photon,
                         const std::array<std::map<std::string,double>,2>& isolation,double vertex,double centrality) {
  Input input;input.et=photon.at("pho_et");input.eta=photon.at("pho_eta");input.phi=photon.at("pho_phi");
  input.vertex=vertex;input.centrality=centrality;input.acceptedEvent=true;input.fields=photon;
  for(size_t r=0;r<2;++r){const auto& q=isolation[r];
    input.nominalThreshold[r]=q.at("isow_threshold");input.nominalSideband[r]=q.at("isow_sideband_threshold");
    input.prescribedCone[r]=q.at("isow_cone_sum");input.method[r]=q.at("isow_subtraction_method");input.role[r]=q.at("isow_reconstructed_or_truth");}
  return input;
}
} }
