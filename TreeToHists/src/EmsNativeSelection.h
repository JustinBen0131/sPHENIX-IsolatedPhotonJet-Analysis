/* Native EMS selection on original eventTree inputs.
   Reuses the captured-input reader and source-bound NPB/ID evaluator. Only
   EMS enters selection ET; named EMR/JER remain response coordinates.
   The required contract declares retained analysis isolation witnesses as a
   mechanical, unadopted choice. Raw AuAu isolation is never substituted. */
#pragma once
#include "EnergyVariation.h"
#include "PhotonEnergyInputs.h"

namespace photonjet { namespace ems_native {
using Json=SelectionJson;
inline void require(bool test,const char* why) {
  if(!test)throw std::runtime_error(std::string("native EMS: ")+why);
}

inline void validateContract(const Json& contract) {
  require(contract.is_object() && contract.value("schema","")=="EmsNativeSelectionContractV1","explicit EMS contract is required");
  require(contract.at("source").at("sha256")=="cabc3a5eb081ce77e0e1d83936e04bcc427bed6f5be5066183983edd4dc32105",
          "PPG18 source identity differs");
  require(contract.at("cone").at("rule")=="retain_captured_analysis_witness" &&
          contract.at("cone").at("status")=="MECHANICAL_UNADOPTED","no adopted cone prescription is installed");
  require(contract.at("cone").at("raw_auau_isolation")=="ARCHIVAL_ONLY","raw AuAu isolation cannot enter EMS selection");
  require(contract.contains("selection") && contract.contains("models") && contract.contains("pp_witness_line"),
          "selection/model/witness bindings are missing");
}

// A view remains exactly nominal until assigned for the current event.
// Clearing before every reader advance prevents rejected events seeing the
// previous event's varied arrays. Original preselection flags stay archival.
template<class Original,class T> class Overlay {
  const Original& original_;
  std::vector<T> changed_;
  bool varied_=false;
public:
  explicit Overlay(const Original& original):original_(original){}
  size_t GetSize()const{return varied_?changed_.size():original_.GetSize();}
  T operator[](size_t i)const{return varied_?changed_.at(i):original_[i];}
  int saved(size_t i)const{return original_.saved(i);}
  void clear(){changed_.clear();varied_=false;}
  template<class Values>void assign(const Values& values) {
    changed_.resize(values.GetSize());for(size_t i=0;i<changed_.size();++i)changed_[i]=values[i];varied_=true;
  }
};

struct Evaluation {
  energy_selection::Values<double> et,score,npb,cone,threshold,sideband;
  energy_selection::Values<int> state,preselection;
};

class Runtime {
  Profile profile_;
  Json contract_;
  std::unique_ptr<energy_selection::Evaluator> evaluator_;
  static energy_selection::ModelBinding binding(const Json& value) {
    return {value.at("file"),value.at("key"),value.at("sha256")};
  }
public:
  Runtime(const Profile& profile,const Json& contract,
      const std::function<std::string(const std::string&)>& digest):profile_(profile),contract_(contract) {
    require(profile.simulation && profile.configured,"resolved SIM profile is required");
    validateContract(contract_);
    require(Json::parse(profile.selectionJson)==contract_.at("selection"),"resolved selection differs from contract");
    require(profile.system!=System::auau || profile.isolationSubtraction==2,"AuAu requires the existing SUB1 witness");
    const auto& line=contract_.at("pp_witness_line");
    energy_selection::IsolationWitnessLine witness{{line.at("intercept"),line.at("slope")},line.at("gap")};
    const auto pp=profile.system==System::pp?binding(contract_.at("models").at("pp")):energy_selection::ModelBinding{};
    evaluator_=std::make_unique<energy_selection::Evaluator>(profile,pp,binding(contract_.at("models").at("npb")),witness,digest);
  }

  template<class Photons>
  Evaluation evaluate(const energy_selection::NativeInputs& inputs,const Photons& original,
      double vertex,double centrality,double emScale) const {
    Evaluation out;const auto n=original.GetSize();
    out.et.resize(n);out.score.assign(n,std::numeric_limits<double>::quiet_NaN());out.npb.assign(n,-1.);
    out.state.assign(n,0);out.preselection.assign(n,0);out.cone.resize(2*n);out.threshold.resize(2*n);out.sideband.resize(2*n);
    for(size_t i=0;i<n;++i) {
      out.et[i]=emScale==1.?original[i]:emScale*original[i];
      auto input=inputs.photon(i,vertex,centrality);
      require(input.et==original[i] || (!std::isfinite(input.et)&&!std::isfinite(original[i])),"captured photon identity differs");
      // Deliberately retain the actual prescribed analysis witness. There is
      // no inference of EM content or UE scaling from one cone-sum scalar.
      for(size_t r=0;r<2;++r){out.cone[2*i+r]=input.prescribedCone[r];out.threshold[2*i+r]=input.nominalThreshold[r];out.sideband[2*i+r]=input.nominalSideband[r];}
      if(!std::isfinite(original[i]))continue;
      const auto result=evaluator_->evaluate(input,out.et[i]);
      out.score[i]=result.score;out.npb[i]=result.npbScore;out.state[i]=result.scoreState;out.preselection[i]=result.preselection;
      for(size_t r=0;r<2;++r){out.threshold[2*i+r]=result.threshold[r];out.sideband[2*i+r]=result.sideband[r];}
    }
    return out;
  }
  const Json& contract()const{return contract_;}
};

inline Json parseContract(const char* text) {
  require(text && *text,"explicit EMS contract is required");auto value=Json::parse(text);validateContract(value);return value;
}
inline bool active(const Profile& profile,const response::Variation& variation) {
  return profile.simulation && response::emsSelectionResponse(variation);
}
} }
