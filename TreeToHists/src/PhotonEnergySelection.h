/* PhotonEnergySelection: shared photon selection after an energy variation.
   Native fills and record replay supply the same captured evaluator inputs.
   This owner recomputes NPB and ID, routes on varied ET, and evaluates the
   existing PhotonJet selection. It never reconstructs towers or invents a
   cone-energy variation: the caller must supply its prescribed cone values.
   Model files must be verified against their bindings before construction. */
#pragma once
#include "PhotonJet.h"
#include "Preselection.h"
#include <functional>

namespace photonjet { namespace energy_selection {

struct ModelBinding {
  std::string file, key, sha256;
};

// The witness line is explicit because a stored threshold alone does not
// determine its dependence on ET. The source configuration owns these units.
struct IsolationWitnessLine {
  Line threshold;
  double gap;
};

struct Input {
  double et, eta, phi, vertex, centrality;
  bool acceptedEvent;
  std::map<std::string,double> fields;
  std::array<double,2> nominalThreshold, nominalSideband, prescribedCone;
  std::array<int,2> method, role;
  double get(const std::string& key) const {
    const auto found=fields.find(key);
    if(found==fields.end()) throw std::runtime_error("energy selection missing captured field: "+key);
    return found->second;
  }
};

struct Result {
  double score=std::numeric_limits<double>::quiet_NaN(), npbScore=-1.;
  int scoreState=0, preselection=0;
  bool idApplicable=false, npbApplicable=false, kinematics=false;
  PhotonClass classification;
  std::array<double,2> threshold, sideband;
};

class Evaluator {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
  using BDT=TMVA::Experimental::RBDT;
#else
  using BDT=TMVA::Experimental::RBDT<>;
#endif
  struct Cell {
    double low, high, centLow, centHigh;
    std::unique_ptr<BDT> model;
  };
  Profile profile_;
  IsolationWitnessLine witness_;
  SelectionJson domain_;
  std::unique_ptr<BDT> npb_,pp_;
  std::vector<Cell> cells_;

  static void require(bool test,const std::string& why) {
    if(!test)throw std::runtime_error("energy selection: "+why);
  }
  static std::vector<float> inputs(const Input& in,const std::string& prefix,size_t count,double variedEt) {
    std::vector<float> values(count);
    for(size_t j=0;j<count;++j) {
      const double value=in.get(prefix+std::to_string(j));
      values[j]=static_cast<float>(value);
      require(std::isfinite(value) && static_cast<double>(values[j])==value,
              "model input is not an exact finite captured float");
    }
    require(values[0]==static_cast<float>(in.et),"captured ET differs from original photon ET");
    values[0]=static_cast<float>(variedEt);
    require(std::isfinite(values[0]),"varied ET overflows model input");
    return values;
  }
  static double predict(const BDT& model,const std::vector<float>& values) {
    const auto result=model.Compute(values);
    require(result.size()==1 && std::isfinite(result[0]) && result[0]>=0 && result[0]<=1,
            "model returned an invalid probability");
    return result[0];
  }
  static bool close(double x,double y) {
    return std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=1e-12*std::max(1.,std::abs(y));
  }

public:
  // digest is the existing caller's SHA256 implementation. Every cell is
  // checked here, including the NPB model used before final photon ID.
  Evaluator(const Profile& profile,const ModelBinding& pp,const ModelBinding& npb,
            const IsolationWitnessLine& witness,
            const std::function<std::string(const std::string&)>& digest)
      :profile_(profile),witness_(witness) {
    require(profile_.configured,"energy variation requires a resolved selection");
    require(bool(digest),"model digest verifier is required");
    const auto load=[&](const ModelBinding& binding) {
      require(binding.sha256.size()==64 && digest(binding.file)==binding.sha256,"model SHA256 differs");
      return std::make_unique<BDT>(binding.key,binding.file);
    };
    require(npb.sha256==canonicalPreselection().at("npb_model_sha256").get<std::string>(),"NPB model is not the selected preselection model");
    npb_=load(npb);
    const auto selection=SelectionJson::parse(profile_.selectionJson);
    if(profile_.system==System::pp) {
      require(!profile_.optionalID,"pp optional-model routing is not implemented");
      require(pp.sha256==selection.at("systems").at("pp").at("model_sha256").get<std::string>(),"pp model differs from selection");
      pp_=load(pp);
      require(std::isfinite(witness_.threshold.intercept)&&std::isfinite(witness_.threshold.slope)&&
              std::isfinite(witness_.gap)&&witness_.gap>=0,"invalid pp witness line");
    } else {
      require(profile_.nativeID,"AuAu energy variation requires native selected-model inputs");
      const auto& native=selection.at("systems").at("auau").at("optional_model").at("native_evaluation");
      const auto& contract=native.at("feature_contract");
      require(contract.at("count")==14 && contract.at("order")[0]=="cluster_Et",
              "unsupported AuAu feature contract");
      domain_=native.at("domain");
      for(const auto& cell:native.at("cells")) {
        ModelBinding binding{cell.at("file"),cell.at("key"),cell.at("sha256")};
        cells_.push_back({cell.at("et")[0],cell.at("et")[1],
                          cell.at("centrality")[0],cell.at("centrality")[1],load(binding)});
      }
      require(!cells_.empty(),"native model has no routing cells");
    }
  }

  Result evaluate(const Input& in,double variedEt) const {
    Result result;
    result.threshold=in.nominalThreshold;result.sideband=in.nominalSideband;
    require(std::isfinite(in.et)&&std::isfinite(variedEt),"nonfinite original or varied photon ET");
    if(!in.acceptedEvent)return result;

    // The NPB producer's domain compares float inputs and has an inclusive
    // upper edge and a 1e-6 GeV lower-edge tolerance. Preserve it exactly.
    require(in.get("pho_npbcap_input_order_is_npb25")==1 &&
            in.get("pho_npbcap_model_available")==1 && in.get("pho_npbcap_domain_mode")==1,
            "NPB capture has no qualified bounded model");
    const float npbEt=static_cast<float>(variedEt),npbEta=static_cast<float>(in.get("npb_model_input_1"));
    const double low=in.get("pho_npbcap_min_et"),high=in.get("pho_npbcap_max_et"),maxEta=in.get("pho_npbcap_max_abs_eta");
    result.npbApplicable=(!std::isfinite(low)||npbEt>low||std::fabs(npbEt-low)<1e-6F) &&
      (!std::isfinite(high)||npbEt<=high) && (!std::isfinite(maxEta)||std::fabs(npbEta)<=maxEta);
    if(result.npbApplicable) {
      require(in.get("pho_npbcap_inputs_captured")==1 && in.get("pho_npbcap_inputs_finite")==1,
              "NPB inputs invalid inside varied domain");
      for(size_t j=0;j<25;++j)require(in.get("npb_model_input_finite_"+std::to_string(j))==1,"NPB finite bit differs");
      result.npbScore=predict(*npb_,inputs(in,"npb_model_input_",25,variedEt));
    }
    result.preselection=newPPG12(result.npbScore,in.get("pho_et1"),in.get("pho_e11_over_e33"),in.get("pho_e32_over_e35"));

    if(profile_.system==System::pp) {
      // The supported analysis is 15 <= ET < 35 GeV. Out-of-analysis stored
      // producer scores remain archival; this evaluator never extrapolates.
      require(profile_.photonEtMin>=15 && profile_.photonEtMax<=35,"pp energy selection exceeds qualified analysis domain");
      result.idApplicable=variedEt>=profile_.photonEtMin && variedEt<profile_.photonEtMax &&
        std::isfinite(in.eta)&&std::abs(in.eta)<profile_.photonAbsEtaMax &&
        std::isfinite(in.vertex)&&std::abs(in.vertex)<profile_.vertexAbsZMax;
      if(result.idApplicable) {
        require(in.get("producer_model_input_count")==11,"pp model input count differs");
        result.score=predict(*pp_,inputs(in,"producer_model_input_",11,variedEt));
      }
      if(profile_.isolationFromWitness)for(size_t r=0;r<2;++r) {
        const double nominal=witness_.threshold.at(in.et);
        require(close(in.nominalThreshold[r],nominal)&&close(in.nominalSideband[r],nominal+witness_.gap),
                "stored isolation witness does not match supplied source line");
        result.threshold[r]=witness_.threshold.at(variedEt);
        result.sideband[r]=result.threshold[r]+witness_.gap;
      }
    } else {
      require(!std::isfinite(in.centrality)||in.centrality==std::round(in.centrality),"noninteger native centrality");
      result.idApplicable=variedEt>=domain_.at("et_min").get<double>() && variedEt<domain_.at("et_max_exclusive").get<double>() &&
        std::abs(in.eta)<domain_.at("eta_abs_max").get<double>() && std::abs(in.vertex)<domain_.at("vertex_abs_max_cm").get<double>() &&
        in.centrality>=domain_.at("centrality_min").get<double>() && in.centrality<domain_.at("centrality_max_exclusive").get<double>();
      if(result.idApplicable) {
        require(in.get("native_model_feature_state")==1,"invalid native features inside varied domain");
        const Cell* selected=nullptr;
        for(const auto& cell:cells_)if(variedEt>=cell.low && variedEt<cell.high && in.centrality>=cell.centLow && in.centrality<cell.centHigh) {
          require(!selected,"native model routing overlaps");selected=&cell;
        }
        require(selected,"native model routing leaves a gap");
        result.score=predict(*selected->model,inputs(in,"native_model_input_",14,variedEt));
      }
    }
    result.scoreState=result.idApplicable?1:0;
    result.kinematics=photonKinematics(profile_,in.vertex,variedEt,in.eta,in.phi);
    // Reuse PhotonJet's ID lines, sidebands, cone parity and region semantics.
    // Raw AuAu diagnostic isolation is never an input to this interface.
    const std::array<double,2> scores{result.score,result.score};
    const std::array<int,2> state{result.scoreState,result.scoreState},unused{0,0};
    result.classification=photonClass(profile_,0,variedEt,in.centrality,scores,state,in.prescribedCone,
        result.threshold,result.sideband,in.method,in.role,unused,unused);
    return result;
  }
};
} }
