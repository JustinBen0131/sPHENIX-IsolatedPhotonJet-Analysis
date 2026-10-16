#pragma once
#include <cmath>
#include <stdexcept>
#include <TTreeReaderArray.h>
#include <nlohmann/json.hpp>
namespace photonjet {
inline nlohmann::json canonicalPreselection(){return nlohmann::json::parse(R"npb({"schema":"CanonicalNPBPreselectionV1","name":"newPPG12","score_branch":"pho_npb_score","npb_model_sha256":"d6086dadac534013cda15cdfb69c1683776d3456d9e439903589653e8ac19eab","finite_fields":["pho_npb_score","pho_et1","pho_e11_over_e33","pho_e32_over_e35"],"score_gt":0.5,"et1_gt":0.6,"et1_lt":1.0,"e11_over_e33_lt":0.98,"e32_over_e35_gt":0.8,"e32_over_e35_lt":1.0,"weta_cut":null,"tree_flag":"pho_presel_active","tree_flag_role":"cross_check_only","applies_to":"pp_and_auau_data_photonjet_inclusivejet"})npb");}
inline bool newPPG12(double score,double et1,double e11e33,double e32e35){
 return std::isfinite(score)&&std::isfinite(et1)&&std::isfinite(e11e33)&&std::isfinite(e32e35)&&
        score>0.5 && et1>0.6 && et1<1.0 && e11e33<0.98 && e32e35>0.8 && e32e35<1.0;
}
/* The producer flag is a cross-check only. Candidate ordering is unchanged. */
class PhotonPreselection {
 TTreeReaderArray<double> score_,et1_,e11e33_,e32e35_;
 TTreeReaderArray<int> saved_;
 public:
 explicit PhotonPreselection(TTreeReader&r):score_(r,"pho_npb_score"),et1_(r,"pho_et1"),
 e11e33_(r,"pho_e11_over_e33"),e32e35_(r,"pho_e32_over_e35"),saved_(r,"pho_presel_active"){}
 size_t GetSize()const{return score_.GetSize();}
 int saved(size_t i)const{if(i>=saved_.GetSize())throw std::runtime_error("preselection flag size differs");return saved_[i];}
 int operator[](size_t i)const{
  const auto n=score_.GetSize();
  if(i>=n||et1_.GetSize()!=n||e11e33_.GetSize()!=n||e32e35_.GetSize()!=n||saved_.GetSize()!=n)
   throw std::runtime_error("canonical preselection primitives are not aligned");
  return newPPG12(score_[i],et1_[i],e11e33_[i],e32e35_[i])?1:0;
 }
};
}
