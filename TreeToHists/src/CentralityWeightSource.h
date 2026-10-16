#pragma once
// One alternative: remove only the AuAu centrality factor from the response.
// The frozen vertex correction, source ownership and source normalization stay.
// This is a model-dependence comparison, not a finite-statistics covariance band.
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <string>
namespace photonjet { namespace centrality_source {
using Json=nlohmann::json;
inline void require(bool good,const std::string& why){if(!good)throw std::runtime_error("centrality source: "+why);}
inline double responseWeight(double nominal,double centrality,bool off,bool auau){
 if(!off)return nominal; // Preserve nominal bytes and operation order.
 require(std::isfinite(nominal)&&nominal>=0&&std::isfinite(centrality)&&centrality>0,"invalid finite weight factors");
 require(auau,"centrality off is AuAu-only");
 const double result=nominal/centrality;
 require(std::isfinite(result),"off weight overflow");return result;
}
inline double recordedWeight(double stored,double base,double vertex,double centrality,double normalization,bool off){
 require(std::isfinite(base)&&base>=0&&std::isfinite(vertex)&&vertex>=0&&std::isfinite(normalization)&&normalization>0,"invalid source/vertex factors");
 require(stored==base*vertex*centrality,"stored weight is not exactly one factor product");
 // The same existing responseWeight operation; normalization is applied once.
 const double result=responseWeight(stored,centrality,off,true)*normalization;
 require(std::isfinite(result),"normalized weight overflow");return result;
}
inline void validateTarget(const Json& photon,const Json& inclusive,const Json& selected,const Json& dataSelection,const Json& counts,const Json& qualification){
 const std::string identity="6f324813586b5e1e0357f7feb82caa439d15d1dec1c50baf3aec67ee0a249871";
 require(selected.at("selection_sha256")==identity,"wrong code_p selection");
 require(photon.at("schema")=="AuAuSimWeightFreezeV1"&&inclusive.at("schema")=="AuAuSimWeightFreezeV1","wrong freeze schema");
 require(photon.at("family")=="photonjet"&&inclusive.at("family")=="inclusivejet","wrong SIM families");
 for(const auto* f:{&photon,&inclusive}){
  require(f->at("status")=="QUALIFIED"&&!f->at("test_only").get<bool>(),"unqualified freeze");
  require(f->at("selection_sha256")==identity&&f->at("data_package_sha256")==qualification.at("data").at("sha256"),"freeze target differs");
  require(f->at("centrality").at("target")=="tight_isolated_data_photons_15_le_et_lt_35_gev"&&f->at("centrality").at("weights").size()==16,"wrong target population");
  require(f->at("vertex").at("coordinate")=="truth_vtx_z","wrong vertex coordinate");
  for(const auto&w:f->at("centrality").at("weights"))require(std::isfinite(w.get<double>())&&w.get<double>()>0,"invalid frozen map");
 }
 require(photon.at("vertex")==inclusive.at("vertex")&&photon.at("data_trigger")==inclusive.at("data_trigger"),"different shared vertex/trigger");
 require(photon.at("centrality").at("edges_percent")==inclusive.at("centrality").at("edges_percent"),"different centrality axes");
 for(int i=0;i<=16;++i)require(photon.at("centrality").at("edges_percent").at(i)==5*i,"wrong centrality edge");
 for(const char* k:{"schema","version","preselection","timing_reference"})require(selected.at(k)==dataSelection.at(k),std::string("DATA selection differs: ")+k);
 require(selected.at("systems").at("auau")==dataSelection.at("systems").at("auau"),"DATA AuAu selection differs");
 require(counts==qualification.at("target_counts")&&counts.size()==16,"selected DATA target counts changed");
 double sum=0;for(const auto&x:counts){require(std::isfinite(x.get<double>())&&x.get<double>()>0,"invalid target cell");sum+=x.get<double>();}
 require(sum==29985&&qualification.at("target_candidates")==29985,"target total differs");
}
} }
