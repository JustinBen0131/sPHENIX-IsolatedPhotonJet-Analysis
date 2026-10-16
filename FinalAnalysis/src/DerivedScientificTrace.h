#pragma once
#include "FinalAnalysis.h"
#include "RootScientificTrace.h"

namespace fa::derived_science {
using photonjet_science::Trace;
using photonjet_science::Json;
inline Trace join(const std::string &operation, std::initializer_list<Trace> inputs, const Json &parameters=Json::object()) {
  return Trace::combine(operation,std::vector<Trace>(inputs),parameters);
}
inline std::vector<double> effectiveScale(const std::vector<double> &values,int n) {
  return values.empty() ? std::vector<double>(n,1.0) : values;
}
inline ScientificLineage corrections(const Measured &m,const Window &w,const SystemSettings &sys) {
  ScientificLineage out;
  const auto measured=[&](const char *name){return lineageOf(m.lineage,name);};
  const auto simulated=[&](const char *name){return lineageOf(w.lineage,name);};
  const int n=m.abcd.GetNcols();
  std::vector<Trace> signalInputs{measured("abcd"),simulated("leakage")};
  if (sys.backgroundCorrelationFromQA) signalInputs.push_back(lineageOf(sys.lineage,"r_bkg"));
  if (sys.purityScaleFromQA) signalInputs.push_back(lineageOf(sys.lineage,"purity_scale"));
  auto signal=Trace::combine("leakage_aware_abcd_physical_root",signalInputs,
    {{"background_correlation",effectiveScale(sys.backgroundCorrelation,n)},
     {"purity_scale",effectiveScale(sys.purityScale,n)},
     {"root_rule","physical_root_nearest_naive"}});
  signal=withSettings(signal,sys.parameterOrigins,{"BackgroundCorrelation","PurityScale"});
  out["purity"]=signal.then("signal_over_region_A");
  out["transfer"]=signal.then("background_A_over_C");
  out["photons"]=signal.then("divide_efficiency",{{"scale",effectiveScale(sys.efficiencyScale,n)}});
  out["photons"]=withSettings(out["photons"],sys.parameterOrigins,{"EfficiencyScale"});
  auto pairsA=measured("pairsA"),pairsC=measured("pairsC");
  if (sys.subtractSideband) {
    const auto transfer=join("combinatoric_window_over_sideband",{simulated("combinatoric"),simulated("sidebandJets")},{{"empty_denominator",0.0}});
    pairsA=join("subtract_sideband_times_transfer",{pairsA,measured("sidebandA"),transfer});
    pairsC=join("subtract_sideband_times_transfer",{pairsC,measured("sidebandC"),transfer});
  }
  out["scaledC"]=join("scale_region_C",{pairsC,out.at("transfer")});
  out["purityCorrected"]=join("subtract_C_divide_leakage_denominator",{pairsA,pairsC,out.at("transfer"),simulated("leakage")});
  out["pairs"]=out.at("purityCorrected");
  out["combinatoric"]=signal.then("no_simulation_template_subtracted");
  if (sys.subtractTemplate) {
    out["combinatoric"]=join("simulation_template_per_matched_photon_times_signal",{simulated("combinatoric"),simulated("leakage"),signal});
    out["pairs"]=join("subtract_simulation_template",{out.at("pairs"),out.at("combinatoric")});
  }
  out["pairs"]=out.at("pairs").then("divide_efficiency",{{"scale",effectiveScale(sys.efficiencyScale,n)}});
  for(const auto *key:{"pairs","combinatoric"}) {
    out[key]=withSettings(out.at(key),sys.parameterOrigins,{"SubtractCombinatoric"});
    if(sys.subtractCombinatoric)out[key]=withSettings(out.at(key),sys.parameterOrigins,{"Combinatoric.Method"});
    if(sys.subtractSideband)out[key]=withSettings(out.at(key),sys.parameterOrigins,{"Combinatoric.Shape"});
  }
  out["pairs"]=withSettings(out.at("pairs"),sys.parameterOrigins,{"EfficiencyScale"});
  return out;
}
inline Trace response(const Window &w,bool photon) {
  const auto t=[&](const char *name){return lineageOf(w.lineage,name);};
  return photon ? join("photon_bayes_problem",{t("photonResponse"),t("photonMisses"),t("photonTruth"),t("photonFakes")})
                : join("joint_bayes_problem",{t("pairResponse"),t("pairMisses"),t("pairTruth"),t("pairFakes"),t("detector"),t("boundaryByXj"),t("boundaryFixed")});
}
// The automatic rule always scans k=1..12 (Unfolding.cc); the separately
// configured display scan length is not an input to that choice.
inline Trace iterationRule(Trace inputs,const Settings &s,bool photon,bool fittedScan) {
  Json parameters={{"data_model",s.dataToyModel},{"data_toys",s.ruleDataToys},
    {"response_toys",s.ruleResponseToys},{"data_seed",s.ruleDataSeed},
    {"response_seed",s.ruleResponseSeed},{"maximum_iterations",12}};
  if(!photon) {
    parameters["prior"]=s.rulePrior;
    parameters["reported_high"]=s.reportedHigh;
    if(s.rulePrior=="fitted")parameters["prior_form"]=s.priorForm;
    if(fittedScan) {
      parameters["prior_fraction"]=s.priorFraction;parameters["fake_model"]=s.fakeModel;
      parameters["fit_range"]=s.fitRange;
    }
  }
  auto trace=inputs.then(photon?"photon_iteration_rule":"joint_iteration_rule",parameters);
  trace=withSettings(trace,s.parameterOrigins,{"Uncertainty.DataToyModel","IterationRule.DataToys","IterationRule.DataSeed","IterationRule.ResponseToys","IterationRule.ResponseSeed"});
  if(!photon) {
    trace=withSettings(trace,s.parameterOrigins,{"IterationRule.Prior","Measurement.ReportedXjMax"});
    if(s.rulePrior=="fitted")trace=withSettings(trace,s.parameterOrigins,{"Prior.Form"});
    if(fittedScan)trace=withSettings(trace,s.parameterOrigins,{"Prior.Fraction","Prior.FakeModel","Prior.FitRange"});
  }
  return trace;
}
inline Trace iterationChoice(Trace inputs,const SystemSettings &sys,const char *key,int requested,const Trace &scan) {
  auto trace=withSettings(inputs,sys.parameterOrigins,{key});
  return requested==0?join("scan_selected_iterations",{trace,scan}):trace;
}
/* A covariance computed from the same factors multiplied in another order (the ratio's b_i b_i b_j b_j) can differ
   from its transpose in the last bit; that rounding is evened out before the strict record check. A larger
   difference is still an error. */
inline Json covarianceRecord(const TMatrixD &matrix,const std::string &method) {
  const int n=matrix.GetNrows();
  if(matrix.GetNcols()!=n)throw std::runtime_error("covariance record of a non-square matrix");
  std::vector<double> values(matrix.GetMatrixArray(),matrix.GetMatrixArray()+matrix.GetNoElements());
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j) {
    double &a=values[i*n+j],&b=values[j*n+i];
    if(a==b)continue;
    const double scale=std::max({std::fabs(a),std::fabs(b),std::sqrt(std::fabs(values[i*n+i]*values[j*n+j]))});
    if(!(std::fabs(a-b)<=1e-12*scale))throw std::runtime_error("nonsymmetric covariance beyond rounding");
    a=b=0.5*(a+b);
  }
  return photonjet_science::covariance("measured",values,n,method);
}
inline void write(TObject *object,const ScientificLineage &lineage,const std::map<std::string,Json> &covariances={}) {
  const std::string name=object->GetName();
  photonjet_science::attachTrace(object,lineageOf(lineage,name));
  const auto covariance=covariances.find(name);
  photonjet_science::setTag(object,"physics_covariance",covariance==covariances.end()
    ? photonjet_science::covariance("not_computed",{},0,"this producer exports bin errors but no full covariance for this quantity")
    : covariance->second);
  object->Write();
}
}  // namespace fa::derived_science
