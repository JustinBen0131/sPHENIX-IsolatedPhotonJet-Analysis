// MbdCaloCoreCut.h: derive or apply the Au+Au MBD-calorimeter core band.
// loadCalibration() applies the frozen sample-specific cubic definition.
// derive() creates a new calibration from an unweighted pre-cut TH2.
// Input: one sample's unit-count TH2, x=MBD charge, y=total calo energy [GeV],
// native 20 x 20 bins, filled BEFORE any MBD-calo cut, AFTER MB classifier,
// |z_vtx|<10 cm, stored centrality [0,80), sample ownership/quality gates,
// and DATA scaled trigger: Photon 12 GeV + MBD NS >= 2, vtx < 150 cm.
// No physics, centrality or vertex weights. Call independently for each family.
// Derive first, freeze the resulting functions, then apply in the event loop.
// Core sigma differs from the original full-population RMS prescription.
// Numerical convergence alone does not qualify model adequacy or physics closure.
#pragma once
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <Math/MinimizerOptions.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <nlohmann/json.hpp>

namespace mbdcalo {
struct Band {
  std::shared_ptr<TF1> mean, width;
  double fitLow, fitHigh, domainHigh;
  bool pass(double Q, double E) const {
    if (!std::isfinite(Q)||!std::isfinite(E)||Q<0||Q>=domainHigh) return false;
    const double x=Q/1000., mu=mean->Eval(x), sigma=width->Eval(x);
    return sigma>0 && E>=mu-3.5*sigma && E<=mu+3.5*sigma;
  }
};
inline void converged(const TFitResultPtr& r) {
  if (int(r)!=0 || !r->IsValid() || r->CovMatrixStatus()!=3)
    throw std::runtime_error("MBD-calo fit did not converge with accurate covariance");
}
inline Band derive(TH2& h, const std::string& name) {
  constexpr int degree=3;
  for (auto* axis : {h.GetXaxis(),h.GetYaxis()}) {
    if (axis->GetXmin()!=0 || axis->GetXmax()!=3000 || axis->GetNbins()!=150)
      throw std::invalid_argument("Expected the validated 150-bin [0,3000] axes");
  }
  std::vector<double> x,mean,width,meanError,widthError,zero;
  double domainHigh=0.;
  for (int j=1;j<=h.GetNbinsX();++j) {
    std::unique_ptr<TH1D> s(h.ProjectionY((name+"_slice_"+std::to_string(j)).c_str(),j,j,"e"));
    s->SetDirectory(nullptr);
    const double count=s->Integral(), peak=s->GetMaximum();
    if (count>0) domainHigh=h.GetXaxis()->GetBinUpEdge(j);
    if (count<1000 || peak<100) continue;
    for (int b=1;b<=s->GetNbinsX();++b) s->SetBinError(b,std::sqrt(s->GetBinContent(b)));
    const int mode=s->GetMaximumBin(); int left=mode,right=mode;
    while (left>1 && s->GetBinContent(left-1)>=peak/2) --left;
    while (right<s->GetNbinsX() && s->GetBinContent(right+1)>=peak/2) ++right;
    double xl=s->GetBinCenter(left), xr=s->GetBinCenter(right);
    if (left>1) xl-=20*(s->GetBinContent(left)-peak/2)/(s->GetBinContent(left)-s->GetBinContent(left-1));
    if (right<s->GetNbinsX()) xr+=20*(s->GetBinContent(right)-peak/2)/(s->GetBinContent(right)-s->GetBinContent(right+1));
    const double seed=std::max(10.,(xr-xl)/2.354820045), center=s->GetBinCenter(mode);
    int lo=std::max(1,int(std::ceil((center-1.3*seed)/20.)));
    int hi=std::min(150,int(std::ceil((center+1.3*seed)/20.)));
    while (hi-lo+1<5) { if (lo>1) --lo; if (hi<150) ++hi; }
    TF1 core((name+"_core").c_str(),"gaus",s->GetXaxis()->GetBinLowEdge(lo),s->GetXaxis()->GetBinUpEdge(hi));
    core.SetParameters(peak,center,seed);
    core.SetParLimits(0,.1*peak,3*peak);
    core.SetParLimits(1,center-seed,center+seed);
    core.SetParLimits(2,std::max(3.,.35*seed),2.5*seed);
    converged(s->Fit(&core,"LISQRN0")); // Bin-integrated Poisson-likelihood core.
    x.push_back(h.GetXaxis()->GetBinCenter(j)/1000.);
    mean.push_back(core.GetParameter(1)); width.push_back(core.GetParameter(2));
    meanError.push_back(core.GetParError(1)); widthError.push_back(core.GetParError(2)); zero.push_back(0.);
  }
  if (x.size()<=unsigned(degree+1)) throw std::runtime_error("Too few qualified core slices");
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
  auto polynomial=[&](const std::vector<double>& y,const std::vector<double>& error,const char* suffix) {
    TGraphErrors g(x.size(),x.data(),y.data(),zero.data(),error.data());
    std::shared_ptr<TF1> f(new TF1((name+suffix).c_str(),("pol"+std::to_string(degree)).c_str(),x.front(),x.back()));
    converged(g.Fit(f.get(),"SQNR0"));  // Automatic weighted linear least-squares seed.
    converged(g.Fit(f.get(),"SQNR0F")); // All coefficients free, Minuit2/Migrad.
    return f;
  };
  auto mu=polynomial(mean,meanError,"_mean");
  auto sigma=polynomial(width,widthError,"_sigma");
  return Band{std::move(mu),std::move(sigma),1000*x.front(),1000*x.back(),domainHigh};
}
// Load once per sample before the event loop. Unknown samples, malformed
// calibrations and non-finite parameters throw; invalid event values fail pass().
// Sample keys: data, phosim (photon+jet), jetsim (inclusive jet).
inline Band loadCalibration(const std::string& path, const std::string& sample) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open MBD-calo calibration: "+path);
  nlohmann::json calibration; input >> calibration;
  if (calibration.at("schema")!="MbdCaloCubicCalibrationV1" ||
      calibration.at("degree")!=3 || calibration.at("multiplier")!=3.5 ||
      calibration.at("charge_scale")!=1000. ||
      calibration.at("valid_charge_domain")!=nlohmann::json::array({0.,2120.}))
    throw std::invalid_argument("Unsupported MBD-calo calibration contract");
  const auto& family=calibration.at("samples").at(sample);
  auto function=[&](const char* key) {
    const auto coefficients=family.at(key).get<std::vector<double>>();
    if (coefficients.size()!=4 || !std::all_of(coefficients.begin(),coefficients.end(),
                                             [](double p){return std::isfinite(p);}))
      throw std::invalid_argument("Expected four finite cubic coefficients");
    auto f=std::make_shared<TF1>(("mbdcalo_"+sample+"_"+key).c_str(),"pol3",0.,2.12);
    f->SetParameters(coefficients.data());
    return f;
  };
  const auto support=family.at("fit_charge_support").get<std::vector<double>>();
  if (support.size()!=2 || !std::isfinite(support[0]) || !std::isfinite(support[1]) ||
      support[0]<0 || support[0]>=support[1] || support[1]>=2120)
    throw std::invalid_argument("Invalid fitted charge support");
  return Band{function("mean"),function("sigma"),support[0],support[1],2120.};
}
} // namespace mbdcalo
