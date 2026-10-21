// MbdCaloCoreCut.h: derive or apply the Au+Au MBD-calorimeter core band.
// Self-contained: loadCalibration("data"), "phosim" or "jetsim" loads the
// frozen coefficients below. No calibration file or JSON library is needed.
// Load once, then require cut.pass(mbdChargeSum, totalCaloEnergyGeV) per event.
// derive() creates a new calibration from an unweighted pre-cut TH2.
// Input: one sample's unit-count TH2, x=MBD charge, y=total calo energy [GeV],
// native 20 x 20 bins, filled BEFORE any MBD-calo cut, AFTER MB classifier,
// |z_vtx|<10 cm, stored centrality [0,80), sample ownership/quality gates,
// and DATA scaled trigger: Photon 12 GeV + MBD NS >= 2, vtx < 150 cm.
// No physics, centrality or vertex weights. Call independently for each family.
// derive() is for a new calibration; loadCalibration() reproduces the frozen cut.
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
#include <array>

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
// Frozen mbd_calo_cubic_20261020: Q/1000 cubics, mean +/- 3.5 core sigma.
// Load once before the event loop; each family uses its own fitted coefficients.
// Unknown sample keys throw. pass() rejects non-finite or out-of-domain events.
inline Band loadCalibration(const std::string& sample) {
  std::array<double,4> mean, sigma;
  double fitHigh;
  if (sample=="data") {
    mean = {{37.74396713, 422.12296297, 509.88784983, -147.54815561}};
    sigma = {{17.49664184, 58.92212451, 56.44770643, -31.74978374}};
    fitHigh = 2110.0;
  }
  else if (sample=="phosim") {
    mean = {{30.07645580, 591.17568846, 328.09305475, -78.92053693}};
    sigma = {{18.33594272, 93.79209401, -7.17806247, -7.86362368}};
    fitHigh = 1970.0;
  }
  else if (sample=="jetsim") {
    mean = {{41.95544342, 597.35734209, 319.05521941, -75.23903719}};
    sigma = {{25.91486124, 74.04522062, 12.02248690, -13.63297638}};
    fitHigh = 1970.0;
  }
  else throw std::invalid_argument("MBD-calo sample must be data, phosim or jetsim");
  auto function=[&](const std::array<double,4>& coefficients,const char* suffix) {
    auto f=std::make_shared<TF1>(("mbdcalo_"+sample+suffix).c_str(),"pol3",0.,2.12);
    f->SetParameters(coefficients.data());
    return f;
  };
  return Band{function(mean,"_mean"),function(sigma,"_sigma"),50.,fitHigh,2120.};
}
} // namespace mbdcalo
