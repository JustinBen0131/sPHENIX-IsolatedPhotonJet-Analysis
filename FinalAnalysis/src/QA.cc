#include "NativeProducerEvidence.h"
#include "RootScientificTrace.h"
/* Analysis QA and correction diagnostics.

   Reads source-aware histogram packages, applies their normalization once
   and saves candidate purity, leakage, stitching, working-point, trigger
   and centrality checks for the maintained IAN renderer.
   Purity study values share the measurement's correction implementation.
   Saved model diagnostics can also be imported independently of the histogram
   packages, using the same release checks and drawing code.
   Missing or empty populations stay unavailable; figures invent no values.
   The inventory below documents the quantities written to qa.root. */
#include "FinalAnalysis.h"
#include "HistogramErrors.h"
#include "Phases.h"
#include "RenderMetadata.h"
#include "../../TreeToHists/src/Configuration.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>

#include <nlohmann/json.hpp>

#include "TDecompSVD.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

/* QA: the inputs of the photon and event QA figures, written once under qa/<package>/ of the QA file (Output.QA). */

/* final_analysis --qa config/analysis.yaml */

/* Per package (PP, AuAu): every histogram under qa/ of the DATA package (data/), and of the simulated samples summed
   with their cross section per generated event, per family (phosim/: the photon samples of the results, jetsim/: the
   QA.<package>.Inclusive samples); the generator and region-A spectra also per sample (samples/<family>/<sample>/) for
   the stitching figures. And the numerics the figures show, computed here and never in drawing:
   abcd/<class>           the DATA event-leading candidates per region (A-D) on the ABCD photon bins
   purity/<class>/        DATA raw and leakage-corrected ABCD purity (the simulated leakage), the leakage fractions
   closure/<class>/       the purity closure on simulation: raw, corrected and truth purity of the simulated
   pseudo-data (kSignalLabel and kBackgroundLabel below: prompt photons from the photon+jet
   samples, the other candidates from the inclusive-jet samples) and ratio, corrected over
   truth; corrected_r_bkg and ratio_r_bkg, the same with the simulated background correlation
   of the result's photon bins in the estimate; each with its statistical error
   purity_study/<result>/ on the result's photon bins, from the same pseudo-data: r_bkg, the
   background correlation bA bD / (bB bC) of its background candidates;
   purity_scale, the purity non-closure (truth over leakage-corrected ABCD purity on 2 GeV
   bins, nonclosure_ratio, fitted with a line in photon pT weighted by the region-A candidates,
   nonclosure_line, and averaged over each photon bin with the region-A spectrum); the transfer
   terms of the xJ-level ABCD correction per photon and unit xJ (transfer_background: (bA/sA)
   (gA - gC)/(1 - k lC), transfer_leakage: k lC (fA - fC)/(1 - k lC), with f and g the recoil
   jets per signal and background candidate); factorisation_data and _simulation, R of the
   background in the score sub-bands 0.1-0.4, 0.4-0.6 and isolation sub-bands 6-12, 12-40 GeV.
   The analysis takes r_bkg and purity_scale for BackgroundCorrelation and PurityScale "qa".
   working_points/<what>/ the photon-ID score and isolation thresholds at 60, 70, 80, 90% signal efficiency per 5% slice
   and photon-pT group, their per-slice medians and a line in centrality (checks of the packages)
   stitching/<family>_<kind>/  the stitched spectrum, its descriptive fit exp[a - (b + c ln x + d x) ln x] and RMS
   triggers/              the maximum-cluster-energy spectra and turn-ons of the runs with every study trigger enabled
   centrality_target      Au+Au: the DATA centrality of the package's trigger (the centrality weights' target)
   A package without a histogram, or with an empty one, leaves its figure and numbers out; nothing is invented. */

namespace fa {
namespace {

using Json = nlohmann::json;
namespace science=photonjet_science;
using science::Trace;

template<class T> T *traced(T *object,const Trace &trace,const std::string &covarianceState="not_computed",
                           const std::string &method="producer stores per-bin errors; no full covariance") {
  science::attachTrace(object,trace);
  science::setTag(object,"physics_covariance",science::covariance(covarianceState,{},0,method));
  return object;
}
Trace grouped(const TH1 *h,const Edges &edges) {
  if (!h) return Trace::unavailable("missing QA input histogram");
  return science::projectionOf(h,{{{"axis",1},{"low",edges.front()},{"high",edges.back()}}})
    .then("regroup_by_native_lower_edge",{{"edges",edges}});
}
void scalarProof(TDirectory *directory,const std::string &name,double value,const Trace &trace,
                 const std::string &method) {
  const Json record={{"schema","NativeScientificScalarV1"},{"trace",trace.json()},
    {"measurement",std::isfinite(value) ? science::scalar(value,"1",method) : Json{{"state","undefined"},{"reason","no finite supporting points"}}},
    {"covariance",science::covariance("not_computed",{},0,method)}};
  directory->WriteTObject(new TObjString(record.dump().c_str()),("physics_"+name).c_str());
}

std::unique_ptr<TFile> openIfPresent(const std::string &path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) return nullptr;
  observeNativeInput(*f);
  return f;
}

double weightOf(TFile &f) {
  auto *w = dynamic_cast<TParameter<double> *>(f.Get("meta/pb_per_generated_event"));
  if (!w || !(w->GetVal() > 0)) throw std::runtime_error(std::string(f.GetName()) + ": no positive meta/pb_per_generated_event");
  return w->GetVal();
}

/* Sums every histogram under a directory into sums[path], scaled by w (Sumw2: w^2). */
void sumTree(TDirectory *dir, const std::string &prefix, double w, std::map<std::string, std::unique_ptr<TObject>> &sums) {
  TIter next(dir->GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    const std::string path = prefix + key->GetName();
    TObject *o = key->ReadObj();
    if (auto *d = dynamic_cast<TDirectory *>(o)) { sumTree(d, path + "/", w, sums); continue; }
    auto &sum = sums[path];
    if (auto *h = dynamic_cast<TH1 *>(o)) {
      if (!sum) { auto *c = static_cast<TH1 *>(h->Clone()); c->SetDirectory(nullptr); c->Scale(w);
        photonjet_science::sumDefinition(c,h,w,true); sum.reset(c); }
      else {
        auto* target=static_cast<TH1 *>(sum.get());
        if (!target->Add(h, w)) throw std::runtime_error("QA histograms disagree in binning: " + path);
        photonjet_science::sumDefinition(target,h,w,false);
      }
    } else if (auto *s = dynamic_cast<THnSparse *>(o)) {
      if (!sum) { auto *c = static_cast<THnSparse *>(s->Clone()); c->Scale(w); sum.reset(c); }
      else static_cast<THnSparse *>(sum.get())->Add(s, w);
    }
    delete o;
  }
}

void writeTree(TDirectory *out, const std::map<std::string, std::unique_ptr<TObject>> &sums) {
  for (const auto &[path, obj] : sums) {
    if (!obj) continue;
    const auto slash = path.rfind('/');
    TDirectory *d = slash == std::string::npos ? out : out->GetDirectory(path.substr(0, slash).c_str());
    if (!d) out->mkdir(path.substr(0, slash).c_str(), "", true), d = out->GetDirectory(path.substr(0, slash).c_str());
    d->WriteTObject(obj.get(), path.substr(slash + 1).c_str());
  }
}

/* Native bins summed into target bins by their lower edge (the target edges are native edges). */
std::vector<double> regroup(const TH1 *h, int row, const Edges &target, bool variance = false) {
  const TAxis *a = h->GetYaxis();
  std::vector<double> out(target.size() - 1);
  for (int b = 1; b <= a->GetNbins(); ++b) {
    const double low = a->GetBinLowEdge(b) + 1e-9;
    const int i = std::upper_bound(target.begin(), target.end(), low) - target.begin() - 1;
    if (i < 0 || i >= static_cast<int>(out.size())) continue;
    const double e = h->GetBinError(row + 1, b);
    out[i] += variance ? e * e : h->GetBinContent(row + 1, b);
  }
  return out;
}

/* The directory at a path below parent, created when needed (TDirectory::mkdir returns the top of a new hierarchy). */
TDirectory *dirOf(TDirectory *parent, const std::string &path) {
  if (!parent->GetDirectory(path.c_str())) parent->mkdir(path.c_str(), "", true);
  return parent->GetDirectory(path.c_str());
}

void put(TDirectory *d, TObject *o, const std::string &name = "") {  /* written once, then released */
  d->WriteTObject(o, name.empty() ? o->GetName() : name.c_str());
  delete o;
}

/* Null variance denotes a table. Physics quantities pass measured or propagated
   variances; an unavailable variance is NaN rather than a fabricated error. */
TH1D *hist(const std::string &name, const Edges &edges, const std::vector<double> &v,
           const std::vector<double>* variance = nullptr) {
  if (v.size() + 1 != edges.size() || (variance && variance->size() != v.size()))
    throw std::runtime_error("QA histogram values/variances disagree with their edges: " + name);
  auto *h = new TH1D(name.c_str(), "", edges.size() - 1, edges.data());
  h->SetDirectory(nullptr); h->Sumw2();
  for (size_t i = 0; i < v.size(); ++i) {
    h->SetBinContent(i + 1, v[i]);
    h->SetBinError(i + 1, variance ? std::sqrt((*variance)[i]) : 0.0);
  }
  return h;
}

/* Raw (no leakage) and leakage-corrected ABCD purity per bin; NaN where undefined. */
std::pair<std::vector<double>, std::vector<double>> purityRows(const std::vector<std::vector<double>> &counts,
                                                               const std::vector<std::vector<double>> &leakage) {
  const size_t n = counts[0].size();
  std::vector<double> raw(n, std::nan("")), corrected(n, std::nan(""));
  for (size_t k = 0; k < n; ++k) {
    double c[4], l[4];
    for (int r = 0; r < 4; ++r) c[r] = std::max(counts[r][k], 0.0), l[r] = leakage[r][k];
    const double none[4] = {1, 0, 0, 0};
    if (c[0] > 0) raw[k] = solveLeakageABCD(c, none).purity;
    if (c[0] > 0 && l[0] > 0) corrected[k] = solveLeakageABCD(c, l).purity;
  }
  return {raw, corrected};
}

/* numpy.interp(x, xp, fp) for a non-decreasing xp. */
double interp(double x, const std::vector<double> &xp, const std::vector<double> &fp) {
  if (x <= xp.front()) return fp.front();
  if (x >= xp.back()) return fp.back();
  size_t j = std::upper_bound(xp.begin(), xp.end(), x) - xp.begin() - 1;
  return fp[j] + (fp[j + 1] - fp[j]) * (x - xp[j]) / (xp[j + 1] - xp[j]);
}

/* The value at which the requested fraction of a distribution passes (above: from the top, the photon-ID score). */
double threshold(const std::vector<double> &v, const std::vector<double> &edges, double efficiency, bool above) {
  double total = 0;
  for (double x : v) total += x;
  if (!(total > 0)) return std::nan("");
  std::vector<double> cdf(v.size() + 1);
  for (size_t i = 0; i < v.size(); ++i) cdf[i + 1] = cdf[i] + v[i] / total;
  return interp(above ? 1 - efficiency : efficiency, cdf, edges);
}

double median(std::vector<double> v) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  if (v.empty()) return std::nan("");
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* Least squares through the SVD (numpy.linalg.lstsq): the coefficients minimising |A x - b|. */
TVectorD leastSquares(const TMatrixD &a, const TVectorD &b) {
  TDecompSVD svd(a);
  bool ok = false;
  const TVectorD x = svd.Solve(b, ok);  /* the first ncols elements: the least-squares solution */
  if (!ok) return TVectorD(a.GetNcols());
  return x;
}

/* The photon-ID and isolation working points of a signal histogram (value x photon pT x 5% slice). */
void workingPoints(TDirectory *out, const TH3 *h, bool above) {
  const std::vector<std::pair<double, double>> groups = {{15, 17}, {17, 19}, {19, 21}, {21, 23}, {23, 26}, {26, 30}, {30, 35}};
  std::vector<double> edges;
  for (int b = 1; b <= h->GetNbinsX() + 1; ++b) edges.push_back(h->GetXaxis()->GetBinLowEdge(b));
  const double efficiencies[] = {0.6, 0.7, 0.8, 0.9};
  for (double eff : efficiencies) {
    Edges slices(17), groupEdges;
    for (int s = 0; s <= 16; ++s) slices[s] = s;
    for (const auto &g : groups) groupEdges.push_back(g.first);
    groupEdges.push_back(35);
    auto *grid = new TH2D(Form("thresholds_%d", int(std::round(100 * eff))), "", 16, slices.data(), groupEdges.size() - 1, groupEdges.data());
    std::vector<double> perSlice(16);
    for (int s = 0; s < 16; ++s) {
      std::vector<double> ys;
      for (size_t g = 0; g < groups.size(); ++g) {
        std::vector<double> v(edges.size() - 1);
        for (int bx = 1; bx <= h->GetNbinsX(); ++bx)
          for (int by = 1; by <= h->GetNbinsY(); ++by) {
            const double lo = h->GetYaxis()->GetBinLowEdge(by), hi = h->GetYaxis()->GetBinUpEdge(by);
            if (lo < groups[g].first - 1e-9 || hi > groups[g].second + 1e-9) continue;
            for (int bz = 1; bz <= h->GetNbinsZ(); ++bz)
              if (h->GetZaxis()->GetBinLowEdge(bz) >= s - 1e-9 && h->GetZaxis()->GetBinUpEdge(bz) <= s + 1 + 1e-9) v[bx - 1] += h->GetBinContent(bx, by, bz);
          }
        const double t = threshold(v, edges, eff, above);
        grid->SetBinContent(s + 1, g + 1, t);
        grid->SetBinError(s + 1, g + 1, 0.0); // Working-point table, not an event-count histogram.
        if (std::isfinite(t)) ys.push_back(t);
      }
      perSlice[s] = median(ys);
    }
    Edges sliceEdges(17);
    for (int s = 0; s <= 16; ++s) sliceEdges[s] = 5 * s;
    const auto gridTrace=grouped(h,groupEdges).then("interpolated_signal_efficiency_threshold",
      {{"efficiency",eff},{"above",above},{"group_edges",groupEdges},{"slice_edges",slices}});
    const auto medianTrace=gridTrace.then("median_of_finite_group_thresholds");
    put(out, traced(hist(Form("median_%d", int(std::round(100 * eff))), sliceEdges, perSlice),medianTrace));
    grid->SetDirectory(nullptr);
    put(out, traced(grid,gridTrace));
    /* the line in centrality over the slice centres (numpy.polyfit, degree 1) */
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    for (int s = 0; s < 16; ++s)
      if (std::isfinite(perSlice[s])) { const double x = 5 * s + 2.5; sx += x, sy += perSlice[s], sxx += x * x, sxy += x * perSlice[s], ++n; }
    const double slope = n >= 2 ? (n * sxy - sx * sy) / (n * sxx - sx * sx) : std::nan(""), intercept = n >= 2 ? (sy - slope * sx) / n : std::nan("");
    Edges two = {0, 1, 2};
    put(out, traced(hist(Form("line_%d", int(std::round(100 * eff))), two, {intercept, slope}),
      medianTrace.then("ordinary_least_squares_line_finite_slice_centres",{{"centres_percent",{2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5}}})));
  }
}

/* The stitched spectrum of a family: per sample and summed, with the descriptive fit on [lo, hi]. A sample stitched by
   ownership (windows: its [low, high) of the spectrum's variable) enters the sum only inside its window; its own
   spectrum keeps its whole range, so neighbouring samples overlap at every boundary. When the spectrum's variable is
   not the one the sample is stitched on, a window cannot be cut on the axis: ownedParts then gives each sample's
   spectrum of the events it owns, and those are summed. */
using Windows = std::map<std::string, std::pair<double, double>>;
void stitching(TDirectory *out, const std::map<std::string, TH1D *> &parts, double lo, double hi, const Windows &windows = {}, bool generator = true, bool relativeErrors = true,
               const std::map<std::string, TH1D *> &ownedParts = {}) {
  if (parts.empty()) return;
  TH1D *total = nullptr;
  std::vector<Trace> summedTraces;
  for (const auto &[sample, h] : parts) {
    out->WriteTObject(h, ("sample_" + sample).c_str());
    std::unique_ptr<TH1D> owned(static_cast<TH1D *>((ownedParts.count(sample) ? ownedParts.at(sample) : h)->Clone()));
    owned->SetDirectory(nullptr);
    auto ownedTrace=science::traceOf(owned.get());
    if (!ownedParts.count(sample) && windows.count(sample))
      ownedTrace=ownedTrace.then("sample_ownership_window",{{"low",windows.at(sample).first},
        {"high",std::isfinite(windows.at(sample).second) ? Json(windows.at(sample).second) : Json("unbounded")},{"flow","excluded"}});
    summedTraces.push_back(ownedTrace);
    if (!ownedParts.count(sample) && windows.count(sample))
      for (int b = 0; b <= owned->GetNbinsX() + 1; ++b) {
        const double a = owned->GetXaxis()->GetBinLowEdge(b), c = owned->GetXaxis()->GetBinUpEdge(b);
        if (b == 0 || b > owned->GetNbinsX() || a < windows.at(sample).first - 1e-9 || c > windows.at(sample).second + 1e-9)
          owned->SetBinContent(b, 0), owned->SetBinError(b, 0);
      }
    if (!total) { total = static_cast<TH1D *>(owned->Clone("total")); total->SetDirectory(nullptr); }
    else total->Add(owned.get());
  }
  std::vector<double> x, y, e;
  const auto totalTrace=Trace::combine("sum_owned_samples",summedTraces);
  traced(total,totalTrace);
  for (int b = 1; b <= total->GetNbinsX(); ++b) {
    const double a = total->GetXaxis()->GetBinLowEdge(b), c = total->GetXaxis()->GetBinUpEdge(b), w = c - a;
    if (a < lo - 1e-9 || c > hi + 1e-9) continue;
    const double yy = total->GetBinContent(b) / w, ee = total->GetBinError(b) / w;
    if (yy > 0 && std::isfinite(yy) && ee > 0) x.push_back(0.5 * (a + c)), y.push_back(yy), e.push_back(ee);
  }
  std::vector<double> coefficients(4, std::nan("")), coefficientVariance(4, std::nan(""));
  if (x.size() >= 5) {
    TMatrixD a(x.size(), generator ? 4 : 3);
    TVectorD rhs(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      const double lx = std::log(x[i] / (generator ? 1.0 : 20.0)), s = relativeErrors ? y[i] / e[i] : 1.0;
      a(i, 0) = s;
      a(i, 1) = (generator ? -lx : lx) * s;
      a(i, 2) = (generator ? -lx * lx : lx * lx) * s;
      if (generator) a(i, 3) = -x[i] * lx * s;
      rhs[i] = std::log(y[i]) * s;
    }
    const TVectorD c = leastSquares(a, rhs);
    for (int k = 0; k < (generator ? 4 : 3); ++k) coefficients[k] = c[k];
    // Conditional fit covariance: the descriptive fit's weights are fixed.
    // Propagate Var(log y)=(error/y)^2 through the same SVD pseudoinverse.
    TDecompSVD covarianceSvd(a);
    bool invertible = false;
    const TMatrixD inverse = covarianceSvd.Invert(invertible);
    if (invertible) for (int k = 0; k < a.GetNcols(); ++k) {
      coefficientVariance[k] = 0;
      for (size_t i = 0; i < x.size(); ++i) {
        const double scaledError = (relativeErrors ? y[i] / e[i] : 1.0) * e[i] / y[i];
        coefficientVariance[k] += std::pow(inverse(k, i) * scaledError, 2);
      }
    }
    if (!generator) coefficients[3] = 0, coefficientVariance[3] = 0;
  }
  /* the RMS of spectrum / fit - 1 over the fitted points */
  double squares = 0;
  int n = 0;
  for (size_t i = 0; i < x.size() && std::isfinite(coefficients[0]); ++i) {
    const double lx = std::log(x[i] / (generator ? 1.0 : 20.0));
    const double exponent = generator
      ? coefficients[0] - (coefficients[1] + coefficients[2] * lx + coefficients[3] * x[i]) * lx
      : coefficients[0] + coefficients[1] * lx + coefficients[2] * lx * lx;
    const double r = y[i] / std::exp(exponent);
    if (std::isfinite(r)) squares += (r - 1) * (r - 1), ++n;
  }
  Edges four = {0, 1, 2, 3, 4};
  const auto fitTrace=totalTrace.then("svd_descriptive_log_fit",{{"form",generator ? "modified_power_law" : "log_parabola_et_over_20"},
    {"low",lo},{"high",hi},{"relative_errors",relativeErrors},{"bin_values","content_per_width"},{"fit_points","positive_finite_content_and_positive_error"}});
  put(out, traced(hist("fit", four, coefficients, &coefficientVariance),fitTrace));
  out->WriteTObject(new TObjString(generator ? "modified_power_law" : "log_parabola_et_over_20"), "fit_kind");
  put(out, new TParameter<double>("rms", n ? std::sqrt(squares / n) : std::nan("")));
  scalarProof(out,"rms",n ? std::sqrt(squares/n) : std::nan(""),fitTrace.then("rms_spectrum_over_fit_minus_one"),"RMS over finite fitted points; uncertainty not computed");
  put(out, total);
}

using Sums = std::map<std::string, std::unique_ptr<TObject>>;
/* Sparse ROOT histograms cannot retain function-list tags. Keep their source
   trace, including unavailable provenance, in a neighboring record. */
void bindSampleSelection(Sums& sums,TDirectory* output,const std::function<Trace(Trace)>& selection) {
  for(auto& [path,object]:sums) {
    if(!object)continue; // Non-histogram keys have no accumulated object.
    const auto trace=selection(science::traceOf(object.get()));
    if(science::tags(object.get())) {
      science::attachTrace(object.get(),trace);
    } else if(dynamic_cast<THnSparse*>(object.get())) {
      const auto slash=path.rfind('/');
      auto* directory=slash==std::string::npos?output:dirOf(output,path.substr(0,slash));
      const auto name=slash==std::string::npos?path:path.substr(slash+1);
      const Json record={{"schema","NativeScientificObjectV1"},{"root_class",object->ClassName()},
        {"trace",trace.json()},{"covariance",science::covariance("not_computed",{},0,"sparse producer retains bin errors; full covariance unavailable")}};
      TObjString proof(record.dump().c_str());
      directory->WriteTObject(&proof,("physics_"+name).c_str());
    } else throw std::runtime_error("unsupported accumulated QA object: "+path+" ("+object->ClassName()+")");
  }
}
const Edges kFineEdges = {15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35};  /* the purity non-closure's 2 GeV bins */

TH1 *sumOf(Sums &sums, const std::string &path) {
  auto it = sums.find(path);
  return it == sums.end() ? nullptr : static_cast<TH1 *>(it->second.get());
}

/* A histogram with any content or error inside its axes; an empty one is absent support. */
bool filled(const TH1 *h) {
  for (int b = 0; h && b < h->GetNcells(); ++b)
    if (!h->IsBinUnderflow(b) && !h->IsBinOverflow(b) && (h->GetBinContent(b) != 0 || h->GetBinError(b) != 0)) return true;
  return false;
}

/* Region (A-D) x bin counts of the candidates with the given family labels on the given ET bins; complete turns false
   when a label has no fill. */
using Labels = std::vector<std::pair<Sums *, std::string>>;
Trace leadersTrace(const Labels &labels,const std::string &centrality,const Edges &edges,const std::string &observable="abcd") {
  std::vector<Trace> inputs;
  for (const auto &label:labels) inputs.push_back(grouped(sumOf(*label.first,"photons/"+centrality+"/"+label.second+"/"+observable),edges));
  return Trace::combine("sum_source_role_labels",inputs);
}

/* The pseudo-data of every simulation-only purity quantity: the analysis-signal candidates of the photon+jet samples
   and the other candidates of the inclusive-jet samples. Every simulated sample is generated with both HardQCD:all and
   PromptPhoton:all and differs only in its truth filter (a photon or a jet above threshold; PPG12 note, Monte Carlo
   samples), so the inclusive-jet samples hold the photon+jet samples' prompt photons again (in the delivered packages
   their analysis-signal candidates are 0.94 +- 0.04 of the photon samples' in every ABCD region and photon E_T bin,
   p+p and Au+Au alike) and the photon+jet samples' other candidates are inclusive-jet candidates too. Summing all four
   labels counts both twice. This is the source-role rule of the photon-ID training: prompt from photon+jet, nonprompt
   from inclusive+jet, the cross-role candidates (jetprompt, phoother) left out. Both stay filled as their own labels. */
const char *const kSignalLabel = "sig";       /* of the photon-family simulation */
const char *const kBackgroundLabel = "bkg";   /* of the inclusive-family simulation */
std::vector<std::vector<double>> leaders(const Labels &labels, const std::string &c, const Edges &edges, bool &complete,
                                         bool variance = false) {
  std::vector<std::vector<double>> out(4, std::vector<double>(edges.size() - 1));
  for (const auto &label : labels) {
    auto *h = dynamic_cast<TH2 *>(sumOf(*label.first, "photons/" + c + "/" + label.second + "/abcd"));
    if (!filled(h)) { complete = false; continue; }
    for (int r = 0; r < 4; ++r) {
      const auto v = regroup(h, r, edges, variance);
      for (size_t i = 0; i < v.size(); ++i) out[r][i] += v[i];
    }
  }
  return out;
}

/* Recoil jets of the region-A (region 0) or region-C (region 1) candidates with the given labels: photon bin x xJ. */
std::vector<std::vector<double>> recoilJets(const Labels &labels, const std::string &c, int region, const Edges &edges, Edges &xj, bool &complete) {
  std::vector<std::vector<double>> out(edges.size() - 1);
  for (const auto &label : labels) {
    auto *h = dynamic_cast<TH3 *>(sumOf(*label.first, "photons/" + c + "/" + label.second + "/recoil"));
    if (!filled(h)) { complete = false; continue; }
    xj.clear();
    for (int z = 1; z <= h->GetNbinsZ() + 1; ++z) xj.push_back(h->GetZaxis()->GetBinLowEdge(z));
    for (auto &row : out) row.resize(h->GetNbinsZ());
    for (int y = 1; y <= h->GetNbinsY(); ++y) {
      const int i = std::upper_bound(edges.begin(), edges.end(), h->GetYaxis()->GetBinLowEdge(y) + 1e-9) - edges.begin() - 1;
      if (i < 0 || i >= static_cast<int>(out.size())) continue;
      for (int z = 1; z <= h->GetNbinsZ(); ++z) out[i][z - 1] += h->GetBinContent(region + 1, y, z);
    }
  }
  return out;
}

/* R of the background in two score and two isolation sub-bands (every ET): a d / (b c) and its statistical error, with
   a, b, c, d laid out as the regions A-D (a: the higher scores and the lower isolation energies), so it reads as R_bkg
   does. The error is that of the weighted contents (for DATA the Poisson error). */
TH1D *factorisation(const std::string &name, const std::vector<TH1 *> &parts) {
  double n[4] = {0, 0, 0, 0}, e2[4] = {0, 0, 0, 0};
  bool any = false;
  std::vector<Trace> ancestry;
  for (TH1 *part : parts) {
    auto *h = dynamic_cast<TH3 *>(part);
    if (!h) continue;
    ancestry.push_back(science::projectionOf(h,{{{"axis",0},{"low",0.1},{"high",0.6}},{{"axis",1},{"low",6.0},{"high",40.0}}}));
    any = true;
    for (int x = 1; x <= h->GetNbinsX(); ++x)
      for (int y = 1; y <= h->GetNbinsY(); ++y) {
        const double s0 = h->GetXaxis()->GetBinLowEdge(x), s1 = h->GetXaxis()->GetBinUpEdge(x);
        const double i0 = h->GetYaxis()->GetBinLowEdge(y), i1 = h->GetYaxis()->GetBinUpEdge(y);
        const int score = s0 >= 0.4 - 1e-9 && s1 <= 0.6 + 1e-9 ? 0 : s0 >= 0.1 - 1e-9 && s1 <= 0.4 + 1e-9 ? 1 : -1;
        const int iso = i0 >= 6 - 1e-9 && i1 <= 12 + 1e-9 ? 0 : i0 >= 12 - 1e-9 && i1 <= 40 + 1e-9 ? 1 : -1;
        if (score < 0 || iso < 0) continue;
        for (int z = 1; z <= h->GetNbinsZ(); ++z) {
          n[2 * score + iso] += h->GetBinContent(x, y, z);
          e2[2 * score + iso] += h->GetBinError(x, y, z) * h->GetBinError(x, y, z);
        }
      }
  }
  if (!any) return nullptr;
  const double r = n[1] > 0 && n[2] > 0 ? n[0] * n[3] / (n[1] * n[2]) : std::nan("");
  double relative2 = 0;
  for (int q = 0; q < 4; ++q) relative2 += n[q] > 0 ? e2[q] / (n[q] * n[q]) : 0;
  const double error = std::min({n[0], n[1], n[2], n[3]}) > 0 ? r * std::sqrt(relative2) : std::nan("");
  auto *h = new TH1D(name.c_str(), "score 0.4-0.6 and 0.1-0.4 x isolation 6-12 and 12-40 GeV", 6, 0, 6);
  h->SetDirectory(nullptr);
  const char *labels[] = {"a", "b", "c", "d", "R", "R_error"};
  const double values[] = {n[0], n[1], n[2], n[3], r, error};
  h->Sumw2();
  for (int b = 0; b < 6; ++b) {
    h->GetXaxis()->SetBinLabel(b + 1, labels[b]); h->SetBinContent(b + 1, values[b]);
    h->SetBinError(b + 1, b < 4 ? std::sqrt(e2[b]) : b == 4 ? error : 0.0);
  }
  return traced(h,Trace::combine("factorisation_abcd_ratio",ancestry,{{"score_bands",{0.1,0.4,0.6}},
    {"isolation_bands",{6.0,12.0,40.0}},{"variance","independent_weighted_cells_delta_method"}}));
}

/* The purity study of one result (see the header). Each quantity needs its candidate labels filled in both simulated
   families (the non-closure also the response's leakage); one that lacks them is left out. */
/* Association QA uses the producer's global-leading jet: choose the highest-pT jet after threshold and photon
   overlap, then apply the fiducial/recoil veto without promoting a softer jet. The dense aliases share that same
   stored definition; the inclusive per-photon measurement is unchanged. */
void leadingData(TFile &f, const Settings &s, const SystemSettings &sys, TDirectory *out) {
  const std::string path="measurement/"+sys.centrality+"/"+s.jet+"/leading_global";
  auto *source=dynamic_cast<TH2*>(f.Get(path.c_str()));
  if(!source)throw std::runtime_error(std::string(f.GetName())+": association QA requires "+path);
  out->WriteTObject(new TObjString("global_leading_v1"),"leading_definition");
  auto *leading=static_cast<TH2*>(source->Clone("leading_data"));leading->SetDirectory(nullptr);put(out,leading);
  if(sys.package=="PP"){
    auto *dense=static_cast<TH2*>(source->Clone("leading_data_u"));dense->SetDirectory(nullptr);put(out,dense);
  }
}

void purityStudy(TDirectory *out, const SystemSettings &sys, const TH2 *leakage, Sums &phosim, Sums &jetsim, Sums &data, bool diagnostics=true) {
  const std::string &c = sys.centrality;
  const Edges &pt = sys.photonEdges;
  const Labels signal = {{&phosim, kSignalLabel}}, background = {{&jetsim, kBackgroundLabel}};
  bool withBackground = true, withSignal = true;
  const auto b = leaders(background, c, pt, withBackground), s = leaders(signal, c, pt, withSignal);
  if (withBackground) {
    const auto bVariance = leaders(background, c, pt, withBackground, true);
    std::vector<double> r(pt.size() - 1, std::nan("")), rVariance = r;
    for (size_t i = 0; i < r.size(); ++i) {
      const auto estimate = abcdRatio({b[0][i], b[1][i], b[2][i], b[3][i]},
                                      {bVariance[0][i], bVariance[1][i], bVariance[2][i], bVariance[3][i]});
      // Preserve the established availability domain of the central value.
      if (b[0][i] > 0 && b[1][i] > 0 && b[2][i] > 0 && b[3][i] > 0)
        r[i] = estimate.value, rVariance[i] = estimate.variance;
    }
    put(out, traced(hist("r_bkg", pt, r, &rVariance),leadersTrace(background,c,pt).then("ABCD_background_correlation",
      {{"formula","A*D/(B*C)"},{"availability","all four cells positive"},{"variance","delta_method_independent_weighted_cells"}})));
    put(out, new TObjString("delta_method_independent_weighted_ABCD_cells; variance=sumw2; event_cross_covariance_not_included"), "r_bkg_uncertainty");
  }
  /* Nonclosure uses the same underlying signal in truth purity and leakage.
     Propagate the native weighted cells through the ratio, fitted line and
     spectrum average together. Matching the leakage moments to the projected
     signal moments is required before treating them as the same population.
     This is a cell-Poisson delta model; cross-region event covariance requires
     event records and is not claimed by these QA errors. */
  if (leakage && withBackground && withSignal) {
    auto* signalHistogram = dynamic_cast<TH2*>(sumOf(phosim, "photons/" + c + "/" + kSignalLabel + "/abcd"));
    auto* backgroundHistogram = dynamic_cast<TH2*>(sumOf(jetsim, "photons/" + c + "/" + kBackgroundLabel + "/abcd"));
    using Counts = std::vector<std::vector<double>>;
    auto edgesOf = [](const TH2* h) {
      Edges edges(h->GetNbinsY() + 1);
      for (int y = 0; y <= h->GetNbinsY(); ++y) edges[y] = h->GetYaxis()->GetBinLowEdge(y + 1);
      return edges;
    };
    const Edges native = edgesOf(signalHistogram), backgroundEdges = edgesOf(backgroundHistogram), leakEdges = edgesOf(leakage);
    if (native != backgroundEdges) throw std::runtime_error("nonclosure signal/background native axes differ");
    bool complete = true;
    const Counts signalNative = leaders(signal, c, native, complete), backgroundNative = leaders(background, c, native, complete);
    const Counts signalVariance = leaders(signal, c, native, complete, true), backgroundVariance = leaders(background, c, native, complete, true);
    auto into = [](const Counts& counts, const Edges& from, const Edges& target) {
      Counts rows(4, std::vector<double>(target.size() - 1));
      for (int q = 0; q < 4; ++q) for (size_t j = 0; j + 1 < from.size(); ++j) {
        const int k = std::upper_bound(target.begin(), target.end(), from[j] + 1e-9) - target.begin() - 1;
        if (k >= 0 && k < static_cast<int>(rows[q].size())) rows[q][k] += counts[q][j];
      }
      return rows;
    };
    Counts frozenLeak(4), frozenLeakVariance(4);
    for (int q = 0; q < 4; ++q) {
      frozenLeak[q] = regroup(leakage, q, leakEdges);
      frozenLeakVariance[q] = regroup(leakage, q, leakEdges, true);
    }
    const Counts projectedSignal = into(signalNative, native, leakEdges), projectedVariance = into(signalVariance, native, leakEdges);
    bool sharedLeakage = complete;
    for (int q = 0; q < 4; ++q) for (size_t j = 0; j < frozenLeak[q].size(); ++j) {
      if (leakEdges[j] < kFineEdges.front() || leakEdges[j] >= kFineEdges.back()) continue;
      // QA contents may have single-precision storage; Sumw2 is double.
      const auto close = [](double a, double b, double relative) { return std::abs(a-b) <= 1e-10 + relative * std::max(std::abs(a), std::abs(b)); };
      sharedLeakage = sharedLeakage && close(frozenLeak[q][j], projectedSignal[q][j], 5e-6) &&
                                      close(frozenLeakVariance[q][j], projectedVariance[q][j], 1e-9);
    }
    struct Estimate { std::vector<double> ratio, line, scale; int used = 0; };
    auto estimate = [&](const Counts& sn, const Counts& bn, bool moved) {
      const Counts sf = into(sn, native, kFineEdges), bf = into(bn, native, kFineEdges);
      // Move leakage by the same source-cell displacement while preserving the
      // original rounded central value exactly at the expansion point.
      Counts rawLeak = frozenLeak;
      if (moved) {
        const Counts changed = into(sn, native, leakEdges);
        for (int q=0;q<4;++q) for (size_t j=0;j<rawLeak[q].size();++j) rawLeak[q][j] += changed[q][j]-projectedSignal[q][j];
      }
      const Counts l = into(rawLeak, leakEdges, kFineEdges);
      Counts total(4, std::vector<double>(kFineEdges.size() - 1));
      for (int q=0;q<4;++q) for (size_t i=0;i<total[q].size();++i) total[q][i]=sf[q][i]+bf[q][i];
      const auto corrected=purityRows(total,l).second;
      Estimate value; value.ratio.assign(total[0].size(),std::nan(""));
      double w=0,wx=0,wy=0,wxx=0,wxy=0;
      for (size_t i=0;i<value.ratio.size();++i) {
        const double truth=total[0][i]>0 ? sf[0][i]/total[0][i] : std::nan(""), x=.5*(kFineEdges[i]+kFineEdges[i+1]);
        value.ratio[i]=truth/corrected[i];
        if(!std::isfinite(value.ratio[i]) || !(total[0][i]>0) || !(corrected[i]>.05) || !(truth>0))continue;
        w+=total[0][i];wx+=total[0][i]*x;wy+=total[0][i]*value.ratio[i];wxx+=total[0][i]*x*x;wxy+=total[0][i]*x*value.ratio[i];++value.used;
      }
      if(value.used<2)return value;
      const double slope=(w*wxy-wx*wy)/(w*wxx-wx*wx), intercept=(wy-slope*wx)/w;
      value.line={intercept,slope}; value.scale.resize(pt.size()-1);
      for(size_t k=0;k<value.scale.size();++k) {
        double sum=0,weighted=0;
        for(const Counts* population:{&sn,&bn})for(size_t j=0;j+1<native.size();++j) {
          const double centre=.5*(native[j]+native[j+1]);
          if(centre<pt[k] || centre>=pt[k+1])continue;
          sum+=(*population)[0][j];weighted+=(*population)[0][j]*(intercept+slope*centre);
        }
        value.scale[k]=sum>0 ? weighted/sum : intercept+slope*.5*(pt[k]+pt[k+1]);
      }
      return value;
    };
    const Estimate central=estimate(signalNative,backgroundNative,false);
    Estimate variance=central;
    for(auto* values:{&variance.ratio,&variance.line,&variance.scale})
      values->assign(values->size(),sharedLeakage ? 0.0 : std::nan(""));
    if(sharedLeakage)for(bool ofSignal:{true,false})for(int q=0;q<4;++q)for(size_t j=0;j+1<native.size();++j) {
      const double cellVariance=(ofSignal ? signalVariance : backgroundVariance)[q][j];
      if(!(cellVariance>0))continue;
      const double n=(ofSignal ? signalNative : backgroundNative)[q][j], step=1e-5*std::max(1.0,std::abs(n));
      Estimate sides[2];
      for(int side=0;side<2;++side) {
        Counts sn=signalNative,bn=backgroundNative;
        (ofSignal ? sn : bn)[q][j]+=side ? -step : step;
        sides[side]=estimate(sn,bn,true);
      }
      auto add=[&](std::vector<double>& result,const std::vector<double>& up,const std::vector<double>& down) {
        if(up.size()!=result.size() || down.size()!=result.size()) {std::fill(result.begin(),result.end(),std::nan(""));return;}
        for(size_t k=0;k<result.size();++k) result[k]+=std::pow((up[k]-down[k])/(2*step),2)*cellVariance;
      };
      add(variance.ratio,sides[0].ratio,sides[1].ratio);add(variance.line,sides[0].line,sides[1].line);add(variance.scale,sides[0].scale,sides[1].scale);
    }
    const auto nonclosure=Trace::combine("purity_nonclosure_native_cell_delta_method",
      {leadersTrace(signal,c,kFineEdges),leadersTrace(background,c,kFineEdges),grouped(leakage,kFineEdges)},
      {{"shared_signal_leakage",sharedLeakage},{"ratio","truth_over_leakage_corrected_purity"},
       {"step_relative",1e-5},{"step_floor",1.0},{"event_cross_covariance","not_included"}});
    const auto nonclosureFit=nonclosure.then("region_A_weighted_line_fit",{{"minimum_corrected_purity",0.05},{"finite_positive_truth",true}});
    put(out,traced(hist("nonclosure_ratio",kFineEdges,central.ratio,&variance.ratio),nonclosure,
      sharedLeakage ? "not_computed" : "unknown",sharedLeakage ? "per-bin delta variances; no full covariance" : "leakage and signal moments differ"));
    if(central.used>=2) {
      put(out,traced(hist("nonclosure_line",{0,1,2},central.line,&variance.line),nonclosureFit,
        sharedLeakage ? "not_computed" : "unknown","per-parameter variances only; unavailable when leakage and signal differ"));
      put(out,traced(hist("purity_scale",pt,central.scale,&variance.scale),Trace::combine("region_A_spectrum_average_line",
        {nonclosureFit,leadersTrace(signal,c,pt),leadersTrace(background,c,pt)},{{"photon_edges",pt},{"empty_bin","line_at_bin_centre"}}),
        sharedLeakage ? "not_computed" : "unknown","per-bin variances only; unavailable when leakage and signal differ"));
    }
    put(out,new TObjString(sharedLeakage
        ? "native_weighted_cell_delta_method; shared_signal_leakage; fit_and_spectrum_weights_propagated; event_cross_covariance_not_included"
        : "unavailable: leakage and QA signal moments differ; joint source/event moments are required; NaN errors"),"nonclosure_uncertainty");
  }
  if(!diagnostics)return;
  /* the transfer terms of the xJ-level correction */
  Edges xj;
  bool complete = withBackground && withSignal;
  const auto fA = recoilJets(signal, c, 0, pt, xj, complete), fC = recoilJets(signal, c, 1, pt, xj, complete);
  const auto gA = recoilJets(background, c, 0, pt, xj, complete), gC = recoilJets(background, c, 1, pt, xj, complete);
  if (complete && !xj.empty()) {
    auto *terms = new TH2D("transfer_background", "(bA/sA)(gA - gC)/(1 - k lC) per signal photon and unit xJ", pt.size() - 1, pt.data(), xj.size() - 1, xj.data());
    auto *leak = new TH2D("transfer_leakage", "k lC (fA - fC)/(1 - k lC) per signal photon and unit xJ", pt.size() - 1, pt.data(), xj.size() - 1, xj.data());
    for (TH2D *h : {terms, leak}) {h->SetDirectory(nullptr); h->Sumw2();}
    for (size_t i = 0; i + 1 < pt.size(); ++i) {
      const double k = b[0][i] / b[2][i], lC = s[2][i] / s[0][i], denominator = 1 - k * lC;
      for (size_t x = 0; x + 1 < xj.size(); ++x) {
        const double width = xj[x + 1] - xj[x];
        terms->SetBinContent(i + 1, x + 1, (b[0][i] / s[0][i]) * (gA[i][x] / b[0][i] - gC[i][x] / b[2][i]) / denominator / width);
        leak->SetBinContent(i + 1, x + 1, k * lC * (fA[i][x] / s[0][i] - fC[i][x] / s[2][i]) / denominator / width);
        // These per-photon all-jet ratios need joint event moments, including
        // multiple jets in one bin. Marginal Sumw2 cannot supply that covariance.
        terms->SetBinError(i + 1, x + 1, std::nan(""));
        leak->SetBinError(i + 1, x + 1, std::nan(""));
      }
    }
    const auto transferTrace=Trace::combine("per_signal_photon_recoil_transfer",
      {leadersTrace(signal,c,pt),leadersTrace(background,c,pt),leadersTrace(signal,c,pt,"recoil"),leadersTrace(background,c,pt,"recoil")},
      {{"xj_edges",xj},{"normalization","per_photon_per_xj_width"},{"cross_moments","unavailable"}});
    put(out, traced(terms,transferTrace.then("background_term"),"unknown","joint photon and recoil moments unavailable"));
    put(out, traced(leak,transferTrace.then("leakage_term"),"unknown","joint photon and recoil moments unavailable"));
    put(out, new TObjString("unavailable: per-event joint photon/recoil moments are required; NaN errors must not be used as uncertainties"), "transfer_uncertainty");
  }
  /* the factorisation of the background: DATA and the simulated non-signal candidates */
  if (TH1D *h = factorisation("factorisation_data", {sumOf(data, "photons/" + c + "/data/scoreiso")})) put(out, h);
  if (TH1D *h = factorisation("factorisation_simulation", {sumOf(jetsim, "photons/" + c + "/" + kBackgroundLabel + "/scoreiso")}))
    put(out, h);
}

}  /* namespace */

/* Save QA quantities from the declared packages for later rendering. Missing QA inputs stay absent;
   the renderer reports unavailable figures rather than treating an empty plot as a measured zero. */
const char *const kSourceStudies = "PhotonJetSourceStudiesV2";  /* V2: the photon live bit among reference-scaled events */

void copyDirectory(TDirectory *from, TDirectory *to, const std::set<std::string> &skip, const std::string &prefix) {
  std::set<std::string> seen;
  TIter next(from->GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    const std::string name = key->GetName();
    if (skip.count(prefix + name) || !seen.insert(name).second) continue;  /* the newest cycle of a name is listed first */
    TClass *cls = TClass::GetClass(key->GetClassName());
    TObject *object = key->ReadObj();
    if (!cls || !object) throw std::runtime_error("cannot read " + prefix + name);
    if (cls->InheritsFrom(TDirectory::Class())) {
      copyDirectory(static_cast<TDirectory *>(object), to->mkdir(name.c_str()), skip, prefix + name + "/");
    } else if (cls->InheritsFrom(TTree::Class())) {
      to->cd();
      TTree *copy = static_cast<TTree *>(object)->CloneTree(-1, "fast");
      science::attachTrace(copy,science::traceOf(object));
      copy->Write(name.c_str(), TObject::kOverwrite);
    } else {
      to->WriteTObject(object, name.c_str());
      if (auto *histogram = dynamic_cast<TH1 *>(object)) histogram->SetDirectory(nullptr);
      delete object;
    }
  }
}

void sourceStudies(TFile *data, const std::string &package, TDirectory *out) {
  auto *sources = dynamic_cast<TTree *>(data->Get("records/sources"));
  auto *cells = dynamic_cast<TTree *>(data->Get("records/source_contributions"));
  auto *bitsRaw = dynamic_cast<TH3 *>(data->Get("qa/events/trigger_bits_raw"));
  const bool perRunStudy=data->Get("qa/events/trigger_study_bits_raw")!=nullptr;
  if (perRunStudy) bitsRaw=dynamic_cast<TH3 *>(data->Get("qa/events/trigger_study_bits_raw"));
  /* The study population (|vz| < 30 cm, a cluster of 1 GeV or more) when the package has it. The efficiency is the
     IAN's: among the events recorded by the scaled reference trigger, the share whose photon trigger was live, a
     binomial ratio inside one triggered sample that needs no prescale counter. */
  const std::string studyCells=perRunStudy ? "qa/events/trigger_study_bits_raw" : "qa/events/trigger_bits_raw";
  if (sources && cells && bitsRaw) {
    const bool auau = package == "AuAu";
    const std::string reference = auau ? "MBD N&S >= 2, vtx < 150 cm" : "MBD N&S >= 1";
    const std::vector<std::pair<std::string, std::string>> photonsStudy = auau
        ? std::vector<std::pair<std::string, std::string>>{{"Photon 10", "Photon 10 GeV + MBD NS >= 2"}, {"Photon 12", "Photon 12 GeV + MBD NS >= 2"}}
        : std::vector<std::pair<std::string, std::string>>{{"Photon 3", "Photon 3 GeV + MBD NS >= 1"}, {"Photon 4", "Photon 4 GeV + MBD NS >= 1"}, {"Photon 5", "Photon 5 GeV + MBD NS >= 1"}};
    const int firstRun = auau ? 0 : 47289;
    auto bitOf = [](const Json &names, const std::string &prefix) -> int {
      for (size_t i = 0; i < names.size(); ++i) if (names[i] == prefix) return i;
      int first = -1, preferred = -1;
      for (size_t i = 0; i < names.size(); ++i) {
        const std::string n = names[i];
        if (n.rfind(prefix, 0) != 0) continue;
        if (first < 0) first = i;
        if (preferred < 0 && n.find("vtx < 150 cm") != std::string::npos) preferred = i;
      }
      return preferred >= 0 ? preferred : first;
    };
    const int nE = bitsRaw->GetNbinsX();
    std::vector<double> refSpectrum(nE), den(nE);
    std::vector<double> refVariance(nE), denVariance(nE);
    std::map<std::string, std::vector<double>> photonSpectrum, num;
    std::map<std::string, std::vector<double>> photonVariance, numVariance;
    for (const auto &p : photonsStudy) photonSpectrum[p.first].assign(nE, 0), num[p.first].assign(nE, 0);
    for (const auto &p : photonsStudy) photonVariance[p.first].assign(nE, 0), numVariance[p.first].assign(nE, 0);
    int runs = 0, firstSeen = 0, lastSeen = 0;
    std::set<int> studiedRuns;
    Json triggerBindings=Json::array(),centralityBindings=Json::array();
    TTreeReader rs(sources), rc(cells);
    TTreeReaderValue<std::string> json(rs, "json"), histogram(rc, "histogram");
    TTreeReaderValue<Long64_t> begin(rs, "contribution_begin"), end(rs, "contribution_end"), bin(rc, "bin");
    TTreeReaderValue<double> value(rc, "value"), variance(rc,"variance");
    std::unique_ptr<TH3> one(static_cast<TH3 *>(bitsRaw->Clone("one_source")));
    one->SetDirectory(nullptr);
    std::vector<double> target, targetVariance;
    std::unique_ptr<TH2> centralityOne;
    auto *centralityBits = dynamic_cast<TH2 *>(data->Get("qa/events/centrality_scaled_bits"));
    if (auau && centralityBits) {
      centralityOne.reset(static_cast<TH2 *>(centralityBits->Clone("centrality_one")));
      centralityOne->SetDirectory(nullptr);
      target.assign(centralityBits->GetNbinsX(), 0.0);
      targetVariance.assign(centralityBits->GetNbinsX(), 0.0);
    }
    const Json selection = Json::parse(static_cast<TObjString *>(data->Get("meta/package"))->GetString().Data()).at("selection");
    const std::string trigger = selection.at("systems").at(auau ? "auau" : "pp").at("event_selection").at("trigger_name");
    while (rs.Next()) {
      const Json source = Json::parse(*json);
      const Json &m = source.at("metadata");
      const Json info = m.contains("trigger_run_info") ? m.at("trigger_run_info") : m.value("trigger_selection", Json());
      const int run = m.at("source").value("run", -1);
      if (info.is_null()) continue;
      one->Reset();
      if (centralityOne) centralityOne->Reset();
      for (Long64_t i = *begin; i < *end; ++i) {
        if (rc.SetEntry(i) != TTreeReader::kEntryValid) break;
        if (*histogram == studyCells) {one->SetBinContent(static_cast<int>(*bin), *value);one->SetBinError(static_cast<int>(*bin),std::sqrt(*variance));}
        else if (centralityOne && *histogram == "qa/events/centrality_scaled_bits") {
          centralityOne->SetBinContent(static_cast<int>(*bin), *value);
          centralityOne->SetBinError(static_cast<int>(*bin), std::sqrt(*variance));
        }
      }
      /* the centrality target: the selected trigger, enabled in this source */
      if (centralityOne) {
        int b = -1;
        for (size_t i = 0; i < info.at("name").size(); ++i) if (info.at("name")[i].get<std::string>() == trigger) b = i;
        if (b >= 0 && info.at("initial_prescale")[b].get<int>() >= 0) {
          centralityBindings.push_back({{"run",run},{"bit",b},{"trigger",trigger}});
          for (int x = 1; x <= centralityOne->GetNbinsX(); ++x) {
            target[x - 1] += centralityOne->GetBinContent(x, b + 1);
            targetVariance[x - 1] += std::pow(centralityOne->GetBinError(x, b + 1), 2);
          }
        }
      }
      if (run < firstRun) continue;
      std::vector<int> bits = {bitOf(info.at("name"), reference)};
      for (const auto &p : photonsStudy) bits.push_back(bitOf(info.at("name"), p.second));
      bool usable = true;
      for (int b : bits) usable = usable && b >= 0 && info.at("initial_prescale")[b].get<int>() >= 0;
      if (!usable) continue;
      double referenceEvents = 0;
      for (int x = 1; x <= nE; ++x) referenceEvents += one->GetBinContent(x, bits[0] + 1, 65);
      if (!(referenceEvents > 0)) continue;
      triggerBindings.push_back({{"run",run},{"bits",bits}});
      studiedRuns.insert(run); runs=studiedRuns.size(); firstSeen = firstSeen ? std::min(firstSeen, run) : run; lastSeen = std::max(lastSeen, run);
      for (int x = 1; x <= nE; ++x) {
        refSpectrum[x - 1] += one->GetBinContent(x, bits[0] + 1, 65);
        den[x - 1] += one->GetBinContent(x, bits[0] + 1, 65);
        const double referenceError=one->GetBinError(x,bits[0]+1,65);
        refVariance[x-1]+=referenceError*referenceError; denVariance[x-1]+=referenceError*referenceError;
        for (size_t p = 0; p < photonsStudy.size(); ++p) {
          photonSpectrum[photonsStudy[p].first][x - 1] += one->GetBinContent(x, bits[p + 1] + 1, 65);
          const int y=bits[0]+1, z=bits[p+1]+1;  /* reference trigger scaled, photon trigger live */
          num[photonsStudy[p].first][x - 1] += one->GetBinContent(x,y,z);
          const double e=one->GetBinError(x,bits[p+1]+1,65), en=one->GetBinError(x,y,z);
          photonVariance[photonsStudy[p].first][x-1]+=e*e; numVariance[photonsStudy[p].first][x-1]+=en*en;
        }
      }
    }
    Edges energy;
    for (int b = 1; b <= nE + 1; ++b) energy.push_back(bitsRaw->GetXaxis()->GetBinLowEdge(b));
    TDirectory *d = dirOf(out, "triggers");
    auto uniqueBindings=[](Json values) {std::sort(values.begin(),values.end());values.erase(std::unique(values.begin(),values.end()),values.end());return values;};
    const auto triggerTrace=science::traceOf(bitsRaw).then("enabled_trigger_source_study",
      {{"reference",reference},{"first_run",firstRun},{"selected_run_bits",uniqueBindings(triggerBindings)},
       {"population",perRunStudy ? "per_run_study" : "legacy_raw"},{"energy_edges",energy}});
    auto studyHist=[&](const std::string &name,const std::vector<double> &v,const std::vector<double> &e) {
      auto *h=hist(name,energy,v);
      for (int b=1;b<=h->GetNbinsX();++b) h->SetBinError(b,std::sqrt(e[b-1]));
      return traced(h,triggerTrace.then("trigger_marginal",{{"quantity",name},{"reference","scaled"},{"photon","live"}}));
    };
    put(d, studyHist("reference",refSpectrum,refVariance));
    put(d, studyHist("reference_denominator",den,denVariance));
    for (const auto &p : photonsStudy) {
      put(d, studyHist("spectrum_"+p.first,photonSpectrum[p.first],photonVariance[p.first]));
      put(d, studyHist("numerator_"+p.first,num[p.first],numVariance[p.first]));
    }
    put(d, new TObjString(Form("runs=%d;first=%d;last=%d;reference=%s;correction=%s", runs, firstSeen, lastSeen, reference.c_str(),perRunStudy ? "reference_live_bit" : "legacy_raw")), "study");
    if (centralityBits) {
      Edges centrality;
      for (int b = 1; b <= centralityBits->GetNbinsX() + 1; ++b) centrality.push_back(centralityBits->GetXaxis()->GetBinLowEdge(b));
      put(out, traced(hist("centrality_target", centrality, target, &targetVariance),science::traceOf(centralityBits).then("enabled_trigger_centrality_target",
        {{"trigger",trigger},{"selected_run_bits",uniqueBindings(centralityBindings)},{"unit","DATA event counts"}})));
    }
  }
}

void writeQA(const Settings &s, TDirectory *results) {
  const auto &kAbcdEdges=s.qaAbcdEdges;
  std::map<std::string, std::vector<std::string>> photonSamples;
  std::map<std::string, std::vector<std::string>> classes;
  for (const auto &sys : s.systems) {
    auto &v = photonSamples[sys.package];
    for (const auto &x : sys.simulation) if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
    classes[sys.package].push_back(sys.centrality);
  }
  for (const auto &[package, photons] : photonSamples) {
    Phases phases("qa_" + std::filesystem::path(s.qaOutput).stem().string() + "_" + package,
                  s.name, 0, package);
    phases.start("qa/" + package + "/data_package");
    const std::string base = s.inputDirectory + "/" + package + "/";
    auto data = openIfPresent(base + "data.root");
    TDirectory *out = dirOf(results, "qa/" + package);
    Sums dataSums;
    // Simulation QA is usable before the DATA package is available.
    // Missing DATA stays absent, including its trigger and purity quantities.
    if (data && data->Get("qa")) {
      sumTree(static_cast<TDirectory *>(data->Get("qa")), "", 1.0, dataSums);
      writeTree(out->mkdir("data"), dataSums);
    }
    const std::map<std::string, std::vector<std::string>> families = {{"phosim", photons}, {"jetsim", s.inclusive.count(package) ? s.inclusive.at(package) : std::vector<std::string>()}};
    size_t sampleCount=0;for(const auto& [family,samples]:families)sampleCount+=samples.size();
    phases.start("qa/" + package + "/simulation_packages",sampleCount,"sample_package");
    std::map<std::string, Sums> familySums;
    std::map<std::string, std::map<std::string, std::map<std::string, TH1D *>>> generator, reco;  /* family -> class -> sample */
    std::map<std::string, std::unique_ptr<TH2>> leakageNative;  /* per class: the simulated signal leakage (A-D x reco ET) */
    /* per class: the leading reconstructed recoil jet of the photon+jet simulation, every one and those matched to the
       leading truth recoil jet (reco ET x xJ on the packages' bins), for the leading-jet association of IAN Fig. 25.
       Packages since the second refill hold them. */
    std::map<std::string, std::array<std::unique_ptr<TH2>, 2>> leadingNative;
    std::map<std::string,std::map<std::string,std::unique_ptr<TH2>>> leadingUniform;
    Windows windows;  /* per sample stitched by ownership: its window of the generator spectrum */
    const bool ppReference = package == "PP" && s.ppStitching == "reference";
    /* The stitching figure shows each sample over its whole generated range, so the turn-on below its window is
       seen: unit counts before the stitching windows (qa/generator/raw_leading_pt) times the cross section per
       generated event. p+p in reference mode requires them. An Au+Au package holds them since the refill of
       2026-10-12; an older one only has the analysis-weighted spectrum inside the producer's own windows. */
    std::map<std::string, bool> rawCounts;  /* per family */
    /* Au+Au photon+jet samples are stitched on the generator-filter photon, a sharp cut, while the IAN draws the
       maximum signal truth-photon pT, where a sample has a turn-on. Packages since the second refill hold that
       coordinate (kind 3 of raw_leading_pt: every event; kind 4: the events the sample owns). */
    std::map<std::string, std::map<std::string, TH1D *>> ownedGenerator;  /* family -> sample */
    std::map<std::string, std::string> generatorCoordinate;
    std::map<std::string, double> interactionMix;  /* per family: the double-interaction fraction of a p+p nominal generator spectrum */
    for (const auto &[family, samples] : families)
      for (const auto &sample : samples) {
        const auto progress=phases.attempt();
        auto f = openIfPresent(base + sample + ".root");
        if (!f || !f->Get("qa")) continue;
        const double w = weightOf(*f);
        /* a sample stitched by ownership (Au+Au embedded simulation): its window of the generator spectrum */
        const Json packageMetadata = Json::parse(static_cast<TObjString *>(f->Get("meta/package"))->GetString().Data());
        const Json normalization = packageMetadata.value("normalization", Json::object());
        const bool referenceSource = !ppReference || packageMetadata.value("si_di_role", "") == "SI";
        if (normalization.contains("ownership")) {
          const Json &window = normalization.at("ownership").at("window_gev");
          windows[sample] = {window[0].get<double>(), window[1].is_null() ? std::numeric_limits<double>::infinity() : window[1].get<double>()};
        }
        sumTree(static_cast<TDirectory *>(f->Get("qa")), "", w, familySums[family]);
        if (referenceSource) {
          auto *raw = dynamic_cast<TH2 *>(f->Get("qa/generator/raw_leading_pt"));
          /* p+p in nominal mode, where the package holds the raw generator counts: drawn as the reference mode draws
             them (each sample over its whole generated range, in its generator window of the sum), with the
             single- and the double-interaction part of every sample, weighted by the simulation's period mix */
          const bool ppCounts = package == "PP" && !ppReference && raw && normalization.contains("generator_diagnostic");
          const bool counts = ppReference || ppCounts || (package == "AuAu" && raw);
          auto *g = counts ? raw : dynamic_cast<TH2 *>(f->Get("qa/generator/leading_pt"));
          if (ppReference && (!g || !normalization.contains("generator_diagnostic")))
            throw std::runtime_error(sample + ": reference stitching needs raw generator counts; analysis weights cannot substitute");
          if (g) {
            if (rawCounts.count(family) && rawCounts.at(family) != counts)
              throw std::runtime_error(package + " " + family + ": only some samples hold raw generator counts");
            rawCounts[family] = counts;
            const bool truthPhoton = counts && package == "AuAu" && family == "phosim" && g->GetNbinsY() >= 4;
            const int kind = truthPhoton ? 3 : family == "phosim" ? 1 : 2;
            generatorCoordinate[family] = truthPhoton ? "maximum signal truth-photon pT"
                                          : family == "phosim" ? "leading prompt truth-photon pT" : "leading truth-jet pT";
            auto *p = g->ProjectionX(("gen_" + sample).c_str(), kind, kind, "e");
            p->SetDirectory(nullptr);
            auto spectrumTrace=science::projectionOf(g,{{{"axis",1},{"low",g->GetYaxis()->GetBinLowEdge(kind)},{"high",g->GetYaxis()->GetBinUpEdge(kind)}}});
            if (truthPhoton) {
              auto *o = g->ProjectionX(("gen_owned_" + sample).c_str(), 4, 4, "e");
              o->SetDirectory(nullptr);
              o->Scale(w);
              traced(o,science::projectionOf(g,{{{"axis",1},{"low",g->GetYaxis()->GetBinLowEdge(4)},{"high",g->GetYaxis()->GetBinUpEdge(4)}}}).then("sample_normalization",{{"factor",w}}));
              ownedGenerator[family][sample] = o;
            }
            double scale = w;
            Json mix;
            if (ppReference || ppCounts) {
              const auto &diagnostic = normalization.at("generator_diagnostic");
              if (diagnostic.at("source_moment") != "raw_counts") throw std::runtime_error("generator stitching requires raw counts");
              scale = diagnostic.at("cross_section_pb").get<double>() / diagnostic.at("generated_events").get<double>();
              const auto &window = diagnostic.at("window_gev");
              windows[sample] = {window[0].get<double>(), window[1].is_null() ? std::numeric_limits<double>::infinity() : window[1].get<double>()};
            }
            if (ppCounts) {
              /* the double-interaction fraction of the simulation: the period shares times their double fractions,
                 as the producer weights its events (the truth-vertex weights are not applied at generator level) */
              const Json &periods = packageMetadata.at("render_selection").at("photon_id").at("simulation_periods").at("periods");
              double shares = 0, doubled = 0;
              for (const auto &[period, p] : periods.items())
                shares += p.at("share").get<double>(), doubled += p.at("share").get<double>() * p.at("double_fraction").get<double>();
              const std::string role = packageMetadata.value("si_di_role", "");
              if (!(shares > 0) || (role != "SI" && role != "DI")) throw std::runtime_error(sample + ": no SI/DI period mix");
              const double fraction = doubled / shares;
              scale *= role == "DI" ? fraction : 1 - fraction;
              mix = {{"si_di_role",role},{"double_fraction",fraction},{"from","simulation period shares x double_fraction"},
                     {"truth_vertex_weights","not applied"}};
              interactionMix[family] = fraction;
            }
            p->Scale(scale);
            traced(p,spectrumTrace.then("sample_normalization",mix.is_null() ? Json{{"factor",scale}} : Json{{"factor",scale},{"interaction_mix",mix}}));
            generator[family]["all"][sample] = p;
          }
        }
        for (const auto &c : classes[package]) {
          TH1D *regionA = nullptr;
          for (const char *label : family == "phosim" ? std::vector<const char *>{"sig", "phoother"} : std::vector<const char *>{"bkg", "jetprompt"})
            if (auto *a = dynamic_cast<TH2 *>(f->Get(("qa/photons/" + c + "/" + label + "/abcd").c_str()))) {
              auto *p = a->ProjectionY(("reco_" + sample + label).c_str(), 1, 1, "e");
              p->SetDirectory(nullptr);
              const auto part=science::projectionOf(a,{{{"axis",0},{"low",a->GetXaxis()->GetBinLowEdge(1)},{"high",a->GetXaxis()->GetBinUpEdge(1)}}});
              traced(p,part);
              if (!regionA) regionA = p; else {
                const auto combined=Trace::combine("sum_candidate_labels",{science::traceOf(regionA),part});
                regionA->Add(p); traced(regionA,combined); delete p;
              }
            }
          if (ppReference) {
            delete regionA;
            regionA = nullptr;
            // The historical reconstructed-tag comparison starts at Jet12;
            // Jet8 remains part of the generator stitching diagnostic.
            if (referenceSource && normalization.at("generator_diagnostic").at("sample_slice") != "jet8") {
              auto *saved = dynamic_cast<TH1 *>(f->Get("qa/stitching/reco_tag_vz30"));
              if (!saved) throw std::runtime_error(sample + ": missing |vz| < 30 cm tag comparison");
              regionA = new TH1D(("reference_"+sample).c_str(), "", saved->GetNbinsX(), saved->GetXaxis()->GetXmin(), saved->GetXaxis()->GetXmax());
              regionA->SetDirectory(nullptr); regionA->Sumw2();
              for (int b=0;b<=saved->GetNbinsX()+1;++b) regionA->SetBinContent(b,saved->GetBinContent(b)), regionA->SetBinError(b,saved->GetBinError(b));
              traced(regionA,science::traceOf(saved).then("copy_reference_tag_bins"));
            }
          }
          if (regionA) { regionA->Scale(w); traced(regionA,science::traceOf(regionA).then("sample_normalization",{{"factor",w}})); reco[family][c][sample] = regionA; }
          if (family == "phosim")
            if (auto *l = dynamic_cast<TH2 *>(f->Get(("response/" + c + "/" + s.jet + "/signal_leakage_truth_matched_abcd").c_str()))) {
              auto &L = leakageNative[c];
              const bool first=!L;
              if (!L) L.reset(static_cast<TH2 *>(l->Clone(("leakage_" + c).c_str()))), L->SetDirectory(nullptr), L->Scale(w);
              else L->Add(l, w);
              science::sumDefinition(L.get(),l,w,first);
            }
          if (family == "phosim")
            if(package=="PP")for(const std::string name:{"leading_reco_u","leading_matched_reco_u","leading_matched_any_reco_u"})
              if(auto *l=dynamic_cast<TH2*>(f->Get(("response/"+c+"/"+s.jet+"/"+(name=="leading_reco_u" ? "xj_leading_global_reco" : name=="leading_matched_reco_u" ? "xj_leading_global_matched_reco" : "xj_"+name)).c_str()))){
                auto &L=leadingUniform[c][name];
                const bool first=!L;
                if(!L)L.reset(static_cast<TH2*>(l->Clone((name+"_"+c).c_str()))),L->SetDirectory(nullptr),L->Scale(w);
                else if(!L->Add(l,w))throw std::runtime_error("uniform association QA binning differs");
                science::sumDefinition(L.get(),l,w,first);
              }
          if (family == "phosim")
            for (int kind = 0; kind < 2; ++kind)
              if (auto *l = dynamic_cast<TH2 *>(f->Get(("response/" + c + "/" + s.jet + (kind ? "/xj_leading_global_matched_reco" : "/xj_leading_global_reco")).c_str()))) {
                auto &L = leadingNative[c][kind];
                const bool first=!L;
                if (!L) L.reset(static_cast<TH2 *>(l->Clone((std::string(kind ? "leading_matched_" : "leading_") + c).c_str()))), L->SetDirectory(nullptr), L->Scale(w);
                else L->Add(l, w);
                science::sumDefinition(L.get(),l,w,first);
              }
        }
      }
    phases.start("qa/" + package + "/derived_quantities");
    const auto sampleSelection=[&](Trace trace,const std::string &family) {
      if(family=="jetsim")return withSettings(trace,s.parameterOrigins,{"QA."+package+".Inclusive"});
      for(const auto &sys:s.systems)if(sys.package==package)trace=withSettings(trace,sys.parameterOrigins,{"Simulation"});
      return trace;
    };
    for(auto &[family,sums]:familySums)
      bindSampleSelection(sums,dirOf(out,family),[&](Trace trace){return sampleSelection(trace,family);});
    for(auto *spectra:{&generator,&reco})for(auto &[family,byClass]:*spectra)for(auto &[c,parts]:byClass)for(auto &[sample,h]:parts) {
      auto trace=sampleSelection(science::traceOf(h),family);
      if(package=="PP")trace=withSettings(trace,s.parameterOrigins,{"QA.PP.Stitching"});
      science::attachTrace(h,trace);
    }
    for(auto &[family,parts]:ownedGenerator)for(auto &[sample,h]:parts) {
      auto trace=sampleSelection(science::traceOf(h),family);
      if(package=="PP")trace=withSettings(trace,s.parameterOrigins,{"QA.PP.Stitching"});
      science::attachTrace(h,trace);
    }
    for(auto &[c,h]:leakageNative)science::attachTrace(h.get(),sampleSelection(science::traceOf(h.get()),"phosim"));
    for(auto &[c,rows]:leadingUniform)for(auto &[name,h]:rows)science::attachTrace(h.get(),sampleSelection(science::traceOf(h.get()),"phosim"));
    for(auto &[c,rows]:leadingNative)for(auto &h:rows)if(h)science::attachTrace(h.get(),sampleSelection(science::traceOf(h.get()),"phosim"));
    for (auto &[family, sums] : familySums) writeTree(dirOf(out,family), sums);
    /* stitching: generator truth and the region-A tag, per family */
    for (const auto &[family, byClass] : generator)
      for (const auto &[c, parts] : byClass) {
        auto *directory = dirOf(out, "stitching/" + family + "_generator");
        /* the range of the descriptive fit: where the stitched spectrum is populated by more than one sample's tail
           (photons 15 to 40 GeV; inclusive jets to 50 GeV, from 10 GeV in p+p and 15 GeV in Au+Au) */
        const bool ppFull = ppReference || interactionMix.count(family);
        const double low = ppFull ? (family == "phosim" ? 5 : 9) : family == "phosim" ? 15 : package == "PP" ? 10 : 15;
        const double high = family == "phosim" ? 40 : 50;
        const bool counts = rawCounts.count(family) && rawCounts.at(family);
        stitching(directory, parts, low, high, windows, true, !counts, ownedGenerator.count(family) ? ownedGenerator.at(family) : std::map<std::string, TH1D *>());
        Json recipe={{"mode",ppReference ? "reference" : "nominal"},{"range_gev",{low,high}},{"coordinate",generatorCoordinate[family]},
          {"weight_channel",counts ? "raw_counts_times_cross_section_over_generated" : "downstream_analysis_weighted"}};
        /* The events of the spectrum, for the figure's header. p+p counts: every recorded event. Au+Au: the events
           the producer classes as analysis events before the stitching windows, 0-80% centrality inside its own
           vertex window (TreeProduction vertex_cut_cm, 60 cm), which is wider than the analysis one. */
        recipe["population"] = package == "PP" ? (counts ? "all generated events" : "events admitted by the source stitching") : "0-80%, |z_vtx| < 60 cm";
        if (interactionMix.count(family))
          recipe["interactions"] = {{"parts","single and double interaction"},{"double_fraction",interactionMix.at(family)}};
        directory->WriteTObject(new TObjString(recipe.dump().c_str()),"recipe");
      }
    const std::string central = *std::min_element(classes[package].begin(), classes[package].end(), [](const std::string &a, const std::string &b) {
      return (a.rfind("cent", 0) == 0 ? std::atof(a.c_str() + 4) : 0) < (b.rfind("cent", 0) == 0 ? std::atof(b.c_str() + 4) : 0);
    });
    for (const auto &[family, byClass] : reco)
      if (byClass.count(central)) {
        auto *directory=dirOf(out, "stitching/" + family + "_reco");
        stitching(directory, byClass.at(central), 15, 35, {}, package != "PP", !ppReference);
        Json recipe={{"mode",ppReference ? "reference" : "nominal"},{"range_gev",{15,35}},
          {"weight_channel","downstream_analysis_weighted"}};
        if (ppReference) recipe["si_di_role"]="SI", recipe["vertex_abs_z_max_cm"]=30;
        directory->WriteTObject(new TObjString(recipe.dump().c_str()),"recipe");
      }
    /* DATA purity and the simulation closure, per class, on the ABCD bins */
    std::map<std::string, std::vector<std::vector<double>>> leakageByClass, leakageVarianceByClass;
    for(const auto &[c,names]:leadingUniform)for(const auto &[name,h]:names)dirOf(out,"association/"+c)->WriteTObject(h.get(),name.c_str());
    for (auto &[c, pair] : leadingNative)
      if (pair[0] && pair[1]) {
        TDirectory *d = dirOf(out, "association/" + c);
        d->WriteTObject(pair[0].get(), "leading_reco");
        d->WriteTObject(pair[1].get(), "leading_matched_reco");
      }
    for (const auto &[c, h] : leakageNative)
      for (int r = 0; r < 4; ++r) {
        leakageByClass[c].push_back(regroup(h.get(), r, kAbcdEdges));
        leakageVarianceByClass[c].push_back(regroup(h.get(), r, kAbcdEdges, true));
      }
    for (const auto &c : classes[package]) {
      auto *abcd = data ? dynamic_cast<TH2 *>(data->Get(("measurement/" + c + "/abcd").c_str())) : nullptr;
      if (abcd) {
        auto *counts = new TH2D(c.c_str(), "event-leading candidates;region (A,B,C,D);photon E_{T} [GeV]", 4, 0, 4, kAbcdEdges.size() - 1, kAbcdEdges.data());
        for (int r = 0; r < 4; ++r) {
          const auto v = regroup(abcd, r, kAbcdEdges), e2 = regroup(abcd, r, kAbcdEdges, true);
          for (size_t i = 0; i < v.size(); ++i) counts->SetBinContent(r + 1, i + 1, v[i]), counts->SetBinError(r + 1, i + 1, std::sqrt(e2[i]));
        }
        counts->SetDirectory(nullptr);
        traced(counts,withSettings(grouped(abcd,kAbcdEdges),s.parameterOrigins,{"QA.AbcdEdges"}));
        put(dirOf(out, "abcd"), counts);
      }
      if (abcd && leakageByClass.count(c)) {
        std::vector<std::vector<double>> counts(4), variance(4);
        for (int r = 0; r < 4; ++r) counts[r] = regroup(abcd, r, kAbcdEdges), variance[r] = regroup(abcd, r, kAbcdEdges, true);
        const auto &leakage = leakageByClass[c], &leakageVariance = leakageVarianceByClass[c];
        const auto [raw, corrected] = purityRows(counts, leakage);
        /* The statistical errors of the DATA purity: each region's count of a bin moved up and down by its error in
           turn, half the difference added in quadrature. The leakage is the simulation's and is held. */
        const size_t n = kAbcdEdges.size() - 1;
        std::vector<double> rawVariance(n), correctedVariance(n);
        for (int r = 0; r < 4; ++r)
          for (size_t i = 0; i < n; ++i) {
            auto up = counts, down = counts;
            up[r][i] += std::sqrt(variance[r][i]), down[r][i] -= std::sqrt(variance[r][i]);
            const auto [rawUp, correctedUp] = purityRows(up, leakage);
            const auto [rawDown, correctedDown] = purityRows(down, leakage);
            const double dRaw = 0.5 * (rawUp[i] - rawDown[i]), dCorrected = 0.5 * (correctedUp[i] - correctedDown[i]);
            if (std::isfinite(dRaw)) rawVariance[i] += dRaw * dRaw;
            if (std::isfinite(dCorrected)) correctedVariance[i] += dCorrected * dCorrected;
          }
        auto withErrors = [&](const char *name, const std::vector<double> &values, const std::vector<double> &squared) {
          TH1D *h = hist(name, kAbcdEdges, values);
          for (size_t i = 0; i < n; ++i) h->SetBinError(i + 1, std::isfinite(values[i]) ? std::sqrt(squared[i]) : 0.0);
          return h;
        };
        TDirectory *d = dirOf(out, "purity/" + c);
        const auto countTrace=withSettings(grouped(abcd,kAbcdEdges),s.parameterOrigins,{"QA.AbcdEdges"}),leakageTrace=withSettings(grouped(leakageNative.at(c).get(),kAbcdEdges),s.parameterOrigins,{"QA.AbcdEdges"});
        put(d, traced(withErrors("raw", raw, rawVariance),countTrace.then("raw_abcd_purity",{{"variance","half_difference_each_DATA_cell_error"}})));
        put(d, traced(withErrors("corrected", corrected, correctedVariance),Trace::combine("leakage_corrected_abcd_purity",{countTrace,leakageTrace},
          {{"variance","half_difference_each_DATA_cell_error"},{"leakage_uncertainty","held_fixed"}})));
        /* the signal leaking into B, C and D over the signal in A, with the simulation's statistical error */
        for (int r = 1; r < 4; ++r) {
          std::vector<double> fraction(n, std::nan("")), squared(n);
          for (size_t i = 0; i < n; ++i) {
            if (!(leakage[0][i] > 0)) continue;
            const auto estimate=signalLeakageRatio(leakage[r][i],leakage[0][i],
                                                    leakageVariance[r][i],leakageVariance[0][i],r==0);
            fraction[i]=estimate.value;squared[i]=estimate.variance;
          }
          put(d, traced(withErrors(Form("leakage_%c", "ABCD"[r]), fraction, squared),leakageTrace.then("signal_region_over_A",{{"region",r},{"variance","independent_regions_delta_method"}})));
        }
      }
      /* closure: the simulated pseudo-data (kSignalLabel, kBackgroundLabel: one source of prompt photons) on the closure
         bins, estimated as the data are (R_bkg = 1) and with the simulated background correlation of the result's photon
         bins (the value BackgroundCorrelation: qa gives the analysis), each over the true purity. Statistical errors by
         linear propagation: every weighted count (signal and background, each region and native E_T bin) moved up and
         down by its error in turn, half the difference added in quadrature; the leakage is built from the same signal
         candidates and moves with them, and R_bkg is recomputed from the moved background. */
      auto labelled = [&](const char *family, const char *label) -> TH2 * {
        if (!familySums.count(family)) return nullptr;
        auto it = familySums[family].find("photons/" + c + "/" + label + "/abcd");
        if (it == familySums[family].end() || !filled(static_cast<TH1 *>(it->second.get()))) return nullptr;
        return static_cast<TH2 *>(it->second.get());
      };
      TH2 *signalCounts = labelled("phosim", kSignalLabel), *backgroundCounts = labelled("jetsim", kBackgroundLabel);
      if (leakageByClass.count(c) && signalCounts && backgroundCounts && signalCounts->GetNbinsY() == backgroundCounts->GetNbinsY()) {
        const SystemSettings *result = nullptr;
        for (const auto &sys : s.systems) if (sys.package == package && sys.centrality == c) result = &sys;
        const int nNative = signalCounts->GetNbinsY();
        Edges native(nNative + 1);
        for (int j = 0; j <= nNative; ++j) native[j] = signalCounts->GetYaxis()->GetBinLowEdge(j + 1);
        using Counts = std::vector<std::vector<double>>;  /* region x bin */
        Counts signal(4, std::vector<double>(nNative)), background = signal, signalError = signal, backgroundError = signal;
        for (int r = 0; r < 4; ++r)
          for (int j = 0; j < nNative; ++j) {
            signal[r][j] = signalCounts->GetBinContent(r + 1, j + 1), signalError[r][j] = signalCounts->GetBinError(r + 1, j + 1);
            background[r][j] = backgroundCounts->GetBinContent(r + 1, j + 1), backgroundError[r][j] = backgroundCounts->GetBinError(r + 1, j + 1);
          }
        auto into = [&](const Counts &fine, const Edges &target) {  /* native bins summed into target bins */
          Counts out(4, std::vector<double>(target.size() - 1));
          for (int r = 0; r < 4; ++r)
            for (int j = 0; j < nNative; ++j) {
              const int i = std::upper_bound(target.begin(), target.end(), native[j] + 1e-9) - target.begin() - 1;
              if (i >= 0 && i < static_cast<int>(target.size()) - 1) out[r][i] += fine[r][j];
            }
          return out;
        };
        const size_t n = kAbcdEdges.size() - 1;
        const Counts centralSignal = into(signal, kAbcdEdges);
        auto closure = [&](const Counts &movedSignal, const Counts &movedBackground) {
          const Counts sg = into(movedSignal, kAbcdEdges), bg = into(movedBackground, kAbcdEdges);
          Counts sum(4, std::vector<double>(n)), leakage = leakageByClass[c];
          for (int r = 0; r < 4; ++r)
            for (size_t i = 0; i < n; ++i) {
              sum[r][i] = sg[r][i] + bg[r][i];
              if (centralSignal[r][i] > 0) leakage[r][i] *= sg[r][i] / centralSignal[r][i];
            }
          const auto [raw, corrected] = purityRows(sum, leakage);
          std::vector<double> truth(n, std::nan("")), ratio = truth, correctedR = truth, ratioR = truth;
          for (size_t i = 0; i < n; ++i)
            if (sum[0][i] > 0) truth[i] = sg[0][i] / sum[0][i], ratio[i] = corrected[i] / truth[i];
          if (result) {
            const Counts b = into(movedBackground, result->photonEdges);
            const int nPhoton = result->photonEdges.size() - 1;
            for (size_t i = 0; i < n; ++i) {
              const double centre = 0.5 * (kAbcdEdges[i] + kAbcdEdges[i + 1]);  /* a closure bin outside the photon bins takes the nearest */
              int p = std::upper_bound(result->photonEdges.begin(), result->photonEdges.end(), centre) - result->photonEdges.begin() - 1;
              p = std::min(std::max(p, 0), nPhoton - 1);
              if (!(b[0][p] > 0 && b[1][p] > 0 && b[2][p] > 0 && b[3][p] > 0) || !(sum[0][i] > 0) || !(leakage[0][i] > 0)) continue;
              const double counts[4] = {std::max(sum[0][i], 0.0), std::max(sum[1][i], 0.0), std::max(sum[2][i], 0.0), std::max(sum[3][i], 0.0)};
              const double leak[4] = {leakage[0][i], leakage[1][i], leakage[2][i], leakage[3][i]};
              correctedR[i] = solveLeakageABCD(counts, leak, b[0][p] * b[3][p] / (b[1][p] * b[2][p])).purity;
              ratioR[i] = correctedR[i] / truth[i];
            }
          }
          return std::array<std::vector<double>, 6>{raw, corrected, truth, ratio, correctedR, ratioR};
        };
        const auto central = closure(signal, background);
        std::array<std::vector<double>, 6> error2;
        for (auto &e : error2) e.assign(n, 0.0);
        for (const bool ofSignal : {true, false})
          for (int r = 0; r < 4; ++r)
            for (int j = 0; j < nNative; ++j) {
              const double error = (ofSignal ? signalError : backgroundError)[r][j];
              if (!(error > 0)) continue;
              std::array<std::array<std::vector<double>, 6>, 2> moved;
              for (int side = 0; side < 2; ++side) {
                Counts movedSignal = signal, movedBackground = background;
                double &x = (ofSignal ? movedSignal : movedBackground)[r][j];
                x = std::max(0.0, x + (side == 0 ? error : -error));
                moved[side] = closure(movedSignal, movedBackground);
              }
              for (int q = 0; q < 6; ++q)
                for (size_t i = 0; i < n; ++i) {
                  const double half = 0.5 * (moved[0][q][i] - moved[1][q][i]);
                  if (std::isfinite(half)) error2[q][i] += half * half;
                }
            }
        TDirectory *d = dirOf(out, "closure/" + c);
        const char *names[] = {"raw", "corrected", "truth", "ratio", "corrected_r_bkg", "ratio_r_bkg"};
        for (int q = 0; q < (result ? 6 : 4); ++q) {
          TH1D *h = hist(names[q], kAbcdEdges, central[q]);
          for (size_t i = 0; i < n; ++i) h->SetBinError(i + 1, std::isfinite(central[q][i]) ? std::sqrt(error2[q][i]) : 0.0);
          std::vector<Trace> ancestry{grouped(signalCounts,kAbcdEdges),grouped(backgroundCounts,kAbcdEdges)};
          if (q!=0 && q!=2) ancestry.push_back(grouped(leakageNative.at(c).get(),kAbcdEdges));
          Json method={{"quantity",names[q]},{"variance","half_difference_each_native_weighted_cell"},{"signal_leakage","shared_displacement"}};
          if (q>=4) { ancestry.push_back(grouped(backgroundCounts,result->photonEdges));method["correlation_bins"]=result->photonEdges; }
          traced(h,withSettings(Trace::combine("simulation_purity_closure",ancestry,method),s.parameterOrigins,{"QA.AbcdEdges"}));
          if(q>=4)science::attachTrace(h,withSettings(science::traceOf(h),result->parameterOrigins,{"PhotonEdges"}));
          put(d, h);
        }
      }
    }
    if (data)
      for (const auto &sys : s.systems)
        if (sys.package == package) leadingData(*data, s, sys, dirOf(out, "association/" + sys.centrality));
    /* the purity study of every result of the package */
    for (const auto &sys : s.systems)
      if (sys.package == package)
        purityStudy(dirOf(out, "purity_study/" + sys.name), sys, leakageNative.count(sys.centrality) ? leakageNative[sys.centrality].get() : nullptr,
                    familySums["phosim"], familySums["jetsim"], dataSums);
    /* working points of the photon samples (Au+Au: per 5% slice) */
    if (package == "AuAu" && familySums.count("phosim")) {
      for (const auto &[what, path, above] : std::vector<std::tuple<std::string, std::string, bool>>{{"bdt", "slices/sig/bdt", true}, {"isolation", "slices/sig/isoslice", false}}) {
        auto it = familySums["phosim"].find(path);
        if (it != familySums["phosim"].end()) workingPoints(dirOf(out, "working_points/" + what), static_cast<TH3 *>(it->second.get()), above);
      }
    }
    /* trigger study and the centrality target: per DATA source, with its own trigger table. A compact package
       (final_analysis.sh compact) carries them in derived/, made by this same function. */
    if (!data) continue;
    if (auto *derived = data->GetDirectory("derived")) {
      auto *version = dynamic_cast<TObjString *>(derived->Get("version"));
      if (!version || version->GetString() != kSourceStudies)
        throw std::runtime_error(base + "data.root: its derived studies come from another version of the analysis; run final_analysis.sh compact again");
      copyDirectory(derived, out, {"version"});
      // V2's target is a sum of unit DATA event fills (cent_bits_data), before
      // any source normalization. Old compact packages saved its counts only.
      // Upgrade this one documented input; the general output guard still
      // refuses every other histogram with implicit errors.
      if (auto* target=dynamic_cast<TH1*>(out->Get("centrality_target")); target && !target->GetSumw2N()) {
        std::unique_ptr<TH1> owned(target); target->SetDirectory(nullptr);
        materializeUnitCountErrors(*target);
        out->WriteTObject(target,"centrality_target","Overwrite");
        put(out,new TObjString("PhotonJetSourceStudiesV2 unit DATA event counts; variance=count; upgraded legacy compact representation"),"centrality_target_uncertainty");
      }
    } else sourceStudies(data.get(), package, out);
    phases.complete();
  }
  results->cd();
}

/* The QA file: written to Output.QA.partial and moved into place; an existing file is never replaced. */
/* Copy the producer's saved model diagnostics without recomputing training
   statistics. A model release, its diagnostics and its collection must agree. */
void savedModelQA(const std::string &collection,TDirectory *output) {
  namespace fs=std::filesystem;
  const auto index=photonjet_config::readYaml(collection);
  Json body=index;body.erase("collection_sha256");
  if (index.value("schema","")!="PhotonIDModelQACollectionV1" || index.at("models").empty() ||
      photonjet_render::sha256(body.dump(-1,' ',true))!=index.at("collection_sha256").get<std::string>())
    throw std::runtime_error("invalid or changed saved model-QA collection");
  const fs::path base=fs::canonical(fs::path(collection).parent_path());
  const fs::path modelRoot=base.parent_path().filename()==".campaigns" ? base.parent_path().parent_path() : base;
  auto readJson=[](const fs::path &path) {std::ifstream file(path);if(!file)throw std::runtime_error("missing model-QA metadata");return Json::parse(file);};
  std::set<std::string> names;
  std::map<std::string,Json> recipes, metadataByName;
  const std::string system=index.at("system");
  if(system!="AuAu" && system!="PP")throw std::runtime_error("unknown model-QA collision system");
  auto *models=output->mkdir("model_qa");
  // The pp comparison imports its separately sealed AuAu reference into a
  // distinct directory: both systems may use the same nominal model name.
  auto import=[&](const Json &entry,const std::string &expectedSystem,TDirectory *parent,bool reference=false) {
    const std::string name=entry.at("name");
    if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
      throw std::runtime_error("invalid or duplicate model-QA name");
    const fs::path folder=fs::canonical(base/entry.at("path").get<std::string>());
    if (folder.lexically_relative(modelRoot).empty() || *folder.lexically_relative(modelRoot).begin()=="..") throw std::runtime_error("model-QA release escapes models root");
    const Json release=readJson(folder/"manifest.json"), diagnostics=readJson(folder/"diagnostics.json");
    Json releaseBody=release;releaseBody.erase("release_sha256");
    if (release.value("schema","")!="PhotonIDReleaseV1" || release.value("status","")!="PASS" ||
        release.at("release_sha256")!=entry.at("release_sha256") ||
        photonjet_render::sha256(releaseBody.dump(-1,' ',true))!=release.at("release_sha256").get<std::string>() ||
        photonjet_render::sha256((folder/"diagnostics.json").string(),true)!=entry.at("diagnostics_json_sha256").get<std::string>() ||
        photonjet_render::sha256((folder/"diagnostics.root").string(),true)!=entry.at("diagnostics_root_sha256").get<std::string>() ||
        diagnostics.at("root_sha256")!=entry.at("diagnostics_root_sha256") ||
        release.at("system")!=expectedSystem || diagnostics.at("system")!=expectedSystem ||
        diagnostics.at("et_ranges_gev")!=release.at("recipe").at("et_ranges_gev") ||
        diagnostics.at("identity").at("request_sha256")!=release.at("request_sha256") ||
        diagnostics.at("identity").at("recipe_sha256")!=photonjet_render::sha256(release.at("recipe").dump(-1,' ',true)))
      throw std::runtime_error("model-QA inputs are not bound to their sealed release");
    for (const auto &[name,expected] : release.at("artifacts").items()) {
      const fs::path artifact=fs::canonical(folder/name);
      if (artifact.lexically_relative(folder).empty() || *artifact.lexically_relative(folder).begin()==".." ||
          photonjet_render::sha256(artifact.string(),true)!=expected.get<std::string>()) throw std::runtime_error("sealed model-QA artifact changed");
    }
    if (entry.at("name")!=release.at("name") || entry.at("name")!=diagnostics.at("name")) throw std::runtime_error("model-QA names differ");
    if (!reference && (entry.at("axes")!=release.at("recipe").at("axes") ||
        entry.at("et_ranges_gev")!=diagnostics.at("et_ranges_gev") ||
        entry.at("max_depth")!=release.at("recipe").at("learner").at("max_depth")))
      throw std::runtime_error("model-QA routing/depth settings differ from release");
    if(!reference) {
      recipes[name]=release.at("recipe");metadataByName[name]=diagnostics;
      for (const char *field : {"population_sha256","split_sha256","event_split_sha256","training_weights_sha256"})
      if (diagnostics.at("identity").at(field)!=entry.at("pins").at(field)) throw std::runtime_error("model-QA population/split/weight pin differs");
    // Retained fits keep the training implementation pin; subsequent reporting
    // code has its own pin and does not change which learner was trained.
    const Json fitCode=release.contains("retained_fit") && release.at("retained_fit").contains("fit_code_sha256")
      ? release.at("retained_fit").at("fit_code_sha256") : release.at("code_sha256");
    if(entry.contains("reporting_code_sha256") && entry.at("reporting_code_sha256")!=release.at("code_sha256"))
      throw std::runtime_error("model-QA reporting implementation pin differs");
    const Json implementation={{"code",fitCode},{"environment",release.at("environment")},{"ROOT",diagnostics.at("ROOT")}};
    if (entry.at("pins").at("sources_sha256")!=photonjet_render::sha256(diagnostics.at("sources").dump(-1,' ',true)) ||
        entry.at("pins").at("implementation_sha256")!=photonjet_render::sha256(implementation.dump(-1,' ',true))) throw std::runtime_error("model-QA source or implementation pin differs");
    }
    TFile input((folder/"diagnostics.root").c_str(),"READ");if(input.IsZombie()) throw std::runtime_error("cannot read model diagnostics");
    observeNativeInput(input);
    auto *saved=dynamic_cast<TObjString*>(input.Get("meta/diagnostics"));if(!saved)throw std::runtime_error("missing saved model-QA metadata");
    const Json metadata=Json::parse(saved->GetString().Data());
    for(const char *field:{"identity","name","system","feature_order","et_ranges_gev","max_depth"})
      if(metadata.at(field)!=diagnostics.at(field))throw std::runtime_error("model-QA ROOT and JSON metadata differ");
    for(const char *table:{"features","performance","holdout_et_c3"})if(!dynamic_cast<TTree*>(input.Get(table)))throw std::runtime_error("missing saved model-QA table");
    auto *dest=parent->mkdir(reference ? "auau" : name.c_str());
    std::function<void(TDirectory*,TDirectory*)> copy=[&](TDirectory *from,TDirectory *to) {
      TIter next(from->GetListOfKeys());while(auto *key=static_cast<TKey*>(next())) {
        std::unique_ptr<TObject> object(key->ReadObj());
        if(auto *dir=dynamic_cast<TDirectory*>(object.get())) copy(dir,to->mkdir(key->GetName()));
        else if(auto *tree=dynamic_cast<TTree*>(object.get())) {
          to->cd();auto *clone=tree->CloneTree(-1,"fast");if(!clone)throw std::runtime_error("cannot copy model-QA records");
          science::attachTrace(clone,science::traceOf(tree));
          const auto columns=science::tag(tree,"physics_column_traces");
          if(!columns.is_null())science::setTag(clone,"physics_column_traces",columns);
          clone->Write(key->GetName());clone->SetDirectory(nullptr);delete clone;
        }
        else {to->cd();object->Write(key->GetName());}
      }
    };
    copy(&input,dest);
    // Keep the released calibration beside the saved diagnostics. Rendering
    // its thresholds must not derive a different set of working points.
    if (!reference && expectedSystem == "AuAu") {
      dest->GetDirectory("meta")->cd();
      TObjString(readJson(folder/"working_points.json").dump().c_str()).Write("working_points");
    }
    if (photonjet_render::sha256((folder/"diagnostics.root").string(),true)!=entry.at("diagnostics_root_sha256").get<std::string>()) throw std::runtime_error("model diagnostics changed during QA import");
    return std::make_pair(release.at("recipe"),diagnostics);
  };
  for(const auto &entry:index.at("models")) {
    if(!names.insert(entry.at("name")).second)throw std::runtime_error("duplicate model-QA name");
    if(entry.at("system")!=system)throw std::runtime_error("model-QA collection systems differ");
    import(entry,system,models);
  }
  if (!names.count(index.at("default_model"))) throw std::runtime_error("model-QA default is absent");
  if(system=="AuAu")for (const std::string kind : {"routing","depth"}) {
    const auto &group=index.at("comparisons").at(kind);
    if (!names.count(group.at("reference"))) throw std::runtime_error("model-QA comparison reference is absent");
    auto entryFor=[&](const Json &name)->const Json& {
      for(const auto &entry:index.at("models")) if(entry.at("name")==name)return entry;
      throw std::runtime_error("model-QA comparison member is absent");
    };
    const auto &reference=entryFor(group.at("reference"));
    for(const auto &name:group.at("models")) {
      const auto &entry=entryFor(name);
      if(entry.at("pins")!=reference.at("pins") || entry.at("et_ranges_gev")!=reference.at("et_ranges_gev"))
        throw std::runtime_error("model-QA comparison populations differ");
      Json a=recipes.at(reference.at("name")),b=recipes.at(entry.at("name"));
      for(const char *field:{"name","score_branch"}) {a.erase(field);b.erase(field);}
      if(kind=="routing") {a.erase("axes");b.erase("axes");}
      else {a.at("learner").erase("max_depth");b.at("learner").erase("max_depth");}
      if(a!=b)throw std::runtime_error("model-QA comparison changes additional settings");
    }
  }
  if(system=="PP") {
    const auto &groups=index.at("comparisons"), &paired=groups.at("ppg12"), &comparison=groups.at("auau_settings");
    const std::string pp= comparison.at("pp_model"), ppg12=paired.at("model");
    if(!names.count(pp) || !names.count(ppg12) || pp!=index.at("default_model").get<std::string>() ||
       ppg12.rfind("ppg12_retrained_",0)!=0 || paired.at("table")!="original_comparison/paired_scores" ||
       comparison.at("status")!="BOUND" || comparison.at("et_range_gev")!=Json::array({15.,35.}) ||
       comparison.at("auau_centrality_range")!=Json::array({50.,80.}))
      throw std::runtime_error("pp model-QA comparison is absent or unbound");
    const auto &ppRecipe=recipes.at(pp), &ppMeta=metadataByName.at(pp), &pairedMeta=metadataByName.at(ppg12);
    if(ppRecipe.at("et_ranges_gev").at("application")!=comparison.at("et_range_gev") ||
       ppMeta.at("heldout_scores").at("valid").get<int>()!=1 ||
       pairedMeta.at("original_comparison").at("reference_model_sha256")!=recipes.at(ppg12).at("reference_model_sha256"))
      throw std::runtime_error("pp saved comparison scope or frozen reference differs");
    const auto &reference= comparison.at("auau_reference");
    if(reference.at("name").get<std::string>().rfind("nominal_retrained_",0)!=0)
      throw std::runtime_error("pp comparison requires the canonical AuAu reference");
    const auto [auauRecipe,auauMeta]=import(reference,"AuAu",models->mkdir("references"),true);
    if(auauMeta.at("peripheral_15to35_scores").at("valid").get<int>()!=1)
      throw std::runtime_error("AuAu reference lacks valid peripheral held-out scores");
    for(const char *field:{"learner","training","weights"})
      if(auauRecipe.at(field)!=ppRecipe.at(field))throw std::runtime_error("pp/AuAu training settings differ");
    Json sharedFeatures=Json::array();
    for(const auto &feature:auauRecipe.at("features"))
      if(feature!="cluster_weta33_cogx" && feature!="cluster_wphi33_cogx" && feature!="centrality")sharedFeatures.push_back(feature);
    if(sharedFeatures!=ppRecipe.at("features"))throw std::runtime_error("pp/AuAu shared feature order differs");
    for(const char *field:{"et_min","et_max_exclusive","eta_abs_max","vertex_abs_max_cm"})
      if(auauRecipe.at("domain").at(field)!=ppRecipe.at("domain").at(field))throw std::runtime_error("pp/AuAu comparison domains differ");
    for(const std::string &key:std::vector<std::string>{pp+"/holdout_scores/roc",ppg12+"/original_comparison/paired_scores", "references/auau/peripheral_15to35_scores/roc"})
      if(!dynamic_cast<TTree*>(models->Get(key.c_str())))throw std::runtime_error("missing saved pp comparison table");
  }
  models->cd();TObjString(index.dump().c_str()).Write("collection");
  const Json validation={{"schema","PhotonIDModelQAImportV1"},{"status","PASS"},
    {"collection_sha256",index.at("collection_sha256")},{"scientific_acceptance",false}};
  TObjString(validation.dump().c_str()).Write("validation");output->cd();
}

int runModelQA(const std::string &collection, const std::string &output) {
  resetNativeInputs();
  namespace fs=std::filesystem;
  const auto input=fs::canonical(collection);
  const auto inputHash=photonjet_render::sha256(input.string(),true);
  const std::string partial=output+".partial";
  for(const auto &path:{output,partial,output+".render.json"})
    if(fs::exists(path) || fs::is_symlink(path)) throw std::runtime_error(path+" exists; choose a new output");
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(),"CREATE"));
  if(!out || out->IsZombie()) throw std::runtime_error("cannot write "+partial);
  savedModelQA(input.string(),out.get());
  if(photonjet_render::sha256(input.string(),true)!=inputHash)
    throw std::runtime_error("model-QA collection changed during import");
  // Empty measurement settings record the executable sources without implying
  // that any DATA, simulation histogram or measurement was processed.
  Json metadata=photonjet_render::verify(Json::parse(renderMetadata(Settings{})));
  metadata["model_qa"]={{"path",input.string()},{"sha256",inputHash}};
  const std::string render=photonjet_render::seal(metadata).dump();
  out->cd();
  TObjString("saved-model-qa").Write("analysis");
  const Json coverage={{"status","SAVED_MODEL_DIAGNOSTICS"},{"histogram_packages_processed",0},
                       {"scientific_acceptance",false}};
  TObjString(coverage.dump().c_str()).Write("input_coverage");
  TObjString(render.c_str()).Write("render_metadata");
  writeNativeProducerWitness(*out);
  requireExplicitHistogramErrors(out.get());
  out->Close();
  if(gSystem->Rename(partial.c_str(),output.c_str())!=0) throw std::runtime_error("cannot publish model QA output");
  renderProductReceipt(output,render);
  std::cout << "[FinalAnalysis] wrote " << output << std::endl;
  return 0;
}

/* Only the inputs that the measurement consumes are recomputed for a saved
   variation. All sums use the same order and one normalization factor as QA.
   Missing or statistically undefined fitted values remain absent/NaN. */
void writeMeasurementInputs(const Settings& s,TDirectory* results) {
  std::map<std::string,std::vector<std::string>> photonSamples;
  std::map<std::string,std::set<std::string>> classes;
  for(const auto& sys:s.systems) {
    auto& samples=photonSamples[sys.package];
    for(const auto& name:sys.simulation)
      if(std::find(samples.begin(),samples.end(),name)==samples.end())samples.push_back(name);
    classes[sys.package].insert(sys.centrality);
  }
  for(const auto& [package,photons]:photonSamples) {
    const std::string base=s.inputDirectory+"/"+package+"/";
    std::map<std::string,Sums> families;
    std::map<std::string,std::unique_ptr<TH2>> leakage;
    const std::map<std::string,std::vector<std::string>> samples={
      {"phosim",photons},{"jetsim",s.inclusive.count(package)?s.inclusive.at(package):std::vector<std::string>()}};
    for(const auto& [family,names]:samples)for(const auto& name:names) {
      auto input=openIfPresent(base+name+".root");
      if(!input)throw std::runtime_error("missing purity input package "+name);
      const double weight=weightOf(*input);
      for(const auto& c:classes.at(package)) {
        const std::string path="photons/"+c+"/"+(family=="phosim"?kSignalLabel:kBackgroundLabel)+"/abcd";
        auto* h=input->Get<TH2>(("qa/"+path).c_str());
        if(!h)throw std::runtime_error("missing physics purity histogram "+name+"/"+path);
        auto& sum=families[family][path];
        if(!sum) {
          auto* clone=static_cast<TH2*>(h->Clone());clone->SetDirectory(nullptr);clone->Scale(weight);sum.reset(clone);
        } else if(!static_cast<TH2*>(sum.get())->Add(h,weight))throw std::runtime_error("physics purity binning differs");
        if(family=="phosim") {
          auto* h=input->Get<TH2>(("response/"+c+"/"+s.jet+"/signal_leakage_truth_matched_abcd").c_str());
          if(!h)throw std::runtime_error("missing physics signal leakage "+name);
          auto& sum=leakage[c];
          if(!sum){sum.reset(static_cast<TH2*>(h->Clone()));sum->SetDirectory(nullptr);sum->Scale(weight);}
          else if(!sum->Add(h,weight))throw std::runtime_error("physics leakage binning differs");
        }
      }
    }
    Sums unusedData;
    for(const auto& sys:s.systems)if(sys.package==package)
      purityStudy(dirOf(results,"qa/"+package+"/purity_study/"+sys.name),sys,leakage.at(sys.centrality).get(),
                  families["phosim"],families["jetsim"],unusedData,false);
  }
}

int runQA(const Settings &s) {
  Phases phases("qa_" + std::filesystem::path(s.qaOutput).stem().string(),s.name,0,"all");
  try {
  phases.start("qa/input_validation");
  resetNativeInputs();
  const std::string coverage = inputCoverage(s, true);
  const std::string render = renderMetadata(s, true);
  const bool measurement=s.qaMode=="measurement_inputs";
  if(Json::parse(coverage).contains("measurement_replay") && !measurement)
    throw std::runtime_error("replayed packages require QA.Mode measurement_inputs; diagnostics remain nominal");
  if(measurement && !s.photonIDModels.empty())throw std::runtime_error("saved model diagnostics are nominal-only");
  if (!gSystem->AccessPathName(s.qaOutput.c_str())) throw std::runtime_error(s.qaOutput + " exists; choose a new output");
  const std::string partial = s.qaOutput + ".partial";
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot write " + partial);
  phases.start("qa/quantities");
  if(measurement)writeMeasurementInputs(s,out.get());else writeQA(s, out.get());
  if (!s.photonIDModels.empty()) savedModelQA(s.photonIDModels,out.get());
  phases.start("qa/write_and_verify");
  out->cd();
  TObjString(s.name.c_str()).Write("analysis");
  TObjString(s.nativeParameterExports.dump().c_str()).Write("native_parameter_exports");
  writeNativeProducerWitness(*out);
  TObjString(scientificParameterBindings(s).dump().c_str()).Write("scientific_parameter_bindings");
  TObjString(coverage.c_str()).Write("input_coverage");
  TObjString(render.c_str()).Write("render_metadata");
  if(measurement)TObjString("measurement_physics_inputs_only; diagnostic_QA=nominal_parent").Write("scope");
  requireExplicitHistogramErrors(out.get());
  out->Close();
  if (gSystem->Rename(partial.c_str(), s.qaOutput.c_str()) != 0) throw std::runtime_error("cannot move the QA file into place");
  renderProductReceipt(s.qaOutput, render);
  std::cout << "[FinalAnalysis] wrote " << s.qaOutput << std::endl;
  phases.complete();
  return 0;
  } catch(const std::exception& error) { phases.fail(error.what());throw; }
}

}  /* namespace fa */
