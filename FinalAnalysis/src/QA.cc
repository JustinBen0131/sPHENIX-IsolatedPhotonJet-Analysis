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
#include "RenderMetadata.h"
#include "../../TreeToHists/src/Configuration.h"

#include <algorithm>
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
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TKey.h"
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
   closure/<class>/       the purity closure on simulation: raw, corrected and truth purity of all simulated candidates
   purity_study/<result>/ on the result's photon bins, from the simulated candidates of both families: r_bkg, the
   background correlation bA bD / (bB bC) of the non-signal candidates (phoother, bkg);
   purity_scale, the purity non-closure (truth over leakage-corrected ABCD purity on 2 GeV
   bins, nonclosure_ratio, fitted with a line in photon pT weighted by the region-A candidates,
   nonclosure_line, and averaged over each photon bin with the region-A spectrum); the transfer
   terms of the xJ-level ABCD correction per photon and unit xJ (transfer_background: (bA/sA)
   (gA - gC)/(1 - k lC), transfer_leakage: k lC (fA - fC)/(1 - k lC), with f and g the recoil
   jets per signal and background candidate); factorisation_data and _simulation, R of the
   background in the score sub-bands 0.1-0.4, 0.4-0.6 and isolation sub-bands 6-12, 12-40 GeV.
   The analysis takes r_bkg and purity_scale for BackgroundCorrelation and PurityScale "qa".
   working_points/<what>/ the photon-ID score and isolation thresholds at 70, 80, 90% signal efficiency per 5% slice
   and photon-pT group, their per-slice medians and a line in centrality (checks of the packages)
   stitching/<family>_<kind>/  the stitched spectrum, its descriptive fit exp[a - (b + c ln x + d x) ln x] and RMS
   triggers/              the maximum-cluster-energy spectra and turn-ons of the runs with every study trigger enabled
   centrality_target      Au+Au: the DATA centrality of the package's trigger (the centrality weights' target)
   A package without a histogram, or with an empty one, leaves its figure and numbers out; nothing is invented. */

namespace fa {
namespace {

using Json = nlohmann::json;
const Edges kAbcdEdges = {15, 17, 19, 21, 23, 25, 35};

std::unique_ptr<TFile> openIfPresent(const std::string &path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) return nullptr;
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
      if (!sum) { auto *c = static_cast<TH1 *>(h->Clone()); c->SetDirectory(nullptr); c->Scale(w); sum.reset(c); }
      else if (!static_cast<TH1 *>(sum.get())->Add(h, w)) throw std::runtime_error("QA histograms disagree in binning: " + path);
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

TH1D *hist(const std::string &name, const Edges &edges, const std::vector<double> &v) {
  auto *h = new TH1D(name.c_str(), "", edges.size() - 1, edges.data());
  h->SetDirectory(nullptr);
  for (size_t i = 0; i < v.size(); ++i) h->SetBinContent(i + 1, v[i]);
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
  const double efficiencies[] = {0.7, 0.8, 0.9};
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
        if (std::isfinite(t)) ys.push_back(t);
      }
      perSlice[s] = median(ys);
    }
    Edges sliceEdges(17);
    for (int s = 0; s <= 16; ++s) sliceEdges[s] = 5 * s;
    put(out, hist(Form("median_%d", int(std::round(100 * eff))), sliceEdges, perSlice));
    grid->SetDirectory(nullptr);
    put(out, grid);
    /* the line in centrality over the slice centres (numpy.polyfit, degree 1) */
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    for (int s = 0; s < 16; ++s)
      if (std::isfinite(perSlice[s])) { const double x = 5 * s + 2.5; sx += x, sy += perSlice[s], sxx += x * x, sxy += x * perSlice[s], ++n; }
    const double slope = n >= 2 ? (n * sxy - sx * sy) / (n * sxx - sx * sx) : std::nan(""), intercept = n >= 2 ? (sy - slope * sx) / n : std::nan("");
    Edges two = {0, 1, 2};
    put(out, hist(Form("line_%d", int(std::round(100 * eff))), two, {intercept, slope}));
  }
}

/* The stitched spectrum of a family: per sample and summed, with the descriptive fit on [lo, hi]. A sample stitched by
   ownership (windows: its [low, high) of the spectrum's variable) enters the sum only inside its window; its own
   spectrum keeps its whole range, so neighbouring samples overlap at every boundary. */
using Windows = std::map<std::string, std::pair<double, double>>;
void stitching(TDirectory *out, const std::map<std::string, TH1D *> &parts, double lo, double hi, const Windows &windows = {}) {
  if (parts.empty()) return;
  TH1D *total = nullptr;
  for (const auto &[sample, h] : parts) {
    out->WriteTObject(h, ("sample_" + sample).c_str());
    std::unique_ptr<TH1D> owned(static_cast<TH1D *>(h->Clone()));
    owned->SetDirectory(nullptr);
    if (windows.count(sample))
      for (int b = 0; b <= owned->GetNbinsX() + 1; ++b) {
        const double a = owned->GetXaxis()->GetBinLowEdge(b), c = owned->GetXaxis()->GetBinUpEdge(b);
        if (b == 0 || b > owned->GetNbinsX() || a < windows.at(sample).first - 1e-9 || c > windows.at(sample).second + 1e-9)
          owned->SetBinContent(b, 0), owned->SetBinError(b, 0);
      }
    if (!total) { total = static_cast<TH1D *>(owned->Clone("total")); total->SetDirectory(nullptr); }
    else total->Add(owned.get());
  }
  std::vector<double> x, y, e;
  for (int b = 1; b <= total->GetNbinsX(); ++b) {
    const double a = total->GetXaxis()->GetBinLowEdge(b), c = total->GetXaxis()->GetBinUpEdge(b), w = c - a;
    if (a < lo - 1e-9 || c > hi + 1e-9) continue;
    const double yy = total->GetBinContent(b) / w, ee = total->GetBinError(b) / w;
    if (yy > 0 && std::isfinite(yy) && ee > 0) x.push_back(0.5 * (a + c)), y.push_back(yy), e.push_back(ee);
  }
  std::vector<double> coefficients(4, std::nan(""));
  if (x.size() >= 5) {
    TMatrixD a(x.size(), 4);
    TVectorD rhs(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      const double lx = std::log(x[i]), s = y[i] / e[i];
      a(i, 0) = s, a(i, 1) = -lx * s, a(i, 2) = -lx * lx * s, a(i, 3) = -x[i] * lx * s;
      rhs[i] = std::log(y[i]) * s;
    }
    const TVectorD c = leastSquares(a, rhs);
    for (int k = 0; k < 4; ++k) coefficients[k] = c[k];
  }
  /* the RMS of spectrum / fit - 1 over the fitted points */
  double squares = 0;
  int n = 0;
  for (size_t i = 0; i < x.size() && std::isfinite(coefficients[0]); ++i) {
    const double lx = std::log(x[i]);
    const double r = y[i] / std::exp(coefficients[0] - (coefficients[1] + coefficients[2] * lx + coefficients[3] * x[i]) * lx);
    if (std::isfinite(r)) squares += (r - 1) * (r - 1), ++n;
  }
  Edges four = {0, 1, 2, 3, 4};
  put(out, hist("fit", four, coefficients));
  put(out, new TParameter<double>("rms", n ? std::sqrt(squares / n) : std::nan("")));
  put(out, total);
}

using Sums = std::map<std::string, std::unique_ptr<TObject>>;
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
std::vector<std::vector<double>> leaders(const Labels &labels, const std::string &c, const Edges &edges, bool &complete) {
  std::vector<std::vector<double>> out(4, std::vector<double>(edges.size() - 1));
  for (const auto &label : labels) {
    auto *h = dynamic_cast<TH2 *>(sumOf(*label.first, "photons/" + c + "/" + label.second + "/abcd"));
    if (!filled(h)) { complete = false; continue; }
    for (int r = 0; r < 4; ++r) {
      const auto v = regroup(h, r, edges);
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

/* R of the background in two score and two isolation sub-bands (every ET): a d / (b c) and its statistical error. */
TH1D *factorisation(const std::string &name, const std::vector<TH1 *> &parts) {
  double n[4] = {0, 0, 0, 0};
  bool any = false;
  for (TH1 *part : parts) {
    auto *h = dynamic_cast<TH3 *>(part);
    if (!h) continue;
    any = true;
    for (int x = 1; x <= h->GetNbinsX(); ++x)
      for (int y = 1; y <= h->GetNbinsY(); ++y) {
        const double s0 = h->GetXaxis()->GetBinLowEdge(x), s1 = h->GetXaxis()->GetBinUpEdge(x);
        const double i0 = h->GetYaxis()->GetBinLowEdge(y), i1 = h->GetYaxis()->GetBinUpEdge(y);
        const int score = s0 >= 0.1 - 1e-9 && s1 <= 0.4 + 1e-9 ? 0 : s0 >= 0.4 - 1e-9 && s1 <= 0.6 + 1e-9 ? 1 : -1;
        const int iso = i0 >= 6 - 1e-9 && i1 <= 12 + 1e-9 ? 0 : i0 >= 12 - 1e-9 && i1 <= 40 + 1e-9 ? 1 : -1;
        if (score < 0 || iso < 0) continue;
        for (int z = 1; z <= h->GetNbinsZ(); ++z) n[2 * score + iso] += h->GetBinContent(x, y, z);
      }
  }
  if (!any) return nullptr;
  const double r = n[1] > 0 && n[2] > 0 ? n[0] * n[3] / (n[1] * n[2]) : std::nan("");
  const double error = std::min({n[0], n[1], n[2], n[3]}) > 0 ? r * std::sqrt(1 / n[0] + 1 / n[1] + 1 / n[2] + 1 / n[3]) : std::nan("");
  auto *h = new TH1D(name.c_str(), "score 0.1-0.4 and 0.4-0.6 x isolation 6-12 and 12-40 GeV", 6, 0, 6);
  h->SetDirectory(nullptr);
  const char *labels[] = {"a", "b", "c", "d", "R", "R_error"};
  const double values[] = {n[0], n[1], n[2], n[3], r, error};
  for (int b = 0; b < 6; ++b) h->GetXaxis()->SetBinLabel(b + 1, labels[b]), h->SetBinContent(b + 1, values[b]);
  return h;
}

/* The purity study of one result (see the header). Each quantity needs its candidate labels filled in both simulated
   families (the non-closure also the response's leakage); one that lacks them is left out. */
void purityStudy(TDirectory *out, const SystemSettings &sys, const TH2 *leakage, Sums &phosim, Sums &jetsim, Sums &data) {
  const std::string &c = sys.centrality;
  const Edges &pt = sys.photonEdges;
  const Labels signal = {{&phosim, "sig"}, {&jetsim, "jetprompt"}}, background = {{&phosim, "phoother"}, {&jetsim, "bkg"}};
  bool withBackground = true, withSignal = true;
  const auto b = leaders(background, c, pt, withBackground), s = leaders(signal, c, pt, withSignal);
  if (withBackground) {
    std::vector<double> r(pt.size() - 1, std::nan(""));
    for (size_t i = 0; i < r.size(); ++i)
      if (b[0][i] > 0 && b[1][i] > 0 && b[2][i] > 0 && b[3][i] > 0) r[i] = b[0][i] * b[3][i] / (b[1][i] * b[2][i]);
    put(out, hist("r_bkg", pt, r));
  }
  /* the purity non-closure on 2 GeV bins and its line in photon pT */
  if (leakage && withBackground && withSignal) {
    bool complete = true;
    const auto sf = leaders(signal, c, kFineEdges, complete), bf = leaders(background, c, kFineEdges, complete);
    std::vector<std::vector<double>> total(4, std::vector<double>(kFineEdges.size() - 1)), l(4);
    for (int q = 0; q < 4; ++q) {
      for (size_t i = 0; i < total[q].size(); ++i) total[q][i] = sf[q][i] + bf[q][i];
      l[q] = regroup(leakage, q, kFineEdges);
    }
    const std::vector<double> corrected = purityRows(total, l).second;
    std::vector<double> ratio(total[0].size(), std::nan(""));
    double w = 0, wx = 0, wy = 0, wxx = 0, wxy = 0;
    int used = 0;
    for (size_t i = 0; i < ratio.size(); ++i) {
      const double truth = total[0][i] > 0 ? sf[0][i] / total[0][i] : std::nan(""), x = 0.5 * (kFineEdges[i] + kFineEdges[i + 1]);
      ratio[i] = truth / corrected[i];
      if (!std::isfinite(ratio[i]) || !(total[0][i] > 0) || !(corrected[i] > 0.05) || !(truth > 0)) continue;
      w += total[0][i], wx += total[0][i] * x, wy += total[0][i] * ratio[i], wxx += total[0][i] * x * x, wxy += total[0][i] * x * ratio[i];
      ++used;
    }
    put(out, hist("nonclosure_ratio", kFineEdges, ratio));
    if (used >= 2) {
      const double slope = (w * wxy - wx * wy) / (w * wxx - wx * wx), intercept = (wy - slope * wx) / w;
      Edges two = {0, 1, 2};
      put(out, hist("nonclosure_line", two, {intercept, slope}));
      /* the line averaged over each photon bin with the region-A candidates of the native bins as weights */
      std::vector<double> scale(pt.size() - 1);
      for (size_t k = 0; k < scale.size(); ++k) {
        double sum = 0, weighted = 0;
        for (const auto &label : Labels{{&phosim, "sig"}, {&phosim, "phoother"}, {&jetsim, "bkg"}, {&jetsim, "jetprompt"}})
          if (auto *h = dynamic_cast<TH2 *>(sumOf(*label.first, "photons/" + c + "/" + label.second + "/abcd")))
            for (int y = 1; y <= h->GetNbinsY(); ++y) {
              const double centre = h->GetYaxis()->GetBinCenter(y);
              if (centre < pt[k] || centre >= pt[k + 1]) continue;
              sum += h->GetBinContent(1, y), weighted += h->GetBinContent(1, y) * (intercept + slope * centre);
            }
        scale[k] = sum > 0 ? weighted / sum : intercept + slope * 0.5 * (pt[k] + pt[k + 1]);
      }
      put(out, hist("purity_scale", pt, scale));
    }
  }
  /* the transfer terms of the xJ-level correction */
  Edges xj;
  bool complete = withBackground && withSignal;
  const auto fA = recoilJets(signal, c, 0, pt, xj, complete), fC = recoilJets(signal, c, 1, pt, xj, complete);
  const auto gA = recoilJets(background, c, 0, pt, xj, complete), gC = recoilJets(background, c, 1, pt, xj, complete);
  if (complete && !xj.empty()) {
    auto *terms = new TH2D("transfer_background", "(bA/sA)(gA - gC)/(1 - k lC) per signal photon and unit xJ", pt.size() - 1, pt.data(), xj.size() - 1, xj.data());
    auto *leak = new TH2D("transfer_leakage", "k lC (fA - fC)/(1 - k lC) per signal photon and unit xJ", pt.size() - 1, pt.data(), xj.size() - 1, xj.data());
    for (TH2D *h : {terms, leak}) h->SetDirectory(nullptr);
    for (size_t i = 0; i + 1 < pt.size(); ++i) {
      const double k = b[0][i] / b[2][i], lC = s[2][i] / s[0][i], denominator = 1 - k * lC;
      for (size_t x = 0; x + 1 < xj.size(); ++x) {
        const double width = xj[x + 1] - xj[x];
        terms->SetBinContent(i + 1, x + 1, (b[0][i] / s[0][i]) * (gA[i][x] / b[0][i] - gC[i][x] / b[2][i]) / denominator / width);
        leak->SetBinContent(i + 1, x + 1, k * lC * (fA[i][x] / s[0][i] - fC[i][x] / s[2][i]) / denominator / width);
      }
    }
    put(out, terms);
    put(out, leak);
  }
  /* the factorisation of the background: DATA and the simulated non-signal candidates */
  if (TH1D *h = factorisation("factorisation_data", {sumOf(data, "photons/" + c + "/data/scoreiso")})) put(out, h);
  if (TH1D *h = factorisation("factorisation_simulation", {sumOf(phosim, "photons/" + c + "/phoother/scoreiso"), sumOf(jetsim, "photons/" + c + "/bkg/scoreiso")}))
    put(out, h);
}

}  /* namespace */

/* Save QA quantities from the declared packages for later rendering. Missing QA inputs stay absent;
   the renderer reports unavailable figures rather than treating an empty plot as a measured zero. */
void writeQA(const Settings &s, TDirectory *results) {
  std::map<std::string, std::vector<std::string>> photonSamples;
  std::map<std::string, std::vector<std::string>> classes;
  for (const auto &sys : s.systems) {
    auto &v = photonSamples[sys.package];
    for (const auto &x : sys.simulation) if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
    classes[sys.package].push_back(sys.centrality);
  }
  for (const auto &[package, photons] : photonSamples) {
    const std::string base = s.inputDirectory + "/" + package + "/";
    auto data = openIfPresent(base + "data.root");
    if (!data || !data->Get("qa")) continue;  /* e.g. converted legacy packages: no QA fills */
    TDirectory *out = dirOf(results, "qa/" + package);
    Sums dataSums;
    sumTree(static_cast<TDirectory *>(data->Get("qa")), "", 1.0, dataSums);
    writeTree(out->mkdir("data"), dataSums);
    const std::map<std::string, std::vector<std::string>> families = {{"phosim", photons}, {"jetsim", s.inclusive.count(package) ? s.inclusive.at(package) : std::vector<std::string>()}};
    std::map<std::string, Sums> familySums;
    std::map<std::string, std::map<std::string, std::map<std::string, TH1D *>>> generator, reco;  /* family -> class -> sample */
    std::map<std::string, std::unique_ptr<TH2>> leakageNative;  /* per class: the simulated signal leakage (A-D x reco ET) */
    Windows windows;  /* per sample stitched by ownership: its window of the generator spectrum */
    for (const auto &[family, samples] : families)
      for (const auto &sample : samples) {
        auto f = openIfPresent(base + sample + ".root");
        if (!f || !f->Get("qa")) continue;
        const double w = weightOf(*f);
        /* a sample stitched by ownership (Au+Au embedded simulation): its window of the generator spectrum */
        const Json normalization = Json::parse(static_cast<TObjString *>(f->Get("meta/package"))->GetString().Data()).value("normalization", Json::object());
        if (normalization.contains("ownership")) {
          const Json &window = normalization.at("ownership").at("window_gev");
          windows[sample] = {window[0].get<double>(), window[1].is_null() ? std::numeric_limits<double>::infinity() : window[1].get<double>()};
        }
        sumTree(static_cast<TDirectory *>(f->Get("qa")), "", w, familySums[family]);
        if (auto *g = dynamic_cast<TH2 *>(f->Get("qa/generator/leading_pt"))) {
          auto *p = g->ProjectionX(("gen_" + sample).c_str(), family == "phosim" ? 1 : 2, family == "phosim" ? 1 : 2, "e");
          p->SetDirectory(nullptr), p->Scale(w);
          generator[family]["all"][sample] = p;
        }
        for (const auto &c : classes[package]) {
          TH1D *regionA = nullptr;
          for (const char *label : family == "phosim" ? std::vector<const char *>{"sig", "phoother"} : std::vector<const char *>{"bkg", "jetprompt"})
            if (auto *a = dynamic_cast<TH2 *>(f->Get(("qa/photons/" + c + "/" + label + "/abcd").c_str()))) {
              auto *p = a->ProjectionY(("reco_" + sample + label).c_str(), 1, 1, "e");
              p->SetDirectory(nullptr);
              if (!regionA) regionA = p; else regionA->Add(p), delete p;
            }
          if (regionA) regionA->Scale(w), reco[family][c][sample] = regionA;
          if (family == "phosim")
            if (auto *l = dynamic_cast<TH2 *>(f->Get(("response/" + c + "/" + s.jet + "/signal_leakage_truth_matched_abcd").c_str()))) {
              auto &L = leakageNative[c];
              if (!L) L.reset(static_cast<TH2 *>(l->Clone(("leakage_" + c).c_str()))), L->SetDirectory(nullptr), L->Scale(w);
              else L->Add(l, w);
            }
        }
      }
    for (auto &[family, sums] : familySums) writeTree(out->mkdir(family.c_str()), sums);
    /* stitching: generator truth and the region-A tag, per family */
    for (const auto &[family, byClass] : generator)
      for (const auto &[c, parts] : byClass) stitching(dirOf(out, "stitching/" + family + "_generator"), parts, family == "phosim" ? 15 : 10, family == "phosim" ? 40 : 60, windows);
    const std::string central = *std::min_element(classes[package].begin(), classes[package].end(), [](const std::string &a, const std::string &b) {
      return (a.rfind("cent", 0) == 0 ? std::atof(a.c_str() + 4) : 0) < (b.rfind("cent", 0) == 0 ? std::atof(b.c_str() + 4) : 0);
    });
    for (const auto &[family, byClass] : reco)
      if (byClass.count(central)) stitching(dirOf(out, "stitching/" + family + "_reco"), byClass.at(central), 15, 35);
    /* DATA purity and the simulation closure, per class, on the ABCD bins */
    std::map<std::string, std::vector<std::vector<double>>> leakageByClass;
    for (const auto &[c, h] : leakageNative)
      for (int r = 0; r < 4; ++r) leakageByClass[c].push_back(regroup(h.get(), r, kAbcdEdges));
    for (const auto &c : classes[package]) {
      auto *abcd = dynamic_cast<TH2 *>(data->Get(("measurement/" + c + "/abcd").c_str()));
      if (abcd) {
        auto *counts = new TH2D(c.c_str(), "event-leading candidates;region (A,B,C,D);photon E_{T} [GeV]", 4, 0, 4, kAbcdEdges.size() - 1, kAbcdEdges.data());
        for (int r = 0; r < 4; ++r) {
          const auto v = regroup(abcd, r, kAbcdEdges), e2 = regroup(abcd, r, kAbcdEdges, true);
          for (size_t i = 0; i < v.size(); ++i) counts->SetBinContent(r + 1, i + 1, v[i]), counts->SetBinError(r + 1, i + 1, std::sqrt(e2[i]));
        }
        counts->SetDirectory(nullptr);
        put(dirOf(out, "abcd"), counts);
      }
      if (abcd && leakageByClass.count(c)) {
        std::vector<std::vector<double>> counts(4);
        for (int r = 0; r < 4; ++r) counts[r] = regroup(abcd, r, kAbcdEdges);
        const auto [raw, corrected] = purityRows(counts, leakageByClass[c]);
        TDirectory *d = dirOf(out, "purity/" + c);
        put(d, hist("raw", kAbcdEdges, raw));
        put(d, hist("corrected", kAbcdEdges, corrected));
        for (int r = 1; r < 4; ++r) {
          std::vector<double> fraction(kAbcdEdges.size() - 1, std::nan(""));
          for (size_t i = 0; i < fraction.size(); ++i) if (leakageByClass[c][0][i] > 0) fraction[i] = leakageByClass[c][r][i] / leakageByClass[c][0][i];
          put(d, hist(Form("leakage_%c", "ABCD"[r]), kAbcdEdges, fraction));
        }
      }
      /* closure: every simulated candidate (both families, all labels) as pseudo-data */
      std::vector<std::vector<double>> total(4, std::vector<double>(kAbcdEdges.size() - 1)), prompt = total;
      bool complete = leakageByClass.count(c) && familySums.count("phosim") && familySums.count("jetsim");
      for (const auto &[family, labels] : std::map<std::string, std::vector<std::string>>{{"phosim", {"sig", "phoother"}}, {"jetsim", {"bkg", "jetprompt"}}})
        for (const auto &label : labels) {
          auto it = familySums[family].find("photons/" + c + "/" + label + "/abcd");
          if (it == familySums[family].end() || !filled(static_cast<TH1 *>(it->second.get()))) { complete = false; continue; }
          auto *h = static_cast<TH2 *>(it->second.get());
          for (int r = 0; r < 4; ++r) {
            const auto v = regroup(h, r, kAbcdEdges);
            for (size_t i = 0; i < v.size(); ++i) {
              total[r][i] += v[i];
              if (label == "sig" || label == "jetprompt") prompt[r][i] += v[i];
            }
          }
        }
      if (complete) {
        const auto [raw, corrected] = purityRows(total, leakageByClass[c]);
        std::vector<double> truth(total[0].size(), std::nan(""));
        for (size_t i = 0; i < truth.size(); ++i) if (total[0][i] > 0) truth[i] = prompt[0][i] / total[0][i];
        TDirectory *d = dirOf(out, "closure/" + c);
        put(d, hist("raw", kAbcdEdges, raw));
        put(d, hist("corrected", kAbcdEdges, corrected));
        put(d, hist("truth", kAbcdEdges, truth));
      }
    }
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
    /* trigger study and the centrality target: per DATA source, with its own trigger table */
    auto *sources = dynamic_cast<TTree *>(data->Get("records/sources"));
    auto *cells = dynamic_cast<TTree *>(data->Get("records/source_contributions"));
    auto *bitsRaw = dynamic_cast<TH3 *>(data->Get("qa/events/trigger_bits_raw"));
    const bool perRunStudy=data->Get("qa/events/trigger_study_bits_raw")!=nullptr;
    if (perRunStudy) bitsRaw=dynamic_cast<TH3 *>(data->Get("qa/events/trigger_study_bits_raw"));
    const std::string studyCells=perRunStudy ? "qa/events/trigger_study_bits_corrected" : "qa/events/trigger_bits_raw";
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
      TTreeReader rs(sources), rc(cells);
      TTreeReaderValue<std::string> json(rs, "json"), histogram(rc, "histogram");
      TTreeReaderValue<Long64_t> begin(rs, "contribution_begin"), end(rs, "contribution_end"), bin(rc, "bin");
      TTreeReaderValue<double> value(rc, "value"), variance(rc,"variance");
      std::unique_ptr<TH3> one(static_cast<TH3 *>(bitsRaw->Clone("one_source")));
      one->SetDirectory(nullptr);
      std::vector<double> target;
      std::unique_ptr<TH2> centralityOne;
      auto *centralityBits = dynamic_cast<TH2 *>(data->Get("qa/events/centrality_scaled_bits"));
      if (auau && centralityBits) {
        centralityOne.reset(static_cast<TH2 *>(centralityBits->Clone("centrality_one")));
        centralityOne->SetDirectory(nullptr);
        target.assign(centralityBits->GetNbinsX(), 0.0);
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
          else if (centralityOne && *histogram == "qa/events/centrality_scaled_bits") centralityOne->SetBinContent(static_cast<int>(*bin), *value);
        }
        /* the centrality target: the selected trigger, enabled in this source */
        if (centralityOne) {
          int b = -1;
          for (size_t i = 0; i < info.at("name").size(); ++i) if (info.at("name")[i].get<std::string>() == trigger) b = i;
          if (b >= 0 && info.at("initial_prescale")[b].get<int>() >= 0)
            for (int x = 1; x <= centralityOne->GetNbinsX(); ++x) target[x - 1] += centralityOne->GetBinContent(x, b + 1);
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
        studiedRuns.insert(run); runs=studiedRuns.size(); firstSeen = firstSeen ? std::min(firstSeen, run) : run; lastSeen = std::max(lastSeen, run);
        for (int x = 1; x <= nE; ++x) {
          refSpectrum[x - 1] += one->GetBinContent(x, bits[0] + 1, 65);
          den[x - 1] += one->GetBinContent(x, bits[0] + 1, 65);
          const double referenceError=one->GetBinError(x,bits[0]+1,65);
          refVariance[x-1]+=referenceError*referenceError; denVariance[x-1]+=referenceError*referenceError;
          for (size_t p = 0; p < photonsStudy.size(); ++p) {
            photonSpectrum[photonsStudy[p].first][x - 1] += one->GetBinContent(x, bits[p + 1] + 1, 65);
            const int y=perRunStudy ? bits[p+1]+1 : bits[0]+1, z=perRunStudy ? 65 : bits[p+1]+1;
            num[photonsStudy[p].first][x - 1] += one->GetBinContent(x,y,z);
            const double e=one->GetBinError(x,bits[p+1]+1,65), en=one->GetBinError(x,y,z);
            photonVariance[photonsStudy[p].first][x-1]+=e*e; numVariance[photonsStudy[p].first][x-1]+=en*en;
          }
        }
      }
      Edges energy;
      for (int b = 1; b <= nE + 1; ++b) energy.push_back(bitsRaw->GetXaxis()->GetBinLowEdge(b));
      TDirectory *d = dirOf(out, "triggers");
      auto studyHist=[&](const std::string &name,const std::vector<double> &v,const std::vector<double> &e) {
        auto *h=hist(name,energy,v);
        for (int b=1;b<=h->GetNbinsX();++b) h->SetBinError(b,std::sqrt(e[b-1]));
        return h;
      };
      put(d, studyHist("reference",refSpectrum,refVariance));
      put(d, studyHist("reference_denominator",den,denVariance));
      for (const auto &p : photonsStudy) {
        put(d, studyHist("spectrum_"+p.first,photonSpectrum[p.first],photonVariance[p.first]));
        put(d, studyHist("numerator_"+p.first,num[p.first],numVariance[p.first]));
      }
      put(d, new TObjString(Form("runs=%d;first=%d;last=%d;reference=%s;correction=%s", runs, firstSeen, lastSeen, reference.c_str(),perRunStudy ? "run_live_over_scaled" : "legacy_raw")), "study");
      if (centralityBits) {
        Edges centrality;
        for (int b = 1; b <= centralityBits->GetNbinsX() + 1; ++b) centrality.push_back(centralityBits->GetXaxis()->GetBinLowEdge(b));
        put(out, hist("centrality_target", centrality, target));
      }
    }
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
          clone->Write(key->GetName());clone->SetDirectory(nullptr);delete clone;
        }
        else {to->cd();object->Write(key->GetName());}
      }
    };
    copy(&input,dest);
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
  out->Close();
  if(gSystem->Rename(partial.c_str(),output.c_str())!=0) throw std::runtime_error("cannot publish model QA output");
  renderProductReceipt(output,render);
  std::cout << "[FinalAnalysis] wrote " << output << std::endl;
  return 0;
}

int runQA(const Settings &s) {
  const std::string coverage = inputCoverage(s);
  const std::string render = renderMetadata(s);
  if (!gSystem->AccessPathName(s.qaOutput.c_str())) throw std::runtime_error(s.qaOutput + " exists; choose a new output");
  const std::string partial = s.qaOutput + ".partial";
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot write " + partial);
  writeQA(s, out.get());
  if (!s.photonIDModels.empty()) savedModelQA(s.photonIDModels,out.get());
  out->cd();
  TObjString(s.name.c_str()).Write("analysis");
  TObjString(coverage.c_str()).Write("input_coverage");
  TObjString(render.c_str()).Write("render_metadata");
  out->Close();
  if (gSystem->Rename(partial.c_str(), s.qaOutput.c_str()) != 0) throw std::runtime_error("cannot move the QA file into place");
  renderProductReceipt(s.qaOutput, render);
  std::cout << "[FinalAnalysis] wrote " << s.qaOutput << std::endl;
  return 0;
}

}  /* namespace fa */
