/* Analysis: the ordered chain of every configured result, its data toys, and results.root. */

/* final_analysis --config config/analysis.cfg */

/* For each result: measured data and response -> purity and combinatoric correction -> photon-yield unfolding and the
   joint (photon pT, xJ) unfolding from the data-fitted prior -> per-photon xJ. Data toys (independent Poisson counts of
   the measured histograms, or a Poisson(1) multiplier per event of the event records), corrected like the data (a toy
   whose ABCD solution changes branch is rejected), give the covariance of the corrected pairs that weights the prior
   fit; every toy is refitted and unfolded, and the statistical uncertainty is the central 68% half width of the toy
   results. results.root holds every quantity the figures draw, so drawing never recomputes. */
#include "FinalAnalysis.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TSystem.h"

namespace fa {
namespace {

TH1D *hist1(const char *name, const Edges &edges, const TVectorD &values, const TVectorD *variance = nullptr) {
  auto *h = new TH1D(name, "", edges.size() - 1, edges.data());
  for (int i = 0; i < values.GetNrows(); ++i) {
    h->SetBinContent(i + 1, values[i]);
    h->SetBinError(i + 1, variance ? std::sqrt(std::max((*variance)[i], 0.0)) : 0.0);
  }
  return h;
}

TH2D *hist2(const char *name, const Edges &x, const Edges &y, const TMatrixD &values, const TMatrixD *variance = nullptr) {
  auto *h = new TH2D(name, "", x.size() - 1, x.data(), y.size() - 1, y.data());
  for (int i = 0; i < values.GetNrows(); ++i)
    for (int j = 0; j < values.GetNcols(); ++j) {
      h->SetBinContent(i + 1, j + 1, values(i, j));
      h->SetBinError(i + 1, j + 1, variance ? std::sqrt(std::max((*variance)(i, j), 0.0)) : 0.0);
    }
  return h;
}

Edges indexEdges(int n) {
  Edges e(n + 1);
  for (int i = 0; i <= n; ++i) e[i] = i;
  return e;
}

TMatrixD cellsToMatrix(const TVectorD &v, int nPt, int nX) {
  TMatrixD m(nPt, nX);
  for (int i = 0; i < nPt; ++i)
    for (int x = 0; x < nX; ++x) m(i, x) = v[i * nX + x];
  return m;
}

TVectorD matrixToCells(const TMatrixD &m) {
  TVectorD v(m.GetNrows() * m.GetNcols());
  for (int i = 0; i < m.GetNrows(); ++i)
    for (int x = 0; x < m.GetNcols(); ++x) v[i * m.GetNcols() + x] = m(i, x);
  return v;
}

TVectorD columnSums(const TMatrixD &m) {
  TVectorD out(m.GetNcols());
  for (int j = 0; j < m.GetNcols(); ++j)
    for (int i = 0; i < m.GetNrows(); ++i) out[j] += m(i, j);
  return out;
}

TVectorD ones(int n) {
  TVectorD v(n);
  for (int i = 0; i < n; ++i) v[i] = 1;
  return v;
}

/* numpy.quantile with linear interpolation */
double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double at = q * (v.size() - 1);
  const size_t i = static_cast<size_t>(std::floor(at));
  return i + 1 < v.size() ? v[i] + (at - i) * (v[i + 1] - v[i]) : v.back();
}

std::string join(const TVectorD &v) {
  std::ostringstream out;
  out.precision(10);
  for (int i = 0; i < v.GetNrows(); ++i) out << (i ? "," : "") << v[i];
  return out.str();
}
std::string join(const std::vector<double> &v) { return v.empty() ? "1" : join(TVectorD(v.size(), v.data())); }

/* BackgroundCorrelation or PurityScale "qa": the values the qa step derived from the simulation for this result
   (qa/<package>/purity_study/<result>/ of the QA file). An implausible purity scale (outside 0.5-1.5) is refused. */
void resolveFromQA(const Settings &s, SystemSettings &sys) {
  if (!sys.backgroundCorrelationFromQA && !sys.purityScaleFromQA) return;
  std::unique_ptr<TFile> f(TFile::Open(s.qaOutput.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error(sys.name + ": the qa values need the QA file " + s.qaOutput + " (final_analysis.sh qa)");
  auto values = [&](const std::string &name) {
    const std::string key = "qa/" + sys.package + "/purity_study/" + sys.name + "/" + name;
    auto *h = dynamic_cast<TH1 *>(f->Get(key.c_str()));
    if (!h || h->GetNbinsX() + 1 != static_cast<int>(sys.photonEdges.size())) throw std::runtime_error(s.qaOutput + ": no " + key + " on the photon bins");
    std::vector<double> v;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      if (std::fabs(h->GetXaxis()->GetBinLowEdge(b) - sys.photonEdges[b - 1]) > 1e-9) throw std::runtime_error(key + ": other photon bins");
      v.push_back(h->GetBinContent(b));
      if (!(v.back() > 0) || !std::isfinite(v.back())) throw std::runtime_error(key + ": no positive value in photon bin " + std::to_string(b - 1));
    }
    return v;
  };
  if (sys.backgroundCorrelationFromQA) sys.backgroundCorrelation = values("r_bkg");
  if (sys.purityScaleFromQA) {
    sys.purityScale = values("purity_scale");
    for (double v : sys.purityScale)
      if (v <= 0.5 || v >= 1.5) throw std::runtime_error(sys.name + ": the purity non-closure scale " + std::to_string(v) + " is outside 0.5-1.5; review it before use");
  }
}

/* A Poisson toy of the measured histograms; its variance is its count. */
Measured poissonToy(const Measured &m, TRandom3 &rng) {
  Measured t = m;
  const std::pair<TMatrixD *, TMatrixD *> parts[] = {
      {&t.abcd, &t.abcdVariance}, {&t.pairsA, &t.pairsAVariance}, {&t.pairsC, &t.pairsCVariance}};
  for (const auto &[x, v] : parts)
    for (int i = 0; i < x->GetNrows(); ++i)
      for (int j = 0; j < x->GetNcols(); ++j) (*x)(i, j) = (*v)(i, j) = rng.Poisson((*x)(i, j));
  return t;
}

/* Corrected photons and pairs of the accepted data toys. */
struct Toys {
  std::vector<TVectorD> photons, pairs;
  int exceptions = 0, branchChanges = 0;
};

template <class Draw> Toys dataToys(Draw draw, const Window &w, const SystemSettings &sys, const Corrected &central, int n) {
  Toys t;
  for (int i = 0; i < n; ++i) {
    try {
      const Corrected c = correct(draw(), w, sys);
      if (c.physical != central.physical) { ++t.branchChanges; continue; }
      t.photons.push_back(c.photons);
      t.pairs.push_back(matrixToCells(c.pairs));
    } catch (const std::exception &) {
      ++t.exceptions;
    }
  }
  return t;
}

TMatrixD covariance(const std::vector<TVectorD> &rows);

/* The covariance of the per-photon results of a toy ensemble, without the toys that lie beyond five robust sigma
   (1.4826 MAD) of the ensemble median in any reported bin: a result-independent trim of the rare toys in which the
   iterative unfolding of signed inputs diverges (about one toy in a million for a Gaussian ensemble). */
TMatrixD robustCovariance(const std::vector<TVectorD> &rows, const std::vector<bool> &reported) {
  const int n = rows.front().GetNrows();
  std::vector<double> median(n), mad(n);
  for (int x = 0; x < n; ++x) {
    std::vector<double> v;
    for (const auto &r : rows) v.push_back(r[x]);
    median[x] = quantile(v, 0.5);
    for (auto &e : v) e = std::fabs(e - median[x]);
    mad[x] = 1.4826 * quantile(v, 0.5);
  }
  std::vector<TVectorD> kept;
  for (const auto &r : rows) {
    bool keep = true;
    for (int x = 0; x < n; ++x) keep = keep && (!reported[x] || std::fabs(r[x] - median[x]) <= 5.0 * std::max(mad[x], 1e-300));
    if (keep) kept.push_back(r);
  }
  return kept.size() >= 2 ? covariance(kept) : TMatrixD(n, n);
}

TMatrixD covariance(const std::vector<TVectorD> &rows) {  /* sample covariance, n - 1 */
  const int n = rows.size(), d = rows[0].GetNrows();
  TVectorD mean(d);
  for (const auto &r : rows) mean += r;
  mean *= 1.0 / n;
  TMatrixD cov(d, d);
  for (const auto &r : rows)
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j) cov(i, j) += (r[i] - mean[i]) * (r[j] - mean[j]);
  cov *= 1.0 / (n - 1);
  return cov;
}

/* Bootstrap multiplicities of the simulated sources: n draws with replacement within each sample (sources ordered by sample). */
std::vector<double> occurrences(const std::vector<SourceWindow> &sources, TRandom3 &draw) {
  std::vector<int> count;
  for (const auto &src : sources) {
    if (src.sample >= static_cast<int>(count.size())) count.resize(src.sample + 1);
    ++count[src.sample];
  }
  std::vector<double> multiplicity(sources.size());
  for (size_t k = 0, first = 0; k < count.size(); first += count[k], ++k)
    for (int d = 0; d < count[k]; ++d) multiplicity[first + draw.Integer(count[k])] += 1;
  return multiplicity;
}

/* Per-photon views of pairs (photon bin x xJ cells) and photons: every photon bin, then their integral (rows) x xJ. */
TMatrixD views(const TVectorD &pairs, const TVectorD &photons, const Edges &xj) {
  const int nX = xj.size() - 1, nPt = photons.GetNrows();
  TMatrixD v(nPt + 1, nX);
  const double all = photons.Sum();
  for (int x = 0; x < nX; ++x) {
    const double width = xj[x + 1] - xj[x];
    double sum = 0;
    for (int i = 0; i < nPt; ++i) v(i, x) = pairs[i * nX + x] / photons[i] / width, sum += pairs[i * nX + x];
    v(nPt, x) = sum / all / width;
  }
  return v;
}

/* One mode of the estimator closure: the central unfolded and truth views, the spread of their difference over source
   resamples and the pulls; per view row the largest |pull| and the mean pull^2 over the reported xJ bins. */
struct ClosureMode {
  TMatrixD unfolded, truth, sigma, pull;
  int toys = 0, failed = 0;
  std::vector<double> maxPull, chi2;
};

/* Prior-fit and unfolding closure, after the correction stage.

   Pseudo-data are response-column sums plus detector and boundary fakes.
   They do not exercise the ABCD or data-driven combinatoric correction.
   The training response supplies the prior fit and both unfoldings, with
   the data covariance converted to simulation units. Full uses one shared
   source resample; A_to_B and B_to_A use source-disjoint halves, ranked by
   source hash within each sample and resampled independently. Uncertainty
   is the spread of unfolded minus truth across those resamples. */
std::map<std::string, ClosureMode> estimatorClosure(const Settings &s, const SystemSettings &sys, const std::vector<SourceWindow> &sources,
                                                    const TMatrixD &dataCovariance, const TVectorD &dataPhotons, const Window &w) {
  const Edges &xj = s.xjEdges;
  const int nX = xj.size() - 1, nPt = sys.photonEdges.size() - 1;
  std::map<int, std::vector<size_t>> bySample;
  for (size_t i = 0; i < sources.size(); ++i) bySample[sources[i].sample].push_back(i);
  std::vector<double> all(sources.size(), 1.0), a(sources.size()), b(sources.size());
  for (const auto &entry : bySample) {
    std::vector<size_t> ranked = entry.second;
    std::sort(ranked.begin(), ranked.end(), [&](size_t p, size_t q) { return sources[p].identity < sources[q].identity; });
    for (size_t r = 0; r < ranked.size(); ++r) (r % 2 ? b : a)[ranked[r]] = 1;
  }
  TMatrixD covariance(dataCovariance);
  for (int i = 0; i < covariance.GetNrows(); ++i)
    for (int j = 0; j < covariance.GetNcols(); ++j)
      covariance(i, j) /= (dataPhotons[i / nX] / w.simPhotons[i / nX]) * (dataPhotons[j / nX] / w.simPhotons[j / nX]);
  auto onePass = [&](const std::vector<double> &train, const std::vector<double> &test) {
    const Window tr = combine(sources, train), te = combine(sources, test);
    TVectorD pairs = columnSums(te.pairResponse), photons = columnSums(te.photonResponse);
    pairs += te.pairFakes, photons += te.photonFakes;
    TVectorD weights = ones(nX);
    if (s.priorForm != "simulation") {
      const PriorFit f = fitPrior(tr, pairs, photons, covariance, s, sys);
      if (!f.ok) throw std::runtime_error("the prior fit did not converge");
      for (int x = 0; x < nX; ++x) weights[x] = std::pow(f.weights[x], s.priorPower);
    }
    return std::make_pair(views(unfoldPairsWithPrior(pairs, tr, sys.pairIterations, weights, s.fakeModel),
                                unfoldPhotons(photons, tr, sys.photonIterations), xj),
                          views(te.pairTruth, te.photonTruth, xj));
  };
  TRandom3 draw(s.closureSeed);
  auto resample = [&](const std::vector<double> &base) {  /* n draws with replacement among a sample's sources in the base */
    std::vector<double> m(sources.size());
    for (const auto &entry : bySample) {
      std::vector<size_t> in;
      for (size_t i : entry.second) if (base[i] > 0) in.push_back(i);
      for (size_t d = 0; d < in.size(); ++d) m[in[draw.Integer(in.size())]] += 1;
    }
    return m;
  };
  std::map<std::string, ClosureMode> out;
  const std::vector<std::tuple<std::string, const std::vector<double> *, const std::vector<double> *>> modes = {
      {"full", &all, &all}, {"A_to_B", &a, &b}, {"B_to_A", &b, &a}};
  for (const auto &[name, train, test] : modes) {
    ClosureMode m;
    const auto central = onePass(*train, *test);
    m.unfolded.ResizeTo(central.first), m.unfolded = central.first;
    m.truth.ResizeTo(central.second), m.truth = central.second;
    std::vector<TMatrixD> differences;
    for (int k = 0; k < s.closureToys; ++k) {
      const std::vector<double> tr = resample(*train), te = name == "full" ? tr : resample(*test);
      try {
        auto [u, t] = onePass(tr, te);
        differences.push_back(u - t);
      } catch (const std::exception &) {  /* an unconverged fit or a resample the unfolding cannot take */
        ++m.failed;
      }
    }
    m.toys = differences.size();
    m.sigma.ResizeTo(nPt + 1, nX), m.pull.ResizeTo(nPt + 1, nX);
    for (int r = 0; r <= nPt; ++r) {
      double largest = 0, squares = 0;
      int reported = 0;
      for (int x = 0; x < nX; ++x) {
        double mean = 0, sum2 = 0;
        for (const auto &d : differences) mean += d(r, x);
        mean /= differences.size();
        for (const auto &d : differences) sum2 += (d(r, x) - mean) * (d(r, x) - mean);
        m.sigma(r, x) = differences.size() >= 2 ? std::sqrt(sum2 / (differences.size() - 1)) : std::nan("");
        /* a spread at rounding level (every resample the same, e.g. a half of one source) gives no pull */
        const double rounding = 1e-12 * (std::fabs(m.unfolded(r, x)) + std::fabs(m.truth(r, x)));
        m.pull(r, x) = m.sigma(r, x) > rounding ? (m.unfolded(r, x) - m.truth(r, x)) / m.sigma(r, x) : 0.0;
        if (xj[x] >= sys.reportedLow - 1e-9 && xj[x + 1] <= s.reportedHigh + 1e-9)
          largest = std::max(largest, std::fabs(m.pull(r, x))), squares += m.pull(r, x) * m.pull(r, x), ++reported;
      }
      m.maxPull.push_back(largest), m.chi2.push_back(reported ? squares / reported : 0.0);
    }
    out.emplace(name, m);
  }
  return out;
}

/* The rows of an iteration scan: k along x, the metrics along y. */
TH2D *scanHistogram(const char *name, const IterationScan &scan, int selected) {
  const char *metrics[] = {"Q", "S_data", "S_response", "D_k", "refold chi2", "MC closure", "MC refold", "negative fraction",
                           "data toy success", "response toy success", "passes gates", "selected"};
  auto *h = new TH2D(name, "", 12, 0.5, 12.5, 12, 0, 12);
  for (int y = 0; y < 12; ++y) h->GetYaxis()->SetBinLabel(y + 1, metrics[y]);
  for (const auto &r : scan.rows) {
    const double v[] = {r.objective, r.dataStat, r.responseStat, r.deviation, r.refold, r.mcClosure, r.mcRefold, r.negativeFraction,
                        r.dataSuccess, r.responseSuccess, r.passes ? 1.0 : 0.0, r.k == selected ? 1.0 : 0.0};
    for (int y = 0; y < 12; ++y) h->SetBinContent(r.k, y + 1, v[y]);
  }
  return h;
}

void runSystem(const Settings &s, SystemSettings sys, TDirectory *out, TRandom3 &rng) {
  resolveFromQA(s, sys);
  const bool events = s.dataToyModel == "event_bootstrap", dataDriven = s.combinatoricMethod == "data_dphi";
  const EventTable table = events ? readEvents(s, sys) : EventTable();
  const Measured m = events ? fromEvents(table, nullptr, s) : readMeasured(s, sys);
  const Window w = readWindow(s, sys);
  const Corrected c = correct(m, w, sys);
  /* The result uses the physical root of the leakage-aware ABCD equation in every photon bin; the clamped estimate
     solveLeakageABCD falls back to is never a result. Toys of the other branch are rejected below. */
  for (size_t i = 0; i < c.physical.size(); ++i)
    if (!c.physical[i])
      throw std::runtime_error(sys.name + ": no physical leakage-corrected ABCD solution in photon bin " + std::to_string(i));
  /* for the figures: the recoil jets and their correction before any combinatoric subtraction */
  Measured raw = m;
  if (dataDriven) { raw.pairsA += m.backgroundA; raw.pairsC += m.backgroundC; }
  const Corrected before = dataDriven ? correct(raw, w, sys) : c;
  TMatrixD combinatoricShown(c.combinatoric);
  if (dataDriven) { combinatoricShown = before.pairs; combinatoricShown -= c.pairs; }
  const Edges &pt = sys.photonEdges, &xj = s.xjEdges;
  const int nPt = pt.size() - 1, nX = xj.size() - 1;
  const Edges regions = indexEdges(4), cellEdges = indexEdges(nPt * nX);
  const TVectorD cells = matrixToCells(c.pairs);
  std::vector<double> multipliers(table.weight.size());
  auto dataToy = [&](TRandom3 &random) {  /* a data toy: Poisson counts, or one Poisson(1) multiplier per event */
    if (!events) return poissonToy(m, random);
    for (auto &x : multipliers) x = random.Poisson(1.0);
    return fromEvents(table, &multipliers, s);
  };
  const bool rule = sys.pairIterations == 0 || sys.photonIterations == 0;
  const std::vector<SourceWindow> sources = s.responseToys > 0 || rule || s.closureToys > 0 ? readSourceWindows(s, sys) : std::vector<SourceWindow>();
  if (!sources.empty()) {  /* the per-source windows add up to the package's */
    const Window sum = combine(sources, std::vector<double>(sources.size(), 1.0));
    double worst = 0, scale = 0;
    for (int i = 0; i < w.pairResponse.GetNrows(); ++i)
      for (int j = 0; j < w.pairResponse.GetNcols(); ++j)
        worst = std::max(worst, std::fabs(sum.pairResponse(i, j) - w.pairResponse(i, j))), scale = std::max(scale, std::fabs(w.pairResponse(i, j)));
    if (worst > 1e-9 * scale) throw std::runtime_error(sys.name + ": the per-source contributions do not add up to the package response");
  }

  /* the frozen iteration rule, on its own data and response toys, with the simulated prior */
  std::unique_ptr<TH2D> rulePairs, rulePhotons;
  if (rule) {
    TRandom3 dataDraw(s.ruleDataSeed), responseDraw(s.ruleResponseSeed);
    const Toys scanToys = dataToys([&]() { return dataToy(dataDraw); }, w, sys, c, s.ruleDataToys);
    std::vector<ResponseToy> pairToys, photonToys;
    std::vector<TVectorD> responsePhotons;
    for (int t = 0; t < s.ruleResponseToys; ++t) {
      try {
        const Window wt = combine(sources, occurrences(sources, responseDraw));
        const Corrected ct = correct(m, wt, sys);
        if (ct.physical != c.physical) continue;
        pairToys.push_back(ResponseToy{pairProblem(wt), matrixToCells(ct.pairs)});
        photonToys.push_back(ResponseToy{photonProblem(wt), ct.photons});
        responsePhotons.push_back(ct.photons);
      } catch (const std::exception &) {
      }
    }
    std::vector<bool> pairMask(nPt * nX), photonMask(nPt, true);
    for (int i = 0; i < nPt * nX; ++i) pairMask[i] = xj[i % nX] >= sys.reportedLow - 1e-9 && xj[i % nX + 1] <= s.reportedHigh + 1e-9;
    const IterationScan pairScan = scanIterations(pairProblem(w), cells, scanToys.pairs, s.ruleDataToys, pairToys, s.ruleResponseToys, pairMask);
    const IterationScan photonScan = scanIterations(photonProblem(w), c.photons, scanToys.photons, s.ruleDataToys, photonToys, s.ruleResponseToys, photonMask);
    const int pairK = pairScan.selected, photonK = photonRefoldSelection(photonScan, c.photons, scanToys.photons, responsePhotons);
    rulePairs.reset(scanHistogram("iteration_rule_pairs", pairScan, pairK));
    rulePhotons.reset(scanHistogram("iteration_rule_photons", photonScan, photonK));
    if (!pairK || !photonK) {
      std::string why;
      for (const auto *scan : {&pairScan, &photonScan})
        for (const auto &row : scan->rows)
          if (row.k >= 2 && row.k <= 4) why += std::string(scan == &pairScan ? " pairs" : " photons") + " k=" + std::to_string(row.k) + ": " + row.failures + ";";
      throw std::runtime_error(sys.name + ": the frozen iteration rule selects no k (pairs " + std::to_string(pairK) + ", photons " +
                               std::to_string(photonK) + "); the cause is upstream, never absorbed by a larger k. Failing gates:" + why);
    }
    if (sys.pairIterations == 0) sys.pairIterations = pairK;
    if (sys.photonIterations == 0) sys.photonIterations = photonK;
    std::cout << "[FinalAnalysis] " << sys.name << ": the frozen iteration rule selects k = " << pairK << " (pairs), " << photonK << " (photons)" << std::endl;
  }

  /* the data toys and the covariance of the corrected pairs */
  const Toys toys = dataToys([&]() { return dataToy(rng); }, w, sys, c, s.dataToys);
  if (toys.pairs.size() < 2) throw std::runtime_error(sys.name + ": fewer than two accepted data toys");
  const TMatrixD cov = covariance(toys.pairs);

  /* the truth weights the joint unfolding starts from: the fitted prior to Prior.Power, or the simulation (all 1) */
  const bool fitted = s.priorForm != "simulation";
  auto fitOf = [&](const Window &win, const TVectorD &pairs, const TVectorD &photons) {
    return fitted ? fitPrior(win, pairs, photons, cov, s, sys) : PriorFit();
  };
  auto startWeights = [&](const PriorFit &f) {
    if (!fitted) return ones(nX);
    if (!f.ok) throw std::runtime_error("the prior fit did not converge");
    TVectorD start(f.weights);
    for (int x = 0; x < nX; ++x) start[x] = std::pow(start[x], s.priorPower);
    return start;
  };

  const PriorFit fit = fitOf(w, cells, c.photons);
  const TVectorD weights = startWeights(fit);
  TVectorD refolded;
  const TVectorD pairs = unfoldPairsWithPrior(cells, w, sys.pairIterations, weights, s.fakeModel, &refolded);
  const TVectorD photons = unfoldPhotons(c.photons, w, sys.photonIterations);
  const TVectorD result = perPhoton(pairs, photons, xj);
  const TVectorD reco = perPhoton(cells, c.photons, xj);
  const TVectorD simulationPrior = perPhoton(unfoldPairsWithPrior(cells, w, sys.pairIterations, ones(nX), s.fakeModel),
                                             photons, xj);

  /* every accepted toy refitted and unfolded: the central 68% half width */
  std::vector<std::vector<double>> spread(nX);
  std::vector<TVectorD> dataRows, responseRows;
  int finite = 0;
  for (size_t i = 0; i < toys.pairs.size(); ++i) {
    try {
      const TVectorD start = startWeights(fitOf(w, toys.pairs[i], toys.photons[i]));
      const TVectorD r = perPhoton(unfoldPairsWithPrior(toys.pairs[i], w, sys.pairIterations, start, s.fakeModel),
                                   unfoldPhotons(toys.photons[i], w, sys.photonIterations), xj);
      bool ok = true;
      for (int x = 0; x < nX; ++x) ok = ok && std::isfinite(r[x]);
      if (!ok) continue;
      for (int x = 0; x < nX; ++x) spread[x].push_back(r[x]);
      dataRows.push_back(r);
      ++finite;
    } catch (const std::exception &) {
    }
  }
  if (finite < 2) throw std::runtime_error(sys.name + ": fewer than two finite toy results");
  TVectorD dataHalf(nX), responseHalf(nX), statVariance(nX);
  for (int x = 0; x < nX; ++x) dataHalf[x] = 0.5 * (quantile(spread[x], 0.84) - quantile(spread[x], 0.16));

  /* the source bootstrap of the response: every toy resamples the simulated sources with replacement within each sample
     (occurrence multiplicities), and the central data are corrected, refitted and unfolded with it */
  int responseAccepted = 0, responseFinite = 0;
  if (s.responseToys > 0) {
    TRandom3 draw(s.responseSeed);
    std::vector<std::vector<double>> rows(nX);
    for (int t = 0; t < s.responseToys; ++t) {
      try {
        const Window wt = combine(sources, occurrences(sources, draw));
        const Corrected ct = correct(m, wt, sys);
        if (ct.physical != c.physical) continue;
        ++responseAccepted;
        const TVectorD cellsT = matrixToCells(ct.pairs);
        const TVectorD r = perPhoton(unfoldPairsWithPrior(cellsT, wt, sys.pairIterations, startWeights(fitOf(wt, cellsT, ct.photons)), s.fakeModel),
                                     unfoldPhotons(ct.photons, wt, sys.photonIterations), xj);
        bool ok = true;
        for (int x = 0; x < nX; ++x) ok = ok && std::isfinite(r[x]);
        if (!ok) continue;
        for (int x = 0; x < nX; ++x) rows[x].push_back(r[x]);
        responseRows.push_back(r);
        ++responseFinite;
      } catch (const std::exception &) {
      }
    }
    if (responseFinite < 2) throw std::runtime_error(sys.name + ": fewer than two finite response toys");
    for (int x = 0; x < nX; ++x) responseHalf[x] = 0.5 * (quantile(rows[x], 0.84) - quantile(rows[x], 0.16));
  }
  for (int x = 0; x < nX; ++x) statVariance[x] = dataHalf[x] * dataHalf[x] + responseHalf[x] * responseHalf[x];
  std::vector<bool> reportedXj(nX);
  for (int x = 0; x < nX; ++x) reportedXj[x] = xj[x] >= sys.reportedLow - 1e-9 && xj[x + 1] <= s.reportedHigh + 1e-9;
  const TMatrixD dataCovariance = robustCovariance(dataRows, reportedXj);
  const TMatrixD responseCovariance = responseRows.size() >= 2 ? robustCovariance(responseRows, reportedXj) : TMatrixD(nX, nX);
  const std::map<std::string, ClosureMode> closure =
      s.closureToys > 0 ? estimatorClosure(s, sys, sources, cov, c.photons, w) : std::map<std::string, ClosureMode>();

  /* iteration scans (joint: photons at their nominal iterations; photons alone) and the simulation closure */
  TMatrixD scanXj(s.scanIterations, nX), scanPhotons(s.scanIterations, nPt);
  for (int k = 1; k <= s.scanIterations; ++k) {
    const TVectorD xk = perPhoton(unfoldPairsWithPrior(cells, w, k, weights, s.fakeModel), photons, xj);
    const TVectorD gk = unfoldPhotons(c.photons, w, k);
    for (int x = 0; x < nX; ++x) scanXj(k - 1, x) = xk[x];
    for (int i = 0; i < nPt; ++i) scanPhotons(k - 1, i) = gk[i];
  }
  const TVectorD closurePhotonReco = columnSums(w.photonResponse) + w.photonFakes;
  const TVectorD closurePairReco = columnSums(w.pairResponse) + w.pairFakes;
  const TVectorD closurePhotons = bayes(closurePhotonReco, w.photonResponse, w.photonMisses, sys.photonIterations,
                                        w.photonTruth, w.photonFakes);
  const TVectorD closurePairs = bayes(closurePairReco, w.pairResponse, w.pairMisses, sys.pairIterations, w.pairTruth,
                                      w.pairFakes);
  TVectorD reweightedTruth(w.pairTruth);
  for (int t = 0; t < reweightedTruth.GetNrows(); ++t) reweightedTruth[t] *= weights[t % nX];

  out->cd();
  const Edges scanEdges = indexEdges(s.scanIterations);
  for (auto *h : std::vector<TObject *>{
           hist2("abcd", regions, pt, m.abcd, &m.abcdVariance), hist2("leakage", regions, pt, w.leakage, &w.leakageVariance),
           hist1("purity", pt, c.purity), hist1("transfer", pt, c.transfer), hist1("photons_corrected", pt, c.photons),
           hist2("pairs_A", pt, xj, raw.pairsA, &raw.pairsAVariance), hist2("pairs_C_scaled", pt, xj, before.scaledC),
           hist2("pairs_purity_corrected", pt, xj, before.purityCorrected),
           hist2("pairs_combinatoric", pt, xj, combinatoricShown),
           hist2("pairs_corrected", pt, xj, c.pairs, &c.pairsVariance), hist2("pair_covariance", cellEdges, cellEdges, cov),
           hist2("pairs_unfolded", pt, xj, cellsToMatrix(pairs, nPt, nX)),
           hist2("pairs_refolded", pt, xj, cellsToMatrix(refolded, nPt, nX)), hist1("photons_unfolded", pt, photons),
           hist1("xj_per_photon", xj, result, &statVariance), hist1("xj_per_photon_reco", xj, reco),
           hist1("xj_half_width_data", xj, dataHalf), hist1("xj_half_width_response", xj, responseHalf),
           hist2("xj_covariance_data", xj, xj, dataCovariance), hist2("xj_covariance_response", xj, xj, responseCovariance),
           hist1("xj_per_photon_simulation_prior", xj, simulationPrior), hist1("prior_weights", xj, weights),
           hist1("xj_prior_simulation", xj, perPhoton(w.pairTruth, w.photonTruth, xj)),
           hist1("xj_prior_reweighted", xj, perPhoton(reweightedTruth, w.photonTruth, xj)),
           hist2("photon_response", pt, pt, w.photonResponse), hist2("photon_response_grid", w.truthEdges, w.recoEdges, w.responseGridPhoton),
           hist2("pair_response", cellEdges, cellEdges, w.pairResponse), hist1("photon_truth", pt, w.photonTruth),
           hist1("photon_reco", pt, w.photonReco), hist2("pairs_reco_simulation", pt, xj, cellsToMatrix(w.pairReco, nPt, nX)),
           hist1("closure_photons", pt, closurePhotons),
           hist1("closure_pairs", cellEdges, closurePairs), hist1("closure_pairs_truth", cellEdges, w.pairTruth),
           hist2("scan_xj_per_photon", scanEdges, xj, scanXj), hist2("scan_photons", scanEdges, pt, scanPhotons)})
    h->Write();
  if (rulePairs) rulePairs->Write(), rulePhotons->Write();
  const Edges viewEdges = indexEdges(nPt + 1);  /* the photon bins, then their integral */
  for (const auto &[mode, cm] : closure)
    for (const auto &[what, matrix] : std::vector<std::pair<std::string, const TMatrixD *>>{
             {"unfolded", &cm.unfolded}, {"truth", &cm.truth}, {"sigma", &cm.sigma}, {"pull", &cm.pull}})
      hist2(("estimator_closure_" + mode + "_" + what).c_str(), viewEdges, xj, *matrix)->Write();
  std::ostringstream summary;
  summary.precision(10);
  summary << "label=" << sys.label << ";package=" << sys.package << ";centrality=" << sys.centrality
          << ";pair_iterations=" << sys.pairIterations << ";photon_iterations=" << sys.photonIterations
          << ";iteration_rule=" << (rule ? "frozen" : "configured")
          << ";data_toys_accepted=" << toys.pairs.size() << "/" << s.dataToys << ";data_toys_finite=" << finite
          << ";response_toys_accepted=" << responseAccepted << "/" << s.responseToys << ";response_toys_finite=" << responseFinite
          << ";reported_xj=" << sys.reportedLow << "," << s.reportedHigh << ";jet_threshold_gev=" << s.jetThreshold
          << ";prior=" << s.priorForm << ";prior_power=" << s.priorPower << ";data_toy_model=" << s.dataToyModel
          << ";background_correlation=" << join(sys.backgroundCorrelation) << ";purity_scale=" << join(sys.purityScale)
          << ";efficiency_scale=" << join(sys.efficiencyScale)
          << ";combinatoric=" << (!sys.subtractCombinatoric ? "none" : dataDriven ? "data_dphi_" + s.combinatoricShape + "_" + s.combinatoricScope : "simulation");
  if (fitted) summary << ";prior_theta=" << join(fit.theta) << ";prior_chi2=" << fit.chi2 << ";prior_ndf=" << fit.ndf;
  for (const auto &[mode, cm] : closure)  /* toys and failures, then per view row (photon bins, integral) max |pull| and mean pull^2 */
    summary << ";closure_" << mode << "=" << cm.toys << "/" << cm.failed << "," << join(cm.maxPull) << "," << join(cm.chi2);
  TObjString(summary.str().c_str()).Write("settings");
  std::cout << "[FinalAnalysis] " << sys.name << ": purity";
  for (int i = 0; i < nPt; ++i) std::cout << " " << c.purity[i];
  std::cout << ", photons " << photons.Sum() << ", toys " << toys.pairs.size() << "/" << s.dataToys << " (" << finite
            << " finite)";
  if (fitted) std::cout << ", prior theta " << join(fit.theta) << " chi2/ndf " << fit.chi2 << "/" << fit.ndf;
  std::cout << ", reported xJ " << sys.reportedLow << "-" << s.reportedHigh << std::endl;
}

}  /* namespace */

/* Results are written to <output>.partial and renamed only when every system succeeded; an existing output is never
   replaced. */
int run(const Settings &s) {
  if (!gSystem->AccessPathName(s.output.c_str())) throw std::runtime_error(s.output + " exists; choose a new output");
  const std::string coverage = inputCoverage(s);
  const std::string render = renderMetadata(s);
  const std::string partial = s.output + ".partial";
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot write " + partial);
  TRandom3 rng(s.seed);
  for (const auto &sys : s.systems) runSystem(s, sys, out->mkdir(sys.name.c_str()), rng);
  out->cd();
  TObjString(s.name.c_str()).Write("analysis");
  TObjString(coverage.c_str()).Write("input_coverage");
  TObjString(render.c_str()).Write("render_metadata");
  out->Close();
  if (gSystem->Rename(partial.c_str(), s.output.c_str()) != 0) throw std::runtime_error("cannot move results into place");
  renderProductReceipt(s.output, render);
  std::cout << "[FinalAnalysis] wrote " << s.output << std::endl;
  return 0;
}

}  /* namespace fa */

namespace fa { int pipeline(const std::string &config, const std::string &package); }

int main(int argc, char **argv) {
  if (argc == 4 && std::strcmp(argv[1], "--all") == 0) {
    try { return fa::pipeline(argv[2], argv[3]); }
    catch (const std::exception &e) { std::cerr << "FinalAnalysis all: " << e.what() << std::endl; return 1; }
  }
  if (argc != 3 || (std::strcmp(argv[1], "--config") != 0 && std::strcmp(argv[1], "--systematics") != 0 && std::strcmp(argv[1], "--qa") != 0)) {
    std::cerr << "usage: final_analysis --config SETTINGS | --qa SETTINGS | --systematics SYSTEMATICS" << std::endl;
    return 2;
  }
  try {
    if (std::strcmp(argv[1], "--systematics") == 0) return fa::systematics(argv[2]);
    if (std::strcmp(argv[1], "--qa") == 0) return fa::runQA(fa::Settings::read(argv[2]));
    return fa::run(fa::Settings::read(argv[2]));
  } catch (const std::exception &e) {
    std::cerr << "FinalAnalysis: " << e.what() << std::endl;
    return 1;
  }
}
