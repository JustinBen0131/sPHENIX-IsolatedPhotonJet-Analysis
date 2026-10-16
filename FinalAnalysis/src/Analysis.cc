#include "NativeProducerEvidence.h"
/* Analysis: the ordered chain of every configured result, its data toys, and results.root. */

/* final_analysis --config config/analysis.yaml */

/* For each result: measured data and response -> purity and combinatoric correction -> photon-yield unfolding and the
   joint (photon pT, xJ) unfolding from the data-fitted prior -> per-photon xJ. Data toys (independent Poisson counts of
   the measured histograms, or a Poisson(1) multiplier per event of the event records), corrected like the data (a toy
   whose ABCD solution changes branch is rejected), give the covariance of the corrected pairs that weights the prior
   fit; every toy is unfolded from the nominal prior, and the statistical uncertainty is the central 68% half width of
   the toy results. The prior is not refitted per toy (as in RooUnfold and PPG18); its dependence is the prior
   systematic, a variation of its own. results.root holds every quantity the figures draw, so drawing never recomputes. */
#include "FinalAnalysis.h"
#include "HistogramErrors.h"
#include "Phases.h"
#include "DerivedScientificTrace.h"

#include <cerrno>
#include <unistd.h>
#include <sys/wait.h>
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
#include "TTree.h"

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

TVectorD diagonal(const TMatrixD &m) {
  TVectorD out(m.GetNrows());
  for (int i = 0; i < m.GetNrows(); ++i) out[i] = m(i, i);
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
  observeNativeInput(*f);
  auto values = [&](const std::string &name) {
    const std::string key = "qa/" + sys.package + "/purity_study/" + sys.name + "/" + name;
    auto *h = dynamic_cast<TH1 *>(f->Get(key.c_str()));
    if (!h || h->GetNbinsX() + 1 != static_cast<int>(sys.photonEdges.size())) throw std::runtime_error(s.qaOutput + ": no " + key + " on the photon bins");
    sys.lineage[name]=photonjet_science::traceOf(h);
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
  std::vector<int> attempt;
  int exceptions = 0, branchChanges = 0;
};

template <class Draw> Toys dataToys(Draw draw, const Window &w, const SystemSettings &sys, const Corrected &central, int n, Phases &phases) {
  Toys t;
  for (int i = 0; i < n; ++i) {
    const auto progress = phases.attempt();
    try {
      const Corrected c = correct(draw(), w, sys);
      if (c.physical != central.physical) { ++t.branchChanges; continue; }
      t.attempt.push_back(i);
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

/* The covariance saved with a result: the correlations of a toy ensemble with the variances the result quotes, the
   squared central 68% half widths. The trimmed sample variance and the half width agree within a few percent where
   the toys are Gaussian and differ by up to a factor of two in bins with very few events; the quoted error is the half
   width everywhere, so every quantity derived from the covariance (the mean, a ratio) carries the same errors. */
TMatrixD withHalfWidths(const TMatrixD &covariance, const TVectorD &half) {
  const int n = half.GetNrows();
  TMatrixD out(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) {
      const double scale = covariance(i, i) * covariance(j, j);
      // Evaluate once: reversed multiplication can differ by one rounding bit.
      out(j, i) = out(i, j) = scale > 0 ? covariance(i, j) * half[i] * half[j] / std::sqrt(scale) : i == j ? half[i] * half[i] : 0.0;
    }
  return out;
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
std::vector<double> occurrences(const Sources &sources, TRandom3 &draw) {
  std::vector<int> count;
  for (int sample : sources.sample) {
    if (sample >= static_cast<int>(count.size())) count.resize(sample + 1);
    ++count[sample];
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
   resamples and the pulls; per view row the largest |pull| and the mean pull^2 over the reported xJ bins. The photon
   yield of the same pass is kept beside them: unfolded, the test sample's truth, its reconstructed input and the spread
   of unfolded minus truth over the resamples. */
struct ClosureMode {
  TMatrixD unfolded, truth, sigma, pull;
  TVectorD photons, photonTruth, photonInput, photonSigma;
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
   is the spread of unfolded minus truth across those resamples. The photon
   yield unfolded in the same pass is saved per mode (full closure and half
   closure of the 1D unfolding). */
std::map<std::string, ClosureMode> estimatorClosure(const Settings &s, const SystemSettings &sys, const Sources &sources,
                                                    const TMatrixD &dataCovariance, const TVectorD &dataPhotons, const Window &w,
                                                    Phases &phases) {
  const Edges &xj = s.xjEdges;
  const int nX = xj.size() - 1, nPt = sys.photonEdges.size() - 1;
  std::map<int, std::vector<size_t>> bySample;
  for (size_t i = 0; i < sources.size(); ++i) bySample[sources.sample[i]].push_back(i);
  std::vector<double> all(sources.size(), 1.0), a(sources.size()), b(sources.size());
  for (const auto &entry : bySample) {
    std::vector<size_t> ranked = entry.second;
    std::sort(ranked.begin(), ranked.end(), [&](size_t p, size_t q) { return sources.identity[p] < sources.identity[q]; });
    for (size_t r = 0; r < ranked.size(); ++r) (r % 2 ? b : a)[ranked[r]] = 1;
  }
  TMatrixD covariance(dataCovariance);
  for (int i = 0; i < covariance.GetNrows(); ++i)
    for (int j = 0; j < covariance.GetNcols(); ++j)
      covariance(i, j) /= (dataPhotons[i / nX] / w.simPhotons[i / nX]) * (dataPhotons[j / nX] / w.simPhotons[j / nX]);
  struct Pass {
    TMatrixD unfolded, truth;
    TVectorD photons, photonTruth, photonInput;
  };
  auto onePass = [&](const std::vector<double> &train, const std::vector<double> &test) {
    const Window tr = combine(sources, train), te = combine(sources, test);
    TVectorD pairs = columnSums(te.pairResponse), photons = columnSums(te.photonResponse);
    pairs += te.pairFakes, photons += te.photonFakes;
    TVectorD weights = ones(nX);
    if (s.priorForm != "simulation") {
      const PriorFit f = fitPrior(tr, pairs, photons, covariance, s, sys);
      if (!f.ok) throw std::runtime_error("the prior fit did not converge");
      for (int x = 0; x < nX; ++x) weights[x] = 1 + s.priorFraction * (f.weights[x] - 1);
    }
    const TVectorD unfoldedPhotons = unfoldPhotons(photons, tr, sys.photonIterations);
    return Pass{views(unfoldPairsWithPrior(pairs, tr, sys.pairIterations, weights, s.fakeModel), unfoldedPhotons, xj),
                views(te.pairTruth, te.photonTruth, xj), unfoldedPhotons, te.photonTruth, photons};
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
    phases.start(sys.name + "/closure_" + name, s.closureToys, "toy_attempt");
    ClosureMode m;
    const Pass central = onePass(*train, *test);
    m.unfolded.ResizeTo(central.unfolded), m.unfolded = central.unfolded;
    m.truth.ResizeTo(central.truth), m.truth = central.truth;
    m.photons.ResizeTo(nPt), m.photons = central.photons;
    m.photonTruth.ResizeTo(nPt), m.photonTruth = central.photonTruth;
    m.photonInput.ResizeTo(nPt), m.photonInput = central.photonInput;
    std::vector<TMatrixD> differences;
    std::vector<TVectorD> photonDifferences;
    for (int k = 0; k < s.closureToys; ++k) {
      const auto progress = phases.attempt();
      const std::vector<double> tr = resample(*train), te = name == "full" ? tr : resample(*test);
      try {
        const Pass p = onePass(tr, te);
        differences.push_back(p.unfolded - p.truth);
        photonDifferences.push_back(p.photons - p.photonTruth);
      } catch (const std::exception &) {  /* an unconverged fit or a resample the unfolding cannot take */
        ++m.failed;
      }
    }
    m.toys = differences.size();
    m.photonSigma.ResizeTo(nPt);
    for (int i = 0; i < nPt; ++i) {
      double mean = 0, sum2 = 0;
      for (const auto &d : photonDifferences) mean += d[i];
      mean /= photonDifferences.size();
      for (const auto &d : photonDifferences) sum2 += (d[i] - mean) * (d[i] - mean);
      m.photonSigma[i] = photonDifferences.size() >= 2 ? std::sqrt(sum2 / (photonDifferences.size() - 1)) : std::nan("");
    }
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
  h->Sumw2(); // Iteration-rule table: no count-statistics interpretation of bin errors.
  for (int y = 0; y < 12; ++y) h->GetYaxis()->SetBinLabel(y + 1, metrics[y]);
  for (const auto &r : scan.rows) {
    const double v[] = {r.objective, r.dataStat, r.responseStat, r.deviation, r.refold, r.mcClosure, r.mcRefold, r.negativeFraction,
                        r.dataSuccess, r.responseSuccess, r.passes ? 1.0 : 0.0, r.k == selected ? 1.0 : 0.0};
    for (int y = 0; y < 12; ++y) h->SetBinContent(r.k, y + 1, v[y]);
  }
  return h;
}

/* One system's ordered calculation: checked inputs, corrections, prior/iteration selection, unfolding and toys.
   Event and source resampling preserve the correlations needed by the final per-photon uncertainty. */
/* The reconstructed-to-truth association figure on its own bins (Association.PhotonEdges and XjEdges, unions of the
   package bins like the measurement's): the DATA region-A recoil jets before any correction, every simulated
   reconstructed pair, and the simulated pairs matched to their truth pair. The packages are read once more on that
   grid; the response grids keep their edges outside the photon range. Nothing of the measurement uses these. */
struct Association {
  Edges photon, xj;
  TMatrixD data, dataVariance, simulation, matched;
  ScientificLineage lineage;
};
std::unique_ptr<Association> readAssociation(const Settings &s, const SystemSettings &sys) {
  if (std::find(s.associationSystems.begin(), s.associationSystems.end(), sys.name) == s.associationSystems.end()) return nullptr;
  Settings fine = s;
  fine.xjEdges = s.associationXjEdges;
  fine.parameterOrigins["Measurement.XjEdges"]=s.parameterOrigins.at("Association.XjEdges");
  SystemSettings bins = sys;
  bins.photonEdges = s.associationPhotonEdges;
  bins.parameterOrigins["PhotonEdges"]=s.parameterOrigins.at("Association.PhotonEdges");
  auto around = [&](const Edges &grid) {
    Edges edges;
    for (double edge : grid) if (edge < bins.photonEdges.front() - 1e-9) edges.push_back(edge);
    edges.insert(edges.end(), bins.photonEdges.begin(), bins.photonEdges.end());
    for (double edge : grid) if (edge > bins.photonEdges.back() + 1e-9) edges.push_back(edge);
    return edges;
  };
  bins.recoEdges = around(sys.recoEdges);
  bins.truthEdges = around(sys.truthEdges);
  const Measured m = readMeasured(fine, bins);
  const Window w = readWindow(fine, bins);
  const int nPt = bins.photonEdges.size() - 1, nX = fine.xjEdges.size() - 1;
  auto a = std::make_unique<Association>(Association{bins.photonEdges, fine.xjEdges, m.pairsA, m.pairsAVariance,
      cellsToMatrix(w.pairReco, nPt, nX), cellsToMatrix(columnSums(w.pairResponse) + w.pairFakes - w.detector, nPt, nX)});
  if (s.combinatoricMethod == "data_dphi") a->data += m.backgroundA;  /* before the combinatoric subtraction, as pairs_A */
  a->lineage["association_pairs_A"]=lineageOf(m.lineage,"rawPairsA");
  a->lineage["association_pairs_reco_simulation"]=lineageOf(w.lineage,"pairReco");
  a->lineage["association_pairs_reco_matched_simulation"]=derived_science::join("column_sum_plus_fakes_minus_detector",
    {lineageOf(w.lineage,"pairResponse"),lineageOf(w.lineage,"pairFakes"),lineageOf(w.lineage,"detector")});
  return a;
}

void runSystem(const Settings &s, SystemSettings sys, TDirectory *out, TRandom3 &rng) {
  // Concurrent radius calculations share system names but have distinct outputs.
  Phases phases(sys.name + "_" + std::filesystem::path(s.output).stem().string(),
                s.name, sys.photonEdges.size() - 1, sys.centrality);
  try {
  phases.start(sys.name + "/inputs");
  const int requestedPhotonIterations=sys.photonIterations,requestedPairIterations=sys.pairIterations;
  resolveFromQA(s, sys);
  const bool events = s.dataToyModel == "event_bootstrap", dataDriven = s.combinatoricMethod == "data_dphi";
  const EventTable table = events ? readEvents(s, sys) : EventTable();
  const Measured m = events ? fromEvents(table, nullptr, s) : readMeasured(s, sys);
  const Window w = readWindow(s, sys);
  const SimulatedDphi simulated = readSimulatedDphi(s, sys);  /* read here: closing input files resets ROOT's current directory */
  const auto association = readAssociation(s, sys);
  const Corrected c = correct(m, w, sys);
  /* The result uses the physical root of the leakage-aware ABCD equation in every photon bin; the clamped estimate
     solveLeakageABCD falls back to is never a result. Toys of the other branch are rejected below. */
  for (size_t i = 0; i < c.physical.size(); ++i)
    if (!c.physical[i])
      throw std::runtime_error(sys.name + ": no physical leakage-corrected ABCD solution in photon bin " + std::to_string(i));
  /* for the figures: the recoil jets and their correction before any combinatoric subtraction */
  Measured raw = m;
  if (dataDriven) { raw.pairsA += m.backgroundA; raw.pairsC += m.backgroundC; }
  SystemSettings unsubtracted = sys;
  unsubtracted.subtractSideband = false;
  const Corrected before = dataDriven ? correct(raw, w, unsubtracted) : c;
  raw.lineage["pairsA"]=lineageOf(m.lineage,"rawPairsA");
  raw.lineage["pairsC"]=lineageOf(m.lineage,"rawPairsC");
  namespace ds=derived_science;
  const ScientificLineage corrected=ds::corrections(m,w,sys), beforeTrace=dataDriven ? ds::corrections(raw,w,unsubtracted) : corrected;
  ScientificLineage outputTrace;
  const auto wt=[&](const char *name){return lineageOf(w.lineage,name);};
  const auto ct=[&](const char *name){return lineageOf(corrected,name);};
  outputTrace["abcd"]=lineageOf(m.lineage,"abcd"); outputTrace["leakage"]=wt("leakage");
  outputTrace["purity"]=ct("purity"); outputTrace["transfer"]=ct("transfer"); outputTrace["photons_corrected"]=ct("photons");
  outputTrace["pairs_A"]=lineageOf(m.lineage,"rawPairsA");
  outputTrace["pairs_C_scaled"]=lineageOf(beforeTrace,"scaledC");
  outputTrace["pairs_purity_corrected"]=lineageOf(beforeTrace,"purityCorrected");
  outputTrace["pairs_corrected"]=ct("pairs");
  outputTrace["pairs_combinatoric"]=dataDriven ? ds::join("before_minus_corrected_pairs",{lineageOf(beforeTrace,"pairs"),ct("pairs")}) : ct("combinatoric");
  const auto dataSeed=rng.GetSeed();
  TMatrixD combinatoricShown(c.combinatoric);
  if (dataDriven) { combinatoricShown = before.pairs; combinatoricShown -= c.pairs; }
  const Edges &pt = sys.photonEdges, &xj = s.xjEdges;
  const int nPt = pt.size() - 1, nX = xj.size() - 1;
  const Edges regions = indexEdges(4), cellEdges = indexEdges(nPt * nX);
  if (s.correctionsOnly) {
    /* These are the identical correction objects serialized by the full calculation below. No prior fit, toy,
       iteration selection or unfolding is run; the product is explicitly insufficient for any unfolded figure. */
    phases.start(sys.name + "/write_pre_unfolding_diagnostic");
    out->cd();
    for (auto *h : std::vector<TObject *>{
             hist1("purity", pt, c.purity), hist1("photons_corrected", pt, c.photons),
             hist2("pairs_A", pt, xj, raw.pairsA, &raw.pairsAVariance),
             hist2("pairs_purity_corrected", pt, xj, before.purityCorrected, &before.purityCorrectedVariance),
             hist2("pairs_combinatoric", pt, xj, combinatoricShown),
             hist2("pairs_corrected", pt, xj, c.pairs, &c.pairsVariance)}) ds::write(h,outputTrace);
    std::ostringstream settings;
    settings << std::setprecision(10) << "label=" << sys.label << ";package=" << sys.package << ";centrality=" << sys.centrality
             << ";reported_xj=" << sys.reportedLow << "," << s.reportedHigh << ";jet_threshold_gev=" << std::atof(s.jet.c_str() + 3)
             << ";combinatoric=" << (dataDriven ? "data_dphi_" + s.combinatoricShape + "_" + s.combinatoricScope : "simulation")
             << ";dphi_sideband=" << s.sidebandLow << "," << s.sidebandHigh
             << ";product_stage=pre_unfolding_diagnostic;unfolding=NOT_RUN;prior_fit=NOT_RUN";
    TObjString(settings.str().c_str()).Write("settings");
    phases.complete();
    return;
  }
  /* Every jet above the threshold around the region-A leading photons per Delta-phi bin, summed over the photon and xJ
     bins: the reconstructed photon-jet separation before the purity and combinatoric corrections (event records only). */
  TVectorD dphiA(dphiEdges().size() - 1);
  if (events) {
    const int nD = dphiA.GetNrows(), perRegion = nPt * nX * nD * 2;
    for (size_t i = 0; i < table.dphiEvent.size(); ++i)
      if (table.dphiCell[i] / perRegion == 0) dphiA[table.dphiCell[i] % perRegion / 2 % nD] += table.weight[table.dphiEvent[i]];
  }
  const TVectorD cells = matrixToCells(c.pairs);
  std::vector<double> multipliers(table.weight.size());
  auto dataToy = [&](TRandom3 &random) {  /* a data toy: Poisson counts, or one Poisson(1) multiplier per event */
    if (!events) return poissonToy(m, random);
    for (auto &x : multipliers) x = random.Poisson(1.0);
    return fromEvents(table, &multipliers, s);
  };
  /* The iteration scan (Unfolding.cc): for every k the statistical terms, the change from k - 1 and their quadrature
     sum Q_k, for the photon yield and for the joint unfolding, on the scan's own data and response toys. It is the
     stability figure of a result, and where a result's settings say `scan` it gives the iterations: the k with the
     smallest Q_k. Configured iterations are used as given and marked in the scan. A run without response toys (a
     variation, run for its central values) has no scan. */
  const bool chooses = sys.pairIterations == 0 || sys.photonIterations == 0;
  const bool scan = chooses || s.responseToys > 0;
  /* the simulated sources, for whatever resamples them: the iteration scan, the response toys, the closure */
  const Sources sources = scan || s.closureToys > 0 ? readSources(s, sys, w) : Sources();

  std::unique_ptr<TH2D> pairTable, photonTable;  /* the two scans as saved */
  bool usedFittedScan=false;
  if (scan) {
    TRandom3 dataDraw(s.ruleDataSeed), responseDraw(s.ruleResponseSeed);
    phases.start(sys.name + "/scan_data_toys", s.ruleDataToys, "toy_attempt");
    const Toys scanToys = dataToys([&]() { return dataToy(dataDraw); }, w, sys, c, s.ruleDataToys, phases);
    /* IterationRule.Prior fitted: the joint scan starts from the prior the measurement uses (the simulated truth
       reweighted to the data, fitted here with the scan's own data-toy covariance), so it describes the unfolding
       that is run; simulation: from the unreweighted simulated truth. */
    TVectorD scanWeights = ones(nX);
    const bool fittedScan = s.rulePrior == "fitted" && s.priorForm != "simulation" && scanToys.pairs.size() >= 2;
    usedFittedScan=fittedScan;
    if (fittedScan) {
      const PriorFit f = fitPrior(w, cells, c.photons, covariance(scanToys.pairs), s, sys);
      if (!f.ok) throw std::runtime_error(sys.name + ": the iteration scan's prior fit did not converge");
      for (int x = 0; x < nX; ++x) scanWeights[x] = 1 + s.priorFraction * (f.weights[x] - 1);
    }
    auto pairsOf = [&](const Window &window) {
      return fittedScan ? pairProblemWithPrior(window, scanWeights, s.fakeModel) : pairProblem(window);
    };
    std::vector<ResponseToy> pairToys, photonToys;
    phases.start(sys.name + "/scan_response_toys", s.ruleResponseToys, "toy_attempt");
    for (int t = 0; t < s.ruleResponseToys; ++t) {
      const auto progress = phases.attempt();
      try {
        const Window wt = combine(sources, occurrences(sources, responseDraw));
        const Corrected ct = correct(m, wt, sys);
        if (ct.physical != c.physical) continue;
        pairToys.push_back(ResponseToy{pairsOf(wt), matrixToCells(ct.pairs)});
        photonToys.push_back(ResponseToy{photonProblem(wt), ct.photons});
      } catch (const std::exception &) {
      }
    }
    phases.start(sys.name + "/scan");
    std::vector<bool> pairMask(nPt * nX), photonMask(nPt, true);
    for (int i = 0; i < nPt * nX; ++i) pairMask[i] = xj[i % nX] >= sys.reportedLow - 1e-9 && xj[i % nX + 1] <= s.reportedHigh + 1e-9;
    const IterationScan pairScan = scanIterations(pairsOf(w), cells, scanToys.pairs, s.ruleDataToys, pairToys, s.ruleResponseToys, pairMask);
    const IterationScan photonScan =
        scanIterations(photonProblem(w), c.photons, scanToys.photons, s.ruleDataToys, photonToys, s.ruleResponseToys, photonMask);
    if ((sys.pairIterations == 0 && !pairScan.selected) || (sys.photonIterations == 0 && !photonScan.selected))
      throw std::runtime_error(sys.name + ": the iteration scan has no finite Q_k to choose from");
    if (sys.pairIterations == 0) sys.pairIterations = pairScan.selected;
    if (sys.photonIterations == 0) sys.photonIterations = photonScan.selected;
    std::cout << "[FinalAnalysis] " << sys.name << ": smallest Q_k at k = " << pairScan.selected << " (pairs), " << photonScan.selected
              << " (photons); iterations used: " << sys.pairIterations << ", " << sys.photonIterations << std::endl;
    /* a quality check that the used k does not pass is said, never hidden and never a reason to move k */
    for (const auto &[what, table, used] : {std::tuple{"pairs", &pairScan, sys.pairIterations}, std::tuple{"photons", &photonScan, sys.photonIterations}})
      for (const auto &row : table->rows)
        if (row.k == used && !row.failures.empty())
          std::cout << "[FinalAnalysis] " << sys.name << ": " << what << " at k = " << used << " do not pass " << row.failures
                    << " (data refold chi2/rank " << row.refold << ", simulation closure " << row.mcClosure << ")" << std::endl;
    /* the scan's "selected" row marks the iterations the result uses */
    pairTable.reset(scanHistogram("iteration_rule_pairs", pairScan, sys.pairIterations));
    photonTable.reset(scanHistogram("iteration_rule_photons", photonScan, sys.photonIterations));
  }

  /* the data toys and the covariance of the corrected pairs */
  phases.start(sys.name + "/data_toys", s.dataToys, "toy_attempt");
  const Toys toys = dataToys([&]() { return dataToy(rng); }, w, sys, c, s.dataToys, phases);
  if (toys.pairs.size() < 2) throw std::runtime_error(sys.name + ": fewer than two accepted data toys");
  const TMatrixD cov = covariance(toys.pairs);

  /* the truth weights the joint unfolding starts from: 1 + f (w - 1) of the fitted prior w (Prior.Fraction f), or the
     simulation (all 1) */
  const bool fitted = s.priorForm != "simulation";
  auto fitOf = [&](const Window &win, const TVectorD &pairs, const TVectorD &photons) {
    return fitted ? fitPrior(win, pairs, photons, cov, s, sys) : PriorFit();
  };
  auto startWeights = [&](const PriorFit &f) {
    if (!fitted) return ones(nX);
    if (!f.ok) throw std::runtime_error("the prior fit did not converge");
    TVectorD start(f.weights);
    for (int x = 0; x < nX; ++x) start[x] = 1 + s.priorFraction * (start[x] - 1);
    return start;
  };

  phases.start(sys.name + "/central_result");
  const PriorFit fit = fitOf(w, cells, c.photons);
  const TVectorD weights = startWeights(fit);
  TVectorD refolded;
  const TVectorD pairs = unfoldPairsWithPrior(cells, w, sys.pairIterations, weights, s.fakeModel, &refolded);
  const TVectorD photons = unfoldPhotons(c.photons, w, sys.photonIterations);
  const TVectorD result = perPhoton(pairs, photons, xj);
  const TVectorD reco = perPhoton(cells, c.photons, xj);
  const TVectorD simulationPrior = perPhoton(unfoldPairsWithPrior(cells, w, sys.pairIterations, ones(nX), s.fakeModel),
                                             photons, xj);

  /* every accepted toy unfolded from the nominal prior: the central 68% half width. The photon yields and the
     reconstructed per-photon xJ of every toy are kept too: their spread is the statistical error those spectra are
     drawn with. */
  std::vector<std::vector<double>> spread(nX);
  std::vector<TVectorD> dataRows, responseRows, unfoldedPhotonRows, recoRows;
  std::vector<int> dataAttempt, responseAttempt;
  int finite = 0;
  phases.start(sys.name + "/data_unfolding", toys.pairs.size(), "toy_attempt");
  for (size_t i = 0; i < toys.pairs.size(); ++i) {
    const auto progress = phases.attempt();
    try {
      const TVectorD unfoldedPhotons = unfoldPhotons(toys.photons[i], w, sys.photonIterations);
      const TVectorD r = perPhoton(unfoldPairsWithPrior(toys.pairs[i], w, sys.pairIterations, weights, s.fakeModel), unfoldedPhotons, xj);
      bool ok = true;
      for (int x = 0; x < nX; ++x) ok = ok && std::isfinite(r[x]);
      if (!ok) continue;
      const TVectorD recoRow = perPhoton(toys.pairs[i], toys.photons[i], xj);  /* before any row is kept: a toy counts in all or none */
      for (int x = 0; x < nX; ++x) spread[x].push_back(r[x]);
      dataRows.push_back(r);
      dataAttempt.push_back(toys.attempt.at(i));
      unfoldedPhotonRows.push_back(unfoldedPhotons);
      recoRows.push_back(recoRow);
      ++finite;
    } catch (const std::exception &) {
    }
  }
  if (finite < 2) throw std::runtime_error(sys.name + ": fewer than two finite toy results");
  TVectorD dataHalf(nX), responseHalf(nX), statVariance(nX);
  for (int x = 0; x < nX; ++x) dataHalf[x] = 0.5 * (quantile(spread[x], 0.84) - quantile(spread[x], 0.16));

  /* the source bootstrap of the response: every toy resamples the simulated sources with replacement within each sample
     (occurrence multiplicities), and the central data are corrected and unfolded with it from the nominal prior */
  int responseAccepted = 0, responseFinite = 0;
  if (s.responseToys > 0) {
    TRandom3 draw(s.responseSeed);
    std::vector<std::vector<double>> rows(nX);
    phases.start(sys.name + "/response_toys", s.responseToys, "toy_attempt");
    for (int t = 0; t < s.responseToys; ++t) {
      const auto progress = phases.attempt();
      try {
        const Window wt = combine(sources, occurrences(sources, draw));
        const Corrected ct = correct(m, wt, sys);
        if (ct.physical != c.physical) continue;
        ++responseAccepted;
        const TVectorD cellsT = matrixToCells(ct.pairs);
        const TVectorD r = perPhoton(unfoldPairsWithPrior(cellsT, wt, sys.pairIterations, weights, s.fakeModel),
                                     unfoldPhotons(ct.photons, wt, sys.photonIterations), xj);
        bool ok = true;
        for (int x = 0; x < nX; ++x) ok = ok && std::isfinite(r[x]);
        if (!ok) continue;
        for (int x = 0; x < nX; ++x) rows[x].push_back(r[x]);
        responseRows.push_back(r);
        responseAttempt.push_back(t);
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
  const TMatrixD dataCovariance = withHalfWidths(robustCovariance(dataRows, reportedXj), dataHalf);
  const TMatrixD responseCovariance =
      responseRows.size() >= 2 ? withHalfWidths(robustCovariance(responseRows, reportedXj), responseHalf) : TMatrixD(nX, nX);
  const std::map<std::string, ClosureMode> closure =
      s.closureToys > 0 ? estimatorClosure(s, sys, sources, cov, c.photons, w, phases) : std::map<std::string, ClosureMode>();

  /* iteration scans (joint: photons at their nominal iterations; photons alone) and the simulation closure */
  TMatrixD scanXj(s.scanIterations, nX), scanPhotons(s.scanIterations, nPt);
  phases.start(sys.name + "/iteration_scan", s.scanIterations, "iteration_attempt");
  for (int k = 1; k <= s.scanIterations; ++k) {
    const auto progress = phases.attempt();
    const TVectorD xk = perPhoton(unfoldPairsWithPrior(cells, w, k, weights, s.fakeModel), photons, xj);
    const TVectorD gk = unfoldPhotons(c.photons, w, k);
    for (int x = 0; x < nX; ++x) scanXj(k - 1, x) = xk[x];
    for (int i = 0; i < nPt; ++i) scanPhotons(k - 1, i) = gk[i];
  }
  phases.start(sys.name + "/write_result");
  const TVectorD closurePhotonReco = columnSums(w.photonResponse) + w.photonFakes;
  const TVectorD closurePairReco = columnSums(w.pairResponse) + w.pairFakes;
  const TVectorD closurePhotons = bayes(closurePhotonReco, w.photonResponse, w.photonMisses, sys.photonIterations,
                                        w.photonTruth, w.photonFakes);
  const TVectorD closurePairs = bayes(closurePairReco, w.pairResponse, w.pairMisses, sys.pairIterations, w.pairTruth,
                                      w.pairFakes);
  const TVectorD matchedPairReco = columnSums(w.pairResponse) + w.pairFakes - w.detector;  /* boundary pairs are matched pairs */
  TVectorD reweightedTruth(w.pairTruth);
  for (int t = 0; t < reweightedTruth.GetNrows(); ++t) reweightedTruth[t] *= weights[t % nX];

  const TVectorD photonVariance = diagonal(covariance(toys.photons));
  const TVectorD unfoldedPhotonVariance = diagonal(covariance(unfoldedPhotonRows)), recoVariance = diagonal(covariance(recoRows));

  auto dataToysTrace=ds::join("data_toy_covariance",{ct("pairs")},
    {{"model",s.dataToyModel},{"attempts",s.dataToys},{"seed",dataSeed},{"branch","same_physical_abcd_root"}});
  auto photonToysTrace=ct("photons").then("data_toy_covariance",
    {{"model",s.dataToyModel},{"attempts",s.dataToys},{"seed",dataSeed},{"branch","same_physical_abcd_root"}});
  dataToysTrace=withSettings(dataToysTrace,s.parameterOrigins,{"Uncertainty.DataToyModel","Uncertainty.DataToys","Uncertainty.Seed"});
  photonToysTrace=withSettings(photonToysTrace,s.parameterOrigins,{"Uncertainty.DataToyModel","Uncertainty.DataToys","Uncertainty.Seed"});
  const auto photonProblemTrace=ds::response(w,true),pairProblemTrace=ds::response(w,false);
  auto prior=wt("pairTruth").then("simulation_prior");
  if (fitted) prior=ds::join("fit_joint_prior",{pairProblemTrace,ct("pairs"),ct("photons"),dataToysTrace},
    {{"form",s.priorForm},{"fraction",s.priorFraction},{"fake_model",s.fakeModel},{"fit_range",s.fitRange},
     {"reported_low",s.fitRange=="reported" ? sys.reportedLow : xj.front()},
     {"reported_high",s.fitRange=="reported" ? s.reportedHigh : xj.back()}});
  prior=withSettings(prior,s.parameterOrigins,{"Prior.Form"});
  if(fitted)prior=withSettings(prior,s.parameterOrigins,{"Prior.Fraction","Prior.FakeModel","Prior.FitRange"});
  if(fitted && s.fitRange=="reported")prior=withSettings(prior,s.parameterOrigins,{"Measurement.ReportedXjMax"});
  auto pairScanInputs=ds::join("joint_iteration_inputs",{pairProblemTrace,ct("pairs")},{{"reported_low",sys.reportedLow}});
  if(usedFittedScan)pairScanInputs=ds::join("fitted_scan_photon_normalization",{pairScanInputs,ct("photons")});
  const auto pairScanTrace=ds::iterationRule(pairScanInputs,s,false,usedFittedScan);
  const auto photonScanTrace=ds::iterationRule(ds::join("photon_iteration_inputs",{photonProblemTrace,ct("photons")}),s,true,false);
  auto photonsTrace=ds::join("bayes_photons",{ct("photons"),photonProblemTrace},{{"iterations",sys.photonIterations}});
  auto pairsTrace=ds::join("bayes_joint_with_prior",{ct("pairs"),pairProblemTrace,prior},{{"iterations",sys.pairIterations},{"fake_model",s.fakeModel}});
  photonsTrace=ds::iterationChoice(photonsTrace,sys,"PhotonIterations",requestedPhotonIterations,photonScanTrace);
  pairsTrace=ds::iterationChoice(withSettings(pairsTrace,s.parameterOrigins,{"Prior.FakeModel"}),sys,"PairIterations",requestedPairIterations,pairScanTrace);
  const auto recoTrace=ds::join("per_photon_density",{ct("pairs"),ct("photons")},{{"xj_edges",xj}});
  const auto resultTrace=ds::join("per_photon_density",{pairsTrace,photonsTrace},{{"xj_edges",xj}});
  auto dataResultTrace=ds::join("data_toy_quantile_and_robust_covariance",{resultTrace,dataToysTrace},
    {{"quantiles",{0.16,0.84}},{"refit_prior",false},{"reported_low",sys.reportedLow},{"reported_high",s.reportedHigh}});
  dataResultTrace=withSettings(dataResultTrace,s.parameterOrigins,{"Measurement.ReportedXjMax"});
  auto responseResultTrace=resultTrace.then("response_source_bootstrap_quantile_and_robust_covariance",
    {{"attempts",s.responseToys},{"seed",s.responseToys>0?s.responseSeed:0},{"within_sample",true},{"refit_prior",false},
     {"quantiles",{0.16,0.84}},{"reported_low",sys.reportedLow},{"reported_high",s.reportedHigh}});
  responseResultTrace=withSettings(responseResultTrace,s.parameterOrigins,{"Uncertainty.ResponseToys"});
  if(s.responseToys>0)responseResultTrace=withSettings(responseResultTrace,s.parameterOrigins,{"Uncertainty.ResponseSeed","Measurement.ReportedXjMax"});
  outputTrace["pair_covariance"]=dataToysTrace;
  outputTrace["photons_corrected"]=ds::join("central_value_with_toy_variance",{ct("photons"),photonToysTrace});
  outputTrace["pairs_unfolded"]=pairsTrace;
  outputTrace["pairs_refolded"]=ds::join("refold_joint",{pairsTrace,pairProblemTrace,prior});
  outputTrace["photons_unfolded"]=ds::join("central_value_with_toy_variance",{photonsTrace,photonToysTrace});
  outputTrace["xj_per_photon"]=ds::join("central_with_quadrature_half_widths",{resultTrace,dataResultTrace,responseResultTrace});
  outputTrace["xj_per_photon_reco"]=ds::join("central_value_with_toy_variance",{recoTrace,dataToysTrace,photonToysTrace});
  outputTrace["xj_half_width_data"]=dataResultTrace; outputTrace["xj_covariance_data"]=dataResultTrace;
  outputTrace["xj_half_width_response"]=responseResultTrace; outputTrace["xj_covariance_response"]=responseResultTrace;
  outputTrace["xj_per_photon_simulation_prior"]=ds::join("bayes_simulation_prior_per_photon",{ct("pairs"),pairProblemTrace,photonsTrace},
    {{"iterations",sys.pairIterations},{"fake_model",s.fakeModel},{"xj_edges",xj}});
  outputTrace["xj_per_photon_simulation_prior"]=ds::iterationChoice(withSettings(outputTrace["xj_per_photon_simulation_prior"],s.parameterOrigins,{"Prior.FakeModel"}),sys,"PairIterations",requestedPairIterations,pairScanTrace);
  outputTrace["prior_weights"]=prior;
  outputTrace["xj_prior_simulation"]=ds::join("per_photon_density",{wt("pairTruth"),wt("photonTruth")},{{"xj_edges",xj}});
  outputTrace["xj_prior_reweighted"]=ds::join("per_photon_density_reweighted",{wt("pairTruth"),wt("photonTruth"),prior},{{"xj_edges",xj}});
  for (const auto &entry : std::map<std::string,std::string>{{"photon_response","photonResponse"},{"photon_response_grid","responseGridPhoton"},
       {"pair_response","pairResponse"},{"photon_truth","photonTruth"},{"photon_reco","photonReco"},{"pairs_reco_simulation","pairReco"},
       {"closure_pairs_truth","pairTruth"}}) outputTrace[entry.first]=lineageOf(w.lineage,entry.second);
  outputTrace["pairs_reco_expected"]=ds::join("column_sum_plus_fakes",{wt("pairResponse"),wt("pairFakes")});
  outputTrace["pairs_reco_matched_simulation"]=ds::join("column_sum_plus_fakes_minus_detector",{wt("pairResponse"),wt("pairFakes"),wt("detector")});
  outputTrace["closure_photons"]=ds::iterationChoice(photonProblemTrace.then("bayes_self_closure",{{"iterations",sys.photonIterations}}),sys,"PhotonIterations",requestedPhotonIterations,photonScanTrace);
  outputTrace["closure_pairs"]=ds::iterationChoice(pairProblemTrace.then("bayes_self_closure",{{"iterations",sys.pairIterations}}),sys,"PairIterations",requestedPairIterations,pairScanTrace);
  outputTrace["scan_xj_per_photon"]=ds::join("joint_iterations_per_photon",{ct("pairs"),pairProblemTrace,prior,photonsTrace},
    {{"maximum_iterations",s.scanIterations},{"fake_model",s.fakeModel},{"xj_edges",xj}});
  outputTrace["scan_xj_per_photon"]=withSettings(outputTrace["scan_xj_per_photon"],s.parameterOrigins,{"Unfolding.ScanIterations","Prior.FakeModel"});
  outputTrace["scan_photons"]=ds::join("photon_iterations",{ct("photons"),photonProblemTrace},{{"maximum_iterations",s.scanIterations}});
  outputTrace["scan_photons"]=withSettings(outputTrace["scan_photons"],s.parameterOrigins,{"Unfolding.ScanIterations"});
  outputTrace["iteration_rule_pairs"]=withSettings(pairScanTrace.then("mark_used_iteration",{{"used",sys.pairIterations}}),sys.parameterOrigins,{"PairIterations"});
  outputTrace["iteration_rule_photons"]=withSettings(photonScanTrace.then("mark_used_iteration",{{"used",sys.photonIterations}}),sys.parameterOrigins,{"PhotonIterations"});
  outputTrace["dphi_region_A"]=ds::join("count_region_A_dphi",{lineageOf(table.lineage,"jets"),lineageOf(table.lineage,"events")},{{"edges",dphiEdges()}});
  outputTrace["dphi_simulation_truthjet"]=lineageOf(simulated.lineage,"truthJet");
  outputTrace["dphi_simulation_combinatoric"]=lineageOf(simulated.lineage,"combinatoric");
  if (association) outputTrace.insert(association->lineage.begin(),association->lineage.end());
  for (const auto &[mode,cm] : closure) {
    auto trace=ds::join("estimator_source_closure",{pairProblemTrace,photonProblemTrace},
      {{"mode",mode},{"prior",s.priorForm},{"fake_model",s.fakeModel},
       {"pair_iterations",sys.pairIterations},{"photon_iterations",sys.photonIterations}});
    trace=withSettings(withSettings(trace,sys.parameterOrigins,{"PairIterations","PhotonIterations"}),s.parameterOrigins,{"Prior.Form","Prior.FakeModel"});
    if(fitted) {
      trace=ds::join("closure_prior_fit_with_data_covariance",{trace,ct("photons"),dataToysTrace},
        {{"fraction",s.priorFraction},{"fit_range",s.fitRange}});
      trace=withSettings(trace,s.parameterOrigins,{"Prior.FitRange","Prior.Fraction"});
      if(s.fitRange=="reported")trace=withSettings(trace,s.parameterOrigins,{"Measurement.ReportedXjMax"});
    }
    if(requestedPairIterations==0)trace=ds::join("scan_selected_iterations",{trace,pairScanTrace});
    if(requestedPhotonIterations==0)trace=ds::join("scan_selected_iterations",{trace,photonScanTrace});
    const auto uncertainty=withSettings(trace.then("source_resampling_uncertainty",{{"toys",s.closureToys},{"seed",s.closureSeed}}),
      s.parameterOrigins,{"Closure.Toys","Closure.Seed"});
    for(const auto *what:{"unfolded","truth","sigma","pull","photons_unfolded","photons_truth","photons_input","photons_sigma"}) {
      const std::string quantity=what;
      const bool sampled=quantity=="sigma" || quantity=="pull" || quantity=="photons_sigma";
      outputTrace["estimator_closure_"+mode+"_"+what]=(sampled?uncertainty:trace).then("closure_quantity",{{"quantity",what}});
    }
  }
  std::map<std::string,ds::Json> outputCovariance;
  outputCovariance["photons_corrected"]=ds::covarianceRecord(covariance(toys.photons),"data toy covariance");
  outputCovariance["pairs_corrected"]=ds::covarianceRecord(cov,"data toy covariance; photon-major xj-minor ordering");
  outputCovariance["photons_unfolded"]=ds::covarianceRecord(covariance(unfoldedPhotonRows),"data toy covariance with fixed nominal prior");
  outputCovariance["xj_per_photon_reco"]=ds::covarianceRecord(covariance(recoRows),"data toy covariance");
  outputCovariance["xj_per_photon"]=ds::covarianceRecord(dataCovariance+responseCovariance,"sum of independent robust data and response covariance with central-68-percent half widths");

  out->cd();
  const Edges scanEdges = indexEdges(s.scanIterations);
  for (auto *h : std::vector<TObject *>{
           hist2("abcd", regions, pt, m.abcd, &m.abcdVariance), hist2("leakage", regions, pt, w.leakage, &w.leakageVariance),
           hist1("purity", pt, c.purity), hist1("transfer", pt, c.transfer), hist1("photons_corrected", pt, c.photons, &photonVariance),
           hist2("pairs_A", pt, xj, raw.pairsA, &raw.pairsAVariance), hist2("pairs_C_scaled", pt, xj, before.scaledC, &before.scaledCVariance),
           hist2("pairs_purity_corrected", pt, xj, before.purityCorrected, &before.purityCorrectedVariance),
           hist2("pairs_combinatoric", pt, xj, combinatoricShown),
           hist2("pairs_corrected", pt, xj, c.pairs, &c.pairsVariance), hist2("pair_covariance", cellEdges, cellEdges, cov),
           hist2("pairs_unfolded", pt, xj, cellsToMatrix(pairs, nPt, nX)),
           hist2("pairs_refolded", pt, xj, cellsToMatrix(refolded, nPt, nX)), hist1("photons_unfolded", pt, photons, &unfoldedPhotonVariance),
           hist1("xj_per_photon", xj, result, &statVariance), hist1("xj_per_photon_reco", xj, reco, &recoVariance),
           hist1("xj_half_width_data", xj, dataHalf), hist1("xj_half_width_response", xj, responseHalf),
           hist2("xj_covariance_data", xj, xj, dataCovariance), hist2("xj_covariance_response", xj, xj, responseCovariance),
           hist1("xj_per_photon_simulation_prior", xj, simulationPrior), hist1("prior_weights", xj, weights),
           hist1("xj_prior_simulation", xj, perPhoton(w.pairTruth, w.photonTruth, xj)),
           hist1("xj_prior_reweighted", xj, perPhoton(reweightedTruth, w.photonTruth, xj)),
           hist2("photon_response", pt, pt, w.photonResponse), hist2("photon_response_grid", w.truthEdges, w.recoEdges, w.responseGridPhoton),
           hist2("pair_response", cellEdges, cellEdges, w.pairResponse), hist1("photon_truth", pt, w.photonTruth),
           hist1("photon_reco", pt, w.photonReco), hist2("pairs_reco_simulation", pt, xj, cellsToMatrix(w.pairReco, nPt, nX)),
           hist2("pairs_reco_expected", pt, xj, cellsToMatrix(closurePairReco, nPt, nX)),
           hist2("pairs_reco_matched_simulation", pt, xj, cellsToMatrix(matchedPairReco, nPt, nX)),
           hist1("closure_photons", pt, closurePhotons),
           hist1("closure_pairs", cellEdges, closurePairs), hist1("closure_pairs_truth", cellEdges, w.pairTruth),
           hist2("scan_xj_per_photon", scanEdges, xj, scanXj), hist2("scan_photons", scanEdges, pt, scanPhotons)})
    ds::write(h,outputTrace,outputCovariance);
  if (association)  /* the same three histograms as pairs_A, pairs_reco_simulation and pairs_reco_matched_simulation */
    for (auto *h : {hist2("association_pairs_A", association->photon, association->xj, association->data, &association->dataVariance),
                    hist2("association_pairs_reco_simulation", association->photon, association->xj, association->simulation),
                    hist2("association_pairs_reco_matched_simulation", association->photon, association->xj, association->matched)})
      ds::write(h,outputTrace,outputCovariance);
  if (events) ds::write(hist1("dphi_region_A", dphiEdges(), dphiA, &dphiA),outputTrace);
  /* the simulation's two components of the same distribution, for the Delta-phi figure */
  ds::write(hist1("dphi_simulation_truthjet", dphiEdges(), simulated.truthJet, &simulated.truthJetVariance),outputTrace);
  ds::write(hist1("dphi_simulation_combinatoric", dphiEdges(), simulated.combinatoric, &simulated.combinatoricVariance),outputTrace);
  if (pairTable) { ds::write(pairTable.get(),outputTrace); ds::write(photonTable.get(),outputTrace); }
  const Edges viewEdges = indexEdges(nPt + 1);  /* the photon bins, then their integral */
  for (const auto &[mode, cm] : closure)
    for (const auto &[what, matrix] : std::vector<std::pair<std::string, const TMatrixD *>>{
             {"unfolded", &cm.unfolded}, {"truth", &cm.truth}, {"sigma", &cm.sigma}, {"pull", &cm.pull}})
      ds::write(hist2(("estimator_closure_" + mode + "_" + what).c_str(), viewEdges, xj, *matrix),outputTrace);
  for (const auto &[mode, cm] : closure)  /* the photon yield of the same pass: full closure and the two half closures */
    for (const auto &[what, vector] : std::vector<std::pair<std::string, const TVectorD *>>{
             {"photons_unfolded", &cm.photons}, {"photons_truth", &cm.photonTruth}, {"photons_input", &cm.photonInput},
             {"photons_sigma", &cm.photonSigma}})
      ds::write(hist1(("estimator_closure_" + mode + "_" + what).c_str(), pt, *vector),outputTrace);
  std::ostringstream summary;
  summary.precision(10);
  summary << "label=" << sys.label << ";package=" << sys.package << ";centrality=" << sys.centrality
          << ";pair_iterations=" << sys.pairIterations << ";photon_iterations=" << sys.photonIterations
          << ";iteration_rule=" << (chooses ? "smallest_Q" : "configured")
          << ";data_toys_accepted=" << toys.pairs.size() << "/" << s.dataToys << ";data_toys_finite=" << finite
          << ";response_toys_accepted=" << responseAccepted << "/" << s.responseToys << ";response_toys_finite=" << responseFinite
          << ";reported_xj=" << sys.reportedLow << "," << s.reportedHigh << ";jet_threshold_gev=" << s.jetThreshold
          << ";prior=" << s.priorForm << ";prior_fraction=" << s.priorFraction << ";data_toy_model=" << s.dataToyModel
          << ";background_correlation=" << join(sys.backgroundCorrelation) << ";purity_scale=" << join(sys.purityScale)
          << ";efficiency_scale=" << join(sys.efficiencyScale)
          << ";combinatoric=" << (!sys.subtractCombinatoric ? "none" : dataDriven ? "data_dphi_" + s.combinatoricShape + "_" + s.combinatoricScope : "simulation")
          << ";dphi_sideband=" << s.sidebandLow << "," << s.sidebandHigh;
  if (!sys.responseComponent.empty())summary<<";response_component="<<sys.responseComponent<<";response_recipe="<<sys.responseRecipe<<";response_grid_sha256="<<sys.responseGridSha256;
  if (fitted) summary << ";prior_theta=" << join(fit.theta) << ";prior_chi2=" << fit.chi2 << ";prior_ndf=" << fit.ndf;
  for (const auto &[mode, cm] : closure)  /* toys and failures, then per view row (photon bins, integral) max |pull| and mean pull^2 */
    summary << ";closure_" << mode << "=" << cm.toys << "/" << cm.failed << "," << join(cm.maxPull) << "," << join(cm.chi2);
  // Numeric annotations are exported from the same native calculation. A
  // renderer never reconstructs their lineage from semicolon-delimited prose.
  auto *metrics=out->mkdir("physics_settings");
  auto metric=[&](const std::string &name,double value,ds::Trace trace,const std::string &method) {
    const ds::Json record={{"schema","NativeScientificScalarV1"},{"trace",trace.json()},
      {"measurement",photonjet_science::scalar(value,"iterations",method)},
      {"covariance",photonjet_science::covariance("not_applicable",{},0,"deterministic selected iteration")}};
    metrics->WriteTObject(new TObjString(record.dump().c_str()),name.c_str());
  };
  auto selectedPair=withSettings(ds::Trace{},sys.parameterOrigins,{"PairIterations"});
  auto selectedPhoton=withSettings(ds::Trace{},sys.parameterOrigins,{"PhotonIterations"});
  if(requestedPairIterations==0)selectedPair=ds::join("smallest_Q_selected_iteration",{selectedPair,pairScanTrace});
  if(requestedPhotonIterations==0)selectedPhoton=ds::join("smallest_Q_selected_iteration",{selectedPhoton,photonScanTrace});
  metric("pair_iterations",sys.pairIterations,selectedPair,requestedPairIterations==0?"smallest_Q":"configured");
  metric("photon_iterations",sys.photonIterations,selectedPhoton,requestedPhotonIterations==0?"smallest_Q":"configured");
  out->cd();
  TObjString(summary.str().c_str()).Write("settings");
  if(!sys.excludedDataRuns.empty()) {
    const photonjet_science::Json mask={{"schema","DataRunExclusionV1"},{"runs",sys.excludedDataRuns},
      {"selected_event_positions",table.weight.size()},{"excluded_events",table.excludedEvents},
      {"algorithm","zero_event_weight_preserve_bootstrap_positions"},{"simulation_changed",false}};
    TObjString(mask.dump().c_str()).Write("data_run_exclusion");
  }
  if(s.saveToyRecords) {
    const auto save=[&](const char* name,const std::vector<TVectorD>& rows,const std::vector<int>& ids,const photonjet_science::Trace& trace) {
      if(rows.size()!=ids.size())throw std::runtime_error("toy attempt identity mismatch");
      TTree tree(name,"Finite native unfolded replicas, identified by original attempted draw");
      tree.SetDirectory(nullptr);
      int attempt=0;std::vector<double> values;
      tree.Branch("attempt",&attempt);tree.Branch("values",&values);
      for(size_t i=0;i<rows.size();++i) {
        attempt=ids[i];values.assign(rows[i].GetMatrixArray(),rows[i].GetMatrixArray()+rows[i].GetNrows());tree.Fill();
      }
      photonjet_science::attachTrace(&tree,trace);tree.Write();
    };
    save("data_toy_records",dataRows,dataAttempt,dataToysTrace);
    save("response_toy_records",responseRows,responseAttempt,outputTrace.at("xj_covariance_response"));
    const photonjet_science::Json contract={{"schema","PairedToyRecordsV1"},{"data_seed",dataSeed},{"response_seed",s.responseSeed},
      {"data_attempts",s.dataToys},{"response_attempts",s.responseToys},{"event_positions",table.weight.size()},
      {"xj_edges",xj},{"matching","intersect original attempt IDs within each ensemble"},
      {"comparison_requires","same event ordering, parent package, settings except declared run veto, and fixed seeds"}};
    TObjString(contract.dump().c_str()).Write("toy_record_contract");
  }
  std::cout << "[FinalAnalysis] " << sys.name << ": purity";
  for (int i = 0; i < nPt; ++i) std::cout << " " << c.purity[i];
  std::cout << ", photons " << photons.Sum() << ", toys " << toys.pairs.size() << "/" << s.dataToys << " (" << finite
            << " finite)";
  if (fitted) std::cout << ", prior theta " << join(fit.theta) << " chi2/ndf " << fit.chi2 << "/" << fit.ndf;
  std::cout << ", reported xJ " << sys.reportedLow << "-" << s.reportedHigh << std::endl;
  phases.complete();
  } catch (const std::exception &error) { phases.fail(error.what()); throw; }
}

}  /* namespace */

/* Results are written to <output>.partial and renamed only when every system succeeded; an existing output is never
   replaced. */
int run(const Settings &s) {
  resetNativeInputs();
  if (!gSystem->AccessPathName(s.output.c_str())) throw std::runtime_error(s.output + " exists; choose a new output");
  Phases phases;
  phases.start("checking the inputs");
  const std::string coverage = inputCoverage(s);
  const std::string render = renderMetadata(s);
  phases.start("the results");
  const std::string partial = s.output + ".partial";
  /* One process per result. The results are independent, and each has its own random stream (the seed and its
     position in Systems), so its numbers do not depend on which results are run together or at once. */
  std::vector<pid_t> children;
  std::cout.flush(); std::cerr.flush();
  for (size_t i = 0; i < s.systems.size(); ++i) {
    const pid_t pid = fork();
    if (pid < 0) throw std::runtime_error("cannot start the process of " + s.systems[i].name);
    if (pid == 0) {
      int code = 0;
      try {
        const std::string part = partial + "." + s.systems[i].name;
        std::unique_ptr<TFile> file(TFile::Open(part.c_str(), "RECREATE"));
        if (!file || file->IsZombie()) throw std::runtime_error("cannot write " + part);
        TRandom3 rng(s.seed + 7919U * static_cast<unsigned>(i));
        runSystem(s, s.systems[i], file->mkdir(s.systems[i].name.c_str()), rng);
        writeNativeProducerWitness(*file);
        requireExplicitHistogramErrors(file.get());
        file->Close();
      } catch (const std::exception &e) {
        std::cerr << "FinalAnalysis: " << e.what() << std::endl;
        code = 1;
      }
      std::cout.flush(); std::cerr.flush();
      _exit(code);  /* the parent owns the output and every exit handler */
    }
    children.push_back(pid);
  }
  bool failed = false;
  for (const pid_t pid : children) {
    int status = 0;
    while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
    failed = failed || !WIFEXITED(status) || WEXITSTATUS(status) != 0;
  }
  if (failed) {
    for (const auto &sys : s.systems) gSystem->Unlink((partial + "." + sys.name).c_str());
    throw std::runtime_error("a result failed; its message is above");
  }
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot write " + partial);
  for (const auto &sys : s.systems) {
    const std::string part = partial + "." + sys.name;
    std::unique_ptr<TFile> file(TFile::Open(part.c_str(), "READ"));
    auto *saved = file && !file->IsZombie() ? file->GetDirectory(sys.name.c_str()) : nullptr;
    if (!saved) throw std::runtime_error("the process of " + sys.name + " left no result");
    observeNativeInput(*file);
    copyDirectory(saved, out->mkdir(sys.name.c_str()));
    file.reset();
    gSystem->Unlink(part.c_str());
  }
  phases.start("writing the results");
  out->cd();
  TObjString(s.name.c_str()).Write("analysis");
  if (s.correctionsOnly) TObjString("PRE_UNFOLDING_DIAGNOSTIC_ONLY").Write("product_stage");
  TObjString(coverage.c_str()).Write("input_coverage");
  TObjString(render.c_str()).Write("render_metadata");
  TObjString(s.nativeParameterExports.dump().c_str()).Write("native_parameter_exports");
  writeNativeProducerWitness(*out);
  TObjString(scientificParameterBindings(s).dump().c_str()).Write("scientific_parameter_bindings");
  requireExplicitHistogramErrors(out.get());
  out->Close();
  if (gSystem->Rename(partial.c_str(), s.output.c_str()) != 0) throw std::runtime_error("cannot move results into place");
  renderProductReceipt(s.output, render);
  std::cout << "[FinalAnalysis] wrote " << s.output << std::endl;
  return 0;
}

}  /* namespace fa */

namespace fa { int pipeline(const std::string &config, const std::string &package); }

int main(int argc, char **argv) {
  if (argc == 4 && std::strcmp(argv[1], "--compact") == 0) {
    try { return fa::compactCollection(argv[2], argv[3]); }
    catch (const std::exception &e) { std::cerr << "FinalAnalysis compact: " << e.what() << std::endl; return 1; }
  }
  if (argc == 4 && std::strcmp(argv[1], "--model-qa") == 0) {
    try { return fa::runModelQA(argv[2], argv[3]); }
    catch (const std::exception &e) { std::cerr << "FinalAnalysis model QA: " << e.what() << std::endl; return 1; }
  }
  if (argc == 4 && std::strcmp(argv[1], "--all") == 0) {
    try { return fa::pipeline(argv[2], argv[3]); }
    catch (const std::exception &e) { std::cerr << "FinalAnalysis all: " << e.what() << std::endl; return 1; }
  }
  if (argc != 3 || (std::strcmp(argv[1], "--config") != 0 && std::strcmp(argv[1], "--corrections-only") != 0 &&
                   std::strcmp(argv[1], "--systematics") != 0 && std::strcmp(argv[1], "--qa") != 0)) {
    std::cerr << "usage: final_analysis --config SETTINGS | --corrections-only SETTINGS | --qa SETTINGS | --model-qa MODELS OUTPUT | --systematics SYSTEMATICS | --compact COLLECTION DIRECTORY" << std::endl;
    return 2;
  }
  try {
    if (std::strcmp(argv[1], "--systematics") == 0) return fa::systematics(argv[2]);
    if (std::strcmp(argv[1], "--qa") == 0) return fa::runQA(fa::Settings::read(argv[2]));
    auto settings = fa::Settings::read(argv[2]);
    settings.correctionsOnly = std::strcmp(argv[1], "--corrections-only") == 0;
    return fa::run(settings);
  } catch (const std::exception &e) {
    std::cerr << "FinalAnalysis: " << e.what() << std::endl;
    return 1;
  }
}
