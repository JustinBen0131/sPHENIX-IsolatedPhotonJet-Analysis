/* Unfolding: iterative Bayes (D'Agostini, as RooUnfoldBayes) with one fake cause, the per-photon normalisation, and the
   frozen iteration rule. */

/* The iteration scan, as the analysis note defines it, unfolds one problem for k = 1..12. For every k: S_data and
   S_response, the median over the reported bins of the central-68% half width of the unfolded data and response toys
   relative to the central value; D_k, the median symmetric relative change from k - 1; and the stability quantity
   Q_k = sqrt(S_data^2 + S_response^2 + D_k^2). The selected k is the one with the smallest Q_k (k >= 2, ties to the
   lower k). The photon yield and the joint (photon pT, xJ) unfolding are scanned and selected separately.
   The scan also records quality checks for every k: a finite result, a symmetric positive semi-definite toy
   covariance, no reported truth bin without efficiency, data and response toy success >= 95%, a data refold chi2/rank
   < 1.5 against the data plus response residual-toy covariance, simulation closure and refold chi2/ndf < 5, a negative
   fraction < 50%, a condition number < 1e8. They are saved with the scan and reported when the selected k fails one;
   they do not change the selection. */
#include "FinalAnalysis.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

namespace fa {

/* measured: reco bins; responseTruthXReco: matched simulation, truth x reco; misses: truth not reconstructed in the
   window; prior: truth distribution the iterations start from; fakes: reco entries without a truth cause in the window. */
TVectorD bayes(const TVectorD &measured, const TMatrixD &response, const TVectorD &misses, int iterations,
               const TVectorD &prior, const TVectorD &fakes, TVectorD *refolded) {
  const int nT = response.GetNrows(), nR = response.GetNcols();
  if (measured.GetNrows() != nR || misses.GetNrows() != nT || prior.GetNrows() != nT || fakes.GetNrows() != nR)
    throw std::runtime_error("unfolding dimensions differ");
  if (iterations < 1) throw std::runtime_error("iterations must be positive");
  double fakeTotal = 0;
  for (int j = 0; j < nR; ++j) fakeTotal += std::max(fakes[j], 0.0);
  const int nC = nT + (fakeTotal > 0 ? 1 : 0);  /* causes: truth bins and one fake cause */
  TMatrixD conditional(nC, nR);                  /* P(reco j | cause i) */
  TVectorD efficiency(nC), state(nC);
  for (int i = 0; i < nT; ++i) {
    double total = misses[i];
    for (int j = 0; j < nR; ++j) total += response(i, j);
    for (int j = 0; j < nR; ++j) conditional(i, j) = total > 0 ? response(i, j) / total : 0.0;
    for (int j = 0; j < nR; ++j) efficiency[i] += conditional(i, j);
    state[i] = prior[i];
  }
  if (nC > nT) {
    for (int j = 0; j < nR; ++j) conditional(nT, j) = std::max(fakes[j], 0.0) / fakeTotal;
    efficiency[nT] = 1.0;
    state[nT] = fakeTotal;
  }
  double sum = state.Sum();
  if (!(sum > 0)) throw std::runtime_error("empty prior");
  state *= 1.0 / sum;
  TVectorD unfolded(nC);
  for (int it = 0; it < iterations; ++it) {
    TVectorD denominator(nR);
    for (int j = 0; j < nR; ++j)
      for (int i = 0; i < nC; ++i) denominator[j] += conditional(i, j) * state[i];
    unfolded.Zero();
    for (int i = 0; i < nC; ++i) {
      if (!(efficiency[i] > 0)) continue;
      double x = 0;
      for (int j = 0; j < nR; ++j)
        if (denominator[j] > 0) x += conditional(i, j) * state[i] / denominator[j] * measured[j];
      unfolded[i] = x / efficiency[i];
    }
    if (!std::isfinite(unfolded.Sum())) throw std::runtime_error("empty signed estimate");
    /* the next iteration's prior: the estimate with negative causes at zero (a probability cannot be negative; with
       negative entries a reco denominator can change sign and the iteration collapses). The estimate itself stays signed. */
    double total = 0;
    for (int i = 0; i < nC; ++i) state[i] = std::max(unfolded[i], 0.0), total += state[i];
    if (!(total > 1e-15)) throw std::runtime_error("empty signed estimate");
    state *= 1.0 / total;
  }
  if (refolded) {  /* conditional^T unfolded: the reco spectrum the result predicts, fake cause included */
    refolded->ResizeTo(nR);
    refolded->Zero();
    for (int j = 0; j < nR; ++j)
      for (int i = 0; i < nC; ++i) (*refolded)[j] += conditional(i, j) * unfolded[i];
  }
  TVectorD out(nT);
  for (int i = 0; i < nT; ++i) out[i] = unfolded[i];
  return out;
}

TVectorD unfoldPhotons(const TVectorD &photons, const Window &w, int iterations) {
  return bayes(photons, w.photonResponse, w.photonMisses, iterations, w.photonTruth, w.photonFakes);
}

/* Pairs summed over photon pT, per photon and per unit xJ; every photon bin must hold a positive yield. */
TVectorD perPhoton(const TVectorD &pairCells, const TVectorD &photons, const Edges &xjEdges) {
  const int nX = xjEdges.size() - 1, nPt = photons.GetNrows();
  for (int i = 0; i < nPt; ++i)
    if (!(std::isfinite(photons[i]) && photons[i] > 0))
      throw std::runtime_error("the per-photon normalisation needs a positive yield in every photon bin");
  const double total = photons.Sum();
  TVectorD out(nX);
  for (int x = 0; x < nX; ++x) {
    double pairs = 0;
    for (int i = 0; i < nPt; ++i) pairs += pairCells[i * nX + x];
    out[x] = pairs / (total * (xjEdges[x + 1] - xjEdges[x]));
  }
  return out;
}

namespace {

const double kInf = std::numeric_limits<double>::infinity();

double median(std::vector<double> v) {  /* numpy.median; inf when empty */
  if (v.empty()) return kInf;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double quantile(std::vector<double> v, double q) {  /* numpy.quantile, linear */
  std::sort(v.begin(), v.end());
  const double at = q * (v.size() - 1);
  const size_t i = static_cast<size_t>(std::floor(at));
  return i + 1 < v.size() ? v[i] + (at - i) * (v[i + 1] - v[i]) : v.back();
}

TMatrixDSym sampleCovariance(const std::vector<TVectorD> &rows) {
  const int n = rows.size(), d = rows.front().GetNrows();
  TVectorD mean(d);
  for (const auto &r : rows) mean += r;
  mean *= 1.0 / n;
  TMatrixDSym c(d);
  for (const auto &r : rows)
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j) c(i, j) += (r[i] - mean[i]) * (r[j] - mean[j]);
  c *= 1.0 / (n - 1);
  return c;
}

double relativeHalfWidth(const TVectorD &central, const std::vector<TVectorD> &toys, const std::vector<bool> &mask) {
  if (toys.empty()) return kInf;
  std::vector<double> values;
  for (int i = 0; i < central.GetNrows(); ++i) {
    std::vector<double> column;
    for (const auto &t : toys) column.push_back(t[i]);
    const double half = 0.5 * (quantile(column, 0.84) - quantile(column, 0.16));
    if (mask[i] && std::isfinite(central[i]) && std::isfinite(half) && central[i] > 0) values.push_back(half / std::max(std::fabs(central[i]), 1e-12));
  }
  return median(values);
}

double deviation(const TVectorD &current, const TVectorD &previous, const std::vector<bool> &mask) {
  std::vector<double> values;
  for (int i = 0; i < current.GetNrows(); ++i) {
    const double change = current[i] - previous[i], scale = std::max(0.5 * (std::fabs(current[i]) + std::fabs(previous[i])), 1e-12);
    if (mask[i] && std::isfinite(change) && current[i] > 0 && previous[i] > 0) values.push_back(std::fabs(change) / scale);
  }
  return median(values);
}

double chi2ndf(const TVectorD &observed, const TVectorD &expected, const TVectorD &variance) {
  double sum = 0;
  int n = 0;
  bool close = true;
  for (int i = 0; i < observed.GetNrows(); ++i) {
    close = close && std::fabs(observed[i] - expected[i]) <= 1e-9 + 1e-9 * std::fabs(expected[i]);
    if (!std::isfinite(observed[i]) || !std::isfinite(expected[i]) || !(variance[i] > 0)) continue;
    sum += (observed[i] - expected[i]) * (observed[i] - expected[i]) / variance[i];
    ++n;
  }
  return n ? sum / n : close ? 0.0 : kInf;
}

/* Generalised chi2/rank of a residual on the supported subspace of a covariance (eigenmodes above 1e-10 of the largest);
   a residual outside that subspace fails closed. */
double covarianceChi2(const TVectorD &residual, const TMatrixDSym &covariance) {
  TMatrixDSym c(covariance);
  for (int i = 0; i < c.GetNrows(); ++i)
    for (int j = 0; j < i; ++j) c(i, j) = c(j, i) = 0.5 * (covariance(i, j) + covariance(j, i));
  for (int i = 0; i < c.GetNrows(); ++i) for (int j = 0; j < c.GetNrows(); ++j) if (!std::isfinite(c(i, j))) return kInf;
  const TMatrixDSymEigen eigen(c);
  const TVectorD lambda = eigen.GetEigenValues();
  const TMatrixD vectors = eigen.GetEigenVectors();
  double maximum = 0;
  for (int i = 0; i < lambda.GetNrows(); ++i) maximum = std::max(maximum, lambda[i]);
  double norm = 0;
  for (int i = 0; i < residual.GetNrows(); ++i) norm += residual[i] * residual[i];
  norm = std::sqrt(norm);
  if (maximum <= 0) return norm <= 1e-9 ? 0.0 : kInf;
  double value = 0, nullNorm = 0;
  int rank = 0;
  for (int k = 0; k < lambda.GetNrows(); ++k) {
    double coordinate = 0;
    for (int i = 0; i < residual.GetNrows(); ++i) coordinate += vectors(i, k) * residual[i];
    if (lambda[k] > 1e-10 * maximum) value += coordinate * coordinate / lambda[k], ++rank;
    else nullNorm += coordinate * coordinate;
  }
  if (rank == 0 || std::sqrt(nullNorm) > 1e-9 * std::max(norm, 1.0)) return kInf;
  return value / rank;
}

/* Largest over smallest nonzero singular value of P(reco | truth) on its occupied rows and columns (0: none). */
double conditionNumber(const TMatrixD &response, const TVectorD &misses) {
  std::vector<int> rows, cols;
  TMatrixD p(response);
  for (int i = 0; i < p.GetNrows(); ++i) {
    double total = misses[i];
    for (int j = 0; j < p.GetNcols(); ++j) total += response(i, j);
    for (int j = 0; j < p.GetNcols(); ++j) p(i, j) = total > 0 ? response(i, j) / total : 0.0;
  }
  for (int i = 0; i < p.GetNrows(); ++i) { bool any = false; for (int j = 0; j < p.GetNcols(); ++j) any = any || p(i, j) > 0; if (any) rows.push_back(i); }
  for (int j = 0; j < p.GetNcols(); ++j) { bool any = false; for (int i = 0; i < p.GetNrows(); ++i) any = any || p(i, j) > 0; if (any) cols.push_back(j); }
  if (rows.empty() || cols.empty()) return 0;
  const int nr = rows.size(), nc = cols.size(), n = std::max(nr, nc);
  TMatrixD a(nr, nc);
  for (int i = 0; i < nr; ++i) for (int j = 0; j < nc; ++j) a(i, j) = p(rows[i], cols[j]);
  const TMatrixD tall = nr >= nc ? a : TMatrixD(TMatrixD::kTransposed, a);  /* TDecompSVD needs rows >= columns */
  TDecompSVD svd(tall);
  if (!svd.Decompose()) return 0;
  const TVectorD sv = svd.GetSig();  /* descending */
  const double floor = n * std::numeric_limits<double>::epsilon() * sv[0];
  double smallest = 0;
  for (int i = 0; i < sv.GetNrows(); ++i) if (sv[i] > floor) smallest = sv[i];
  return smallest > 0 ? sv[0] / smallest : 0;
}

/* The masked toy covariance exists (two finite toys: then it is symmetric) and is positive semi-definite
   (lambda_min >= -1e-10 max |lambda|). */
std::pair<bool, bool> covarianceStatus(const std::vector<TVectorD> &toys, const std::vector<bool> &mask) {
  std::vector<int> keep;
  for (size_t i = 0; i < mask.size(); ++i) if (mask[i]) keep.push_back(i);
  if (toys.size() < 2 || keep.empty()) return {false, false};
  std::vector<TVectorD> masked;
  for (const auto &t : toys) {
    TVectorD m(keep.size());
    for (size_t i = 0; i < keep.size(); ++i) { m[i] = t[keep[i]]; if (!std::isfinite(m[i])) return {false, false}; }
    masked.push_back(m);
  }
  const TVectorD lambda = TMatrixDSymEigen(sampleCovariance(masked)).GetEigenValues();
  double scale = 0, minimum = 0;
  for (int i = 0; i < lambda.GetNrows(); ++i) scale = std::max(scale, std::fabs(lambda[i])), minimum = std::min(minimum, lambda[i]);
  return {true, minimum >= (scale > 0 ? -1e-10 * scale : 0.0)};
}

/* Unfolded values and refold residuals of toys; a failed or non-finite unfolding is dropped. */
void unfoldToys(const std::vector<std::pair<const Problem *, const TVectorD *>> &toys, int k, std::vector<TVectorD> &unfolded,
                std::vector<TVectorD> &residuals) {
  for (const auto &[p, measured] : toys) {
    try {
      TVectorD refolded;
      const TVectorD u = bayes(*measured, p->response, p->misses, k, p->prior, p->fakes, &refolded);
      bool finite = true;
      for (int i = 0; i < u.GetNrows(); ++i) finite = finite && std::isfinite(u[i]);
      if (!finite) continue;
      unfolded.push_back(u);
      residuals.push_back(*measured - refolded);
    } catch (const std::exception &) {
    }
  }
}

const double kToySuccessMin = 0.95, kRefoldMax = 1.5, kMcMax = 5.0, kNegativeMax = 0.5, kConditionMax = 1e8;

}  /* namespace */

Problem photonProblem(const Window &w) { return {w.photonResponse, w.photonMisses, w.photonTruth, w.photonFakes}; }
Problem pairProblem(const Window &w) { return {w.pairResponse, w.pairMisses, w.pairTruth, w.pairFakes}; }

/* The scan above on accepted toys, retaining the requested totals for success rates. selected = 0 means no k >= 2
   has a finite Q_k. */
IterationScan scanIterations(const Problem &p, const TVectorD &measured, const std::vector<TVectorD> &dataToys, int dataRequested,
                             const std::vector<ResponseToy> &responseToys, int responseRequested, const std::vector<bool> &mask) {
  if (responseRequested < 2) throw std::runtime_error("the iteration rule needs at least two response toys");
  IterationScan scan;
  const int nT = p.response.GetNrows();
  int zeroEfficiency = 0;
  for (int i = 0; i < nT; ++i) {
    double matched = 0;
    for (int j = 0; j < p.response.GetNcols(); ++j) matched += p.response(i, j);
    const double total = matched + p.misses[i];
    if (mask[i] && total > 0 && matched / total == 0) ++zeroEfficiency;
  }
  const double condition = conditionNumber(p.response, p.misses);
  TVectorD pseudo(p.response.GetNcols()), pseudoVariance(p.response.GetNcols());
  for (int j = 0; j < pseudo.GetNrows(); ++j) {
    for (int i = 0; i < nT; ++i) pseudo[j] += p.response(i, j);
    pseudo[j] += p.fakes[j];
    pseudoVariance[j] = std::max(pseudo[j], 1e-12);
  }
  std::vector<std::pair<const Problem *, const TVectorD *>> data, response;
  for (const auto &t : dataToys) data.emplace_back(&p, &t);
  for (const auto &t : responseToys) response.emplace_back(&t.problem, &t.measured);
  for (int k = 1; k <= 12; ++k) {
    ScanRow row;
    row.k = k;
    TVectorD refolded;
    const TVectorD central = bayes(measured, p.response, p.misses, k, p.prior, p.fakes, &refolded);
    std::vector<TVectorD> dataUnfolded, dataResiduals, responseUnfolded, responseResiduals;
    unfoldToys(data, k, dataUnfolded, dataResiduals);
    unfoldToys(response, k, responseUnfolded, responseResiduals);
    row.dataSuccess = double(dataUnfolded.size()) / std::max(1, dataRequested);
    row.responseSuccess = double(responseUnfolded.size()) / responseRequested;
    row.dataStat = relativeHalfWidth(central, dataUnfolded, mask);
    row.responseStat = relativeHalfWidth(central, responseUnfolded, mask);
    row.deviation = k == 1 ? std::nan("") : deviation(central, scan.unfolded.back(), mask);
    row.objective = k == 1 ? std::nan("") : std::sqrt(row.dataStat * row.dataStat + row.responseStat * row.responseStat + row.deviation * row.deviation);
    /* the data and response residual-toy covariances are independent and add */
    TMatrixDSym covariance(measured.GetNrows());
    bool responseCovariance = false;
    if (dataResiduals.size() >= 2) covariance += sampleCovariance(dataResiduals);
    if (responseResiduals.size() >= 2) {
      const TMatrixDSym component = sampleCovariance(responseResiduals);
      covariance += component;
      responseCovariance = component.E2Norm() > 0;
    }
    row.refold = responseCovariance ? covarianceChi2(measured - refolded, covariance) : kInf;
    TVectorD closureRefolded;
    const TVectorD closure = bayes(pseudo, p.response, p.misses, k, p.prior, p.fakes, &closureRefolded);
    TVectorD closureVariance(nT);
    for (int i = 0; i < nT; ++i) closureVariance[i] = std::max(p.prior[i] + std::fabs(closure[i]), 1e-12);
    row.mcClosure = chi2ndf(p.prior, closure, closureVariance);
    row.mcRefold = chi2ndf(pseudo, closureRefolded, pseudoVariance);
    double negative = 0, total = 0;
    bool finite = true;
    for (int i = 0; i < nT; ++i) {
      finite = finite && std::isfinite(central[i]);
      if (!mask[i]) continue;
      total += std::fabs(central[i]);
      if (central[i] < 0) negative += std::fabs(central[i]);
    }
    row.negativeFraction = negative / std::max(total, 1e-15);
    const auto [symmetric, psd] = covarianceStatus(dataUnfolded, mask);
    const std::pair<const char *, bool> gates[] = {
        {"finite_unfolded", finite}, {"covariance_symmetric", symmetric}, {"covariance_psd", psd},
        {"reported_zero_efficiency_states", zeroEfficiency == 0}, {"data_toy_success", row.dataSuccess >= kToySuccessMin},
        {"response_toy_success", row.responseSuccess >= kToySuccessMin}, {"response_residual_covariance", responseCovariance},
        {"data_refold", row.refold < kRefoldMax}, {"mc_closure", row.mcClosure < kMcMax}, {"mc_refold", row.mcRefold < kMcMax},
        {"negative_fraction", row.negativeFraction < kNegativeMax}, {"conditioning", condition > 0 && condition < kConditionMax}};
    for (const auto &[name, passed] : gates) if (!passed) row.failures += (row.failures.empty() ? "" : "|") + std::string(name);
    row.passes = row.failures.empty();
    scan.rows.push_back(row);
    scan.unfolded.push_back(central);
    scan.refolded.push_back(refolded);
  }
  double best = kInf;
  for (const auto &row : scan.rows)
    if (row.k >= 2 && std::isfinite(row.objective) && row.objective < best) best = row.objective, scan.selected = row.k;
  return scan;
}

}  /* namespace fa */
