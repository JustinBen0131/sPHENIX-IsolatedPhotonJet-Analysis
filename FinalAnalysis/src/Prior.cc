/* Prior: the data-driven prior of the joint unfolding. */

/* The simulated truth is reweighted by a smooth function of xJ,
   w(xJ) = exp(t0 + t1 (u - 0.8) [+ t2 (u - 0.8)^2]),  u = truth xJ bin centre, held constant outside the reported
   range, until the folded simulation, scaled to the measured photons of each reconstructed photon bin, describes the
   corrected pairs (generalized least squares with the data-toy covariance). P(reco | truth) is untouched: the weights
   change the prior, the boundary migration and, with the "scaled" model, the detector fakes in proportion to the
   signal. The unfolding starts from w^power. */
#include "FinalAnalysis.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "TDecompChol.h"
#include "TMatrixDSym.h"

namespace fa {
namespace {

bool reported(const Edges &xj, int x, double low, double high) { return xj[x] >= low - 1e-9 && xj[x + 1] <= high + 1e-9; }

/* Matched plus boundary-migrated simulated pairs in every reco cell for truth weights w. */
TVectorD signal(const Window &win, const TVectorD &w) {
  const int nCells = win.pairResponse.GetNcols(), nX = w.GetNrows();
  TVectorD out(win.boundaryFixed);
  for (int t = 0; t < win.pairResponse.GetNrows(); ++t)
    for (int r = 0; r < nCells; ++r) out[r] += win.pairResponse(t, r) * w[t % nX];
  for (int x = 0; x < nX; ++x)
    for (int r = 0; r < nCells; ++r) out[r] += w[x] * win.boundaryByXj(x, r);
  return out;
}

TVectorD detectorFakes(const Window &win, const TVectorD &w, const std::string &model) {
  if (model == "fixed") return win.detector;
  const TVectorD one = [&] { TVectorD u(w.GetNrows()); for (int i = 0; i < u.GetNrows(); ++i) u[i] = 1; return u; }();
  const TVectorD nominal = signal(win, one), now = signal(win, w);
  TVectorD out(win.detector.GetNrows());
  for (int r = 0; r < out.GetNrows(); ++r) out[r] = nominal[r] > 0 ? win.detector[r] / nominal[r] * now[r] : win.detector[r];
  return out;
}

TVectorD folded(const Window &win, const TVectorD &w, const std::string &model) {
  return signal(win, w) + detectorFakes(win, w, model);
}

/* Levenberg-Marquardt on whitened residuals: the minimum of the generalized chi2. */
struct FitProblem {
  const Window &win;
  const TVectorD &cells, &photons;
  const Settings &s;
  double reportedLow;
  std::vector<int> fitted;       /* fitted reco cells */
  TMatrixD lower;                /* Cholesky factor L of the covariance of the fitted cells (V = L L^T) */
  TVectorD residual(const TVectorD &theta) const {
    const TVectorD w = priorWeights(theta, s.xjEdges, reportedLow, s.reportedHigh);
    const TVectorD f = folded(win, w, s.fakeModel);
    const int nX = s.xjEdges.size() - 1, n = fitted.size();
    TVectorD d(n), y(n);
    for (int i = 0; i < n; ++i) {
      const int r = fitted[i], pt = r / nX;
      d[i] = cells[r] - photons[pt] / win.simPhotons[pt] * f[r];
    }
    for (int i = 0; i < n; ++i) {  /* y = L^-1 d */
      double v = d[i];
      for (int j = 0; j < i; ++j) v -= lower(i, j) * y[j];
      y[i] = v / lower(i, i);
    }
    return y;
  }
};

}  /* namespace */

TVectorD priorWeights(const TVectorD &theta, const Edges &xj, double low, double high) {
  const int nX = xj.size() - 1;
  int first = -1, last = -1;
  for (int x = 0; x < nX; ++x)
    if (reported(xj, x, low, high)) { if (first < 0) first = x; last = x; }
  const double uLow = 0.5 * (xj[first] + xj[first + 1]), uHigh = 0.5 * (xj[last] + xj[last + 1]);
  TVectorD w(nX);
  for (int x = 0; x < nX; ++x) {
    const double u = std::min(std::max(0.5 * (xj[x] + xj[x + 1]), uLow), uHigh) - 0.8;
    double t = 0;
    for (int p = 0; p < theta.GetNrows(); ++p) t += theta[p] * std::pow(u, p);
    w[x] = std::exp(t);
  }
  return w;
}

PriorFit fitPrior(const Window &win, const TVectorD &cells, const TVectorD &photons, const TMatrixD &covariance,
                  const Settings &s, const SystemSettings &sys) {
  const int nX = s.xjEdges.size() - 1, nCells = cells.GetNrows(), nPar = s.priorForm == "exp2" ? 3 : 2;
  FitProblem prob{win, cells, photons, s, sys.reportedLow, {}, TMatrixD()};
  for (int r = 0; r < nCells; ++r)
    if (s.fitRange == "all" ? covariance(r, r) > 0 : reported(s.xjEdges, r % nX, sys.reportedLow, s.reportedHigh))
      prob.fitted.push_back(r);
  const int n = prob.fitted.size();
  /* Fewer populated cells than parameters (e.g. no corrected pairs with a toy variance) has no fit: refuse it here,
     before an empty covariance reaches the matrix algebra. */
  if (n < nPar)
    throw std::runtime_error("the prior fit has " + std::to_string(n) + " corrected-pair cells with a data variance; it needs at least " +
                             std::to_string(nPar));
  TMatrixDSym v(n);
  double trace = 0;
  for (int i = 0; i < n; ++i) trace += covariance(prob.fitted[i], prob.fitted[i]);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) v(i, j) = covariance(prob.fitted[i], prob.fitted[j]) + (i == j ? 1e-12 * trace / n : 0.0);
  TDecompChol chol(v);
  PriorFit out;
  if (!chol.Decompose()) return out;
  prob.lower.ResizeTo(n, n);
  prob.lower.Transpose(chol.GetU());
  TVectorD theta(nPar);
  TVectorD r = prob.residual(theta);
  double chi2 = r.Norm2Sqr(), lambda = 1e-3;
  for (int it = 0; it < 500; ++it) {
    TMatrixD J(n, nPar);
    for (int p = 0; p < nPar; ++p) {
      const double h = 1e-7 * std::max(1.0, std::fabs(theta[p]));
      TVectorD up(theta), down(theta);
      up[p] += h; down[p] -= h;
      const TVectorD ru = prob.residual(up), rd = prob.residual(down);
      for (int i = 0; i < n; ++i) J(i, p) = (ru[i] - rd[i]) / (2 * h);
    }
    const TMatrixD Jt(TMatrixD::kTransposed, J);
    const TVectorD g = Jt * r;
    const TMatrixD JtJ(Jt, TMatrixD::kMult, J);
    bool improved = false;
    for (int tries = 0; tries < 30 && !improved; ++tries) {
      TMatrixD a(JtJ);
      for (int p = 0; p < nPar; ++p) a(p, p) *= 1 + lambda;
      TVectorD step = -1.0 * g;
      TDecompChol solve(TMatrixDSym(nPar, a.GetMatrixArray()));
      if (!solve.Decompose()) { lambda *= 10; continue; }
      solve.Solve(step);
      const TVectorD next = theta + step;
      const TVectorD rn = prob.residual(next);
      if (rn.Norm2Sqr() < chi2) {
        const double change = chi2 - rn.Norm2Sqr();
        theta = next; r = rn; improved = true;
        lambda = std::max(lambda / 10, 1e-12);
        if (change <= 1e-14 * std::max(1.0, chi2)) { chi2 = rn.Norm2Sqr(); out.ok = true; break; }
        chi2 = rn.Norm2Sqr();
      } else {
        lambda *= 10;
      }
    }
    if (out.ok) break;
    if (!improved) { out.ok = true; break; }  /* no step lowers chi2: the minimum within numerical precision */
  }
  out.theta.ResizeTo(nPar);
  out.theta = theta;
  out.weights.ResizeTo(nX);
  out.weights = priorWeights(theta, s.xjEdges, sys.reportedLow, s.reportedHigh);
  out.chi2 = chi2;
  out.ndf = n - nPar;
  return out;
}

TVectorD unfoldPairsWithPrior(const TVectorD &cells, const Window &w, int iterations, const TVectorD &weights,
                              const std::string &fakeModel, TVectorD *refolded) {
  const int nX = weights.GetNrows();
  TVectorD prior(w.pairTruth), fakes = detectorFakes(w, weights, fakeModel) + w.boundaryFixed;
  for (int t = 0; t < prior.GetNrows(); ++t) prior[t] *= weights[t % nX];
  for (int x = 0; x < nX; ++x)
    for (int r = 0; r < fakes.GetNrows(); ++r) fakes[r] += weights[x] * w.boundaryByXj(x, r);
  return bayes(cells, w.pairResponse, w.pairMisses, iterations, prior, fakes, refolded);
}

}  /* namespace fa */
