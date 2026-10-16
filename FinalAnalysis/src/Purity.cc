/* Purity: the leakage-aware ABCD background with its correlation R_bkg, and the combinatoric-jet subtraction. */

/* ABCD regions: A tight isolated, B tight non-isolated, C non-tight isolated, D non-tight non-isolated. With the signal
   s in A and the prompt-photon leakage l_r = s_r / s_A of the simulation,
   (A - s) = R (B - l_B s)(C - l_C s) / (D - l_D s),
   where R = 1 when the background factorises between identification and isolation (ATLAS, PRC 93 (2016) 034914,
   Eq. 4). The region-A recoil jets are corrected with the transfer k = (A - s)/(C - l_C s):
   corrected = (A_pairs - k C_pairs) / (1 - k l_C),
   then the simulated combinatoric template, scaled to the corrected photons, is subtracted. A purity scale (the
   non-closure systematic) multiplies the signal and re-derives k from it; an efficiency scale (DATA/simulation photon
   efficiency) divides the corrected photons and jets of its bin. Both are 1 in the nominal. */

/* The data-driven estimate measures the level of the combinatoric jets in the data, in a Delta-phi sideband away from
   the photon and from its recoil (Combinatoric.Sideband, nominal 0.8 < Delta phi < 2.5), per photon bin and xJ bin.
   Combinatoric.Shape simulation (nominal) carries that level into the recoil window with the simulation:
     combinatoric jets in the window = (DATA jets in the sideband) x T,
     T = (simulated jets without a truth jet nearby in the window) / (every simulated jet in the sideband),
   both of matched region-A photons (Window::combinatoric and Window::sidebandJets). T is not the ratio of the two
   Delta-phi ranges: where the real recoil jet stands, fewer jets of the underlying event are found, so in the embedded
   simulation of central Au+Au the density of jets without a truth jet nearby falls to 0.6-0.7 of the sideband inside
   the window at low xJ. An estimate that takes the sideband density there removes 1.2-1.4 times the true combinatoric
   jets and most of the real recoil jets at low xJ. With T the estimate returns the simulated combinatoric jets exactly
   when it is applied to the simulation, in p+p (where T is near zero) as in Au+Au, and the response toys carry its
   statistical uncertainty. The same T is used for the region A and region C photons.
   Shape flow and flat (combinatoricFromDphi) are the earlier estimate, which assumes the sideband density under the
   recoil peak: jets that do not come
   from the photon's hard scattering populate every Delta phi with a flow-modulated shape. The Delta-phi distribution of
   the jets at |Delta eta| > 0.8, free of the photon's own scattering, is fitted over 0 < Delta phi < 2.5 with
   f = a0 + a2 cos(2 Delta phi) + a3 cos(3 Delta phi)   (flow; flat: a constant)
   by weighted least squares on the bin integrals (weight 1 / max(variance, 1)); per xJ bin it is normalised to the
   jets at every Delta eta in the sideband 0.8 < Delta phi < 2.5, and its integral over the recoil window
   Delta phi > 7 pi / 8 is the combinatoric count, removed from the region A and C recoil jets before the purity
   correction. */
#include "FinalAnalysis.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "TDecompLU.h"

namespace fa {

/* Both arrays are ordered A, B, C, D; leakage contains signal counts and is normalised to A here.
   ok=false marks an unphysical or absent root, even if the returned diagnostic estimate is finite. */
ABCDSolution solveLeakageABCD(const double n[4], const double leakage[4], double r) {
  ABCDSolution out;
  if (leakage[0] <= 0 || !(r > 0)) return out;
  double l[4];
  for (int i = 0; i < 4; ++i) l[i] = leakage[i] / leakage[0];
  const double a = l[3] - r * l[1] * l[2];
  const double b = -n[0] * l[3] - n[3] + r * n[1] * l[2] + r * n[2] * l[1];
  const double c = n[0] * n[3] - r * n[1] * n[2];
  double upper = n[0];
  for (int i = 1; i < 4; ++i)
    if (l[i] > 0) upper = std::min(upper, n[i] / l[i]);
  const double naive = n[3] > 0 ? std::min(std::max(n[0] - r * n[1] * n[2] / n[3], 0.0), upper) : upper;
  std::vector<double> roots;
  if (std::fabs(a) < 1e-14) {
    if (std::fabs(b) > 1e-14) roots.push_back(-c / b);
  } else {
    const double disc = b * b - 4 * a * c;
    if (disc >= 0) {
      roots.push_back((-b - std::sqrt(disc)) / (2 * a));
      roots.push_back((-b + std::sqrt(disc)) / (2 * a));
    }
  }
  double signal = naive;
  bool physical = false;
  for (double x : roots)
    if (x >= -1e-9 && x <= upper + 1e-9 && (!physical || std::fabs(x - naive) < std::fabs(signal - naive))) {
      signal = x;
      physical = true;
    }
  const double backgroundA = n[0] - l[0] * signal, backgroundC = n[2] - l[2] * signal;
  out.signal = std::max(0.0, signal);
  out.purity = n[0] > 0 ? std::min(1.0, std::max(0.0, signal / n[0])) : 0.0;
  out.transfer = backgroundC > 0 ? std::max(0.0, backgroundA / backgroundC) : 0.0;
  out.ok = physical;
  return out;
}

Corrected correct(const Measured &m, const Window &w, const SystemSettings &sys) {
  const int nPt = m.abcd.GetNcols(), nX = m.pairsA.GetNcols();
  Corrected c;
  c.photons.ResizeTo(nPt); c.purity.ResizeTo(nPt); c.transfer.ResizeTo(nPt);
  c.physical.assign(nPt, false);
  for (TMatrixD *x : {&c.pairs, &c.pairsVariance, &c.scaledC, &c.scaledCVariance, &c.purityCorrected, &c.purityCorrectedVariance,
                      &c.combinatoric})
    x->ResizeTo(nPt, nX);
  for (int pt = 0; pt < nPt; ++pt) {
    const double counts[4] = {m.abcd(0, pt), m.abcd(1, pt), m.abcd(2, pt), m.abcd(3, pt)};
    const double leak[4] = {w.leakage(0, pt), w.leakage(1, pt), w.leakage(2, pt), w.leakage(3, pt)};
    const double r = sys.backgroundCorrelation.empty() ? 1.0 : sys.backgroundCorrelation[pt];
    ABCDSolution s = solveLeakageABCD(counts, leak, r);
    c.physical[pt] = s.ok;
    const double relativeC = leak[2] / leak[0];
    /* a purity scale multiplies the signal; the transfer is re-derived so the corrected jets integrate to it */
    const double purityScale = sys.purityScale.empty() ? 1.0 : sys.purityScale[pt];
    if (purityScale != 1.0) {
      const double signal = purityScale * s.signal, backgroundC = counts[2] - relativeC * signal;
      if (!(signal > 0 && signal < counts[0]) || !(backgroundC > 0))
        throw std::runtime_error(sys.name + ": the purity scale leaves no physical ABCD solution in photon bin " + std::to_string(pt));
      s.signal = signal, s.purity = signal / counts[0], s.transfer = (counts[0] - signal) / backgroundC;
    }
    const double denominator = 1.0 - s.transfer * relativeC;
    if (!(denominator > 1e-6)) throw std::runtime_error(sys.name + ": nonphysical region-C leakage denominator");
    c.photons[pt] = s.signal;
    c.purity[pt] = s.purity;
    c.transfer[pt] = s.transfer;
    const double normalisation = w.leakage(0, pt);  /* truth-matched region-A photons of the template */
    for (int x = 0; x < nX; ++x) {
      double jetsA = m.pairsA(pt, x), varianceA = m.pairsAVariance(pt, x), jetsC = m.pairsC(pt, x), varianceC = m.pairsCVariance(pt, x);
      if (sys.subtractSideband) {
        /* the sideband jets of the data times the simulated transfer (see the header); the variance of a count is the count */
        const double transfer = w.sidebandJets(pt, x) > 0 ? w.combinatoric(pt, x) / w.sidebandJets(pt, x) : 0.0;
        jetsA -= transfer * m.sidebandA(pt, x), varianceA += transfer * transfer * m.sidebandA(pt, x);
        jetsC -= transfer * m.sidebandC(pt, x), varianceC += transfer * transfer * m.sidebandC(pt, x);
      }
      c.scaledC(pt, x) = s.transfer * jetsC;
      c.scaledCVariance(pt, x) = s.transfer * s.transfer * varianceC;
      c.pairs(pt, x) = (jetsA - s.transfer * jetsC) / denominator;
      c.pairsVariance(pt, x) = (varianceA + s.transfer * s.transfer * varianceC) / (denominator * denominator);
      c.purityCorrected(pt, x) = c.pairs(pt, x);
      c.purityCorrectedVariance(pt, x) = c.pairsVariance(pt, x);
      if (sys.subtractTemplate && normalisation > 0) {
        const double scale = s.signal / normalisation;
        const double slope = -s.signal * w.combinatoric(pt, x) / (normalisation * normalisation);
        c.combinatoric(pt, x) = w.combinatoric(pt, x) * scale;
        c.pairs(pt, x) -= c.combinatoric(pt, x);
        c.pairsVariance(pt, x) += std::max(0.0, scale * scale * w.combinatoricVariance(pt, x) +
                                                    slope * slope * w.leakageVariance(0, pt));
      }
    }
    /* an efficiency scale divides the bin's photons and jets as the last step (the displays stay before it) */
    const double efficiency = sys.efficiencyScale.empty() ? 1.0 : sys.efficiencyScale[pt];
    if (efficiency != 1.0) {
      c.photons[pt] /= efficiency;
      for (int x = 0; x < nX; ++x) c.pairs(pt, x) /= efficiency, c.pairsVariance(pt, x) /= efficiency * efficiency;
    }
  }
  return c;
}

namespace {

const double kWindowLow = 7 * M_PI / 8, kWindowHigh = M_PI;

double integral(const double a[3], double lo, double hi) {
  return a[0] * (hi - lo) + a[1] * (std::sin(2 * hi) - std::sin(2 * lo)) / 2 + a[2] * (std::sin(3 * hi) - std::sin(3 * lo)) / 3;
}

bool inside(const Edges &e, int d, double lo, double hi) { return e[d] >= lo - 1e-12 && e[d + 1] <= hi + 1e-12; }

/* The flow form fitted to far-side counts over 0 < Delta phi < 2.5; returns the window / sideband integral ratio. */
double windowOverSideband(const std::vector<double> &counts, const std::vector<double> &variance, const std::string &shape,
                          double sidebandLow, double sidebandHigh) {
  double a[3] = {1, 0, 0};
  if (shape == "flow") {
    const Edges &e = dphiEdges();
    TMatrixD normal(3, 3);
    TVectorD rhs(3);
    double total = 0;
    for (size_t d = 0; d + 1 < e.size(); ++d) {
      if (!inside(e, d, 0.0, sidebandHigh)) continue;
      const double basis[3] = {e[d + 1] - e[d], (std::sin(2 * e[d + 1]) - std::sin(2 * e[d])) / 2,
                               (std::sin(3 * e[d + 1]) - std::sin(3 * e[d])) / 3};
      const double weight = 1.0 / std::max(variance[d], 1.0);
      for (int i = 0; i < 3; ++i) {
        rhs[i] += weight * basis[i] * counts[d];
        for (int j = 0; j < 3; ++j) normal(i, j) += weight * basis[i] * basis[j];
      }
      total += counts[d];
    }
    if (!(total > 0)) throw std::runtime_error("no jets in the Delta-phi fit range");
    TDecompLU lu(normal);
    if (!lu.Solve(rhs)) throw std::runtime_error("singular Delta-phi fit");
    for (int i = 0; i < 3; ++i) a[i] = rhs[i];
  }
  const double sideband = integral(a, sidebandLow, sidebandHigh);
  if (!(sideband > 0)) throw std::runtime_error("the fitted shape is not positive in the sideband");
  return integral(a, kWindowLow, kWindowHigh) / sideband;
}

}  /* namespace */

const Edges &dphiEdges() {  /* pi/32 steps with 0.8 and 2.5 rad added (TreeToHists FillPhotonJet.C) */
  static const Edges edges = [] {
    Edges e;
    for (int i = 0; i <= 32; ++i) e.push_back(i * M_PI / 32);
    e.push_back(0.8);
    e.push_back(2.5);
    std::sort(e.begin(), e.end());
    return e;
  }();
  return edges;
}

Background combinatoricFromDphi(const std::vector<double> dphi[2], const std::vector<double> variance[2], int nPt, int nXj,
                                const std::string &shape, const std::string &scope, double sidebandLow, double sidebandHigh) {
  const Edges &e = dphiEdges();
  const int nD = e.size() - 1;
  auto at = [&](int p, int x, int d, int side) { return ((p * nXj + x) * nD + d) * 2 + side; };
  Background out;
  for (int r = 0; r < 2; ++r) {
    if (static_cast<int>(dphi[r].size()) != nPt * nXj * nD * 2 || variance[r].size() != dphi[r].size())
      throw std::runtime_error("Delta-phi counts do not match the photon, xJ and Delta-phi bins");
    out.counts[r].ResizeTo(nPt, nXj);
    out.variance[r].ResizeTo(nPt, nXj);
    for (int p = 0; p < nPt; ++p) {
      std::vector<double> sideband(nXj), sidebandVariance(nXj);
      for (int x = 0; x < nXj; ++x)
        for (int d = 0; d < nD; ++d)
          if (inside(e, d, sidebandLow, sidebandHigh))
            for (int side = 0; side < 2; ++side) {
              sideband[x] += dphi[r][at(p, x, d, side)];
              sidebandVariance[x] += variance[r][at(p, x, d, side)];
            }
      for (int x = 0; x < nXj; ++x) {
        if (scope == "pt" && x > 0) break;
        std::vector<double> far(nD), farVariance(nD);
        for (int y = scope == "pt" ? 0 : x; y < (scope == "pt" ? nXj : x + 1); ++y)
          for (int d = 0; d < nD; ++d) {
            far[d] += dphi[r][at(p, y, d, 1)];
            farVariance[d] += variance[r][at(p, y, d, 1)];
          }
        const double fraction = windowOverSideband(far, farVariance, shape, sidebandLow, sidebandHigh);
        for (int y = scope == "pt" ? 0 : x; y < (scope == "pt" ? nXj : x + 1); ++y) {
          out.counts[r](p, y) = fraction * sideband[y];
          out.variance[r](p, y) = fraction * fraction * sidebandVariance[y];
        }
      }
    }
  }
  return out;
}

}  /* namespace fa */
