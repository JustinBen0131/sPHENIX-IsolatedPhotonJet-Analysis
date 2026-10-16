/* Explicit statistical errors of saved analysis histograms.
   Counts supply their measured sum of squared weights. Derived quantities
   supply propagated variances; tables and covariance matrices use explicit
   zero bin errors. The publication guard rejects ROOT's implicit sqrt(content)
   convention, which has no meaning for a ratio or a matrix element. */
#ifndef FINAL_ANALYSIS_HISTOGRAM_ERRORS_H
#define FINAL_ANALYSIS_HISTOGRAM_ERRORS_H

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include "TDirectory.h"
#include "TH1.h"
#include "TKey.h"

namespace fa {

struct RatioEstimate { double value, variance; };

/* Only for a source contract that explicitly records unweighted event counts.
   Integer contents alone do not establish that contract. This upgrades the
   legacy SourceStudiesV2 centrality target without inventing ratio errors or
   rerunning a DATA reduction whose source records were deliberately compacted. */
inline void materializeUnitCountErrors(TH1& histogram) {
  if (histogram.GetSumw2N()) return;
  for (int bin=0;bin<histogram.GetNcells();++bin) {
    const double n=histogram.GetBinContent(bin);
    if (!std::isfinite(n) || n<0 || n!=std::trunc(n))
      throw std::runtime_error("unit-count variance contract violated: " + std::string(histogram.GetName()));
  }
  histogram.Sumw2();
  for (int bin=0;bin<histogram.GetNcells();++bin)
    histogram.SetBinError(bin,std::sqrt(histogram.GetBinContent(bin)));
}

/* Delta-method variance of A*D/(B*C) for independent weighted cells.
   This is a count-statistics model, not a calibration or nonclosure nuisance.
   A zero numerator can have a nonzero variance; only the denominators must
   be positive. Unknown or invalid input moments remain explicitly undefined. */
inline RatioEstimate abcdRatio(const std::array<double, 4>& n,
                               const std::array<double, 4>& sumw2) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (!(n[1] > 0 && n[2] > 0)) return {nan, nan};
  for (unsigned q = 0; q < 4; ++q)
    if (!std::isfinite(n[q]) || !(sumw2[q] >= 0) || !std::isfinite(sumw2[q])) return {nan, nan};
  const double value = n[0] * n[3] / (n[1] * n[2]);
  const std::array<double, 4> derivative = {n[3] / (n[1] * n[2]), -value / n[1],
                                          -value / n[2], n[0] / (n[1] * n[2])};
  double variance = 0;
  for (unsigned q = 0; q < 4; ++q) variance += derivative[q] * derivative[q] * sumw2[q];
  return {value, variance};
}

/* A leakage fraction uses region A as reference. A/A is the same integral
   twice and is exactly one, with zero variance; it is not two independent
   Poisson measurements. Other regions use the declared weighted-cell model. */
inline RatioEstimate signalLeakageRatio(double numerator, double reference,
                                        double numeratorVariance, double referenceVariance,
                                        bool sameIntegral) {
  if (!(reference > 0)) {
    const double nan=std::numeric_limits<double>::quiet_NaN(); return {nan,nan};
  }
  if (sameIntegral) {
    if (numerator!=reference || numeratorVariance!=referenceVariance)
      throw std::runtime_error("shared leakage reference has different moments");
    return {1.,0.};
  }
  return abcdRatio({numerator,reference,1.,1.},{numeratorVariance,referenceVariance,0.,0.});
}

/* Inspect serialized keys before an output is promoted from its partial name.
   Include under/overflow, negative and undefined contents. A table's explicit
   zero Sumw2 is valid; merely allocating a ROOT histogram is not sufficient. */
inline void requireExplicitHistogramErrors(TDirectory* directory, const std::string& prefix = "") {
  // Newly created child directories may not have a parent key until Close.
  // Inspect their saved histogram keys through the live directory hierarchy.
  std::set<std::string> visited;
  TIter live(directory->GetList());
  while (auto* object = live()) if (auto* child = dynamic_cast<TDirectory*>(object)) {
    visited.insert(child->GetName());
    requireExplicitHistogramErrors(child, prefix + child->GetName() + "/");
  }
  TIter next(directory->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next())) {
    if (visited.count(key->GetName())) continue;
    std::unique_ptr<TObject> object(key->ReadObj());
    const std::string path = prefix + key->GetName();
    if (auto* child = dynamic_cast<TDirectory*>(object.get())) {
      requireExplicitHistogramErrors(child, path + "/");
    } else if (auto* histogram = dynamic_cast<TH1*>(object.get())) {
      if (histogram->GetSumw2N() != 0) continue;
      for (int bin = 0; bin < histogram->GetNcells(); ++bin)
        if (histogram->GetBinContent(bin) != 0)
          throw std::runtime_error("histogram has implicit sqrt(content) errors: " + path);
    }
  }
}

} // namespace fa
#endif
