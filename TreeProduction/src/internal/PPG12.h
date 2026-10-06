#ifndef PHOTONJETTREE_INTERNAL_PPG12_H
#define PHOTONJETTREE_INTERNAL_PPG12_H

/**
 * @file PPG12.h
 * The established producer's fixed PPG12 photon-yield definitions, shared by
 * its DATA bookkeeping (PPG12.cc) and its simulation bookkeeping
 * (PPG12Simulation.cc). Values are the producer's own constants.
 */

#include "Types.h"

#include <array>
#include <cstdint>
#include <vector>

class TH1;

namespace photonjet::ppg12
{
// Reconstructed and truth ET bins of the yield histograms, GeV.
inline constexpr std::array<double, 12> kRecoEtBins{{10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36}};
inline constexpr std::array<double, 14> kTruthEtBins{{8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45}};

// Isolation range of the A/B/C/D regions, GeV.
inline constexpr double kIsolatedMin = -20.0;
inline constexpr double kNonIsolatedMax = 20.0;

// The ownership-window test: low <= value <= high (open above).
bool windowContains(double value, double low, double high);

// The truth-vertex weight at z: interpolated, 1 for a missing or sentinel
// vertex, 0 outside the histogram.
double truthVertexWeight(const TH1* weights, double z);

// Both ends of the response window: reconstructed and truth ET inside the
// outer edges of their yield bins.
bool inResponseWindow(double recoEt, double truthEt);

// The shower variables the yield selection reads, from the builder's primary
// representation of one candidate.
struct YieldShapes
{
  double probability, weta, wphi, et1, et2, et3, et4, e11OverE33, e32OverE35, npb, score;

  explicit YieldShapes(const PhotonRecord& photon);

  // Common selection: cluster probability and E11/E33 in open windows, a
  // positive width ratio, w_eta below 2 and an NPB score above 0.5.
  bool common() const;
};

// The PPG12 tightness tag of a candidate at the given ET.
enum class Tightness { Fail, Tight, NonTight, Neither };
Tightness tightness(const YieldShapes& shapes, double et);

// The R = 0.4 topocluster isolation of one candidate (the builder's stored
// value) and the audit stages of one producer lookup of it: 1 stored value
// usable, 2 accepted, 3 rejected, 4 recomputation fallback, 6 reconstructed
// vertex used.
struct IsolationLookup
{
  double value = 1e9;  // the producer's fail-safe: not isolated
  std::vector<std::int32_t> stages;
};
IsolationLookup isolationLookup(const PhotonRecord& photon, double eventVertexZ);
}  // namespace photonjet::ppg12

#endif  // PHOTONJETTREE_INTERNAL_PPG12_H
