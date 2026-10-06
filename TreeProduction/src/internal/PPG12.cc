/**
 * @file PPG12.cc
 * The established pp producer's PPG12 photon-yield bookkeeping.
 *
 * While it produced trees, the established pp producer also filled its own
 * PPG12 photon-yield histograms and logged an audit occurrence for the fills
 * it checked: the vertex/period contract of each classified event, the
 * simulation sample ownership and, in DATA, for every photon in the analysis
 * ET bins, its isolation lookups, common and tight selection and
 * isolation-stack entry. The tree keeps those occurrences (pdiag_*) in their
 * original order. This file reproduces them from the captured measurements;
 * it fills no histogram and changes no event, object or label. The
 * simulation histogram bookkeeping after the weight configuration is
 * PPG12Simulation.cc.
 */
#include "../PhotonJetTree.h"
#include "PPG12.h"
#include "Production.h"

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace
{
// The PPG12 period label of a DATA run, "other" outside every period.
std::string periodLabel(const photonjet::EventClassificationConfig& c, const int run)
{
  for (const auto& period : c.periods)
    if (run >= period.runMin && run < period.runMaxExclusive) return period.key;
  return "other";
}


// DATA yield weight: the inverse bit-30 trigger efficiency, a Gumbel CDF in ET.
double dataTriggerEfficiencyWeight(const double et)
{
  if (!std::isfinite(et)) return 0.0;
  const double efficiency = std::exp(-std::exp(-(et + 2.04) / 3.33));
  return std::isfinite(efficiency) && efficiency > 1e-6 ? 1.0 / efficiency : 0.0;
}

// The producer's fixed shower windows. Common selection: cluster probability
// and E11/E33 in open windows, a positive width ratio, w_eta below 2 and an
// NPB score above 0.5. Tight and non-tight windows; every other window is the
// open unit interval and the tight score also stays below 1.
constexpr double kCommonE11E33Max = 0.98;
constexpr double kCommonWetaMax = 2.0;
constexpr double kCommonNpbMin = 0.5;
constexpr double kTightEt1Min = 0.5;
constexpr double kNonTightEt1Min = 0.6;
constexpr double kTightE32E35Min = 0.8;
constexpr double kTightScoreMax = 1.0;

bool inOpenInterval(const double x, const double low, const double high)
{
  return std::isfinite(x) && x > low && x < high;
}

double ratio(const double numerator, const double denominator)
{
  return std::isfinite(numerator) && std::isfinite(denominator) && denominator > 0.0
             ? numerator / denominator
             : std::numeric_limits<double>::quiet_NaN();
}
}  // namespace


// ============================================================================
// The fixed yield definitions (PPG12.h)
// ============================================================================

double photonjet::ppg12::truthVertexWeight(const TH1* weights, const double z)
{
  if (!weights) return 0.0;
  if (!std::isfinite(z) || z <= -9000.0) return 1.0;
  if (z < weights->GetXaxis()->GetXmin() || z > weights->GetXaxis()->GetXmax()) return 0.0;
  const double w = const_cast<TH1*>(weights)->Interpolate(z);
  return std::isfinite(w) && w > 0.0 ? w : 0.0;
}


bool photonjet::ppg12::windowContains(const double value, const double low, const double high)
{
  if (!std::isfinite(value) || !std::isfinite(low) || value < low) return false;
  return std::isinf(high) || value <= high;
}

bool photonjet::ppg12::inResponseWindow(const double recoEt, const double truthEt)
{
  return std::isfinite(recoEt) && std::isfinite(truthEt) && recoEt > kRecoEtBins.front() &&
         recoEt < kRecoEtBins.back() && truthEt > kTruthEtBins.front() && truthEt < kTruthEtBins.back();
}

photonjet::ppg12::YieldShapes::YieldShapes(const PhotonRecord& photon)
    : probability(photon.clusterProbability),
      weta(photon.primaryShapes[0]),
      wphi(photon.primaryShapes[1]),
      et1(photon.primaryShapes[6]),
      et2(photon.primaryShapes[7]),
      et3(photon.primaryShapes[8]),
      et4(photon.primaryShapes[9]),
      e11OverE33(ratio(photon.primaryShapes[10], photon.primaryShapes[11])),
      e32OverE35(ratio(photon.primaryShapes[12], photon.primaryShapes[13])),
      npb(photon.npbCapture.rawScore),
      score(photon.productionScore)
{
}

bool photonjet::ppg12::YieldShapes::common() const
{
  const double widthRatio = std::isfinite(wphi) && std::isfinite(weta) && weta != 0.0
                                ? wphi / weta
                                : std::numeric_limits<double>::quiet_NaN();
  return inOpenInterval(probability, 0.0, 1.0) && inOpenInterval(e11OverE33, 0.0, kCommonE11E33Max) &&
         std::isfinite(widthRatio) && widthRatio > 0.0 && std::isfinite(weta) && weta < kCommonWetaMax &&
         std::isfinite(npb) && npb > kCommonNpbMin;
}

photonjet::ppg12::Tightness photonjet::ppg12::tightness(const YieldShapes& v, const double et)
{
  if (!v.common()) return Tightness::Fail;

  const bool prob = inOpenInterval(v.probability, 0.0, 1.0);
  const bool weta = inOpenInterval(v.weta, 0.0, 1.0);
  const bool wphi = inOpenInterval(v.wphi, 0.0, 1.0);
  const bool et1 = inOpenInterval(v.et1, kTightEt1Min, 1.0);
  const bool et2 = inOpenInterval(v.et2, 0.0, 1.0);
  const bool et3 = inOpenInterval(v.et3, 0.0, 1.0);
  const bool et4 = inOpenInterval(v.et4, 0.0, 1.0);
  const bool e11e33 = inOpenInterval(v.e11OverE33, 0.0, 1.0);
  const bool e32e35 = inOpenInterval(v.e32OverE35, kTightE32E35Min, 1.0);
  const double tightLine = PPG12IdLines::tight(et);
  const bool score =
      std::isfinite(v.score) && std::isfinite(tightLine) && v.score > tightLine && v.score < kTightScoreMax;
  if (prob && weta && wphi && et1 && et2 && et3 && et4 && e11e33 && e32e35 && score) return Tightness::Tight;

  // Non-tight: the loose shower windows, the non-tight score band, and at
  // least one failed tight criterion among w_eta, the probability and the score.
  const bool shape = prob && weta && wphi && inOpenInterval(v.et1, kNonTightEt1Min, 1.0) && et4 && e11e33 && e32e35;
  const double low = PPG12IdLines::nonTightLow(et), high = PPG12IdLines::nonTightHigh(et);
  const bool band =
      std::isfinite(v.score) && std::isfinite(low) && std::isfinite(high) && v.score > low && v.score < high;
  const int failures = int(!weta) + int(!prob) + int(!score);
  return shape && band && failures > 0 ? Tightness::NonTight : Tightness::Neither;
}

photonjet::ppg12::IsolationLookup photonjet::ppg12::isolationLookup(const PhotonRecord& photon,
                                                                    const double eventVertexZ)
{
  IsolationLookup lookup;
  const photonjet::IsolationRecord* cone = nullptr;
  for (const auto& candidate : photon.isolation)
    if (candidate.method == photonjet::IsolationMethod::Topocluster && candidate.radius == 0.4) cone = &candidate;
  if (cone && cone->valid && std::isfinite(cone->coneSum) && cone->coneSum < 1e8)
  {
    lookup.value = cone->coneSum;
    lookup.stages = {1, 2};
    if (std::isfinite(photon.producerVertexZ)) lookup.stages.push_back(6);
    return lookup;
  }
  // The recomputation reads the same topocluster cone and fails the same way.
  lookup.stages = {3, 4};
  if (cone && std::isfinite(cone->axisEta) && std::isfinite(cone->axisPhi) && std::isfinite(eventVertexZ))
    lookup.stages.push_back(6);
  return lookup;
}


int photonjet::production::PPG12WeightLaneCode(const WeightConfig& weights)
{
  if (!weights.ppg12Nominal) return 0;
  return 2 * (weights.ppg12LaneCode - 1) + (weights.doubleInteraction ? 2 : 1);
}


using namespace photonjet::ppg12;


// ============================================================================
// Event vertex/period audit
// ============================================================================

/**
 * Record one PPG12 vertex/period audit block for every active stream, before
 * (preVertexCut) or after the |vz| cut, exactly as the producer logged it.
 */
void PhotonJetTree::recordPPG12VertexAudit(const bool preVertexCut)
{
  const photonjet::EventClassificationConfig& c = m_config.classification;
  const bool periodContract = !c.period.key.empty();
  const std::string runPeriod =
      isData() && m_event.run > 0 ? periodLabel(c, m_event.run) : (periodContract ? c.period.key : "unset");
  const bool runPass = !c.filterDataByPeriod ||
                       (periodContract && m_event.run >= c.period.runMin && m_event.run < c.period.runMaxExclusive);

  double vertexWeight = 1.0;
  bool usedTruthWeight = false;
  const bool missingMinimumBiasTruth =
      isSimulation() && c.doubleInteraction &&
      (!std::isfinite(m_event.truthMinimumBiasVertexZ) || m_event.truthMinimumBiasVertexZ <= -9000.0);
  if (isSimulation() && m_config.weights.applyVertexReweight && m_vertexReweight)
  {
    vertexWeight = truthVertexWeight(m_vertexReweight, m_event.truthHardVertexZ);
    if (c.doubleInteraction) vertexWeight *= truthVertexWeight(m_vertexReweight, m_event.truthMinimumBiasVertexZ);
    usedTruthWeight = true;
  }
  // Without a period contract the producer's unset values are recorded.
  const double luminosityWeight =
      c.useLuminosityWeight && periodContract ? c.period.luminosityPbInv / c.period.targetLuminosityPbInv : 1.0;
  const double singleFraction = periodContract ? 1.0 - c.period.doubleFraction : -1.0;
  const double periodEventWeight = isSimulation() ? vertexWeight * c.mixWeight * luminosityWeight : 1.0;
  // The classification precedes the weight configuration: every row carries
  // the record's unit fill weight.

  for (const std::string& stream : m_event.activeTriggers)
  {
    auto audit = [&](const std::int32_t stage)
    {
      photonjet::PPG12DiagnosticRecord row;
      row.trigger = stream;
      row.kind = "VERTEX_CONTRACT_AUDIT";
      row.stageCode = stage;
      m_event.ppg12Diagnostics.push_back(row);
    };
    audit(1);
    if (periodContract) audit(2);
    if (c.filterDataByPeriod) audit(3);
    audit(runPeriod == "0mrad" ? 4 : runPeriod == "1p5mrad" ? 5 : 6);
    if (isSimulation()) audit(c.doubleInteraction ? 8 : 7);
    if (c.mixWeightFromPeriod)
      audit(9);
    else if (isSimulation())
      audit(10);
    if (c.vertexFileFromPeriod)
      audit(11);
    else if (isSimulation())
      audit(12);
    if (isSimulation() && c.useLuminosityWeight) audit(13);
    if (usedTruthWeight) audit(14);
    if (missingMinimumBiasTruth) audit(15);
    if (c.filterDataByPeriod) audit(runPass ? 16 : 17);
    audit(preVertexCut ? 18 : 19);
    if (!periodContract) audit(20);

    const std::array<double, 13> values{
        {static_cast<double>(c.period.runMin), static_cast<double>(c.period.runMaxExclusive), c.period.luminosityPbInv,
         c.period.targetLuminosityPbInv, luminosityWeight, singleFraction, c.period.doubleFraction,
         c.period.closureZCutCm, c.doubleInteraction ? 1.0 : 0.0, c.mixWeight, c.mixWeightFromPeriod ? 1.0 : 0.0,
         c.vertexFileFromPeriod ? 1.0 : 0.0, c.filterDataByPeriod ? 1.0 : 0.0}};
    for (std::size_t index = 0; index < values.size(); ++index)
    {
      photonjet::PPG12DiagnosticRecord row;
      row.trigger = stream;
      row.kind = "PERIOD_CONTRACT_VALUE";
      row.stageCode = static_cast<std::int32_t>(index) + 1;
      row.flowValue = values[index];
      m_event.ppg12Diagnostics.push_back(row);
    }

    photonjet::PPG12DiagnosticRecord row;
    row.trigger = stream;
    row.kind = "VERTEX_QA_STAGE";
    row.stageCode = preVertexCut ? 18 : 19;
    row.recoVertexZ = m_event.recoVertexZ;
    row.hardTruthVertexZ = m_event.truthHardVertexZ;
    row.minimumBiasTruthVertexZ = m_event.truthMinimumBiasVertexZ;
    row.vertexWeight = vertexWeight;
    row.periodEventWeight = periodEventWeight;
    m_event.ppg12Diagnostics.push_back(row);
  }
}


// ============================================================================
// pp simulation sample ownership
//
// A pp simulation sample owns the events of its generator slice. Inclusive-jet
// samples decide before the event cuts, from the leading truth R=0.4 jet (an
// event outside the window ends with terminal status 4). Photon+jet samples log
// their stitch source before the event cuts and, for an accepted event, every
// established stitch definition; the stored Geant4 photon decides.
// ============================================================================

/**
 * Before the event cuts (beforeEventCuts) or after them: the ownership rows and
 * decision. Returns false only for an inclusive-jet event outside its window.
 */
bool PhotonJetTree::recordPPG12SampleOwnership(const bool beforeEventCuts)
{
  const photonjet::WeightConfig& w = m_config.weights;
  if (!w.ppg12Nominal) return true;
  const double low = w.ownershipLowGeV, high = w.ownershipHighGeV;
  auto source = [&](const std::int32_t code, const double pt, const bool pass)
  {
    photonjet::PPG12DiagnosticRecord row;
    row.trigger = "SIM";
    row.kind = "PP_STITCH_SOURCE";
    row.sourceCode = code;
    row.sampleCode = w.ppg12SliceCode;
    row.decision = pass ? 1 : pt >= 0.0 ? 2 : 3;
    if (pt >= 0.0) row.sourcePhotonPt = pt;
    m_event.ppg12Diagnostics.push_back(row);
  };

  if (w.ppg12LaneCode == 2)
  {
    if (!beforeEventCuts) return true;
    m_event.ppg12Owned =
        m_event.ppg12MaxTruthJetPtR04 >= 0.0 && windowContains(m_event.ppg12MaxTruthJetPtR04, low, high);
    return m_event.ppg12Owned;
  }

  if (beforeEventCuts)
  {
    source(7, m_event.ppg12SourcePhotonPt, windowContains(m_event.ppg12SourcePhotonPt, low, high));
    return true;
  }
  // The stitch definitions keep low <= pT < high; the |eta| < 0.7 source
  // includes the upper edge.
  for (std::size_t variant = 0; variant < m_event.ppg12StitchPhotonPt.size(); ++variant)
  {
    const double pt = m_event.ppg12StitchPhotonPt[variant];
    source(static_cast<std::int32_t>(variant) + 1, pt, pt >= 0.0 && pt >= low && pt < high);
  }
  source(6, m_event.ppg12SourcePhotonPt, windowContains(m_event.ppg12SourcePhotonPt, low, high));
  const double canonical = m_event.ppg12StitchPhotonPt[0];
  m_event.ppg12Owned = canonical >= 0.0 && canonical >= low && canonical < high;
  return true;
}


// ============================================================================
// Photon-yield audit
//
// Every photon of an accepted event inside the analysis ET bins entered the
// yield histograms, in container order. For each one the producer logged:
//
//   EISO_VERTEX_AUDIT        the path of each isolation lookup: one for the
//                            isolation QA, one for the yield and, for a photon
//                            passing the active preselection, one per stream
//                            (DATA rows carry the scaled-bit audit stream)
//   PHOTON_ALL_COMMON        per stream, its ET and common-selection bits
//   PHOTON_TIGHT_ABCD        per stream, a tight or non-tight photon with its
//                            A/B/C/D region bits
//   PP_ISOLATION_STACK_FILL  per stream, a tight or non-tight photon passing
//                            the common selection, in the reconstructed ET bins
// ============================================================================

void PhotonJetTree::recordPPG12PhotonYield()
{
  const photonjet::EventClassificationConfig& c = m_config.classification;
  if (!c.enabled || !isPP() || !m_event.nodesReady || m_event.terminalStatus != 0 || m_event.activeTriggers.empty() ||
      !isData())
    return;

  const std::vector<double>& bins = m_config.ppg12.photonEtBins;
  const std::string& auditStream = c.scaledBitStream;
  const auto record = [&](photonjet::PPG12DiagnosticRecord row) { m_event.ppg12Diagnostics.push_back(std::move(row)); };
  const auto isolationAudit = [&](const IsolationLookup& lookup)
  {
    for (const std::int32_t stage : lookup.stages)
    {
      photonjet::PPG12DiagnosticRecord row;
      row.trigger = auditStream;
      row.kind = "EISO_VERTEX_AUDIT";
      row.stageCode = stage;
      record(row);
    }
  };

  for (const PhotonRecord& photon : m_photons)
  {
    const double et = photon.et;
    if (!std::isfinite(et) || !std::isfinite(photon.eta) || !std::isfinite(photon.phi)) continue;
    if (photon.ppg12TowerMaskState == 1) continue;
    if (et < 5.0 || !(std::fabs(photon.eta) < m_config.photon.maxAbsEta)) continue;
    bool inBin = false;
    for (std::size_t i = 0; i + 1 < bins.size(); ++i) inBin |= et >= bins[i] && et < bins[i + 1];
    if (!inBin) continue;

    const IsolationLookup isolation = isolationLookup(photon, m_event.recoVertexZ);
    isolationAudit(isolation);  // isolation QA
    isolationAudit(isolation);  // yield

    const YieldShapes shapes(photon);
    const bool common = shapes.common();
    const double weight = dataTriggerEfficiencyWeight(et);
    for (const std::string& stream : m_event.activeTriggers)
    {
      photonjet::PPG12DiagnosticRecord row;
      row.trigger = stream;
      row.kind = "PHOTON_ALL_COMMON";
      row.recoPhotonPt = et;
      row.occurrenceWeight = weight;
      // common | signal response window (no truth class in DATA) | inclusive
      row.selectionBits = (common ? 1ULL : 0ULL) | 8ULL | 16ULL;
      record(row);
    }

    // The isolation regions use the stored R = 0.4 witness lines.
    const photonjet::IsolationWitnessRecord* witness = nullptr;
    for (const auto& candidate : photon.isolationWitnesses)
      if (candidate.radius == 0.4) witness = &candidate;
    const double eiso = isolation.value;
    const bool valid = std::isfinite(eiso) && eiso < 1e8 && witness;
    const bool isolated = valid && eiso > kIsolatedMin && eiso < witness->threshold;
    const bool nonIsolated = valid && eiso > witness->sidebandThreshold && eiso < kNonIsolatedMax;

    const Tightness tag = tightness(shapes, et);
    const bool tight = tag == Tightness::Tight, nonTight = tag == Tightness::NonTight;
    if (tight || nonTight)
    {
      const double truthEt = std::numeric_limits<double>::quiet_NaN();
      for (const std::string& stream : m_event.activeTriggers)
      {
        photonjet::PPG12DiagnosticRecord row;
        row.trigger = stream;
        row.kind = "PHOTON_TIGHT_ABCD";
        row.recoPhotonPt = et;
        row.responsePhotonPt = et;
        row.truthPhotonPt = truthEt;
        row.occurrenceWeight = weight;
        row.selectionBits = (tight ? 1ULL : 0ULL) | (nonTight ? 2ULL : 0ULL) | (isolated ? 4ULL : 0ULL) |
                            (nonIsolated ? 8ULL : 0ULL) | (inResponseWindow(et, truthEt) ? 64ULL | 128ULL : 0ULL) |
                            256ULL;
        record(row);
      }
    }

    if (common && (tight || nonTight) && std::isfinite(eiso) && eiso <= 1e8)
    {
      std::int32_t stackBin = -1;
      for (std::size_t i = 0; i + 1 < kRecoEtBins.size() && stackBin < 0; ++i)
        if (et > kRecoEtBins[i] && et < kRecoEtBins[i + 1]) stackBin = static_cast<std::int32_t>(i);
      if (stackBin >= 0)
        for (const std::string& stream : m_event.activeTriggers)
        {
          photonjet::PPG12DiagnosticRecord row;
          row.trigger = stream;
          row.kind = "PP_ISOLATION_STACK_FILL";
          row.sampleCode = stackBin + 1;
          row.decision = tight ? 1 : 2;
          row.flowValue = eiso;
          row.occurrenceWeight = weight;
          record(row);
        }
    }

    // A photon passing the active preselection is classified per stream.
    if (photonjet::production::ActivePreselection(m_config.photon, photon) == 1)
      for (std::size_t i = 0; i < m_event.activeTriggers.size(); ++i) isolationAudit(isolation);
  }
}
