/**
 * @file PPG12Simulation.cc
 * The established pp simulation producer's histogram bookkeeping.
 *
 * Once a simulated event's PPG12 weight was configured, the established
 * producer filled its analysis histograms: jet QA, photon identification and
 * isolation spectra, ABCD purity counters, truth efficiencies, response and
 * unfolding inputs. Every fill went through a weighted or a raw counting
 * histogram, and each event recorded how many fills of each kind it made per
 * weighting lane (ppg12_wfill*, ppg12_rfill*), beside an audit occurrence for
 * the fills it checked (pdiag_*).
 *
 * This file walks the producer's decisions over the captured event, in its
 * order, counting those fills and logging those occurrences. It fills no
 * histogram and changes no event, object or label. Section titles name the
 * producer's histogram families.
 */
#include "../PhotonJetTree.h"
#include "PPG12.h"
#include "Production.h"

#include <TH1.h>
#include <TRandom3.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using namespace photonjet;
using ppg12::Tightness;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------------------
// The producer's fixed simulation definitions
// ---------------------------------------------------------------------------

constexpr double kPhotonEtFloor = 5.0;          // reconstructed photon floor, GeV
constexpr double kOverlapDeltaR2 = 0.4 * 0.4;   // photon-jet overlap veto
constexpr double kAcceptanceEtaMax = 0.7;       // truth-flow and Fig. 8/24 acceptance
constexpr double kTruthFigureVertexMax = 30.0;  // Fig. 2/3/8 |z|, cm
constexpr double kFig8RecoEtMin = 5.0;
constexpr std::array<double, 14> kFig24EtEdges{{10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36}};

// Simulated isolation is scaled to data: 1.2 * raw + 0.1 GeV.
constexpr double kSimIsolationScale = 1.2;
constexpr double kSimIsolationShift = 0.1;

// The response arm adds the data-minus-simulation energy resolution.
constexpr double kResolutionData[3] = {0.15, 0.05, 0.05};
constexpr double kResolutionSim[3] = {0.185, 0.0, 0.04};

// The photon-ID shower windows of the "reference" tight variant: widths below
// 0.15 + 0.006 ET, 0.4 < E11/E33 < 0.98, 0.9 < et1 < 1, 0.92 < E32/E35 < 1.
constexpr double kTightWidthIntercept = 0.15, kTightWidthSlope = 0.006;

bool inOpen(const double x, const double low, const double high)
{
  return std::isfinite(x) && std::isfinite(low) && std::isfinite(high) && x > low && x < high;
}

// The pseudorapidity reach of a jet collection, |eta| < 1.1 - R, with R read
// from its key digits ("r04" -> 0.1 * 4) as the producer did.
double jetEtaMax(const int radiusDigits)
{
  const double etaMax = 1.1 - 0.1 * radiusDigits;
  return etaMax > 0.0 ? etaMax : 0.0;
}

double deltaR(const double eta1, const double phi1, const double eta2, const double phi2)
{
  const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

// The shower variables of the producer's reference ID decisions.
struct Shapes
{
  double et, weta, wphi, et1, e11OverE33, e32OverE35;

  Shapes(const PhotonRecord& photon, const double candidateEt)
  {
    const ppg12::YieldShapes v(photon);
    et = candidateEt;
    weta = v.weta;
    wphi = v.wphi;
    et1 = v.et1;
    e11OverE33 = v.e11OverE33;
    e32OverE35 = v.e32OverE35;
  }

  // Reference preselection: finite inputs, E11/E33 < 0.98, 0.6 < et1 < 1,
  // 0.8 < E32/E35 < 1 and w_eta < 0.6.
  bool preselected() const
  {
    return std::isfinite(weta) && std::isfinite(wphi) && std::isfinite(et1) && std::isfinite(e11OverE33) &&
           std::isfinite(e32OverE35) && std::isfinite(et) && e11OverE33 < 0.98 && inOpen(et1, 0.6, 1.0) &&
           inOpen(e32OverE35, 0.8, 1.0) && weta < 0.6;
  }

  // The five reference tight windows; bit i set when window i fails, in the
  // order w_eta, w_phi, E11/E33, et1, E32/E35.
  int tightFailures() const
  {
    const double widthMax = std::isfinite(et) && et > 0.0 ? kTightWidthIntercept + kTightWidthSlope * et : kNaN;
    int mask = 0;
    if (!inOpen(weta, 0.0, widthMax)) mask |= 1;
    if (!inOpen(wphi, 0.0, widthMax)) mask |= 2;
    if (!inOpen(e11OverE33, 0.40, 0.98)) mask |= 4;
    if (!inOpen(et1, 0.90, 1.00)) mask |= 8;
    if (!inOpen(e32OverE35, 0.92, 1.00)) mask |= 16;
    return mask;
  }
};

// Tight with no failed window, non-tight with two or more.
Tightness referenceTag(const int failures)
{
  const int count = std::popcount(static_cast<unsigned>(failures));
  return count == 0 ? Tightness::Tight : count >= 2 ? Tightness::NonTight : Tightness::Neither;
}

// The producer's classifyPhotonTightness: preselection first; an undefined
// width window is neither.
Tightness classifiedTag(const Shapes& v)
{
  if (!v.preselected()) return Tightness::Fail;
  if (!(std::isfinite(v.et) && v.et > 0.0)) return Tightness::Neither;
  return referenceTag(v.tightFailures());
}

// ---------------------------------------------------------------------------
// One event's bookkeeping
// ---------------------------------------------------------------------------

struct SimulationEvent
{
  const Config& config;
  const EventRecord& event;
  const std::vector<PhotonRecord>& photons;                 // every candidate, container order
  const std::vector<PPG12ClusterRecord>& clusters;          // the builder's input clusters
  const std::vector<PPG12ClusterRecord>* responseClusters;  // Fig. 8 node; nullptr: the input clusters
  const std::vector<JetRecord>& jets;
  const std::vector<TruthPhotonRecord>& truthPhotons;
  const std::vector<TruthJetRecord>& truthJets;
  const std::vector<int>& truthJetRadiiFound;  // radius codes of the truth containers read
  bool evaluatorReady;                         // the cluster truth evaluator had its nodes
};

class Bookkeeping
{
 public:
  Bookkeeping(const SimulationEvent& in, std::vector<PPG12DiagnosticRecord>& log);

  void run();

  std::int64_t weighted = 0;
  std::int64_t raw = 0;

 private:
  // A reconstructed jet or truth jet as the histogram code reads it.
  struct Jet
  {
    double pt, eta, phi;
  };
  // The jets of one radius key ("r02", "r03", "r04").
  struct Radius
  {
    int digits;
    std::vector<Jet> reco;
    bool truthPresent = false;
    std::vector<Jet> truth;
  };
  // A truth photon of the producer's truth maps (keyed by Geant4 track).
  struct TruthPhoton
  {
    const TruthPhotonRecord* row;
    bool signal;
  };

  void fill(const std::int64_t n = 1) { weighted += n * streams(); }
  void fillRaw(const std::int64_t n = 1) { raw += n * streams(); }
  std::int64_t streams() const { return static_cast<std::int64_t>(m_in.event.activeTriggers.size()); }
  PPG12DiagnosticRecord occurrence(const std::string& stream, const char* kind) const;
  void log(PPG12DiagnosticRecord row) { m_log.push_back(std::move(row)); }

  // Producer helpers
  double isoThreshold(const double et) const { return m_isoIntercept + m_isoSlope * et; }
  double nonIsoThreshold(const double et) const { return isoThreshold(et) + m_isoGap; }
  int photonBin(double et) const;
  double simIsolation(const PhotonRecord& photon);
  double rawIsolation(const PhotonRecord& photon);
  bool isolated(const PhotonRecord& photon, double et);
  const TruthPhoton* signalByTrack(int trackId) const;
  bool primaryIsSignal(const DominantTruthRecord& truth) const;
  bool primaryIsSignal(const PPG12ClusterRecord& cluster) const;
  const PhotonRecord* matchedReco(int trackId, double& et, double& eta, double& phi) const;
  double responseEt(double recoEt, double truthPt, int candidateIndex) const;

  // Producer stages, in its order
  void jetQuality(const Radius& radius);
  void candidates();
  void truthLeakage();
  void responseTargets();
  void truthRecoil();
  void inclusiveJetSample();
  void photonLoop();
  void yieldAllCommon(double et, bool common, bool haveTruthClass, bool truthSignal, bool window, bool inclusive);
  void yieldTightAndAbcd(double et, double responseEt, double truthPt, Tightness tag, bool iso, bool nonIso,
                         bool haveTruthClass, bool truthSignal, bool inclusive);
  void isolationStack(double et, double isolation, Tightness tag);
  void abcdCounters(const PhotonRecord& photon, const Shapes& v, Tightness tag, double et);
  void leadPhotonJets(int leadBin, double leadEt, double leadEta, double leadPhi, bool haveTruthSignal,
                      double truthSignalPt, bool haveTruthPhoton, double tPt, double tEta, double tPhi);
  void unfoldingResponse(const Radius& radius, double leadEt, double leadEta, double leadPhi, bool haveTruthPhoton,
                         double tPt, double tEta, double tPhi, const std::vector<const Jet*>& fid,
                         const std::vector<bool>& fidRecoil, const Jet* recoil1);
  void truthJes(const Radius& radius, double tPt, double tEta, double tPhi, const Jet* recoil1);
  void sidebandRecoil(double leadEt, double leadEta, double leadPhi);
  void leadTag(const char* name, double value);
  void efficiencyStage(const std::string& name, double truthPt);
  void stitchFlow(int stage, double value);

  const SimulationEvent& m_in;
  std::vector<PPG12DiagnosticRecord>& m_log;
  const PPG12BookkeepingConfig& m_book;
  const WeightConfig& m_weights;

  double m_vertexZ, m_weight;
  double m_isoIntercept, m_isoSlope, m_isoGap;
  bool m_photonSample, m_inclusiveSample;
  double m_windowLow, m_windowHigh;
  std::array<bool, 5> m_stitchKept{};
  std::vector<TruthPhoton> m_truth;  // the witness map, ascending track
  std::vector<Radius> m_radii;
};

Bookkeeping::Bookkeeping(const SimulationEvent& in, std::vector<PPG12DiagnosticRecord>& log)
    : m_in(in),
      m_log(log),
      m_book(in.config.ppg12),
      m_weights(in.config.weights),
      m_vertexZ(in.event.recoVertexZ),
      m_weight(in.event.eventWeight)
{
  // The pp sliding working point of the R = 0.4 cone.
  const auto& radii = in.config.photon.isolationRadii;
  const auto r04 = std::find_if(radii.begin(), radii.end(), [](const double r) { return std::fabs(r - 0.4) < 1e-6; });
  if (r04 == radii.end() || in.config.photon.isolationWitness.lines.size() != radii.size())
    throw std::runtime_error("PPG12 simulation bookkeeping needs the R = 0.4 isolation working point");
  const IsolationLine& line = in.config.photon.isolationWitness.lines[static_cast<std::size_t>(r04 - radii.begin())];
  m_isoIntercept = line.interceptGeV;
  m_isoSlope = line.slope;
  m_isoGap = line.sidebandGapGeV;

  m_photonSample = m_weights.ppg12LaneCode == 1;
  m_inclusiveSample = m_weights.ppg12LaneCode == 2;
  m_windowLow = m_weights.ownershipLowGeV;
  m_windowHigh = m_weights.ownershipHighGeV;
  // Photon samples keep a stitch-flow diagnostic for every stitch definition
  // owning the event (low <= pT < high).
  for (std::size_t v = 0; v < m_stitchKept.size(); ++v)
  {
    const double pt = in.event.ppg12StitchPhotonPt[v];
    m_stitchKept[v] = m_photonSample && pt >= 0.0 && pt >= m_windowLow && pt < m_windowHigh;
  }

  // The producer's truth maps: every Geant4 truth photon in track order; the
  // signal map is the analysis-signal subset.
  for (const TruthPhotonRecord& row : in.truthPhotons) m_truth.push_back({&row, row.analysisSignal});

  // Jet collections by radius key, in the producer's key order.
  for (const int digits : {2, 3, 4})
  {
    Radius radius;
    radius.digits = digits;
    const double r = 0.1 * digits;
    std::vector<const JetRecord*> reco;
    for (const JetRecord& jet : in.jets)
      if (jet.view == JetView::PP && std::fabs(jet.radius - r) < 1e-6) reco.push_back(&jet);
    std::sort(reco.begin(), reco.end(),
              [](const JetRecord* a, const JetRecord* b) { return a->encounterOrdinal < b->encounterOrdinal; });
    for (const JetRecord* jet : reco) radius.reco.push_back({jet->unitScalePt, jet->eta, jet->phi});
    radius.truthPresent = std::find(in.truthJetRadiiFound.begin(), in.truthJetRadiiFound.end(), 10 * digits) !=
                          in.truthJetRadiiFound.end();
    for (const TruthJetRecord& jet : in.truthJets)
      if (std::fabs(jet.radius - r) < 1e-6) radius.truth.push_back({jet.pt, jet.eta, jet.phi});
    m_radii.push_back(std::move(radius));
  }
}

PPG12DiagnosticRecord Bookkeeping::occurrence(const std::string& stream, const char* kind) const
{
  PPG12DiagnosticRecord row;
  row.trigger = stream;
  row.kind = kind;
  return row;
}

int Bookkeeping::photonBin(const double et) const
{
  const auto& bins = m_book.photonEtBins;
  if (bins.size() < 2 || !std::isfinite(et)) return -1;
  for (std::size_t i = 0; i + 1 < bins.size(); ++i)
    if (et >= bins[i] && et < bins[i + 1]) return static_cast<int>(i);
  return -1;
}

// One producer lookup of the topocluster isolation: logs its audit stages
// (raw fills of the audit histogram) and returns the raw value.
double Bookkeeping::rawIsolation(const PhotonRecord& photon)
{
  const ppg12::IsolationLookup lookup = ppg12::isolationLookup(photon, m_vertexZ);
  for (const std::int32_t stage : lookup.stages)
  {
    PPG12DiagnosticRecord row = occurrence("SIM", "EISO_VERTEX_AUDIT");
    row.stageCode = stage;
    log(row);
    raw += 1;  // the audit histogram lives in the SIM directory, once per lookup
  }
  return lookup.value;
}

double Bookkeeping::simIsolation(const PhotonRecord& photon)
{
  const double value = rawIsolation(photon);
  return std::isfinite(value) && value < 1e8 ? kSimIsolationScale * value + kSimIsolationShift : value;
}

bool Bookkeeping::isolated(const PhotonRecord& photon, const double et)
{
  if (!std::isfinite(et) || et <= 0.0) return false;
  return simIsolation(photon) < isoThreshold(et);
}

const Bookkeeping::TruthPhoton* Bookkeeping::signalByTrack(const int trackId) const
{
  for (const TruthPhoton& truth : m_truth)
    if (truth.signal && truth.row->trackId == trackId) return &truth;
  return nullptr;
}

// The producer's truth tag of a cluster: its dominant primary is a signal
// photon. An empty signal map tags nothing.
bool Bookkeeping::primaryIsSignal(const DominantTruthRecord& truth) const
{
  const bool primary =
      truth.state == DominantTruthState::ValidPrimary || truth.state == DominantTruthState::InvalidPrimary;
  return primary && signalByTrack(truth.trackId);
}

bool Bookkeeping::primaryIsSignal(const PPG12ClusterRecord& cluster) const
{
  return cluster.primaryFound && signalByTrack(cluster.truthTrackId);
}

// The reconstructed photon matched to a signal truth photon: one whose
// dominant primary is that photon, with |eta| < 0.7 and ET above 5 GeV; the
// largest energy contribution wins, then the smallest Delta R.
const PhotonRecord* Bookkeeping::matchedReco(const int trackId, double& et, double& eta, double& phi) const
{
  const TruthPhoton* target = signalByTrack(trackId);
  if (!target) return nullptr;
  const PhotonRecord* best = nullptr;
  float bestContribution = std::numeric_limits<float>::lowest();
  double bestDeltaR = std::numeric_limits<double>::infinity();
  for (const PhotonRecord& photon : m_in.photons)
  {
    const double photonEta = photon.eta, photonPhi = TVector2::Phi_mpi_pi(photon.phi), photonEt = photon.et;
    if (!std::isfinite(photonEta) || !std::isfinite(photonPhi) || !std::isfinite(photonEt)) continue;
    if (std::fabs(photonEta) >= kAcceptanceEtaMax || photonEt <= kPhotonEtFloor) continue;
    if (!primaryIsSignal(photon.dominantTruth) || photon.dominantTruth.trackId != trackId) continue;
    const double dphi = TVector2::Phi_mpi_pi(target->row->phi - photonPhi);
    const double deta = target->row->eta - photonEta;
    const double dr = std::sqrt(deta * deta + dphi * dphi);
    const float contribution = static_cast<float>(photon.dominantTruth.energyContribution);
    if (!best || contribution > bestContribution || (contribution == bestContribution && dr < bestDeltaR))
    {
      best = &photon;
      et = photonEt;
      eta = photonEta;
      phi = photonPhi;
      bestDeltaR = dr;
      bestContribution = contribution;
    }
  }
  return best;
}

// The response ET of a signal cluster: one Gaussian draw of the extra
// data-minus-simulation resolution, seeded per physical event and candidate.
double Bookkeeping::responseEt(const double recoEt, const double truthPt, const int candidateIndex) const
{
  if (!std::isfinite(recoEt) || recoEt <= 0.0 || !std::isfinite(truthPt) || truthPt <= 0.0) return recoEt;
  const auto sigma2 = [truthPt](const double* p)
  { return p[0] * p[0] / truthPt + p[1] * p[1] / (truthPt * truthPt) + p[2] * p[2]; };
  const double sigma = std::sqrt(std::max(0.0, sigma2(kResolutionData) - sigma2(kResolutionSim))) * truthPt;
  if (!(sigma > 0.0) || !std::isfinite(sigma)) return recoEt;

  // The established jobs read each simulation file from its first entry, so
  // their event counter is the source entry plus one.
  const auto eventCount = static_cast<std::uint64_t>(m_in.event.sourceEntry + 1);
  std::uint64_t key = m_book.responseSeedRun;
  const auto mix = [&key](const std::uint64_t value)
  { key ^= value + 0x9e3779b97f4a7c15ULL + (key << 6U) + (key >> 2U); };
  mix(m_in.event.physicalEventSequenceValid
          ? static_cast<std::uint64_t>(static_cast<std::uint32_t>(m_in.event.physicalEventSequence))
          : eventCount);
  mix(eventCount);
  mix(static_cast<std::uint64_t>(std::max(candidateIndex, 0)));
  mix(std::bit_cast<std::uint64_t>(recoEt));
  mix(std::bit_cast<std::uint64_t>(truthPt));
  key += 0x9e3779b97f4a7c15ULL;
  key = (key ^ (key >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  key = (key ^ (key >> 27U)) * 0x94d049bb133111ebULL;
  key ^= key >> 31U;
  unsigned int seed = static_cast<unsigned int>((key ^ (key >> 32U)) & 0xffffffffULL);
  if (seed == 0U) seed = 42U;
  TRandom3 random(seed);
  const double smeared = recoEt + random.Gaus(0.0, sigma);
  return std::isfinite(smeared) && smeared > 0.0 ? smeared : recoEt;
}

// ===========================================================================
// The event after weight configuration
// ===========================================================================

void Bookkeeping::run()
{
  // The weight configuration logs its own audit stage.
  PPG12DiagnosticRecord configured = occurrence("SIM", "EISO_VERTEX_AUDIT");
  configured.stageCode = 7;
  log(configured);
  raw += 1;

  // Events outside the vertex window stop here.
  const auto& classification = m_in.config.classification;
  if (!std::isfinite(m_vertexZ) || (classification.useVertexCut && std::fabs(m_vertexZ) >= classification.vertexCutCm))
    return;

  for (const Radius& radius : m_radii) jetQuality(radius);
  candidates();
}

// ===========================================================================
// Inclusive jet QA (h_jet*_all, h_jet*_incl, h_nJets, h_HT, h_leadJet*)
//
// The calibrated jets carry no energy, so their mass is undefined and the
// mass histograms never fill.
// ===========================================================================

void Bookkeeping::jetQuality(const Radius& radius)
{
  const double etaMax = jetEtaMax(radius.digits);
  int fiducial = 0;
  double leadPt = -1.0, subleadPt = -1.0;
  for (const Jet& jet : radius.reco)
  {
    if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta)) continue;
    if (jet.pt < m_book.recoilJetPtMinGeV) continue;
    fill(4);  // pT, eta, phi, eta-phi of all jets
    if (std::fabs(jet.eta) >= etaMax) continue;
    ++fiducial;
    if (jet.pt > leadPt)
    {
      subleadPt = leadPt;
      leadPt = jet.pt;
    }
    else if (jet.pt > subleadPt)
      subleadPt = jet.pt;
    fill(4);  // the same for fiducial jets
  }
  (void)fiducial;
  fill(2);                     // multiplicity, HT
  if (leadPt >= 0.0) fill(3);  // leading pT, eta, phi
  if (subleadPt >= 0.0) fill(1);
}

// ===========================================================================
// Candidate processing
// ===========================================================================

void Bookkeeping::candidates()
{
  truthLeakage();
  truthRecoil();
  if (m_inclusiveSample) inclusiveJetSample();
  photonLoop();
}

// ---------------------------------------------------------------------------
// Ownership audits, truth isolation, efficiency stages and stitch flow
// ---------------------------------------------------------------------------

void Bookkeeping::stitchFlow(const int stage, const double value)
{
  if (!m_photonSample || !std::isfinite(value)) return;
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    PPG12DiagnosticRecord row = occurrence(stream, "PP_STITCH_FLOW");
    row.stageCode = stage;
    row.flowValue = value;
    log(row);
    raw += 1;
  }
}

void Bookkeeping::efficiencyStage(const std::string& name, const double truthPt)
{
  if (!std::isfinite(truthPt) || truthPt <= 0.0) return;
  static const std::map<std::string, int> kFig6Stage = {{"h_photonEffPpg12Fig6TruthDen_pTgamma_0", 1},
                                                        {"h_photonEffPpg12Fig6Reco_pTgamma_0", 2},
                                                        {"h_photonEffPpg12Fig6RecoTight_pTgamma_0", 3},
                                                        {"h_photonEffPpg12Fig6RecoIso_pTgamma_0", 4},
                                                        {"h_photonEffPpg12Fig6RecoTightIso_pTgamma_0", 5}};
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    weighted += 1;
    if (name == "h_photonEffTruthDen_pTgamma_0")
    {
      PPG12DiagnosticRecord row = occurrence(stream, "TRUTH_EFFICIENCY_DEN_FILL");
      row.truthPhotonPt = truthPt;
      row.occurrenceWeight = m_weight;
      log(row);
    }
    if (const auto stage = kFig6Stage.find(name); stage != kFig6Stage.end())
    {
      PPG12DiagnosticRecord row = occurrence(stream, "PP_FIG6_EFFICIENCY_STAGE");
      row.stageCode = stage->second;
      row.truthPhotonPt = truthPt;
      row.occurrenceWeight = m_weight;
      log(row);
    }
  }
}

void Bookkeeping::truthLeakage()
{
  if (!m_in.evaluatorReady) return;

  // Sample ownership for the truth figures: Fig. 24 by the leading stored
  // photon inside |eta| < 0.7, Figs. 3 and 6 by the leading stored photon.
  const bool haveWindow = m_photonSample && std::isfinite(m_windowLow);
  const double fig24Owner = m_photonSample ? m_in.event.ppg12SourcePhotonPt : -1.0;
  const double fig6Owner = m_photonSample ? m_in.event.ppg12StitchPhotonPt[0] : -1.0;
  const bool fig24HaveOwner = fig24Owner >= 0.0;
  const bool fig24Owned = haveWindow && fig24HaveOwner && ppg12::windowContains(fig24Owner, m_windowLow, m_windowHigh);
  const bool fig6HaveOwner = fig6Owner >= 0.0;
  const bool fig6Accepted =
      m_photonSample && haveWindow && fig6HaveOwner && ppg12::windowContains(fig6Owner, m_windowLow, m_windowHigh);

  const auto ownershipAudit = [&](const char* kind, const int stage)
  {
    for (const std::string& stream : m_in.event.activeTriggers)
    {
      weighted += 1;
      PPG12DiagnosticRecord row = occurrence(stream, kind);
      row.stageCode = stage;
      row.occurrenceWeight = m_weight;
      log(row);
    }
  };
  ownershipAudit("FIG3_OWNERSHIP_AUDIT", 1);
  if (m_photonSample) ownershipAudit("FIG3_OWNERSHIP_AUDIT", 2);
  if (fig6HaveOwner) ownershipAudit("FIG3_OWNERSHIP_AUDIT", 3);
  ownershipAudit("FIG3_OWNERSHIP_AUDIT", fig6Accepted ? 4 : 5);

  for (const TruthPhoton& truth : m_truth)
  {
    const TruthPhotonRecord& t = *truth.row;
    const int photonClass = static_cast<int>(t.promptClass);
    const double truthIsolation = t.isolationR03;
    const bool passIsolation = truth.signal;
    const bool directOrFragment = photonClass == 1 || photonClass == 2;
    const bool closeVertex = std::isfinite(m_vertexZ) && std::fabs(m_vertexZ) < kTruthFigureVertexMax;

    // Fig. 3: truth isolation (and, for isolated photons in 10-35 GeV, pT)
    // of direct and fragmentation photons in an accepted sample.
    if (directOrFragment && std::isfinite(t.pt) && t.pt > 0.0 && closeVertex && fig6Accepted)
    {
      if (std::isfinite(truthIsolation))
        for (const std::string& stream : m_in.event.activeTriggers)
        {
          weighted += 1;
          PPG12DiagnosticRecord row = occurrence(stream, "PP_FIG3_TRUTH_ISOLATION");
          row.sourceCode = photonClass;
          row.truthPhotonPt = t.pt;
          row.flowValue = truthIsolation;
          row.occurrenceWeight = m_weight;
          log(row);
        }
      if (passIsolation && t.pt >= 10.0 && t.pt < 35.0) fill(2);
    }

    // Truth isolation QA, and Fig. 2 in three pT windows.
    if (std::isfinite(truthIsolation))
    {
      const bool fig2Window =
          (t.pt >= 10.0 && t.pt < 15.0) || (t.pt >= 15.0 && t.pt < 20.0) || (t.pt >= 25.0 && t.pt < 30.0);
      if (directOrFragment && std::isfinite(t.pt) && t.pt > 0.0 && closeVertex && fig2Window) fill(1);
      for (const std::string& stream : m_in.event.activeTriggers)
      {
        weighted += 2;  // distribution and decision
        PPG12DiagnosticRecord row = occurrence(stream, "TRUTH_ISOLATION_FILL");
        row.decision = passIsolation ? 1 : 2;
        row.generatorBarcode = t.barcode;
        row.flowValue = truthIsolation;
        row.truthPhotonPt = t.pt;
        row.occurrenceWeight = m_weight;
        log(row);
      }
    }

    if (!passIsolation) continue;
    const double truthPt = t.pt;
    if (!std::isfinite(truthPt) || truthPt <= 0.0) continue;

    ownershipAudit("FIG6_OWNERSHIP_AUDIT", 1);
    if (m_photonSample) ownershipAudit("FIG6_OWNERSHIP_AUDIT", 2);
    if (fig6HaveOwner) ownershipAudit("FIG6_OWNERSHIP_AUDIT", 3);
    ownershipAudit("FIG6_OWNERSHIP_AUDIT", fig6Accepted ? 4 : 5);
    efficiencyStage("h_photonEffTruthDen_pTgamma_0", truthPt);
    if (fig6Accepted) efficiencyStage("h_photonEffPpg12Fig6TruthDen_pTgamma_0", truthPt);
    stitchFlow(1, truthPt);
    stitchFlow(2, truthPt);  // the truth photon is a stored Geant4 signal photon

    if (t.trackId >= 0)
    {
      // Raw clusters and photon candidates near the truth photon, and those
      // whose dominant primary is it.
      bool rawNear005 = false, rawNear010 = false, rawTrack = false, rawTrackFiducial = false;
      if (std::isfinite(t.eta) && std::isfinite(t.phi))
        for (const PPG12ClusterRecord& cluster : m_in.clusters)
        {
          if (!std::isfinite(cluster.directionEnergy) || cluster.directionEnergy <= 0.0) continue;
          if (!std::isfinite(cluster.pt) || !std::isfinite(cluster.eta) || !std::isfinite(cluster.phi)) continue;
          const double dr = deltaR(cluster.eta, cluster.phi, t.eta, t.phi);
          if (cluster.pt > 1.0 && dr < 0.05) rawNear005 = true;
          if (cluster.pt > 1.0 && dr < 0.10) rawNear010 = true;
          if (primaryIsSignal(cluster) && cluster.truthTrackId == t.trackId)
          {
            rawTrack = true;
            if (cluster.pt > 5.0 && std::fabs(cluster.eta) < kAcceptanceEtaMax) rawTrackFiducial = true;
          }
        }
      if (rawNear005) stitchFlow(3, truthPt);
      if (rawNear010) stitchFlow(4, truthPt);
      if (rawTrack) stitchFlow(5, truthPt);
      if (rawTrackFiducial) stitchFlow(6, truthPt);

      bool near005 = false, near010 = false, track = false, trackFiducial = false, barcode005 = false;
      if (std::isfinite(t.eta) && std::isfinite(t.phi))
        for (const PhotonRecord& photon : m_in.photons)
        {
          const double eta = photon.eta, phi = TVector2::Phi_mpi_pi(photon.phi), pt = photon.et;
          if (!std::isfinite(eta) || !std::isfinite(phi) || !std::isfinite(pt)) continue;
          const double dr = deltaR(eta, phi, t.eta, t.phi);
          if (pt > 1.0 && dr < 0.05) near005 = true;
          if (pt > 1.0 && dr < 0.10) near010 = true;
          const auto& dominant = photon.dominantTruth;
          const bool primary = dominant.state == DominantTruthState::ValidPrimary ||
                               dominant.state == DominantTruthState::InvalidPrimary;
          if (primary && dominant.pid == 22 && dominant.barcode == t.barcode && dr < 0.05) barcode005 = true;
          if (primaryIsSignal(dominant) && dominant.trackId == t.trackId)
          {
            track = true;
            if (pt > 5.0 && std::fabs(eta) < kAcceptanceEtaMax) trackFiducial = true;
          }
        }
      if (near005) stitchFlow(7, truthPt);
      if (near010) stitchFlow(8, truthPt);
      if (track) stitchFlow(9, truthPt);
      if (trackFiducial) stitchFlow(10, truthPt);
      if (barcode005) stitchFlow(11, truthPt);

      // Fig. 6 efficiency: any unmasked candidate of this truth photon.
      bool fig6Reco = false, fig6Iso = false, fig6Id = false;
      for (const PhotonRecord& photon : m_in.photons)
      {
        if (photon.ppg12TowerMaskState == 1) continue;
        if (!primaryIsSignal(photon.dominantTruth) || photon.dominantTruth.trackId != t.trackId) continue;
        const double eta = photon.eta, et = photon.et;
        if (!std::isfinite(eta) || std::fabs(eta) >= kAcceptanceEtaMax) continue;
        if (!std::isfinite(et) || et < kPhotonEtFloor) continue;
        fig6Reco = true;
        const double isolation = simIsolation(photon);
        if (std::isfinite(isolation) && isolation <= 1e8 && isolation > ppg12::kIsolatedMin &&
            isolation < isoThreshold(et))
          fig6Iso = true;
        if (ppg12::tightness(ppg12::YieldShapes(photon), et) == Tightness::Tight) fig6Id = true;
      }
      if (fig6Accepted && fig6Reco) efficiencyStage("h_photonEffPpg12Fig6Reco_pTgamma_0", truthPt);
      if (fig6Accepted && fig6Reco && fig6Id) efficiencyStage("h_photonEffPpg12Fig6RecoTight_pTgamma_0", truthPt);
      if (fig6Accepted && fig6Reco && fig6Iso) efficiencyStage("h_photonEffPpg12Fig6RecoIso_pTgamma_0", truthPt);
      if (fig6Accepted && fig6Reco && fig6Iso && fig6Id)
        efficiencyStage("h_photonEffPpg12Fig6RecoTightIso_pTgamma_0", truthPt);
    }

    // The reconstructed photon of this truth photon.
    double recoEt = kNaN, recoEta = kNaN, recoPhi = kNaN;
    const PhotonRecord* match = matchedReco(t.trackId, recoEt, recoEta, recoPhi);
    if (!match) continue;
    stitchFlow(12, recoEt);
    if (match->ppg12TowerMaskState == 1) continue;
    const double isolation = simIsolation(*match);
    const bool isolationValid = std::isfinite(isolation) && isolation <= 1e8;

    // Fig. 24: the reconstructed isolation of owned signal photons, with its
    // audit flow (raw fills).
    {
      int flow = 0;
      const auto audit = [&](const int bin)
      {
        (void)bin;
        fillRaw();
        ++flow;
      };
      audit(1);
      int fig24Bin = -1;
      if (std::isfinite(recoEt) && recoEt >= kFig24EtEdges.front() && recoEt < kFig24EtEdges.back())
        for (std::size_t i = 0; i + 1 < kFig24EtEdges.size(); ++i)
          if (recoEt >= kFig24EtEdges[i] && recoEt < kFig24EtEdges[i + 1]) fig24Bin = static_cast<int>(i);
      if (std::isfinite(recoEta) && std::fabs(recoEta) < kAcceptanceEtaMax)
      {
        audit(2);
        if (fig24Bin >= 0)
        {
          audit(3);
          if (std::isfinite(isolation) && isolation < 1e8)
          {
            audit(4);
            if (!m_photonSample)
              audit(12);
            else
            {
              audit(5);
              if (!haveWindow)
                audit(12);
              else
              {
                audit(6);
                if (!fig24HaveOwner)
                  audit(11);
                else if (!fig24Owned)
                  audit(10);
                else
                {
                  audit(7);
                  fill(2);     // compatibility and named ET-isolation maps
                  fillRaw(1);  // audit 8
                  fill(1);     // isolation in its ET bin
                  fillRaw(1);  // audit 9
                }
              }
            }
          }
        }
      }
      (void)flow;
    }

    const int bin = photonBin(recoEt);
    if (bin < 0) continue;
    efficiencyStage("h_photonEffReco_pTgamma_0", truthPt);
    const bool tight = ppg12::tightness(ppg12::YieldShapes(*match), recoEt) == Tightness::Tight;
    if (tight) efficiencyStage("h_photonEffRecoTight_pTgamma_0", truthPt);
    if (!isolationValid) continue;
    const bool iso = isolation > ppg12::kIsolatedMin && isolation < isoThreshold(recoEt);
    if (iso) efficiencyStage("h_photonEffRecoIso_pTgamma_0", truthPt);
    if (tight && iso) efficiencyStage("h_photonEffRecoTightIso_pTgamma_0", truthPt);
    stitchFlow(13, recoEt);

    // Reconstructed isolation of truth-matched signal photons, also per kept
    // stitch definition: weighted, unweighted (raw) and its stitch flow.
    fill(1);
    if (m_photonSample)
      for (std::size_t v = 0; v < m_stitchKept.size(); ++v)
      {
        if (!m_stitchKept[v]) continue;
        fill(1);
        fillRaw(1);
        stitchFlow(14 + static_cast<int>(v), recoEt);
      }
  }

  responseTargets();

  // Cluster-anchored signal leakage: every unmasked signal candidate.
  for (const PhotonRecord& photon : m_in.photons)
  {
    if (photon.ppg12TowerMaskState == 1) continue;
    if (!primaryIsSignal(photon.dominantTruth)) continue;
    const TruthPhoton* truth = signalByTrack(photon.dominantTruth.trackId);
    const double et = photon.et;
    const double response = responseEt(et, truth->row->pt, static_cast<int>(photon.encounterOrdinal));
    if (photonBin(et) < 0) continue;
    const double isolation = simIsolation(photon);
    if (!std::isfinite(isolation) || isolation > 1e8) continue;

    const ppg12::YieldShapes shapes(photon);
    const Tightness tag = ppg12::tightness(shapes, et);
    const bool common = shapes.common();
    const bool iso = isolation > ppg12::kIsolatedMin && isolation < isoThreshold(et);
    const bool nonIso = isolation > nonIsoThreshold(et) && isolation < ppg12::kNonIsolatedMax;
    const bool tagged = tag == Tightness::Tight || tag == Tightness::NonTight;
    const bool window = ppg12::inResponseWindow(et, truth->row->pt);
    yieldAllCommon(et, common, true, true, window, false);
    if (tagged) yieldTightAndAbcd(et, response, truth->row->pt, tag, iso, nonIso, true, true, false);
    if (!tagged || !(iso || nonIso)) continue;
    // Region A only inside the response window; B, C and D by reco ET.
    if (iso && tag == Tightness::Tight && !window) continue;
    fill(1);  // signal ABCD leakage
  }
}

// ---------------------------------------------------------------------------
// Fig. 8 response targets: every input cluster, with an audit flow
// ---------------------------------------------------------------------------

void Bookkeeping::responseTargets()
{
  const auto audit = [&](const int stage)
  {
    for (const std::string& stream : m_in.event.activeTriggers)
    {
      weighted += 1;
      PPG12DiagnosticRecord row = occurrence(stream, "FIG8_AUDIT");
      row.stageCode = stage;
      row.occurrenceWeight = m_weight;
      log(row);
    }
  };
  audit(1);
  if (!std::isfinite(m_vertexZ) || std::fabs(m_vertexZ) >= kTruthFigureVertexMax) return;
  audit(2);
  audit(3);
  audit(m_in.responseClusters ? 4 : 5);  // 5: the unsplit node is absent, the input clusters serve

  const auto& clusters = m_in.responseClusters ? *m_in.responseClusters : m_in.clusters;
  int ordinal = 0;
  for (const PPG12ClusterRecord& cluster : clusters)
  {
    const int clusterOrdinal = ordinal++;
    audit(6);
    if (!primaryIsSignal(cluster)) continue;
    audit(7);
    const TruthPhotonRecord& truth = *signalByTrack(cluster.truthTrackId)->row;
    if (!std::isfinite(truth.pt) || truth.pt <= 0.0 || !std::isfinite(truth.eta) ||
        std::fabs(truth.eta) >= kAcceptanceEtaMax)
      continue;
    const double energy = cluster.energy, eta = cluster.vertexEta;
    const double et = std::isfinite(energy) && std::isfinite(eta) ? energy / std::cosh(eta) : kNaN;
    if (!std::isfinite(et) || et < kFig8RecoEtMin) continue;
    audit(8);
    double truthEnergy = truth.energy;
    if (!std::isfinite(truthEnergy) || truthEnergy <= 0.0) truthEnergy = truth.pt * std::cosh(truth.eta);
    const double responseEt = et / truth.pt, responseE = energy / truthEnergy;
    if (!std::isfinite(responseEt) || !std::isfinite(responseE) || !std::isfinite(truthEnergy) ||
        !std::isfinite(energy))
      continue;
    for (const std::string& stream : m_in.event.activeTriggers)
    {
      weighted += 4;  // ET and E response, E and ET correlations
      PPG12DiagnosticRecord row = occurrence(stream, "FIG8_CLUSTER");
      row.clusterOrdinal = clusterOrdinal;
      row.truthTrackId = cluster.truthTrackId;
      row.generatorBarcode = truth.barcode;
      row.clusterEnergy = energy;
      row.clusterEta = eta;
      row.clusterEt = et;
      row.truthEnergy = truthEnergy;
      row.truthEta = truth.eta;
      row.truthEt = truth.pt;
      row.energyContribution = cluster.energyContribution;
      row.occurrenceWeight = m_weight;
      log(row);
    }
    audit(9);
  }
}

// ---------------------------------------------------------------------------
// Pure truth recoil (h_JES3TruthPure_pT_xJ_alpha): the leading signal truth
// photon and its leading truth recoil jet
// ---------------------------------------------------------------------------

void Bookkeeping::truthRecoil()
{
  const TruthPhoton* lead = nullptr;
  for (const TruthPhoton& truth : m_truth)
    if (truth.signal && (!lead || truth.row->pt > lead->row->pt)) lead = &truth;
  if (!lead || lead->row->pt <= 0.0) return;
  const double tEta = lead->row->eta, tPhi = lead->row->phi;
  for (const Radius& radius : m_radii)
  {
    if (!radius.truthPresent) continue;
    const Jet* first = nullptr;
    for (const Jet& jet : radius.truth)
    {
      if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
      if (jet.pt < m_book.recoilJetPtMinGeV) continue;
      const double dphi = TVector2::Phi_mpi_pi(jet.phi - tPhi), deta = jet.eta - tEta;
      const double dr2 = deta * deta + dphi * dphi;
      if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
      if (!first || jet.pt > first->pt) first = &jet;
    }
    if (!first || first->pt <= 0.0) continue;
    if (!std::isfinite(first->eta) || std::fabs(first->eta) >= jetEtaMax(radius.digits)) continue;
    const double dphi = std::fabs(TVector2::Phi_mpi_pi(first->phi - tPhi));
    if (!std::isfinite(dphi) || dphi < m_book.recoilDeltaPhiMin) continue;
    fill(1);
  }
}

// ---------------------------------------------------------------------------
// Inclusive-jet samples (h_ppInclusiveJetStitch_*)
// ---------------------------------------------------------------------------

void Bookkeeping::inclusiveJetSample()
{
  double maxAny = -1.0, maxR04 = -1.0;
  for (const Radius& radius : m_radii)
    for (const Jet& jet : radius.truth)
      if (std::isfinite(jet.pt))
      {
        if (jet.pt > maxAny) maxAny = jet.pt;
        if (radius.digits == 4 && jet.pt > maxR04) maxR04 = jet.pt;
      }
  const bool haveR04 = maxR04 >= 0.0, haveAny = maxAny >= 0.0;
  const bool passR04 = haveR04 && ppg12::windowContains(maxR04, m_windowLow, m_windowHigh);
  const bool passAny = haveAny && ppg12::windowContains(maxAny, m_windowLow, m_windowHigh);

  PPG12DiagnosticRecord row = occurrence("SIM", "PP_INCLUSIVE_JET_EVENT");
  row.sampleCode = m_weights.ppg12SliceCode;
  if (haveAny) row.sourcePhotonPt = maxAny;
  if (haveR04) row.flowValue = maxR04;
  row.selectionBits =
      (haveAny ? 1ULL : 0ULL) | (passAny ? 2ULL : 0ULL) | (haveR04 ? 4ULL : 0ULL) | (passR04 ? 8ULL : 0ULL);
  log(row);

  // Maximum truth-jet pT (any key, R = 0.4) and the sample code; then the
  // metadata count and, for each maximum that exists, its all and kept or
  // rejected spectra.
  fillRaw((haveAny ? 1 : 0) + (haveR04 ? 1 : 0) + 1);
  fillRaw(1 + (haveR04 ? 2 : 0) + (haveAny ? 2 : 0));
}

// ---------------------------------------------------------------------------
// The photon loop
// ---------------------------------------------------------------------------

void Bookkeeping::yieldAllCommon(const double et, const bool common, const bool haveTruthClass, const bool truthSignal,
                                 const bool window, const bool inclusive)
{
  if (!std::isfinite(et) || et <= 0.0) return;
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    PPG12DiagnosticRecord row = occurrence(stream, "PHOTON_ALL_COMMON");
    row.recoPhotonPt = et;
    row.occurrenceWeight = m_weight;
    row.selectionBits = (common ? 1ULL : 0ULL) | (haveTruthClass ? 2ULL : 0ULL) | (truthSignal ? 4ULL : 0ULL) |
                        (window ? 8ULL : 0ULL) | (inclusive ? 16ULL : 0ULL);
    log(row);
    if (inclusive) weighted += common ? 2 : 1;
    if (haveTruthClass && (!truthSignal || window)) weighted += common ? 2 : 1;
  }
}

void Bookkeeping::yieldTightAndAbcd(const double et, const double response, const double truthPt, const Tightness tag,
                                    const bool iso, const bool nonIso, const bool haveTruthClass,
                                    const bool truthSignal, const bool inclusive)
{
  if (!std::isfinite(et) || et <= 0.0) return;
  const bool tight = tag == Tightness::Tight, nonTight = tag == Tightness::NonTight;
  if (!tight && !nonTight) return;
  const bool analysisWindow = ppg12::inResponseWindow(et, truthPt);
  const bool responseWindow = ppg12::inResponseWindow(response, truthPt);
  const bool truthValid = std::isfinite(truthPt) && truthPt > 0.0;

  for (const std::string& stream : m_in.event.activeTriggers)
  {
    PPG12DiagnosticRecord row = occurrence(stream, "PHOTON_TIGHT_ABCD");
    row.recoPhotonPt = et;
    row.responsePhotonPt = response;
    row.truthPhotonPt = truthPt;
    if (truthValid)
    {
      // PPG12 Fig. 36 truth prior: (p0 + p1 x + p3 x^2) / (1 + p2 x + p4 x^2), x in [10, 30] GeV.
      const double x = std::min(30.0, std::max(10.0, truthPt));
      const double denominator = 1.0 + -0.478818 * x + 0.0522681 * x * x;
      const double prior = (4.05592 + -0.984728 * x + 0.0723232 * x * x) / denominator;
      row.truthPriorWeight = std::isfinite(denominator) && denominator != 0.0 && std::isfinite(prior) ? prior : 1.0;
    }
    row.occurrenceWeight = m_weight;
    row.selectionBits = (tight ? 1ULL : 0ULL) | (nonTight ? 2ULL : 0ULL) | (iso ? 4ULL : 0ULL) |
                        (nonIso ? 8ULL : 0ULL) | (haveTruthClass ? 16ULL : 0ULL) | (truthSignal ? 32ULL : 0ULL) |
                        (analysisWindow ? 64ULL : 0ULL) | (responseWindow ? 128ULL : 0ULL) |
                        (inclusive ? 256ULL : 0ULL);
    log(row);
    if (inclusive && tight) weighted += 1;
    if (tight && haveTruthClass && (!truthSignal || analysisWindow)) weighted += 1;
  }

  const bool region = (tight && iso) || (tight && nonIso) || (nonTight && iso) || (nonTight && nonIso);
  if (!region) return;
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    (void)stream;
    if (inclusive) weighted += 1;
    if (haveTruthClass && (!(truthSignal && tight && iso) || analysisWindow)) weighted += 1;
    if (haveTruthClass && truthSignal && tight && iso && responseWindow)
    {
      // Inside the window both ETs are valid: truth, response-reco and
      // response, their Fig. 36 prior copies and the Fig. 37 inputs; then the
      // half split by the job's event parity (the established jobs read each
      // simulation file from its first entry).
      const bool evenEvent = (m_in.event.sourceEntry + 1) % 2 == 0;
      weighted += 9 + (evenEvent ? 7 : 4);
    }
    else if (haveTruthClass && !truthSignal && tight && iso)
      weighted += 2;  // reco fake, Fig. 37 fake input
  }
}

void Bookkeeping::isolationStack(const double et, const double isolation, const Tightness tag)
{
  if (!std::isfinite(et) || !std::isfinite(isolation) || isolation > 1e8) return;
  if (tag != Tightness::Tight && tag != Tightness::NonTight) return;
  int bin = -1;
  for (std::size_t i = 0; i + 1 < ppg12::kRecoEtBins.size() && bin < 0; ++i)
    if (et > ppg12::kRecoEtBins[i] && et < ppg12::kRecoEtBins[i + 1]) bin = static_cast<int>(i);
  if (bin < 0) return;
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    weighted += 1;
    PPG12DiagnosticRecord row = occurrence(stream, "PP_ISOLATION_STACK_FILL");
    row.sampleCode = bin + 1;
    row.decision = tag == Tightness::Tight ? 1 : 2;
    row.flowValue = isolation;
    row.occurrenceWeight = m_weight;
    log(row);
  }
}

// The ABCD counters of one preselected candidate, per stream.
void Bookkeeping::abcdCounters(const PhotonRecord& photon, const Shapes& v, const Tightness tag, const double et)
{
  if (photonBin(et) < 0) return;
  const double isolation = simIsolation(photon);
  if (!std::isfinite(isolation) || isolation > 1e8) return;
  const bool iso = isolation<isoThreshold(et), nonIso = isolation> nonIsoThreshold(et);
  if (!iso && !nonIso) return;
  if (tag != Tightness::Tight && tag != Tightness::NonTight) return;
  (void)v;
  weighted += 3 + 5;  // count, isolation and pT of the region; five shower variables
}

void Bookkeeping::photonLoop()
{
  const auto& triggers = m_in.event.activeTriggers;
  bool haveLeadIsoTight = false, haveLeadIsoNonTight = false, haveLeadNonIsoTight = false;
  bool haveLeadNonIsoNonTight = false, haveLeadAny = false;
  int leadBin = -1, leadAnyBin = -1;
  double leadEt = -1.0, leadEta = 0.0, leadPhi = 0.0;
  double leadNonTightEt = -1.0, leadNonTightEta = 0.0, leadNonTightPhi = 0.0;
  double leadNonIsoTightEt = -1.0, leadNonIsoNonTightEt = -1.0;
  double leadAnyEt = -1.0, leadAnyPhi = 0.0;
  const PhotonRecord* lead = nullptr;
  const PhotonRecord* leadNonTight = nullptr;
  const PhotonRecord* leadNonIsoTight = nullptr;
  const PhotonRecord* leadNonIsoNonTight = nullptr;
  int isoTightCount = 0;
  std::vector<const PhotonRecord*> isoTight;
  std::vector<double> isoTightEt;

  for (const PhotonRecord& photon : m_in.photons)
  {
    const double eta = photon.eta, phi = photon.phi, et = photon.et;
    if (!std::isfinite(eta) || !std::isfinite(phi) || !std::isfinite(et)) continue;
    if (photon.ppg12TowerMaskState == 1) continue;
    if (et < kPhotonEtFloor) continue;
    if (m_inclusiveSample && std::isfinite(m_weights.clusterEtUpperGeV) && et > m_weights.clusterEtUpperGeV) continue;
    if (std::fabs(eta) >= kAcceptanceEtaMax) continue;
    const int bin = photonBin(et);
    if (bin < 0) continue;
    if (!haveLeadAny || et > leadAnyEt)
    {
      haveLeadAny = true;
      leadAnyBin = bin;
      leadAnyEt = et;
      leadAnyPhi = phi;
    }

    // Isolation QA: spectra and decision.
    const double qaIsolation = simIsolation(photon);
    (void)qaIsolation;
    fill(4);

    const Shapes v(photon, et);
    const bool signal = m_in.evaluatorReady && primaryIsSignal(photon.dominantTruth);
    const double signalTruthPt = signal ? signalByTrack(photon.dominantTruth.trackId)->row->pt : kNaN;

    // The yield path of this candidate.
    const double isolation = simIsolation(photon);
    const ppg12::YieldShapes shapes(photon);
    const bool common = shapes.common();
    const bool haveTruthClass = !signal;
    yieldAllCommon(et, common, haveTruthClass, signal, true, true);
    const bool isolationValid = std::isfinite(isolation) && isolation < 1e8;
    const bool yieldIso = isolationValid && isolation > ppg12::kIsolatedMin && isolation < isoThreshold(et);
    const bool yieldNonIso = isolationValid && isolation > nonIsoThreshold(et) && isolation < ppg12::kNonIsolatedMax;
    const Tightness yieldTag = ppg12::tightness(shapes, et);
    yieldTightAndAbcd(et, et, signalTruthPt, yieldTag, yieldIso, yieldNonIso, haveTruthClass, signal, true);
    if (common) isolationStack(et, isolation, yieldTag);

    // Shower-shape spectra: inclusive, and isolated or non-isolated.
    const int finiteShapes = std::isfinite(v.weta) + std::isfinite(v.wphi) + std::isfinite(v.et1) +
                             std::isfinite(v.e11OverE33) + std::isfinite(v.e32OverE35);
    const bool shapeIso = isolationValid && isolation < isoThreshold(et);
    const bool shapeNonIso = isolationValid && isolation > nonIsoThreshold(et);
    fill(finiteShapes * (1 + (shapeIso ? 1 : 0) + (shapeNonIso ? 1 : 0)));

    // Reference preselection; a failure counts its failed criteria (raw).
    if (!v.preselected())
    {
      const bool failE11 = !(v.e11OverE33 < 0.98), failEt1 = !inOpen(v.et1, 0.6, 1.0);
      const bool failE32 = !inOpen(v.e32OverE35, 0.8, 1.0), failWeta = !(v.weta < 0.6);
      fillRaw(1 + failWeta + failEt1 + failE11 + failE32);
      continue;
    }

    const bool iso = isolationValid && isolation < isoThreshold(et);
    const bool nonIso = isolationValid && isolation > nonIsoThreshold(et);
    const int failures = v.tightFailures();
    const Tightness tag = referenceTag(failures);
    // Tight-failure mask (all, and non-tight or neither) and failed windows (raw).
    fillRaw(1 + (tag != Tightness::Tight ? 1 : 0) + std::popcount(static_cast<unsigned>(failures)));

    // Shower-shape templates: preselected, and tight or non-tight.
    const int templateFills = finiteShapes * (1 + (tag == Tightness::Tight || tag == Tightness::NonTight ? 1 : 0));
    fill(templateFills);

    // Isolation by tightness: total, raw and reconstructed topocluster.
    if (tag == Tightness::Tight || tag == Tightness::NonTight) fill(3);

    for (std::size_t s = 0; s < triggers.size(); ++s) abcdCounters(photon, v, tag, et);

    const bool regionA = iso && tag == Tightness::Tight, regionC = iso && tag == Tightness::NonTight;
    const bool regionB = nonIso && tag == Tightness::Tight, regionD = nonIso && tag == Tightness::NonTight;
    if (regionB && (!haveLeadNonIsoTight || et > leadNonIsoTightEt))
    {
      haveLeadNonIsoTight = true;
      leadNonIsoTightEt = et;
      leadNonIsoTight = &photon;
    }
    if (regionD && (!haveLeadNonIsoNonTight || et > leadNonIsoNonTightEt))
    {
      haveLeadNonIsoNonTight = true;
      leadNonIsoNonTightEt = et;
      leadNonIsoNonTight = &photon;
    }
    if (!regionA && !regionC) continue;
    if (regionA)
    {
      ++isoTightCount;
      fill(1);  // (pT, eta, phi) of isolated tight photons
      isoTight.push_back(&photon);
      isoTightEt.push_back(et);
      if (!haveLeadIsoTight || et > leadEt)
      {
        haveLeadIsoTight = true;
        leadBin = bin;
        leadEt = et;
        leadEta = eta;
        leadPhi = phi;
        lead = &photon;
      }
      continue;
    }
    if (!haveLeadIsoNonTight || et > leadNonTightEt)
    {
      haveLeadIsoNonTight = true;
      leadNonTightEt = et;
      leadNonTightEta = eta;
      leadNonTightPhi = phi;
      leadNonTight = &photon;
    }
  }

  // Event-level photon counts.
  if (haveLeadIsoTight && leadBin >= 0) fillRaw(1);  // isolated tight multiplicity
  const auto leadCount = [&](const bool have, const double et)
  {
    if (have && std::isfinite(et) && et > 0.0 && photonBin(et) >= 0) fill(1);
  };
  leadCount(haveLeadIsoTight, leadEt);
  leadCount(haveLeadNonIsoTight, leadNonIsoTightEt);
  leadCount(haveLeadIsoNonTight, leadNonTightEt);
  leadCount(haveLeadNonIsoNonTight, leadNonIsoNonTightEt);
  if (m_in.evaluatorReady)
  {
    const auto leakage = [&](const bool have, const PhotonRecord* photon, const double et)
    {
      if (!have || !photon || !(std::isfinite(et) && et > 0.0)) return;
      if (!primaryIsSignal(photon->dominantTruth)) return;
      if (photonBin(et) < 0) return;
      fill(1);
    };
    leakage(haveLeadIsoTight, lead, leadEt);
    leakage(haveLeadNonIsoTight, leadNonIsoTight, leadNonIsoTightEt);
    leakage(haveLeadIsoNonTight, leadNonTight, leadNonTightEt);
    leakage(haveLeadNonIsoNonTight, leadNonIsoNonTight, leadNonIsoNonTightEt);
  }

  // The event-leading truth signal photon.
  const TruthPhoton* truthLead = nullptr;
  for (const TruthPhoton& truth : m_truth)
    if (truth.signal && (!truthLead || truth.row->pt > truthLead->row->pt)) truthLead = &truth;
  const bool haveTruthSignal = truthLead != nullptr;
  const double tPtSignal = haveTruthSignal ? truthLead->row->pt : 0.0;
  if (haveTruthSignal)
  {
    leadTag("h_ppg12_fig37_input_leadtag_truth_pT_0", tPtSignal);
    if (!haveLeadIsoTight) leadTag("h_ppg12_fig37_input_leadtag_miss_pT_0", tPtSignal);
    fill(1 + (haveLeadIsoTight ? 0 : 1));  // truth photon spectrum, and its miss
  }
  for (const TruthPhoton& truth : m_truth)
  {
    if (!truth.signal) continue;
    const double pt = truth.row->pt;
    if (!std::isfinite(pt) || pt <= 0.0) continue;
    for (const std::string& stream : triggers)
    {
      weighted += 6;  // truth pT spectra
      PPG12DiagnosticRecord row = occurrence(stream, "PP_TRUTH_SPECTRUM_FILL");
      row.generatorBarcode = truth.row->barcode;
      row.truthPhotonPt = pt;
      row.occurrenceWeight = m_weight;
      log(row);
      weighted += 1 + (haveLeadIsoTight ? 0 : 1);  // per-object truth spectrum, and its miss
    }
  }

  if (haveLeadIsoTight)
  {
    bool haveTruthPhoton = false;
    double tPt = haveTruthSignal ? tPtSignal : 0.0;
    double tEta = haveTruthSignal ? truthLead->row->eta : 0.0;
    double tPhi = haveTruthSignal ? truthLead->row->phi : 0.0;
    std::vector<const PhotonRecord*> matchedObjects;
    if (haveTruthSignal && truthLead->row->trackId >= 0 && m_in.evaluatorReady)
    {
      double matchEt = 0.0, matchEta = 0.0, matchPhi = 0.0;
      // Its tightness and isolation are evaluated; only the isolation lookup
      // leaves a trace.
      if (const PhotonRecord* match = matchedReco(truthLead->row->trackId, matchEt, matchEta, matchPhi))
        isolated(*match, matchEt);
      if (primaryIsSignal(lead->dominantTruth) && lead->dominantTruth.trackId == truthLead->row->trackId)
      {
        const TruthPhotonRecord& t = *signalByTrack(lead->dominantTruth.trackId)->row;
        haveTruthPhoton = true;
        tPt = t.pt;
        tEta = t.eta;
        tPhi = t.phi;
      }
      for (const TruthPhoton& truth : m_truth)
      {
        if (!truth.signal) continue;
        const double pt = truth.row->pt;
        if (!std::isfinite(pt) || pt <= 0.0) continue;
        bool matchedObject = false;
        double et = 0.0, eta = 0.0, phi = 0.0;
        if (const PhotonRecord* match = matchedReco(truth.row->trackId, et, eta, phi))
        {
          const Shapes matched(*match, et);
          const bool tight = classifiedTag(matched) == Tightness::Tight;
          const bool iso = isolated(*match, et);
          if (tight && iso)
          {
            matchedObject = true;
            matchedObjects.push_back(match);
          }
        }
        weighted += static_cast<std::int64_t>(triggers.size()) * (matchedObject ? 2 : 1);
      }
    }
    // Isolated tight candidates matching no signal object are fakes.
    for (std::size_t i = 0; i < isoTight.size(); ++i)
      if (std::find(matchedObjects.begin(), matchedObjects.end(), isoTight[i]) == matchedObjects.end()) fill(1);
    leadPhotonJets(leadBin, leadEt, leadEta, leadPhi, haveTruthSignal, tPtSignal, haveTruthPhoton, tPt, tEta, tPhi);
  }
  else if (haveTruthSignal)
  {
    const std::vector<const Jet*> none;
    const std::vector<bool> noneRecoil;
    for (const Radius& radius : m_radii)
      unfoldingResponse(radius, 0.0, 0.0, 0.0, false, tPtSignal, truthLead->row->eta, truthLead->row->phi, none,
                        noneRecoil, nullptr);
  }

  if (haveLeadIsoNonTight) sidebandRecoil(leadNonTightEt, leadNonTightEta, leadNonTightPhi);

  // Jet-by-jet cutflow against the leading candidate.
  if (haveLeadAny && leadAnyBin >= 0)
    for (const Radius& radius : m_radii)
      for (const Jet& jet : radius.reco)
        if (std::isfinite(jet.pt) && std::isfinite(jet.eta) && std::isfinite(jet.phi)) fill(1);
  (void)leadAnyPhi;
  (void)isoTightCount;
  (void)isoTightEt;
}

void Bookkeeping::leadTag(const char* name, const double value)
{
  if (!std::isfinite(value) || value <= 0.0) return;
  const std::string histogram(name);
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    weighted += 1;
    if (histogram == "h_ppg12_fig37_input_leadtag_truth_pT_0" || histogram == "h_ppg12_fig37_input_leadtag_miss_pT_0")
    {
      PPG12DiagnosticRecord row = occurrence(stream, "PP_LEADTAG_TRUTH_FILL");
      row.decision = histogram == "h_ppg12_fig37_input_leadtag_truth_pT_0" ? 1 : 2;
      row.truthPhotonPt = value;
      row.occurrenceWeight = m_weight;
      log(row);
    }
  }
}

// ---------------------------------------------------------------------------
// The leading isolated tight photon: photon unfolding, then recoil jets per
// radius (matching QA, xJ, alpha, JES and unfolding inputs)
// ---------------------------------------------------------------------------

void Bookkeeping::leadPhotonJets(const int leadBin, const double leadEt, const double leadEta, const double leadPhi,
                                 const bool haveTruthSignal, const double truthSignalPt, const bool haveTruthPhoton,
                                 const double tPt, const double tEta, const double tPhi)
{
  (void)leadBin;
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    (void)stream;
    if (haveTruthPhoton) weighted += 1;                      // reco photon spectrum
    if (!haveTruthPhoton) weighted += 1;                     // reco photon fake
    if (haveTruthSignal && haveTruthPhoton) weighted += 1;   // photon response
    if (haveTruthSignal && !haveTruthPhoton) weighted += 1;  // truth photon miss
  }
  leadTag(haveTruthPhoton ? "h_ppg12_fig37_input_leadtag_reco_pT_0" : "h_ppg12_fig37_input_leadtag_fake_pT_0", leadEt);
  if (haveTruthSignal && haveTruthPhoton)
  {
    if (std::isfinite(truthSignalPt) && truthSignalPt > 0.0 && std::isfinite(leadEt) && leadEt > 0.0) fill(1);
  }
  else if (haveTruthSignal)
    leadTag("h_ppg12_fig37_input_leadtag_miss_pT_0", truthSignalPt);

  for (const Radius& radius : m_radii)
  {
    const double etaMax = jetEtaMax(radius.digits);
    const double jetRadius = 0.1 * radius.digits;
    int passPt = 0, passEta = 0, passDphi = 0;
    std::vector<const Jet*> fid;
    std::vector<bool> fidRecoil;
    const Jet* all1 = nullptr;
    double all1Pt = -1.0, maxDphi = -1.0;
    for (const Jet& jet : radius.reco)
    {
      if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
      if (jet.pt < m_book.recoilJetPtMinGeV) continue;
      ++passPt;
      const bool inEta = std::fabs(jet.eta) < etaMax;
      const double dphiPhoton = TVector2::Phi_mpi_pi(jet.phi - leadPhi), detaPhoton = jet.eta - leadEta;
      const double dr2 = detaPhoton * detaPhoton + dphiPhoton * dphiPhoton;
      if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
      if (jet.pt > all1Pt)
      {
        all1Pt = jet.pt;
        all1 = &jet;
      }
      if (inEta)
      {
        ++passEta;
        const double dphi = std::fabs(TVector2::Phi_mpi_pi(jet.phi - leadPhi));
        if (std::isfinite(dphi) && dphi > maxDphi) maxDphi = dphi;
        const bool recoil = dphi >= m_book.recoilDeltaPhiMin;
        fid.push_back(&jet);
        fidRecoil.push_back(recoil);
        if (recoil) ++passDphi;
      }
    }
    if (all1 && (!std::isfinite(all1->eta) || std::fabs(all1->eta) >= etaMax)) all1 = nullptr;
    const Jet* recoil1 = nullptr;
    if (all1)
    {
      const double dphi = std::fabs(TVector2::Phi_mpi_pi(all1->phi - leadPhi));
      if (std::isfinite(dphi) && dphi >= m_book.recoilDeltaPhiMin) recoil1 = all1;
    }
    (void)passEta;
    (void)passDphi;
    fill(3);               // match status, maximum |dphi|, recoil multiplicity
    if (recoil1) fill(2);  // recoil |dphi|, recoil-is-leading
    if (haveTruthPhoton)
      for (std::size_t i = 0; i < fid.size(); ++i)
        if (fidRecoil[i] && std::isfinite(fid[i]->pt) && fid[i]->pt > 0.0 && std::isfinite(fid[i]->phi)) fill(2);
    unfoldingResponse(radius, leadEt, leadEta, leadPhi, haveTruthPhoton, tPt, tEta, tPhi, fid, fidRecoil, recoil1);
    if (!(recoil1 && recoil1->pt > 0.0)) continue;

    // Selected-jet QA (the subleading fiducial jet away from the photon).
    const Jet* jet2 = nullptr;
    double jet2Pt = 0.0;
    const double minPhotonJetDR = std::max(0.30, jetRadius);
    for (const Jet* jet : fid)
    {
      if (jet == recoil1) continue;
      if (!std::isfinite(jet->pt) || !std::isfinite(jet->eta) || !std::isfinite(jet->phi)) continue;
      if (deltaR(jet->eta, jet->phi, leadEta, leadPhi) <= minPhotonJetDR) continue;
      if (jet->pt > jet2Pt)
      {
        jet2Pt = jet->pt;
        jet2 = jet;
      }
    }
    fill(2 + (jet2 ? 2 : 0));      // selected-jet eta, phi (and of jet 2)
    fill(2);                       // recoil (pT, eta, phi) and its balance profile
    if (haveTruthPhoton) fill(1);  // leading-xJ reco component
    fill(6);                       // xJ, jet 1 pT, jet 2 pT, alpha, JES maps
    if (haveTruthPhoton) truthJes(radius, tPt, tEta, tPhi, recoil1);
  }
}

// h_JES3RecoTruthPhoTagged, h_JES3TruthRecoCondNoJetMatch, h_JES3Truth*, h_JES3RecoTruthTagged
void Bookkeeping::truthJes(const Radius& radius, const double tPt, const double tEta, const double tPhi,
                           const Jet* recoil1)
{
  fill(1);
  if (!radius.truthPresent || tPt <= 0.0) return;
  const Jet* first = nullptr;
  for (const Jet& jet : radius.truth)
  {
    if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
    if (jet.pt < m_book.recoilJetPtMinGeV) continue;
    const double dphi = TVector2::Phi_mpi_pi(jet.phi - tPhi), deta = jet.eta - tEta;
    const double dr2 = deta * deta + dphi * dphi;
    if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
    if (!first || jet.pt > first->pt) first = &jet;
  }
  if (!first || first->pt <= 0.0) return;
  if (!std::isfinite(first->eta) || std::fabs(first->eta) >= jetEtaMax(radius.digits)) return;
  const double dphi = std::fabs(TVector2::Phi_mpi_pi(first->phi - tPhi));
  if (!std::isfinite(dphi) || dphi < m_book.recoilDeltaPhiMin) return;
  fill(1);
  if (deltaR(recoil1->eta, recoil1->phi, first->eta, first->phi) > m_book.jetMatchMaxDeltaR) return;
  fill(3);
}

// Unfolding inputs of one radius: truth recoil jets, jet efficiency, leading
// truth-recoil matching, leading-xJ components and the pairwise response.
void Bookkeeping::unfoldingResponse(const Radius& radius, const double leadEt, const double leadEta,
                                    const double leadPhi, const bool haveTruthPhoton, const double tPt,
                                    const double tEta, const double tPhi, const std::vector<const Jet*>& fid,
                                    const std::vector<bool>& fidRecoil, const Jet* recoil1)
{
  (void)leadEta;
  (void)leadPhi;
  const bool haveRecoPhoton = std::isfinite(leadEt) && leadEt > 0.0;
  const bool haveTruthTarget = std::isfinite(tPt) && tPt > 0.0;
  if (haveTruthTarget && !radius.truthPresent) return;
  const double etaMax = jetEtaMax(radius.digits);
  const double matchDR = m_book.jetMatchMaxDeltaR;
  const auto& xjBins = m_book.unfoldXJBins;
  const double xjMin = xjBins.empty() ? 0.0 : xjBins.front();
  const double xjMax = xjBins.empty() ? std::numeric_limits<double>::infinity() : xjBins.back();

  // Truth recoil candidates.
  std::vector<const Jet*> truthFid;
  std::vector<bool> truthRecoil;
  if (haveTruthTarget)
    for (const Jet& jet : radius.truth)
    {
      if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
      if (jet.pt < m_book.recoilJetPtMinGeV) continue;
      const double dphiPhoton = TVector2::Phi_mpi_pi(jet.phi - tPhi), detaPhoton = jet.eta - tEta;
      const double dr2 = detaPhoton * detaPhoton + dphiPhoton * dphiPhoton;
      if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
      if (std::fabs(jet.eta) >= etaMax) continue;
      const double dphi = std::fabs(TVector2::Phi_mpi_pi(jet.phi - tPhi));
      const bool recoil = dphi >= m_book.recoilDeltaPhiMin;
      truthFid.push_back(&jet);
      truthRecoil.push_back(recoil);
      if (recoil) fill(2);  // truth (pT, xJ) and (pT, |dphi|)
    }

  // Jet efficiency: truth recoil jets reached by an independent reco selection.
  std::vector<bool> effReached(truthFid.size(), false);
  if (haveTruthTarget)
  {
    std::vector<const Jet*> recoEff;
    for (const Jet& jet : radius.reco)
    {
      if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
      if (jet.pt < m_book.recoilJetPtMinGeV) continue;
      const double dphiPhoton = TVector2::Phi_mpi_pi(jet.phi - tPhi), detaPhoton = jet.eta - tEta;
      const double dr2 = detaPhoton * detaPhoton + dphiPhoton * dphiPhoton;
      if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
      if (std::fabs(jet.eta) >= etaMax) continue;
      const double dphi = std::fabs(TVector2::Phi_mpi_pi(jet.phi - tPhi));
      if (!std::isfinite(dphi) || dphi < m_book.recoilDeltaPhiMin) continue;
      const double xj = jet.pt / tPt;
      if (!(xj >= xjMin && xj < xjMax)) continue;
      recoEff.push_back(&jet);
    }
    struct Pair
    {
      double dr;
      int reco, truth;
    };
    std::vector<Pair> pairs;
    for (int r = 0; r < static_cast<int>(recoEff.size()); ++r)
      for (int t = 0; t < static_cast<int>(truthFid.size()); ++t)
      {
        if (!truthRecoil[t]) continue;
        const double dr = deltaR(recoEff[r]->eta, recoEff[r]->phi, truthFid[t]->eta, truthFid[t]->phi);
        if (dr < matchDR) pairs.push_back({dr, r, t});
      }
    std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) { return a.dr < b.dr; });
    std::vector<bool> recoUsed(recoEff.size(), false), truthUsed(truthFid.size(), false);
    for (const Pair& pair : pairs)
    {
      if (recoUsed[pair.reco] || truthUsed[pair.truth]) continue;
      recoUsed[pair.reco] = truthUsed[pair.truth] = true;
      effReached[pair.truth] = true;
    }
  }

  // Leading truth recoil jet and its reconstruction.
  const Jet* truthLead = nullptr;
  bool hasRecoMatch = false, leadMatches = false;
  if (haveTruthTarget && haveRecoPhoton)
  {
    const Jet* all1 = nullptr;
    double all1Pt = -1.0;
    for (const Jet& jet : radius.truth)
    {
      if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
      if (jet.pt < m_book.recoilJetPtMinGeV) continue;
      const double dphiPhoton = TVector2::Phi_mpi_pi(jet.phi - tPhi), detaPhoton = jet.eta - tEta;
      const double dr2 = detaPhoton * detaPhoton + dphiPhoton * dphiPhoton;
      if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
      if (jet.pt > all1Pt)
      {
        all1Pt = jet.pt;
        all1 = &jet;
      }
    }
    if (all1 && all1Pt > 0.0 && std::isfinite(all1->eta) && std::fabs(all1->eta) < etaMax)
    {
      const double dphi = std::fabs(TVector2::Phi_mpi_pi(all1->phi - tPhi));
      if (std::isfinite(dphi) && dphi >= m_book.recoilDeltaPhiMin) truthLead = all1;
    }
    if (truthLead)
    {
      fill(1);  // denominator
      const Jet* best = nullptr;
      double bestDR = 1e9;
      for (const Jet* jet : fid)
      {
        const double dr = deltaR(jet->eta, jet->phi, truthLead->eta, truthLead->phi);
        if (dr < bestDR)
        {
          bestDR = dr;
          best = jet;
        }
      }
      hasRecoMatch = best && bestDR < matchDR;
      bool matchPassesRecoil = false;
      int failCut = 0;
      if (hasRecoMatch)
      {
        const bool okPt = std::isfinite(best->pt) && best->pt >= m_book.recoilJetPtMinGeV;
        const bool okEta = std::isfinite(best->eta) && std::fabs(best->eta) < etaMax;
        const double dphi = std::fabs(TVector2::Phi_mpi_pi(best->phi - tPhi));
        const bool okDphi = std::isfinite(dphi) && dphi >= m_book.recoilDeltaPhiMin;
        matchPassesRecoil = okPt && okEta && okDphi;
        failCut = !okPt ? 1 : !okEta ? 2 : !okDphi ? 3 : 0;
      }
      double leadDR = 1e9;
      if (recoil1)
      {
        leadDR = deltaR(recoil1->eta, recoil1->phi, truthLead->eta, truthLead->phi);
        leadMatches = leadDR < matchDR;
      }
      if (leadMatches)
        fill(1);
      else if (hasRecoMatch)
        fill(2 + (!matchPassesRecoil && failCut > 0 ? 1 : 0));
      else
        fill(1);

      // The selected recoil jet in each class.
      if (recoil1 && std::isfinite(recoil1->pt) && std::isfinite(truthLead->pt) && recoil1->pt > 0.0 &&
          truthLead->pt > 0.0)
      {
        if (leadMatches)
          fill(4 + (hasRecoMatch && best && std::isfinite(best->pt) && best->pt > 0.0 ? 1 : 0));
        else if (hasRecoMatch)
          fill(2 + (best && std::isfinite(best->pt) && best->pt > 0.0 ? 2 : 0) + 3);
        else
          fill(4);
      }
    }
  }

  // Leading-xJ components.
  if (haveRecoPhoton || truthLead)
  {
    const bool validTruthLead =
        truthLead && std::isfinite(tPt) && tPt > 0.0 && std::isfinite(truthLead->pt) && truthLead->pt > 0.0;
    const bool validRecoLead =
        recoil1 && std::isfinite(leadEt) && leadEt > 0.0 && std::isfinite(recoil1->pt) && recoil1->pt > 0.0;
    if (validTruthLead) fill(1);
    if (haveTruthPhoton && validTruthLead && validRecoLead && leadMatches)
      fill(1);
    else
    {
      if (validTruthLead) fill(1);
      if (validRecoLead) fill(1);
      if (haveTruthPhoton && validTruthLead && hasRecoMatch && !leadMatches) fill(1);
    }
  }

  // Pairwise response, one-to-one by increasing Delta R.
  struct Pair
  {
    double dr;
    int reco, truth;
  };
  std::vector<Pair> pairs;
  for (int r = 0; r < static_cast<int>(fid.size()); ++r)
    for (int t = 0; t < static_cast<int>(truthFid.size()); ++t)
    {
      const double dr = deltaR(fid[r]->eta, fid[r]->phi, truthFid[t]->eta, truthFid[t]->phi);
      if (dr < matchDR) pairs.push_back({dr, r, t});
    }
  std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) { return a.dr < b.dr; });
  for (const std::string& stream : m_in.event.activeTriggers)
  {
    (void)stream;
    if (haveTruthTarget)
    {
      const auto& truthBins = m_book.unfoldTruthPhotonEtBins;
      const double ptMin = truthBins.empty() ? 0.0 : truthBins.front();
      const double ptMax = truthBins.empty() ? std::numeric_limits<double>::infinity() : truthBins.back();
      if (tPt >= ptMin && tPt < ptMax)
        for (std::size_t t = 0; t < truthFid.size(); ++t)
        {
          if (!truthRecoil[t]) continue;
          const double xj = truthFid[t]->pt / tPt;
          if (!(xj >= xjMin && xj < xjMax)) continue;
          weighted += 1 + (effReached[t] ? 1 : 0);
        }
    }
    if (!haveTruthPhoton)
    {
      if (haveTruthTarget)
        for (std::size_t t = 0; t < truthFid.size(); ++t)
          if (truthRecoil[t]) weighted += 1;
      if (haveRecoPhoton)
        for (std::size_t r = 0; r < fid.size(); ++r)
          if (fidRecoil[r]) weighted += 1;
      continue;
    }
    std::vector<bool> recoUsed(fid.size(), false), truthUsed(truthFid.size(), false);
    for (const Pair& pair : pairs)
    {
      if (recoUsed[pair.reco] || truthUsed[pair.truth]) continue;
      recoUsed[pair.reco] = truthUsed[pair.truth] = true;
      const Jet* reco = fid[pair.reco];
      const Jet* truth = truthFid[pair.truth];
      const bool recoSelected = fidRecoil[pair.reco], truthSelected = truthRecoil[pair.truth];
      if (truth->pt > 0.0) weighted += 1;  // inclusive matched response
      if (recoil1 && reco == recoil1) weighted += 1 + (truth->pt > 0.0 ? 2 : 0);
      if (truthSelected && recoSelected)
        weighted += 4 + (truth->pt > 0.0 ? 1 : 0);
      else if (truthSelected)
        weighted += 2;
      else if (recoSelected)
        weighted += 2;
    }
    for (std::size_t t = 0; t < truthFid.size(); ++t)
      if (!truthUsed[t] && truthRecoil[t]) weighted += 2;
    for (std::size_t r = 0; r < fid.size(); ++r)
      if (!recoUsed[r] && fidRecoil[r]) weighted += 2;
  }
}

// Sideband C recoil jets (h2_unfoldReco_pTgamma_xJ_incl_sidebandC).
void Bookkeeping::sidebandRecoil(const double leadEt, const double leadEta, const double leadPhi)
{
  (void)leadEt;
  for (const Radius& radius : m_radii)
  {
    const double etaMax = jetEtaMax(radius.digits);
    for (const Jet& jet : radius.reco)
    {
      if (!std::isfinite(jet.pt) || !std::isfinite(jet.eta) || !std::isfinite(jet.phi)) continue;
      if (jet.pt < m_book.recoilJetPtMinGeV) continue;
      const double dphiPhoton = TVector2::Phi_mpi_pi(jet.phi - leadPhi), detaPhoton = jet.eta - leadEta;
      const double dr2 = detaPhoton * detaPhoton + dphiPhoton * dphiPhoton;
      if (!std::isfinite(dr2) || dr2 < kOverlapDeltaR2) continue;
      if (std::fabs(jet.eta) >= etaMax) continue;
      const double dphi = std::fabs(TVector2::Phi_mpi_pi(jet.phi - leadPhi));
      if (!(dphi >= m_book.recoilDeltaPhiMin)) continue;
      if (jet.pt <= 0.0) continue;
      fill(1);
    }
  }
}

// ===========================================================================
// Fig. 6 inclusive-jet truth spectrum
// (h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_*)
//
// Every inclusive-jet event, first of all, before its own weight is configured
// and whatever its reconstructed vertex: the maximum R = 0.4 truth-jet pT,
// weighted by the slice cross section over Jet50's, the truth-vertex weights,
// the SI/DI mixture and the period luminosity. A double-interaction event
// needs its minimum-bias truth vertex. Its occurrence is the event's first;
// its all and kept or rejected fills count in the inclusive-jet lane.
// ===========================================================================

std::int64_t inclusiveJetSpectrum(const SimulationEvent& in, const TH1* vertexReweight,
                                  std::vector<PPG12DiagnosticRecord>& log)
{
  const WeightConfig& w = in.config.weights;
  const EventClassificationConfig& c = in.config.classification;
  double vertexWeight = 1.0;
  if (w.applyVertexReweight)
  {
    const double hard = in.event.truthHardVertexValid ? in.event.truthHardVertexZ : kNaN;
    if (!vertexReweight || !std::isfinite(hard)) return 0;
    vertexWeight *= ppg12::truthVertexWeight(vertexReweight, hard);
    if (w.doubleInteraction)
    {
      const double minimumBias = in.event.truthMinimumBiasVertexValid ? in.event.truthMinimumBiasVertexZ : kNaN;
      if (!std::isfinite(minimumBias) || minimumBias <= -9000.0) return 0;
      vertexWeight *= ppg12::truthVertexWeight(vertexReweight, minimumBias);
    }
  }
  const bool periodContract = !c.period.key.empty();
  const double luminosity =
      c.useLuminosityWeight && periodContract ? c.period.luminosityPbInv / c.period.targetLuminosityPbInv : 1.0;
  const std::array<double, 4> factors{{w.sliceFactor, vertexWeight, c.mixWeight, luminosity}};
  if (!(w.sliceFactor > 0.0) || !std::isfinite(w.sliceFactor)) return 0;
  for (const double factor : factors)
    if (!std::isfinite(factor) || factor < 0.0) return 0;
  const double weight = w.sliceFactor * vertexWeight * c.mixWeight * luminosity;
  if (!std::isfinite(weight) || weight < 0.0) return 0;

  double maxPt = -1.0;
  for (const TruthJetRecord& jet : in.truthJets)
    if (std::fabs(jet.radius - 0.4) < 1e-6 && std::isfinite(jet.pt) && jet.pt > maxPt) maxPt = jet.pt;
  const bool haveJet = maxPt >= 0.0;
  const bool kept = haveJet && ppg12::windowContains(maxPt, w.ownershipLowGeV, w.ownershipHighGeV);

  PPG12DiagnosticRecord row;
  row.trigger = "SIM";
  row.kind = "PP_FIG6_WEIGHTED_SOURCE";
  row.sampleCode = w.ppg12SliceCode;
  row.decision = kept ? 1 : haveJet ? 2 : 3;
  if (haveJet) row.flowValue = maxPt;
  row.occurrenceWeight = weight;
  log.insert(log.begin(), row);

  const bool container =
      std::find(in.truthJetRadiiFound.begin(), in.truthJetRadiiFound.end(), 40) != in.truthJetRadiiFound.end();
  return container && haveJet ? 2 : 0;  // all, and kept or rejected
}
}  // namespace


// ============================================================================
// Entry point
// ============================================================================

/**
 * The pp simulation histogram bookkeeping of one event: an inclusive-jet
 * event's Fig. 6 spectrum, then, once its PPG12 weight is configured, the
 * analysis histograms. Records their occurrences and the weighted and raw fill
 * counts of the lane.
 */
void PhotonJetTree::recordPPG12Simulation()
{
  if (!recordsPPG12Simulation()) return;
  const SimulationEvent in{m_config,
                           m_event,
                           m_ppg12Photons,
                           m_ppg12Clusters,
                           m_ppg12ResponseClustersPresent ? &m_ppg12ResponseClusters : nullptr,
                           m_jets,
                           m_truthPhotons,
                           m_truthJets,
                           m_ppg12TruthJetRadiiFound,
                           m_ppg12EvaluatorReady};
  const int lane = photonjet::production::PPG12WeightLaneCode(m_config.weights);
  if (m_config.weights.ppg12LaneCode == 2)
    m_event.ppg12WeightedFills[lane] += inclusiveJetSpectrum(in, m_vertexReweight, m_event.ppg12Diagnostics);
  if (!m_event.ppg12WeightActive) return;
  Bookkeeping bookkeeping(in, m_event.ppg12Diagnostics);
  bookkeeping.run();
  m_event.ppg12WeightedFills[lane] += bookkeeping.weighted;
  m_event.ppg12RawFills[lane] += bookkeeping.raw;
}
