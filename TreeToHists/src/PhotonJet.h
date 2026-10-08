/* PhotonJet.h: what the TreeToHists fills share, in four parts. */

/* Selection           photon ID and isolation regions, the event-leading photon of each region and its recoil jets
   Input               the eventTree file list, system and metadata, centrality classes, slices and simulation weights,
   the stitching ownership, the recoil-jet thresholds and the native photon E_T and xJ grids
   Resolved selection  the resolved selection of the plan and the event-level validation
   Truth               the truth photon, its truth recoil jets, the reco-truth jet match and the response categories */
#pragma once

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TObjString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <RVersion.h>
#include <TMVA/RBDT.hxx>
/* === Selection =======================================================================================================
   Photon + recoil-jet selection on eventTree v1: photon ID and isolation regions, the event-leading
   photon of each region, and its recoil jets. */

/* One place for every selection decision. The functions take plain values so the ROOT reader
   (FillPhotonJet.C) and any other caller make identical decisions. */
namespace photonjet {

enum class System { pp, auau };
enum Region { kA = 0, kB = 1, kC = 2, kD = 3, kNone = -1 };  /* A tight+iso, B tight+non-iso, C non-tight+iso, D non-tight+non-iso */

struct Line {
  double intercept, slope;
  double at(double x) const { return intercept + slope * x; }
};

/* One half-open routing cell from a saved PhotonID selection definition. */
struct IDCell {
  double etMin, etMax, centMin, centMax;
  Line tight, nonTightHigh;
  bool contains(double et, double cent) const {
    return etMin <= et && et < etMax && centMin <= cent && cent < centMax;
  }
};

/* Every selection parameter. Profiles below are the only sanctioned settings. */
struct Profile {
  std::string name;
  System system = System::pp;
  /* photon */
  double photonEtMin = 15.0, photonEtMax = 35.0, photonAbsEtaMax = 0.7;
  double vertexAbsZMax = 60.0;
  /* photon ID on pho_bdt_score. pp: surfaces in photon ET. Au+Au: surfaces in centrality. */
  Line tight{0.815625, -0.0015625};
  Line nonTightLow{0.7333333333333333, -0.01333333333333333};
  Line nonTightHigh{0.684375, 0.0015625};
  bool nonTightLowIsConstant = false;  /* Au+Au nominal: 0.1 */
  double nonTightLowConstant = 0.1;
  /* isolation witness: index 0 is R = 0.3, 1 is R = 0.4 (meta rules isolation_witnesses) */
  int isolationRadiusIndex = 1;
  int isolationSubtraction = 1;  /* producer isolation method: 1 pp topocluster, 2 Au+Au calorimeter SUB1 */
  /* recoil jets (the analysis collection jet04 = anti-kT R = 0.4) */
  double jetPtMin = 5.0, jetAbsEtaMax = 0.7, overlapDeltaR = 0.4, deltaPhiMin = 7.0 * M_PI / 8.0;
  /* Regions from the augmented classes pho_id_class and pho_iso_class, computed once at conversion from the frozen
     photon-ID and isolation packages, instead of the score surfaces and the stored witness thresholds. */
  bool augmentedClasses = false;
  /* ABCD sideband variations (profile variants, see profileVariant): NaN or 0 is the profile's own definition. */
  double nonTightMin = std::numeric_limits<double>::quiet_NaN();  /* Au+Au: the non-tight lower bound on the score */
  double nonTightLowShift = 0.0;                                  /* pp: shift of the non-tight lower score line */
  double sidebandGap = std::numeric_limits<double>::quiet_NaN();  /* GeV from the isolation threshold to the sideband */
  /* Resolved YAML selections evaluate the primitives, never historical producer labels. */
  bool configured = false, isolationFromWitness = true, verifyNominalID = false, verifyNominalIsolation = false;
  Line isolationThreshold{0.0, 0.0};  /* pp: ET, AuAu: centrality */
  double isolationCentMin = 0.0, isolationCentMax = 80.0;
  std::string scoreBranch = "pho_bdt_score", stateBranch = "pho_bdt_state";
  std::string selectionJson, selectionId;
  bool optionalID = false;  /* finite score is applicability; no extra state branch */
  bool nativeID = false;    /* evaluate the sealed ROOT model during the existing read */
  std::vector<IDCell> idCells;

};

/* Persisted photon-ID and isolation class codes from the eventTree contract. */
constexpr int kIdTight = 3, kIdNonTight = 1, kIdBelowNonTight = 0, kIsoIsolated = 2, kIsoGap = 1, kIsoNonIsolated = 0;

/* The historical reference analysis (schema-10 package v3): R = 0.4 isolation with the stored per-photon
   thresholds, Au+Au non-tight WP90 < s < WP80, pp |z| < 30 cm. Used only to prove this code reproduces it. */
inline Profile historical(System system) {
  Profile p;
  p.name = "historical";
  p.system = system;
  p.isolationRadiusIndex = 1;
  if (system == System::pp) {
    p.vertexAbsZMax = 30.0;
    p.isolationSubtraction = 1;
  } else {
    p.vertexAbsZMax = 10.0;
    p.isolationSubtraction = 2;
    p.tight = {0.6927385742113304, 0.0010476808935335573};         /* WP70 */
    p.nonTightHigh = {0.6071724560222583, 0.0013019695025559648};  /* WP80 */
    p.nonTightLow = {0.4634577166843198, 0.0015466253667456042};   /* WP90 */
  }
  return p;
}

/* The nominal selection. pp is the validated PPG12-parity selection (NPB preselection, PPG12 BDT tight and
   non-tight bands, R = 0.4 isolation below 0.490 + 0.037 ET with the sideband 0.8 GeV above it, |z| < 60 cm):
   the historical selection with the 60 cm vertex window. Au+Au takes its regions from the augmented classes: tight
   s > WP70(c) and non-tight 0.1 < s < WP80(c) (pho_id_class), isolated below the R = 0.3 package line and non-isolated
   above it plus the gap (pho_iso_class); the other photon and event cuts are the historical ones. */
inline Profile nominal(System system) {
  Profile p = historical(system);
  p.name = "nominal";
  if (system == System::pp) {
    p.vertexAbsZMax = 60.0;
  } else {
    p.augmentedClasses = true;
    p.isolationRadiusIndex = 0;  /* R = 0.3, applied through pho_iso_class */
  }
  return p;
}

inline double wrappedAbsDeltaPhi(double a, double b) {
  /* |dphi| in [0, pi], the same arithmetic as the reference: fmod into [-pi, pi) then abs */
  double d = std::fmod(a - b + M_PI, 2.0 * M_PI);
  if (d < 0) d += 2.0 * M_PI;
  return std::fabs(d - M_PI);
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return std::hypot(eta1 - eta2, wrappedAbsDeltaPhi(phi1, phi2));
}

/* The photon-ID and isolation classes of photon i, which already passed the kinematic selection: the one definition of
   the regions, used for the region leaders and the QA stages. score is pho_bdt_score (finite when pho_bdt_state == 1);
   cone, threshold, sideband, method and role are the isolation witnesses, two per photon (R = 0.3, then R = 0.4).
   With augmented classes the codes decide, and a variation reads the stored score (non-tight lower bound) or joins the
   gap code to the non-isolated side (sideband gap 0), so no package edge is recomputed. */
struct PhotonClass {
  bool scored = false;   /* a photon-ID score, or with augmented classes always */
  bool witness = false;  /* a reconstructed witness of the profile's radius and method with a finite cone sum */
  bool tight = false, nonTight = false, isolated = false, nonIsolated = false;
  int idCode = -1, isoCode = -1;
  Region region() const {
    if (tight && isolated) return kA;
    if (tight && nonIsolated) return kB;
    if (nonTight && isolated) return kC;
    if (nonTight && nonIsolated) return kD;
    return kNone;
  }
};

template <class Values, class Ints, class Scores>
PhotonClass photonClass(const Profile& p, std::size_t i, double et, double centrality, const Scores& score, const Ints& state,
                        const Values& cone, const Values& threshold, const Values& sideband, const Ints& method,
                        const Ints& role, const Ints& idClass, const Ints& isoClass) {
  PhotonClass c;
  const std::size_t w = 2 * i + p.isolationRadiusIndex;
  c.witness = role[w] == 0 && method[w] == p.isolationSubtraction && std::isfinite(cone[w]);
  if (p.configured) {
    c.scored = (p.optionalID || state[i] == 1) && std::isfinite(score[i]);
    if (c.scored) {
      const double x = p.system == System::pp ? et : centrality;
      double low = p.nonTightLow.at(x), high = p.nonTightHigh.at(x), tight = p.tight.at(x);
      if (p.optionalID) {
        const IDCell* selected = nullptr;
        for (const auto& cell : p.idCells) if (cell.contains(et, p.system == System::pp ? 0.0 : centrality)) {
          if (selected) throw std::runtime_error("optional photon-ID cells overlap");
          selected = &cell;
        }
        if (!selected) throw std::runtime_error("finite optional score lies outside its routed domain");
        low = 0.1; high = selected->nonTightHigh.at(x); tight = selected->tight.at(x);
      }
      c.idCode = score[i] > tight ? kIdTight : low < score[i] && score[i] < high ? kIdNonTight
                               : score[i] >= high ? 2 : kIdBelowNonTight;
      c.tight = c.idCode == kIdTight;
      c.nonTight = c.idCode == kIdNonTight;
    }
    if (p.system == System::auau && !(std::isfinite(centrality) && p.isolationCentMin <= centrality && centrality < p.isolationCentMax))
      c.witness = false;
    const double cut = p.isolationFromWitness ? threshold[w]
                       : p.isolationThreshold.at(p.system == System::pp ? et : centrality);
    /* When a pp threshold is overridden, move its sideband with it; retaining
       the old absolute upper edge would silently change the requested gap. */
    const double upper = std::isnan(p.sidebandGap)
                             ? (p.isolationFromWitness ? sideband[w] : cut + (sideband[w] - threshold[w]))
                             : cut + p.sidebandGap;
    c.witness = c.witness && std::isfinite(cut) && std::isfinite(upper) && upper >= cut;
    if (c.witness) {
      c.isoCode = cone[w] < cut ? kIsoIsolated : cone[w] > upper ? kIsoNonIsolated : kIsoGap;
      c.isolated = c.isoCode == kIsoIsolated;
      c.nonIsolated = c.isoCode == kIsoNonIsolated;
    }
    return c;
  }
  if (p.augmentedClasses) {
    c.scored = true;
    c.tight = idClass[i] == kIdTight;
    c.nonTight = std::isnan(p.nonTightMin)
                     ? idClass[i] == kIdNonTight
                     : (idClass[i] == kIdNonTight || idClass[i] == kIdBelowNonTight) && std::isfinite(score[i]) &&
                           score[i] > p.nonTightMin;  /* both codes lie below the non-tight upper bound */
    c.isolated = isoClass[i] == kIsoIsolated;
    c.nonIsolated = isoClass[i] == kIsoNonIsolated || (p.sidebandGap == 0.0 && isoClass[i] == kIsoGap);
    return c;
  }
  if (state[i] != 1 || !std::isfinite(score[i])) return c;
  c.scored = true;
  double tightCut, low, high;
  if (p.system == System::pp) {
    tightCut = p.tight.at(et);
    low = p.nonTightLow.at(et) + p.nonTightLowShift;
    high = p.nonTightHigh.at(et);
  } else {
    tightCut = p.tight.at(centrality);
    high = p.nonTightHigh.at(centrality);
    low = !std::isnan(p.nonTightMin) ? p.nonTightMin
                                      : (p.nonTightLowIsConstant ? p.nonTightLowConstant : p.nonTightLow.at(centrality));
  }
  c.tight = score[i] > tightCut;
  c.nonTight = low < score[i] && score[i] < high;
  if (c.witness && std::isfinite(threshold[w]) && std::isfinite(sideband[w])) {
    c.isolated = cone[w] < threshold[w];
    c.nonIsolated = cone[w] > (std::isnan(p.sidebandGap) ? sideband[w] : threshold[w] + p.sidebandGap);
  }
  return c;
}

inline bool photonKinematics(const Profile& p, double vz, double et, double eta, double phi) {
  if (!(std::isfinite(vz) && std::isfinite(et) && std::isfinite(eta) && std::isfinite(phi))) return false;
  if (std::fabs(vz) >= p.vertexAbsZMax) return false;
  if (!(p.photonEtMin <= et && et < p.photonEtMax)) return false;
  return std::fabs(eta) < p.photonAbsEtaMax;
}

/* The event-leading photon of each region (A, B, C, D), -1 when the region has none. A photon enters when it has
   active preselection, passes the kinematic selection, has a photon-ID score (pho_bdt_state 1) and a reconstructed
   isolation witness of the profile's radius and method; with augmented classes, when its classes give a region. The
   leader is the highest-E_T photon of its region; an exact tie goes to the later photon, as the reference's (hi, lo)
   identity order does. Photon arrays hold one value per photon; isolation arrays hold two per photon (R = 0.3, then
   R = 0.4). idClass and isoClass (pho_id_class, pho_iso_class) are read only with augmented classes. et may be of
   another indexable type than the other arrays (a varied energy scale). */
/* One event-local photon index per ABCD region; -1 means no eligible leader.
   Equal ET keeps the later encounter, matching the established >= rule.
   All candidate arrays must share the same photon ordering. */
template <class Energies, class Values, class Ints, class Scores>
std::array<long, 4> regionLeaders(const Profile& p, double vz, double centrality, std::size_t n, const Energies& et,
                                  const Values& eta, const Values& phi, const Scores& score, const Ints& presel,
                                  const Ints& state, const Values& cone, const Values& threshold, const Values& sideband,
                                  const Ints& method, const Ints& role, const Ints& idClass, const Ints& isoClass) {
  std::array<long, 4> leader{-1, -1, -1, -1};
  for (std::size_t i = 0; i < n; ++i) {
    if (presel[i] != 1) continue;
    if (!photonKinematics(p, vz, et[i], eta[i], phi[i])) continue;
    const Region r =
        photonClass(p, i, et[i], centrality, score, state, cone, threshold, sideband, method, role, idClass, isoClass).region();
    if (r == kNone) continue;
    if (leader[r] < 0 || et[i] >= et[leader[r]]) leader[r] = static_cast<long>(i);
  }
  return leader;
}

/* Kinematics use GeV and radians. The same recoil decision feeds measurement,
   response and QA, so changing a cut belongs in the resolved Profile. */
inline bool recoilJet(const Profile& p, double photonEta, double photonPhi, double pt, double eta, double phi) {
  if (!(std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi))) return false;
  if (pt < p.jetPtMin) return false;
  if (deltaR(photonEta, photonPhi, eta, phi) < p.overlapDeltaR) return false;
  if (std::fabs(eta) >= p.jetAbsEtaMax) return false;
  return wrappedAbsDeltaPhi(phi, photonPhi) >= p.deltaPhiMin;
}

}  /* namespace photonjet */

/* === Input ===========================================================================================================
   What every TreeToHists reader shares: the input file list, the system and metadata of an eventTree file, the
   centrality classes, slices and simulation weights, the recoil-jet thresholds and the native photon E_T and xJ grids. */
namespace photonjet {

/* Set only by the versioned public lifecycle, before booking any histograms.
   Empty keeps the historical macro interface and its reference axes unchanged. */
inline std::vector<double> configuredXjEdges;
inline std::vector<double> configuredMeasurementEtEdges;
inline bool strictReduction = false;
inline std::string progressPath;
inline long long progressTotal = 0;
/* A production row reads a half-open range of the original tree. Entry numbers
   remain global, including response-smearing seeds and saved event identities. */
inline Long64_t entryBegin = 0, entryEnd = -1;
inline void applyEntryRange(TTreeReader &reader) {
  if (entryEnd < 0) return;
  if (entryBegin < 0 || entryEnd < entryBegin || entryEnd > reader.GetEntries() ||
      (entryBegin==entryEnd && reader.GetEntries()!=0))
    throw std::runtime_error("invalid eventTree entry range");
  if (reader.GetEntries()==0) return;
  reader.SetEntriesRange(entryBegin, entryEnd);
}
inline nlohmann::json progressIdentity;
inline std::string progressPhase;
inline long long progressOffset = 0, progressRowTotal = 0, progressRowProcessed = 0;
inline auto progressStart = std::chrono::steady_clock::now();
inline auto progressLast = progressStart;
/* A grouped worker mirrors the same in-process counts into its row record;
   it does not start a separate polling process. */
inline std::function<void(const nlohmann::json &)> progressObserver;
/* The event loop owns progress. Count all passes and refresh at most every
   two seconds or 10,000 entries; phase boundaries are written immediately. */
inline void reportProgress(long long processed) {
  if (progressPath.empty()) return;
  processed -= entryBegin;
  if (processed < 0 || processed > progressRowTotal)
    throw std::runtime_error("progress exceeds the planned entry range");
  progressRowProcessed = processed;
  const auto now=std::chrono::steady_clock::now();
  if (processed != 0 && processed != progressRowTotal && processed % 10000 != 0 &&
      std::chrono::duration<double>(now-progressLast).count() < 2) return;
  progressLast = now;
  const long long completed = progressOffset + processed;
  const double seconds=std::chrono::duration<double>(now-progressStart).count();
  const double rate=seconds>0 ? completed/seconds : 0;
  nlohmann::json record = progressIdentity;
  record.update({{"phase",progressPhase}, {"processed",completed}, {"total",progressTotal},
    {"remaining",progressTotal-completed}, {"row_processed",processed}, {"row_total",progressRowTotal},
    {"row_remaining",progressRowTotal-processed}, {"deterministic_skip_processed",0},
    {"deterministic_skip_total",0}, {"deterministic_skip_remaining",0},
    {"elapsed_seconds",seconds}, {"rate",rate}, {"eta_seconds",rate>0 ? (progressTotal-completed)/rate : -1},
    {"updated_at_unix",std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},
    {"latest_error",nullptr}});
  std::ofstream out(progressPath+".next");
  out << record.dump();
  if (!out) throw std::runtime_error("cannot write fill progress");
  out.close(); std::filesystem::rename(progressPath+".next",progressPath);
  if (progressObserver) progressObserver(record);
}

/* The stitching ownership of the job's sample (Au+Au embedded simulation, config/auau_embedded_stitching.json), set by
   the driver before the fills: an inclusive-jet sample owns the events whose leading R = 0.4 truth jet pT lies in its
   window, a photon+jet sample the events whose generator-filter photon pT lies in it ([low, high)). An event it does
   not own belongs to another sample and is left out of every fill. Owner::none: every event is the sample's. */
struct Ownership {
  enum class Owner { none, leadingTruthJet, filterPhoton } owner = Owner::none;
  double low = 0, high = std::numeric_limits<double>::infinity();
  const char* branch() const { return owner == Owner::leadingTruthJet ? "truthjet04_pt" : "embdiag_stitch_photon_pt"; }
  /* The event's owner variable from its branch (NaN when the event has none): the largest finite positive truth-jet
     pT, or the generator-filter photon pT of the event's stitch record (the other diagnostic records carry none). */
  double value(const TTreeReaderArray<double>& v) const {
    double out = std::numeric_limits<double>::quiet_NaN();
    for (size_t i = 0; i < v.GetSize(); ++i) {
      if (!std::isfinite(v[i])) continue;
      if (owner == Owner::leadingTruthJet) {
        if (v[i] > 0 && !(v[i] <= out)) out = v[i];
      } else if (std::isnan(out)) {
        out = v[i];
      } else if (v[i] != out) {
        throw std::runtime_error("one event carries two different generator-filter photons");
      }
    }
    return out;
  }
  bool contains(double x) const { return x >= low && x < high; }
};
inline Ownership configuredOwnership;

/* The ownership test of one fill's reader, made beside its other branches (before the first Next()). */
class OwnershipGate {
 public:
  explicit OwnershipGate(TTreeReader& reader, const Ownership& ownership = configuredOwnership,
                         bool strict = strictReduction) : ownership_(ownership), strict_(strict) {
    if (ownership_.owner != Ownership::Owner::none)
      values_.reset(new TTreeReaderArray<double>(reader, ownership_.branch()));
    if (strict_ && ownership_.owner == Ownership::Owner::filterPhoton) {
      record_.reset(new TTreeReaderValue<int>(reader, "evt_record"));
      nodes_.reset(new TTreeReaderValue<int>(reader, "evt_nodes_ready"));
      terminal_.reset(new TTreeReaderValue<int>(reader, "evt_terminal_status"));
    }
  }
  bool owned() const { return !values_ || ownership_.contains(value()); }
  /* A photon source rejected before ownership still contributes one event
     to its all-events denominator, but has no photon witness. */
  bool rejectedBeforeOwnership() const {
    return record_ && (**record_ == 1 || **nodes_ == 0 || (**terminal_ >= 1 && **terminal_ <= 3));
  }
  double value() const {
    const double x = values_ ? ownership_.value(*values_) : std::numeric_limits<double>::quiet_NaN();
    if (strict_ && ownership_.owner == Ownership::Owner::filterPhoton && !std::isfinite(x)) {
      /* Archived AuAu trigger/vz/minimum-bias rejection occurs before stitching.
         These events legitimately have no occurrence. A native ungated retained
         event must instead carry a finite witness, including the explicit -1
         sentinel when the generator has no qualifying photon. */
      if (!rejectedBeforeOwnership())
        throw std::runtime_error("retained event lacks its generator-filter photon witness");
    }
    return x;
  }

 private:
  Ownership ownership_;
  bool strict_;
  std::unique_ptr<TTreeReaderArray<double>> values_;
  std::unique_ptr<TTreeReaderValue<int>> record_, nodes_, terminal_;
};

struct CentralityClass {
  std::string name;
  double low, high;  /* low <= centrality < high; ignored in pp */
};

inline std::vector<CentralityClass> centralityClasses(System system) {
  if (system == System::auau) return {{"cent0_20", 0, 20}, {"cent20_50", 20, 50}, {"cent50_80", 50, 80}};
  return {{"all", 0, 0}};
}

/* The 5% centrality slice of an Au+Au event, floor(c / 5): 0..15 inside the classes (0-80%). */
inline int centralitySlice(double centrality) { return static_cast<int>(std::floor(centrality / 5.0)); }

/* Au+Au simulation centrality weights, one per 5% slice of 0-80% and 1 outside: w(c) = p_DATA(c) / p_SIM(c), which
   weights the simulated centrality distribution to the data's. They are derived for each production from its data
   and simulation centrality distributions, so they are an input. Every simulation fill is weighted by w_event x w(c). */
struct CentralityWeights {
  std::vector<double> perSlice;  /* empty: no centrality weight (pp) */
  double operator()(double centrality) const {
    if (perSlice.empty() || !(centrality >= 0.0 && centrality < 80.0)) return 1.0;
    return perSlice[centralitySlice(centrality)];
  }
};

/* Parse the weights of a system: 16 finite positive numbers for Au+Au, none for pp. */
inline CentralityWeights centralityWeights(const std::string& text, System system) {
  CentralityWeights w;
  std::istringstream stream(text);
  for (double v; stream >> v;) {
    if (!(std::isfinite(v) && v > 0.0)) throw std::runtime_error("centrality weights must be finite and positive");
    w.perSlice.push_back(v);
  }
  if (!stream.eof()) throw std::runtime_error("centrality weights are not numbers: " + text);
  if (system == System::pp && !w.perSlice.empty()) throw std::runtime_error("pp takes no centrality weights");
  if (system == System::auau && w.perSlice.size() != 16)
    throw std::runtime_error("Au+Au simulation needs 16 centrality weights, one per 5% slice of 0-80%");
  return w;
}

const std::vector<double> kJetThresholds = {5.0, 7.0, 10.0, 12.0};

inline std::string thresholdName(double t) { return "jet" + std::to_string(static_cast<int>(t)); }

inline std::vector<double> photonEtEdges() {
  if (!configuredMeasurementEtEdges.empty()) return configuredMeasurementEtEdges;
  std::vector<double> e;
  for (int v = 15; v <= 35; ++v) e.push_back(v);  /* native 1 GeV grid */
  return e;
}

inline std::vector<double> xjEdges() {
  if (!configuredXjEdges.empty()) return configuredXjEdges;
  std::vector<double> e;
  for (int i = 0; i <= 40; ++i) e.push_back(std::round(i * 0.05 * 1e6) / 1e6);  /* 0.05 steps to 2.0 */
  e.push_back(3.0);
  return e;
}

/* |Delta eta| above which a jet is on the far side of the photon: the shape of the data-driven combinatoric estimate. */
const double kFarDeltaEta = 0.8;

/* pi/32 steps with 0.8 and 2.5 rad added, so the recoil window (dphi >= 7pi/8) and the sideband and fit
   ranges of the combinatorial subtraction are bin edges. */
inline std::vector<double> deltaPhiEdges() {
  std::vector<double> e;
  for (int k = 0; k <= 32; ++k) e.push_back(k * M_PI / 32.0);
  e.push_back(0.8);
  e.push_back(2.5);
  std::sort(e.begin(), e.end());
  return e;
}

inline std::vector<std::string> inputPaths(const char* inputs) {
  std::vector<std::string> paths;
  std::map<std::string, std::string> pinnedGroups;
  std::istringstream stream(inputs);
  for (std::string path; stream >> path;) {
    const auto at = path.find("/CURRENT/");
    if (at != std::string::npos) {
      const std::string group = path.substr(0, at);
      auto found = pinnedGroups.find(group);
      if (found == pinnedGroups.end())
        found = pinnedGroups.emplace(group, std::filesystem::canonical(group + "/CURRENT").string()).first;
      path = found->second + path.substr(at + 8);
    }
    paths.push_back(path);
  }
  if (paths.empty()) throw std::runtime_error("no input files");
  return paths;
}

inline std::string metaText(TFile& f, const std::string& key) {
  auto* text = dynamic_cast<TObjString*>(f.Get(key.c_str()));
  if (!text) throw std::runtime_error("no " + key + " in " + f.GetName());
  return text->GetString().Data();
}

/* The system of an eventTree file (meta/schema is compact JSON, Format tables.dumps). */
inline System fileSystem(TFile& f) {
  const std::string text = metaText(f, "meta/schema");
  if (text.find("\"system\":\"auau\"") != std::string::npos) return System::auau;
  if (text.find("\"system\":\"pp\"") != std::string::npos) return System::pp;
  throw std::runtime_error(std::string("unknown system in ") + f.GetName());
}

inline System fileSystem(const std::string& path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  return fileSystem(*f);
}

/* True when the file is simulation: meta/source dataset isSim, isSimInclusive, isSimEmbedded or isSimEmbeddedInclusive. */
inline bool isSimulation(TFile& f) {
  const std::string source = metaText(f, "meta/source");
  for (const char* dataset : {"isSim", "isSimInclusive", "isSimEmbedded", "isSimEmbeddedInclusive"})
    if (source.find(std::string("\"dataset\":\"") + dataset + "\"") != std::string::npos) return true;
  return false;
}

/* The identity of a file's production input: meta/source input_file_sha256. */
inline std::string sourceSha(TFile& f) {
  const std::string source = metaText(f, "meta/source");
  /* Consolidated/augmented versions share one original logical source identity.
     Their container hash changes with augmentation and must not become a new
     normalization or resampling unit. Original chunk hashes remain in provenance. */
  const std::string logical = "\"logical_source_id\":\"";
  const std::string key = source.find(logical) != std::string::npos ? logical : "\"input_file_sha256\":\"";
  const size_t at = source.find(key);
  if (at == std::string::npos) throw std::runtime_error(std::string("no input_file_sha256 in meta/source of ") + f.GetName());
  const size_t begin = at + key.size(), end = source.find('"', begin);
  return source.substr(begin, end - begin);
}

/* {"path":...,"input_file_sha256":...} of one input. Opening the file moves the current ROOT directory, so outputs
   call this before they open their own file. */
inline std::string sourceRecord(const std::string& path) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  return "{\"path\":\"" + path + "\",\"input_file_sha256\":\"" + sourceSha(*file) + "\"}";
}

/* The JSON list of sourceRecord for every input. */
inline std::string sourceRecords(const std::vector<std::string>& paths) {
  std::string out = "[";
  for (size_t i = 0; i < paths.size(); ++i) out += (i ? "," : "") + sourceRecord(paths[i]);
  return out + "]";
}

/* A profile by name: "nominal", "historical", or a variant of one for the ABCD sideband systematic,
   "<base>:key=value[,key=value]", e.g. "nominal:nontight_min=0.2" or "nominal:sideband_gap=0". Keys:
   nontight_min        Au+Au: the non-tight lower bound on the score (nominal 0.1); with augmented classes it is read
   from pho_id_class and the stored score, so the non-tight upper bound stays the package's
   nontight_low_shift  pp: a shift of the non-tight lower score line
   sideband_gap        GeV from the isolation threshold to the non-isolated sideband (pp nominal: the stored sideband);
   with augmented classes only 0, which joins the gap class to the non-isolated side
   The name, variant included, is written to every output, and FinalAnalysis requires every file of a result to carry
   the configured one. */
inline Profile namedProfile(const std::string& name, System system) {
  const size_t colon = name.find(':');
  const std::string base = name.substr(0, colon);
  Profile p;
  if (base == "nominal") {
    p = nominal(system);
  } else if (base == "historical") {
    p = historical(system);
  } else {
    throw std::runtime_error("unknown profile " + name);
  }
  if (colon == std::string::npos) return p;
  const std::string variant = name.substr(colon + 1);
  std::vector<std::string> seen;
  size_t start = 0;
  while (true) {
    const size_t comma = variant.find(',', start);
    const std::string item = variant.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
    const size_t at = item.find('=');
    if (at == std::string::npos) throw std::runtime_error("a profile variant item is key=value: '" + item + "' in " + name);
    const std::string key = item.substr(0, at), value = item.substr(at + 1);
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) throw std::runtime_error("a profile variant key appears twice: " + name);
    seen.push_back(key);
    size_t used = 0;
    double x = 0.0;
    try {
      x = std::stod(value, &used);
    } catch (const std::exception&) {
      used = 0;
    }
    if (used == 0 || used != value.size() || !std::isfinite(x))
      throw std::runtime_error("a profile variant value is one finite number: " + item);
    if (key == "nontight_min") {
      if (system != System::auau) throw std::runtime_error("nontight_min is the Au+Au bound; pp shifts its line (nontight_low_shift)");
      if (!(0.0 <= x && x < 1.0)) throw std::runtime_error("nontight_min is a score in [0, 1): " + item);
      p.nonTightMin = x;
    } else if (key == "nontight_low_shift") {
      if (system != System::pp) throw std::runtime_error("nontight_low_shift is the pp line shift; Au+Au sets nontight_min");
      if (!(std::fabs(x) < 1.0)) throw std::runtime_error("nontight_low_shift is a score shift below 1 in size: " + item);
      p.nonTightLowShift = x;
    } else if (key == "sideband_gap") {
      if (x < 0.0) throw std::runtime_error("sideband_gap is not negative: " + item);
      if (p.augmentedClasses && x != 0.0)
        throw std::runtime_error("with augmented isolation classes the sideband gap can only be 0; another gap needs its own isolation package");
      p.sidebandGap = x;
    } else {
      throw std::runtime_error("unknown profile variant key " + key + " in " + name);
    }
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  p.name = name;
  return p;
}

}  /* namespace photonjet */

/* === Resolved selection ==============================================================================================
   The resolved selection (selection.json, resolved by the plan) and the event-level validation of the fills. */
namespace photonjet {
using SelectionJson = nlohmann::json;

inline Profile configuredProfile(const char* legacyName, System system, const char* resolvedPath) {
  if (!resolvedPath || !*resolvedPath) return namedProfile(legacyName, system);
  if (std::string(legacyName) != "nominal")
    throw std::runtime_error("use one resolved selection; do not combine it with a legacy profile variant");
  std::ifstream stream(resolvedPath);
  if (!stream) throw std::runtime_error("cannot open resolved selection");
  SelectionJson j; stream >> j;
  if (j.at("schema") != "PhotonJetSelectionV1" || j.at("version") != 1)
    throw std::runtime_error("unsupported resolved selection");
  Profile p = nominal(system);
  const auto& s = j.at("systems").at(system == System::pp ? "pp" : "auau");
  auto line = [](const SelectionJson& v) {
    const double a = v.at("intercept").get<double>(), b = v.at("slope").get<double>();
    if (!std::isfinite(a) || !std::isfinite(b)) throw std::runtime_error("non-finite selection line");
    return Line{a,b};
  };
  p.configured = true;
  p.augmentedClasses = false;
  p.name = j.at("profile").get<std::string>();
  p.selectionId = j.at("selection_sha256").get<std::string>();
  if (p.selectionId.size()!=64 || p.name != "selection:"+p.selectionId)
    throw std::runtime_error("resolved selection identity differs");
  p.selectionJson = j.dump();
  p.vertexAbsZMax = s.at("event_selection").at("vertex_abs_z_max_cm").get<double>();
  if (!std::isfinite(p.vertexAbsZMax) || p.vertexAbsZMax <= 0)
    throw std::runtime_error("invalid event vertex window");
  p.tight = line(s.at("tight"));
  p.nonTightLow = line(s.at("non_tight_low"));
  p.nonTightHigh = line(s.at("non_tight_high"));
  p.scoreBranch = s.at("score_branch").get<std::string>();
  p.optionalID = s.contains("optional_model");
  p.nativeID=p.optionalID && s.at("optional_model").contains("native_evaluation");
  if (p.optionalID) {
    const auto &domain=s.at("optional_model").at("domain");
    p.photonEtMin=domain.at("et_min"); p.photonEtMax=domain.at("et_max_exclusive");
    p.photonAbsEtaMax=domain.at("eta_abs_max");
  }
  /* The nominal state remains part of every eventTree, but is never used to
     infer optional applicability. Optional scores encode that with NaN. */
  p.stateBranch = p.optionalID ? "pho_bdt_state" : s.at("state_branch").get<std::string>();
  if (p.optionalID) for (const auto& cell : s.at("optional_model").at("cells")) {
    const auto& et = cell.at("et"); const auto& cent = cell.at("centrality");
    p.idCells.push_back(IDCell{et.at(0), et.at(1), cent.at(0), cent.at(1),
                             line(cell.at("tight")), line(cell.at("non_tight_high"))});
  }
  if (p.optionalID && p.idCells.empty()) throw std::runtime_error("optional photon-ID has no routing cells");
  p.verifyNominalID = s.at("verify_nominal_id").get<bool>();
  p.verifyNominalIsolation = s.at("verify_nominal_isolation").get<bool>();
  const auto& iso = s.at("isolation");
  const double radius = iso.at("radius").get<double>();
  if (radius != 0.3 && radius != 0.4) throw std::runtime_error("isolation radius is not a retained witness");
  p.isolationRadiusIndex = radius == 0.3 ? 0 : 1;
  p.isolationSubtraction = iso.at("subtraction_method").get<int>();
  p.isolationFromWitness = iso.at("threshold").is_string() && iso.at("threshold") == "stored";
  if (!p.isolationFromWitness) p.isolationThreshold = line(iso.at("threshold"));
  p.sidebandGap = iso.at("sideband_gap_gev").is_null() ? std::numeric_limits<double>::quiet_NaN()
                 : iso.at("sideband_gap_gev").get<double>();
  if (!std::isnan(p.sidebandGap) && (!std::isfinite(p.sidebandGap) || p.sidebandGap < 0))
    throw std::runtime_error("invalid isolation sideband gap");
  p.isolationCentMin = iso.at("centrality_min").get<double>();
  p.isolationCentMax = iso.at("centrality_max").get<double>();
  for (double x : (system==System::pp ? std::vector<double>{15,35} : std::vector<double>{0,80}))
    if (!(p.nonTightLow.at(x) < p.nonTightHigh.at(x) && p.nonTightHigh.at(x) <= p.tight.at(x)))
      throw std::runtime_error("selection bands overlap or are empty");
  return p;
}

inline Profile validateSelectionInput(TFile& file, Profile p) {
  if (!p.configured) return p;
  const auto resolved = SelectionJson::parse(p.selectionJson);
  const auto& s = resolved.at("systems").at(p.system==System::pp ? "pp" : "auau");
  const std::string key = s.at("score_metadata").get<std::string>();
  const auto id = p.nativeID ? SelectionJson::object() : SelectionJson::parse(metaText(file,key));
  if (p.optionalID && !p.nativeID) {
    const auto& requested = s.at("optional_model");
    const auto& receipt = id.at(requested.at("model_name").get<std::string>());
    if (receipt.value("status", std::string()) != "PASS" ||
        receipt.at("release_sha256") != requested.at("release_sha256") ||
        receipt.at("score_branch") != requested.at("score_branch") || receipt.at("domain") != requested.at("domain"))
      throw std::runtime_error("optional score receipt differs from selected release/domain");
  } else if (!p.nativeID) {
    if (id.at("model_sha256") != s.at("model_sha256") ||
        id.value("freeze_sha256", std::string()) != s.at("photon_id_package_sha256").get<std::string>())
      throw std::runtime_error("saved score differs from selected feature/model package");
    if (p.scoreBranch != "pho_bdt_score" || p.stateBranch != "pho_bdt_state")
      if (id.value("score_branch", std::string()) != p.scoreBranch || id.value("state_branch", std::string()) != p.stateBranch)
        throw std::runtime_error("augmentation receipt does not bind requested score/state branches");
  }
  if (p.verifyNominalIsolation && p.system==System::auau) {
    const auto iso = SelectionJson::parse(metaText(file,"meta/isolation"));
    if (iso.at("freeze_sha256") != s.at("isolation").at("package_sha256"))
      throw std::runtime_error("saved nominal isolation definition differs from the configured package");
  }
  /* Older pp deliveries have no filled isolation class; their witness remains
     sufficient to evaluate PPG12. New pp deliveries additionally prove label parity. */
  if (p.system==System::pp && !file.Get("meta/isolation")) p.verifyNominalIsolation=false;
  return p;
}

/* The score has the same photon-aligned interface whether it is retained in
   eventTree or evaluated from a sealed ROOT model. H70 float32 inputs are
   already stored by TreeProduction; this reader adds no feature algorithm. */
class PhotonScores {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
  using BDT=TMVA::Experimental::RBDT;
#else
  using BDT=TMVA::Experimental::RBDT<>;
#endif
  struct Cell { double etLow,etHigh,centLow,centHigh; std::unique_ptr<BDT> model; };
  TTreeReader &reader_;
  std::unique_ptr<TTreeReaderArray<double>> saved_,et_,eta_;
  std::unique_ptr<TTreeReaderArray<float>> features_;
  std::unique_ptr<TTreeReaderArray<int>> finite_;
  std::unique_ptr<TTreeReaderValue<double>> vertex_,centrality_;
  std::unique_ptr<TTreeReaderValue<int>> record_,ready_,vertexValid_,terminal_;
  std::vector<Cell> cells_;
  SelectionJson domain_;
  size_t view_=0,stride_=0;
  bool auau_=false;
  mutable Long64_t entry_=-2;
  mutable std::vector<double> values_;

  void evaluate() const {
    if (entry_==reader_.GetCurrentEntry()) return;
    entry_=reader_.GetCurrentEntry(); const size_t n=et_->GetSize();
    if (eta_->GetSize()!=n || features_->GetSize()!=n*7*stride_ || finite_->GetSize()!=n*7)
      throw std::runtime_error("native H70 photon/features alignment differs");
    const double z=**vertex_,cent=**centrality_;
    if (auau_ && std::isfinite(cent) && cent!=std::round(cent)) throw std::runtime_error("native model centrality is not an integer bin");
    values_.assign(n,std::numeric_limits<double>::quiet_NaN());
    const bool event=**record_==0 && **ready_==1 && **vertexValid_==1 && **terminal_==0;
    std::vector<float> input(stride_);
    for (size_t i=0;i<n;++i) {
      const double et=(*et_)[i],eta=(*eta_)[i];
      if (!event || !(et>=domain_.at("et_min").get<double>() && et<domain_.at("et_max_exclusive").get<double>() &&
          std::abs(eta)<domain_.at("eta_abs_max").get<double>() && std::abs(z)<domain_.at("vertex_abs_max_cm").get<double>()) ||
          (auau_ && !(cent>=domain_.at("centrality_min").get<double>() && cent<domain_.at("centrality_max_exclusive").get<double>()))) continue;
      if ((*finite_)[7*i+view_]!=1) throw std::runtime_error("native model has invalid features inside its domain");
      for (size_t j=0;j<stride_;++j) {
        input[j]=(*features_)[(7*i+view_)*stride_+j];
        if (!std::isfinite(input[j])) throw std::runtime_error("nonfinite native model feature");
      }
      const Cell *selected=nullptr;
      for (const auto &cell : cells_) if (et>=cell.etLow && et<cell.etHigh && (!auau_ || (cent>=cell.centLow && cent<cell.centHigh))) {
        if (selected) throw std::runtime_error("native model cells overlap"); selected=&cell;
      }
      if (!selected) throw std::runtime_error("native model leaves an in-domain routing gap");
      const auto prediction=selected->model->Compute(input);
      if (prediction.size()!=1 || !std::isfinite(prediction[0]) || prediction[0]<0 || prediction[0]>1)
        throw std::runtime_error("native model returned an invalid score");
      values_[i]=prediction[0];
    }
  }
 public:
  PhotonScores(TFile &file,TTreeReader &reader,const Profile &profile) : reader_(reader) {
    if (!profile.nativeID) {saved_=std::make_unique<TTreeReaderArray<double>>(reader,profile.scoreBranch.c_str());return;}
    const auto selection=SelectionJson::parse(profile.selectionJson);
    const auto &native=selection.at("systems").at(profile.system==System::pp ? "pp" : "auau").at("optional_model").at("native_evaluation");
    et_=std::make_unique<TTreeReaderArray<double>>(reader,"pho_et");
    if (reader.GetTree()->GetEntries()==0) return;  // no photon feature dictionary is needed for an empty source
    const auto &contract=native.at("feature_contract"); domain_=native.at("domain"); auau_=profile.system==System::auau;
    const auto constants=SelectionJson::parse(metaText(file,"meta/constants"));
    const auto schema=SelectionJson::parse(metaText(file,"meta/schema"));
    const auto names=constants.at("showerViews.definition_name").get<std::vector<std::string>>();
    const std::string view=contract.at("view");
    if (names.size()!=7 || std::count(names.begin(),names.end(),view)!=1) throw std::runtime_error("native model shower-view dictionary differs");
    view_=std::find(names.begin(),names.end(),view)-names.begin();stride_=contract.at("count");
    if (schema.at("strides").at("shv_features")!=stride_ || contract.at("order").size()!=stride_ ||
        constants.at("showerViews.semantic_sha256").at(view_)!=contract.at("view_semantic_sha256"))
      throw std::runtime_error("native model feature count or semantics differ");
    const uint64_t bits=constants.at("showerViews.floor_gev").at(view_).at("f64").get<uint64_t>();double floor=0;
    std::memcpy(&floor,&bits,sizeof floor);
    if (floor!=contract.at("floor_gev").get<double>()) throw std::runtime_error("native model H70 tower floor differs");
    eta_=std::make_unique<TTreeReaderArray<double>>(reader,"pho_eta");
    features_=std::make_unique<TTreeReaderArray<float>>(reader,"shv_features");finite_=std::make_unique<TTreeReaderArray<int>>(reader,"shv_finite_feature_state");
    vertex_=std::make_unique<TTreeReaderValue<double>>(reader,"vtx_z");centrality_=std::make_unique<TTreeReaderValue<double>>(reader,"cent_value");
    record_=std::make_unique<TTreeReaderValue<int>>(reader,"evt_record");ready_=std::make_unique<TTreeReaderValue<int>>(reader,"evt_nodes_ready");
    vertexValid_=std::make_unique<TTreeReaderValue<int>>(reader,"vtx_valid");terminal_=std::make_unique<TTreeReaderValue<int>>(reader,"evt_terminal_status");
    for (const auto &cell : native.at("cells")) cells_.push_back({cell.at("et")[0],cell.at("et")[1],cell.at("centrality")[0],cell.at("centrality")[1],
      std::make_unique<BDT>(cell.at("key").get<std::string>(),cell.at("file").get<std::string>())});
    if (cells_.empty()) throw std::runtime_error("native model has no cells");
  }
  size_t GetSize() const {return saved_ ? saved_->GetSize() : et_->GetSize();}
  double operator[](size_t i) const {if (saved_) return (*saved_)[i];evaluate();return values_.at(i);}
};

template <class Values, class Ints, class Scores>
void validateSelectionEvent(const Profile& p, double centrality, const Values& et, const Scores& score, const Ints& state,
                            const Values& cone, const Values& threshold, const Values& sideband, const Ints& method,
                            const Ints& role, const Ints& idClass, const Ints& isoClass) {
  if (!p.configured) return;
  const std::size_t n=et.GetSize();
  if (score.GetSize()!=n || (!p.optionalID && state.GetSize()!=n) || cone.GetSize()!=2*n || threshold.GetSize()!=2*n ||
      sideband.GetSize()!=2*n || method.GetSize()!=2*n || role.GetSize()!=2*n ||
      (p.verifyNominalID && idClass.GetSize()!=n) || (p.verifyNominalIsolation && isoClass.GetSize()!=n))
    throw std::runtime_error("selection branches are not aligned with photons");
  for (std::size_t i=0;i<n;++i) {
    if (p.optionalID && (!std::isnan(score[i]) && (!std::isfinite(score[i]) || score[i] < 0. || score[i] > 1.)))
      throw std::runtime_error("optional score must be finite in [0,1] or NaN outside its domain");
    if (!p.optionalID && (state[i]<0 || state[i]>2 || (state[i]==1 && (!std::isfinite(score[i]) ||
        (p.system==System::auau && !std::isfinite(centrality))))))
      throw std::runtime_error("invalid score applicability/finite-value witness");
    const auto c=photonClass(p,i,et[i],centrality,score,state,cone,threshold,sideband,method,role,idClass,isoClass);
    if (p.verifyNominalID && idClass[i]!=c.idCode)
      throw std::runtime_error("saved nominal photon-ID class disagrees with its pinned score/definition");
    if (p.verifyNominalIsolation && isoClass[i]!=c.isoCode)
      throw std::runtime_error("saved nominal isolation class disagrees with its pinned witness/definition");
  }
}

inline void writeSelection(const Profile& p) {
  if (p.configured) TObjString(p.selectionJson.c_str()).Write("selection");
}

/* The public job binds one qualified AuAu fit to its exact selection and
   centrality weights. pp already contains its vertex factor in w_event. */
inline std::unique_ptr<TH1> auauVertexFit;
inline void configureVertexFit(const SelectionJson& spec) {
  auauVertexFit.reset();
  if(spec.is_null() || spec.empty()) return;  /* historical macro/reference mode */
  if(spec.at("coordinate")!="truth_vtx_z") throw std::runtime_error("unsupported vertex-fit coordinate");
  TDirectory::TContext context;
  std::unique_ptr<TFile> file(TFile::Open(spec.at("path").get<std::string>().c_str(),"READ"));
  if(!file || file->IsZombie()) throw std::runtime_error("cannot open AuAu vertex fit");
  auto* h=dynamic_cast<TH1*>(file->Get(spec.at("histogram").get<std::string>().c_str()));
  if(!h || h->GetDimension()!=1 || h->GetNbinsX()<2) throw std::runtime_error("invalid AuAu vertex histogram");
  for(int i=1;i<=h->GetNbinsX();++i)
    if(!std::isfinite(h->GetBinContent(i)) || h->GetBinContent(i)<0.)
      throw std::runtime_error("invalid AuAu vertex-fit weight");
  auauVertexFit.reset(static_cast<TH1*>(h->Clone("auau_vertex_fit")));auauVertexFit->SetDirectory(nullptr);
}

/* Stage-3 trees carry diagnostic copies of the configured AuAu factors.
   Their hashes determine whether to check the copies or ignore them. */
inline std::string auauReceiptSha, auauStitchingSha;
inline bool declaredAuAuCopies(const SelectionJson& w) {
  auto sha=[](const std::string& x) {return x.size()==64 && x.find_first_not_of("0123456789abcdef")==std::string::npos;};
  return w.value("composition", "")=="slice*vertex*centrality" && sha(w.value("receipt_sha256", "")) && sha(w.value("stitching_sha256", ""));
}
inline void checkAuAuCopies(double slice,double vertex,double centrality,double event,
                           double owned,double expectedVertex,double expectedCentrality) {
  auto same=[](double a,double b) {return std::isfinite(a) && std::isfinite(b) && (a==b || std::abs(a-b)<=1e-12*std::abs(b));};
  if (slice!=owned || !same(vertex,expectedVertex) || centrality!=expectedCentrality || !same(event,owned*expectedVertex*expectedCentrality))
    throw std::runtime_error("tree AuAu weight differs from the configured weight");
}
class AuAuTreeWeightCopies {
  bool auau_=false;
  CentralityWeights weights_;
  std::unique_ptr<OwnershipGate> ownership_;
  std::unique_ptr<TTreeReaderValue<int>> record_,ready_,terminal_,sampleValid_,vertexValid_;
  std::unique_ptr<TTreeReaderValue<double>> z_,cent_,slice_,vertex_,centrality_,event_,xsec_;
 public:
  AuAuTreeWeightCopies(TFile& file,TTreeReader& reader,const CentralityWeights& weights) : weights_(weights) {
    const bool simulation=fileSystem(file)==System::auau && isSimulation(file);
    const SelectionJson w=file.Get("meta/weighting") ? SelectionJson::parse(metaText(file,"meta/weighting")) : SelectionJson::object();
    const bool declared=simulation && declaredAuAuCopies(w);
    auau_=simulation && (strictReduction || declared);
    if (!declared) return;
    auto *tree=file.Get<TTree>("eventTree");
    for (const char *branch : {"w_xsec","w_sample_valid","w_vertex_valid"})
      if (!tree || !tree->GetBranch(branch)) throw std::runtime_error(std::string("declared AuAu copies lack ")+branch);
    if (!declaredAuAuCopies(w) || w.at("receipt_sha256")!=auauReceiptSha || w.at("stitching_sha256")!=auauStitchingSha) return;
    if (!auauVertexFit) throw std::runtime_error("matching AuAu copies need the configured truth-vertex fit");
    ownership_=std::make_unique<OwnershipGate>(reader);
    z_=std::make_unique<TTreeReaderValue<double>>(reader,"truth_vtx_z");
    cent_=std::make_unique<TTreeReaderValue<double>>(reader,"cent_value");
    slice_=std::make_unique<TTreeReaderValue<double>>(reader,"w_slice");
    vertex_=std::make_unique<TTreeReaderValue<double>>(reader,"w_vertex");
    centrality_=std::make_unique<TTreeReaderValue<double>>(reader,"w_centrality");
    event_=std::make_unique<TTreeReaderValue<double>>(reader,"w_event");
    xsec_=std::make_unique<TTreeReaderValue<double>>(reader,"w_xsec");
    record_=std::make_unique<TTreeReaderValue<int>>(reader,"evt_record");
    ready_=std::make_unique<TTreeReaderValue<int>>(reader,"evt_nodes_ready");
    terminal_=std::make_unique<TTreeReaderValue<int>>(reader,"evt_terminal_status");
    sampleValid_=std::make_unique<TTreeReaderValue<int>>(reader,"w_sample_valid");
    vertexValid_=std::make_unique<TTreeReaderValue<int>>(reader,"w_vertex_valid");
  }
  void check() const {
    if (!z_) return;
    const bool weighted=**record_==0 && **ready_==1 && **terminal_==0 && ownership_->owned();
    if (**sampleValid_!=1 || **vertexValid_!=(weighted ? 1 : 0) || **xsec_!=**slice_)
      throw std::runtime_error("tree AuAu weight state differs from the configured weight");
    if (!weighted) {
      checkAuAuCopies(**slice_,**vertex_,**centrality_,**event_,0.,1.,1.);
      return;
    }
    const double z=**z_;
    if (!std::isfinite(z) || z<=-9000.) throw std::runtime_error("matching AuAu tree weight lacks a truth vertex");
    const double v=(z<auauVertexFit->GetXaxis()->GetXmin() || z>auauVertexFit->GetXaxis()->GetXmax()) ? 0. : auauVertexFit->Interpolate(z);
    checkAuAuCopies(**slice_,**vertex_,**centrality_,**event_,1.,v,weights_(**cent_));
    if (xsec_ && **xsec_!=**slice_) throw std::runtime_error("tree AuAu w_xsec differs from w_slice");
  }
  double baseWeight(double existing) const {return auau_ ? 1. : existing;}
};

class VertexEventWeight {
  bool copied_=false;
  std::unique_ptr<TTreeReaderValue<double>> z_, previous_;
 public:
  VertexEventWeight(TFile& file,TTreeReader& reader) {
    if(!auauVertexFit) return;
    if(fileSystem(file)!=System::auau || !isSimulation(file))
      throw std::runtime_error("additional vertex fit applies only to AuAu SIM");
    const SelectionJson copies=file.Get("meta/weighting") ? SelectionJson::parse(metaText(file,"meta/weighting")) : SelectionJson::object();
    copied_=declaredAuAuCopies(copies);
    if(copies.value("vertex_enabled",false) && !copied_)
      throw std::runtime_error("input already includes a vertex weight");
    z_=std::make_unique<TTreeReaderValue<double>>(reader,"truth_vtx_z");
    previous_=std::make_unique<TTreeReaderValue<double>>(reader,"w_vertex");
  }
  double operator()() const {
    if(!z_) return 1.;
    if(!copied_ && **previous_!=1.) throw std::runtime_error("input vertex factor is not unity; refusing double weighting");
    const double z=**z_;
    if(!std::isfinite(z) || z<=-9000.) throw std::runtime_error("missing truth vertex for nominal AuAu weight");
    if(z<auauVertexFit->GetXaxis()->GetXmin() || z>auauVertexFit->GetXaxis()->GetXmax()) return 0.;
    const double w=auauVertexFit->Interpolate(z);
    if(!std::isfinite(w) || w<0.) throw std::runtime_error("invalid interpolated AuAu vertex weight");
    return w;
  }
};

/* The p+p simulation crossing-angle periods (0 and 1.5 mrad): the PPG12 two-period sum. PPG12, and our PPG12 parity
   production (the97, 2026-07-09), run the simulation once per period, each pass with that period's luminosity share,
   SI/DI fraction and truth-vertex fit, and add the passes. TreeProduction writes each event's factors for the one period
   its source is catalogued in; with the periods configured (selection pp.simulation_periods, mode all, read from
   TreeProduction's ppg12_weights) TreeToHists adds the passes per event from the same trees:
     w = w_slice x sum_p share_p x mix_p x v_p(z),  mix_p = 1 - f_DI,p (single) or f_DI,p (double interaction),
   v_p the period's fit at the hard-scatter truth vertex, times at the minimum-bias truth vertex for a double
   interaction. The term of the source's own period must equal its stored w_period x w_sidi x w_vertex, which proves the
   fits and constants are the producer's. */
struct PpPeriod {
  std::string name;
  double share = 0, doubleFraction = 0;
  std::unique_ptr<TH1> fit;
};
inline std::vector<PpPeriod> ppPeriods;  /* empty: each source's own period, w_event as written */
inline void configurePpPeriods(const SelectionJson& spec) {
  ppPeriods.clear();
  if (spec.is_null() || spec.empty() || spec.value("mode", std::string()) != "all") return;
  for (const auto& [name, p] : spec.at("periods").items()) {
    PpPeriod period;
    period.name = name;
    period.share = p.at("share").get<double>();
    period.doubleFraction = p.at("double_fraction").get<double>();
    if (!(period.share > 0 && period.share <= 1) || !(period.doubleFraction >= 0 && period.doubleFraction <= 1))
      throw std::runtime_error("invalid p+p period share or double-interaction fraction for " + name);
    TDirectory::TContext context;
    const auto& fit = p.at("fit");
    std::unique_ptr<TFile> file(TFile::Open(fit.at("path").get<std::string>().c_str(), "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("cannot open the " + name + " truth-vertex fit");
    auto* h = dynamic_cast<TH1*>(file->Get(fit.at("histogram").get<std::string>().c_str()));
    if (!h || h->GetDimension() != 1 || h->GetNbinsX() < 2) throw std::runtime_error("invalid " + name + " truth-vertex histogram");
    period.fit.reset(static_cast<TH1*>(h->Clone(("pp_vertex_fit_" + name).c_str())));
    period.fit->SetDirectory(nullptr);
    ppPeriods.push_back(std::move(period));
  }
}

/* The producer's truth-vertex lookup: interpolated, 0 outside the fit or for a non-positive value, NaN without a
   vertex (TreeProduction Event.cc, as PPG12). */
inline double truthVertexWeight(TH1& fit, double z) {
  if (!std::isfinite(z) || z <= -9000.) return std::numeric_limits<double>::quiet_NaN();
  if (z < fit.GetXaxis()->GetXmin() || z > fit.GetXaxis()->GetXmax()) return 0.;
  const double w = fit.Interpolate(z);
  return std::isfinite(w) && w > 0. ? w : 0.;
}

/* The lane component code of a pp simulation source. The eventTree branch ppg12_weight_lane_code carries the
   producer's ppg12_weight_lane_component_code (eventTree_v1: events.ppg12_weight_lane_component_code), not its raw
   lane: 1 photon+jet SI, 2 photon+jet DI, 3 inclusive-jet SI, 4 inclusive-jet DI, and 0 for an event whose canonical
   weight is inactive. The family follows from the meta/source dataset (isSim or isSimInclusive, checked against its
   lane when one is declared) and the SI/DI component from its si_di_role. */
inline int ppWeightComponent(const SelectionJson& source) {
  const std::string dataset = source.at("dataset").get<std::string>();
  if (dataset != "isSim" && dataset != "isSimInclusive")
    throw std::runtime_error("pp simulation dataset is neither isSim nor isSimInclusive: " + dataset);
  const bool inclusive = dataset == "isSimInclusive";
  const std::string role = source.value("si_di_role", std::string());
  if (role != "SI" && role != "DI")
    throw std::runtime_error("pp simulation source does not declare its SI/DI component (meta/source si_di_role)");
  // Older DI metadata names the SIM profile; the explicit role and the
  // numeric event component still identify the double-interaction weight.
  const std::string simLane = inclusive ? "pp_inclusive_sim" : "pp_photon_sim";
  const std::string diLane = inclusive ? "pp_inclusive_di" : "pp_photon_di";
  if (source.contains("lane") && source.at("lane") != simLane && !(role == "DI" && source.at("lane") == diLane))
    throw std::runtime_error("pp simulation source lane disagrees with its dataset/role " + dataset + "/" + role);
  return (inclusive ? 3 : 1) + (role == "DI" ? 1 : 0);
}

/* Canonical pp simulation selection.

   A structurally retained event can have been rejected by source stitching (ppg12_weight_active_state 0, component
   code 0). Only active canonical weights enter measurement, response and selected QA; an active event must carry its
   source's component code (ppWeightComponent) and one weight application. The saved slice, truth-vertex, SI/DI and
   period factors are checked against w_event and used once; weight() gives the event weight to fill with (w_event, or
   with the periods configured the PPG12 two-period sum above). Historical direct macros remain a separate reference
   interface. */
class PpSimulationGate {
  int expectedLane_ = 0, ownPeriod_ = -1;
  bool doubleInteraction_ = false;
  std::unique_ptr<TTreeReaderValue<int>> active_, lane_, applications_;
  std::unique_ptr<TTreeReaderValue<double>> event_, slice_, vertex_, mix_, period_, hardZ_, minimumBiasZ_;

 public:
  PpSimulationGate(TFile& file, TTreeReader& reader) {
    if (!strictReduction || fileSystem(file) != System::pp || !isSimulation(file)) return;
    const SelectionJson source = SelectionJson::parse(metaText(file, "meta/source"));
    expectedLane_ = ppWeightComponent(source);
    doubleInteraction_ = expectedLane_ == 2 || expectedLane_ == 4;
    if (!ppPeriods.empty()) {
      const std::string own = source.value("period", std::string());
      for (size_t i = 0; i < ppPeriods.size(); ++i) if (ppPeriods[i].name == own) ownPeriod_ = static_cast<int>(i);
      if (ownPeriod_ < 0) throw std::runtime_error("pp simulation source period '" + own + "' is not a configured period");
      hardZ_ = std::make_unique<TTreeReaderValue<double>>(reader, "truth_vtx_z");
      if (doubleInteraction_) minimumBiasZ_ = std::make_unique<TTreeReaderValue<double>>(reader, "truth_mbvtx_z");
    }
    active_ = std::make_unique<TTreeReaderValue<int>>(reader, "ppg12_weight_active_state");
    lane_ = std::make_unique<TTreeReaderValue<int>>(reader, "ppg12_weight_lane_code");
    applications_ = std::make_unique<TTreeReaderValue<int>>(reader, "w_application_count");
    event_ = std::make_unique<TTreeReaderValue<double>>(reader, "w_event");
    slice_ = std::make_unique<TTreeReaderValue<double>>(reader, "w_slice");
    vertex_ = std::make_unique<TTreeReaderValue<double>>(reader, "w_vertex");
    mix_ = std::make_unique<TTreeReaderValue<double>>(reader, "w_sidi");
    period_ = std::make_unique<TTreeReaderValue<double>>(reader, "w_period");
  }

  bool accepts() const {
    if (!active_) return true;
    if (**active_ == 0) {
      if (**lane_ != 0) throw std::runtime_error("inactive canonical pp simulation weight carries a lane component code");
      return false;
    }
    if (**active_ != 1 || **lane_ != expectedLane_ || **applications_ != 1)
      throw std::runtime_error("invalid canonical pp simulation weight state/lane component/application count");
    double product = 1.;
    for (double value : {**slice_, **vertex_, **mix_, **period_}) {
      if (!std::isfinite(value) || value < 0.)
        throw std::runtime_error("invalid canonical pp simulation weight factor");
      product *= value;
    }
    if (!std::isfinite(product) || !std::isfinite(**event_) ||
        std::abs(**event_ - product) > 1e-12 * std::max(1., std::abs(product)))
      throw std::runtime_error("pp event weight differs from its canonical factor product");
    return true;
  }

  /* The weight to fill an accepted event with: w_event as written, or w_slice x the sum over the configured periods. */
  double weight(double wEvent) const {
    if (!active_ || ppPeriods.empty() || **active_ != 1) return wEvent;
    double sum = 0., own = std::numeric_limits<double>::quiet_NaN();
    for (size_t i = 0; i < ppPeriods.size(); ++i) {
      const PpPeriod& p = ppPeriods[i];
      double vertex = truthVertexWeight(*p.fit, **hardZ_);
      if (doubleInteraction_) vertex *= truthVertexWeight(*p.fit, **minimumBiasZ_);
      const double term = p.share * (doubleInteraction_ ? p.doubleFraction : 1. - p.doubleFraction) * vertex;
      if (static_cast<int>(i) == ownPeriod_) own = term;
      sum += term;
    }
    const double stored = **period_ * **mix_ * **vertex_;
    if (!std::isfinite(sum) || !std::isfinite(own) || std::abs(own - stored) > 1e-12 * std::max(1., std::abs(stored)))
      throw std::runtime_error("the source period's term differs from the producer's w_period x w_sidi x w_vertex");
    return **slice_ * sum;
  }
};

/* Exact trigger names and enabled prescales are resolved per source, never by
   a fixed bit index. Absence/disablement is a valid empty stream; missing or
   ambiguous metadata is an error. The mask is applied to scaled event bits. */
inline ULong64_t dataTriggerMask(TFile& file, const std::string& key) {
  const std::string name = key=="photon4" ? "Photon 4 GeV + MBD NS >= 1" :
                          key=="photon12" ? "Photon 12 GeV + MBD NS >= 2, vtx < 150 cm" : "";
  if(name.empty()) throw std::runtime_error("unknown nominal DATA trigger: "+key);
  const auto table = SelectionJson::parse(metaText(file,file.Get("meta/trigger_selection") ?
                                                   "meta/trigger_selection" : "meta/trigger_run_info"));
  const auto& names=table.at("name");const auto& prescales=table.at("initial_prescale");
  if(!names.is_array()||names.size()!=64||prescales.size()!=64)
    throw std::runtime_error("trigger table must contain 64 names/prescales");
  int bit=-1;
  for(int i=0;i<64;++i) if(names[i]==name) {
    if(bit!=-1) throw std::runtime_error("ambiguous exact trigger name");bit=i;
  }
  return bit<0||prescales[bit].get<int>()<0 ? 0 : ULong64_t(1)<<bit;
}

/* DATA admission is shared by measurement and selected-object QA. Raw trigger
   QA remains inclusive. Run windows are half-open; a declared source period
   restricts the file to that period, otherwise the two supported periods apply. */
struct DataEventSelection {
  bool active=false, minimumBias=false;
  ULong64_t mask=0;
  std::vector<std::pair<int,int>> runWindows;
  DataEventSelection(TFile& file, const Profile& p) {
    active=p.configured && !isSimulation(file);
    if(!active) return;
    const auto definition=SelectionJson::parse(p.selectionJson);
    const auto& s=definition.at("systems").at(p.system==System::pp?"pp":"auau");
    const auto& e=s.at("event_selection");
    minimumBias=e.at("require_minimum_bias").get<bool>();
    mask=dataTriggerMask(file,s.at("data_trigger").get<std::string>());
    const auto source=SelectionJson::parse(metaText(file,"meta/source"));
    const auto period=source.value("period",std::string());
    const auto& periods=e.at("data_run_periods");
    if(p.system==System::pp && !period.empty() && !periods.contains(period))
      throw std::runtime_error("unrecognized pp source run period: "+period);
    for(auto it=periods.begin();it!=periods.end();++it)
      if(period.empty() || it.key()==period) runWindows.emplace_back(it.value().at(0),it.value().at(1));
  }
  bool accepts(int run, ULong64_t scaled, int mb) const {
    if(!active) return true;
    if(!(scaled&mask) || (minimumBias && mb!=1)) return false;
    if(runWindows.empty()) return true;
    return std::any_of(runWindows.begin(),runWindows.end(),[&](auto w){return run>=w.first && run<w.second;});
  }
};
}  /* namespace photonjet */

/* === Truth ===========================================================================================================
   Truth side of the photon + recoil-jet response on eventTree v1 simulation: the truth photon, its truth recoil
   jets, the reco-truth jet match, the nearby-truth-jet test of combinatoric reco jets, and the response categories. */

/* Companion of the selection above, which holds every reco decision. The truth recoil cuts are the reco recoil
   cuts of the same Profile. Functions take plain values, or indexable arrays (TTreeReaderArray, std::vector) with
   their length, so the ROOT reader (FillPhotonJetResponse.C) and any other caller make identical decisions. */

/* Order: eventTree collections keep the row order of the production tables, and a row identity is (domain,
   row + 1), so within an event "identity order" is "index order". The reference breaks ties by identity; the
   functions below break them by index, which is the same decision. */
namespace photonjet {

/* Truth photon acceptance of the reference response: |eta| below this. The reference keeps it apart from the reco
   photon window (Profile::photonAbsEtaMax); both are 0.7. */
constexpr double kTruthPhotonAbsEtaMax = 0.7;
/* Analysis jet radius (collections jet04 and truthjet04). The nearby-truth-jet test uses it as its Delta R. */
constexpr double kAnalysisJetRadius = 0.4;
/* truthpho_source_role of an analysis-signal truth photon. */
constexpr int kSignalSourceRole = 1;
/* The analysis truth signal, recorded by the producer as source role 1: a valid Geant4 photon of a prompt class below
   hadronic, |eta| < kTruthPhotonAbsEtaMax, and truth isolation in R = 0.3 below 4 GeV. */
constexpr int kHadronicPromptClass = 3;
constexpr double kTruthIsolationMaxGeV = 4.0;
inline bool analysisTruthSignal(int g4Valid, int promptClass, double eta, double isolationR03) {
  return g4Valid == 1 && promptClass < kHadronicPromptClass && std::fabs(eta) < kTruthPhotonAbsEtaMax &&
         isolationR03 < kTruthIsolationMaxGeV;
}
/* plink_class / jlink_class of a final reco-truth match. */
constexpr int kLinkMatch = 0;
/* jlink_class of a candidate edge: a reco and a truth jet of the same radius within the producer's match distance. */
constexpr int kLinkMatchCandidate = 5;

/* An analysis event enters the response when its vertex is finite and inside the profile window. This applies to
   every event, with or without a reco photon, before the truth side is filled. */
inline bool responseEvent(const Profile& p, double vz) { return std::isfinite(vz) && std::fabs(vz) < p.vertexAbsZMax; }

/* The truth photon of the event: the highest-pT truth photon with source role 1, finite kinematics, pT > 0 and
   |eta| < kTruthPhotonAbsEtaMax. Returns its index, or -1 when there is none. */
template <class Values, class Ints>
long truthPhoton(const Values& pt, const Values& eta, const Values& phi, const Ints& sourceRole, std::size_t n) {
  long best = -1;
  for (std::size_t t = 0; t < n; ++t) {
    if (sourceRole[t] != kSignalSourceRole) continue;
    if (!(std::isfinite(pt[t]) && std::isfinite(eta[t]) && std::isfinite(phi[t]))) continue;
    if (!(pt[t] > 0.0 && std::fabs(eta[t]) < kTruthPhotonAbsEtaMax)) continue;
    /* an exact pT tie goes to the later photon: the reference takes the maximum of (pT, identity) */
    if (best < 0 || pt[t] >= pt[best]) best = static_cast<long>(t);
  }
  return best;
}

/* A truth recoil jet of the truth photon: the reco recoil cuts of the profile (recoilJet) applied to a truth jet.
   Both angles take the jet first, as the reference does; the wrapped Delta phi is not bitwise symmetric in its
   arguments, so the order is kept. */
inline bool truthRecoilJet(const Profile& p, double photonEta, double photonPhi, double pt, double eta, double phi) {
  if (!(std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi))) return false;
  if (pt < p.jetPtMin) return false;
  if (deltaR(eta, phi, photonEta, photonPhi) < p.overlapDeltaR) return false;
  if (std::fabs(eta) >= p.jetAbsEtaMax) return false;
  return wrappedAbsDeltaPhi(phi, photonPhi) >= p.deltaPhiMin;
}

/* The indices i < n with accept(i), by descending pt[i]; equal pT keep index order. This is the reference's order of
   recoil jets: (-pT, identity) for truth jets, and (-pT, deterministic_order, identity) for reco jets, whose
   deterministic_order is the index within the collection. */
template <class Accept, class Values>
std::vector<long> orderedByPt(std::size_t n, Accept accept, const Values& pt) {
  std::vector<long> out;
  for (std::size_t i = 0; i < n; ++i)
    if (accept(i)) out.push_back(static_cast<long>(i));
  std::stable_sort(out.begin(), out.end(), [&pt](long a, long b) { return pt[a] > pt[b]; });
  return out;
}

/* A jet candidate edge: event-local reco and truth jet indices, its Delta R and the reco jet's corrected pT. */
struct JetEdge {
  std::size_t reco, truth;
  double deltaR, recoPt;
};

/* The producer's one-to-one jet match over a set of candidate edges: edges in order of Delta R, then reco pT
   descending, reco index, truth index (index order is identity order); an edge is taken when neither of its jets is
   taken yet. Returns the (reco, truth) pairs taken, sorted. */
inline std::vector<std::pair<std::size_t, std::size_t>> greedyJetMatch(std::vector<JetEdge> edges) {
  std::sort(edges.begin(), edges.end(), [](const JetEdge& a, const JetEdge& b) {
    if (a.deltaR != b.deltaR) return a.deltaR < b.deltaR;
    if (a.recoPt != b.recoPt) return a.recoPt > b.recoPt;
    if (a.reco != b.reco) return a.reco < b.reco;
    return a.truth < b.truth;
  });
  std::vector<std::pair<std::size_t, std::size_t>> taken;
  auto free = [&taken](std::size_t reco, std::size_t truth) {
    for (const auto& pair : taken)
      if (pair.first == reco || pair.second == truth) return false;
    return true;
  };
  for (const JetEdge& e : edges)
    if (free(e.reco, e.truth)) taken.emplace_back(e.reco, e.truth);
  std::sort(taken.begin(), taken.end());
  return taken;
}

/* True when a truth jet (truthjet04) above the threshold lies within matchDeltaR of the reco jet. A selected reco
   recoil of a matched photon without such a jet is combinatoric. Fiducial and recoil cuts are not applied to the
   truth jet: they belong to the observable, not to the label. */
template <class Values>
bool nearbyTruthJet(double threshold, double recoEta, double recoPhi, const Values& pt, const Values& eta,
                    const Values& phi, std::size_t n, double matchDeltaR = kAnalysisJetRadius) {
  for (std::size_t t = 0; t < n; ++t) {
    if (!(std::isfinite(pt[t]) && std::isfinite(eta[t]) && std::isfinite(phi[t]))) continue;
    if (pt[t] >= threshold && deltaR(eta[t], phi[t], recoEta, recoPhi) < matchDeltaR) return true;
  }
  return false;
}

/* Photon-level response category. A matched pair (the region-A photon's final match is the truth photon) fills the
   response when both photons are inside their pT axes. Otherwise the reco photon is a fake and the truth photon a
   miss, each only when inside its axis; so a matched photon whose truth pT is outside the truth axis is a reco fake. */
struct PhotonResponse {
  bool matched, fake, miss;
};

inline PhotonResponse photonResponse(bool linked, bool recoInAxis, bool truthInAxis) {
  if (linked && recoInAxis && truthInAxis) return {true, false, false};
  return {false, recoInAxis, truthInAxis};
}

/* Category of one selected reco recoil jet (inside the xJ axis) of the region-A photon:
   kMatched       photon matched, and the jet's final match is a selected truth recoil not yet taken
   kFakePhoton    photon not matched
   kCombinatoric  photon matched, no truth jet nearby (nearbyTruthJet)
   kDetectorFake  photon matched, a truth jet nearby, but no free linked truth recoil
   Reco recoils are classified in orderedByPt order and each truth recoil is taken by the first jet linked to it; a
   later jet linked to the same truth recoil is a fake. The result therefore depends on that order. truthJetNearby
   is only called for a matched photon without a match, as in the reference. */
enum class RecoilClass { kMatched, kFakePhoton, kCombinatoric, kDetectorFake };

template <class NearbyTest>
RecoilClass recoilClass(bool photonMatched, bool linkedToFreeTruthRecoil, NearbyTest truthJetNearby) {
  if (photonMatched && linkedToFreeTruthRecoil) return RecoilClass::kMatched;
  if (!photonMatched) return RecoilClass::kFakePhoton;
  return truthJetNearby() ? RecoilClass::kDetectorFake : RecoilClass::kCombinatoric;
}

/* Signal leakage: one region per signal truth photon. The leader of region r (A, B, C, D) is a witness of the truth
   photon of its final match when that photon has source role 1 (any truth photon of the event; no acceptance cut).
   A truth photon keeps its smallest witness by (match metric, -ET, region, photon index). Returns which regions'
   leaders are kept. photonTruth / photonMetric: final match of each reco photon (-1 / NaN when none).
   Throws when a witness has no finite match metric, as the reference does. */
template <class Values, class Ints>
std::array<bool, 4> leakageWitnesses(const std::array<long, 4>& leader, const std::vector<long>& photonTruth,
                                     const std::vector<double>& photonMetric, const Values& et,
                                     const Ints& truthSourceRole) {
  using Order = std::tuple<double, double, int, long>;
  std::array<long, 4> truth{-1, -1, -1, -1};
  std::array<Order, 4> order{};
  for (int r = 0; r < 4; ++r) {
    const long i = leader[r];
    if (i < 0 || photonTruth[i] < 0 || truthSourceRole[photonTruth[i]] != kSignalSourceRole) continue;
    if (!std::isfinite(photonMetric[i])) throw std::runtime_error("final photon match without a finite metric");
    truth[r] = photonTruth[i];
    order[r] = Order(photonMetric[i], -et[i], r, i);
  }
  std::array<bool, 4> keep{false, false, false, false};
  for (int r = 0; r < 4; ++r) {
    if (truth[r] < 0) continue;
    keep[r] = true;
    for (int s = 0; s < 4; ++s)
      if (s != r && truth[s] == truth[r] && order[s] < order[r]) keep[r] = false;
  }
  return keep;
}

}  /* namespace photonjet */
