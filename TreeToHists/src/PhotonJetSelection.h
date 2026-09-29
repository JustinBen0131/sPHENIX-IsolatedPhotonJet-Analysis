// Photon + recoil-jet selection on eventTree v1: photon ID and isolation regions, the event-leading
// photon of each region, and its recoil jets.
//
// One place for every selection decision. The functions take plain values so the ROOT reader
// (FillPhotonJet.C) and any other caller make identical decisions.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace photonjet {

enum class System { pp, auau };
enum Region { kA = 0, kB = 1, kC = 2, kD = 3, kNone = -1 };  // A tight+iso, B tight+non-iso, C non-tight+iso, D non-tight+non-iso

struct Line {
  double intercept, slope;
  double at(double x) const { return intercept + slope * x; }
};

// One half-open routing cell from a saved PhotonID selection definition.
struct IDCell {
  double etMin, etMax, centMin, centMax;
  Line tight, nonTightHigh;
  bool contains(double et, double cent) const {
    return etMin <= et && et < etMax && centMin <= cent && cent < centMax;
  }
};

// Every selection parameter. Profiles below are the only sanctioned settings.
struct Profile {
  std::string name;
  System system = System::pp;
  // photon
  double photonEtMin = 15.0, photonEtMax = 35.0, photonAbsEtaMax = 0.7;
  double vertexAbsZMax = 60.0;
  // photon ID on pho_bdt_score. pp: surfaces in photon ET. Au+Au: surfaces in centrality.
  Line tight{0.815625, -0.0015625};
  Line nonTightLow{0.7333333333333333, -0.01333333333333333};
  Line nonTightHigh{0.684375, 0.0015625};
  bool nonTightLowIsConstant = false;  // Au+Au nominal: 0.1
  double nonTightLowConstant = 0.1;
  // isolation witness: index 0 is R = 0.3, 1 is R = 0.4 (meta rules isolation_witnesses)
  int isolationRadiusIndex = 1;
  int isolationSubtraction = 1;  // producer isolation method: 1 pp topocluster, 2 Au+Au calorimeter SUB1
  // recoil jets (the analysis collection jet04 = anti-kT R = 0.4)
  double jetPtMin = 5.0, jetAbsEtaMax = 0.7, overlapDeltaR = 0.4, deltaPhiMin = 7.0 * M_PI / 8.0;
  // Regions from the augmented classes pho_id_class and pho_iso_class, computed once at conversion from the frozen
  // photon-ID and isolation packages, instead of the score surfaces and the stored witness thresholds.
  bool augmentedClasses = false;
  // ABCD sideband variations (profile variants, see profileVariant): NaN or 0 is the profile's own definition.
  double nonTightMin = std::numeric_limits<double>::quiet_NaN();  // Au+Au: the non-tight lower bound on the score
  double nonTightLowShift = 0.0;                                  // pp: shift of the non-tight lower score line
  double sidebandGap = std::numeric_limits<double>::quiet_NaN();  // GeV from the isolation threshold to the sideband
  // Resolved YAML selections evaluate the primitives, never historical producer labels.
  bool configured = false, isolationFromWitness = true, verifyNominalID = false, verifyNominalIsolation = false;
  Line isolationThreshold{0.0, 0.0};  // pp: ET, AuAu: centrality
  double isolationCentMin = 0.0, isolationCentMax = 80.0;
  std::string scoreBranch = "pho_bdt_score", stateBranch = "pho_bdt_state";
  std::string selectionJson, selectionId;
  bool optionalID = false;  // finite saved score is applicability; no extra state branch
  std::vector<IDCell> idCells;

};

// Persisted photon-ID and isolation class codes from the eventTree contract.
constexpr int kIdTight = 3, kIdNonTight = 1, kIdBelowNonTight = 0, kIsoIsolated = 2, kIsoGap = 1, kIsoNonIsolated = 0;

// The historical reference analysis (schema-10 package v3): R = 0.4 isolation with the stored per-photon
// thresholds, Au+Au non-tight WP90 < s < WP80, pp |z| < 30 cm. Used only to prove this code reproduces it.
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
    p.tight = {0.6927385742113304, 0.0010476808935335573};         // WP70
    p.nonTightHigh = {0.6071724560222583, 0.0013019695025559648};  // WP80
    p.nonTightLow = {0.4634577166843198, 0.0015466253667456042};   // WP90
  }
  return p;
}

// The nominal selection. pp is the validated PPG12-parity selection (NPB preselection, PPG12 BDT tight and
// non-tight bands, R = 0.4 isolation below 0.490 + 0.037 ET with the sideband 0.8 GeV above it, |z| < 60 cm):
// the historical selection with the 60 cm vertex window. Au+Au takes its regions from the augmented classes: tight
// s > WP70(c) and non-tight 0.1 < s < WP80(c) (pho_id_class), isolated below the R = 0.3 package line and non-isolated
// above it plus the gap (pho_iso_class); the other photon and event cuts are the historical ones.
inline Profile nominal(System system) {
  Profile p = historical(system);
  p.name = "nominal";
  if (system == System::pp) {
    p.vertexAbsZMax = 60.0;
  } else {
    p.augmentedClasses = true;
    p.isolationRadiusIndex = 0;  // R = 0.3, applied through pho_iso_class
  }
  return p;
}

inline double wrappedAbsDeltaPhi(double a, double b) {
  // |dphi| in [0, pi], the same arithmetic as the reference: fmod into [-pi, pi) then abs
  double d = std::fmod(a - b + M_PI, 2.0 * M_PI);
  if (d < 0) d += 2.0 * M_PI;
  return std::fabs(d - M_PI);
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return std::hypot(eta1 - eta2, wrappedAbsDeltaPhi(phi1, phi2));
}

// The photon-ID and isolation classes of photon i, which already passed the kinematic selection: the one definition of
// the regions, used for the region leaders and the QA stages. score is pho_bdt_score (finite when pho_bdt_state == 1);
// cone, threshold, sideband, method and role are the isolation witnesses, two per photon (R = 0.3, then R = 0.4).
// With augmented classes the codes decide, and a variation reads the stored score (non-tight lower bound) or joins the
// gap code to the non-isolated side (sideband gap 0), so no package edge is recomputed.
struct PhotonClass {
  bool scored = false;   // a photon-ID score, or with augmented classes always
  bool witness = false;  // a reconstructed witness of the profile's radius and method with a finite cone sum
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

template <class Values, class Ints>
PhotonClass photonClass(const Profile& p, std::size_t i, double et, double centrality, const Values& score, const Ints& state,
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
    // When a pp threshold is overridden, move its sideband with it; retaining
    // the old absolute upper edge would silently change the requested gap.
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
                           score[i] > p.nonTightMin;  // both codes lie below the non-tight upper bound
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

// The event-leading photon of each region (A, B, C, D), -1 when the region has none. A photon enters when it has
// active preselection, passes the kinematic selection, has a photon-ID score (pho_bdt_state 1) and a reconstructed
// isolation witness of the profile's radius and method; with augmented classes, when its classes give a region. The
// leader is the highest-E_T photon of its region; an exact tie goes to the later photon, as the reference's (hi, lo)
// identity order does. Photon arrays hold one value per photon; isolation arrays hold two per photon (R = 0.3, then
// R = 0.4). idClass and isoClass (pho_id_class, pho_iso_class) are read only with augmented classes. et may be of
// another indexable type than the other arrays (a varied energy scale).
template <class Energies, class Values, class Ints>
std::array<long, 4> regionLeaders(const Profile& p, double vz, double centrality, std::size_t n, const Energies& et,
                                  const Values& eta, const Values& phi, const Values& score, const Ints& presel,
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

inline bool recoilJet(const Profile& p, double photonEta, double photonPhi, double pt, double eta, double phi) {
  if (!(std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi))) return false;
  if (pt < p.jetPtMin) return false;
  if (deltaR(photonEta, photonPhi, eta, phi) < p.overlapDeltaR) return false;
  if (std::fabs(eta) >= p.jetAbsEtaMax) return false;
  return wrappedAbsDeltaPhi(phi, photonPhi) >= p.deltaPhiMin;
}

}  // namespace photonjet
