// Fill the photon-selection and event-property histograms of the analysis note (the figures that are not the
// photon + recoil-jet measurement) from eventTree v1 files, data or simulation, in one pass.
//
//   root -l -q 'TreeToHists/src/FillPhotonQA.C("a.eventtree.root b.eventtree.root", "out.root", "nominal", "<weights>")'
//
// Every input shares one kind: DATA, photon-family simulation (isSim, isSimEmbedded) or inclusive-family simulation
// (isSimInclusive, isSimEmbeddedInclusive), written as "family" (data, phosim, jetsim). Photon candidates are labelled
// by that kind and their final reco-truth match (plink_class 0, the match FillPhotonJetResponse.C uses):
//   data        every candidate of a DATA file
//   sig         photon-family simulation, matched to an analysis-signal truth photon (source role 1: Geant4-valid,
//               prompt, |eta| < 0.7, truth isolation R = 0.3 below 4 GeV; PhotonJetTruth.h), the signal the
//               measurement unfolds and the ABCD leakage is normalised to
//   phoother    photon-family simulation, any other candidate
//   bkg         inclusive-family simulation, not matched to an analysis-signal truth photon
//   jetprompt   inclusive-family simulation, matched to an analysis-signal truth photon (kept apart)
// Selection stages of a candidate (PhotonJetSelection.h profile; analysis events, |vz| inside the profile window,
// finite kinematics, |eta| < 0.7):
//   0 every candidate   1 active preselection   2 and tight photon ID   3 and isolated (Region A)
// Tight, non-tight and isolated are the profile's: Au+Au from the augmented classes, pp from the score surfaces and
// the isolation witness, exactly as regionLeaders decides them. The isolation value is the profile's witness (Au+Au
// R = 0.3 SUB1, pp R = 0.4 topocluster). "No UE subtraction" (Au+Au) is that witness plus the energy the SUB1
// subtraction removed from its cone constituents, sum over isoc_ of (calibrated - sub1) / cosh(eta_photon + deta).
//
// Histograms (TH*F; Sumw2 in simulation). <c>: centrality class (Au+Au cent0_20, cent20_50, cent50_80; pp all);
// <f>: candidate label; <F>: file family; slices: the 16 5% slices of 0-80% (one slice in pp); ET axis 10-40 GeV in
// 1 GeV bins unless stated:
//   ss_<var>_<f>_<c>        var x ET x stage      shower shapes and the photon-ID score: e11_over_e33, e32_over_e35,
//                                                 bdt, weta_cogx, wphi_cogx, weta33_cogx, wphi33_cogx, et1..et4
//                                                 (et1..et4: the NPB-input values)
//   ncb_<f>_<c>             NPB score x ET        stage 0 (the NPB score enters the preselection)
//   iso_<f>_<c>             iso x ET x ID         stage 1; ID 0 tight, 1 non-tight, 2 neither
//   isonoue_<f>_<c>         the same without UE subtraction (Au+Au)
//   abcd_<f>_<c>            region x ET           the event-leading photon of each region (regionLeaders, as in
//                                                 FillPhotonJet.C) by its label; ET 15-35 GeV in 1 GeV bins
//   recoil_<f>_<c>          region (A,C) x ET x xJ  every recoil jet (5 GeV) of the A and C leaders by the leader's
//                                                 label, on FillPhotonJet.C's axes: summed over labels it is
//                                                 FillPhotonJet's recoil_<c>_jet5 (the xJ-level purity closure)
//   abcdown_<f>_<c>         region x ET x class   simulation: abcd by the event's ownership class, the family's
//   recoilown_<f>_<c>       (region x class) x ET x xJ   generator variable (leading prompt truth photon pT for the
//                                                 photon family, leading truth jet pT for the inclusive family) on
//                                                 edges holding every sample window, so a stitched mixture keeps
//                                                 each sample to its owned window (region-major flattened axis)
//   scoreiso_<f>_<c>        score x iso x ET      stage 1 with a score and a witness: 20 score bins, iso 1 GeV bins
//                                                 over -10..40, ET 15-35 in 2 GeV (factorisation tests in sub-bands)
//   clu_et_<f>_<c>          ET (0-50 GeV) x stage   (every stage axis is cumulative: a candidate fills each stage
//   clu_etaphi_<f>_<c>      eta x phi x stage       it reached); ET >= 10 GeV for eta x phi
//   bdt_<f>                 score x ET x slice    stage 1, no isolation requirement
//   isoslice_<f>            iso x ET x slice      stage 1, any ID
//   isonoueslice_<f>        the same without UE subtraction (Au+Au)
//   width_<var>_<f>         var x ET x slice      stage 1; weta_cogx, wphi_cogx, weta33_cogx, wphi33_cogx in 0.01 bins
//                                                 over 0-1 (the shape histograms keep the full range)
//   e11e33_dt_data_<c>      e11/e33 x ET x dt     DATA, stage 0; dt = cluster time - MBD t0 [ns]
//   etatime_data_<c>        eta x dt              DATA, stage 0, ET >= 10 GeV
//   evt_vz_<F>_<c>          vertex z              events before the vertex window
//   jet_pt_<F>_<c>          jet pT x selection    analysis-view R = 0.4 jets; selection 0 every analysis event,
//   jet_etaphi_<F>_<c>      eta x phi x selection 1 events with a Region-A leader (15 <= ET < 35 GeV); pT >= 5 GeV
//   evt_cent_<F>            centrality x weight   Au+Au, inside the vertex window; weight 0 w_event, 1 the full weight
//   cent_bits_data          centrality x scaled bit   Au+Au DATA inside the vertex window, minimum-bias classified
//                                                 (cent_mb_decision == 1): one entry per fired scaled bit. The data target
//                                                 of the simulation centrality weights for a trigger chosen by name from
//                                                 the run's own table (centrality.py; the schema-10 package used Photon 10)
//   evt_calo_<F>            total calo E x centrality x MB state   state = (cent_mb_decision == 1) +
//                                                 2 (evt_embedded_mb_decision == 1); every analysis event
//   trig_<F>                max cluster E x scaled bit x live bit  DATA, every recorded event: for each scaled bit s,
//                                                 once at live bit 64 (the reference count) and once per live bit set
//   truth_iso_<c>           truth iso x truth pT x prompt class    photon-family simulation: Geant4-valid truth photons
//                                                 with |eta| < 0.7 and a valid truth isolation (R = 0.3)
//   eff_<c>                 truth pT x stage x reco ET   photon-family simulation, analysis-signal truth photons (source
//                                                 role 1): stage 0 every one (reco ET bin [-2.5, 0)), 1 matched to a
//                                                 candidate, 2 and preselected,
//                                                 3 and isolated, 4 and tight, 5 and tight and isolated
//   gen_<F>                 pT x kind             simulation, every recorded event: kind 0 the leading prompt truth
//                                                 photon, 1 the leading truth jet (R = 0.4)
// Analysis events: evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && evt_terminal_status == 0 (recorded
// events: the first two and the last). Weights are w_event (1 when not finite) times, in Au+Au simulation, the
// centrality weight of the event's 5% slice (an input, as in FillPhotonJetResponse.C; 1 when the product is not
// finite). Written with "profile", "family", "centrality_weights", "sources" and "qa_contract" (version and axes), and for
// DATA "trigger_run_info" (the inputs' meta/trigger_run_info: bit names, prescales, run totals) with
// "trigger_run_info_consistent" (1 when every input carried the same table: one run per output).
#include "SelectionInput.h"
#include "PhotonJetTruth.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libTreePlayer)
#endif

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TObjString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace photonjet {
namespace qa {

constexpr int kQaVersion = 1;
constexpr double kQaEtMin = 10.0, kQaEtMax = 40.0, kStoredAbsEtaMax = 0.7;
const std::vector<std::string> kShapes = {"e11_over_e33", "e32_over_e35", "bdt", "weta_cogx", "wphi_cogx",
                                          "weta33_cogx", "wphi33_cogx", "et1", "et2", "et3", "et4"};
const std::vector<std::string> kWidths = {"weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx"};

struct Range {
  int bins;
  double low, high;
};

inline Range shapeRange(const std::string& var) {
  if (var == "e11_over_e33" || var == "bdt") return {100, 0.0, 1.0};
  if (var == "e32_over_e35" || var == "et1") return {105, 0.0, 1.05};
  if (var == "et2" || var == "et3") return {80, -0.5, 1.1};
  if (var == "et4") return {100, 0.0, 0.5};
  if (var == "weta_cogx" || var == "wphi_cogx") return {100, 0.0, 2.5};
  return {100, 0.0, 1.5};  // 3x3 widths
}

enum class Family { data, phosim, jetsim };

inline const char* familyName(Family f) { return f == Family::data ? "data" : f == Family::phosim ? "phosim" : "jetsim"; }

inline Family fileFamily(TFile& f) {
  if (!isSimulation(f)) return Family::data;
  const std::string source = metaText(f, "meta/source");
  for (const char* dataset : {"isSimInclusive", "isSimEmbeddedInclusive"})
    if (source.find(std::string("\"dataset\":\"") + dataset + "\"") != std::string::npos) return Family::jetsim;
  return Family::phosim;
}

inline std::vector<std::string> labels(Family f) {
  if (f == Family::data) return {"data"};
  if (f == Family::phosim) return {"sig", "phoother"};
  return {"bkg", "jetprompt"};
}

class Book {
 public:
  explicit Book(bool sumw2) : sumw2_(sumw2) {}
  TH1F& h1(const std::string& n, const char* t, Range x) { return get(h1_, n, [&] { return new TH1F(n.c_str(), t, x.bins, x.low, x.high); }); }
  TH2F& h2(const std::string& n, const char* t, Range x, Range y) {
    return get(h2_, n, [&] { return new TH2F(n.c_str(), t, x.bins, x.low, x.high, y.bins, y.low, y.high); });
  }
  TH3F& h3(const std::string& n, const char* t, Range x, Range y, Range z) {
    return get(h3_, n, [&] { return new TH3F(n.c_str(), t, x.bins, x.low, x.high, y.bins, y.low, y.high, z.bins, z.low, z.high); });
  }
  TH3F& h3v(const std::string& n, const char* t, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z) {
    return get(h3_, n, [&] {
      return new TH3F(n.c_str(), t, x.size() - 1, x.data(), y.size() - 1, y.data(), z.size() - 1, z.data());
    });
  }
  std::vector<std::string> names() const {
    std::vector<std::string> out;
    for (const auto& [n, h] : h1_) out.push_back(n);
    for (const auto& [n, h] : h2_) out.push_back(n);
    for (const auto& [n, h] : h3_) out.push_back(n);
    return out;
  }
  void write() {
    for (auto& [n, h] : h1_) h->Write();
    for (auto& [n, h] : h2_) h->Write();
    for (auto& [n, h] : h3_) h->Write();
  }

 private:
  template <class Map, class Make>
  typename Map::mapped_type::element_type& get(Map& map, const std::string& name, Make make) {
    auto it = map.find(name);
    if (it != map.end()) return *it->second;
    auto* h = make();
    h->SetDirectory(nullptr);
    if (sumw2_) h->Sumw2();
    map[name].reset(h);
    return *h;
  }
  bool sumw2_;
  std::map<std::string, std::unique_ptr<TH1F>> h1_;
  std::map<std::string, std::unique_ptr<TH2F>> h2_;
  std::map<std::string, std::unique_ptr<TH3F>> h3_;
};

const std::vector<double> kRecoilRegions = {0.0, 1.0, 2.0}, kRegions4 = {0.0, 1.0, 2.0, 3.0, 4.0};
// ownership classes of the generator variable: every sample window of config/samples.yaml and PPG18's photon cuts
const std::vector<double> kOwnershipEdges = {-1.0, 0.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 21.0, 24.0, 25.0, 30.0, 31.0, 35.0, 40.0,
                                             41.0, 50.0, 60.0, 1.0e4};
inline std::vector<double> recoilOwnershipAxis() {
  std::vector<double> e;
  for (std::size_t k = 0; k <= 2 * (kOwnershipEdges.size() - 1); ++k) e.push_back(static_cast<double>(k));
  return e;
}
const std::vector<double> kRecoilOwnership = recoilOwnershipAxis();
const Range kEt{30, kQaEtMin, kQaEtMax}, kStage{4, 0.0, 4.0}, kIso{500, -10.0, 40.0}, kId{3, 0.0, 3.0}, kSlice16{16, 0.0, 16.0},
    kEta{58, -kStoredAbsEtaMax, kStoredAbsEtaMax}, kPhi{256, -M_PI, M_PI}, kDt{160, -40.0, 40.0}, kDtCoarse{40, -20.0, 20.0}, kWidth{100, 0.0, 1.0},
    kJetEta{44, -1.1, 1.1}, kJetPhi{64, -M_PI, M_PI}, kSel{2, 0.0, 2.0};

// Public packages always contain the same objects, including empty selections.
inline void prebook(Book& b, Family family, System system, const std::vector<CentralityClass>& classes) {
  const bool auau=system==System::auau, sim=family!=Family::data;
  const std::string F=familyName(family);
  const Range slices=auau?kSlice16:Range{1,0,1};
  for (const auto& cls:classes) {
    const std::string c=cls.name;
    for (const auto& label:labels(family)) {
      const std::string t=label+"_"+c;
      for (const auto& v:kShapes) b.h3("ss_"+v+"_"+t,"",shapeRange(v),kEt,kStage);
      b.h2("ncb_"+t,"",{110,-1.1,1.1},kEt);
      b.h3("iso_"+t,"",kIso,kEt,kId);
      if (auau) b.h3("isonoue_"+t,"",kIso,kEt,kId);
      b.h2("abcd_"+t,"",{4,0,4},{20,15,35});
      b.h3v("recoil_"+t,"",kRecoilRegions,photonEtEdges(),xjEdges());
      b.h3("scoreiso_"+t,"",{20,0,1},{50,-10,40},{10,15,35});
      b.h2("clu_et_"+t,"",{100,0,50},kStage);
      b.h3("clu_etaphi_"+t,"",kEta,kPhi,kStage);
      if (sim) {
        b.h3v("abcdown_"+t,"",kRegions4,photonEtEdges(),kOwnershipEdges);
        b.h3v("recoilown_"+t,"",kRecoilOwnership,photonEtEdges(),xjEdges());
      }
    }
    b.h1("evt_vz_"+F+"_"+c,"",{240,-60,60});
    b.h2("jet_pt_"+F+"_"+c,"",{120,0,60},kSel);
    b.h3("jet_etaphi_"+F+"_"+c,"",kJetEta,kJetPhi,kSel);
    if (!sim) {
      b.h3("e11e33_dt_data_"+c,"",shapeRange("e11_over_e33"),kEt,kDtCoarse);
      b.h2("etatime_data_"+c,"",kEta,kDt);
    }
    if (family==Family::phosim) {
      b.h3("truth_iso_"+c,"",{300,0,30},{35,5,40},{4,0,4});
      b.h3("eff_"+c,"",{35,5,40},{6,0,6},{21,-2.5,50});
    }
  }
  for (const auto& label:labels(family)) {
    b.h3("bdt_"+label,"",{200,0,1},kEt,slices);
    b.h3("isoslice_"+label,"",kIso,kEt,slices);
    if (auau) b.h3("isonoueslice_"+label,"",kIso,kEt,slices);
    for (const auto& v:kWidths) b.h3("width_"+v+"_"+label,"",kWidth,kEt,slices);
  }
  b.h3("evt_calo_"+F,"",auau?Range{300,0,3000}:Range{200,0,200},auau?Range{100,0,100}:Range{1,0,100},{4,0,4});
  if (auau) {
    b.h2("evt_cent_"+F,"",{100,0,100},kSel);
    if (!sim) b.h2("cent_bits_data","",{100,0,100},{64,0,64});
  }
  if (sim) b.h2("gen_"+F,"",{160,0,80},{2,0,2});
  else b.h3("trig_data","",{100,0,40},{64,0,64},{65,0,65});
}

// The candidate's photon-ID and isolation decisions under the profile (the per-candidate form of regionLeaders).
struct Decision {
  bool scored = false, tight = false, nonTight = false, isolated = false, nonIsolated = false, witness = false;
  double iso = std::numeric_limits<double>::quiet_NaN();
};

template <class Values, class Ints>
Decision decide(const Profile& p, std::size_t i, double et, double centrality, const Values& score, const Ints& state,
                const Values& cone, const Values& threshold, const Values& sideband, const Ints& method, const Ints& role,
                const Ints& idClass, const Ints& isoClass) {
  // the selection's own classes (PhotonJetSelection.h photonClass), plus the witness cone sum for the isolation plots
  const PhotonClass c = photonClass(p, i, et, centrality, score, state, cone, threshold, sideband, method, role, idClass, isoClass);
  Decision d;
  d.scored = c.scored;
  d.witness = c.witness;
  d.tight = c.tight;
  d.nonTight = c.nonTight;
  d.isolated = c.isolated;
  d.nonIsolated = c.nonIsolated;
  if (c.witness) d.iso = cone[2 * i + p.isolationRadiusIndex];
  return d;
}

inline bool candidateInWindow(const Profile& p, double vz, double et, double eta, double phi) {
  if (!(std::isfinite(vz) && std::isfinite(et) && std::isfinite(eta) && std::isfinite(phi))) return false;
  return std::fabs(vz) < p.vertexAbsZMax && std::fabs(eta) < p.photonAbsEtaMax;
}

inline void fillFile(const std::string& path, const Profile& requested, Family family, const CentralityWeights& centralityWeight,
                     const std::vector<CentralityClass>& classes, Book& b) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  const Profile profile = validateSelectionInput(*file, requested);
  if (fileSystem(*file) != profile.system) throw std::runtime_error("profile system differs from " + path);
  if (fileFamily(*file) != family) throw std::runtime_error("inputs mix data and simulation families: " + path);
  const bool auau = profile.system == System::auau, sim = family != Family::data;
  const DataEventSelection dataSelection(*file,profile);
  const std::string F = familyName(family);
  TTreeReader reader("eventTree", file.get());
  const VertexEventWeight vertexWeight(*file,reader);
  TTreeReaderValue<int> record(reader, "evt_record"), ready(reader, "evt_nodes_ready"), vtxValid(reader, "vtx_valid"),
      terminal(reader, "evt_terminal_status"), centMb(reader, "cent_mb_decision"), embMb(reader, "evt_embedded_mb_decision"), run(reader,"evt_run");
  TTreeReaderValue<double> vz(reader, "vtx_z"), centrality(reader, "cent_value"), weightValue(reader, "w_event"),
      caloTotal(reader, "calo_total"), maxCluster(reader, "calo_max_cluster_e"), mbdT0(reader, "mbd_t0");
  TTreeReaderValue<ULong64_t> scaled(reader, "trig_scaled"), live(reader, "trig_live");
  TTreeReaderArray<double> et(reader, "pho_et"), eta(reader, "pho_eta"), phi(reader, "pho_phi"), score(reader, profile.scoreBranch.c_str()),
      time(reader, "pho_time_ns"), npb(reader, "pho_npb_score");
  TTreeReaderArray<int> presel(reader, "pho_presel_active"), state(reader, profile.stateBranch.c_str()), idClass(reader, "pho_id_class"),
      isoClass(reader, "pho_iso_class");
  std::vector<std::unique_ptr<TTreeReaderArray<double>>> shape;
  for (const std::string& var : kShapes) {
    std::string branch = "pho_" + var;
    if (var == "bdt") branch = profile.scoreBranch;
    if (var == "et1" || var == "et2" || var == "et3" || var == "et4") branch = "pho_npb_cluster_" + var;
    shape.emplace_back(new TTreeReaderArray<double>(reader, branch.c_str()));
  }
  TTreeReaderArray<double> cone(reader, "isow_cone_sum"), threshold(reader, "isow_threshold"),
      sideband(reader, "isow_sideband_threshold");
  TTreeReaderArray<int> method(reader, "isow_subtraction_method"), role(reader, "isow_reconstructed_or_truth");
  TTreeReaderArray<int> isocPho(reader, "isoc_pho");
  TTreeReaderArray<double> isocDr(reader, "isoc_delta_r"), isocDeta(reader, "isoc_delta_eta"), isocCal(reader, "isoc_calibrated_energy"),
      isocSub(reader, "isoc_sub1_energy");
  TTreeReaderArray<double> jetPt(reader, "jet04_pt"), jetEta(reader, "jet04_eta"), jetPhi(reader, "jet04_phi");
  // simulation only
  std::unique_ptr<TTreeReaderArray<double>> truthPt, truthEta, truthIso, truthJetPt;
  std::unique_ptr<TTreeReaderArray<int>> truthRole, truthValid, truthPrompt, truthIsoValid, linkReco, linkTruth, linkClass;
  if (sim) {
    truthPt.reset(new TTreeReaderArray<double>(reader, "truthpho_pt"));
    truthEta.reset(new TTreeReaderArray<double>(reader, "truthpho_eta"));
    truthIso.reset(new TTreeReaderArray<double>(reader, "truthpho_truth_isolation_r03"));
    truthJetPt.reset(new TTreeReaderArray<double>(reader, "truthjet04_pt"));
    truthRole.reset(new TTreeReaderArray<int>(reader, "truthpho_source_role"));
    truthValid.reset(new TTreeReaderArray<int>(reader, "truthpho_g4_photon_valid"));
    truthPrompt.reset(new TTreeReaderArray<int>(reader, "truthpho_prompt_class"));
    truthIsoValid.reset(new TTreeReaderArray<int>(reader, "truthpho_truth_isolation_valid"));
    linkReco.reset(new TTreeReaderArray<int>(reader, "plink_reco"));
    linkTruth.reset(new TTreeReaderArray<int>(reader, "plink_truth"));
    linkClass.reset(new TTreeReaderArray<int>(reader, "plink_class"));
  }
  const double isoRadius = profile.isolationRadiusIndex == 0 ? 0.3 : 0.4;
  const Range caloRange = auau ? Range{300, 0.0, 3000.0} : Range{200, 0.0, 200.0};
  const Range centRange = auau ? Range{100, 0.0, 100.0} : Range{1, 0.0, 100.0};

  while (reader.Next()) {
    reportProgress(reader.GetCurrentEntry()+1);
    validateSelectionEvent(profile, *centrality, et, score, state, cone, threshold, sideband, method, role, idClass, isoClass);
    if (*record != 0 || *ready != 1 || *terminal != 0) continue;
    if (strictReduction && !std::isfinite(*weightValue)) throw std::runtime_error("nonfinite QA event weight");
    const double wEvent = (std::isfinite(*weightValue) ? *weightValue : 1.0) * vertexWeight();
    double w = wEvent;
    if (sim && auau) {
      w = wEvent * centralityWeight(*centrality);
      if (strictReduction && !std::isfinite(w)) throw std::runtime_error("nonfinite QA centrality weight");
      if (!std::isfinite(w)) w = 1.0;
    }
    // every recorded event: trigger characterisation (DATA) and generator-level spectra (simulation)
    double ownership = -1.0;  // simulation: the family's generator variable (leading prompt truth photon or truth jet pT)
    if (!sim) {
      TH3F& trig = b.h3("trig_" + F, "recorded events;max cluster E [GeV];scaled bit;live bit (64: any)", {100, 0.0, 40.0},
                        {64, 0.0, 64.0}, {65, 0.0, 65.0});
      const double e = std::isfinite(*maxCluster) ? *maxCluster : -1.0;
      for (int s = 0; s < 64; ++s) {
        if (!((*scaled >> s) & 1ULL)) continue;
        trig.Fill(e, s + 0.5, 64.5);
        for (int l = 0; l < 64; ++l)
          if ((*live >> l) & 1ULL) trig.Fill(e, s + 0.5, l + 0.5);
      }
    } else {
      TH2F& gen = b.h2("gen_" + F, "recorded events;p_{T} [GeV];0 leading prompt truth photon, 1 leading truth jet",
                       {160, 0.0, 80.0}, {2, 0.0, 2.0});
      double photon = -1.0, jet = -1.0;
      for (std::size_t t = 0; t < truthPt->GetSize(); ++t)
        if ((*truthValid)[t] == 1 && (*truthPrompt)[t] < kHadronicPromptClass && std::isfinite((*truthPt)[t]))
          photon = std::max(photon, (*truthPt)[t]);
      for (std::size_t j = 0; j < truthJetPt->GetSize(); ++j)
        if (std::isfinite((*truthJetPt)[j])) jet = std::max(jet, (*truthJetPt)[j]);
      if (photon >= 0) gen.Fill(photon, 0.5, wEvent);
      if (jet >= 0) gen.Fill(jet, 1.5, wEvent);
      ownership = family == Family::phosim ? photon : jet;
    }
    if (*vtxValid != 1) continue;
    // the event's class and 5% slice
    int k = 0;
    if (auau) {
      while (k < static_cast<int>(classes.size()) && !(classes[k].low <= *centrality && *centrality < classes[k].high)) ++k;
    }
    const bool inClass = k < static_cast<int>(classes.size());
    const std::string c = inClass ? classes[k].name : "";
    const int slice = auau ? ((*centrality >= 0.0 && *centrality < 80.0) ? centralitySlice(*centrality) : -1) : 0;
    const bool inWindow = std::isfinite(*vz) && std::fabs(*vz) < profile.vertexAbsZMax;
    if (inClass) b.h1("evt_vz_" + F + "_" + c, "analysis events;v_{z} [cm]", {240, -60.0, 60.0}).Fill(*vz, w);
    {
      const int mb = ((*centMb == 1) ? 1 : 0) + ((*embMb == 1) ? 2 : 0);
      b.h3("evt_calo_" + F, "analysis events;total calorimeter E [GeV];centrality [%];MB state", caloRange, centRange,
           {4, 0.0, 4.0}).Fill(*caloTotal, *centrality, mb + 0.5, w);
    }
    if (!inWindow) continue;
    if (auau) {
      TH2F& cent = b.h2("evt_cent_" + F, "analysis events;centrality [%];0 w_event, 1 full weight", {100, 0.0, 100.0}, kSel);
      cent.Fill(*centrality, 0.5, wEvent);
      cent.Fill(*centrality, 1.5, w);
    }
    if (auau && !sim && *centMb == 1) {
      TH2F& bits = b.h2("cent_bits_" + F, "analysis events, minimum-bias classified;centrality [%];scaled bit", {100, 0.0, 100.0},
                        {64, 0.0, 64.0});
      for (int s = 0; s < 64; ++s)
        if ((*scaled >> s) & 1ULL) bits.Fill(*centrality, s + 0.5);
    }
    if (!inClass || !dataSelection.accepts(*run,*scaled,*centMb)) continue;
    const std::size_t n = et.GetSize();
    if (cone.GetSize() != 2 * n) throw std::runtime_error("isolation witnesses are not 2 per photon in " + path);
    if (profile.augmentedClasses && (idClass.GetSize() != n || isoClass.GetSize() != n))
      throw std::runtime_error("the augmented classes pho_id_class, pho_iso_class are not filled in " + path);
    // the UE energy SUB1 removed from each candidate's cone (Au+Au)
    std::vector<double> ue(n, 0.0);
    if (auau)
      for (std::size_t q = 0; q < isocPho.GetSize(); ++q) {
        const int i = isocPho[q];
        if (i < 0 || static_cast<std::size_t>(i) >= n || !(isocDr[q] < isoRadius)) continue;
        const double ch = std::cosh(eta[i] + isocDeta[q]);
        if (std::isfinite(isocCal[q]) && std::isfinite(isocSub[q]) && ch > 0) ue[i] += (isocCal[q] - isocSub[q]) / ch;
      }
    // the prompt truth photon matched to each candidate (simulation)
    std::vector<long> truthOf(n, -1);
    if (sim)
      for (std::size_t l = 0; l < linkReco->GetSize(); ++l) {
        if ((*linkClass)[l] != 0 || (*linkReco)[l] < 0 || (*linkTruth)[l] < 0) continue;
        const std::size_t r = (*linkReco)[l], t = (*linkTruth)[l];
        if (r >= n || t >= truthPt->GetSize()) throw std::runtime_error("photon link outside its collections in " + path);
        truthOf[r] = static_cast<long>(t);
      }
    // the candidate's label: the analysis truth signal (source role 1) decides sig / jetprompt
    auto labelOf = [&](std::size_t i) -> std::string {
      if (!sim) return "data";
      const long t = truthOf[i];
      const bool signal = t >= 0 && (*truthRole)[t] == kSignalSourceRole;
      return family == Family::phosim ? (signal ? "sig" : "phoother") : (signal ? "jetprompt" : "bkg");
    };
    std::vector<int> stageOf(n, -1);
    std::vector<Decision> decisionOf(n);
    long leaderA = -1;
    for (std::size_t i = 0; i < n; ++i) {
      if (!candidateInWindow(profile, *vz, et[i], eta[i], phi[i])) continue;
      const Decision d = decide(profile, i, et[i], *centrality, score, state, cone, threshold, sideband, method, role, idClass, isoClass);
      decisionOf[i] = d;
      int stage = 0;
      if (presel[i] == 1) {
        stage = 1;
        if (d.tight) {
          stage = 2;
          if (d.isolated) stage = 3;
        }
      }
      stageOf[i] = stage;
      if (stage == 3 && profile.photonEtMin <= et[i] && et[i] < profile.photonEtMax && (leaderA < 0 || et[i] >= et[leaderA]))
        leaderA = static_cast<long>(i);
      const std::string f = labelOf(i);
      const std::string tag = f + "_" + c;
      // stage axes are cumulative: a candidate fills every stage it reached
      TH2F& cluEt = b.h2("clu_et_" + tag, "candidates;E_{T} [GeV];stage (reached)", {100, 0.0, 50.0}, kStage);
      for (int s = 0; s <= stage; ++s) cluEt.Fill(et[i], s + 0.5, w);
      if (et[i] >= 10.0) {
        TH3F& cluEtaPhi = b.h3("clu_etaphi_" + tag, "candidates E_{T} >= 10 GeV;#eta;#phi;stage (reached)", kEta, kPhi, kStage);
        for (int s = 0; s <= stage; ++s) cluEtaPhi.Fill(eta[i], phi[i], s + 0.5, w);
      }
      if (!(kQaEtMin <= et[i] && et[i] < kQaEtMax)) continue;
      for (std::size_t v = 0; v < kShapes.size(); ++v) {
        const double x = (*shape[v])[i];
        if (!std::isfinite(x)) continue;
        TH3F& h = b.h3("ss_" + kShapes[v] + "_" + tag, (";" + kShapes[v] + ";E_{T} [GeV];stage").c_str(), shapeRange(kShapes[v]),
                       kEt, kStage);
        for (int s = 0; s <= stage; ++s) h.Fill(x, et[i], s + 0.5, w);
      }
      if (std::isfinite(npb[i])) b.h2("ncb_" + tag, ";NPB score;E_{T} [GeV]", {110, -1.1, 1.1}, kEt).Fill(npb[i], et[i], w);
      const double dt = time[i] - *mbdT0;
      if (!sim && std::isfinite(dt)) {
        const double e11 = (*shape[0])[i];
        if (std::isfinite(e11))
          b.h3("e11e33_dt_data_" + c, ";E_{11}/E_{33};E_{T} [GeV];t_{cluster} - t_{MBD} [ns]", shapeRange("e11_over_e33"), kEt,
               kDtCoarse).Fill(e11, et[i], dt, w);
        b.h2("etatime_data_" + c, "E_{T} >= 10 GeV;#eta;t_{cluster} - t_{MBD} [ns]", kEta, kDt).Fill(eta[i], dt, w);
      }
      if (stage < 1) continue;
      const int id = d.tight ? 0 : d.nonTight ? 1 : 2;
      if (d.witness && state[i] == 1 && std::isfinite(score[i]) && et[i] >= 15.0 && et[i] < 35.0)
        b.h3("scoreiso_" + tag, "preselected;photon-ID score;E_{T}^{iso} [GeV];E_{T} [GeV]", {20, 0.0, 1.0}, {50, -10.0, 40.0},
             {10, 15.0, 35.0}).Fill(score[i], d.iso, et[i], w);
      if (d.witness) {
        b.h3("iso_" + tag, ";E_{T}^{iso} [GeV];E_{T} [GeV];ID (tight, non-tight, neither)", kIso, kEt, kId).Fill(d.iso, et[i], id + 0.5, w);
        if (auau)
          b.h3("isonoue_" + tag, "no UE subtraction;E_{T}^{iso} [GeV];E_{T} [GeV];ID (tight, non-tight, neither)", kIso, kEt, kId)
              .Fill(d.iso + ue[i], et[i], id + 0.5, w);
      }
      if (slice < 0) continue;
      const Range slices = auau ? kSlice16 : Range{1, 0.0, 1.0};
      if (state[i] == 1 && std::isfinite(score[i]))
        b.h3("bdt_" + f, "preselected;photon-ID score;E_{T} [GeV];5% slice", {200, 0.0, 1.0}, kEt, slices).Fill(score[i], et[i], slice + 0.5, w);
      if (d.witness) {
        b.h3("isoslice_" + f, "preselected;E_{T}^{iso} [GeV];E_{T} [GeV];5% slice", kIso, kEt, slices).Fill(d.iso, et[i], slice + 0.5, w);
        if (auau)
          b.h3("isonoueslice_" + f, "preselected, no UE subtraction;E_{T}^{iso} [GeV];E_{T} [GeV];5% slice", kIso, kEt, slices)
              .Fill(d.iso + ue[i], et[i], slice + 0.5, w);
      }
      for (const std::string& var : kWidths) {
        std::size_t v = 0;
        while (kShapes[v] != var) ++v;
        const double x = (*shape[v])[i];
        if (std::isfinite(x))
          b.h3("width_" + var + "_" + f, ("preselected;" + var + ";E_{T} [GeV];5% slice").c_str(), kWidth, kEt, slices)
              .Fill(x, et[i], slice + 0.5, w);
      }
    }
    // the event-leading photon of each region, as FillPhotonJet.C selects it, by its label
    const std::array<long, 4> leader = regionLeaders(profile, *vz, *centrality, n, et, eta, phi, score, presel, state, cone, threshold,
                                                     sideband, method, role, idClass, isoClass);
    for (int r = 0; r < 4; ++r) {
      if (leader[r] < 0) continue;
      const long i = leader[r];
      const std::string f = labelOf(i);
      b.h2("abcd_" + f + "_" + c, "event-leading photons;region (A,B,C,D);E_{T} [GeV]", {4, 0.0, 4.0}, {20, 15.0, 35.0}).Fill(r + 0.5, et[i], w);
      if (sim)  // the same by the ownership class of the event (a sample's owned window of its generator variable)
        b.h3v("abcdown_" + f + "_" + c, "event-leading photons;region (A,B,C,D);E_{T} [GeV];generator p_{T} class", kRegions4,
              photonEtEdges(), kOwnershipEdges).Fill(r + 0.5, et[i], std::min(ownership, kOwnershipEdges.back() - 1.0), w);
      if (r != kA && r != kC) continue;
      // its recoil jets at the nominal threshold, exactly as FillPhotonJet.C fills recoil_<c>_jet5
      TH3F& recoil = b.h3v("recoil_" + f + "_" + c, "recoil jets (5 GeV);region (A,C);photon E_{T} [GeV];x_{J#gamma}", kRecoilRegions,
                          photonEtEdges(), xjEdges());
      Profile p = profile;
      p.jetPtMin = kJetThresholds[0];
      for (std::size_t j = 0; j < jetPt.GetSize(); ++j) {
        const double xj = jetPt[j] / et[i];
        if (!(xj >= 0.0 && xj < 3.0)) continue;
        if (!recoilJet(p, eta[i], phi[i], jetPt[j], jetEta[j], jetPhi[j])) continue;
        recoil.Fill(r == kA ? 0.5 : 1.5, et[i], xj, w);
        if (sim) {  // region (A, C) x ownership class flattened: bin (region) * n_class + class
          const int cls = static_cast<int>(std::upper_bound(kOwnershipEdges.begin(), kOwnershipEdges.end(),
                                                            std::min(ownership, kOwnershipEdges.back() - 1.0)) - kOwnershipEdges.begin()) - 1;
          const int nCls = static_cast<int>(kOwnershipEdges.size()) - 1;
          b.h3v("recoilown_" + f + "_" + c, "recoil jets (5 GeV);region (A,C) x generator p_{T} class;photon E_{T} [GeV];x_{J#gamma}",
                kRecoilOwnership, photonEtEdges(), xjEdges()).Fill((r == kA ? 0 : 1) * nCls + std::max(cls, 0) + 0.5, et[i], xj, w);
        }
      }
    }
    // jets of every analysis event, and of events with a Region-A leader
    for (std::size_t j = 0; j < jetPt.GetSize(); ++j) {
      if (!(std::isfinite(jetPt[j]) && std::isfinite(jetEta[j]) && std::isfinite(jetPhi[j]))) continue;
      for (int s = 0; s < (leaderA >= 0 ? 2 : 1); ++s) {
        b.h2("jet_pt_" + F + "_" + c, "R = 0.4 jets;p_{T} [GeV];0 every event, 1 with a Region-A leader", {120, 0.0, 60.0}, kSel)
            .Fill(jetPt[j], s + 0.5, w);
        if (jetPt[j] >= 5.0)
          b.h3("jet_etaphi_" + F + "_" + c, "R = 0.4 jets p_{T} >= 5 GeV;#eta;#phi;0 every event, 1 with a Region-A leader", kJetEta,
               kJetPhi, kSel).Fill(jetEta[j], jetPhi[j], s + 0.5, w);
      }
    }
    if (!sim) continue;
    // truth photons: isolation and prompt class (photon family), and the selection efficiency of the analysis signal
    for (std::size_t t = 0; t < truthPt->GetSize(); ++t) {
      if (!(std::isfinite((*truthPt)[t]) && std::isfinite((*truthEta)[t]))) continue;
      if (family == Family::phosim && (*truthValid)[t] == 1 && std::fabs((*truthEta)[t]) < kStoredAbsEtaMax &&
          (*truthIsoValid)[t] == 1 && std::isfinite((*truthIso)[t]))
        b.h3("truth_iso_" + c, "truth photons |#eta| < 0.7;E_{T}^{iso,truth} (R = 0.3) [GeV];p_{T}^{truth} [GeV];prompt class",
             {300, 0.0, 30.0}, {35, 5.0, 40.0}, {4, 0.0, 4.0}).Fill((*truthIso)[t], (*truthPt)[t], (*truthPrompt)[t] + 0.5, w);
      if (family != Family::phosim || (*truthRole)[t] != kSignalSourceRole) continue;
      // reco E_T axis: first bin [-2.5, 0) holds the truth photons of stage 0 (no reco candidate)
      TH3F& eff = b.h3("eff_" + c, "analysis-signal truth photons;p_{T}^{truth} [GeV];stage;matched reco E_{T} [GeV]", {35, 5.0, 40.0},
                       {6, 0.0, 6.0}, {21, -2.5, 50.0});
      eff.Fill((*truthPt)[t], 0.5, -1.0, w);
      long r = -1;
      for (std::size_t i = 0; i < n; ++i)
        if (truthOf[i] == static_cast<long>(t) && stageOf[i] >= 0) r = static_cast<long>(i);
      if (r < 0) continue;
      const Decision& d = decisionOf[r];
      const bool pre = presel[r] == 1;
      const std::array<bool, 5> reached{true, pre, pre && d.isolated, pre && d.tight, pre && d.tight && d.isolated};
      for (int s = 0; s < 5; ++s)
        if (reached[s]) eff.Fill((*truthPt)[t], s + 1.5, et[r], w);
    }
  }
  // a missing branch or a read error ends the loop like the end of the tree; only the end is accepted
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
}

}  // namespace qa
}  // namespace photonjet

// inputs: whitespace-separated eventTree paths of one family. profile: "nominal" or "historical". centralityWeights:
// Au+Au simulation only and required there, 16 weights for the 5% slices of 0-80% (CentralityWeights).
int FillPhotonQA(const char* inputs, const char* output, const char* profileName = "nominal", const char* centralityWeightText = "", const char* selection = "") {
  using namespace photonjet;
  using namespace photonjet::qa;
  const std::vector<std::string> paths = inputPaths(inputs);
  std::unique_ptr<TFile> first(TFile::Open(paths.front().c_str(), "READ"));
  if (!first || first->IsZombie()) throw std::runtime_error("cannot open " + paths.front());
  const System system = fileSystem(*first);
  const Family family = fileFamily(*first);
  first.reset();
  const bool needWeights = family != Family::data && system == System::auau;
  if (needWeights == (std::string(centralityWeightText).find_first_not_of(" \t") == std::string::npos))
    throw std::runtime_error(needWeights ? "Au+Au simulation needs its 16 centrality weights" : "centrality weights apply to Au+Au simulation only");
  const CentralityWeights weights = needWeights ? photonjet::centralityWeights(centralityWeightText, system) : CentralityWeights{};
  const Profile profile = configuredProfile(profileName, system, selection);
  const std::vector<CentralityClass> classes = centralityClasses(system);
  Book b(strictReduction || family != Family::data);
  if (strictReduction) prebook(b, family, system, classes);
  std::ostringstream sources;
  sources << "[";
  std::string triggerInfo;  // DATA: the run's trigger table (bit names, prescales, run totals), identical in every input
  bool triggerConsistent = true;
  for (std::size_t f = 0; f < paths.size(); ++f) {
    fillFile(paths[f], profile, family, weights, classes, b);
    std::unique_ptr<TFile> in(TFile::Open(paths[f].c_str(), "READ"));
    sources << (f ? "," : "") << "{\"path\":\"" << paths[f] << "\",\"input_file_sha256\":\"" << sourceSha(*in) << "\"}";
    const char* triggerKey = in->Get("meta/trigger_selection") ? "meta/trigger_selection" : "meta/trigger_run_info";
    if (family == Family::data && in->Get(triggerKey)) {
      const std::string info = metaText(*in, triggerKey);
      if (f == 0) triggerInfo = info;
      else if (info != triggerInfo) triggerConsistent = false;
    }
  }
  sources << "]";
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  b.write();
  std::ostringstream contract;
  contract << "{\"schema\":\"RJPhotonQAV1\",\"version\":" << kQaVersion << ",\"et\":[" << kQaEtMin << "," << kQaEtMax
           << "],\"stages\":[\"every candidate\",\"preselected\",\"tight\",\"tight and isolated\"],\"histograms\":" << b.names().size() << "}";
  TObjString(profile.name.c_str()).Write("profile");
  writeSelection(profile);
  TObjString(familyName(family)).Write("family");
  std::ostringstream w;  // the record FillPhotonJetResponse.C writes: {"per_slice":[...]} at full precision
  w.precision(17);
  w << "{\"per_slice\":[";
  for (std::size_t s = 0; s < weights.perSlice.size(); ++s) w << (s ? "," : "") << weights.perSlice[s];
  w << "]}";
  TObjString(w.str().c_str()).Write("centrality_weights");
  TObjString(sources.str().c_str()).Write("sources");
  TObjString(contract.str().c_str()).Write("qa_contract");
  if (family == Family::data) {
    TObjString(triggerInfo.c_str()).Write("trigger_run_info");
    TObjString(triggerConsistent ? "1" : "0").Write("trigger_run_info_consistent");
  }
  out->Close();
  return 0;
}
