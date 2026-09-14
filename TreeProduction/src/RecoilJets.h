// Tell Emacs this is C++   -*- C++ -*-
#ifndef RECOILJETS_H
#define RECOILJETS_H

// ============================================================================
// RecoilJets
//   - Photon (tight/iso) selection + recoil jet matching + xJ/alpha/JES outputs
//   - Supports pp, Au+Au, and simulation
//   - IMPORTANT: jet-dependent outputs are radius-tagged (r02, r04, ...)
// ============================================================================

// ---------------------------------------------------------------------------
// Fun4All / PHOOL
// ---------------------------------------------------------------------------
#include <fun4all/SubsysReco.h>
#include "RJGl1TriggerContract.h"
#include <phool/PHCompositeNode.h>

// ---------------------------------------------------------------------------
// ROOT
// ---------------------------------------------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile3D.h>
#include <TLorentzVector.h>
#include <TVector2.h>

// ---------------------------------------------------------------------------
// sPHENIX object headers (used as pointer types throughout)
// NOTE: These are kept here for "drop-in" builds since the .cc relies on
//       this header to provide these definitions in multiple places.
// ---------------------------------------------------------------------------
#include <calobase/RawClusterContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <globalvertex/GlobalVertexMap.h>
#include <jetbase/JetContainer.h>

// You currently include PhotonClusterBuilder in the header; keep for drop-in.
// (This is not ideal for portability, but preserves your existing build.)
#include "PhotonClusterBuilder.h"

// ---------------------------------------------------------------------------
// STL
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Helper alias: maps histogram name -> ROOT object*
// ---------------------------------------------------------------------------
using HistMap = std::map<std::string, TObject*>;

// Forward declarations (pointers only in the interface)
class CentralityInfo;
class Gl1Packet;
class EventHeader;
class MbdOut;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHG4VtxPoint;

class GlobalVertex;
class RawCluster;
class PhotonClusterv1;
class Jet;
namespace RJReplayRuntimeV1 { class Runtime; }
namespace RJTriggerScalersV1 { class Writer; }
namespace RJTriggerRunInfoV1 { class Writer; }
namespace RJPhotonTrainingViewV1 { class Runtime; }
namespace RJReplayFoundationV1 { struct TruthPhotonInventoryV1; }

// g4eval: used for truth↔reco association of EMCal clusters
class CaloRawClusterEval;

// HepMC forward declarations (used by unified truth signal classifier)
namespace HepMC
{
  class GenEvent;
  class GenParticle;
}

// ============================================================================
// Photon ID working points / helper functions  (PPG12 Table 4)
// ============================================================================
// NOTE:
//   .cc uses:  using namespace PhoIDCuts;
//   PPG12 defines:
//     - Preselection cuts applied first (fail => rejected; does NOT enter A–B–C–D)
//     - Preselection variants:
//       reference = hard-cut SS preselection
//       newPPG12 = PPG12 common SS preselection plus NPB score, no weta cut
//       noPreCriteria = no preselection
//       onlyNPB = NPB score only
//       refPlusNPB = reference hard-cut SS preselection plus NPB score
//       auauOnlyNPB = AuAu-trained NPB score only (AuAu pipeline)
//     - Tight cuts applied only after preselection
//     - Non-tight: passes preselection AND fails >=2 of the 5 tight cuts
//   Also: for photons in this context, pT^gamma ≡ ET^gamma (use pt_gamma everywhere).
// ============================================================================
namespace PhoIDCuts
{
  // ------------------- Pre-selection (Table 4) -------------------
  // Applied to EMCal photon-cluster candidates BEFORE tight/non-tight classification.
  inline constexpr double PRE_E11E33_MAX = 0.98;  // E11/E33 < 0.98
  inline constexpr double PRE_ET1_MIN    = 0.60;  // 0.6 < et1 < 1.0
  inline constexpr double PRE_ET1_MAX    = 1.00;
  inline constexpr double PRE_E32E35_MIN = 0.80;  // 0.8 < E32/E35 < 1.0
  inline constexpr double PRE_E32E35_MAX = 1.00;
  inline constexpr double PRE_WETA_MAX   = 0.60;  // weta_cogx < 0.6  (only η width capped in preselection)

  // ------------------- Tight ID (Table 4) -------------------
  // Applied ONLY to candidates that PASS preselection.
  // Widths use a pT(=ET) dependent upper bound:
  inline constexpr double TIGHT_W_LO = 0.0;       // 0 < w < (0.15 + 0.006 * ET^gamma)

  inline constexpr double TIGHT_E11E33_MIN = 0.40; // 0.4 < E11/E33 < 0.98
  inline constexpr double TIGHT_E11E33_MAX = 0.98;

  inline constexpr double TIGHT_ET1_MIN    = 0.90; // 0.9 < et1 < 1.0
  inline constexpr double TIGHT_ET1_MAX    = 1.00;

  inline constexpr double TIGHT_E32E35_MIN = 0.92; // 0.92 < E32/E35 < 1.0
  inline constexpr double TIGHT_E32E35_MAX = 1.00;

  // Open interval helper: (lo < x < hi)
  inline bool in_open_interval(double x, double lo, double hi)
  {
    return (std::isfinite(x) && std::isfinite(lo) && std::isfinite(hi) && (x > lo) && (x < hi));
  }

  // Tight upper width bound: w_hi(ET^gamma) = 0.15 + 0.006 * ET^gamma
  // In your code, pass pt_gamma (since pT^gamma == ET^gamma here).
  inline double tight_w_hi(double pt_gamma)
  {
    if (!std::isfinite(pt_gamma) || pt_gamma <= 0.0)
      return std::numeric_limits<double>::quiet_NaN();
    return 0.15 + 0.006 * pt_gamma;
  }

} // namespace PhoIDCuts

// ============================================================================
// RecoilJets analysis module
// ============================================================================
class RecoilJets : public SubsysReco
{
public:
  // -------------------------------------------------------------------------
  // Lightweight categorization used by photon ID QA
  // -------------------------------------------------------------------------
  enum class TightTag : std::uint8_t
  {
    kPreselectionFail = 0,
    kTight            = 1,
    kNonTight         = 2,
    kNeither          = 3
  };

  enum class EventReject : std::uint8_t
  {
    None    = 0,
    Trigger = 1,
    Vz      = 2,
    Period  = 3,
    Ownership = 4
  };

  // Shower-shape variables extracted from PhotonClusterv1
    struct SSVars
    {
      double pt_gamma       = 0.0;
      double weta_cogx      = 0.0;
      double wphi_cogx      = 0.0;
      double weta33_cogx    = std::numeric_limits<double>::quiet_NaN();
      double wphi33_cogx    = std::numeric_limits<double>::quiet_NaN();
      double weta35_cogx    = std::numeric_limits<double>::quiet_NaN();
      double wphi53_cogx    = std::numeric_limits<double>::quiet_NaN();
      double et1            = 0.0;
      double et2            = std::numeric_limits<double>::quiet_NaN();
      double et3            = std::numeric_limits<double>::quiet_NaN();
      double et4            = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e33   = 0.0;
      double e17_over_e77   = std::numeric_limits<double>::quiet_NaN();
      double e32_over_e35   = 0.0;
      double e11_over_e22   = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e13   = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e15   = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e17   = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e31   = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e51   = std::numeric_limits<double>::quiet_NaN();
      double e11_over_e71   = std::numeric_limits<double>::quiet_NaN();
      double e22_over_e33   = std::numeric_limits<double>::quiet_NaN();
      double e22_over_e35   = std::numeric_limits<double>::quiet_NaN();
      double e22_over_e37   = std::numeric_limits<double>::quiet_NaN();
      double e22_over_e53   = std::numeric_limits<double>::quiet_NaN();
      double w32            = std::numeric_limits<double>::quiet_NaN();
      double w52            = std::numeric_limits<double>::quiet_NaN();
      double w72            = std::numeric_limits<double>::quiet_NaN();
      double mean_time      = std::numeric_limits<double>::quiet_NaN();
      double cluster_prob   = std::numeric_limits<double>::quiet_NaN();
      double npb_score      = std::numeric_limits<double>::quiet_NaN();
      double tight_bdt_score = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_n_owned = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_owned_e = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_all_e = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_den_e = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_weta_cogx_num = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_wphi_cogx_num = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_cog_eta = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_cog_phi = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_center_ieta = std::numeric_limits<double>::quiet_NaN();
      double ppg12_shape_center_iphi = std::numeric_limits<double>::quiet_NaN();
    };

  // Per-(trigger,slice) category counters printed in End()
  struct CatStat
  {
    long long n_iso_tight        = 0;
    long long n_nonIso_tight     = 0;
    long long n_iso_nonTight     = 0;
    long long n_nonIso_nonTight  = 0;
    long long seen               = 0;
  };

  // Job-wide bookkeeping counters printed in printCutSummary()
  struct Bookkeeping
  {
    // Event-level
    long long evt_seen         = 0;
    long long evt_fail_trigger = 0;
    long long evt_fail_vz      = 0;
    long long evt_accepted     = 0;

    // Photon-level
    long long pho_total            = 0;
    long long pho_early_E          = 0;
    long long pho_eta_fail         = 0;
    long long pho_etbin_out        = 0;
    long long pho_noRC             = 0;
    long long pho_reached_pre_iso  = 0;

    // Preselection breakdown
    long long pre_pass              = 0;
    long long pre_fail_weta         = 0;
    long long pre_fail_et1_low      = 0;
    long long pre_fail_et1_high     = 0;
    long long pre_fail_e11e33_high  = 0;
    long long pre_fail_e32e35_low   = 0;
    long long pre_fail_e32e35_high  = 0;

    // Tight classification breakdown
    long long tight_tight       = 0;
    long long tight_neither     = 0;
    long long tight_nonTight    = 0;

    // Tight per-cut fails
    long long tight_fail_weta   = 0;
    long long tight_fail_wphi   = 0;
    long long tight_fail_et1    = 0;
    long long tight_fail_e11e33 = 0;
    long long tight_fail_e32e35 = 0;

    // Isolation
    long long iso_pass          = 0;
    long long iso_fail          = 0;
  };

  // -------------------------------------------------------------------------
  // Jet radius mapping (key + node names per dataset)
  // -------------------------------------------------------------------------
  struct JetRadiusDef
  {
    std::string key;      // e.g. "r02"
    std::string pp_node;  // jet node name used in pp
    std::string aa_node;  // jet node name used in Au+Au
  };

  // Radii to run in parallel. These are used in fetchNodes(), jet QA, matching.
  //
  // NOTE: Replace node names here if your DST uses different names.
    inline static const std::array<JetRadiusDef, 2> kJetRadii = {{
      // key   pp_node              aa_node
      { "r02", "AntiKt_Tower_r02",  "AntiKt_Tower_r02_Sub1"  },
      { "r04", "AntiKt_Tower_r04",  "AntiKt_Tower_r04_Sub1"  }
    }};

  // Trigger maps:
  //  - pp:   DB/GL1 name -> short key (directory name)
  //  - AuAu: bit index   -> short key (directory name)
  //
  // NOTE: These are placeholders if you don't already have your project-specific
  //       maps. Keep/restore your exact trigger mappings here.
  inline static const std::vector<std::pair<std::string, std::string>> triggerNameMap_pp = {
      //     {"MBD N&S >= 1",          "MBD_NandS_geq_1"},
      //      {"Photon 3 GeV + MBD NS >= 1","Photon_3_GeV_plus_MBD_NS_geq_1"},
      {"Photon 4 GeV + MBD NS >= 1","Photon_4_GeV_plus_MBD_NS_geq_1"}
      //      {"Photon 5 GeV + MBD NS >= 1","Photon_5_GeV_plus_MBD_NS_geq_1"}
  };

  // doNotScale turn-on bookkeeping is configured in RecoilJets.cc as
  // run-resolved (baseline raw bit, probe live bit) pairs.
  // Pair-specific histograms are filled as:
  //   h_maxEnergyClus_NewTriggerFilling_doNotScale_<probeHistKey>
  //   h_maxEnergyClus_NewTriggerFilling_doNotScale_baseline_<probeHistKey>
      
  inline static const std::vector<std::pair<int, std::string>> triggerNameMapAuAu = {
    //        {10, "MBD_NS_geq_2"},
    //        {11, "MBD_NS_geq_1"},
    //        {12, "MBD_NS_geq_2_vtx_lt_10"},
    //        {13, "MBD_NS_geq_2_vtx_lt_30"},
            {14, "MBD_NS_geq_2_vtx_lt_150"}
    //        {15, "MBD_NS_geq_1_vtx_lt_10"},
    //        {16, "photon_6_plus_MBD_NS_geq_2_vtx_lt_10"},
    //        {17, "photon_8_plus_MBD_NS_geq_2_vtx_lt_10"},
    //        {18, "photon_10_plus_MBD_NS_geq_2_vtx_lt_10"},
    //        {19, "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"},
    //        {20, "photon_6_plus_MBD_NS_geq_2_vtx_lt_150"},
    //        {21, "photon_8_plus_MBD_NS_geq_2_vtx_lt_150"},
    //        {22, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"},
    //        {23, "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"}
  };

  // -------------------------------------------------------------------------
  // Construction / Fun4All hooks
  // -------------------------------------------------------------------------
  explicit RecoilJets(const std::string& outFile,
                      const std::string& moduleName = "RecoilJets");
  ~RecoilJets() override;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // -------------------------------------------------------------------------
  // Configuration helpers (keep names stable; safe defaults)
  // -------------------------------------------------------------------------
  void setDataType(const std::string& s)
  {
      // Supported user-facing keys (also overridden by env RJ_DATASET / RJ_IS_SIM)
      //   "isSim"/"sim"  -> sim
      //   "isAuAu"/"auau"-> Au+Au
      //   "ispp"/"pp"    -> pp
      const std::string t = s;
      if (t == "isSim" || t == "sim" || t == "isSimJet5" || t == "isSimInclusive" || t == "isSimMB") { m_isSim = true;  m_isAuAu = false; }
      if (t == "isAuAu" || t == "auau"|| t == "aa") { m_isSim = false; m_isAuAu = true; }
      if (t == "ispp"  || t == "pp")  { m_isSim = false; m_isAuAu = false; }
  }

  // Optional: limit active jet radii (keys like "r02","r04"). Empty => all kJetRadii.
  void setActiveJetRKeys(const std::vector<std::string>& keys) { m_activeJetRKeys = keys; }

  void setPhotonEtaAbsMax(double v) { m_etaAbsMax = v; }
  void setMinJetPt(double v)        { m_minJetPt = v; }
  void setInternalJetPtScan(const std::vector<double>& values);
  void setMinBackToBack(double v)   { m_minBackToBack = v; }
  void setInternalBackToBackScan(const std::vector<double>& values);

  void setUseVzCut(bool on, double vz = 60.0)
  {
      m_useVzCut = on;
      m_vzCut    = static_cast<float>(vz);
  }

  void setVertexReweighting(bool on, const std::string& filePath, const std::string& histPath)
  {
      m_vertexReweightOn   = on;
      m_vertexReweightFile = filePath;
      m_vertexReweightHist = histPath;
  }

  void setGammaPtBins(const std::vector<double>& bins)
          {
            // Canonical pT^gamma binning for ALL photon-binned histograms (reco + truth):
            //   [15,17,19,21,23,26,35] GeV  (6 bins, start at 15 GeV)
            //
            // NOTE: Unfolding histograms book their own extended pT^gamma binning:
            //   reco : [10,15] + (canonical bins) + [35,40]
            //   truth: [5,10] + [10,15] + (canonical bins) + [35,40]
            if (!bins.empty())
            {
              m_gammaPtBins = bins;
            }
  }

  // Phase-1 YAML knobs (matching thresholds)
  void setPhoMatchDRMax(double v) { m_phoMatchDRMax = v; }
  void setJetMatchDRMax(double v) { m_jetMatchDRMax = v; }

  // Phase-1 YAML knobs (unfolding explicit bin edges)
  void setUnfoldRecoPhotonPtBins(const std::vector<double>& bins)  { m_unfoldRecoPhotonPtBins = bins; }
  void setUnfoldTruthPhotonPtBins(const std::vector<double>& bins) { m_unfoldTruthPhotonPtBins = bins; }
  void setUnfoldJetPtBins(const std::vector<double>& bins)         { m_unfoldJetPtBins = bins; }
  void setUnfoldXJBins(const std::vector<double>& bins)            { m_unfoldXJBins = bins; }
  void setLeadingResponseFamilyLabel(const std::string& label)     { m_leadingResponseFamilyLabel = label; }
  void setPPG12PhotonYieldEnabled(bool on = true)                  { m_ppg12PhotonYieldEnabled = on; }
  void setPPG12PhotonYieldApplyBinning(bool on = true)             { m_ppg12PhotonYieldApplyBinning = on; }
  void setPPG12Fig8ClusterNode(const std::string& node)            { m_ppg12Fig8ClusterNode = node; }
  void usePPG12PhotonYieldBinning()
  {
      m_gammaPtBins = m_ppg12PhotonYieldRecoPtBins;
      m_unfoldRecoPhotonPtBins = m_ppg12PhotonYieldRecoPtBins;
      m_unfoldTruthPhotonPtBins = m_ppg12PhotonYieldTruthPtBins;
  }

  // Analysis provenance stamping (written once into output ROOT by RecoilJets.cc)
  void setAnalysisConfigYAML(const std::string& yamlText, const std::string& tag)
  {
          m_analysisConfigYAMLText = yamlText;
          m_analysisConfigTag      = tag;
  }

  void enablePi0Analysis(bool on = true) { m_doPi0Analysis = on; }

  void setCentEdges(const std::vector<int>& edges)     { m_centEdges = edges; }

  // Isolation WP (implemented in .cc)
  void setIsolationWP(double aGeV, double bPerGeV,
	                                  double sideGapGeV, double coneR, double towerMin,
	                                  double fixedGeV = 2.0);
  void setPPIsoWPForCone(double coneR, double aGeV, double bPerGeV, double sideGapGeV);
  void setIsSlidingIso(bool on) { m_isSlidingIso = on; }
  struct IsoView
  {
      std::string label;
      double coneR = 0.3;
      bool isSliding = false;
      double fixedGeV = 2.0;
  };
  void setInternalIsoViews(const std::vector<IsoView>& views);

  // Truth isolation max (independent of sliding/fixed reco mode)
  void setTruthIsoMaxGeV(double isoGeV)
  {
          if (!std::isfinite(isoGeV)) return;
          m_truthIsoMaxGeV = isoGeV;
          if (m_truthIsoMaxGeV < 0.0) m_truthIsoMaxGeV = 0.0;
  }

  // Photon ID cuts (PPG12 Table 4): allow YAML override while preserving baseline defaults
  void setPhotonIDCuts(double pre_e11e33_max,
                           double pre_et1_min,
                           double pre_et1_max,
                           double pre_e32e35_min,
                           double pre_e32e35_max,
                           double pre_weta_max,
                           double tight_w_lo,
                           double tight_w_hi_intercept,
                           double tight_w_hi_slope,
                           double tight_e11e33_min,
                           double tight_e11e33_max,
                           double tight_et1_min,
                           double tight_et1_max,
                           double tight_e32e35_min,
                           double tight_e32e35_max)
  {
        m_phoid_pre_e11e33_max = pre_e11e33_max;
        m_phoid_pre_et1_min    = pre_et1_min;
        m_phoid_pre_et1_max    = pre_et1_max;
        m_phoid_pre_e32e35_min = pre_e32e35_min;
        m_phoid_pre_e32e35_max = pre_e32e35_max;
        m_phoid_pre_weta_max   = pre_weta_max;

        m_phoid_tight_w_lo           = tight_w_lo;
        m_phoid_tight_w_hi_intercept = tight_w_hi_intercept;
        m_phoid_tight_w_hi_slope     = tight_w_hi_slope;

        m_phoid_tight_e11e33_min = tight_e11e33_min;
        m_phoid_tight_e11e33_max = tight_e11e33_max;

        m_phoid_tight_et1_min    = tight_et1_min;
        m_phoid_tight_et1_max    = tight_et1_max;

        m_phoid_tight_e32e35_min = tight_e32e35_min;
        m_phoid_tight_e32e35_max = tight_e32e35_max;
  }

  void setPhotonIDVariants(const std::string& preselection,
                           const std::string& tight,
                           const std::string& nonTight,
                           const std::string& preselectionPhotonNode = "PHOTONCLUSTER_CEMC",
                           const std::string& tightPhotonNode = "PHOTONCLUSTER_CEMC");

  void setAuAuTightBDTRuntimeConfig(const std::string&,
                                    const std::vector<std::string>&,
                                    const std::vector<int>& = {},
                                    const std::vector<std::string>& = {},
                                    const std::vector<double>& = {},
                                    const std::vector<std::string>& = {},
                                    const std::vector<std::string>& = {},
                                    const std::string& = "",
                                    const std::vector<std::string>& = {},
                                    double = 35.0,
                                    double = 40.0,
                                    double = std::numeric_limits<double>::quiet_NaN(),
                                    double = std::numeric_limits<double>::quiet_NaN(),
                                    const std::vector<std::string>& = {},
                                    const std::string& = "") {}

  void setAuAuTightMLPRuntimeConfig(const std::string&,
                                    double = 0.80,
                                    double = 0.0,
                                    double = 1.0,
                                    double = 0.20,
                                    double = 0.0,
                                    double = 0.80,
                                    double = 0.0,
                                    double = std::numeric_limits<double>::quiet_NaN(),
                                    double = std::numeric_limits<double>::quiet_NaN(),
                                    const std::vector<std::string>& = {},
                                    const std::string& = "") {}

  void setAuAuTightBDTMLPStackRuntimeConfig(const std::string&,
                                            double = 0.80,
                                            double = 0.0,
                                            double = 1.0,
                                            double = 0.20,
                                            double = 0.0,
                                            double = 0.80,
                                            double = 0.0,
                                            double = std::numeric_limits<double>::quiet_NaN(),
                                            double = std::numeric_limits<double>::quiet_NaN(),
                                            const std::vector<std::string>& = {}) {}

  void setAuAuTightLogRegRuntimeConfig(const std::string&,
                                       double = 0.80,
                                       double = 0.0,
                                       double = 1.0,
                                       double = 0.20,
                                       double = 0.0,
                                       double = 0.80,
                                       double = 0.0,
                                       double = std::numeric_limits<double>::quiet_NaN(),
                                       double = std::numeric_limits<double>::quiet_NaN(),
                                       const std::vector<std::string>& = {}) {}

    // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
    // When enabled, a compact per-event-per-radius TTree ("EventDisplayTree") is written to the output ROOT file.
    void enableEventDisplayDiagnostics(bool on = true) { m_evtDiagEnabled = on; }
    void setEventDisplayDiagnosticsMaxPerBin(int n)    { m_evtDiagMaxPerBin = n; }

  // -------------------------------------------------------------------------
  // Read-only state access
  // -------------------------------------------------------------------------
  const Bookkeeping& bookkeeping() const { return m_bk; }
  EventReject lastReject() const         { return m_lastReject; }

private:
  // -------------------------------------------------------------------------
  // Internal helpers (nodes, event selection)
  // -------------------------------------------------------------------------
  bool jetRKeyActive(const std::string& rKey) const
  {
    if (m_activeJetRKeys.empty()) return true;
    return (std::find(m_activeJetRKeys.begin(), m_activeJetRKeys.end(), rKey) != m_activeJetRKeys.end());
  }

  bool fetchNodes(PHCompositeNode* topNode);
  bool firstEventCuts(PHCompositeNode* topNode, std::vector<std::string>& activeTrig);
  void fillPPG12Fig7TriggerQA(PHCompositeNode* topNode);
  void createHistos_Data();
  void fillReplayFoundationCaptureWitness();

  void fillUnfoldResponseMatrixAndTruthDistributions(
            const std::vector<std::string>& activeTrig,
            const std::string& rKey,
            const int effCentIdx_M,
            const double leadPtGamma,
            const double leadEtaGamma,
            const double leadPhiGamma,
            const bool haveTruthPho,
            const double tPt,
            const double tEta,
            const double tPhi,
            const std::vector<const Jet*>& recoJetsFid,
            const std::vector<char>& recoJetsFidIsRecoil,
            const Jet* recoil1Jet);

      
  void fillRecoTruthJES3MatchingQA(const std::vector<std::string>& activeTrig,
                                      const std::string& rKey,
                                      const int effCentIdx_M,
                                      const double leadPtGamma,
                                      const double xJ,
                                      const double alpha,
                                      const double tPt,
                                      const double tEta,
                                      const double tPhi,
                                      const Jet* recoil1Jet);

  bool runLeadIsoTightPhotonJetLoopAllRadii(
              const std::vector<std::string>& activeTrig,
              const int effCentIdx_M,
              const int centIdx,
              const int leadPhoIndex,
              const int leadPtIdx,
              const double leadPtGamma,
              const double leadEtaGamma,
              const double leadPhiGamma,
              const bool haveTruthPho,
              const double tPt,
              const double tEta,
              const double tPhi,
              PHG4TruthInfoContainer* truth);

  void runLeadIsoNonTightPhotonJetLoopAllRadii_SidebandC(
              const std::vector<std::string>& activeTrig,
              const int effCentIdx_M,
              const double leadPtGamma,
              const double leadEtaGamma,
              const double leadPhiGamma);

        
  bool runLeadIsoTightPhotonJetMatchingAndUnfolding(
                const std::vector<std::string>& activeTrig,
                const int effCentIdx_M,
                const int centIdx,
                const int leadPhoIndex,
                const int leadPtIdx,
                const double leadPtGamma,
                const double leadEtaGamma,
                const double leadPhiGamma,
                const bool haveTruthSigPho,
                const double tPtSig,
                const bool haveTruthPhoPPG12,
                const double recoPtTruthMatchPPG12,
                const bool haveTruthPho,
                const double tPt,
                const double tEta,
                const double tPhi,
                PHG4TruthInfoContainer* truth);

  void fillPureIsolationQA(PHCompositeNode* topNode,
                               const std::vector<std::string>& activeTrig,
                               const PhotonClusterv1* pho,
                               const RawCluster* rc,
                               const int ptIdx,
                               const int centIdx,
                               const double pt_gamma);

  void fillTruthSigABCDLeakageCounters(PHCompositeNode* topNode,
                                           const std::vector<std::string>& activeTrig,
                                           const int centIdx);

  void processCandidates(PHCompositeNode* topNode, const std::vector<std::string>& activeTrig);
  void initPPPhotonIDTrainingTree();
  void fillPPPhotonIDTrainingTree(const SSVars& v,
                                  double eta,
                                  double phi,
                                  int clusterIndex,
                                  double scoreInputEt,
                                  double responseEt,
                                  double eiso,
                                  double ppg12RawEiso,
                                  double ppg12RecoEiso,
                                  double ppg12IsoThreshold,
                                  double ppg12NonIsoThreshold,
                                  bool ppg12Iso,
                                  bool ppg12NonIso,
                                  bool ppg12CommonPass,
                                  int ppg12TightTag,
                                  int ptIdx,
                                  bool isSignal,
                                  int truthClass,
                                  int ppg12AnalysisWindowPass,
                                  int ppg12ResponseWindowPass,
                                  int ppg12LogicalABCDRegion,
                                  int ppg12SignalFillA,
                                  int ppg12SignalFillB,
                                  int ppg12SignalFillC,
                                  int ppg12SignalFillD,
                                  int ppg12SignalFillMultiplicity,
                                  int truthTrackId,
                                  int truthBarcode,
                                  float truthEnergyContribution,
                                  int ppg12SampleBin = 0,
                                  float ppg12XsecPb = -999.0f,
                                  float ppg12XsecWeight = 1.0f,
                                  float ppg12WindowLow = -999.0f,
                                  float ppg12WindowHigh = -999.0f,
                                  float maxTruthJetPtR04 = -999.0f,
                                  int ppg12TruthWindowPassR04 = -1,
                                  double truthIsoEtR03 = std::numeric_limits<double>::quiet_NaN(),
                                  double truthIsoEtR04 = std::numeric_limits<double>::quiet_NaN(),
                                  int truthIsoValid = -1,
                                  int truthHepMCAssociationValid = -1);
  void fillPi0MassVsPtHistograms(const std::string& trig, RawClusterContainer* clusterContainer, bool useCorr);

  bool getCentralitySlice(int& lo, int& hi, std::string& tag) const;

  // Au+Au scaled trigger helper (safe default behavior)
  static std::bitset<64> extractTriggerBits(std::uint64_t scaledVec, int /*eventNumber*/)
  {
    return std::bitset<64>(scaledVec);
  }
  static bool checkTriggerCondition(const std::bitset<64>& bits, int bit)
  {
    if (bit < 0 || bit >= 64) return false;
    return bits.test(static_cast<std::size_t>(bit));
  }

  // -------------------------------------------------------------------------
  // Photon ID helpers
  // -------------------------------------------------------------------------
    SSVars makeSSFromPhoton(const PhotonClusterv1* pho, double pt_gamma) const;
    void attachVariantScoresToSSVars(const PhotonClusterv1* pho, SSVars& v) const;
    const PhotonClusterv1* findMatchedPhotonByKinematics(const RawClusterContainer* container,
                                                         const PhotonClusterv1* ref) const;
    bool   passesPhotonPreselection(const SSVars& v);
    TightTag classifyPhotonTightness(const SSVars& v);
    bool passesPPG12PhotonYieldCommon(const SSVars& v) const;
    TightTag classifyPPG12PhotonYieldTightness(const SSVars& v) const;
    double configuredTightBDTMin(double et) const;
    double configuredNonTightBDTMin(double et) const;
    double configuredNonTightBDTMax(double et) const;
    bool configuredTightBDTPass(double score, double et) const;
    bool configuredNonTightBDTPass(double score, double et) const;
    TightTag configuredVariantABDTTag(double score, double et) const;


    // Isolation helpers
    struct PPG12EisoConeMemberAudit
    {
      std::uint64_t key = 0;
      double energy = std::numeric_limits<double>::quiet_NaN();
      double eta = std::numeric_limits<double>::quiet_NaN();
      double phi = std::numeric_limits<double>::quiet_NaN();
      double et = std::numeric_limits<double>::quiet_NaN();
      double dEta = std::numeric_limits<double>::quiet_NaN();
      double dPhi = std::numeric_limits<double>::quiet_NaN();
      double dR = std::numeric_limits<double>::quiet_NaN();
    };
    struct PPG12EisoPathAudit
    {
      double storedRaw = std::numeric_limits<double>::quiet_NaN();
      double legacyPositiveOnlyRaw = std::numeric_limits<double>::quiet_NaN();
      double storedTopoSumEt = std::numeric_limits<double>::quiet_NaN();
      double storedPositiveOnlyTopoSumEt = std::numeric_limits<double>::quiet_NaN();
      double storedValid = std::numeric_limits<double>::quiet_NaN();
      double storedVertexZ = std::numeric_limits<double>::quiet_NaN();
      double expectedVertexZ = std::numeric_limits<double>::quiet_NaN();
      double isoAxisEta = std::numeric_limits<double>::quiet_NaN();
      double isoAxisPhi = std::numeric_limits<double>::quiet_NaN();
      double recomputedTopoSumEt = std::numeric_limits<double>::quiet_NaN();
      double recomputedPositiveOnlyTopoSumEt = std::numeric_limits<double>::quiet_NaN();
      double recomputedRaw = std::numeric_limits<double>::quiet_NaN();
      double finalRaw = std::numeric_limits<double>::quiet_NaN();
      double candidateEt = std::numeric_limits<double>::quiet_NaN();
      long long topoNodeCount = 0;
      long long inConeMemberCount = 0;
      bool vertexCompatible = false;
      bool storedUsable = false;
      bool storedAccepted = false;
      bool storedRejected = false;
      bool recomputeFallbackUsed = false;
      bool vertexMismatchRejected = false;
      bool recoVertexUsedForEiso = false;
      int pathCode = 0;  // 0=invalid, 1=stored, 2=recomputed
      std::vector<PPG12EisoConeMemberAudit> coneMembers;
    };
    double eiso(const RawCluster* clus, PHCompositeNode* topNode) const;
    double ppg12PhotonYieldRawEiso(const RawCluster* clus, PHCompositeNode* topNode);
    double ppg12PhotonYieldRawEisoDetailed(const RawCluster* clus,
                                           PHCompositeNode* topNode,
                                           PPG12EisoPathAudit* audit);
    double ppg12PhotonYieldEiso(double eisoEt) const;
    double ppg12PhotonYieldClusterEtForCuts(double ptGamma, int candidateIndex) const;
    double ppg12PhotonYieldClusterEtForResponse(double recoEt,
                                                double truthPt,
                                                int candidateIndex) const;
    bool ppg12PhotonYieldTowerMasked(const PhotonClusterv1* pho) const;
    double ppg12PhotonYieldKinematicVertexZ() const;
    bool ppg12PhotonYieldInResponseWindow(double recoPt, double truthPt) const;
    double ppg12Fig36TruthPriorWeight(double truthPt) const;
    bool loadPPG12PhotonYieldTowerMask();
    int currentRunNumber() const;
    bool ppg12PeriodRunContains(int runNumber) const;
    void fillPPG12VertexContractQA(const std::vector<std::string>& activeTrig,
                                   bool preVzCut);
    void recordPPG12EisoVertexContractAudit(const PPG12EisoPathAudit& audit);
    bool configurePPG12SimEventWeight(double sliceFactor, int laneCode);
    void fillPPG12SimWeightAudit();
    void finalizePPG12SimWeightAudit();
    void recordPPG12EisoPathCanary(const PPG12EisoPathAudit& audit,
                                   int candidateIndex,
                                   double recoEt,
                                   double eta,
                                   double phi,
                                   bool isTruthSignal,
                                   const SSVars& vars,
                                   double correctedEiso,
                                   double isoThreshold,
                                   double nonIsoLower,
                                   TightTag tightTag);
  bool   isIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode);
  bool   isNonIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode);

    // Unified truth-MC signal definition for "isolated prompt photon" (SIM only)
    // Definition
    //   |eta| < 0.7, PID=22, embedded G4 primary photon; HepMC association uses (embed ID, barcode),
    //   prompt classification via CaloAna photon_type logic:
    //     - walk back photon-in/photon-out vertices
    //     - direct=1 if 2->2 with |pdg|<=22 on all legs
    //     - frag  =2 if 1->2 with |incoming pdg|<=11 and outgoing contains incoming pid (and photon)
    //   valid unknown classes -1/0 are also accepted (class < 3),
    //   and truth isolation ETiso_truth < 4 GeV where (Blair/Shuhang CaloAna truth-iso):
    //     ETiso = sum_{ΔR<0.3} Et(G4 primary particles with embed>=1)
    //           - sum_{ΔR<0.001} Et(G4 primary particles with embed>=1)
    //     (the ΔR<0.001 subtraction removes the photon itself, and any ultra-merged pieces).
      
  struct TruthSignalPhotonInfo
  {
    int trackId = -1;
    int barcode = -1;
    double pt = std::numeric_limits<double>::quiet_NaN();
    double eta = std::numeric_limits<double>::quiet_NaN();
    double phi = std::numeric_limits<double>::quiet_NaN();
    double isoEt = std::numeric_limits<double>::quiet_NaN();
    double isoEtR03 = std::numeric_limits<double>::quiet_NaN();
    double isoEtR04 = std::numeric_limits<double>::quiet_NaN();
    bool truthIsolationValid = false;
    bool g4PhotonValid = false;
    bool hepmcAssociationValid = false;
    int photonClass = -999;
    int embedId = -999;
    const HepMC::GenParticle* hep = nullptr;
    const PHG4Particle* g4 = nullptr;
  };
  using TruthSignalPhotonMap = std::map<int, TruthSignalPhotonInfo>;

  // Selection-neutral truth capture. One embedded-primary scan computes both
  // R=0.3 and R=0.4 raw cone sums; the signal map below applies the nominal
  // PPG12 R=0.3/class/eta cuts to this retained witness population.
  TruthSignalPhotonMap buildPPG12TruthPhotonWitnessMap(
      PHCompositeNode* topNode,
      RJReplayFoundationV1::TruthPhotonInventoryV1* inventory=nullptr) const;

  // PPG12 photon truth association:
  //   cluster_truthtrkID = clustereval.max_truth_primary_particle_by_energy(cluster)->get_track_id()
  //   signal if that track id maps to a truth photon in the nominal G4-driven signal map.
  // This deliberately does not impose an additional truth-reco ΔR cut.
  TruthSignalPhotonMap buildPPG12TruthSignalPhotonMap(PHCompositeNode* topNode) const;
  bool classifyRecoPhotonWithPPG12TruthTrack(const RawCluster* rc,
                                             CaloRawClusterEval& clustereval,
                                             const TruthSignalPhotonMap& truthSignalByTrackId,
                                             TruthSignalPhotonInfo& matchedTruth,
                                             int& clusterTruthTrackId,
                                             float& eContrib) const;

  // Compatibility helper for callers that start from a truth photon. Uses the
  // same PPG12 track-id association above and returns the reco photon with the
  // largest contribution for this truth photon. drBest is diagnostic only.
  bool findRecoPhotonMatchedToTruthSignal(const TruthSignalPhotonMap& truthSignalByTrackId,
                                             int targetTrackId,
                                             CaloRawClusterEval& clustereval,
                                             const RawCluster*& recoPho,
                                             double& recoPt,
                                             double& recoEta,
                                             double& recoPhi,
                                             double& drBest,
                                             float& eContribBest) const;

  // PPG12 truth-tagging for SS template overlays (SIM only)
  //   - signal: reco cluster best-matched to a truth photon in the nominal G4-driven signal map
  //   - background: everything else (complement of signal)
  const HepMC::GenParticle* findHepMCParticleByBarcode(const HepMC::GenEvent* evt, int bc) const;


  // -------------------------------------------------------------------------
  // Binning helpers
  // -------------------------------------------------------------------------
  int         findPtBin(double pt) const;
  int         findCentBin(int cent) const;
  std::string suffixForBins(int ptIdx, int centIdx) const;
  void        configureInternalJetPtScanFromEnv();
  std::vector<double> activeJetPtCuts() const;
  std::string jetPtKeyForCut(double ptGeV) const;
  std::string histRKeyForJetPt(const std::string& rKey, double ptGeV) const;
  std::string histRKeyForJetPtAndDphi(const std::string& rKey, double ptGeV, double cutRad) const;
  void        configureInternalIsoViewsFromEnv();
  bool        getPPIsoParamsForCone(double coneR, double& aGeV, double& bPerGeV, double& sideGapGeV) const;
  void        getActiveIsoParams(double& aGeV, double& bPerGeV, double& sideGapGeV) const;
  double      recoIsoThreshold(double ptGamma) const;
  double      recoNonIsoThreshold(double ptGamma) const;
  enum class HistViewScope { Canonical, IsoCone, IsoView };
  bool        fillCanonicalThisView() const;
  bool        fillConeThisView() const;
  std::string histBaseForScope(const std::string& base, HistViewScope scope) const;
  std::string recoilRKeyForScope(const std::string& rKey, double ptGeV, double cutRad, HistViewScope scope) const;
  std::string withIsoConeSuffix(const std::string& base) const;
  std::string withIsoViewSuffix(const std::string& base) const;
  std::string stripIsoViewSuffix(const std::string& key) const;
  void        configureInternalBackToBackScanFromEnv();
  std::vector<double> activeBackToBackCuts() const;
  std::string dphiKeyForCut(double cutRad) const;
  std::string histRKeyForDphi(const std::string& rKey, double cutRad) const;
  std::string baseRKeyFromDphiHistKey(const std::string& rKey) const;

  // EventDisplay categories (used by the diagnostics payload stored in EventDisplayTree)
  enum class EventDisplayCat : int { NUM = 0, MissA = 1, MissB = 2 };

  // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
  void        initEventDisplayDiagnosticsTree();
  void        resetEventDisplayDiagnosticsBuffers();
  bool        eventDisplayDiagnosticsNeed(const std::string& rKey, int ptBin, EventDisplayCat cat);
  void        appendEventDisplayDiagnosticsFromJet(const Jet* jet,
                                                      std::vector<int>& calo,
                                                      std::vector<int>& ieta,
                                                      std::vector<int>& iphi,
                                                      std::vector<float>& eta,
                                                      std::vector<float>& phi,
                                                      std::vector<float>& et,
                                                      std::vector<float>& e) const;
    
    void        fillEventDisplayDiagnostics(const std::string& rKey,
                                               int ptBin,
                                               EventDisplayCat cat,
                                               double truthGammaPt,
                                               double truthGammaPhi,
                                               double recoGammaPt,
                                               double recoGammaEta,
                                               double recoGammaPhi,
                                               const Jet* selectedRecoilJet,
                                               const Jet* recoTruthBest,
                                               const Jet* truthLeadRecoilJet);
  // -------------------------------------------------------------------------
  // Histogram utilities (bookers)
  // -------------------------------------------------------------------------
  TH1I* getOrBookCountHist(const std::string& trig,
                             const std::string& base,
                             int ptIdx, int centIdx,
                             HistViewScope scope = HistViewScope::IsoView);

  // Isolation spectra (reco)
  TH1F* getOrBookIsoHist(const std::string& trig, int ptIdx, int centIdx);
  TH1F* getOrBookIsoPartHist(const std::string& trig,
                                     const std::string& base,
                                     const std::string& xAxisTitle,
                                     int ptIdx, int centIdx,
                                     HistViewScope scope = HistViewScope::IsoCone);
  TH1F* getOrBookPtGammaHist(const std::string& trig,
                                 const std::string& base,
                                 int centIdx);
  TH1F* getOrBookPPG12PhotonYield1D(const std::string& trig,
                                    const std::string& name,
                                    const std::vector<double>& bins,
                                    const std::string& xAxisTitle,
                                    const std::string& yAxisTitle);
  TH2F* getOrBookPPG12PhotonYieldResponse2D(const std::string& trig,
                                            const std::string& name);
  TH2F* getOrBookPPG12Fig8Response2D(const std::string& trig,
                                     const std::string& name,
                                     const std::string& yAxisTitle);
  TH2F* getOrBookPPG12Fig8Primitive2D(const std::string& trig,
                                      const std::string& name,
                                      const std::string& xAxisTitle,
                                      const std::string& yAxisTitle);
  TH1I* getOrBookPPG12Fig8AuditHist(const std::string& trig);
  TH1I* getOrBookPPG12EisoVertexContractAuditHist(const std::string& trig);
  TH1F* getOrBookPPG12VertexQA1D(const std::string& trig,
                                 const std::string& name,
                                 const std::string& title,
                                 int nbins,
                                 double xmin,
                                 double xmax);
  TH1I* getOrBookPPG12VertexQAAuditHist(const std::string& trig);
  TH1D* getOrBookPPG12PeriodContractValues(const std::string& trig);
  void fillPPG12Fig81DataHDReferenceVertexQA(PHCompositeNode* topNode);
  void fillPPG12Fig81SlimtreeVertexQA(PHCompositeNode* topNode);
  void bookPPG12PhotonYieldSchema(const std::vector<std::string>& activeTrig);
  void fillPPG12Fig8ResponseTargets(const std::vector<std::string>& activeTrig,
                                    CaloRawClusterEval& clustereval,
                                    const TruthSignalPhotonMap& truthSignalByTrackId);
  void fillPPG12PhotonYieldAllCommon(const std::vector<std::string>& activeTrig,
                                     double ptGamma,
                                     bool commonPass,
                                     bool haveTruthClass,
                                     bool isTruthSignal,
                                     double weight,
                                     bool signalResponseWindow = true,
                                     bool fillInclusive = true);
  void fillPPG12PhotonYieldTightAndABCD(const std::vector<std::string>& activeTrig,
                                        double ptGamma,
                                        double responsePtGamma,
                                        double truthPt,
                                        TightTag tightTag,
                                        bool iso,
                                        bool nonIso,
                                        bool haveTruthClass,
                                        bool isTruthSignal,
                                        double weight,
                                        bool fillInclusive = true);
  TH1I* getOrBookIsoDecisionHist(const std::string& trig, int ptIdx, int centIdx);

  // Event-level photon multiplicity diagnostic:
  // N = number of reco photon candidates passing the SAME tight+iso+fiducial cuts
  // used for the leading photon selection (filled once per event; binned by leading-photon pT bin).
  TH1I* getOrBookNIsoTightPhoCandHist(const std::string& trig, int ptIdx, int centIdx);

  TH2F* getOrBookIsoCompareHist(const std::string& trig, int ptIdx, int centIdx);

  // -------------------------------------------------------------------------
  // SIM ONLY: matched truth-signal → reco ABCD leakage counters (per pT[/cent] slice)
  //   bins: 1=A, 2=B, 3=C, 4=D
  // -------------------------------------------------------------------------
  TH1* getOrBookSigABCDLeakageHist(const std::string& trig, int ptIdx, int centIdx,
                                   const std::string& base = "h_sigABCD_MC");

  // Isolation QA (truth, SIM only)
  // These are NOT sliced; they live in the trigger directory (SIM => /SIM/).
  TH1F* getOrBookTruthIsoHist(const std::string& trig,
                                const std::string& name,
                                int nbins, double xmin, double xmax);

  TH1I* getOrBookTruthIsoDecisionHist(const std::string& trig,
                                        const std::string& name);

  // Shower-shape distributions
  TH1F* getOrBookSSHist(const std::string& trig,
                        const std::string& varKey,
                        const std::string& tagKey,
                        int ptIdx, int centIdx,
                        HistViewScope scope = HistViewScope::IsoView);
  TH1F* getOrBookPPG12TableQA1DHist(const std::string& trig,
                                    const std::string& varKey,
                                    const std::string& ptToken,
                                    const std::string& centToken,
                                    int cutIdx);
  TH2F* getOrBookPPG12TableQAHist(const std::string& trig,
                                  const std::string& varKey,
                                  const std::string& ptToken,
                                  const std::string& centToken,
                                  int cutIdx);
  TH2D* getOrBookPPG12TimingHist(const std::string& trig,
                                 const std::string& name,
                                 const std::string& title,
                                 int xBins, double xMin, double xMax,
                                 int yBins, double yMin, double yMax);
  TH1F* getOrBookPPG12Fig13E11E33Hist(const std::string& trig,
                                      const std::string& sampleKey,
                                      const std::string& ptToken,
                                      int cutIdx);
  TH2F* getOrBookPPG12Fig13E11E33IsoHist(const std::string& trig,
                                         const std::string& sampleKey,
                                         const std::string& ptToken,
                                         int cutIdx);
  TH1F* getOrBookPPG12IsoStackHist(const std::string& trig,
                                   const std::string& tightness,
                                   int ptIdx);
  TH2F* getOrBookPPG12Fig11ETIsoHist(const std::string& trig,
                                     const std::string& sampleKey);
  TH1D* getOrBookPPG12Fig11AuditHist(const std::string& trig);
  int findPPG12Fig24RecoPtBin(double recoEt) const;
  std::string ppg12Fig24RecoPtSuffix(int ptIdx) const;
  TH2F* getOrBookPPG12Fig24ETIsoHist(const std::string& trig,
                                     const std::string& name);
  TH1F* getOrBookPPG12Fig24IsoHist(const std::string& trig,
                                   int ptIdx);
  TH1D* getOrBookPPG12Fig24AuditHist(const std::string& trig);
  void fillPPG12Fig11SBDiagnostic(const std::vector<std::string>& activeTrig,
                                  const SSVars& v,
                                  double eta,
                                  double recoEtForBinning,
                                  double rawTopoEisoEt,
                                  bool ppPhotonSignalContext,
                                  bool isPPG12Signal,
                                  double photonSampleWeight,
                                  bool ppInclusiveJetContext,
                                  bool ppInclusiveJetPassR04,
                                  double inclusiveJetSampleWeight);
  void fillPPG12Fig13E11E33(const std::vector<std::string>& activeTrig,
                            const SSVars& v,
                            double eisoEt,
                            const std::string& sampleKey,
                            int cutIdx,
                            double ptForBinning = std::numeric_limits<double>::quiet_NaN());
  void fillPPG12IsoStackQA(const std::vector<std::string>& activeTrig,
                           double ptGamma,
                           double eisoEt,
                           TightTag tightTag,
                           double weight = 1.0);
  double ppg12DataTriggerEfficiencyWeight(double recoEt) const;
  void fillPPG12DataTriggerEfficiencyAudit(const std::vector<std::string>& activeTrig,
                                           double recoEt,
                                           double correctionWeight);
  void bookPPG12TableQASchema(const std::vector<std::string>& activeTrig);
  void fillPPG12TableQA(const std::vector<std::string>& activeTrig,
                        const SSVars& v,
                        double eisoEt,
                        int centIdx,
                        int cutIdx,
                        double rowWeight = 1.0,
                        double ptForBinning = std::numeric_limits<double>::quiet_NaN());
  bool computePPG12ClusterMbdTiming(const SSVars& v,
                                    double& clusterMbdDeltaT,
                                    double& clusterTime,
                                    double& mbdTime) const;
  bool passPPG12Fig14NpbMbdTimingSelection(const SSVars& v,
                                           double phi) const;
  void fillPPG12Fig14Fig15TimingQA(const std::vector<std::string>& activeTrig,
                                   const SSVars& v,
                                   double eta,
                                   double phi,
                                   double ptForBinning = std::numeric_limits<double>::quiet_NaN(),
                                   double rowWeight = 1.0);
  double ppg12TableQABDTScore(const SSVars& v) const;
  bool isPPG12TableQADataNPBTaggedCluster(const SSVars& v,
                                          double phi,
                                          double& clusterMbdDeltaT,
                                          double& mbdTime,
                                          bool& hasAwayJet,
                                          bool requireTaggingEnabled = true) const;
  void loadPPG12TableQAMbdT0Corrections(const std::string& path);
  double ppg12TableQAMbdT0OffsetForRun(int runNumber) const;
  TH1F* getOrBookBDTScoreHist(const std::string& trig,
                              const std::string& base,
                              int ptIdx, int centIdx);

  // Physics outputs (already radius-tagged in your .cc)
  TH1F* getOrBookXJHist(const std::string& trig,
                          const std::string& rKey,
                          int ptIdx, int centIdx);
  TH1F* getOrBookJet1PtHist(const std::string& trig,
                              const std::string& rKey,
                              int ptIdx, int centIdx);
  TH1F* getOrBookJet2PtHist(const std::string& trig,
                              const std::string& rKey,
                              int ptIdx, int centIdx);
  TH1F* getOrBookAlphaHist(const std::string& trig,
                             const std::string& rKey,
                             int ptIdx, int centIdx);

  // -------------------------------------------------------------------------
  // inclusive γ–jet unfolding histograms (ATLAS-style 2D in (pT^γ, xJγ))
  //
  //  - Reco (DATA + MC):
  //      h2_unfoldReco_pTgamma_xJ_incl_<rKey><centSuffix>
  //
  //  - Truth (MC):
  //      h2_unfoldTruth_pTgamma_xJ_incl_<rKey><centSuffix>
  //
  //  - Response (MC, matched truth↔reco pairs):
  //      h2_unfoldResponse_pTgamma_xJ_incl_<rKey><centSuffix>
  //    stored as global-bin(truth) vs global-bin(reco), using TH2::FindBin().
  //
  //  - Fakes/Misses (MC, for closure / response completeness):
  //      h2_unfoldRecoFakes_pTgamma_xJ_incl_<rKey><centSuffix>
  //      h2_unfoldTruthMisses_pTgamma_xJ_incl_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  TH2F* getOrBookUnfoldRecoPtXJIncl      (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoPtXJInclSidebandC(const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthPtXJIncl     (const std::string& trig, const std::string& rKey, int centIdx);

  // inclusive |Δphi(gamma,jet)| per recoil jet that passes pT+eta+recoil Δphi cuts
  TH2F* getOrBookUnfoldRecoPtDphiIncl    (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthPtDphiIncl   (const std::string& trig, const std::string& rKey, int centIdx);

  TH2F* getOrBookUnfoldResponsePtXJIncl  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoFakesPtXJIncl (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthMissesPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);

  // Event-leading recoil-jet unfolding family, separate from the inclusive
  // per-photon recoil family. `component` is reco, truth, response, fake,
  // miss, or selectionLoss.
  TH2F* getOrBookUnfoldLeadPtXJComponent(const std::string& trig,
                                         const std::string& rKey,
                                         int centIdx,
                                         const std::string& component);

  // photon-only unfolding (for N_gamma normalization): DATA+SIM reco, SIM truth
  TH1F* getOrBookUnfoldRecoPhoPtGamma        (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldTruthPhoPtGamma       (const std::string& trig, int centIdx);
  TH2F* getOrBookUnfoldResponsePhoPtGamma    (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldRecoPhoFakesPtGamma   (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldTruthPhoMissesPtGamma (const std::string& trig, int centIdx);

  // exploratory 1D photon-only family:
  // same reco-leading photon spectrum, but response/miss/fake bookkeeping uses the
  // PPG12-style object match (matched reco iso∧tight photon, not necessarily leading).
  TH1F* getOrBookUnfoldRecoPhoPtGammaPPG12Obj        (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldTruthPhoPtGammaPPG12Obj       (const std::string& trig, int centIdx);
  TH2F* getOrBookUnfoldResponsePhoPtGammaPPG12Obj    (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldRecoPhoFakesPtGammaPPG12Obj   (const std::string& trig, int centIdx);
  TH1F* getOrBookUnfoldTruthPhoMissesPtGammaPPG12Obj (const std::string& trig, int centIdx);

  // Explicit Fig.37-style response inputs using the xJgamma leading-photon tag.
  // These mirror the existing photon-only unfolding bookkeeping, but with stable
  // names for downstream PPG12 iteration-stability reproduction checks.
  TH1F* getOrBookPPG12Fig37LeadTagInput1D(const std::string& trig,
                                          const std::string& name,
                                          const std::vector<double>& bins,
                                          const std::string& xAxisTitle,
                                          const std::string& yAxisTitle);
  TH2F* getOrBookPPG12Fig37LeadTagResponse2D(const std::string& trig,
                                             const std::string& name);
  void fillPPG12Fig37LeadTagInput1D(const std::vector<std::string>& activeTrig,
                                    const std::string& name,
                                    const std::vector<double>& bins,
                                    const std::string& xAxisTitle,
                                    const std::string& yAxisTitle,
                                    double value);
  void fillPPG12Fig37LeadTagResponse(const std::vector<std::string>& activeTrig,
                                     const std::string& name,
                                     double truthPt,
                                     double recoPt);

  // unfolding QA helpers (SIM only): explicit matched distributions + type-split fakes/misses
  TH2F* getOrBookUnfoldTruthMatchedPtXJIncl      (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoMatchedPtXJIncl       (const std::string& trig, const std::string& rKey, int centIdx);

  // dedicated post-unfold jet-efficiency inputs (SIM only)
  //   Den: truth recoil jets in truth unfolding phase space
  //   Num: subset with a selection-aware reco jet match, built independently
  //        of the geometry-first response MissA/MissB taxonomy
  TH2F* getOrBookUnfoldJetEffDenPtXJIncl         (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldJetEffNumPtXJIncl         (const std::string& trig, const std::string& rKey, int centIdx);

  TH2F* getOrBookUnfoldRecoFakesPtXJIncl_typeA   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldRecoFakesPtXJIncl_typeB   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthMissesPtXJIncl_typeA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldTruthMissesPtXJIncl_typeB (const std::string& trig, const std::string& rKey, int centIdx);

  // jet-match QA (SIM only): match ΔR and pT response for matched recoil jets
  TH1F* getOrBookUnfoldJetMatchDR           (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookUnfoldJetPtResponsePtTruth (const std::string& trig, const std::string& rKey, int centIdx);

  // additional JES diagnostics (SIM only):
  //  - ALL matched fiducial jet pairs (no recoil/truth recoil-selection gating)
  //  - lead recoil jet1-only response + ΔR sanity
  //  - lead recoil jet1 scatter: pT(reco) vs pT(truth)
  TH2F* getOrBookUnfoldJetPtResponseAllPtTruth   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadRecoilJetPtResponsePtTruth  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadRecoilJetPtTruthPtReco      (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadRecoilJetMatchDR            (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // (SIM ONLY): JES3-style *leading truth recoil jet1* match bookkeeping vs truth pT^gamma
  //   Den  : truth leading recoil jet1 exists (truth recoil definition)
  //   Num  : Den + reco recoil jet1 matches truth jet1 (ΔR < m_jetMatchDRMax)
  //   MissA: Den + some reco fid jet within ΔR < m_jetMatchDRMax of truth jet1, but Num failed
  //   MissB: Den + no reco fid jet within ΔR < m_jetMatchDRMax of truth jet1
  // -------------------------------------------------------------------------
  TH1F* getOrBookLeadTruthRecoilMatchDenPtGammaTruth    (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchNumPtGammaTruth    (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchMissA_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchMissB_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);

  // NEW: MissA subtypes (SIM-only)
  //   MissA1: truth-matched reco jet passes recoil definition → competitor/ordering
  //   MissA2: truth-matched reco jet fails recoil definition  → gate-exclusion
  TH1F* getOrBookLeadTruthRecoilMatchMissA1_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);
  TH1F* getOrBookLeadTruthRecoilMatchMissA2_PtGammaTruth (const std::string& trig, const std::string& rKey, int centIdx);
  TH1I* getOrBookLeadTruthRecoilMatchMissA2_Cutflow      (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // (SIM ONLY): diagnostics for why reco xJ differs from truth-conditioned xJ
  //   Filled inside the same DEN / NUM / MissA / MissB classification block.
  //
  //   NOTE: all are radius-tagged and use centrality-only suffix (like the
  //   existing LeadTruthRecoilMatch bookkeeping).
  // -------------------------------------------------------------------------

  // (A1) pT(recoilJet1^reco) vs pT(truth-leading recoil jet), split by class
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_num    (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA1 (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missA2 (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtTruthLead_missB  (const std::string& trig, const std::string& rKey, int centIdx);

  // (A2) pT(recoilJet1^reco) vs pT(reco jet matched to truth-leading recoil jet), for NUM / MissA
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_num    (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA  (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA1 (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchPtRecoJet1VsPtRecoTruthMatch_missA2 (const std::string& trig, const std::string& rKey, int centIdx);

  // (B3) |Δphi(γ^truth, recoilJet1^reco)| vs pTγ,truth, split by class
  TH2F* getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_num   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDphiRecoJet1VsPtGammaTruth_missB (const std::string& trig, const std::string& rKey, int centIdx);

  // (B4) ΔR(recoilJet1^reco, truth-leading recoil jet) vs pTγ,truth, split by class
  TH2F* getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_num   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchDRRecoJet1VsPtGammaTruth_missB (const std::string& trig, const std::string& rKey, int centIdx);

  // (C5) xJ(recoilJet1^reco) vs |Δphi(γ^truth, recoilJet1^reco)|, split by class
  TH2F* getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_num   (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missA (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F* getOrBookLeadTruthRecoilMatchXJRecoJet1VsDphiRecoJet1_missB (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  //  radius-tagged matching-QA bookers
  //   - centrality suffix only (Au+Au)
  //   - pT^gamma is an axis
  //   - name pattern: <base>_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  TH2F*     getOrBookMatchStatusVsPtGamma     (const std::string& trig, const std::string& rKey, int centIdx);

  // per-jet cutflow status vs pT^gamma (NOT event-level; not iso/tight-conditioned)
  // name: h_jetcutflow_status_vs_pTgamma_<rKey><centSuffix>
  // y bins:
  //   1 = FailJetPt, 2 = FailJetEta, 3 = FailBackToBack, 4 = PassAll
  TH2F*     getOrBookJetCutflowStatusVsPtGamma(const std::string& trig, const std::string& rKey, int centIdx);

  TH2F*     getOrBookMatchMaxDphiVsPtGamma    (const std::string& trig, const std::string& rKey, int centIdx);
  TProfile* getOrBookNRecoilJetsVsPtGamma     (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F*     getOrBookMatchDphiVsPtGamma       (const std::string& trig, const std::string& rKey, int centIdx);
  TH2F*     getOrBookRecoilIsLeadingVsPtGamma (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // radius-tagged JES3 bookers (data)
  //   - centrality suffix only (Au+Au)
  //   - name pattern: <base>_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  TH3F* getOrBookJES3_xJ_alphaHist     (const std::string& trig, const std::string& rKey, int centIdx);
  TH3F* getOrBookJES3_jet1Pt_alphaHist (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // radius-tagged Jet13 + Balance3 bookers
  //   - no centrality suffix
  //   - name pattern: <base>_<rKey>
  // -------------------------------------------------------------------------
  TH3F*       getOrBookPho3TightIso       (const std::string& trig);               // photon baseline
  TH3F*       getOrBookJet13RecoilJet1    (const std::string& trig, const std::string& rKey);
  TProfile3D* getOrBookBalance3           (const std::string& trig, const std::string& rKey);

  // -------------------------------------------------------------------------
  // radius-tagged JES3 bookers (truth)
  //   - centrality suffix only (Au+Au)
  //   - name pattern: <base>_<rKey><centSuffix>
  // -------------------------------------------------------------------------
  // TRUTH reco-conditioned + reco↔truth jet1 match (your existing "works too well" truth)
  TH3F* getOrBookJES3Truth_xJ_alphaHist         (const std::string& trig, const std::string& rKey, int centIdx);

  // TRUTH reco-conditioned, but WITHOUT requiring reco jet1 ↔ truth jet1 ΔR match
  TH3F* getOrBookJES3TruthRecoCondNoJetMatch_xJ_alphaHist(const std::string& trig, const std::string& rKey, int centIdx);

  TH3F* getOrBookJES3Truth_jet1Pt_alphaHist     (const std::string& trig, const std::string& rKey, int centIdx);

  // PURE truth xJgamma distribution (no reco gating, no reco↔truth jet matching)
  TH3F* getOrBookJES3TruthPure_xJ_alphaHist     (const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // reco-side JES3 subsets to diagnose fakes / wrong-jet assignment
  // -------------------------------------------------------------------------
  // RECO xJ,alpha for events where the LEADING reco iso∧tight photon is truth-signal tagged (photon match only)
  TH3F* getOrBookJES3RecoTruthPhoTagged_xJ_alphaHist(const std::string& trig, const std::string& rKey, int centIdx);

  // RECO xJ,alpha for truth-tagged pairs (photon match + reco jet1 matched to truth jet1)
  TH3F* getOrBookJES3RecoTruthTagged_xJ_alphaHist(const std::string& trig, const std::string& rKey, int centIdx);

  // -------------------------------------------------------------------------
  // Jet QA: generic bookers + fillers (already radius-tagged in your .cc)
  // -------------------------------------------------------------------------
  TH1I* getOrBookJetQA1I(const std::string& trig,
                         const std::string& base,
                         const std::string& xAxisTitle,
                         const std::string& rKey,
                         int ptIdx, int centIdx,
                         int nbins, double xmin, double xmax);

  TH1F* getOrBookJetQA1F(const std::string& trig,
                         const std::string& base,
                         const std::string& xAxisTitle,
                         const std::string& rKey,
                         int ptIdx, int centIdx,
                         int nbins, double xmin, double xmax);

  TH2F* getOrBookJetQA2F(const std::string& trig,
                         const std::string& base,
                         const std::string& xAxisTitle,
                         const std::string& yAxisTitle,
                         const std::string& rKey,
                         int ptIdx, int centIdx,
                         int nxbins, double xmin, double xmax,
                         int nybins, double ymin, double ymax);

  void fillInclusiveJetQA(const std::vector<std::string>& activeTrig,
                          int centIdx,
                          const std::string& rKey);
  bool fillPPG12InclusiveJetTruthSpectrumQA(PHCompositeNode* topNode = nullptr);

  void fillSelectedJetQA(const std::vector<std::string>& activeTrig,
                         int ptIdx, int centIdx,
                         const std::string& rKey,
                         const Jet* jet1,
                         const Jet* jet2);

  // -------------------------------------------------------------------------
  // SS + Iso category accounting
  // -------------------------------------------------------------------------
    void processCandidatesForCurrentIsoView(PHCompositeNode* topNode,
                                            const std::vector<std::string>& activeTrig);
    bool initReplayFoundation();
    void writeReplayFoundationEvent(PHCompositeNode* topNode, int terminalStatus);
  void fillIsoSSTagCounters(const std::string& trig,
                            const RawCluster* clus,
                            const SSVars& v,
                            TightTag tightTag,
                            double pt_gamma,
                            int centIdx,
                            PHCompositeNode* topNode);

  void bumpHistFill(const std::string& trig, const std::string& hnameWithSuffix);
  void printCutSummary() const;

  // -------------------------------------------------------------------------
  // Member data
  // -------------------------------------------------------------------------
  std::string Outfile;
  TFile*      out    = nullptr;

  bool captureReplayTrigger(PHCompositeNode* topNode);
  RJGl1TriggerContract::Decisions m_replayGl1Decisions;
  std::uint64_t m_replayTriggerSnapshotId=0;

  // Dataset flags (these may be overridden at runtime by env in fetchNodes())
  bool m_isSim  = false;
  bool m_isAuAu = false;
  bool m_doPi0Analysis = false;

  // Vertex / event selection
  GlobalVertex* m_vtx = nullptr;
  float m_vx = 0.0f;
  float m_vy = 0.0f;
  float m_vz = 0.0f;
  double m_truthVz = std::numeric_limits<double>::quiet_NaN();
  double m_truthVzMB = std::numeric_limits<double>::quiet_NaN();

  bool  m_useVzCut = true;
  float m_vzCut    = 30.0f;

  bool        m_vertexReweightOn = false;
  std::string m_vertexReweightFile;
  std::string m_vertexReweightHist = "h_w_iterative";
  TH1*        m_vertexReweightH = nullptr;
  double      m_mcVertexWeight = 1.0;
  double      m_mcEventWeight = 1.0;
  bool        m_ppg12SimWeightActive = false;
  int         m_ppg12SimWeightLaneCode = 0;  // 1=photon+jet, 2=inclusive-jet
  int         m_ppg12SimWeightLaneComponentCode = 0;
  double      m_ppg12SimWeightFactorSlice = 1.0;
  double      m_ppg12SimWeightFactorVertex = 1.0;
  double      m_ppg12SimWeightFactorMix = 1.0;
  double      m_ppg12SimWeightFactorPeriod = 1.0;
  double      m_ppg12SimWeightFactorFinal = 1.0;

  // pp PhotonJet5/10/20 stitching diagnostics.  These are intentionally
  // parallel to the nominal event weight so a debug production can compare
  // stitch-variable choices and vertex weighting from the same output files.
  static constexpr std::size_t kPPStitchDiagVariantCount = 5;
  bool m_ppStitchDiagContext = false;
  std::array<bool, kPPStitchDiagVariantCount>   m_ppStitchDiagKeep{};
  std::array<double, kPPStitchDiagVariantCount> m_ppStitchDiagPhotonPt{};

  // Selection-neutral source witnesses for direct PPG12 diagnostics that
  // cannot be reconstructed from the candidate/jet rows alone.  These are
  // typed scientific occurrences, never histogram names or bin coordinates.
  struct ReplayPPG12DiagnosticOccurrence
  {
    std::string trigger;
    std::string occurrence_kind;
    int stage_code = 0;
    int source_code = 0;
    int sample_code = 0;
    int decision = 0;
    int cluster_ordinal = -1;
    int truth_track_id = -1;
    int generator_barcode = -1;
    double source_photon_pt = std::numeric_limits<double>::quiet_NaN();
    double flow_value = std::numeric_limits<double>::quiet_NaN();
    double reco_photon_pt = std::numeric_limits<double>::quiet_NaN();
    double response_photon_pt = std::numeric_limits<double>::quiet_NaN();
    double truth_photon_pt = std::numeric_limits<double>::quiet_NaN();
    double truth_prior_weight = std::numeric_limits<double>::quiet_NaN();
    double reco_vertex_z = std::numeric_limits<double>::quiet_NaN();
    double hard_truth_vertex_z = std::numeric_limits<double>::quiet_NaN();
    double mb_truth_vertex_z = std::numeric_limits<double>::quiet_NaN();
    double vertex_weight = std::numeric_limits<double>::quiet_NaN();
    double period_event_weight = std::numeric_limits<double>::quiet_NaN();
    double cluster_energy = std::numeric_limits<double>::quiet_NaN();
    double cluster_eta = std::numeric_limits<double>::quiet_NaN();
    double cluster_et = std::numeric_limits<double>::quiet_NaN();
    double truth_energy = std::numeric_limits<double>::quiet_NaN();
    double truth_eta = std::numeric_limits<double>::quiet_NaN();
    double truth_et = std::numeric_limits<double>::quiet_NaN();
    double energy_contribution = std::numeric_limits<double>::quiet_NaN();
    double occurrence_weight = 1.0;
    std::uint64_t selection_bitmask = 0;
  };
  std::vector<ReplayPPG12DiagnosticOccurrence> m_replayPPG12DiagnosticOccurrences;
  std::array<unsigned long long, 5> m_replayFullFillCountsAtEventStart{};
  std::array<unsigned long long, 5> m_replayRawFillCountsAtEventStart{};

  // Centrality
  int m_centBin = -1;                 // 0..99 (Au+Au), or -1 in pp
  std::vector<int> m_centEdges;       // centrality bin edges, e.g. {0,10,20,...,100}

  // Photon fiducial + binning
  double m_etaAbsMax = 0.7;           // photon |eta| cut

  std::vector<double> m_gammaPtBins = {15,17,19,21,23,26,35};  // canonical photon pT bin edges (6 bins, start at 15 GeV)

    // Photon ID cuts (PPG12 Table 4) defaults (match PhoIDCuts namespace baseline)
    double m_phoid_pre_e11e33_max = 0.98;
    double m_phoid_pre_et1_min    = 0.60;
    double m_phoid_pre_et1_max    = 1.00;
    double m_phoid_pre_e32e35_min = 0.80;
    double m_phoid_pre_e32e35_max = 1.00;
    double m_phoid_pre_weta_max   = 0.60;

    std::string m_preselectionVariant = "reference";
    std::string m_tightVariant = "reference";
    std::string m_nonTightVariant = "reference";
    std::string m_preselectionPhotonNode = "PHOTONCLUSTER_CEMC";
    std::string m_tightPhotonNode = "PHOTONCLUSTER_CEMC";
    bool m_explicitPhotonIDVariants = false;
    double m_npbCut = 0.5;
    double m_tightBDTMinIntercept = 0.815625;
    double m_tightBDTMinSlope = -0.0015625;
    double m_tightBDTMax = 1.0;
    double m_nonTightBDTMinIntercept = 0.7333333333333333;
    double m_nonTightBDTMinSlope = -0.01333333333333333;
    double m_nonTightBDTMaxIntercept = 0.684375;
    double m_nonTightBDTMaxSlope = 0.0015625;

    double m_phoid_tight_w_lo           = 0.0;
    double m_phoid_tight_w_hi_intercept = 0.15;
    double m_phoid_tight_w_hi_slope     = 0.006;

    double m_phoid_tight_e11e33_min = 0.40;
    double m_phoid_tight_e11e33_max = 0.98;

    double m_phoid_tight_et1_min    = 0.90;
    double m_phoid_tight_et1_max    = 1.00;

    double m_phoid_tight_e32e35_min = 0.92;
    double m_phoid_tight_e32e35_max = 1.00;


  // Phase-1 YAML knobs (matching thresholds)
  double m_phoMatchDRMax = 0.05;
  double m_jetMatchDRMax = 0.3;

  // Phase-1 YAML knobs (explicit unfolding bin edges)
  std::vector<double> m_unfoldRecoPhotonPtBins  = {10,15,17,19,21,23,26,35,40};
  std::vector<double> m_unfoldTruthPhotonPtBins = {5,10,15,17,19,21,23,26,35,40};
  std::vector<double> m_unfoldJetPtBins;
  std::vector<double> m_unfoldXJBins = {0.0,0.20,0.24,0.29,0.35,0.41,0.50,0.60,0.72,0.86,1.03,1.24,1.49,1.78,2.14,3.0};
  std::string m_leadingResponseFamilyLabel = "";

  // Analysis provenance stamping (written once into output ROOT by RecoilJets.cc)
  std::string m_analysisConfigYAMLText = "";
  std::string m_analysisConfigTag      = "";
  bool        m_analysisConfigStamped  = false;


  // Generic/default isolation WP. PPG12 photon-yield pp mode requires:
  // sliding R=0.4, Eiso < 0.490 + 0.037*pT, non-iso gap=0.8.
  double m_isoA      = 0.490;
  double m_isoB      = 0.037;
  double m_isoGap    = 0.8;
  double m_isoFixed  = 2.0;          // used ONLY when m_isSlidingIso==false (RECO)
  double m_truthIsoMaxGeV = 4.0;     // TRUTH isolation max (independent of sliding/fixed mode)
  double m_isoConeR  = 0.3;
  double m_isoTowMin = 0.0;
  bool   m_isSlidingIso = true;
  struct PPIsoWP
  {
      double aGeV = 0.0;
      double bPerGeV = 0.0;
      double sideGapGeV = 1.0;
      bool enabled = false;
  };
  PPIsoWP m_ppIsoWPR30;
  PPIsoWP m_ppIsoWPR40;
  std::vector<IsoView> m_internalIsoViews;
  std::string m_activeIsoViewSuffix;

  // Jet selection WP
  //  double m_minJetPt      = 5.0;
  double m_minJetPt      = 10.0;
  std::vector<double> m_internalJetPtCuts;
  double m_internalJetPtBaseKeyCut = -1.0;
  double m_minBackToBack = 7.0 * M_PI / 8.0;    // radians
  std::vector<double> m_internalBackToBackCuts;
  double m_internalBackToBackBaseKeyCut = -1.0;

  // Legacy "primary" reco jet key (still used for printing/overrides only)
  std::string m_xjRecoJetKey = "r04";

  // Nodes: truth / photons / clusters
  PHG4TruthInfoContainer* m_truthInfo      = nullptr;
  RawClusterContainer* m_clus              = nullptr;
  RawClusterContainer* m_clus_nocorr       = nullptr;
  RawClusterContainer* m_photons           = nullptr;
  RawClusterContainer* m_photons_npb       = nullptr;
  RawClusterContainer* m_photons_tightbdt  = nullptr;
  RawClusterContainer* m_ppg12TopoClusters = nullptr;
  RawClusterContainer* m_ppg12Fig8Clusters = nullptr;
  MbdOut* m_mbdout                         = nullptr;

  // Calo tower bundles (node cache)
  struct CaloBundle
  {
    TowerInfoContainer*   towers = nullptr;
    RawTowerGeomContainer* geom  = nullptr;
    double                sumEt  = 0.0;
  };

  std::vector<std::tuple<std::string, std::string, std::string>> m_caloInfo;
  std::map<std::string, CaloBundle> m_calo;

    // -------------------------------------------------------------------------
    // parallel jet containers by radius key
    // -------------------------------------------------------------------------
    std::map<std::string, JetContainer*> m_jets;
    std::map<std::string, JetContainer*> m_jetsRaw;
    // Au+Au-only companion view built from the retowered, un-subtracted
    // calorimeter towers.  The canonical m_jets/m_jetsRaw maps continue to
    // own the SUB1 view used by the established analysis.
    std::map<std::string, JetContainer*> m_jetsNoSub;
    std::map<std::string, JetContainer*> m_jetsNoSubRaw;

  // Optional: limit active jet radii (keys like "r02","r04"). Empty => all kJetRadii.
  std::vector<std::string> m_activeJetRKeys;

  // -------------------------------------------------------------------------
  // NEW: truth jet containers by radius key (SIM only)
  // -------------------------------------------------------------------------
  std::map<std::string, JetContainer*> m_truthJetsByRKey;
  std::map<std::string, std::string>   m_truthJetsNodeByRKey;

  // -------------------------------------------------------------------------
  // EventDisplay diagnostics payload support nodes (optional; never affects physics)
  // -------------------------------------------------------------------------
  EventHeader*           m_evtHeader          = nullptr;
  TowerInfoContainer*    m_evtDispTowersCEMC  = nullptr;
  TowerInfoContainer*    m_evtDispTowersIHCal = nullptr;
  TowerInfoContainer*    m_evtDispTowersOHCal = nullptr;
  RawTowerGeomContainer* m_evtDispGeomCEMC    = nullptr;
  RawTowerGeomContainer* m_evtDispGeomIHCal   = nullptr;
  RawTowerGeomContainer* m_evtDispGeomOHCal   = nullptr;

  // -------------------------------------------------------------------------
  // EventDisplay diagnostics payload (offline rendering; independent of Verbosity()).
  //
  //  - One TTree entry per (event, rKey) when enabled.
  //  - Stores compact "ingredients" (jet kinematics + sparse tower constituent lists)
  //    so that single-event displays can be rendered OFFLINE without any online scanning
  //    or Verbosity() hacks.
  // -------------------------------------------------------------------------
  bool m_evtDiagEnabled    = false;  // user-controlled switch (macro setter or env)
  bool m_evtDiagNodesReady = false;  // per-event: true if required tower/geom nodes are available

  int  m_evtDiagMaxPerBin = 0;       // 0 => unlimited; else limit entries per (rKey,ptBin,cat)
  std::unordered_map<std::string, int> m_evtDiagSavedPerBin;

  long long m_evtDiagNFill = 0;
  long long m_evtDiagNFillWithSelTowers = 0;
  long long m_evtDiagNFillWithBestTowers = 0;
  long long m_evtDiagNFillWithAnyTowers = 0;
  long long m_evtDiagNFillByCat[3] = {0, 0, 0};
  long long m_evtDiagNFillWithAnyTowersByCat[3] = {0, 0, 0};

  TTree* m_evtDiagTree = nullptr;

  int       m_evtDiag_run        = 0;
  int       m_evtDiag_evt        = 0;
  long long m_evtDiag_eventCount = 0;
  float     m_evtDiag_vz         = 0.0f;

    std::string m_evtDiag_rKey;
    int         m_evtDiag_ptBin = -1;
    int         m_evtDiag_cat   = -1;
    int         m_evtDiag_isSim = 0;

    float m_evtDiag_ptGammaTruth  = 0.0f;
    float m_evtDiag_phiGammaTruth = 0.0f;
    float m_evtDiag_ptGammaReco   = 0.0f;
    float m_evtDiag_etaGammaReco  = 0.0f;
    float m_evtDiag_phiGammaReco  = 0.0f;

  float m_evtDiag_sel_pt  = 0.0f;
  float m_evtDiag_sel_eta = 0.0f;
  float m_evtDiag_sel_phi = 0.0f;

  float m_evtDiag_best_pt  = 0.0f;
  float m_evtDiag_best_eta = 0.0f;
  float m_evtDiag_best_phi = 0.0f;

  float m_evtDiag_truthLead_pt  = 0.0f;
  float m_evtDiag_truthLead_eta = 0.0f;
  float m_evtDiag_truthLead_phi = 0.0f;

  float m_evtDiag_drSelToTruthLead  = -1.0f;
  float m_evtDiag_drBestToTruthLead = -1.0f;

  std::vector<int>   m_evtDiag_sel_calo;
  std::vector<int>   m_evtDiag_sel_ieta;
  std::vector<int>   m_evtDiag_sel_iphi;
  std::vector<float> m_evtDiag_sel_etaTower;
  std::vector<float> m_evtDiag_sel_phiTower;
  std::vector<float> m_evtDiag_sel_etTower;
  std::vector<float> m_evtDiag_sel_eTower;

  std::vector<int>   m_evtDiag_best_calo;
  std::vector<int>   m_evtDiag_best_ieta;
  std::vector<int>   m_evtDiag_best_iphi;
  std::vector<float> m_evtDiag_best_etaTower;
  std::vector<float> m_evtDiag_best_phiTower;
  std::vector<float> m_evtDiag_best_etTower;
  std::vector<float> m_evtDiag_best_eTower;

  // Compatibility photon-ID training tree for pp ML studies.
  // The tree name intentionally matches the Au+Au trainer contract.
  bool m_ppPhotonIDTrainingTreeEnabled = false;
  bool m_ppPhotonIDExtractOnly = false;
  bool m_ppPhotonIDPPG12Filter = true;
  bool m_ppPhotonIDRequirePreselection = false;
  long long m_ppPhotonIDTrainingTreeMaxEntries = 0;
  long long m_ppPhotonIDTrainingTreeEntries = 0;
  TTree* m_ppPhotonIDTrainingTree = nullptr;
  std::string m_ppPhotonIDSourceRole = "auto";  // auto, signal, background, all

  bool m_ppg12TableQAEnabled = false;
  // THE-119 canary-only, selection-neutral witness for the loose replay
  // capture population. It never sets a tag or source-ownership decision.
  bool m_replayFoundationCaptureWitnessEnabled = false;
  bool m_ppg12TableQANPBDataTaggingEnabled = false;
  bool m_ppg12Fig7TriggerDiagnostic = false;
  bool m_ppg12Fig11SBDiagnostic = false;
  bool m_ppg12Fig13Bit30Diagnostic = false;
  bool m_ppg12Fig13ParityQA = false;
  double m_ppg12TableQANPBTagTimeSampleNs = 17.6;
  double m_ppg12TableQANPBDeltaTCut = -5.0;
  double m_ppg12TableQANPBWetaMin = 0.4;
  double m_ppg12TableQANPBAwayJetPtMin = 5.0;
  double m_ppg12TableQANPBAwayJetDPhiMin = 1.5707963267948966;
  double m_ppg12Fig14NPBMbdWetaMin = 0.6;
  double m_ppg12TableQAMcIsoScale = 1.2;
  double m_ppg12TableQAMcIsoShift = 0.2;
  std::string m_ppg12TableQAMbdT0CorrectionFile;
  std::map<int, double> m_ppg12TableQAMbdT0Correction;
  std::set<std::string> m_ppg12TableQASchemaBookedTriggers;
  bool m_ppg12PhotonYieldEnabled = false;
  bool m_ppg12PhotonYieldUseTopoIso = true;
  bool m_ppg12PhotonYieldApplyTowerMask = true;
  bool m_ppg12PhotonYieldExcludeCandidateTopoCluster = false;
  bool m_ppg12PhotonYieldDoubleInteraction = false;
  bool m_ppg12PhotonYieldDiagFeatures = false;
  bool m_ppg12PhotonYieldApplyBinning = false;
  bool m_ppg12Fig8UseFallbackClusterNode = true;
  bool m_ppg12PeriodContractEnabled = false;
  bool m_ppg12PeriodUseLumiWeight = true;
  bool m_ppg12PeriodFilterData = false;
  bool m_ppg12PeriodStrictDoubleMB = true;
  bool m_ppg12PeriodMixWeightAuto = false;
  bool m_ppg12PeriodVertexFileAuto = false;
  std::string m_ppg12PhotonYieldTowerMaskFile;
  std::string m_ppg12PhotonYieldTowerMaskName;
  std::string m_ppg12Fig8ClusterNode = "CLUSTERINFO_CEMC_NO_SPLIT";
  std::string m_ppg12PeriodKey = "unset";
  std::string m_ppg12PeriodLabel = "unset";
  std::string m_ppg12PeriodExpectedVertexFile;
  TH2* m_ppg12PhotonYieldTowerMask = nullptr;
  double m_ppg12PhotonYieldMcIsoScale = 1.2;
  double m_ppg12PhotonYieldMcIsoShift = 0.1;
  double m_ppg12PhotonYieldMixWeight = 1.0;
  int m_ppg12PeriodRunMin = -1;
  int m_ppg12PeriodRunMaxExclusive = -1;
  double m_ppg12PeriodLumi = 1.0;
  double m_ppg12PeriodLumiTarget = 1.0;
  double m_ppg12PeriodLumiWeight = 1.0;
  double m_ppg12PeriodFDouble = -1.0;
  double m_ppg12PeriodFSingle = -1.0;
  double m_ppg12PeriodClosureZCut = -1.0;
  bool m_ppg12Fig6EventCanaryEnabled = false;
  bool m_ppg12Fig6FilledThisEvent = false;
  long long m_ppg12Fig6EventCanaryMaxRows = 6000;
  long long m_ppg12Fig6EventCanaryRowsWritten = 0;
  std::string m_ppg12Fig6EventCanaryPath;
  std::string m_ppg12Fig6EventCanarySummaryPath;
  std::ofstream m_ppg12Fig6EventCanaryOut;
  bool m_ppg12EisoPathCanaryEnabled = false;
  long long m_ppg12EisoPathCanaryMaxRows = 5000;
  long long m_ppg12EisoPathCanaryRowsWritten = 0;
  double m_ppg12EisoPathCanaryMinEt = 10.0;
  double m_ppg12EisoPathCanaryMaxEt = 36.0;
  std::string m_ppg12EisoPathCanaryPath;
  std::ofstream m_ppg12EisoPathCanaryOut;
  long long m_ppg12EisoConeMemberCanaryMaxRows = 100000;
  long long m_ppg12EisoConeMemberCanaryRowsWritten = 0;
  std::string m_ppg12EisoConeMemberCanaryPath;
  std::ofstream m_ppg12EisoConeMemberCanaryOut;
  std::string m_ppg12EisoVertexSource = "unset";
  long long m_ppg12Fig6CanaryEventsSeen = 0;
  long long m_ppg12Fig6CanaryValidTruthJets = 0;
  long long m_ppg12Fig6CanaryOwnedWindowEvents = 0;
  long long m_ppg12Fig6CanaryRawOwnedEntries = 0;
  double m_ppg12Fig6CanaryWeightedOwnedIntegral = 0.0;
  double m_ppg12Fig6CanarySumXsecFactor = 0.0;
  double m_ppg12Fig6CanarySumLumiFactor = 0.0;
  double m_ppg12Fig6CanarySumMixFactor = 0.0;
  double m_ppg12Fig6CanarySumVertexWeight = 0.0;
  double m_ppg12Fig6CanarySumFinalFillWeight = 0.0;
  double m_ppg12Fig6CanaryIANSuffixIntegral = std::numeric_limits<double>::quiet_NaN();
  double m_ppg12Fig6CanaryTruePeriodIntegral = std::numeric_limits<double>::quiet_NaN();
  std::set<std::string> m_ppg12PhotonYieldSchemaBookedTriggers;
  std::vector<double> m_ppg12PhotonYieldRecoPtBins =
      {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};
  std::vector<double> m_ppg12PhotonYieldTruthPtBins =
      {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};

  int m_bdtTrain_run = 0;
  long long m_bdtTrain_evt = 0;
  long long m_bdtTrain_eventnumber = 0;
  int m_bdtTrain_cluster_index = -1;
  int m_bdtTrain_is_signal = 0;
  int m_bdtTrain_pt_bin = -1;
  int m_bdtTrain_cent_bin = -1;
  float m_bdtTrain_pt = 0.0f;
  float m_bdtTrain_score_input_et = 0.0f;
  float m_bdtTrain_response_et = 0.0f;
  float m_bdtTrain_eta = 0.0f;
  float m_bdtTrain_phi = 0.0f;
  float m_bdtTrain_cent = -1.0f;
  float m_bdtTrain_vz = 0.0f;
  float m_bdtTrain_ppg12_kin_vertexz = 0.0f;
  float m_bdtTrain_weight = 1.0f;
  float m_bdtTrain_eiso = 0.0f;
  float m_bdtTrain_ppg12_raw_eiso = 0.0f;
  float m_bdtTrain_ppg12_reco_eiso = 0.0f;
  float m_bdtTrain_ppg12_iso_threshold = 0.0f;
  float m_bdtTrain_ppg12_noniso_threshold = 0.0f;
  int m_bdtTrain_ppg12_is_iso = 0;
  int m_bdtTrain_ppg12_is_noniso = 0;
  int m_bdtTrain_ppg12_common_pass = 0;
  int m_bdtTrain_ppg12_tight_tag = -1;
  int m_bdtTrain_ppg12_logical_abcd_region = 0;
  int m_bdtTrain_ppg12_analysis_window_pass = -1;
  int m_bdtTrain_ppg12_response_window_pass = -1;
  int m_bdtTrain_ppg12_signal_fill_a = 0;
  int m_bdtTrain_ppg12_signal_fill_b = 0;
  int m_bdtTrain_ppg12_signal_fill_c = 0;
  int m_bdtTrain_ppg12_signal_fill_d = 0;
  int m_bdtTrain_ppg12_signal_fill_multiplicity = 0;
  int m_bdtTrain_ppg12_sample_bin = 0;
  float m_bdtTrain_ppg12_xsec_pb = -999.0f;
  float m_bdtTrain_ppg12_xsec_weight = 1.0f;
  float m_bdtTrain_ppg12_window_low = -999.0f;
  float m_bdtTrain_ppg12_window_high = -999.0f;
  float m_bdtTrain_max_truth_jet_pt_r04 = -999.0f;
  int m_bdtTrain_ppg12_truth_window_pass_r04 = -1;
  int m_bdtTrain_truth_track_id = -1;
  int m_bdtTrain_truth_barcode = -1;
  int m_bdtTrain_truth_class = -999;
  float m_bdtTrain_truth_energy_contribution = -999.0f;
  int m_bdtTrain_ppg12_weight_lane_code = 0;
  int m_bdtTrain_ppg12_weight_component_code = 0;
  float m_bdtTrain_ppg12_weight_slice = 1.0f;
  float m_bdtTrain_ppg12_weight_vertex = 1.0f;
  float m_bdtTrain_ppg12_weight_mix = 1.0f;
  float m_bdtTrain_ppg12_weight_period = 1.0f;
  float m_bdtTrain_ppg12_weight_final = 1.0f;
  float m_bdtTrain_weta = 0.0f;
  float m_bdtTrain_wphi = 0.0f;
  float m_bdtTrain_weta33 = 0.0f;
  float m_bdtTrain_wphi33 = 0.0f;
  float m_bdtTrain_weta35 = 0.0f;
  float m_bdtTrain_wphi53 = 0.0f;
  float m_bdtTrain_et1 = 0.0f;
  float m_bdtTrain_et2 = 0.0f;
  float m_bdtTrain_et3 = 0.0f;
  float m_bdtTrain_et4 = 0.0f;
  float m_bdtTrain_e11e33 = 0.0f;
  float m_bdtTrain_e32e35 = 0.0f;
  float m_bdtTrain_e11e22 = 0.0f;
  float m_bdtTrain_e11e13 = 0.0f;
  float m_bdtTrain_e11e15 = 0.0f;
  float m_bdtTrain_e11e17 = 0.0f;
  float m_bdtTrain_e11e31 = 0.0f;
  float m_bdtTrain_e11e51 = 0.0f;
  float m_bdtTrain_e11e71 = 0.0f;
  float m_bdtTrain_e22e33 = 0.0f;
  float m_bdtTrain_e22e35 = 0.0f;
  float m_bdtTrain_e22e37 = 0.0f;
  float m_bdtTrain_e22e53 = 0.0f;
  float m_bdtTrain_w32 = 0.0f;
  float m_bdtTrain_w52 = 0.0f;
  float m_bdtTrain_w72 = 0.0f;
  float m_bdtTrain_cluster_prob = -2.0f;
  float m_bdtTrain_npb_score = -2.0f;
  float m_bdtTrain_tight_bdt_score = -2.0f;
  float m_bdtTrain_ppg12_shape_n_owned = -999.0f;
  float m_bdtTrain_ppg12_shape_owned_e = -999.0f;
  float m_bdtTrain_ppg12_shape_all_e = -999.0f;
  float m_bdtTrain_ppg12_shape_den_e = -999.0f;
  float m_bdtTrain_ppg12_shape_weta_cogx_num = -999.0f;
  float m_bdtTrain_ppg12_shape_wphi_cogx_num = -999.0f;
  float m_bdtTrain_ppg12_shape_cog_eta = -999.0f;
  float m_bdtTrain_ppg12_shape_cog_phi = -999.0f;
  float m_bdtTrain_ppg12_shape_center_ieta = -999.0f;
  float m_bdtTrain_ppg12_shape_center_iphi = -999.0f;

  // -------------------------------------------------------------------------
  // Diagnostics / accounting
  // -------------------------------------------------------------------------
  EventReject m_lastReject = EventReject::None;

  long long event_count  = 0;
  long long m_evtNoTrig  = 0;

  Bookkeeping m_bk{};

  // For End() histogram summary (counts how many times each histogram was filled)
  std::unordered_map<std::string, long long> m_histFill;

  // -------------------------------------------------------------------------
  // NEW: per-pT-bin negative-isolation bookkeeping (independent of SS classification)
  //
  //  - "Builder"    = RecoilJets::eiso() using PhotonClusterBuilder iso_* layer sums
  //  - "ClusterIso" = RawCluster::get_et_iso(radiusx10,false,true) (UNSUBTRACTED)
  // -------------------------------------------------------------------------
  std::vector<unsigned long long> m_nIsoBuilderByPt;
  std::vector<unsigned long long> m_nIsoBuilderNegByPt;
  std::vector<unsigned long long> m_nIsoClusterIsoByPt;
  std::vector<unsigned long long> m_nIsoClusterIsoNegByPt;

  // Histogram storage: trigger -> (name -> object*)
  std::map<std::string, HistMap> qaHistogramsByTrigger;

  // Per-trigger slice counters printed in End()
  std::map<std::string, std::map<std::string, CatStat>> m_catByTrig;

  bool m_replayFoundationEnabled = false;
  bool m_replayNodesReady = false;
  bool m_collaboratorRecoVertexValid = false;
  bool m_replayWriteFailed = false;
  bool m_the134MultiviewSidecarOnly = false;
  bool m_the134FastExtraction = false;
  std::unique_ptr<RJReplayRuntimeV1::Runtime> m_replayRuntime;
  std::unique_ptr<RJTriggerScalersV1::Writer> m_triggerScalerWriter;
  std::unique_ptr<RJTriggerRunInfoV1::Writer> m_triggerRunInfoWriter;
  std::unique_ptr<RJPhotonTrainingViewV1::Runtime> m_photonTrainingViewRuntime;
};

#endif // RECOILJETS_H
