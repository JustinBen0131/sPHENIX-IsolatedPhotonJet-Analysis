// Tell Emacs this is C++   -*- C++ -*-
#ifndef RECOILJETS_AuAu_H
#define RECOILJETS_AuAu_H

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
#include <TH1F.h>
#include <TH1I.h>
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

// Replay schema types used by value in the class definition below must be
// complete before this header is parsed independently.  Use the installed
// public-header name: a relative source-tree include becomes invalid after
// makeProject installs RecoilJets_AuAu.h under include/caloana/.
#include "RJReplayFoundationV1.h"

// ---------------------------------------------------------------------------
// Helper alias: maps histogram name -> ROOT object*
// ---------------------------------------------------------------------------
using HistMap = std::map<std::string, TObject*>;

// Forward declarations (pointers only in the interface)
class CentralityInfo;
class Gl1Packet;
class EventHeader;
class PHG4TruthInfoContainer;
namespace TMVA { namespace Experimental { class RBDT; } }
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
class MbdPmtContainer;
class MbdOut;
class MbdGeom;
class EpdGeom;

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
//       auauOnlyNPB = AuAu-trained NPB score only
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
        kNonTight         = 2,  // >=2 fails
        kNeither          = 3   // exactly 1 fail
    };
    
    enum class EventReject : std::uint8_t
    {
        None       = 0,
        Trigger    = 1,
        Vz         = 2,
        MinBias    = 3,
        Stitch     = 4,
        Centrality = 5
    };
    
    // Shower-shape variables extracted from PhotonClusterv1
    struct SSVars
    {
        double pt_gamma       = 0.0;
        double weta_cogx      = 0.0;
        double wphi_cogx      = 0.0;
        double weta33_cogx    = 0.0;
        double wphi33_cogx    = 0.0;
        double weta35_cogx    = 0.0;
        double wphi53_cogx    = 0.0;
        double et1            = 0.0;
        double et2            = 0.0;
        double et3            = 0.0;
        double et4            = 0.0;
        double e11_over_e33   = 0.0;
        double e17_over_e77   = std::numeric_limits<double>::quiet_NaN();
        double e32_over_e35   = 0.0;
        double e11_over_e22   = 0.0;
        double e11_over_e13   = 0.0;
        double e11_over_e15   = 0.0;
        double e11_over_e17   = 0.0;
        double e11_over_e31   = 0.0;
        double e11_over_e51   = 0.0;
        double e11_over_e71   = 0.0;
        double e22_over_e33   = 0.0;
        double e22_over_e35   = 0.0;
        double e22_over_e37   = 0.0;
        double e22_over_e53   = 0.0;
        double w32            = 0.0;
        double w52            = 0.0;
        double w72            = 0.0;
        double mean_time      = std::numeric_limits<double>::quiet_NaN();
        double npb_score      = std::numeric_limits<double>::quiet_NaN();
        double tight_bdt_score = std::numeric_limits<double>::quiet_NaN();
        double auau_npb_score = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_bdt_score = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_mlp_score = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_bdt_mlp_score = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_logreg_score = std::numeric_limits<double>::quiet_NaN();
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
        long long evt_fail_minbias = 0;
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
    
    // NOTE: bit 18 changed meaning across the Run-3 Au+Au GRL.
    //   Epoch 1 (runs 67599–68155, 21 runs): bit 18 = "Photon 10 GeV + MBD NS >= 2" (no vtx cut)
    //   Epoch 2–6 (runs 68208–78954, 1065 runs): bit 18 = "Photon 10 GeV + MBD NS >= 2, vtx < 10 cm"
    // The photon_10_plus_MBD_NS_geq_2_vtx_lt_10 folder will therefore include some events
    // with |vtx| > 10 cm from those 21 epoch-1 runs. Given 21/1086 runs the contamination
    // is negligible for a first look. If needed, add a software vtx < 10 cm cut in
    // RecoilJets_AuAu.cc specifically for the bit-18 path.
    inline static const std::vector<std::pair<int, std::string>> triggerNameMapAuAu = {
        //        {10, "MBD_NS_geq_2"},
        //        {11, "MBD_NS_geq_1"},
        //        {12, "MBD_NS_geq_2_vtx_lt_10"},
        //        {13, "MBD_NS_geq_2_vtx_lt_30"},
        {14, "MBD_NS_geq_2_vtx_lt_150"},
        //        {15, "MBD_NS_geq_1_vtx_lt_10"},
        //        {16, "photon_6_plus_MBD_NS_geq_2_vtx_lt_10"},
        //        {17, "photon_8_plus_MBD_NS_geq_2_vtx_lt_10"},
        //                {18, "photon_10_plus_MBD_NS_geq_2_vtx_lt_10"},
        //        {19, "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"},
        //        {20, "photon_6_plus_MBD_NS_geq_2_vtx_lt_150"},
        //        {21, "photon_8_plus_MBD_NS_geq_2_vtx_lt_150"},
        {22, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"},
        {23, "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"}
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
        //   "isSim"/"sim"                 -> sim
        //   "isSimEmbedded"/"simembedded" -> sim input + AuAu-like reconstruction
        //   "isSimEmbeddedInclusive"      -> inclusive sim input + AuAu-like reconstruction
        //   "isAuAu"/"auau"               -> Au+Au
        //   "ispp"/"pp"                   -> pp
        const std::string t = s;
        if (t == "isSimEmbedded" || t == "simembedded" ||
            t == "isSimEmbeddedInclusive" || t == "simembeddedinclusive") { m_isSim = true;  m_isAuAu = true;  m_isSimEmbedded = true;  }
        if (t == "isSim" || t == "sim")                 { m_isSim = true;  m_isAuAu = false; m_isSimEmbedded = false; }
        if (t == "isAuAu" || t == "auau"|| t == "aa")  { m_isSim = false; m_isAuAu = true;  m_isSimEmbedded = false; }
        if (t == "ispp"  || t == "pp")                 { m_isSim = false; m_isAuAu = false; m_isSimEmbedded = false; }
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
    void setMinBiasClassifier(bool on) { m_useMinBiasClassifier = on; }
    // pp-only unified-macro compatibility; AuAu vertex handling is unchanged.
    void setPPG12PhotonYieldUseTruthVertexInPPSim(bool = true) {}
    void setPPG12PhotonYieldUseTruthVertexForRecoObjectsInPPSim(bool = true) {}
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
    void setRequireTowerInfoTruthMatching(bool on = true)            { m_requireTowerInfoTruthMatching = on; }
    
    // Analysis provenance stamping (written once into output ROOT by RecoilJets.cc)
    void setAnalysisConfigYAML(const std::string& yamlText, const std::string& tag)
    {
        m_analysisConfigYAMLText = yamlText;
        m_analysisConfigTag      = tag;
    }
    
    void enablePi0Analysis(bool on = true) { m_doPi0Analysis = on; }

    void setCentEdges(const std::vector<int>& edges)     { m_centEdges = edges; }
    std::uint64_t validCentralityObservedEvents() const { return m_validCentralityObservedEvents; }
    std::uint64_t invalidCentralityObservedEvents() const { return m_invalidCentralityObservedEvents; }
    std::uint64_t staleCentralityGuardedEvents() const { return m_staleCentralityGuardedEvents; }

    void setVertexReweighting(bool on,
                              const std::string& filePath,
                              const std::string& histPath = "data_over_MC_ratios/h_zvtx_ratio_data_over_photonJet")
    {
        m_vertexReweightOn   = on;
        m_vertexReweightFile = filePath;
        m_vertexReweightHist = histPath;
    }

    void setCentralityReweighting(bool on,
                                  const std::string& filePath,
                                  const std::string& histPath = "nom_cent_rw_hist")
    {
        m_centralityReweightOn   = on;
        m_centralityReweightFile = filePath;
        m_centralityReweightHist = histPath;
    }
    
    // Isolation WP (implemented in .cc)
    void setIsolationWP(double aGeV, double bPerGeV,
                        double sideGapGeV, double coneR, double towerMin,
                        double fixedGeV = 2.0);
    void setIsSlidingIso(bool on) { m_isSlidingIso = on; }
    struct IsoView
    {
        std::string label;
        double coneR = 0.3;
        bool isSliding = false;
        double fixedGeV = 4.0;
    };
    void setInternalIsoViews(const std::vector<IsoView>& views);
    
    // AuAu / embedded-SIM isolation WPs.
    // If m_isSlidingIso is true and exactly ONE entry is provided, RecoilJets_AuAu
    // interprets it as an event-centrality fit:
    //   thrReco(cent) = aGeV + bPerGeV * centrality
    // and applies that event's threshold as the reco isolation cut.
    // If multiple entries are provided, the entry indexed by findCentBin() is used.
    struct CentIsoWP { double aGeV; double bPerGeV; double sideGapGeV; };
    void setCentIsoWPs(const std::vector<CentIsoWP>& wps) { m_centIsoWPs = wps; }
    void setCentIsoWPsForCone(double coneR, const std::vector<CentIsoWP>& wps);
    
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

    void setAuAuTightBDTRuntimeConfig(const std::string& modelFile,
                                      const std::vector<std::string>& features,
                                      const std::vector<int>& centDepEdges = {},
                                      const std::vector<std::string>& centDepModelFiles = {},
                                      const std::vector<double>& ptBinEdges = {},
                                      const std::vector<std::string>& ptBinModelFiles = {},
                                      const std::vector<std::string>& ptCentModelFiles = {},
                                      const std::string& ptFallbackModelFile = "",
                                      const std::vector<std::string>& ptFallbackCentModelFiles = {},
                                      double ptFallbackMin = 35.0,
                                      double ptFallbackMax = 40.0,
                                      double applyPtMin = std::numeric_limits<double>::quiet_NaN(),
                                      double applyPtMax = std::numeric_limits<double>::quiet_NaN(),
                                      const std::vector<std::string>& workingPointEntries = {},
                                      const std::string& scoringMode = "");

    void setAuAuTightMLPRuntimeConfig(const std::string& modelFile,
                                      double minIntercept = 0.80,
                                      double minSlope = 0.0,
                                      double maxScore = 1.0,
                                      double nonTightMinIntercept = 0.20,
                                      double nonTightMinSlope = 0.0,
                                      double nonTightMaxIntercept = 0.80,
                                      double nonTightMaxSlope = 0.0,
                                      double applyPtMin = std::numeric_limits<double>::quiet_NaN(),
                                      double applyPtMax = std::numeric_limits<double>::quiet_NaN(),
                                      const std::vector<std::string>& workingPointEntries = {},
                                      const std::string& scoringMode = "");

    void setAuAuTightBDTMLPStackRuntimeConfig(const std::string& modelFile,
                                              double minIntercept = 0.80,
                                              double minSlope = 0.0,
                                              double maxScore = 1.0,
                                              double nonTightMinIntercept = 0.20,
                                              double nonTightMinSlope = 0.0,
                                              double nonTightMaxIntercept = 0.80,
                                              double nonTightMaxSlope = 0.0,
                                              double applyPtMin = std::numeric_limits<double>::quiet_NaN(),
                                              double applyPtMax = std::numeric_limits<double>::quiet_NaN(),
                                              const std::vector<std::string>& workingPointEntries = {});

    void setAuAuTightLogRegRuntimeConfig(const std::string& modelFile,
                                         double minIntercept = 0.80,
                                         double minSlope = 0.0,
                                         double maxScore = 1.0,
                                         double nonTightMinIntercept = 0.20,
                                         double nonTightMinSlope = 0.0,
                                         double nonTightMaxIntercept = 0.80,
                                         double nonTightMaxSlope = 0.0,
                                         double applyPtMin = std::numeric_limits<double>::quiet_NaN(),
                                         double applyPtMax = std::numeric_limits<double>::quiet_NaN(),
                                         const std::vector<std::string>& workingPointEntries = {});
    
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
    bool firstEventCuts(PHCompositeNode* topNode,
                        std::vector<std::string>& activeTrig,
                        bool applyVzCut = true,
                        bool* outMinimumBiasPass = nullptr,
                        bool* outTriggerPass = nullptr);
    void createHistos_Data();
    
    struct IsoAuditScalarStats
    {
        unsigned long long n = 0;
        double sum = 0.0;
        double sumsq = 0.0;
        std::vector<double> values;
        
        void fill(double x)
        {
            if (!std::isfinite(x)) return;
            ++n;
            sum += x;
            sumsq += x * x;
            values.push_back(x);
        }
    };
    
    struct IsoAuditMeanStats
    {
        unsigned long long n = 0;
        double sum = 0.0;
        
        void fill(double x)
        {
            if (!std::isfinite(x)) return;
            ++n;
            sum += x;
        }
    };
    
    struct IsoAuditEventFlow
    {
        unsigned long long evt_seen = 0;
        unsigned long long mandatory_nodes_ok = 0;
        unsigned long long minimum_bias_pass = 0;
        unsigned long long trigger_pass = 0;
        unsigned long long valid_centrality_info = 0;
        unsigned long long valid_reco_vertex_found = 0;
        unsigned long long vz_pass = 0;
        unsigned long long events_reaching_photon_loop = 0;
    };
    
    struct IsoAuditCell
    {
        unsigned long long n_seen_container = 0;
        unsigned long long n_pass_pt_floor = 0;
        unsigned long long n_pass_eta = 0;
        unsigned long long n_in_pt_bin = 0;
        unsigned long long n_inclusive = 0;
        
        unsigned long long n_total_finite = 0;
        unsigned long long n_fail_safe_overflow = 0;
        unsigned long long n_nonfinite_components = 0;
        unsigned long long n_negative_total = 0;
        unsigned long long n_negative_emcal = 0;
        
        IsoAuditScalarStats emcal;
        IsoAuditScalarStats hcalin;
        IsoAuditScalarStats hcalout;
        IsoAuditScalarStats total;
        
        unsigned long long n_sumdiff = 0;
        double sum_abs_sumdiff = 0.0;
        double max_abs_sumdiff = 0.0;
        
        unsigned long long n_isoPass = 0;
        unsigned long long n_gap = 0;
        unsigned long long n_nonIso = 0;
        unsigned long long n_decision = 0;
        double sum_thrIso = 0.0;
        double sum_eiso_minus_thrIso = 0.0;
        
        IsoAuditMeanStats tight;
        IsoAuditMeanStats nonTight;
        IsoAuditMeanStats neither;
        
        std::vector<std::string> ex_negative_total;
        std::vector<std::string> ex_near_zero_total;
        std::vector<std::string> ex_large_positive_total;
        std::vector<std::string> ex_overflow_total;
        std::vector<std::string> ex_large_sum_mismatch;
    };
    
    struct IsoAuditSample
    {
        int centIdx = -1;
        int ptIdx = -1;
        double ptGamma = 0.0;
        double etaGamma = 0.0;
        double phiGamma = 0.0;
        double eisoTot = 1e9;
        double eisoEmcal = 1e9;
        double eisoHcalIn = 1e9;
        double eisoHcalOut = 1e9;
        double thrIso = 0.0;
        double thrNonIso = 0.0;
        bool supportedCone = false;
        bool componentsFinite = false;
    };
    
    struct IsoAuditSkipSnapshot
    {
        long long evt_seen = 0;
        long long evt_accepted = 0;
        long long pho_total = 0;
        long long pho_early_E = 0;
        long long pho_eta_fail = 0;
        long long pho_etbin_out = 0;
        long long pho_reached_pre_iso = 0;
        long long pre_pass = 0;
        long long pre_fail_weta = 0;
        long long pre_fail_et1_low = 0;
        long long pre_fail_et1_high = 0;
        long long pre_fail_e11e33_high = 0;
        long long pre_fail_e32e35_low = 0;
        long long pre_fail_e32e35_high = 0;
        long long tight_tight = 0;
        long long tight_nonTight = 0;
        long long tight_neither = 0;
        unsigned long long noPhotonContainerEvents = 0;
        unsigned long long nonFiniteKinematics = 0;
        unsigned long long outsideConfiguredCentrality = 0;
        unsigned long long enteredInclusive = 0;
    };
    
    void initIsolationAudit();
    void recordIsolationAuditInclusive(const IsoAuditSample& sample);
    void recordIsolationAuditFollowup(const std::vector<std::string>& activeTrig,
                                      const IsoAuditSample& sample,
                                      TightTag tightTag);
    bool isolationAuditTargetsMet() const;
    std::string isoAuditCentLabel(int centIdx) const;
    std::string isoAuditPtLabel(int ptIdx) const;
    void printIsolationAuditProgress(bool force = false) const;
    IsoAuditSkipSnapshot makeIsolationAuditSkipSnapshot() const;
    void printIsolationAuditSkipSummary(bool force = false);
    void printIsolationAuditSummary() const;
    
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

    void fillAuAuEmbeddedTruthIsolationDiagnostics(
        PHCompositeNode* topNode,
        const std::vector<std::string>& activeTrig);
    
    void processCandidates(PHCompositeNode* topNode, const std::vector<std::string>& activeTrig);
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
    bool initAuAuTightBDTModelsIfNeeded() const;
    bool initAuAuTightMLPModelIfNeeded() const;
    bool initAuAuTightBDTMLPStackModelIfNeeded() const;
    bool initAuAuTightLogRegModelIfNeeded() const;
    double auauTightBDTFeatureValue(const std::string& feature,
                                    const PhotonClusterv1* pho,
                                    const SSVars& v) const;
    double predictAuAuTightBDTScore(const PhotonClusterv1* pho, const SSVars& v) const;
    double predictAuAuTightMLPScore(const PhotonClusterv1* pho, const SSVars& v) const;
    double predictAuAuTightBDTMLPStackScore(const PhotonClusterv1* pho, const SSVars& v) const;
    double predictAuAuTightLogRegScore(const PhotonClusterv1* pho, const SSVars& v) const;
    int auauTightBDTCentBinIndex() const;
    int auauTightBDTPtBinIndex(double pt) const;
    struct AuAuTightBDTWorkingPoint
    {
        std::string variant;
        bool binned = false;
        bool grid2d = false;
        bool centralityLinear = false;
        double intercept = std::numeric_limits<double>::quiet_NaN();
        double slope = 0.0;
        double ptMin = std::numeric_limits<double>::quiet_NaN();
        double ptMax = std::numeric_limits<double>::quiet_NaN();
        double maxScore = 1.0;
        std::vector<double> edges;
        std::vector<double> centEdges;
        std::vector<double> thresholds;
    };
    void parseAuAuTightBDTWorkingPointEntries(const std::vector<std::string>& entries);
    const AuAuTightBDTWorkingPoint* activeAuAuTightBDTWorkingPoint() const;
    double configuredAuAuTightBDTMin(double et) const;
    double configuredAuAuTightBDTMax(double et) const;
    double configuredAuAuNonTightBDTMin(const SSVars& v) const;
    double configuredAuAuNonTightBDTMax(const SSVars& v) const;
    bool configuredAuAuNonTightBDTPass(double score, const SSVars& v) const;
    void parseAuAuTightMLPWorkingPointEntries(const std::vector<std::string>& entries);
    const AuAuTightBDTWorkingPoint* activeAuAuTightMLPWorkingPoint() const;
    double configuredAuAuTightMLPMin(double et) const;
    double configuredAuAuTightMLPMax(double et) const;
    void parseAuAuTightBDTMLPStackWorkingPointEntries(const std::vector<std::string>& entries);
    const AuAuTightBDTWorkingPoint* activeAuAuTightBDTMLPStackWorkingPoint() const;
    double configuredAuAuTightBDTMLPStackMin(double et) const;
    double configuredAuAuTightBDTMLPStackMax(double et) const;
    void parseAuAuTightLogRegWorkingPointEntries(const std::vector<std::string>& entries);
    const AuAuTightBDTWorkingPoint* activeAuAuTightLogRegWorkingPoint() const;
    double configuredAuAuTightLogRegMin(double et) const;
    double configuredAuAuTightLogRegMax(double et) const;
    const PhotonClusterv1* findMatchedPhotonByKinematics(const RawClusterContainer* container,
                                                         const PhotonClusterv1* ref) const;
    bool   passesPhotonPreselection(const SSVars& v);
    TightTag classifyPhotonTightness(const SSVars& v);
    double configuredTightBDTMin(double et) const;
    double configuredNonTightBDTMin(double et) const;
    double configuredNonTightBDTMax(double et) const;
    bool configuredTightBDTPass(double score, double et) const;
    bool configuredNonTightBDTPass(double score, double et) const;
    TightTag configuredVariantABDTTag(double score, double et) const;
    struct TruthSignalPhotonInfo;
    void initAuAuBDTTrainingTree();
    void fillAuAuBDTTrainingTree(const SSVars& v,
                                 double eta,
                                 double phi,
                                 int clusterIndex,
                                 double eiso,
                                 int ptIdx,
                                 int centIdx,
                                 bool isSignal,
                                 int npbLabel = -1,
                                 int isNPB = -1,
                                 double clusterMbdDeltaT = std::numeric_limits<double>::quiet_NaN(),
                                 double mbdTime = std::numeric_limits<double>::quiet_NaN(),
                                 bool hasAwayJet = false,
                                 double eisoR30 = std::numeric_limits<double>::quiet_NaN(),
                                 double eisoR40 = std::numeric_limits<double>::quiet_NaN(),
                                 int truthMatchFound = -1,
                                 int truthPhotonClass = -1,
                                 double truthIsoEt = std::numeric_limits<double>::quiet_NaN(),
                                 int truthIsoPass = -1,
                                 int clusterTruthTrackId = -1,
                                 int clusterTruthPid = 0,
                                 int clusterTruthBarcode = -1,
                                 float truthEnergyContribution = std::numeric_limits<float>::quiet_NaN(),
                                 int sourceRole = 0,
                                 int sourceSampleCode = 0,
                                 int ppg12SourceRoleLabel = -1,
                                 double truthIsoEtR03 = std::numeric_limits<double>::quiet_NaN(),
                                 double truthIsoEtR04 = std::numeric_limits<double>::quiet_NaN(),
                                 int truthIsoValid = -1,
                                 int truthHepMCAssociationValid = -1);
    void initAuAuPhotonCandidateSkimTree();
    void fillAuAuPhotonCandidateSkimTree(PHCompositeNode* topNode,
                                         const PhotonClusterv1* pho,
                                         const std::vector<std::string>& activeTrig,
                                         const SSVars& v,
                                         double eta,
                                         double phi,
                                         double eiso,
                                         double eisoR30,
                                         double eisoR40,
                                         int ptIdx,
                                         int centIdx,
                                         bool preselectionPass,
                                         TightTag activeTightTag,
                                         bool iso,
                                         bool nonIso,
                                         double isoThreshold,
                                         double nonIsoThreshold,
                                         double isoSideGap,
                                         bool haveTruthLabel,
                                         bool isSignal,
                                         int clusterTruthTrackId = -1,
                                         int clusterTruthPid = 0,
                                         int clusterTruthBarcode = -1,
                                         float clusterTruthEContrib = std::numeric_limits<float>::quiet_NaN());
    void initTHE44PythiaAutopsyTree();
    void fillTHE44PythiaAutopsyTree(const SSVars& v,
                                    double eta,
                                    double phi,
                                    double eiso,
                                    int ptIdx,
                                    int centIdx,
                                    bool isSignal,
                                    int clusterTruthTrackId,
                                    int clusterTruthPid,
                                    int clusterTruthBarcode,
                                    float eContrib,
                                    const TruthSignalPhotonInfo& matchedTruth,
                                    const HepMC::GenEvent* evt);
    bool isPPG12DataNPBTaggedCluster(const SSVars& v,
                                     double phi,
                                     double& clusterMbdDeltaT,
                                     double& mbdTime,
                                     bool& hasAwayJet,
                                     bool requireTaggingEnabled = true) const;
    void initJetMLTrainingTree();
    void fillJetMLTrainingTree(const std::string& rKey,
                               double Rjet,
                               double leadPtGamma,
                               double leadEtaGamma,
                               double leadPhiGamma,
                               const Jet* recoJet,
                               bool recoIsRecoil,
                               const Jet* matchedTruthJet,
                               double matchDR);
    bool initJetMLModelIfNeeded();
    double predictJetMLDeltaPt(const std::string& rKey,
                               double Rjet,
                               double leadPtGamma,
                               double leadEtaGamma,
                               double leadPhiGamma,
                               const Jet* recoJet) const;
    double jetMLFeatureValue(const std::string& feature,
                             const std::string& rKey,
                             double Rjet,
                             double leadPtGamma,
                             double leadEtaGamma,
                             double leadPhiGamma,
                             const Jet* recoJet) const;
    
    
    // Isolation helpers
    double eiso(const RawCluster* clus, PHCompositeNode* topNode) const;
    double eisoForCone(const RawCluster* clus, double coneR) const;
    bool   isIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const;
    bool   isNonIsolated(const RawCluster* clus, double et_gamma, PHCompositeNode* topNode) const;
    
    // Unified truth-photon classification and isolation metadata (SIM only).
    // Definition
    //   |eta| < 0.7, PID=22, embedded G4 primary photon matched by HepMC barcode,
    //   prompt classification via CaloAna photon_type logic:
    //     - walk back photon-in/photon-out vertices
    //     - direct=1 if 2->2 with |pdg|<=22 on all legs
    //     - frag  =2 if 1->2 with |incoming pdg|<=11 and outgoing contains incoming pid (and photon)
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
        int embedId = 0;
        int photonClass = 0;
        bool truthIsoPass = false;
        const HepMC::GenParticle* hep = nullptr;
        const PHG4Particle* g4 = nullptr;
    };
    using TruthSignalPhotonMap = std::map<int, TruthSignalPhotonInfo>;

    // PPG12 photon truth association:
    //   cluster_truthtrkID = clustereval.max_truth_primary_particle_by_energy(cluster)->get_track_id()
    //   signal if that track id maps to a truth photon passing isTruthPromptIsolatedSignalPhoton().
    // This deliberately does not impose an additional truth-reco ΔR cut.
    TruthSignalPhotonMap buildPPG12TruthPhotonMap(
        PHCompositeNode* topNode,
        RJReplayFoundationV1::TruthPhotonInventoryV1* inventory=nullptr) const;
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
    //   - signal: reco cluster best-matched to a truth photon that satisfies isTruthPromptIsolatedSignalPhoton()
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
    TH1F* getOrBookRecoClusterEtFineDiagHist(const std::string& trig,
                                             const std::string& base,
                                             int centIdx);
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
    TH1I* getOrBookSigABCDLeakageHist(const std::string& trig, int ptIdx, int centIdx,
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
    void bookPPG12TableQASchema(const std::vector<std::string>& activeTrig);
    void fillPPG12TableQA(const std::vector<std::string>& activeTrig,
                          const SSVars& v,
                          double eisoEt,
                          int centIdx,
                          int cutIdx,
                          double rowWeight = 1.0);
    double ppg12TableQABDTScore(const SSVars& v) const;

    TH3F* getOrBookAuAuDualViewScoreIsoSurface(const std::string& trig,
                                               int centIdx,
                                               const std::string& category);
    void fillAuAuDualViewScoreIsoSurface(const std::vector<std::string>& activeTrig,
                                         const SSVars& v,
                                         double eisoEt,
                                         int centIdx,
                                         const std::string& category);

    TH2F* getOrBookAuAuFig25CorrelationSurface(const std::string& trig,
                                                int centIdx,
                                                const std::string& axisKey);
    void fillAuAuFig25CorrelationSurfaces(const std::vector<std::string>& activeTrig,
                                          const SSVars& v,
                                          double eisoEt,
                                          int centIdx,
                                          bool isBackground);

    // Opt-in embedded-simulation diagnostics for the low-calorimeter-energy
    // centrality band. PMT objects are filled before the optional embedded
    // MinimumBiasInfo gate; downstream physics output may then reject failures.
    void bookMbdPmtLowCaloDiagnostics();
    void fillMbdPmtLowCaloDiagnostics(const std::vector<std::string>& activeTrig,
                                      double emcalEnergy,
                                      double ihcalEnergy,
                                      double ohcalEnergy,
                                      double totalCaloEnergy);
    bool buildSepdChannelMapForDiagnostics();
    
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
    TH2F* getOrBookUnfoldRecoCombinatoricPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldRecoFakesPtXJIncl (const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldTruthMissesPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);

    // Event-leading recoil-jet unfolding family, separate from the inclusive
    // per-photon recoil family. `component` is reco, truth, response, fake,
    // miss, or selectionLoss.
    TH2F* getOrBookUnfoldLeadPtXJComponent(const std::string& trig,
                                           const std::string& rKey,
                                           int centIdx,
                                           const std::string& component);

    TH2F* getOrBookUnfoldRecoJetMLPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldRecoJetMLPtXJInclSidebandC(const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldResponseJetMLPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldRecoCombinatoricJetMLPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldRecoFakesJetMLPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);
    TH2F* getOrBookUnfoldTruthMissesJetMLPtXJIncl(const std::string& trig, const std::string& rKey, int centIdx);
    
    // photon-only unfolding (for N_gamma normalization): DATA+SIM reco, SIM truth
    TH1F* getOrBookUnfoldRecoPhoPtGamma        (const std::string& trig, int centIdx);
    TH1F* getOrBookUnfoldTruthPhoPtGamma       (const std::string& trig, int centIdx);
    TH2F* getOrBookUnfoldResponsePhoPtGamma    (const std::string& trig, int centIdx);
    TH1F* getOrBookUnfoldRecoPhoFakesPtGamma   (const std::string& trig, int centIdx);
    TH1F* getOrBookUnfoldTruthPhoMissesPtGamma (const std::string& trig, int centIdx);

    // photon-efficiency stage QA (SIM only): truth-pT binned denominators and numerators
    TH1F* getOrBookPhotonEffStagePtGamma       (const std::string& trig, int centIdx,
                                                const std::string& base,
                                                const std::string& yTitle,
                                                bool useIsoViewSuffix);
    TH1F* getOrBookPhotonEffTruthDenPtGamma    (const std::string& trig, int centIdx);
    TH1F* getOrBookPhotonEffRecoPtGamma        (const std::string& trig, int centIdx);
    TH1F* getOrBookPhotonEffRecoIsoPtGamma     (const std::string& trig, int centIdx);
    TH1F* getOrBookPhotonEffRecoTightPtGamma   (const std::string& trig, int centIdx);
    TH1F* getOrBookPhotonEffRecoTightIsoPtGamma(const std::string& trig, int centIdx);
    
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
    void recordReplayTruthPhotonMissOccurrence(const std::string& trigger,
                                               int centralityBinIndex,
                                               double truthPt,
                                               const char* decisionPath);
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
    bool m_isSimEmbedded = false;
    bool m_doPi0Analysis = false;
    bool m_requireTowerInfoTruthMatching = false;
    bool m_towerInfoTruthMatchingContractViolation = false;
    std::uint64_t m_towerInfoTruthMatchingChecks = 0;
    
    // Vertex / event selection
    GlobalVertex* m_vtx = nullptr;
    float m_vx = 0.0f;
    float m_vy = 0.0f;
    float m_vz = 0.0f;
    
    bool  m_useVzCut = true;
    float m_vzCut    = 30.0f;
    bool  m_useMinBiasClassifier = false;
    
    // Centrality
    int m_centBin = -1;                 // 0..99 (Au+Au), or -1 in pp
    double m_centPercent = -1.0;        // event centrality percentile (float, used for reweighting)
    std::vector<int> m_centEdges;       // centrality bin edges, e.g. {0,10,20,...,100}
    std::uint64_t m_validCentralityObservedEvents = 0;
    std::uint64_t m_invalidCentralityObservedEvents = 0;
    std::uint64_t m_staleCentralityGuardedEvents = 0;

    // Event-level calorimeter sums cached once per accepted event.
    float m_eventCaloCemcEnergy = std::numeric_limits<float>::quiet_NaN();
    float m_eventCaloIhcalEnergy = std::numeric_limits<float>::quiet_NaN();
    float m_eventCaloOhcalEnergy = std::numeric_limits<float>::quiet_NaN();
    float m_eventCaloTotalEnergy = std::numeric_limits<float>::quiet_NaN();
    float m_eventCaloLog10TotalEnergyPlus1 = std::numeric_limits<float>::quiet_NaN();
    double m_eventTruthVertexZ = std::numeric_limits<double>::quiet_NaN();
    double m_eventMbdTotalCharge = std::numeric_limits<double>::quiet_NaN();
    int m_eventCaloRequireIsGood = -1;
    
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
    double m_auauNPBCut = 0.5;
    double m_auauTightBDTMinIntercept = 0.8333333333333334;
    double m_auauTightBDTMinSlope = -0.003333333333333336;
    double m_auauTightBDTMax = 1.0;
    double m_auauNonTightBDTMinIntercept = 0.7333333333333333;
    double m_auauNonTightBDTMinSlope = -0.01333333333333333;
    double m_auauNonTightBDTMaxIntercept = 0.6666666666666666;
    double m_auauNonTightBDTMaxSlope = 0.003333333333333336;
    std::string m_auauNonTightBDTSidebandMode = "etLinear";
    double m_auauNonTightBDTRelativeMinOffset = -0.20;
    double m_auauNonTightBDTRelativeMaxOffset = -0.03;
    std::vector<int> m_auauTightBDTCentDepEdges;
    std::vector<std::string> m_auauTightBDTCentDepScoreNames;
    std::string m_auauTightBDTModelFile;
    std::string m_auauTightBDTScoringMode;
    std::vector<std::string> m_auauTightBDTFeatures;
    std::vector<std::string> m_auauTightBDTCentDepModelFiles;
    std::vector<double> m_auauTightBDTPtBinEdges;
    std::vector<std::string> m_auauTightBDTPtBinModelFiles;
    std::vector<std::string> m_auauTightBDTPtCentModelFiles;
    std::string m_auauTightBDTPtFallbackModelFile;
    std::vector<std::string> m_auauTightBDTPtFallbackCentModelFiles;
    double m_auauTightBDTPtFallbackMin = 35.0;
    double m_auauTightBDTPtFallbackMax = 40.0;
    double m_auauTightBDTApplyPtMin = std::numeric_limits<double>::quiet_NaN();
    double m_auauTightBDTApplyPtMax = std::numeric_limits<double>::quiet_NaN();
    std::vector<AuAuTightBDTWorkingPoint> m_auauTightBDTWorkingPoints;
    bool m_explicitAuAuTightBDTConfig = false;
    mutable bool m_auauTightBDTModelInitAttempted = false;
    mutable std::unique_ptr<TMVA::Experimental::RBDT> m_auauTightBDTModel;
    mutable std::vector<std::unique_ptr<TMVA::Experimental::RBDT>> m_auauTightBDTCentDepModels;
    mutable std::vector<std::unique_ptr<TMVA::Experimental::RBDT>> m_auauTightBDTPtBinModels;
    mutable std::vector<std::unique_ptr<TMVA::Experimental::RBDT>> m_auauTightBDTPtCentModels;
    mutable std::unique_ptr<TMVA::Experimental::RBDT> m_auauTightBDTPtFallbackModel;
    mutable std::vector<std::unique_ptr<TMVA::Experimental::RBDT>> m_auauTightBDTPtFallbackCentModels;
    struct AuAuTightMLPLayer
    {
        std::vector<std::vector<double>> weightsOutByIn;
        std::vector<double> bias;
    };
    struct AuAuTightMLPModel
    {
        std::string schema;
        std::string variant;
        std::vector<std::string> features;
        std::vector<double> mean;
        std::vector<double> scale;
        double clip = 8.0;
        double outputTemperature = 1.0;
        std::string activation = "tanh";
        std::vector<AuAuTightMLPLayer> layers;
        bool valid = false;
    };
    std::string m_auauTightMLPModelFile;
    std::string m_auauTightMLPScoringMode;
    double m_auauTightMLPMinIntercept = 0.80;
    double m_auauTightMLPMinSlope = 0.0;
    double m_auauTightMLPMax = 1.0;
    double m_auauNonTightMLPMinIntercept = 0.20;
    double m_auauNonTightMLPMinSlope = 0.0;
    double m_auauNonTightMLPMaxIntercept = 0.80;
    double m_auauNonTightMLPMaxSlope = 0.0;
    double m_auauTightMLPApplyPtMin = std::numeric_limits<double>::quiet_NaN();
    double m_auauTightMLPApplyPtMax = std::numeric_limits<double>::quiet_NaN();
    std::vector<AuAuTightBDTWorkingPoint> m_auauTightMLPWorkingPoints;
    bool m_explicitAuAuTightMLPConfig = false;
    mutable bool m_auauTightMLPModelInitAttempted = false;
    mutable AuAuTightMLPModel m_auauTightMLPModel;

    struct AuAuStackTree
    {
        std::vector<int> childrenLeft;
        std::vector<int> childrenRight;
        std::vector<int> feature;
        std::vector<double> threshold;
        std::vector<double> value;
    };
    struct AuAuStackSubModel
    {
        std::string kind;
        std::vector<std::string> featureNames;
        std::vector<double> impute;
        std::vector<double> mean;
        std::vector<double> scale;
        std::vector<double> coef;
        double intercept = 0.0;
        double initialRawScore = 0.0;
        double learningRate = 1.0;
        std::vector<AuAuStackTree> trees;
        std::string activation = "relu";
        std::string output = "sigmoid";
        std::vector<AuAuTightMLPLayer> layers;
        bool valid = false;
    };
    struct AuAuStackRoute
    {
        std::string label;
        double ptLo = std::numeric_limits<double>::quiet_NaN();
        double ptHi = std::numeric_limits<double>::quiet_NaN();
        double centLo = std::numeric_limits<double>::quiet_NaN();
        double centHi = std::numeric_limits<double>::quiet_NaN();
        bool hasCent = false;
        AuAuStackSubModel model;
    };
    struct AuAuTightBDTMLPStackModel
    {
        std::string schema;
        std::string name;
        AuAuStackSubModel globalModel;
        std::vector<AuAuStackRoute> routes;
        bool valid = false;
    };
    std::string m_auauTightBDTMLPStackModelFile;
    double m_auauTightBDTMLPStackMinIntercept = 0.80;
    double m_auauTightBDTMLPStackMinSlope = 0.0;
    double m_auauTightBDTMLPStackMax = 1.0;
    double m_auauNonTightBDTMLPStackMinIntercept = 0.20;
    double m_auauNonTightBDTMLPStackMinSlope = 0.0;
    double m_auauNonTightBDTMLPStackMaxIntercept = 0.80;
    double m_auauNonTightBDTMLPStackMaxSlope = 0.0;
    double m_auauTightBDTMLPStackApplyPtMin = std::numeric_limits<double>::quiet_NaN();
    double m_auauTightBDTMLPStackApplyPtMax = std::numeric_limits<double>::quiet_NaN();
    std::vector<AuAuTightBDTWorkingPoint> m_auauTightBDTMLPStackWorkingPoints;
    bool m_explicitAuAuTightBDTMLPStackConfig = false;
    mutable bool m_auauTightBDTMLPStackModelInitAttempted = false;
    mutable AuAuTightBDTMLPStackModel m_auauTightBDTMLPStackModel;

    std::string m_auauTightLogRegModelFile;
    double m_auauTightLogRegMinIntercept = 0.80;
    double m_auauTightLogRegMinSlope = 0.0;
    double m_auauTightLogRegMax = 1.0;
    double m_auauNonTightLogRegMinIntercept = 0.20;
    double m_auauNonTightLogRegMinSlope = 0.0;
    double m_auauNonTightLogRegMaxIntercept = 0.80;
    double m_auauNonTightLogRegMaxSlope = 0.0;
    double m_auauTightLogRegApplyPtMin = std::numeric_limits<double>::quiet_NaN();
    double m_auauTightLogRegApplyPtMax = std::numeric_limits<double>::quiet_NaN();
    std::vector<AuAuTightBDTWorkingPoint> m_auauTightLogRegWorkingPoints;
    bool m_explicitAuAuTightLogRegConfig = false;
    mutable bool m_auauTightLogRegModelInitAttempted = false;
    mutable AuAuTightBDTMLPStackModel m_auauTightLogRegModel;

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
    
    
    // Isolation WP (PPG12 nominal)
    double m_isoA      = 1.08128;
    double m_isoB      = 0.0299107;
    double m_isoGap    = 1.0;
    double m_isoFixed  = 2.0;          // used ONLY when m_isSlidingIso==false (RECO)
    double m_truthIsoMaxGeV = 4.0;     // TRUTH isolation max (independent of sliding/fixed mode)
    double m_isoConeR  = 0.3;
    double m_isoTowMin = 0.0;
    bool   m_isSlidingIso = true;
    bool   m_useTopoClusterIsolationForEiso = false;
    std::vector<IsoView> m_internalIsoViews;
    std::string m_activeIsoViewSuffix;
    
    // Per-centrality sliding WPs (indexed by findCentBin result; falls back to global if empty)
    std::vector<CentIsoWP> m_centIsoWPs;
    std::vector<CentIsoWP> m_centIsoWPsR30;
    std::vector<CentIsoWP> m_centIsoWPsR40;
    void getIsoParamsForCone(int centIdx, double coneR, double& outA, double& outB, double& outGap) const;
    void getIsoParams(int centIdx, double& outA, double& outB, double& outGap) const;

    // Blair-style SIM event reweighting
    bool        m_vertexReweightOn = false;
    std::string m_vertexReweightFile;
    std::string m_vertexReweightHist = "data_over_MC_ratios/h_zvtx_ratio_data_over_photonJet";
    TH1*        m_vertexReweightH = nullptr;

    bool        m_centralityReweightOn = false;
    std::string m_centralityReweightFile;
    std::string m_centralityReweightHist = "nom_cent_rw_hist";
    TH1*        m_centralityReweightH = nullptr;

    double      m_mcVertexWeight = 1.0;
    double      m_mcCentralityWeight = 1.0;
    double      m_mcEventWeight = 1.0;
    
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
    
    // Calo tower bundles (node cache)
    struct CaloBundle
    {
        TowerInfoContainer*   towers = nullptr;
        RawTowerGeomContainer* geom  = nullptr;
        double                sumEt  = 0.0;
    };
    
    std::vector<std::tuple<std::string, std::string, std::string>> m_caloInfo;
    std::map<std::string, CaloBundle> m_calo;
    MbdPmtContainer* m_mbdpmts = nullptr;
    MbdOut* m_mbdout = nullptr;
    MbdGeom* m_mbdgeom = nullptr;
    TowerInfoContainer* m_sepdTowers = nullptr;
    EpdGeom* m_epdgeom = nullptr;
    std::vector<unsigned int> m_sepdDiagnosticKeys;
    bool m_sepdDiagnosticMapAttempted = false;
    bool m_sepdDiagnosticMapReady = false;
    
    // -------------------------------------------------------------------------
    // parallel jet containers by radius key
    // -------------------------------------------------------------------------
    std::map<std::string, JetContainer*> m_jets;
    std::map<std::string, JetContainer*> m_jetsRaw;
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

    bool m_auauBDTTrainingTreeEnabled = false;
    bool m_auauBDTExtractOnly = false;
    long long m_auauBDTTrainingTreeMaxEntries = 0;
    long long m_auauBDTTrainingTreeEntries = 0;
    TTree* m_auauBDTTrainingTree = nullptr;
    bool m_auauBDTNPBDataTaggingEnabled = false;
    double m_auauNPBTagDeltaTCut = -7.0;
    double m_auauNPBTagWetaMin = 0.0;
    double m_auauNPBTagAwayJetPtMin = 5.0;
    double m_auauNPBTagAwayJetDPhiMin = 1.5707963267948966;
    double m_auauNPBTagTimeSampleNs = 17.6;
    double m_auauNPBMbdT0Offset = 0.0;
    bool m_ppg12TableQAEnabled = false;
    bool m_auauDualViewDiagnosticsEnabled = false;
    bool m_auauFig25CorrelationDiagnosticsEnabled = false;
    bool m_mbdPmtLowCaloDiagnosticsEnabled = false;
    bool m_requireEmbeddedMinBiasClassifier = false;
    // -1: not evaluated, 0: node missing, 1: classifier fail, 2: classifier pass.
    int m_embeddedMinBiasDecision = -1;
    bool m_mbdPmtGeometryFilled = false;
    std::set<std::string> m_ppg12TableQASchemaBookedTriggers;

    int m_bdtTrain_run = 0;
    long long m_bdtTrain_evt = 0;
    int m_bdtTrain_is_signal = 0;
    int m_bdtTrain_truth_match_found = -1;
    int m_bdtTrain_truth_photon_class = -1;
    int m_bdtTrain_truth_is_prompt = -1;
    float m_bdtTrain_truth_iso_et = -999.0f;
    int m_bdtTrain_truth_iso_pass = -1;
    int m_bdtTrain_cluster_truth_track_id = -1;
    int m_bdtTrain_cluster_truth_pid = 0;
    int m_bdtTrain_cluster_truth_barcode = -1;
    int m_bdtTrain_source_role = 0;
    int m_bdtTrain_source_sample_code = 0;
    int m_bdtTrain_ppg12_source_role_label = -1;
    int m_bdtTrain_pt_bin = -1;
    int m_bdtTrain_cent_bin = -1;
    int m_bdtTrain_minbias_decision = -1;
    float m_bdtTrain_pt = 0.0f;
    float m_bdtTrain_eta = 0.0f;
    float m_bdtTrain_phi = 0.0f;
    float m_bdtTrain_cent = -1.0f;
    float m_bdtTrain_vz = 0.0f;
    float m_bdtTrain_weight = 1.0f;
    float m_bdtTrain_event_calo_cemc_energy = 0.0f;
    float m_bdtTrain_event_calo_ihcal_energy = 0.0f;
    float m_bdtTrain_event_calo_ohcal_energy = 0.0f;
    float m_bdtTrain_event_calo_total_energy = 0.0f;
    float m_bdtTrain_event_calo_log10_total_energy_plus1 = 0.0f;
    float m_bdtTrain_eiso = 0.0f;
    float m_bdtTrain_eiso_r30 = 0.0f;
    float m_bdtTrain_eiso_r40 = 0.0f;
    int m_bdtTrain_npb_label = -1;
    int m_bdtTrain_is_npb = -1;
    float m_bdtTrain_cluster_mean_time = -999.0f;
    float m_bdtTrain_mbd_time = -999.0f;
    float m_bdtTrain_cluster_mbd_delta_t = -999.0f;
    int m_bdtTrain_npb_has_away_jet = 0;
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
    float m_bdtTrain_npb_score = -2.0f;
    float m_bdtTrain_auau_npb_score = -2.0f;
    float m_bdtTrain_auau_tight_bdt_score = -2.0f;
    float m_bdtTrain_auau_tight_mlp_score = -2.0f;
    float m_bdtTrain_auau_tight_bdt_mlp_score = -2.0f;
    float m_bdtTrain_auau_tight_logreg_score = -2.0f;

    bool m_auauPhotonCandidateSkimEnabled = false;
    long long m_auauPhotonCandidateSkimMaxEntries = 0;
    long long m_auauPhotonCandidateSkimEntries = 0;
    TTree* m_auauPhotonCandidateSkimTree = nullptr;
    int m_phoSkim_run = 0;
    long long m_phoSkim_evt = 0;
    long long m_phoSkim_event_count = 0;
    int m_phoSkim_is_sim = 0;
    int m_phoSkim_is_sim_embedded = 0;
    int m_phoSkim_sample_code = 0;
    int m_phoSkim_pt_bin = -1;
    int m_phoSkim_cent_bin = -1;
    std::string m_phoSkim_active_triggers;
    unsigned long long m_phoSkim_gl1_scaled_vector = 0ULL;
    float m_phoSkim_cluster_et = 0.0f;
    float m_phoSkim_cluster_eta = 0.0f;
    float m_phoSkim_cluster_phi = 0.0f;
    float m_phoSkim_centrality = -1.0f;
    float m_phoSkim_vertexz = 0.0f;
    float m_phoSkim_event_weight = 1.0f;
    float m_phoSkim_event_calo_cemc_energy = 0.0f;
    float m_phoSkim_event_calo_ihcal_energy = 0.0f;
    float m_phoSkim_event_calo_ohcal_energy = 0.0f;
    float m_phoSkim_event_calo_total_energy = 0.0f;
    float m_phoSkim_weta = 0.0f;
    float m_phoSkim_wphi = 0.0f;
    float m_phoSkim_weta33 = 0.0f;
    float m_phoSkim_wphi33 = 0.0f;
    float m_phoSkim_weta35 = 0.0f;
    float m_phoSkim_wphi53 = 0.0f;
    float m_phoSkim_et1 = 0.0f;
    float m_phoSkim_et2 = 0.0f;
    float m_phoSkim_et3 = 0.0f;
    float m_phoSkim_et4 = 0.0f;
    float m_phoSkim_e11e33 = 0.0f;
    float m_phoSkim_e32e35 = 0.0f;
    float m_phoSkim_e11e22 = 0.0f;
    float m_phoSkim_e11e13 = 0.0f;
    float m_phoSkim_e11e15 = 0.0f;
    float m_phoSkim_e11e17 = 0.0f;
    float m_phoSkim_e11e31 = 0.0f;
    float m_phoSkim_e11e51 = 0.0f;
    float m_phoSkim_e11e71 = 0.0f;
    float m_phoSkim_e22e33 = 0.0f;
    float m_phoSkim_e22e35 = 0.0f;
    float m_phoSkim_e22e37 = 0.0f;
    float m_phoSkim_e22e53 = 0.0f;
    float m_phoSkim_w32 = 0.0f;
    float m_phoSkim_w52 = 0.0f;
    float m_phoSkim_w72 = 0.0f;
    float m_phoSkim_mean_time = -999.0f;
    float m_phoSkim_eiso = 0.0f;
    float m_phoSkim_eiso_r30 = 0.0f;
    float m_phoSkim_eiso_r40 = 0.0f;
    float m_phoSkim_iso_threshold = 0.0f;
    float m_phoSkim_noniso_threshold = 0.0f;
    float m_phoSkim_iso_side_gap = 0.0f;
    int m_phoSkim_iso_pass = 0;
    int m_phoSkim_noniso_pass = 0;
    int m_phoSkim_iso_gap = 0;
    int m_phoSkim_preselection_pass = 0;
    int m_phoSkim_active_tight_tag = 0;
    int m_phoSkim_active_abcd_region = 0;
    float m_phoSkim_npb_score = -2.0f;
    int m_phoSkim_npb_pass = 0;
    float m_phoSkim_auau_npb_score = -2.0f;
    float m_phoSkim_ppg12_tight_bdt_score = -2.0f;
    float m_phoSkim_auau_tight_bdt_score = -2.0f;
    float m_phoSkim_baseline_wp80_threshold = -2.0f;
    int m_phoSkim_baseline_bdt_tight = 0;
    int m_phoSkim_baseline_bdt_nontight = 0;
    int m_phoSkim_baseline_bdt_abcd_region = 0;
    int m_phoSkim_box_tight = 0;
    int m_phoSkim_box_nontight = 0;
    int m_phoSkim_box_fail_count = 0;
    int m_phoSkim_box_abcd_region = 0;
    int m_phoSkim_have_truth_label = 0;
    int m_phoSkim_is_truth_signal = 0;
    int m_phoSkim_cluster_truth_track_id = -1;
    int m_phoSkim_cluster_truth_pid = 0;
    int m_phoSkim_cluster_truth_barcode = -1;
    float m_phoSkim_cluster_truth_econtrib = -999.0f;
    std::string m_phoSkim_lead_recoil_rkey;
    float m_phoSkim_lead_recoil_pt = -1.0f;
    float m_phoSkim_lead_recoil_eta = -999.0f;
    float m_phoSkim_lead_recoil_phi = -999.0f;
    float m_phoSkim_lead_recoil_dphi = -999.0f;
    float m_phoSkim_lead_xj = -1.0f;

    bool m_the44PythiaAutopsyEnabled = false;
    long long m_the44PythiaAutopsyMaxEntries = 200;
    long long m_the44PythiaAutopsyEntries = 0;
    long long m_the44PythiaAutopsyHighBkgEntries = 0;
    long long m_the44PythiaAutopsyMidSigEntries = 0;
    double m_the44PythiaAutopsyHighBDTMin = 0.80;
    double m_the44PythiaAutopsyMidBDTMin = 0.45;
    double m_the44PythiaAutopsyMidBDTMax = 0.65;
    double m_the44PythiaAutopsyMinPt = 15.0;
    double m_the44PythiaAutopsyMaxPt = 35.0;
    double m_the44PythiaAutopsyParticleCone = 0.30;
    double m_the44PythiaAutopsyParticleMinPt = 0.30;
    int m_the44PythiaAutopsyMaxParticles = 64;
    TTree* m_the44PythiaAutopsyTree = nullptr;
    int m_the44_run = 0;
    long long m_the44_evt = 0;
    int m_the44_category = 0;
    int m_the44_is_signal = 0;
    int m_the44_pt_bin = -1;
    int m_the44_cent_bin = -1;
    int m_the44_source_sample_code = 0;
    float m_the44_cluster_et = 0.0f;
    float m_the44_cluster_eta = 0.0f;
    float m_the44_cluster_phi = 0.0f;
    float m_the44_centrality = -1.0f;
    float m_the44_reco_eiso = 0.0f;
    float m_the44_bdt_score = -2.0f;
    float m_the44_npb_score = -2.0f;
    float m_the44_weta = 0.0f;
    float m_the44_wphi = 0.0f;
    float m_the44_e11e33 = 0.0f;
    float m_the44_e32e35 = 0.0f;
    int m_the44_cluster_truth_track_id = -1;
    int m_the44_cluster_truth_pid = 0;
    int m_the44_cluster_truth_barcode = -1;
    float m_the44_cluster_truth_econtrib = -999.0f;
    int m_the44_matched_signal_barcode = -1;
    float m_the44_matched_signal_pt = -1.0f;
    float m_the44_matched_signal_eta = -999.0f;
    float m_the44_matched_signal_phi = -999.0f;
    float m_the44_matched_signal_iso = -999.0f;
    std::vector<int> m_the44_part_pdg;
    std::vector<int> m_the44_part_status;
    std::vector<int> m_the44_part_barcode;
    std::vector<int> m_the44_part_parent0_pdg;
    std::vector<int> m_the44_part_parent1_pdg;
    std::vector<float> m_the44_part_pt;
    std::vector<float> m_the44_part_eta;
    std::vector<float> m_the44_part_phi;
    std::vector<float> m_the44_part_dr;

    bool m_jetMLTrainingTreeEnabled = false;
    long long m_jetMLTrainingTreeMaxEntries = 0;
    long long m_jetMLTrainingTreeEntries = 0;
    TTree* m_jetMLTrainingTree = nullptr;

    bool m_jetMLCorrectionEnabled = false;
    bool m_jetMLModelInitAttempted = false;
    std::string m_jetMLModelFile;
    std::vector<std::string> m_jetMLFeatures;
    std::unique_ptr<TMVA::Experimental::RBDT> m_jetMLModel;

    int m_jetMLTrain_run = 0;
    long long m_jetMLTrain_evt = 0;
    std::string m_jetMLTrain_sample;
    std::string m_jetMLTrain_rKey;
    int m_jetMLTrain_is_sim = 0;
    int m_jetMLTrain_is_recoil = 0;
    int m_jetMLTrain_is_matched = 0;
    float m_jetMLTrain_R = 0.0f;
    float m_jetMLTrain_centrality = -1.0f;
    float m_jetMLTrain_vertexz = 0.0f;
    float m_jetMLTrain_event_weight = 1.0f;
    float m_jetMLTrain_photon_pt = 0.0f;
    float m_jetMLTrain_photon_eta = 0.0f;
    float m_jetMLTrain_photon_phi = 0.0f;
    float m_jetMLTrain_reco_pt = 0.0f;
    float m_jetMLTrain_raw_pt = 0.0f;
    float m_jetMLTrain_jet_eta = 0.0f;
    float m_jetMLTrain_jet_phi = 0.0f;
    float m_jetMLTrain_jet_area = 0.0f;
    float m_jetMLTrain_dphi = 0.0f;
    float m_jetMLTrain_rho = 0.0f;
    float m_jetMLTrain_local_rho = 0.0f;
    float m_jetMLTrain_local_et = 0.0f;
    float m_jetMLTrain_truth_pt = -1.0f;
    float m_jetMLTrain_truth_eta = 0.0f;
    float m_jetMLTrain_truth_phi = 0.0f;
    float m_jetMLTrain_match_dr = -1.0f;
    float m_jetMLTrain_delta_pt = 0.0f;
    float m_jetMLTrain_response_ratio = 0.0f;
    
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
    
    // IsolationAudit (local-only, env-gated)
    bool m_isoAuditMode = false;
    bool m_isoAuditTargetReached = false;
    bool m_isoAuditStopAnnounced = false;
    int m_isoAuditTargetPerCent = 200;
    int m_isoAuditRequiredCentBins = 0;
    int m_isoAuditProgressEveryEvents = 1;
    int m_isoAuditSkipSummaryEveryEvents = 1000;
    int m_isoAuditExemplarsPerCell = 3;
    double m_isoAuditLargePositiveGeV = 15.0;
    double m_isoAuditNearZeroAbsGeV = 1.0;
    double m_isoAuditMismatchGeV = 0.1;
    unsigned long long m_isoAuditStopEvent = 0;
    unsigned long long m_isoAuditRunNumber = 0;
    std::string m_isoAuditStopReason;
    
    unsigned long long m_isoAuditNoPhotonContainerEvents = 0;
    unsigned long long m_isoAuditNonFiniteKinematics = 0;
    unsigned long long m_isoAuditOutsideConfiguredCentrality = 0;
    IsoAuditSkipSnapshot m_isoAuditSkipLastSummary{};
    
    IsoAuditEventFlow m_isoAuditFlowGlobal{};
    std::vector<IsoAuditEventFlow> m_isoAuditFlowByCent;
    std::vector<std::vector<IsoAuditCell>> m_isoAuditCells;
    
    // Histogram storage: trigger -> (name -> object*)
    std::map<std::string, HistMap> qaHistogramsByTrigger;
    
    // Per-trigger slice counters printed in End()
    std::map<std::string, std::map<std::string, CatStat>> m_catByTrig;

    bool m_replayFoundationEnabled = false;
    bool m_replayNodesReady = false;
    bool m_collaboratorRecoVertexValid=false;
    bool m_replayWriteFailed = false;
    bool m_the134MultiviewSidecarOnly = false;
    bool m_the134FastExtraction = false;
    std::vector<RJReplayFoundationV1::TruthPhotonMissOccurrenceRow> m_replayTruthPhotonMissOccurrences;
    std::vector<RJReplayFoundationV1::EmbeddedPhotonDiagnosticOccurrenceRow> m_replayEmbeddedPhotonDiagnosticOccurrences;
    std::unique_ptr<RJReplayRuntimeV1::Runtime> m_replayRuntime;
    std::unique_ptr<RJTriggerScalersV1::Writer> m_triggerScalerWriter;
  std::unique_ptr<RJTriggerRunInfoV1::Writer> m_triggerRunInfoWriter;
    std::unique_ptr<RJPhotonTrainingViewV1::Runtime> m_photonTrainingViewRuntime;
};

#endif // RECOILJETS_AuAu_H
