#ifndef CALORECO_PHOTONCLUSTERBUILDER_H
#define CALORECO_PHOTONCLUSTERBUILDER_H

#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerDefs.h>

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
#include <limits>

class PHCompositeNode;
class RawClusterContainer;
class RawCluster;
class PhotonClusterv1;
class TowerInfo;
class TowerInfoContainer;
class RawTowerGeomContainer;
class RawTowerGeom;

namespace TMVA
{
  namespace Experimental
  {
    class RBDT;
  }
}  // namespace TMVA

// Simple builder that wraps existing RawClusters above an energy threshold
// into PhotonClusterv1 objects and stores them in a PhotonClusterContainer node.
class PhotonClusterBuilder : public SubsysReco
{
 public:
  explicit PhotonClusterBuilder(const std::string& name = "PhotonClusterBuilder");
  ~PhotonClusterBuilder() override = default;

    int InitRun(PHCompositeNode* topNode) override;
    int process_event(PHCompositeNode* topNode) override;
    int End(PHCompositeNode* topNode) override;

    void set_input_cluster_node(const std::string& n) { m_input_cluster_node = n; }
    void set_output_photon_node(const std::string& n) { m_output_photon_node = n; }
    void set_ET_threshold(float e) { m_min_cluster_et = e; }
    void set_shower_shape_min_tower_energy(float e) { m_shape_min_tower_E = e; }
    void set_iso_min_tower_energy(float e) { m_iso_min_tower_E = e; }
    void set_enable_ss_3x3_moments(bool enable) { m_enable_ss_3x3_moments = enable; }
    void set_bdt_model_file(const std::string& path) { m_bdt_model_file = path; }
    void set_bdt_feature_list(const std::vector<std::string>& features) { m_bdt_feature_list = features; }
    void set_do_bdt(bool do_bdt) { m_do_bdt = do_bdt; }
    void add_named_bdt_score(const std::string& score_name,
                             const std::string& model_file,
                             const std::vector<std::string>& features,
                             float min_et = std::numeric_limits<float>::quiet_NaN(),
                             float max_et = std::numeric_limits<float>::quiet_NaN(),
                             float max_abs_eta = std::numeric_limits<float>::quiet_NaN());
    const std::vector<std::string>& get_bdt_feature_list() const { return m_bdt_feature_list; }

    void set_use_ppg12_pp_iso_axis(bool use) { m_use_ppg12_pp_iso_axis = use; }
    void set_isolation_source(const std::string& source);
    void set_use_ppg12_topocluster_isolation(bool use) { m_use_ppg12_topocluster_isolation = use; }
    void set_ppg12_topocluster_node(const std::string& n) { m_ppg12_topocluster_node = n; }
    void set_ppg12_topocluster_iso_radius(float radius);
    void set_ppg12_topocluster_exclude_candidate(bool use) { m_ppg12_topocluster_exclude_candidate = use; }
    void set_use_ppg12_pp_sim_truth_vertex(bool use) { m_use_ppg12_pp_sim_truth_vertex = use; }
    void set_use_ppg12_pp_sim_global_mbd_vertex(bool use) { m_use_ppg12_pp_sim_global_mbd_vertex = use; }
    // Canonical core/PPG12 TowerInfo acceptance: use get_isGood() only. The
    // historical method name is retained for source/ABI compatibility.
    void set_use_ppg12_pp_sim_towerinfo_shapes(bool use) { m_use_ppg12_pp_sim_towerinfo_shapes = use; }
    void set_skip_ppg12_edge_clusters(bool skip) { m_skip_ppg12_edge_clusters = skip; }

    void set_vz_cut(bool use, float vz_cm) { m_use_vz_cut = use; m_vz_cut_cm = vz_cm; }
    void set_use_vz_cut(bool use) { m_use_vz_cut = use; }
    void set_vz_cut_cm(float vz_cm) { m_vz_cut_cm = vz_cm; }
    bool get_use_vz_cut() const { return m_use_vz_cut; }
    float get_vz_cut_cm() const { return m_vz_cut_cm; }

    void set_emc_tower_node(const std::string& n) { m_emc_tower_node = n; }
    void set_ihcal_tower_node(const std::string& n) { m_ihcal_tower_node = n; }
    void set_ohcal_tower_node(const std::string& n) { m_ohcal_tower_node = n; }

    // Au+Au mode: use UE-subtracted tower nodes for isolation cone sums
    void set_is_auau(bool isAuAu) { m_is_auau = isAuAu; }
    bool get_is_auau() const { return m_is_auau; }

    // TowerInfo node prefix used to construct UE-subtracted node names:
    //   <prefix>_CEMC_RETOWER_SUB1, <prefix>_HCALIN_SUB1, <prefix>_HCALOUT_SUB1
    void set_tower_node_prefix(const std::string& p) { m_tower_node_prefix = p; }
    const std::string& get_tower_node_prefix() const { return m_tower_node_prefix; }
    // Diagnostic-only alternative to the canonical complete good-TowerInfo
    // CEMC shower-shape grid. Unified production must leave this false.
    void set_use_raw_cluster_towermap_for_cemc_shapes(bool use) { m_use_raw_cluster_towermap_for_cemc_shapes = use; }

 private:
  void CreateNodes(PHCompositeNode* topNode);
  bool calculate_shower_shapes(RawCluster* rc, PhotonClusterv1* photon, float eta, float phi);
  void calculate_bdt_score(PhotonClusterv1* photon);
  void calculate_named_bdt_scores(PhotonClusterv1* photon);
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
  std::vector<int> find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer* geom, TowerInfoContainer* towerContainer, float vertex_z, bool isihcal);
  double deltaR(double eta1, double phi1, double eta2, double phi2);
  float calculate_layer_et(float seed_eta, float seed_phi, float radius, TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, RawTowerDefs::CalorimeterId calo_id, float vertex_z);
  float calculate_ppg12_topocluster_raw_eiso(float seed_eta, float seed_phi, float candidate_et, float& topo_sum_et);
  float calculate_ppg12_topocluster_raw_eiso_for_radius(float seed_eta,
                                                         float seed_phi,
                                                         float candidate_et,
                                                         float radius,
                                                         float& topo_sum_et,
                                                         float* positive_only_sum_et = nullptr);
  void load_cemc_bad_tower_mask();
  bool is_cemc_tower_good(TowerInfo* tower, unsigned int tower_key) const;
    bool m_do_bdt{false};

      struct AuditSnapshot
      {
        unsigned long long evt_seen = 0;
        unsigned long long evt_missing_rawclusters = 0;
        unsigned long long evt_vertex_from_mbd = 0;
        unsigned long long evt_vertex_from_global = 0;
        unsigned long long evt_no_finite_vertex = 0;
        unsigned long long evt_skip_vz = 0;
        unsigned long long evt_missing_ppg12_topocluster = 0;
        unsigned long long evt_zero_input_clusters = 0;
        unsigned long long evt_zero_pass_et_clusters = 0;
        unsigned long long evt_zero_built_after_pass_et = 0;
        unsigned long long evt_built_photons = 0;
        unsigned long long clusters_seen_total = 0;
        unsigned long long clusters_pass_et_total = 0;
        unsigned long long photons_built_total = 0;
        unsigned long long mean_time_default_total = 0;
        unsigned long long nonfinite_eta_total = 0;
        unsigned long long nonfinite_phi_total = 0;
        unsigned long long nonfinite_et_total = 0;
      };

      AuditSnapshot make_audit_snapshot() const;
      void print_audit_summary(bool force);

    std::string m_input_cluster_node{"CLUSTERINFO_CEMC"};
    std::string m_output_photon_node{"PHOTONCLUSTER_CEMC"};
    float m_min_cluster_et{5.0f};
    float m_shape_min_tower_E{0.070f};
    float m_iso_min_tower_E{0.0f};
    std::string m_bdt_model_file{"myBDT_5.root"};
        std::vector<std::string> m_bdt_feature_list;
        struct NamedBDTScoreConfig
        {
          std::string score_name;
          std::string model_file;
          std::vector<std::string> features;
          float min_et = std::numeric_limits<float>::quiet_NaN();
          float max_et = std::numeric_limits<float>::quiet_NaN();
          float max_abs_eta = std::numeric_limits<float>::quiet_NaN();
          std::unique_ptr<TMVA::Experimental::RBDT> model;
        };
        std::vector<NamedBDTScoreConfig> m_named_bdt_scores;
        float m_vertex{std::numeric_limits<float>::quiet_NaN()};
        bool m_use_vz_cut{true};
        float m_vz_cut_cm{30.0f};
        bool m_use_ppg12_pp_iso_axis{false};
        bool m_use_ppg12_topocluster_isolation{false};
        bool m_use_ppg12_pp_sim_truth_vertex{false};
        bool m_use_ppg12_pp_sim_global_mbd_vertex{false};
        bool m_use_ppg12_pp_sim_towerinfo_shapes{false};
        bool m_skip_ppg12_edge_clusters{false};
        bool m_enable_ss_3x3_moments{false};

        RawClusterContainer* m_rawclusters{nullptr};
      RawClusterContainer* m_photon_container{nullptr};
      std::string m_emc_tower_node{"TOWERINFO_CALIB_CEMC"};
      TowerInfoContainer* m_emc_tower_container{nullptr};
      RawTowerGeomContainer* m_geomEM{nullptr};
      bool m_use_raw_cluster_towermap_for_cemc_shapes{false};
      bool m_apply_cemc_bad_tower_mask{false};
      bool m_cemc_bad_tower_mask_loaded{false};
      std::unordered_set<unsigned int> m_cemc_bad_tower_keys;
      std::string m_ihcal_tower_node{"TOWERINFO_CALIB_HCALIN"};
      TowerInfoContainer* m_ihcal_tower_container{nullptr};
      RawTowerGeomContainer* m_geomIH{nullptr};
      std::string m_ohcal_tower_node{"TOWERINFO_CALIB_HCALOUT"};
      TowerInfoContainer* m_ohcal_tower_container{nullptr};
      RawTowerGeomContainer* m_geomOH{nullptr};

      std::string m_ppg12_topocluster_node{"TOPOCLUSTER_ALLCALO"};
      RawClusterContainer* m_ppg12_topocluster_container{nullptr};
      float m_ppg12_topocluster_iso_radius{0.4f};
      bool m_ppg12_topocluster_exclude_candidate{false};

      // Au+Au UE-subtracted tower nodes (used ONLY for isolation sums)
      bool m_is_auau{false};
      std::string m_tower_node_prefix{"TOWERINFO_CALIB"};
      TowerInfoContainer* m_emc_tower_container_iso{nullptr};
      RawTowerGeomContainer* m_geomEM_iso{nullptr};  // retowered CEMC uses HCALIN geometry
      TowerInfoContainer* m_ihcal_tower_container_iso{nullptr};
      TowerInfoContainer* m_ohcal_tower_container_iso{nullptr};

      bool m_iso_audit_mode{false};
      int m_iso_audit_summary_every_events{1000};

      unsigned long long m_evt_seen{0};
      unsigned long long m_evt_missing_rawclusters{0};
      unsigned long long m_evt_vertex_from_mbd{0};
      unsigned long long m_evt_vertex_from_global{0};
      unsigned long long m_evt_no_finite_vertex{0};
      unsigned long long m_evt_skip_vz{0};
      unsigned long long m_evt_missing_ppg12_topocluster{0};
      unsigned long long m_evt_zero_input_clusters{0};
      unsigned long long m_evt_zero_pass_et_clusters{0};
      unsigned long long m_evt_zero_built_after_pass_et{0};
      unsigned long long m_evt_built_photons{0};

      unsigned long long m_clusters_seen_total{0};
      unsigned long long m_clusters_pass_et_total{0};
      unsigned long long m_photons_built_total{0};
      unsigned long long m_mean_time_default_total{0};
      unsigned long long m_nonfinite_eta_total{0};
      unsigned long long m_nonfinite_phi_total{0};
      unsigned long long m_nonfinite_et_total{0};

      AuditSnapshot m_audit_last_summary{};

    std::unique_ptr<TMVA::Experimental::RBDT> m_bdt;
};

#endif  // CALORECO_PHOTONCLUSTERBUILDER_H
