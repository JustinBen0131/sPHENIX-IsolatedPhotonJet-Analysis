#ifndef CALORECO_PHOTONCLUSTERBUILDER_H
#define CALORECO_PHOTONCLUSTERBUILDER_H

#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerDefs.h>

#include <limits>
#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class RawClusterContainer;
class RawCluster;
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
// into PhotonClusterv1 objects and stores them in a RawClusterContainer node.
//
// Candidate retention depends on the cluster transverse energy only. Shower
// shapes, isolation and any classifier score are attached to the retained
// candidate as named shower-shape parameters and never decide retention.
class PhotonClusterBuilder : public SubsysReco
{
 public:
  explicit PhotonClusterBuilder(const std::string& name = "PhotonClusterBuilder");
  ~PhotonClusterBuilder() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void set_input_cluster_node(const std::string& n) { m_input_cluster_node = n; }
  void set_output_photon_node(const std::string& n) { m_output_photon_node = n; }
  void set_ET_threshold(float e) { m_min_cluster_et = e; }
  void set_shower_shape_min_tower_energy(float e) { m_shape_min_tower_E = e; }

  // Also write the second moments restricted to the 3x3 core about the
  // centre of gravity (weta33_cogx, wphi33_cogx), alongside the full-grid
  // moments. Intended backward-compatible when new options are disabled;
  // runtime equivalence has not been established.
  // Optional isolation axis: the CEMC tower at the shower centre of gravity.
  // Default false preserves the cluster-axis layer-cone convention.
  void set_isolation_axis_cog(bool enable) { m_isolation_axis_cog = enable; }

  void set_enable_3x3_moments(bool enable) { m_enable_3x3_moments = enable; }

  // Skip candidate building when |vertex z| exceeds max_abs_z (cm). The
  // event itself is not aborted. Off by default.
  void set_vertex_cut(bool enable, float max_abs_z)
  {
    m_do_vertex_cut = enable;
    m_vertex_cut_max_abs_z = max_abs_z;
  }

  // Legacy single-model interface: one TMVA RBDT model whose score is written
  // as "bdt_score". Kept for existing macros.
  void set_bdt_model_file(const std::string& path) { m_bdt_model_file = path; }
  void set_bdt_feature_list(const std::vector<std::string>& features) { m_bdt_feature_list = features; }
  void set_do_bdt(bool do_bdt) { m_do_bdt = do_bdt; }
  const std::vector<std::string>& get_bdt_feature_list() const { return m_bdt_feature_list; }

  // Generic interface: any number of named TMVA RBDT models, each with its
  // own ordered feature list, written as shower-shape parameter <score_name>.
  // The model file is either a path or, when cdb_key is non-empty, resolved
  // from the conditions database at InitRun. No model is built in by default.
  // Feature names follow resolve_bdt_feature; "centrality" reads the
  // CentralityInfo node when present.
  void add_bdt_model(const std::string& score_name,
                     const std::string& model_file,
                     const std::vector<std::string>& features,
                     const std::string& cdb_key = "");

  void set_do_subtracted_iso(bool do_subtracted_iso) { m_do_subtracted_iso = do_subtracted_iso; }
  void set_do_topocluster_isolation(bool do_topo_iso) { m_do_topocluster_isolation = do_topo_iso; }
  void set_topocluster_node(const std::string& n) { m_topocluster_node = n; }
  void set_centrality_node(const std::string& n) { m_centrality_node = n; }

 private:
  struct BdtModel
  {
    std::string score_name;
    std::string model_file;
    std::string cdb_key;
    std::vector<std::string> features;
    std::unique_ptr<TMVA::Experimental::RBDT> bdt;
  };

  void CreateNodes(PHCompositeNode* topNode);
  bool calculate_shower_shapes(RawCluster* rc, RawCluster* photon, float eta, float phi);
  void calculate_bdt_score(RawCluster* photon);
  void calculate_bdt_scores(RawCluster* photon);
  static float resolve_bdt_feature(const RawCluster* photon, const std::string& feature, float vertex_z, float centrality);
  bool calculate_topocluster_iso(float eta, float phi, float candidate_et, float& iso03, float& iso04);
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
  std::vector<int> find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer* geom, TowerInfoContainer* towerContainer, float vertex_z, bool isihcal);
  double deltaR(double eta1, double phi1, double eta2, double phi2);
  float calculate_layer_et(float seed_eta, float seed_phi, float radius, TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, RawTowerDefs::CalorimeterId calo_id, float vertex_z);
  bool m_do_bdt{false};
  bool m_do_subtracted_iso{false};
  bool m_do_topocluster_isolation{false};
  bool m_enable_3x3_moments{false};
  bool m_isolation_axis_cog{false};
  bool m_do_vertex_cut{false};
  float m_vertex_cut_max_abs_z{std::numeric_limits<float>::infinity()};

  std::string m_input_cluster_node{"CLUSTERINFO_CEMC"};
  std::string m_output_photon_node{"PHOTONCLUSTER_CEMC"};
  std::string m_topocluster_node{"TOPOCLUSTER_ALLCALO"};
  std::string m_centrality_node{"CentralityInfo"};
  float m_min_cluster_et{5.0f};
  float m_shape_min_tower_E{0.070f};
  std::string m_bdt_model_file;  // explicit file required when legacy evaluation is enabled
  std::vector<std::string> m_bdt_feature_list;
  std::vector<BdtModel> m_bdt_models;
  float m_vertex{std::numeric_limits<float>::quiet_NaN()};
  float m_centrality{std::numeric_limits<float>::quiet_NaN()};
  float m_subtracted_iso_defval{-999};
  float m_topo_iso_defval{-999};

  RawClusterContainer* m_rawclusters{nullptr};
  RawClusterContainer* m_photon_container{nullptr};
  RawClusterContainer* m_topocluster_container{nullptr};
  TowerInfoContainer* m_emc_tower_container{nullptr};
  RawTowerGeomContainer* m_geomEM{nullptr};
  TowerInfoContainer* m_ihcal_tower_container{nullptr};
  RawTowerGeomContainer* m_geomIH{nullptr};
  TowerInfoContainer* m_ohcal_tower_container{nullptr};
  RawTowerGeomContainer* m_geomOH{nullptr};
  TowerInfoContainer* m_emc_sub1_tower_container{nullptr};
  TowerInfoContainer* m_ihcal_sub1_tower_container{nullptr};
  TowerInfoContainer* m_ohcal_sub1_tower_container{nullptr};
  std::unique_ptr<TMVA::Experimental::RBDT> m_bdt;
};

#endif  // CALORECO_PHOTONCLUSTERBUILDER_H
