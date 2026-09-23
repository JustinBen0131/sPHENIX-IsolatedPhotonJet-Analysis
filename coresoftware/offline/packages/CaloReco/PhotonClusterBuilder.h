#ifndef CALORECO_PHOTONCLUSTERBUILDER_H
#define CALORECO_PHOTONCLUSTERBUILDER_H

#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerDefs.h>

#include <limits>
#include <memory>
#include <string>
#include <utility>
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
// Candidates require a reconstructed vertex and the broad ET threshold. An
// optional complete-window requirement restricts geometrical shower support.
// Numerical shower validity, isolation and classifier scores never select them.
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

  // Independent isolation constituent policy, in GeV. Zero means NO energy
  // cut (including signed SUB1 energies), unlike a zero shower floor (E>0).
  // Without this setter, retain the legacy raw-shower-floor/SUB1-no-cut policy.
  void set_isolation_min_tower_energy(float e) { m_iso_min_tower_E = e; }

  enum class VertexSource { Mbd, MbdThenGlobal, GlobalMbd };
  void set_vertex_source(VertexSource source) { m_vertex_source = source; }
  // Current event's selected reconstructed z (cm), or NaN. Read-only event
  // witness, including events with no candidates; never a truth substitution.
  float get_vertex_z() const { return m_vertex; }

  // Optional exact encoded-key-to-channel lookup for cluster-owned towers.
  // Off preserves the existing get_tower_at_key behavior for older users.
  void set_use_explicit_tower_channel_lookup(bool enable) { m_use_explicit_tower_channels = enable; }

  // Optional additional views use the same reconstruction arithmetic at a
  // different tower-energy floor. Named values are stored as
  // "shower_<name>_<quantity>"; the primary, unprefixed values and defaults
  // are unchanged. The name must contain only letters, digits or underscores.
  void add_shower_shape_view(const std::string& name, float min_tower_energy);

  // Optional CEMC COG-tower axis, with cluster-axis fallback for missing
  // geometry. An explicit native floor decouples this axis from the primary
  // shower view. NaN uses the primary floor; default false uses cluster axis.
  void set_isolation_axis_cog(bool enable, float native_floor = std::numeric_limits<float>::quiet_NaN())
  {
    m_isolation_axis_cog = enable;
    m_isolation_axis_floor = native_floor;
  }

  // Optional complete native 7x7 eta support (center bins 3..92). This is a
  // geometrical requirement, never an identification cut. NaN uses the primary
  // shower floor; an explicit floor keeps population support independent of it.
  void set_require_complete_shower_window(bool enable, float native_floor = std::numeric_limits<float>::quiet_NaN())
  {
    m_require_complete_shower_window = enable;
    m_candidate_window_floor = native_floor;
  }

  // Also write center-excluded 3x3 second moments. Off by default.
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
  float select_vertex_z(PHCompositeNode* topNode) const;
  bool calculate_shower_shapes(RawCluster* rc, RawCluster* photon,
                               float tower_floor, const std::string& key_prefix,
                               bool include_ancillary);
  void calculate_hcal_shapes(RawCluster* photon, float eta, float phi);
  void calculate_isolation(RawCluster* rc, RawCluster* photon, float eta, float phi, float candidate_et);
  void calculate_bdt_score(RawCluster* photon);
  void calculate_bdt_scores(RawCluster* photon);
  static float resolve_bdt_feature(const RawCluster* photon, const std::string& feature, float vertex_z, float centrality);
  bool calculate_topocluster_iso(float eta, float phi, float candidate_et, float radius, float& isolation);
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
  bool m_use_explicit_tower_channels{false};
  bool m_require_complete_shower_window{false};
  VertexSource m_vertex_source{VertexSource::Mbd};
  float m_isolation_axis_floor{std::numeric_limits<float>::quiet_NaN()};
  float m_candidate_window_floor{std::numeric_limits<float>::quiet_NaN()};
  float m_vertex_cut_max_abs_z{std::numeric_limits<float>::infinity()};

  std::string m_input_cluster_node{"CLUSTERINFO_CEMC"};
  std::string m_output_photon_node{"PHOTONCLUSTER_CEMC"};
  std::string m_topocluster_node{"TOPOCLUSTER_ALLCALO"};
  std::string m_centrality_node{"CentralityInfo"};
  float m_min_cluster_et{5.0f};
  float m_shape_min_tower_E{0.070f};
  float m_iso_min_tower_E{std::numeric_limits<float>::quiet_NaN()};
  std::vector<std::pair<std::string, float>> m_additional_shower_views;
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
