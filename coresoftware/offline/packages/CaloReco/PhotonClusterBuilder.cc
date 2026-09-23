#include "PhotonClusterBuilder.h"

#include <calobase/PhotonClusterv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// Tower stuff
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfo.h>

// for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <centrality/CentralityInfo.h>

#include <ffamodules/CDBInterface.h>

#include <TMVA/RBDT.hxx>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>

namespace
{
  // Helper function to shift tower indices for wrapping in phi
  void shift_tower_index(int& ieta, int& iphi, int etadiv, int phidiv)
  {
    while (iphi < 0)
    {
      iphi += phidiv;
    }
    while (iphi >= phidiv)
    {
      iphi -= phidiv;
    }
    if (ieta < 0 || ieta >= etadiv)
    {
      ieta = -1;  // invalid
    }
  }

  // TowerInfo's encoded key is detector-specific. Resolve cluster-owned
  // towers through the corresponding channel, as in the native production
  // builder, rather than relying on container-specific get_tower_at_key.
  TowerInfo* tower_at_encoded_key(TowerInfoContainer* container, unsigned int key,
                                  RawTowerDefs::CalorimeterId detector,
                                  bool use_explicit_channel)
  {
    if (!container) return nullptr;
    if (!use_explicit_channel) return container->get_tower_at_key(key);
    if (detector == RawTowerDefs::CalorimeterId::CEMC)
      return container->get_tower_at_channel(TowerInfoDefs::decode_emcal(key));
    return container->get_tower_at_channel(TowerInfoDefs::decode_hcal(key));
  }
}  // namespace

PhotonClusterBuilder::PhotonClusterBuilder(const std::string& name)
  : SubsysReco(name)
{
}

PhotonClusterBuilder::~PhotonClusterBuilder() = default;

void PhotonClusterBuilder::add_shower_shape_view(const std::string& name, float min_tower_energy)
{
  if (name.empty() || !std::isfinite(min_tower_energy) || min_tower_energy < 0.0F ||
      !std::all_of(name.begin(), name.end(), [](unsigned char c) { return std::isalnum(c) || c == '_'; }))
  {
    throw std::invalid_argument("PhotonClusterBuilder: invalid additional shower view");
  }
  if (std::any_of(m_additional_shower_views.begin(), m_additional_shower_views.end(),
          [&name](const auto& view) { return view.first == name; }))
  {
    throw std::invalid_argument("PhotonClusterBuilder: duplicate shower view name: " + name);
  }
  m_additional_shower_views.emplace_back(name, min_tower_energy);
}

void PhotonClusterBuilder::add_bdt_model(const std::string& score_name,
                                         const std::string& model_file,
                                         const std::vector<std::string>& features,
                                         const std::string& cdb_key)
{
  if (score_name.empty() || features.empty() || (model_file.empty() == cdb_key.empty()))
  {
    throw std::runtime_error("PhotonClusterBuilder::add_bdt_model: a score name, features, and exactly one file or CDB key are required");
  }
  for (const auto& model : m_bdt_models)
  {
    if (model.score_name == score_name)
    {
      throw std::runtime_error("PhotonClusterBuilder::add_bdt_model: duplicate score name '" + score_name + "'");
    }
  }
  BdtModel model;
  model.score_name = score_name;
  model.model_file = model_file;
  model.cdb_key = cdb_key;
  model.features = features;
  m_bdt_models.push_back(std::move(model));
}

int PhotonClusterBuilder::InitRun(PHCompositeNode* topNode)
{
  if (m_input_cluster_node == m_output_photon_node || !std::isfinite(m_min_cluster_et) || m_min_cluster_et < 0 ||
      !std::isfinite(m_shape_min_tower_E) || m_shape_min_tower_E < 0)
    return Fun4AllReturnCodes::ABORTRUN;
  // NaN is the documented unset/default marker for optional independent floors.
  for (const float floor : {m_iso_min_tower_E, m_isolation_axis_floor, m_candidate_window_floor})
    if (!std::isnan(floor) && (!std::isfinite(floor) || floor < 0.0F))
      return Fun4AllReturnCodes::ABORTRUN;
  if (m_do_vertex_cut && (!std::isfinite(m_vertex_cut_max_abs_z) || m_vertex_cut_max_abs_z <= 0))
    return Fun4AllReturnCodes::ABORTRUN;
  if (m_do_bdt && (m_bdt_model_file.empty() || m_bdt_feature_list.empty()))
    return Fun4AllReturnCodes::ABORTRUN;
  if (m_do_bdt && std::any_of(m_bdt_models.begin(), m_bdt_models.end(),
      [](const BdtModel& model) { return model.score_name == "bdt_score"; }))
    return Fun4AllReturnCodes::ABORTRUN;
  // BDT
  if (m_do_bdt)
  {
    m_bdt = std::make_unique<TMVA::Experimental::RBDT>("myBDT", m_bdt_model_file);
  }

  for (auto& model : m_bdt_models)
  {
    std::string file = model.model_file;
    if (!model.cdb_key.empty())
    {
      file = CDBInterface::instance()->getUrl(model.cdb_key);
      if (file.empty())
      {
        std::cerr << Name() << ": no conditions-database payload for BDT model key '" << model.cdb_key << "'" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
    if (file.empty())
    {
      std::cerr << Name() << ": BDT model '" << model.score_name << "' has neither a file nor a CDB key" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    model.bdt = std::make_unique<TMVA::Experimental::RBDT>("myBDT", file);
  }

  // locate input raw cluster container
  m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
  if (!m_rawclusters)
  {
    std::cerr << Name() << ": could not find RawClusterContainer node '" << m_input_cluster_node << "'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_emc_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (!m_emc_tower_container)
  {
    std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_CEMC'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!m_geomEM)
  {
    std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_CEMC'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ihcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  if (!m_ihcal_tower_container)
  {
    std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_HCALIN'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!m_geomIH)
  {
    std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALIN'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ohcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if (!m_ohcal_tower_container)
  {
    std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_HCALOUT'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!m_geomOH)
  {
    std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALOUT'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_do_subtracted_iso)
  {
    m_emc_sub1_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    m_ihcal_sub1_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
    m_ohcal_sub1_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");

    if (!m_emc_sub1_tower_container || !m_ihcal_sub1_tower_container || !m_ohcal_sub1_tower_container)
    {
      std::cout << Name() << ": subtracted isolation enabled but one or more SUB1 tower nodes are missing; "
                << "iso_sub_* values will remain at " << m_subtracted_iso_defval << std::endl;
    }
  }

  if (m_do_topocluster_isolation)
  {
    m_topocluster_container = findNode::getClass<RawClusterContainer>(topNode, m_topocluster_node);
    if (!m_topocluster_container)
    {
      std::cout << Name() << ": topo-cluster isolation enabled but node '" << m_topocluster_node
                << "' is missing; iso_topo_* values will remain at " << m_topo_iso_defval << std::endl;
    }
  }

  CreateNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonClusterBuilder::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    throw std::runtime_error("PhotonClusterBuilder: DST node not found");
  }
  m_photon_container = findNode::getClass<RawClusterContainer>(dstNode, m_output_photon_node);
  if (!m_photon_container)
  {
    m_photon_container = new RawClusterContainer();
    auto* photonNode = new PHIODataNode<PHObject>(m_photon_container, m_output_photon_node, "PHObject");
    dstNode->addNode(photonNode);
  }
}

// Select a reconstructed vertex once per event. GlobalMbd deliberately has no
// fallback: the MBD component and the aggregate Global vertex are different facts.
float PhotonClusterBuilder::select_vertex_z(PHCompositeNode* topNode) const
{
  const float unavailable = std::numeric_limits<float>::quiet_NaN();
  if (m_vertex_source != VertexSource::GlobalMbd)
  {
    auto* vertices = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
    if (vertices && !vertices->empty())
    {
      const auto* vertex = vertices->begin()->second;
      if (vertex && std::isfinite(vertex->get_z())) return vertex->get_z();
    }
    if (m_vertex_source == VertexSource::Mbd) return unavailable;
  }

  auto* vertices = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertices || vertices->empty() || !vertices->begin()->second) return unavailable;
  const auto* vertex = vertices->begin()->second;
  if (m_vertex_source == VertexSource::MbdThenGlobal)
    return std::isfinite(vertex->get_z()) ? vertex->get_z() : unavailable;

  float z = unavailable;
  for (auto it = vertex->find_vertexes(GlobalVertex::MBD);
       it != vertex->end_vertexes(); ++it)
  {
    if (it->first != GlobalVertex::MBD) continue;
    for (const auto* component : it->second)
      if (component && std::isfinite(component->get_z())) z = component->get_z();
  }
  return z;
}

int PhotonClusterBuilder::process_event(PHCompositeNode* topNode)
{
  if (!m_rawclusters)
  {
    m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
    if (!m_rawclusters)
    {
      std::cerr << Name() << ": missing RawClusterContainer '" << m_input_cluster_node << "'" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Clear output even on an early return; never reuse the previous event.
  m_photon_container->Reset();

  m_vertex = select_vertex_z(topNode);
  if (!std::isfinite(m_vertex)) return Fun4AllReturnCodes::EVENT_OK;

  if (m_do_vertex_cut && !(std::abs(m_vertex) < m_vertex_cut_max_abs_z))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // centrality is an optional model input; NaN when the node is absent
  m_centrality = std::numeric_limits<float>::quiet_NaN();
  if (!m_bdt_models.empty())
  {
    CentralityInfo* centrality = findNode::getClass<CentralityInfo>(topNode, m_centrality_node);
    if (centrality && centrality->has_centile(CentralityInfo::PROP::mbd_NS))
    {
      m_centrality = centrality->get_centile(CentralityInfo::PROP::mbd_NS);
    }
  }

  // iterate over clusters via map to have access to keys if needed
  const auto& rcmap = m_rawclusters->getClustersMap();
  for (const auto& kv : rcmap)
  {
    RawCluster* rc = kv.second;
    if (!rc)
    {
      continue;
    }

    CLHEP::Hep3Vector vertex_vec(0, 0, m_vertex);



    float eta = RawClusterUtility::GetPseudorapidity(*rc, vertex_vec);
    float phi = RawClusterUtility::GetAzimuthAngle(*rc, vertex_vec);
    float E = rc->get_energy();
    float ET = E / std::cosh(eta);
    if (!std::isfinite(ET) || !std::isfinite(eta) || !std::isfinite(phi) || ET < m_min_cluster_et)
    {
      continue;
    }

    // Preserve the configured native geometrical population independently of
    // the primary output view. A model or a numerical validity flag never vetoes.
    if (m_require_complete_shower_window)
    {
      const float floor = std::isnan(m_candidate_window_floor)
          ? m_shape_min_tower_E : m_candidate_window_floor;
      const auto shape = rc->get_shower_shapes(floor);
      if (shape.size() < 12 || !std::isfinite(shape[4]) || !std::isfinite(shape[5]))
        continue;
      const float center_eta = std::floor(shape[4] + 0.5F);
      if (center_eta < 3.0F || center_eta > 92.0F) continue;
    }

    RawCluster* photon = new PhotonClusterv1(*rc);
    // Publish the exact reconstruction kinematics before any optional quantity
    // can be unavailable. Capture must not silently choose another vertex.
    photon->set_shower_shape_parameter("vertex_z", m_vertex);
    photon->set_shower_shape_parameter("cluster_eta", eta);
    photon->set_shower_shape_parameter("cluster_phi", phi);
    photon->set_shower_shape_parameter("cluster_et", ET);
    photon->set_shower_shape_parameter("cluster_pt", ET);

    const bool shower_shapes_valid = calculate_shower_shapes(
        rc, photon, m_shape_min_tower_E, "", true);
    for (const auto& view : m_additional_shower_views)
    {
      calculate_shower_shapes(rc, photon, view.second,
                              "shower_" + view.first + "_", false);
    }
    calculate_hcal_shapes(photon, eta, phi);
    calculate_isolation(rc, photon, eta, phi, ET);
    //this is defensive coding, if do bdt is set false the bdt object should be nullptr
    //and this method will simply pass
    if (m_do_bdt)
    {
      photon->set_shower_shape_parameter("bdt_score", -1.0F);
      if (shower_shapes_valid)
      {
        calculate_bdt_score(photon);
      }
    }
    if (!m_bdt_models.empty())
    {
      for (const auto& model : m_bdt_models)
      {
        if (photon->get_all_shower_shapes().count(model.score_name))
          throw std::runtime_error("model score name collides with a reconstruction field: " + model.score_name);
        photon->set_shower_shape_parameter(model.score_name, std::numeric_limits<float>::quiet_NaN());
      }
      if (shower_shapes_valid)
      {
        calculate_bdt_scores(photon);
      }
    }

    m_photon_container->AddCluster(photon);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonClusterBuilder::calculate_bdt_score(RawCluster* photon)
{
  if (!m_bdt)
  {
    return;
  }

  std::vector<float> x;
  for (const auto& feature : m_bdt_feature_list)
  {
    x.push_back(resolve_bdt_feature(photon, feature, m_vertex, m_centrality));
    //check if the thing we pushed back is NaN
    if (!std::isfinite(x.back()))
    {
      std::cerr << "PhotonClusterBuilder - feature name: " << feature << " is nonfinite" << std::endl;
      photon->set_shower_shape_parameter("bdt_score", std::numeric_limits<float>::quiet_NaN());
      return;
    }
  }

  float bdt_score = -1;  // default value

  const auto prediction = m_bdt->Compute(x);
  if (prediction.size() != 1) throw std::runtime_error("BDT must return one score");
  bdt_score = prediction[0];

  photon->set_shower_shape_parameter("bdt_score", bdt_score);
}

void PhotonClusterBuilder::calculate_bdt_scores(RawCluster* photon)
{
  for (const auto& model : m_bdt_models)
  {
    if (!model.bdt)
    {
      continue;
    }
    std::vector<float> x;
    x.reserve(model.features.size());
    bool complete = true;
    for (const auto& feature : model.features)
    {
      x.push_back(resolve_bdt_feature(photon, feature, m_vertex, m_centrality));
      complete = complete && std::isfinite(x.back());
    }
    // A model is evaluated only on a complete, finite input vector. The
    // score stays NaN otherwise; downstream code distinguishes "not
    // evaluated" from a low score.
    if (!complete)
    {
      continue;
    }
    const auto prediction = model.bdt->Compute(x);
    if (prediction.size() != 1) throw std::runtime_error("BDT must return one score");
    photon->set_shower_shape_parameter(model.score_name, prediction[0]);
  }
}

bool PhotonClusterBuilder::calculate_shower_shapes(RawCluster* rc, RawCluster* photon,
                                                   float tower_floor, const std::string& key_prefix,
                                                   bool include_ancillary)
{
  const auto set_shape = [photon, &key_prefix](const std::string& key, float value)
  {
    photon->set_shower_shape_parameter(key_prefix + key, value);
  };
  set_shape("shower_shape_valid", 0.0F);
  set_shape("shower_shape_floor_gev", tower_floor);

  const auto leadtowerindex = rc->get_lead_tower();
  const int lead_ieta = leadtowerindex.first;
  const int lead_iphi = leadtowerindex.second;

  // for detamax, dphimax, nsaturated
  int detamax = 0;
  int dphimax = 0;
  int nsaturated = 0;
  float clusteravgtime = 0;
  float cluster_total_e = 0;
  unsigned int cluster_time_tower_count = 0;
  const RawCluster::TowerMap& tower_map = rc->get_towermap();
  std::set<unsigned int> towers_in_cluster;
  for (auto tower_iter : tower_map)
  {
    RawTowerDefs::keytype tower_key = tower_iter.first;
    int ieta = RawTowerDefs::decode_index1(tower_key);
    int iphi = RawTowerDefs::decode_index2(tower_key);


    unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
    towers_in_cluster.insert(towerinfokey);
    TowerInfo* towerinfo = tower_at_encoded_key(
        m_emc_tower_container, towerinfokey,
        RawTowerDefs::CalorimeterId::CEMC, m_use_explicit_tower_channels);
    if (towerinfo)
    {
      if (towerinfo->get_isSaturated())
      {
        nsaturated++;
      }
      if (include_ancillary)
      {
        clusteravgtime += towerinfo->get_time() * towerinfo->get_energy();
        cluster_total_e += towerinfo->get_energy();
        ++cluster_time_tower_count;
      }
    }

    int totalphibins = 256;
    auto dphiwrap = [totalphibins](int towerphi, int maxiphi_arg)
    {
      int idphi = towerphi - maxiphi_arg;
      if (idphi > totalphibins / 2)
      {
        idphi -= totalphibins;
      }
      if (idphi < -totalphibins / 2)
      {
        idphi += totalphibins;
      }
      return idphi;
    };

    int deta = ieta - lead_ieta;
    int dphi_val = dphiwrap(iphi, lead_iphi);

    detamax = std::max(std::abs(deta), detamax);
    dphimax = std::max(std::abs(dphi_val), dphimax);
  }

  const float cluster_time_numerator = clusteravgtime;
  if (include_ancillary && cluster_total_e > 0)
  {
    clusteravgtime /= cluster_total_e;
  }
  else if (include_ancillary)
  {
    clusteravgtime = -999.0F;
  }
  if (include_ancillary)
  {
    photon->set_shower_shape_parameter("time_energy_numerator", cluster_time_numerator);
    photon->set_shower_shape_parameter("time_energy_denominator", cluster_total_e);
    photon->set_shower_shape_parameter("time_contributing_towers", static_cast<float>(cluster_time_tower_count));
    photon->set_shower_shape_parameter("time_valid",
        cluster_total_e > 0.0F && std::isfinite(cluster_time_numerator) &&
        std::isfinite(cluster_total_e) && std::isfinite(clusteravgtime) ? 1.0F : 0.0F);
  }

  if (include_ancillary) photon->set_shower_shape_parameter("mean_time", clusteravgtime);

  std::vector<float> showershape = rc->get_shower_shapes(tower_floor);
  if (showershape.size() < 12 || !std::isfinite(showershape[4]) || !std::isfinite(showershape[5]))
  {
    return false;
  }


  float avg_eta = showershape[4] + 0.5F;
  float avg_phi = showershape[5] + 0.5F;


  int maxieta = std::floor(avg_eta);
  int maxiphi = std::floor(avg_phi);


  float E77[7][7] = {{0.0F}};
  int E77_ownership[7][7] = {{0}};

  for (int ieta = maxieta - 3; ieta < maxieta + 4; ieta++)
  {
    for (int iphi = maxiphi - 3; iphi < maxiphi + 4; iphi++)
    {
      //this is defensive coding, if ieta is out of range, set the energy to 0
      //even without this, the requirement for towerinfo object will take care of it
      if (ieta < 0 || ieta > 95)
      {
        E77[ieta - maxieta + 3][iphi - maxiphi + 3] = 0.0F;
        E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 0;
        continue;
      }

      int temp_ieta = ieta;
      int temp_iphi = iphi;
      shift_tower_index(temp_ieta, temp_iphi, 96, 256);
      if (temp_ieta < 0)
      {
        continue;
      }

      unsigned int towerinfokey = TowerInfoDefs::encode_emcal(temp_ieta, temp_iphi);

      if (towers_in_cluster.count(towerinfokey) != 0)
      {
        E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 1;
      }

      TowerInfo* towerinfo = tower_at_encoded_key(
          m_emc_tower_container, towerinfokey,
          RawTowerDefs::CalorimeterId::CEMC, m_use_explicit_tower_channels);
      if (towerinfo && towerinfo->get_isGood())
      {
        float energy = towerinfo->get_energy();
        if (energy > tower_floor)
        {
          E77[ieta - maxieta + 3][iphi - maxiphi + 3] = energy;
        }
      }
    }
  }

  float e11 = E77[3][3];
  float e33 = 0;
  float e55 = 0;
  float e77 = 0;
  float e13 = 0;
  float e15 = 0;
  float e17 = 0;
  float e31 = 0;
  float e51 = 0;
  float e71 = 0;
  float e35 = 0;
  float e37 = 0;
  float e53 = 0;
  float e73 = 0;
  float e57 = 0;
  float e75 = 0;
  float weta = 0;
  float wphi = 0;
  float weta_cog = 0;
  float wphi_cog = 0;
  float weta_cogx = 0;
  float wphi_cogx = 0;
  float weta33_cogx = 0;
  float wphi33_cogx = 0;
  float Eetaphi = 0;
  float Eetaphi33 = 0;
  float shift_eta = avg_eta - std::floor(avg_eta) - 0.5;
  float shift_phi = avg_phi - std::floor(avg_phi) - 0.5;
  float cog_eta = 3 + shift_eta;
  float cog_phi = 3 + shift_phi;

  float w32 = 0;
  float e32 = 0;
  float w52 = 0;
  float e52 = 0;
  float w72 = 0;
  float e72 = 0;
  float detacog = std::abs(maxieta - avg_eta);
  float dphicog = std::abs(maxiphi - avg_phi);
  float drad = std::sqrt(dphicog * dphicog + detacog * detacog);


  int signphi = (avg_phi - std::floor(avg_phi)) > 0.5 ? 1 : -1;

  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < 7; j++)
    {
      int di = std::abs(i - 3);
      int dj = std::abs(j - 3);
      float di_float = i - cog_eta;
      float dj_float = j - cog_phi;

      if (E77_ownership[i][j] == 1)
      {
        weta += E77[i][j] * di * di;
        wphi += E77[i][j] * dj * dj;
        weta_cog += E77[i][j] * di_float * di_float;
        wphi_cog += E77[i][j] * dj_float * dj_float;
        Eetaphi += E77[i][j];
        if (i != 3 || j != 3)
        {
          weta_cogx += E77[i][j] * di_float * di_float;
          wphi_cogx += E77[i][j] * dj_float * dj_float;
        }
        if (di <= 1 && dj <= 1)
        {
          // 3x3 core: same centre of gravity, centre tower excluded from
          // the numerator and included in the denominator.
          Eetaphi33 += E77[i][j];
          if (i != 3 || j != 3)
          {
            weta33_cogx += E77[i][j] * di_float * di_float;
            wphi33_cogx += E77[i][j] * dj_float * dj_float;
          }
        }
      }

      e77 += E77[i][j];
      if (di <= 1 && (dj == 0 || j == (3 + signphi)))
      {
        w32 += E77[i][j] * (i - 3) * (i - 3);
        e32 += E77[i][j];
      }
      if (di <= 2 && (dj == 0 || j == (3 + signphi)))
      {
        w52 += E77[i][j] * (i - 3) * (i - 3);
        e52 += E77[i][j];
      }
      if (di <= 3 && (dj == 0 || j == (3 + signphi)))
      {
        w72 += E77[i][j] * (i - 3) * (i - 3);
        e72 += E77[i][j];
      }

      if (di <= 0 && dj <= 1)
      {
        e13 += E77[i][j];
      }
      if (di <= 0 && dj <= 2)
      {
        e15 += E77[i][j];
      }
      if (di <= 0 && dj <= 3)
      {
        e17 += E77[i][j];
      }
      if (di <= 1 && dj <= 0)
      {
        e31 += E77[i][j];
      }
      if (di <= 2 && dj <= 0)
      {
        e51 += E77[i][j];
      }
      if (di <= 3 && dj <= 0)
      {
        e71 += E77[i][j];
      }
      if (di <= 1 && dj <= 1)
      {
        e33 += E77[i][j];
      }
      if (di <= 1 && dj <= 2)
      {
        e35 += E77[i][j];
      }
      if (di <= 1 && dj <= 3)
      {
        e37 += E77[i][j];
      }
      if (di <= 2 && dj <= 1)
      {
        e53 += E77[i][j];
      }
      if (di <= 3 && dj <= 1)
      {
        e73 += E77[i][j];
      }
      if (di <= 2 && dj <= 2)
      {
        e55 += E77[i][j];
      }
      if (di <= 2 && dj <= 3)
      {
        e57 += E77[i][j];
      }
      if (di <= 3 && dj <= 2)
      {
        e75 += E77[i][j];
      }
    }
  }

  if (Eetaphi > 0)
  {
    weta /= Eetaphi;
    wphi /= Eetaphi;
    weta_cog /= Eetaphi;
    wphi_cog /= Eetaphi;
    weta_cogx /= Eetaphi;
    wphi_cogx /= Eetaphi;
  }
  if (Eetaphi33 > 0)
  {
    weta33_cogx /= Eetaphi33;
    wphi33_cogx /= Eetaphi33;
  }
  if (e32 > 0)
  {
    w32 /= e32;
  }
  if (e52 > 0)
  {
    w52 /= e52;
  }
  if (e72 > 0)
  {
    w72 /= e72;
  }

  set_shape("et1", showershape[0]);
  set_shape("et2", showershape[1]);
  set_shape("et3", showershape[2]);
  set_shape("et4", showershape[3]);
  set_shape("e11", e11);
  set_shape("e22", showershape[8] + showershape[9] + showershape[10] + showershape[11]);
  set_shape("e33", e33);
  set_shape("e55", e55);
  set_shape("e77", e77);
  set_shape("e13", e13);
  set_shape("e15", e15);
  set_shape("e17", e17);
  set_shape("e31", e31);
  set_shape("e51", e51);
  set_shape("e71", e71);
  set_shape("e35", e35);
  set_shape("e37", e37);
  set_shape("e53", e53);
  set_shape("e73", e73);
  set_shape("e57", e57);
  set_shape("e75", e75);
  set_shape("weta", weta);
  set_shape("wphi", wphi);
  set_shape("weta_cog", weta_cog);
  set_shape("wphi_cog", wphi_cog);
  set_shape("weta_cogx", weta_cogx);
  set_shape("wphi_cogx", wphi_cogx);
  if (m_enable_3x3_moments)
  {
    set_shape("weta33_cogx", weta33_cogx);
    set_shape("wphi33_cogx", wphi33_cogx);
    set_shape("moment33_valid", Eetaphi33 > 0.0F && std::isfinite(Eetaphi33) ? 1.0F : 0.0F);
  }
  set_shape("detamax", detamax);
  set_shape("dphimax", dphimax);
  set_shape("nsaturated", nsaturated);
  set_shape("e32", e32);
  set_shape("e52", e52);
  set_shape("e72", e72);
  set_shape("w32", w32);
  set_shape("w52", w52);
  set_shape("w72", w72);
  set_shape("cluster_ietacent", showershape[4]);
  set_shape("cluster_iphicent", showershape[5]);
  set_shape("detacog", detacog);
  set_shape("dphicog", dphicog);
  set_shape("drad", drad);

  set_shape("center_ieta", static_cast<float>(maxieta));
  set_shape("center_iphi", static_cast<float>((maxiphi % 256 + 256) % 256));
  set_shape("cog_eta_local", cog_eta);
  set_shape("cog_phi_local", cog_phi);
  set_shape("shower_shape_valid", Eetaphi > 0.0F &&
      std::isfinite(Eetaphi) &&
      (!m_enable_3x3_moments || (Eetaphi33 > 0.0F && std::isfinite(Eetaphi33))) &&
      std::isfinite(showershape[0]) && std::isfinite(showershape[1]) &&
      std::isfinite(showershape[2]) && std::isfinite(showershape[3]) ? 1.0F : 0.0F);
  return true;
}

// Local HCAL shower witnesses keep their native detector-centered geometry.
void PhotonClusterBuilder::calculate_hcal_shapes(RawCluster* photon, float cluster_eta, float cluster_phi)
{
  // HCAL info
  std::vector<int> ihcal_tower = find_closest_hcal_tower(cluster_eta, cluster_phi, m_geomIH, m_ihcal_tower_container, 0.0, true);
  std::vector<int> ohcal_tower = find_closest_hcal_tower(cluster_eta, cluster_phi, m_geomOH, m_ohcal_tower_container, 0.0, false);

  float ihcal_et = 0;
  float ohcal_et = 0;
  float ihcal_et22 = 0;
  float ohcal_et22 = 0;
  float ihcal_et33 = 0;
  float ohcal_et33 = 0;

  int ihcal_ieta = ihcal_tower[0];
  int ihcal_iphi = ihcal_tower[1];
  float ihcalEt33[3][3] = {{0.0F}};

  int ohcal_ieta = ohcal_tower[0];
  int ohcal_iphi = ohcal_tower[1];
  float ohcalEt33[3][3] = {{0.0F}};

  for (int ieta_h = ihcal_ieta - 1; ieta_h <= ihcal_ieta + 1; ieta_h++)
  {
    for (int iphi_h = ihcal_iphi - 1; iphi_h <= ihcal_iphi + 1; iphi_h++)
    {
      int temp_ieta = ieta_h;
      int temp_iphi = iphi_h;
      shift_tower_index(temp_ieta, temp_iphi, 24, 64);
      if (temp_ieta < 0)
      {
        continue;
      }

      unsigned int towerinfokey = TowerInfoDefs::encode_hcal(temp_ieta, temp_iphi);
      TowerInfo* towerinfo = tower_at_encoded_key(
          m_ihcal_tower_container, towerinfokey,
          RawTowerDefs::CalorimeterId::HCALIN, m_use_explicit_tower_channels);
      if (towerinfo && towerinfo->get_isGood())
      {
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, temp_ieta, temp_iphi);
        RawTowerGeom* tower_geom = m_geomIH->get_tower_geometry(key);
        if (tower_geom)
        {
          float energy = towerinfo->get_energy();
          float eta = getTowerEta(tower_geom, 0, 0, 0);
          float sintheta = 1.0 / std::cosh(eta);
          float Et = energy * sintheta;
          ihcalEt33[ieta_h - ihcal_ieta + 1][iphi_h - ihcal_iphi + 1] = Et;
        }
      }
    }
  }

  for (int ieta_h = ohcal_ieta - 1; ieta_h <= ohcal_ieta + 1; ieta_h++)
  {
    for (int iphi_h = ohcal_iphi - 1; iphi_h <= ohcal_iphi + 1; iphi_h++)
    {
      int temp_ieta = ieta_h;
      int temp_iphi = iphi_h;
      shift_tower_index(temp_ieta, temp_iphi, 24, 64);
      if (temp_ieta < 0)
      {
        continue;
      }

      unsigned int towerinfokey = TowerInfoDefs::encode_hcal(temp_ieta, temp_iphi);
      TowerInfo* towerinfo = tower_at_encoded_key(
          m_ohcal_tower_container, towerinfokey,
          RawTowerDefs::CalorimeterId::HCALOUT, m_use_explicit_tower_channels);
      if (towerinfo && towerinfo->get_isGood())
      {
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, temp_ieta, temp_iphi);
        RawTowerGeom* tower_geom = m_geomOH->get_tower_geometry(key);
        if (tower_geom)
        {
          float energy = towerinfo->get_energy();
          float eta = getTowerEta(tower_geom, 0, 0, 0);
          float sintheta = 1.0 / std::cosh(eta);
          float Et = energy * sintheta;
          ohcalEt33[ieta_h - ohcal_ieta + 1][iphi_h - ohcal_iphi + 1] = Et;
        }
      }
    }
  }

  ihcal_et = ihcalEt33[1][1];
  ohcal_et = ohcalEt33[1][1];

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      ihcal_et33 += ihcalEt33[i][j];
      ohcal_et33 += ohcalEt33[i][j];
      if (i == 1 || j == 1 + ihcal_tower[2])
      {
        if (j == 1 || i == 1 + ihcal_tower[3])
        {
          ihcal_et22 += ihcalEt33[i][j];
        }
      }
      if (i == 1 || j == 1 + ohcal_tower[2])
      {
        if (j == 1 || i == 1 + ohcal_tower[3])
        {
          ohcal_et22 += ohcalEt33[i][j];
        }
      }
    }
  }

  photon->set_shower_shape_parameter("ihcal_et", ihcal_et);
  photon->set_shower_shape_parameter("ohcal_et", ohcal_et);
  photon->set_shower_shape_parameter("ihcal_et22", ihcal_et22);
  photon->set_shower_shape_parameter("ohcal_et22", ohcal_et22);
  photon->set_shower_shape_parameter("ihcal_et33", ihcal_et33);
  photon->set_shower_shape_parameter("ohcal_et33", ohcal_et33);
  photon->set_shower_shape_parameter("ihcal_ieta", ihcal_ieta);
  photon->set_shower_shape_parameter("ihcal_iphi", ihcal_iphi);
  photon->set_shower_shape_parameter("ohcal_ieta", ohcal_ieta);
  photon->set_shower_shape_parameter("ohcal_iphi", ohcal_iphi);

}

// Isolation is evaluated once, independently of whether either shower view is
// usable. The native COG floor is an explicit axis policy, not a shower selection.
void PhotonClusterBuilder::calculate_isolation(RawCluster* rc, RawCluster* photon,
                                                float cluster_eta, float cluster_phi, float ET)
{
  float cog_eta = cluster_eta, cog_phi = cluster_phi;
  // Topo isolation has always used the COG axis. The optional layer-axis
  // switch must not change that established default for other builder users.
  if (m_isolation_axis_cog || m_do_topocluster_isolation)
  {
    const float floor = std::isnan(m_isolation_axis_floor)
        ? m_shape_min_tower_E : m_isolation_axis_floor;
    const auto shape = rc->get_shower_shapes(floor);
    if (shape.size() >= 6 && std::isfinite(shape[4]) && std::isfinite(shape[5]))
    {
      const int ieta = static_cast<int>(std::floor(shape[4] + 0.5F));
      const int iphi = static_cast<int>(std::floor(shape[5] + 0.5F));
      if (ieta >= 0 && ieta < 96)
      {
        const auto key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC,
            ieta, (iphi % 256 + 256) % 256);
        if (auto* geometry = m_geomEM->get_tower_geometry(key))
        {
          cog_eta = getTowerEta(geometry, 0, 0, m_vertex);
          cog_phi = geometry->get_phi();
        }
      }
    }
  }
  const float isolation_eta = m_isolation_axis_cog ? cog_eta : cluster_eta;
  const float isolation_phi = m_isolation_axis_cog ? cog_phi : cluster_phi;
  photon->set_shower_shape_parameter("isolation_axis_eta", isolation_eta);
  photon->set_shower_shape_parameter("isolation_axis_phi", isolation_phi);

  auto compute_layer_iso = [&](RawTowerDefs::CalorimeterId calo_id, float radius)
  {
    TowerInfoContainer* container = nullptr;
    RawTowerGeomContainer* geom = nullptr;
    if (calo_id == RawTowerDefs::CalorimeterId::CEMC)
    {
      container = m_emc_tower_container;
      geom = m_geomEM;
    }
    else if (calo_id == RawTowerDefs::CalorimeterId::HCALIN)
    {
      container = m_ihcal_tower_container;
      geom = m_geomIH;
    }
    else
    {
      container = m_ohcal_tower_container;
      geom = m_geomOH;
    }
    return calculate_layer_et(isolation_eta, isolation_phi, radius, container, geom, calo_id, m_vertex);
  };

  const float emcal_et_04 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.4);
  const float ihcal_et_04 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALIN, 0.4);
  const float ohcal_et_04 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALOUT, 0.4);

  const float emcal_et_03 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.3);
  const float ihcal_et_03 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALIN, 0.3);
  const float ohcal_et_03 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALOUT, 0.3);

  const float emcal_et_02 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.2);
  const float emcal_et_01 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.1);
  const float emcal_et_005 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.05);

  photon->set_shower_shape_parameter("iso_04_emcal", emcal_et_04 - ET);
  photon->set_shower_shape_parameter("iso_04_hcalin", ihcal_et_04);
  photon->set_shower_shape_parameter("iso_04_hcalout", ohcal_et_04);

  photon->set_shower_shape_parameter("iso_03_emcal", emcal_et_03 - ET);
  photon->set_shower_shape_parameter("iso_03_hcalin", ihcal_et_03);
  photon->set_shower_shape_parameter("iso_03_hcalout", ohcal_et_03);
  photon->set_shower_shape_parameter("iso_04_valid",
      std::isfinite(emcal_et_04) && std::isfinite(ihcal_et_04) &&
      std::isfinite(ohcal_et_04) && std::isfinite(ET) ? 1.0F : 0.0F);
  photon->set_shower_shape_parameter("iso_03_valid",
      std::isfinite(emcal_et_03) && std::isfinite(ihcal_et_03) &&
      std::isfinite(ohcal_et_03) && std::isfinite(ET) ? 1.0F : 0.0F);
  photon->set_shower_shape_parameter("iso_02_emcal", emcal_et_02 - ET);
  photon->set_shower_shape_parameter("iso_01_emcal", emcal_et_01 - ET);
  photon->set_shower_shape_parameter("iso_005_emcal", emcal_et_005 - ET);

  photon->set_shower_shape_parameter("iso_sub_04_emcal", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_04_hcalin", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_04_hcalout", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_03_emcal", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_03_hcalin", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_03_hcalout", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_04_valid", 0.0F);
  photon->set_shower_shape_parameter("iso_sub_03_valid", 0.0F);
  photon->set_shower_shape_parameter("iso_sub_02_emcal", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_01_emcal", m_subtracted_iso_defval);
  photon->set_shower_shape_parameter("iso_sub_005_emcal", m_subtracted_iso_defval);

  if (m_do_subtracted_iso && m_emc_sub1_tower_container && m_ihcal_sub1_tower_container && m_ohcal_sub1_tower_container)
  {
    const float sub_emcal_et_04 = calculate_layer_et(isolation_eta, isolation_phi, 0.4, m_emc_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);
    const float sub_ihcal_et_04 = calculate_layer_et(isolation_eta, isolation_phi, 0.4, m_ihcal_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);
    const float sub_ohcal_et_04 = calculate_layer_et(isolation_eta, isolation_phi, 0.4, m_ohcal_sub1_tower_container, m_geomOH, RawTowerDefs::CalorimeterId::HCALOUT, m_vertex);

    const float sub_emcal_et_03 = calculate_layer_et(isolation_eta, isolation_phi, 0.3, m_emc_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);
    const float sub_ihcal_et_03 = calculate_layer_et(isolation_eta, isolation_phi, 0.3, m_ihcal_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);
    const float sub_ohcal_et_03 = calculate_layer_et(isolation_eta, isolation_phi, 0.3, m_ohcal_sub1_tower_container, m_geomOH, RawTowerDefs::CalorimeterId::HCALOUT, m_vertex);

    const float sub_emcal_et_02 = calculate_layer_et(isolation_eta, isolation_phi, 0.2, m_emc_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);
    const float sub_emcal_et_01 = calculate_layer_et(isolation_eta, isolation_phi, 0.1, m_emc_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);
    const float sub_emcal_et_005 = calculate_layer_et(isolation_eta, isolation_phi, 0.05, m_emc_sub1_tower_container, m_geomIH, RawTowerDefs::CalorimeterId::HCALIN, m_vertex);

    photon->set_shower_shape_parameter("iso_sub_04_emcal", sub_emcal_et_04 - ET);
    photon->set_shower_shape_parameter("iso_sub_04_hcalin", sub_ihcal_et_04);
    photon->set_shower_shape_parameter("iso_sub_04_hcalout", sub_ohcal_et_04);
    photon->set_shower_shape_parameter("iso_sub_03_emcal", sub_emcal_et_03 - ET);
    photon->set_shower_shape_parameter("iso_sub_03_hcalin", sub_ihcal_et_03);
    photon->set_shower_shape_parameter("iso_sub_03_hcalout", sub_ohcal_et_03);
    photon->set_shower_shape_parameter("iso_sub_04_valid",
        std::isfinite(sub_emcal_et_04) && std::isfinite(sub_ihcal_et_04) &&
        std::isfinite(sub_ohcal_et_04) && std::isfinite(ET) ? 1.0F : 0.0F);
    photon->set_shower_shape_parameter("iso_sub_03_valid",
        std::isfinite(sub_emcal_et_03) && std::isfinite(sub_ihcal_et_03) &&
        std::isfinite(sub_ohcal_et_03) && std::isfinite(ET) ? 1.0F : 0.0F);
    photon->set_shower_shape_parameter("iso_sub_02_emcal", sub_emcal_et_02 - ET);
    photon->set_shower_shape_parameter("iso_sub_01_emcal", sub_emcal_et_01 - ET);
    photon->set_shower_shape_parameter("iso_sub_005_emcal", sub_emcal_et_005 - ET);
  }

  if (m_do_topocluster_isolation)
  {
    photon->set_shower_shape_parameter("iso_topo_axis_eta", cog_eta);
    photon->set_shower_shape_parameter("iso_topo_axis_phi", cog_phi);
    float iso03 = m_topo_iso_defval, iso04 = m_topo_iso_defval;
    const bool valid03 = calculate_topocluster_iso(cog_eta, cog_phi, ET, 0.3F, iso03);
    const bool valid04 = calculate_topocluster_iso(cog_eta, cog_phi, ET, 0.4F, iso04);
    photon->set_shower_shape_parameter("iso_topo_03", iso03);
    photon->set_shower_shape_parameter("iso_topo_04", iso04);
    photon->set_shower_shape_parameter("iso_topo_03_valid", valid03 ? 1.0F : 0.0F);
    photon->set_shower_shape_parameter("iso_topo_04_valid", valid04 ? 1.0F : 0.0F);
    photon->set_shower_shape_parameter("iso_topo_valid", valid03 && valid04 ? 1.0F : 0.0F);
  }
}

double PhotonClusterBuilder::getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz)
{
  if (!tower_geom)
  {
    return -9999;
  }
  if (vx == 0 && vy == 0 && vz == 0)
  {
    return tower_geom->get_eta();
  }

  double radius = sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) + (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
  double theta = atan2(radius, tower_geom->get_center_z() - vz);
  return -log(tan(theta / 2.));
}

std::vector<int> PhotonClusterBuilder::find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer* geom, TowerInfoContainer* towerContainer, float vertex_z, bool isihcal)
{
  int matchedieta = -1;
  int matchediphi = -1;
  double matchedeta = -999;
  double matchedphi = -999;

  if (!geom || !towerContainer)
  {
    return {-1, -1, 0, 0};
  }

  unsigned int ntowers = towerContainer->size();
  float minR = 999;

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower)
    {
      continue;
    }

    // Detector-explicit mapping also works with persisted containers whose
    // transient detector id was not restored.
    unsigned int towerkey = m_use_explicit_tower_channels
        ? TowerInfoDefs::encode_hcal(channel) : towerContainer->encode_key(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);

    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(isihcal ? RawTowerDefs::CalorimeterId::HCALIN : RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
    RawTowerGeom* tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom)
    {
      continue;
    }

    double this_phi = tower_geom->get_phi();
    double this_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double dR_val = deltaR(eta, phi, this_eta, this_phi);
    if (dR_val < minR)
    {
      minR = dR_val;
      matchedieta = ieta;
      matchediphi = iphi;
      matchedeta = this_eta;
      matchedphi = this_phi;
    }
  }

  float deta = eta - matchedeta;
  float dphi_val = phi - matchedphi;
  if (dphi_val > M_PI)
  {
    dphi_val -= 2 * M_PI;
  }
  if (dphi_val < -M_PI)
  {
    dphi_val += 2 * M_PI;
  }

  int dphisign = (dphi_val > 0) ? 1 : -1;
  int detasign = (deta > 0) ? 1 : -1;

  return {matchedieta, matchediphi, detasign, dphisign};
}

float PhotonClusterBuilder::calculate_layer_et(float seed_eta, float seed_phi, float radius, TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, RawTowerDefs::CalorimeterId calo_id, float vertex_z)
{
  if (!towerContainer || !geomContainer)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float layer_et = 0.0;
  const unsigned int ntowers = towerContainer->size();
  for (unsigned int channel = 0; channel < ntowers; ++channel)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower || !tower->get_isGood())
    {
      continue;
    }

    unsigned int towerkey = calo_id == RawTowerDefs::CalorimeterId::CEMC
        ? TowerInfoDefs::encode_emcal(channel) : TowerInfoDefs::encode_hcal(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);

    RawTowerDefs::keytype geom_key = RawTowerDefs::encode_towerid(calo_id, ieta, iphi);
    RawTowerGeom* tower_geom = geomContainer->get_tower_geometry(geom_key);
    if (!tower_geom)
    {
      continue;
    }

    double tower_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double tower_phi = tower_geom->get_phi();
    if (!std::isfinite(tower_eta) || !std::isfinite(tower_phi)) continue;

    if (deltaR(seed_eta, seed_phi, tower_eta, tower_phi) >= radius)
    {
      continue;
    }

    float energy = tower->get_energy();
    const bool subtracted = towerContainer == m_emc_sub1_tower_container ||
        towerContainer == m_ihcal_sub1_tower_container || towerContainer == m_ohcal_sub1_tower_container;
    if (!std::isfinite(energy)) continue;
    if (std::isnan(m_iso_min_tower_E))
    {
      // Backward-compatible default only. Production explicitly binds its
      // independent isolation policy rather than inheriting a shower floor.
      if (!subtracted && energy <= m_shape_min_tower_E) continue;
    }
    else if (m_iso_min_tower_E > 0.0F && energy <= m_iso_min_tower_E)
    {
      continue;
    }

    float et = energy / std::cosh(tower_eta);
    layer_et += et;
  }

  return layer_et;
}

bool PhotonClusterBuilder::calculate_topocluster_iso(float eta, float phi, float candidate_et,
                                                     float radius, float& isolation)
{
  isolation = m_topo_iso_defval;
  if (!m_topocluster_container || !std::isfinite(eta) || !std::isfinite(phi) ||
      !std::isfinite(candidate_et) || candidate_et <= 0.0F ||
      !std::isfinite(radius) || radius <= 0.0F || !std::isfinite(m_vertex))
    return false;

  double sum = 0.0;
  const CLHEP::Hep3Vector vertex(0, 0, m_vertex);
  const auto range = m_topocluster_container->getClusters();
  for (auto it = range.first; it != range.second; ++it)
  {
    const RawCluster* topo = it->second;
    if (!topo) continue;
    // A negative-energy topo cluster can fail isValid() while supplying a
    // legitimate signed contribution. Test the values used, not that status.
    const double topo_eta = RawClusterUtility::GetPseudorapidity(*topo, vertex);
    const double topo_phi = RawClusterUtility::GetAzimuthAngle(*topo, vertex);
    if (!std::isfinite(topo_eta) || !std::isfinite(topo_phi)) continue;
    const double topo_et = topo->get_energy() / std::cosh(topo_eta);
    if (!std::isfinite(topo_et)) continue;
    if (deltaR(eta, phi, topo_eta, topo_phi) < radius) sum += topo_et;
  }

  // Subtract once in double precision, then cast, as in the native calculation.
  const float result = static_cast<float>(sum - candidate_et);
  if (!std::isfinite(static_cast<float>(sum)) ||
      !std::isfinite(result) || result >= 1.0e8F) return false;
  isolation = result;
  return true;
}

float PhotonClusterBuilder::resolve_bdt_feature(const RawCluster* photon,
                                                const std::string& feature,
                                                float vertex_z,
                                                float centrality)
{
  if (!photon)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  // Exact stored keys are authoritative. In particular, lower-case
  // cluster_* kinematics must not be rewritten before this lookup.
  const auto& shower_shapes = photon->get_all_shower_shapes();
  const auto exact = shower_shapes.find(feature);
  if (exact != shower_shapes.end())
  {
    return exact->second;
  }

  if (feature == "vertex_z" || feature == "vertexz")
  {
    return vertex_z;
  }
  if (feature == "centrality")
  {
    return centrality;
  }
  if (feature == "ET" || feature == "cluster_Et" || feature == "cluster_et")
  {
    const float eta = photon->get_shower_shape_parameter("cluster_eta");
    const float et = photon->get_energy() / std::cosh(eta);
    return std::isfinite(et) ? et : std::numeric_limits<float>::quiet_NaN();
  }
  if (feature == "cluster_Eta")
  {
    return photon->get_shower_shape_parameter("cluster_eta");
  }
  if (feature == "cluster_Phi")
  {
    return photon->get_shower_shape_parameter("cluster_phi");
  }

  const std::string delim = "_over_";
  const auto pos = feature.find(delim);
  if (pos != std::string::npos)
  {
    const std::string numerator_name = feature.substr(0, pos);
    const std::string denominator_name = feature.substr(pos + delim.size());
    const float numerator = resolve_bdt_feature(photon, numerator_name, vertex_z, centrality);
    const float denominator = resolve_bdt_feature(photon, denominator_name, vertex_z, centrality);
    return (denominator > 0) ? (numerator / denominator) : 0.0F;
  }

  if (feature.rfind("cluster_", 0) == 0)
  {
    return photon->get_shower_shape_parameter(feature.substr(8));
  }

  return photon->get_shower_shape_parameter(feature);
}

double PhotonClusterBuilder::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double dphi = phi1 - phi2;
  while (dphi > M_PI)
  {
    dphi -= 2 * M_PI;
  }
  while (dphi <= -M_PI)
  {
    dphi += 2 * M_PI;
  }
  return sqrt(pow(eta1 - eta2, 2) + pow(dphi, 2));
}
