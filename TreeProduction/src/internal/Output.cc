/**
 * @file Output.cc
 * Write one eventTree directly from bounded in-memory event records.
 * The shared writer owns ROOT storage; capture owners supply measured values.
 * Original source identity and new execution provenance are separate metadata.
 * Unimplemented applicable fields remain listed in native_capture and block
 * delivery; default storage values are never evidence of a successful capture.
 */
#include "../PhotonJetTree.h"
#include "Production.h"
#include "../../scripts/eventtree/writer.hxx"
#include <nlohmann/json.hpp>
#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <tuple>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <ctime>

namespace {
std::set<std::string> outputRules();
const std::set<std::string> sharedNominal={"pho_bdt_score","pho_bdt_state","pho_id_class","pho_iso_class"};
}

namespace photonjet {
class DirectEventWriter {
 public:
  std::vector<eventtree::FieldSpec> specs;
  eventtree::Entry entry;
  eventtree::Writer writer;
  std::set<std::string> covered, missing;
  nlohmann::json meta;
  nlohmann::json codes=nlohmann::json::object();
  nlohmann::json policies=nlohmann::json::object();
  std::uint64_t entries = 0;
  std::unique_ptr<TH1D> centrality;
  std::int64_t nextEntry;
  static std::vector<eventtree::FieldSpec> fields(const Config& c) {
    std::vector<eventtree::FieldSpec> result;
    for(const auto& f:c.output.fields) result.push_back({f.name,f.type,f.title,f.counter,f.length});
    if(result.empty())throw std::runtime_error("eventTree schema is not bound to the production plan");
    return result;
  }
  explicit DirectEventWriter(const Config& c):specs(fields(c)),entry(specs),
      writer(c.outputFile,"eventTree","RJEventTreeV1",404,0),nextEntry(c.source.firstEntry) {
    const auto rules=outputRules();
    for(const auto& field:specs)
      (rules.count(field.name)?covered:missing).insert(field.name);
    if(c.system==CollisionSystem::AuAu && c.dataKind==DataKind::Simulation) {
      centrality=std::make_unique<TH1D>("h_centrality",";Centrality [%];Events",100,0.,100.);
      centrality->SetDirectory(nullptr);centrality->Sumw2();
    }
  }
  template<class T> void put(const std::string& n,T v) {entry.put(n,v);}
  template<class T> void all(const std::string& n,const T& v) {entry.putAll(n,v);}
  void empty(const std::string& n) {entry.empty(n);}
  void begin() {entry.reset();}
  // One dictionary per coded branch (meta/codes/<branch>), the shared
  // consolidation and check contract; -1 encodes "no value" (a valid state).
  int code(const std::string& dictionary,const std::string& value) {
    if(value.empty())return -1;
    auto& values=codes[dictionary];if(values.is_null())values=nlohmann::json::array();
    for(std::size_t i=0;i<values.size();++i)if(values[i]==value)return static_cast<int>(i);
    values.push_back(value);return static_cast<int>(values.size()-1);
  }
  void finish(bool rejected=false,bool bookOnly=false) {
    // Empty collections still have a filling rule. Missing rules remain
    // distinct from an observed empty collection in the terminal metadata.
    for(const auto& f:specs) if(f.length==-1 && covered.count(f.name)) entry.empty(f.name);
    for(const auto& name:entry.missing()) {
      const auto& f=*std::find_if(specs.begin(),specs.end(),[&](const auto& s){return s.name==name;});
      if(!rejected&&!bookOnly&&covered.count(name)&&f.length!=-1)
        throw std::runtime_error("declared native scalar rule did not fill field: "+name);
      if(f.length==-1)entry.empty(name);
      else for(int i=0;i<std::max(1,f.length);++i) {
        if(f.type=="double"||f.type=="float")entry.put(name,std::numeric_limits<double>::quiet_NaN());
        else entry.put(name,0);
      }
    }
    entry.write(writer,bookOnly); if(!bookOnly)++entries;
  }
};
}

namespace {
using json=nlohmann::json;
constexpr double missingValue=std::numeric_limits<double>::quiet_NaN();
using Key=std::pair<std::uint64_t,std::uint64_t>;
Key key(const photonjet::Identity& i){return {i.hi,i.lo};}
template<class Row,class Member> std::map<Key,int> indices(const std::vector<Row>& rows,Member member) {
  std::map<Key,int> out;
  for(std::size_t i=0;i<rows.size();++i)if(!out.emplace(key(rows[i].*member),static_cast<int>(i)).second)
    throw std::runtime_error("duplicate event-local object identity");
  return out;
}
int local(const std::map<Key,int>& objects,const photonjet::Identity& id,bool optional=false) {
  if(!id.valid()&&optional)return -1;
  const auto it=objects.find(key(id));
  if(it==objects.end())throw std::runtime_error("dangling event-local object reference");
  return it->second;
}
json bits(double value){std::uint64_t raw;std::memcpy(&raw,&value,sizeof(raw));return {{"f64",raw}};}
double feature(const std::string& name,const photonjet::PhotonRecord& p,
               const photonjet::EventRecord& e,const photonjet::ShowerShapeRecord* v) {
  double value=missingValue;
  if(name=="cluster_Et")value=p.et;
  else if(name=="cluster_Eta")value=p.eta;
  else if(name=="vertexz")value=e.recoVertexZ;
  else if(name=="centrality")value=e.analysisCentrality;
  else if(v) {
    const std::map<std::string,double> values={
      {"cluster_weta_cogx",v->wetaCogX},{"cluster_wphi_cogx",v->wphiCogX},
      {"cluster_weta33_cogx",v->weta33CogX},{"cluster_wphi33_cogx",v->wphi33CogX},
      {"e11_over_e33",v->e11OverE33},{"e32_over_e35",v->e32OverE35},
      {"cluster_et1",v->et1},{"cluster_et2",v->et2},{"cluster_et3",v->et3},{"cluster_et4",v->et4}};
    value=values.at(name);
  } else if(name=="e11_over_e33"||name=="e32_over_e35") {
    const int a=name=="e11_over_e33"?10:12,b=a+1;
    value=std::isfinite(p.primaryShapes[a])&&std::isfinite(p.primaryShapes[b])&&p.primaryShapes[b]>0?
        p.primaryShapes[a]/p.primaryShapes[b]:missingValue;
  } else {
    const std::string suffix=name.substr(8);
    const auto it=std::find(photonjet::PhotonRecord::primaryShapeNames.begin(),
                          photonjet::PhotonRecord::primaryShapeNames.end(),suffix);
    if(it==photonjet::PhotonRecord::primaryShapeNames.end())throw std::runtime_error("unknown frozen model feature: "+name);
    value=p.primaryShapes[std::distance(photonjet::PhotonRecord::primaryShapeNames.begin(),it)];
  }
  return static_cast<float>(value);  // exact evaluator input, widened for storage
}
template<class T>void appendLine(std::ostringstream& out,const std::string& name,const T& value){out<<name<<'='<<value<<'\n';}
}

namespace {
std::set<std::string> outputRules() {
  // Static native rules remain independent of observed event populations.
  return {
    "calo_available_mask",
    "calo_capture_version",
    "calo_emcal",
    "calo_emcal_native",
    "calo_emcal_reason",
    "calo_emcal_state",
    "calo_ihcal",
    "calo_ihcal_native",
    "calo_ihcal_reason",
    "calo_ihcal_state",
    "calo_max_cluster_e",
    "calo_native_mask",
    "calo_ohcal",
    "calo_ohcal_native",
    "calo_ohcal_reason",
    "calo_ohcal_state",
    "calo_require_isgood",
    "calo_require_isgood_native",
    "calo_total",
    "calo_total_native",
    "calo_total_reason",
    "calo_total_state",
    "calo_valid_mask",
    "cell_denominator_member",
    "cell_e",
    "cell_floor0_member",
    "cell_floor70_member",
    "cell_grid_bits",
    "cell_ieta",
    "cell_iphi",
    "cell_is_good",
    "cell_is_negative",
    "cell_is_nonfinite",
    "cell_is_zero",
    "cell_local_ieta",
    "cell_local_iphi",
    "cell_pho",
    "cell_rawcluster_owned",
    "cell_rawcluster_value",
    "cell_rawcluster_value_present",
    "cell_seed_state",
    "cell_time_available",
    "cell_time_finite",
    "cell_time_samples",
    "cell_time_version",
    "cell_tower_status",
    "cent_inputs_valid",
    "cent_mb_decision",
    "cent_mbd_clock",
    "cent_mbd_event",
    "cent_mbd_femclock",
    "cent_mbd_z",
    "cent_mbd_z_valid",
    "cent_native_bin",
    "cent_native_centile",
    "cent_native_valid",
    "cent_pmt_available",
    "cent_pmt_selected",
    "cent_pmt_time",
    "cent_pmt_valid",
    "cent_replay_version",
    "cent_selected_charge",
    "cent_value",
    "evt_aug_donor_evidence",
    "evt_collab_iface_version",
    "evt_embedded_mb_decision",
    "evt_n_candidates",
    "evt_n_recoils",
    "evt_n_tags",
    "evt_nodes_ready",
    "evt_phys_seq",
    "evt_phys_seq_valid",
    "evt_record",
    "evt_run",
    "evt_seq",
    "evt_src_entry",
    "evt_src_file_ordinal",
    "evt_src_global_entry",
    "evt_terminal_status",
    "gl1_bco",
    "gl1_discontinuity",
    "gl1_live",
    "gl1_raw",
    "gl1_scaled",
    "gl1_state",
    "isoc_calibrated_energy",
    "isoc_candidate_removal_state",
    "isoc_delta_eta",
    "isoc_delta_phi",
    "isoc_delta_r",
    "isoc_mask_state",
    "isoc_native_key",
    "isoc_native_source",
    "isoc_pho",
    "isoc_phosub_residual",
    "isoc_quality_state",
    "isoc_raw_energy",
    "isoc_sub1_energy",
    "isoc_subsystem",
    "isow_cone_sum",
    "isow_constituent_count",
    "isow_pass_state",
    "isow_reconstructed_or_truth",
    "isow_sideband_threshold",
    "isow_subtraction_method",
    "isow_threshold",
    "jcon_constituent_ordinal",
    "jcon_energy",
    "jcon_eta",
    "jcon_native_key",
    "jcon_native_source",
    "jcon_phi",
    "jcon_quality_state",
    "jcon_subsystem",
    "jet02_area",
    "jet02_emcal_frac",
    "jet02_eta",
    "jet02_ihcal_frac",
    "jet02_jes_valid",
    "jet02_ohcal_frac",
    "jet02_phi",
    "jet02_pt",
    "jet02_rawpt",
    "jet03_area",
    "jet03_emcal_frac",
    "jet03_eta",
    "jet03_ihcal_frac",
    "jet03_jes_valid",
    "jet03_ohcal_frac",
    "jet03_phi",
    "jet03_pt",
    "jet03_rawpt",
    "jet04_area",
    "jet04_emcal_frac",
    "jet04_eta",
    "jet04_ihcal_frac",
    "jet04_jes_valid",
    "jet04_ohcal_frac",
    "jet04_phi",
    "jet04_pt",
    "jet04_rawpt",
    "jetunsub02_area",
    "jetunsub02_eta",
    "jetunsub02_phi",
    "jetunsub02_pt",
    "jetunsub02_rawpt",
    "jetunsub03_area",
    "jetunsub03_eta",
    "jetunsub03_phi",
    "jetunsub03_pt",
    "jetunsub03_rawpt",
    "jetunsub04_area",
    "jetunsub04_eta",
    "jetunsub04_phi",
    "jetunsub04_pt",
    "jetunsub04_rawpt",
    "jetview_available",
    "jlink_class",
    "jlink_metric",
    "jlink_reco",
    "jlink_truth",
    "mbd_charge_sum",
    "mbd_charge_sum_origin",
    "mbd_charge_sum_state",
    "mbd_native_available",
    "mbd_native_capture_version",
    "mbd_native_north_npmt",
    "mbd_native_north_time_finite",
    "mbd_native_north_time_ns",
    "mbd_native_south_time_finite",
    "mbd_native_south_time_ns",
    "mbd_native_t0_finite",
    "mbd_native_valid",
    "mbd_npmt_north",
    "mbd_npmt_south",
    "mbd_pmt_available",
    "mbd_pmt_q",
    "mbd_t0",
    "mbd_t0_reason",
    "mbd_t0_state",
    "mbd_time_north",
    "mbd_time_north_reason",
    "mbd_time_north_state",
    "mbd_time_south",
    "mbd_time_south_reason",
    "mbd_time_south_state",
    "mbd_total_charge_native",
    "pair_delta_phi",
    "pair_jet",
    "pair_jet_rank",
    "pair_pho",
    "pair_photon_rank",
    "pair_recoil_state",
    "pair_wrong_photon_class",
    "pair_wrong_recoil_class",
    "pair_xjgamma",
    "pho_assoc_dr",
    "pho_assoc_state",
    "pho_assoc_truthpho",
    "pho_auau_npb_score_native",
    "pho_aug_donor_evidence",
    "pho_aunpbcap_binding_policy",
    "pho_aunpbcap_capture_version",
    "pho_aunpbcap_configured_cut",
    "pho_aunpbcap_domain_mode",
    "pho_aunpbcap_evaluation_state",
    "pho_aunpbcap_in_producer_domain",
    "pho_aunpbcap_input_et_bits",
    "pho_aunpbcap_input_eta_bits",
    "pho_aunpbcap_input_finite",
    "pho_aunpbcap_input_order_is_npb25",
    "pho_aunpbcap_inputs",
    "pho_aunpbcap_inputs_captured",
    "pho_aunpbcap_inputs_finite",
    "pho_aunpbcap_kinematic_input_bits_available",
    "pho_aunpbcap_max_abs_eta",
    "pho_aunpbcap_max_et",
    "pho_aunpbcap_min_et",
    "pho_aunpbcap_model_available",
    "pho_aunpbcap_raw_score",
    "pho_aunpbcap_source_cluster_et",
    "pho_aunpbcap_source_encounter_ordinal",
    "pho_aunpbcap_source_eta",
    "pho_aunpbcap_source_identity_available",
    "pho_aunpbcap_source_native_cluster_key",
    "pho_aunpbcap_source_phi",
    "pho_aunpbcap_source_photon_id",
    "pho_below15_state",
    "pho_dom_barcode",
    "pho_dom_embedding_id",
    "pho_dom_energy",
    "pho_dom_evaluator_mode",
    "pho_dom_pid",
    "pho_dom_state",
    "pho_dom_track_id",
    "pho_dom_vertex_id",
    "pho_e11_over_e33",
    "pho_e32_over_e35",
    "pho_encounter",
    "pho_et",
    "pho_et1",
    "pho_eta",
    "pho_finite_feature_state",
    "pho_key",
    "pho_npb_applicability",
    "pho_npb_cluster_et",
    "pho_npb_cluster_et1",
    "pho_npb_cluster_et1_evidence_index",
    "pho_npb_cluster_et1_origin_index",
    "pho_npb_cluster_et1_reason_index",
    "pho_npb_cluster_et1_state",
    "pho_npb_cluster_et2",
    "pho_npb_cluster_et2_evidence_index",
    "pho_npb_cluster_et2_origin_index",
    "pho_npb_cluster_et2_reason_index",
    "pho_npb_cluster_et2_state",
    "pho_npb_cluster_et3",
    "pho_npb_cluster_et3_evidence_index",
    "pho_npb_cluster_et3_origin_index",
    "pho_npb_cluster_et3_reason_index",
    "pho_npb_cluster_et3_state",
    "pho_npb_cluster_et4",
    "pho_npb_cluster_et4_evidence_index",
    "pho_npb_cluster_et4_origin_index",
    "pho_npb_cluster_et4_reason_index",
    "pho_npb_cluster_et4_state",
    "pho_npb_cluster_et_evidence_index",
    "pho_npb_cluster_et_origin_index",
    "pho_npb_cluster_et_reason_index",
    "pho_npb_cluster_et_state",
    "pho_npb_cluster_eta",
    "pho_npb_cluster_eta_evidence_index",
    "pho_npb_cluster_eta_origin_index",
    "pho_npb_cluster_eta_reason_index",
    "pho_npb_cluster_eta_state",
    "pho_npb_cluster_w32",
    "pho_npb_cluster_w32_evidence_index",
    "pho_npb_cluster_w32_origin_index",
    "pho_npb_cluster_w32_reason_index",
    "pho_npb_cluster_w32_state",
    "pho_npb_cluster_w52",
    "pho_npb_cluster_w52_evidence_index",
    "pho_npb_cluster_w52_origin_index",
    "pho_npb_cluster_w52_reason_index",
    "pho_npb_cluster_w52_state",
    "pho_npb_cluster_w72",
    "pho_npb_cluster_w72_evidence_index",
    "pho_npb_cluster_w72_origin_index",
    "pho_npb_cluster_w72_reason_index",
    "pho_npb_cluster_w72_state",
    "pho_npb_cluster_weta_cogx",
    "pho_npb_cluster_weta_cogx_evidence_index",
    "pho_npb_cluster_weta_cogx_origin_index",
    "pho_npb_cluster_weta_cogx_reason_index",
    "pho_npb_cluster_weta_cogx_state",
    "pho_npb_cluster_wphi_cogx",
    "pho_npb_cluster_wphi_cogx_evidence_index",
    "pho_npb_cluster_wphi_cogx_origin_index",
    "pho_npb_cluster_wphi_cogx_reason_index",
    "pho_npb_cluster_wphi_cogx_state",
    "pho_npb_cut",
    "pho_npb_e11_over_e13",
    "pho_npb_e11_over_e13_evidence_index",
    "pho_npb_e11_over_e13_origin_index",
    "pho_npb_e11_over_e13_reason_index",
    "pho_npb_e11_over_e13_state",
    "pho_npb_e11_over_e15",
    "pho_npb_e11_over_e15_evidence_index",
    "pho_npb_e11_over_e15_origin_index",
    "pho_npb_e11_over_e15_reason_index",
    "pho_npb_e11_over_e15_state",
    "pho_npb_e11_over_e17",
    "pho_npb_e11_over_e17_evidence_index",
    "pho_npb_e11_over_e17_origin_index",
    "pho_npb_e11_over_e17_reason_index",
    "pho_npb_e11_over_e17_state",
    "pho_npb_e11_over_e22",
    "pho_npb_e11_over_e22_evidence_index",
    "pho_npb_e11_over_e22_origin_index",
    "pho_npb_e11_over_e22_reason_index",
    "pho_npb_e11_over_e22_state",
    "pho_npb_e11_over_e31",
    "pho_npb_e11_over_e31_evidence_index",
    "pho_npb_e11_over_e31_origin_index",
    "pho_npb_e11_over_e31_reason_index",
    "pho_npb_e11_over_e31_state",
    "pho_npb_e11_over_e33",
    "pho_npb_e11_over_e33_evidence_index",
    "pho_npb_e11_over_e33_origin_index",
    "pho_npb_e11_over_e33_reason_index",
    "pho_npb_e11_over_e33_state",
    "pho_npb_e11_over_e51",
    "pho_npb_e11_over_e51_evidence_index",
    "pho_npb_e11_over_e51_origin_index",
    "pho_npb_e11_over_e51_reason_index",
    "pho_npb_e11_over_e51_state",
    "pho_npb_e11_over_e71",
    "pho_npb_e11_over_e71_evidence_index",
    "pho_npb_e11_over_e71_origin_index",
    "pho_npb_e11_over_e71_reason_index",
    "pho_npb_e11_over_e71_state",
    "pho_npb_e22_over_e33",
    "pho_npb_e22_over_e33_evidence_index",
    "pho_npb_e22_over_e33_origin_index",
    "pho_npb_e22_over_e33_reason_index",
    "pho_npb_e22_over_e33_state",
    "pho_npb_e22_over_e35",
    "pho_npb_e22_over_e35_evidence_index",
    "pho_npb_e22_over_e35_origin_index",
    "pho_npb_e22_over_e35_reason_index",
    "pho_npb_e22_over_e35_state",
    "pho_npb_e22_over_e37",
    "pho_npb_e22_over_e37_evidence_index",
    "pho_npb_e22_over_e37_origin_index",
    "pho_npb_e22_over_e37_reason_index",
    "pho_npb_e22_over_e37_state",
    "pho_npb_e22_over_e53",
    "pho_npb_e22_over_e53_evidence_index",
    "pho_npb_e22_over_e53_origin_index",
    "pho_npb_e22_over_e53_reason_index",
    "pho_npb_e22_over_e53_state",
    "pho_npb_e32_over_e35",
    "pho_npb_e32_over_e35_evidence_index",
    "pho_npb_e32_over_e35_origin_index",
    "pho_npb_e32_over_e35_reason_index",
    "pho_npb_e32_over_e35_state",
    "pho_npb_evaluation",
    "pho_npb_pass",
    "pho_npb_pass_reason",
    "pho_npb_pass_state",
    "pho_npb_score",
    "pho_npb_score_native",
    "pho_npb_score_reason",
    "pho_npb_score_state",
    "pho_npb_scoring_encounter",
    "pho_npb_scoring_et_bits",
    "pho_npb_scoring_eta",
    "pho_npb_scoring_eta_bits",
    "pho_npb_scoring_evidence",
    "pho_npb_scoring_key",
    "pho_npb_scoring_node",
    "pho_npb_scoring_phi",
    "pho_npb_scoring_policy",
    "pho_npb_scoring_pt",
    "pho_npb_vertexz",
    "pho_npb_vertexz_evidence_index",
    "pho_npb_vertexz_origin_index",
    "pho_npb_vertexz_reason_index",
    "pho_npb_vertexz_state",
    "pho_npbcap_association",
    "pho_npbcap_binding_policy",
    "pho_npbcap_capture_version",
    "pho_npbcap_configured_cut",
    "pho_npbcap_domain_mode",
    "pho_npbcap_evaluation_state",
    "pho_npbcap_in_producer_domain",
    "pho_npbcap_input_et_bits",
    "pho_npbcap_input_eta_bits",
    "pho_npbcap_input_finite",
    "pho_npbcap_input_order_is_npb25",
    "pho_npbcap_inputs",
    "pho_npbcap_inputs_captured",
    "pho_npbcap_inputs_finite",
    "pho_npbcap_kinematic_input_bits_available",
    "pho_npbcap_max_abs_eta",
    "pho_npbcap_max_et",
    "pho_npbcap_min_et",
    "pho_npbcap_model_available",
    "pho_npbcap_raw_score",
    "pho_npbcap_source_cluster_et",
    "pho_npbcap_source_encounter_ordinal",
    "pho_npbcap_source_eta",
    "pho_npbcap_source_identity_available",
    "pho_npbcap_source_native_cluster_key",
    "pho_npbcap_source_phi",
    "pho_npbcap_source_photon_id",
    "pho_ntime_den",
    "pho_ntime_finite",
    "pho_ntime_ntowers",
    "pho_ntime_num",
    "pho_ntime_version",
    "pho_phi",
    "pho_presel_active",
    "pho_presel_bits",
    "pho_presel_reference",
    "pho_prodbdt_applicability",
    "pho_prodbdt_delta_wp70",
    "pho_prodbdt_delta_wp80",
    "pho_prodbdt_delta_wp90",
    "pho_prodbdt_finite",
    "pho_prodbdt_inputs",
    "pho_prodbdt_score",
    "pho_prodbdt_wp70",
    "pho_prodbdt_wp80",
    "pho_prodbdt_wp90",
    "pho_producer_vtx_z",
    "pho_rank_keys",
    "pho_resp_class",
    "pho_time_den_e",
    "pho_time_den_e_reason",
    "pho_time_den_e_state",
    "pho_time_ns",
    "pho_time_ntowers",
    "pho_time_raw",
    "pho_time_raw_reason",
    "pho_time_raw_state",
    "pho_truth_match_barcode",
    "pho_truth_match_state",
    "pho_weta33_cogx",
    "pho_weta35_cogx",
    "pho_weta_cogx",
    "pho_wphi33_cogx",
    "pho_wphi53_cogx",
    "pho_wphi_cogx",
    "plink_class",
    "plink_metric",
    "plink_policy",
    "plink_reco",
    "plink_truth",
    "ppg12_weight_active_state",
    "ppg12_weight_lane_code",
    "recopho_capture_domain",
    "recopho_capture_state",
    "recopho_capture_version",
    "recopho_eligible",
    "recopho_inspected",
    "recopho_unclassifiable",
    "recopho_written",
    "shv_active_moment_cell_count",
    "shv_active_sum_cell_count",
    "shv_center_eta_index",
    "shv_center_phi_index",
    "shv_cog_eta",
    "shv_cog_phi",
    "shv_e11",
    "shv_e11_over_e33",
    "shv_e32",
    "shv_e32_over_e35",
    "shv_e33",
    "shv_e35",
    "shv_energy_source",
    "shv_et1",
    "shv_et2",
    "shv_et3",
    "shv_et4",
    "shv_exact_zero_count",
    "shv_features",
    "shv_finite_feature_state",
    "shv_good_cell_count",
    "shv_moment33_denominator",
    "shv_moment33_eta_numerator",
    "shv_moment33_phi_numerator",
    "shv_moment_denominator",
    "shv_moment_eta_numerator",
    "shv_moment_membership",
    "shv_moment_phi_numerator",
    "shv_negative_count",
    "shv_nonfinite_count",
    "shv_owned_cell_count",
    "shv_raw_center_eta",
    "shv_raw_center_phi",
    "shv_rectangular_membership",
    "shv_weta33_cogx",
    "shv_weta_cogx",
    "shv_wphi33_cogx",
    "shv_wphi_cogx",
    "snap_quota_class",
    "snap_selection_reason",
    "snap_serialized_payload_hash",
    "trig_capture_version",
    "trig_decision_available",
    "trig_decision_state",
    "trig_live",
    "trig_packet_status",
    "trig_packet_version",
    "trig_scaled",
    "truth_denominator_complete",
    "truth_mbvtx_z",
    "truth_vtx_z",
    "truthcen_analysis_signal_count",
    "truthcen_capture_state",
    "truthcen_duplicate_track_count",
    "truthcen_isolation_input_incomplete_count",
    "truthcen_native_embedded_primary_pid22_count",
    "truthcen_rejected_count",
    "truthcen_rejected_invalid_kinematics_count",
    "truthcen_rejected_nonfinite_isolation_count",
    "truthcen_serialized_raw_count",
    "truthjet02_eta",
    "truthjet02_guard_state",
    "truthjet02_phi",
    "truthjet02_pt",
    "truthjet03_eta",
    "truthjet03_guard_state",
    "truthjet03_phi",
    "truthjet03_pt",
    "truthjet04_eta",
    "truthjet04_guard_state",
    "truthjet04_phi",
    "truthjet04_pt",
    "truthjetcap_container_valid",
    "truthjetcap_state",
    "truthpho_analysis_signal_r03",
    "truthpho_embedding_id",
    "truthpho_eta",
    "truthpho_g4_photon_valid",
    "truthpho_gen_embedding_id",
    "truthpho_generator_barcode",
    "truthpho_hepmc_association_valid",
    "truthpho_native_track_id",
    "truthpho_native_vertex_id",
    "truthpho_phi",
    "truthpho_prompt_class",
    "truthpho_pt",
    "truthpho_reporting_guard_state",
    "truthpho_resp_class",
    "truthpho_sample_source_role",
    "truthpho_source_role",
    "truthpho_truth_isolation_r03",
    "truthpho_truth_isolation_r04",
    "truthpho_truth_isolation_valid",
    "truthpho_truth_isolation_witness",
    "truthvtx_capture_state",
    "truthvtx_embedding_id",
    "truthvtx_embedding_valid",
    "truthvtx_id",
    "truthvtx_primary_id",
    "truthvtx_z",
    "truthvtx_z_valid",
    "ue_contract_state",
    "ue_emcal",
    "ue_eta_bins",
    "ue_flow_failure_state",
    "ue_flow_mode",
    "ue_ihcal",
    "ue_n_strips",
    "ue_n_towers",
    "ue_ohcal",
    "ue_psi2",
    "ue_v2",
    "vtx_valid",
    "vtx_z",
    "w_application_count",
    "w_event",
    "w_exposure",
    "w_legacy",
    "w_period",
    "w_sample_valid",
    "w_sidi",
    "w_slice",
    "w_vertex",
    "w_vertex_valid",
    "w_xsec"
  };
}
}

std::string photonjet::production::OutputMappingStatusJson(const Config& config) {
  const auto rules=outputRules();
  nlohmann::json implemented=nlohmann::json::array(),missing=nlohmann::json::array();
  bool ready=true;
  for(const auto& field:config.output.fields) {
    if(rules.count(field.name))implemented.push_back(field.name);
    else {missing.push_back(field.name);ready&=sharedNominal.count(field.name)!=0;}
  }
  ready&=!config.output.fields.empty();
  return nlohmann::json{{"schema","NativeOutputMappingV1"},
      {"schema_sha256",config.output.schemaSha256},{"implemented",implemented},
      {"missing_branches",missing},{"shared_nominal",sharedNominal},{"ready",ready},
      {"replacement_certified",false}}.dump();
}

void PhotonJetTree::initializeOutput() {
  m_directOutput=std::make_shared<photonjet::DirectEventWriter>(m_config);
  writeFileMetadata();m_outputReady=true;
}

void photonjet::production::WriteProgress(const Config& c,const std::string& phase,
    std::int64_t processed,std::int64_t skipped,const std::string& error,bool force) {
  struct State {
    std::chrono::steady_clock::time_point started=std::chrono::steady_clock::now(),written=started;
    std::int64_t processed=0,skipped=0;std::string phase;
  };
  static std::map<std::string,State> states;
  auto& state=states[c.outputFile];
  if(c.source.maxEvents<0)throw std::runtime_error("progress requires an exact finite event budget");
  const auto skipTotal=c.source.maxEvents==0?0:c.source.firstEntry;
  if(processed<0)processed=state.processed;if(skipped<0)skipped=state.skipped;
  if(processed<state.processed||processed>c.source.maxEvents||skipped<state.skipped||skipped>skipTotal)
    throw std::runtime_error("invalid or regressing native progress count");
  state.processed=processed;state.skipped=skipped;
  const auto now=std::chrono::steady_clock::now();
  if(!force&&state.phase==phase&&processed%128!=0&&
      std::chrono::duration<double>(now-state.written).count()<5.)return;
  const double seconds=std::chrono::duration<double>(now-state.started).count();
  const double rate=seconds>0?(processed+skipped)/seconds:0.;
  const auto remaining=c.source.maxEvents+skipTotal-processed-skipped;
  const auto stamp=std::time(nullptr);std::ostringstream utc;
  utc<<std::put_time(std::gmtime(&stamp),"%Y-%m-%dT%H:%M:%SZ");
  json record={{"schema","NativeTreeProductionProgressV1"},{"phase",phase},
    {"dataset",c.source.dataset},{"sample",c.source.sample},{"run",c.source.run},
    {"segment",c.source.segment},{"source_file_ordinal",c.source.sourceFileOrdinal},
    {"source_manifest_sha256",c.source.sourceManifestSha256},{"output",c.outputFile},
    {"processed",processed+skipped},{"total",c.source.maxEvents+skipTotal},
    {"remaining",remaining},{"row_processed",processed},{"row_total",c.source.maxEvents},
    {"row_remaining",c.source.maxEvents-processed},{"deterministic_skip_processed",skipped},
    {"deterministic_skip_total",skipTotal},{"deterministic_skip_remaining",skipTotal-skipped},
    {"rate",rate},{"eta_seconds",rate>0?json(remaining/rate):json(nullptr)},
    {"updated_at",utc.str()},{"latest_error",error}};
  const auto target=c.outputFile+".progress.json",temporary=target+".tmp";
  {std::ofstream out(temporary,std::ios::trunc);out<<record.dump()<<'\n';out.flush();
    if(!out)throw std::runtime_error("cannot write native progress: "+temporary);}
  if(std::rename(temporary.c_str(),target.c_str()))throw std::runtime_error("cannot publish native progress: "+target);
  state.written=now;state.phase=phase;
}

void PhotonJetTree::writeFileMetadata()
{
  std::ostringstream text;

  // --------------------------------------------------------------------------
  // Product contract and production mode
  // --------------------------------------------------------------------------

  appendLine(text, "contract", std::string("PhotonJetTrees"));
  appendLine(text, "schema_version", 1);
  appendLine(text, "collision_system", static_cast<int>(m_config.system));
  appendLine(text, "data_kind", static_cast<int>(m_config.dataKind));
  appendLine(text, "simulation_role", static_cast<int>(m_config.simulationRole));
  appendLine(text, "input_mode", static_cast<int>(m_config.inputMode));

  // --------------------------------------------------------------------------
  // Physical source binding
  // --------------------------------------------------------------------------

  appendLine(text, "lane", m_config.source.lane);
  appendLine(text, "dataset", m_config.source.dataset);
  appendLine(text, "sample", m_config.source.sample);
  appendLine(text, "period", m_config.source.period);
  appendLine(text, "si_di_role", m_config.source.siDiRole);
  appendLine(text, "run", m_config.source.run);
  appendLine(text, "segment", m_config.source.segment);
  appendLine(text, "source_file_ordinal", m_config.source.sourceFileOrdinal);
  appendLine(text, "first_entry", m_config.source.firstEntry);
  appendLine(text, "max_events", m_config.source.maxEvents);
  appendLine(text, "input_uri_sha256", m_config.source.inputUriSha256);
  appendLine(text, "input_file_sha256", m_config.source.inputFileSha256);
  appendLine(text, "source_manifest_sha256", m_config.source.sourceManifestSha256);

  // --------------------------------------------------------------------------
  // Software and configuration provenance
  // --------------------------------------------------------------------------

  appendLine(text, "production_tag", m_config.provenance.productionTag);
  appendLine(text, "software_release", m_config.provenance.softwareRelease);
  appendLine(text, "producer_git_commit", m_config.provenance.producerGitCommit);
  appendLine(text, "producer_source_sha256", m_config.provenance.producerSourceSha256);
  appendLine(text, "producer_library_sha256", m_config.provenance.producerLibrarySha256);
  appendLine(text, "macro_sha256", m_config.provenance.macroSha256);
  appendLine(text, "configuration_sha256", m_config.provenance.configurationSha256);
  appendLine(text, "calibration_manifest_sha256", m_config.provenance.calibrationManifestSha256);

  // --------------------------------------------------------------------------
  // Calibration and producer-owned weighting bindings
  // --------------------------------------------------------------------------

  appendLine(text, "jet_energy_scale_payload", m_config.jets.energyScalePayloadPath);
  appendLine(text, "jet_energy_scale_payload_sha256", m_config.jets.energyScalePayloadSha256);
  appendLine(text, "jet_energy_scale_uses_em_fraction", m_config.jets.energyScaleUsesEmFraction ? 1 : 0);
  appendLine(text, "centrality_calibration_tag", m_config.centrality.calibrationTag);
  appendLine(text, "vertex_reweight_applied", m_config.weights.applyVertexReweight ? 1 : 0);
  appendLine(text, "vertex_reweight_file_sha256", m_config.weights.vertexReweightFileSha256);
  appendLine(text, "vertex_reweight_histogram", m_config.weights.vertexReweightHistogram);

  // --------------------------------------------------------------------------
  // Storage and relation contract
  // --------------------------------------------------------------------------

  appendLine(text, "photon_min_et_gev", m_config.photon.minEtGeV);
  appendLine(text, "photon_max_et_gev", m_config.photon.maxEtGeV);
  appendLine(text, "photon_max_abs_eta", m_config.photon.maxAbsEta);
  appendLine(text, "object_vertex_abs_z_max_cm", m_config.photon.objectVertexAbsZMaxCm);
  appendLine(text, "jet_truth_match_max_delta_r", m_config.jets.truthMatchMaxDeltaR);
  appendLine(text, "recoil_delta_phi_min", photonjet::kRecoilDeltaPhiMin);
  appendLine(text, "truth_isolation_signal_radius", m_config.photon.truthIsolationSignalRadius);
  appendLine(text, "truth_isolation_max_et_gev", m_config.photon.truthIsolationMaxEtGeV);
  appendLine(text, "calorimeter_require_good_towers", m_config.calorimeter.requireGoodTowers ? 1 : 0);
  appendLine(text, "retain_every_producer_encounter", m_config.retainEveryProducerEncounter ? 1 : 0);
  appendLine(text, "write_photon_cells", m_config.output.writePhotonCells ? 1 : 0);
  appendLine(text, "write_isolation_constituents", m_config.output.writeIsolationConstituents ? 1 : 0);
  appendLine(text, "write_photon_jet_pairs", m_config.output.writePhotonJetPairs ? 1 : 0);

  // Serialize the ordered set of shower definitions used for every candidate.
  {
    std::string definitions;
    for (const std::string& name : m_config.showerDefinitions)
    {
      if (!definitions.empty()) definitions += ',';
      definitions += name;
    }
    appendLine(text, "shower_definitions", definitions);
  }

  /*
   * Bind every reconstructed jet collection by view/radius to the exact raw,
   * calibrated, truth, input, and subtraction identities that define it.
   */
  for (const JetNodeConfig& collection : m_config.jets.nodes)
  {
    std::ostringstream key;
    key << "jet_collection_view" << static_cast<int>(collection.view)
        << "_r" << static_cast<int>(photonjet::radiusCode(collection.radius));
    appendLine(text, key.str(),
               collection.rawNode + "|" + collection.correctedNode + "|" + collection.truthNode + "|" +
                   collection.inputIdentity + "|" + collection.subtractionIdentity);
  }

  m_directOutput->writer.text("metadata",text.str());
  m_directOutput->writer.text("configuration",m_config.provenance.configurationText);
}

void PhotonJetTree::fillEventOutput() {
  auto& d=*m_directOutput;d.begin();const auto& r=m_event;
  if(r.sourceEntry!=d.nextEntry++)throw std::runtime_error("source entry gap or duplicate");
  d.put("evt_record",0);
  d.put("evt_run", r.run);
  d.put("evt_seq", r.producerEventOrdinal);
  d.put("evt_phys_seq", r.physicalEventSequence);
  d.put("evt_phys_seq_valid", r.physicalEventSequenceValid);
  d.put("evt_src_file_ordinal", r.sourceFileOrdinal);
  d.put("evt_src_entry", r.sourceEntry);
  d.put("evt_src_global_entry", r.sourceGlobalEntry);
  d.put("evt_collab_iface_version",2);
  d.put("evt_n_candidates",m_photons.size());
  d.put("evt_n_recoils",std::count_if(m_photonJetPairs.begin(),m_photonJetPairs.end(),
      [](const auto& pair){return pair.recoilWitness;}));
  // The historical node check failed exactly when the reconstructed vertex
  // the analysis requires was unavailable.
  d.put("evt_nodes_ready",r.recoVertexValid?1:0);
  // The producer applies no event gate: every encounter is kept with its GL1
  // words, run trigger table, vertex, centrality and minimum-bias decisions,
  // and readers select triggers by exact name and events downstream. 0 is
  // "no producer terminal gate", never a legacy trigger/vertex decision.
  d.put("evt_terminal_status",0);
  d.put("trig_live", r.triggerLiveBits);
  d.put("trig_scaled", r.triggerScaledBits);
  d.put("trig_packet_status", r.triggerPacketStatus);
  d.put("trig_packet_version", r.triggerPacketVersion);
  d.put("trig_decision_available", r.triggerDecisionsAvailable);
  d.put("trig_capture_version",1);
  d.put("trig_decision_state",r.triggerDecisionState);
  d.put("evt_embedded_mb_decision",static_cast<std::int32_t>(r.embeddedMinimumBias));
  d.put("vtx_z", r.recoVertexZ);
  d.put("vtx_valid", r.recoVertexValid);
  d.put("mbd_pmt_available", static_cast<std::int32_t>(r.mbdPmtAvailable));
  d.all("mbd_pmt_q", r.mbdPmtCharge);
  // PMT arms and validity are not stored per channel. The file records the
  // one rule that derives them; a channel that follows none is refused.
  if(!r.mbdPmtId.empty()) {
    if(r.mbdPmtArm.size()!=r.mbdPmtId.size()||r.mbdPmtValid.size()!=r.mbdPmtId.size())
      throw std::runtime_error("MBD PMT witnesses are not aligned");
    bool divide=true,unknown=true,valid=true;
    for(std::size_t i=0;i<r.mbdPmtId.size();++i) {
      divide&=r.mbdPmtArm[i]==r.mbdPmtId[i]/64;unknown&=r.mbdPmtArm[i]==-1;valid&=r.mbdPmtValid[i]==1;
    }
    if(!divide&&!unknown)throw std::runtime_error("MBD PMT arms follow no recordable rule");
    const std::string arm=divide?"id_div_64":"minus_one";
    if(d.policies.contains("mbd_pmt_arm")&&d.policies["mbd_pmt_arm"]!=arm)
      throw std::runtime_error("MBD PMT arm rule changed within the source");
    d.policies["mbd_pmt_arm"]=arm;
    d.policies["mbd_pmt_charge_valid"]=valid&&d.policies.value("mbd_pmt_charge_valid",1)==1?1:0;
  }
  d.put("cent_replay_version", r.centralityReplayVersion);
  d.put("cent_pmt_available", static_cast<std::int32_t>(r.centralityPmtAvailable));
  d.all("cent_pmt_time", r.centralityPmtTimeNs);
  d.all("cent_pmt_valid", r.centralityPmtValid);
  d.all("cent_pmt_selected", r.centralityPmtSelected);
  d.put("cent_inputs_valid", r.centralityInputsValid);
  d.put("cent_mbd_z", r.centralityMbdZ);
  d.put("cent_mbd_z_valid", r.centralityMbdZValid);
  d.put("cent_mbd_event", r.centralityMbdEvent);
  d.put("cent_mbd_clock", r.centralityMbdClock);
  d.put("cent_mbd_femclock", r.centralityMbdFemclock);
  d.put("cent_native_bin", r.centralityNativeBin);
  d.put("cent_native_centile", r.centralityNativeCentile);
  d.put("cent_native_valid", r.centralityNativeValid);
  d.put("cent_selected_charge", r.centralitySelectedCharge);
  d.put("cent_mb_decision", static_cast<std::int32_t>(r.minimumBiasDecision));
  d.put("mbd_native_available", r.mbdAvailable);
  d.put("mbd_native_valid", r.mbdNativeValid);
  d.put("mbd_npmt_south", r.mbdSouthNpmt);
  d.put("mbd_npmt_north", r.mbdNorthNpmt);
  d.put("mbd_native_north_npmt", r.mbdNorthNpmt);
  d.put("mbd_t0", r.mbdT0Ns);
  d.put("mbd_native_south_time_ns", r.mbdSouthTimeNs);
  d.put("mbd_native_north_time_ns", r.mbdNorthTimeNs);
  d.put("cent_value", r.analysisCentrality);
  d.put("calo_emcal_native", r.cemcEnergy);
  d.put("calo_ihcal_native", r.ihcalEnergy);
  d.put("calo_ohcal_native", r.ohcalEnergy);
  d.put("calo_total_native", r.totalCaloEnergy);
  d.put("calo_require_isgood_native", r.caloRequiredGoodTowers);
  d.put("calo_max_cluster_e",r.maximumRawClusterEnergy);
  d.put("calo_native_mask",r.caloAvailableMask);
  d.put("w_event", r.eventWeight);
  d.put("w_legacy",r.eventWeight);
  // Saved collaborator interface convention: SIM sample/vertex normalization
  // remains unverified independently of a finite composed producer weight.
  d.put("w_sample_valid",m_config.dataKind==DataKind::Data?1:-1);
  d.put("w_vertex_valid",m_config.dataKind==DataKind::Data?1:-1);
  d.put("truth_vtx_z", r.truthHardVertexZ);
  d.put("truth_mbvtx_z", r.truthMinimumBiasVertexZ);
  d.put("truth_denominator_complete", r.truthDenominatorComplete);
  d.put("truthcen_native_embedded_primary_pid22_count", r.truthPhotonEmbeddedPrimaryCount);
  d.put("truthcen_serialized_raw_count", r.truthPhotonWrittenCount);
  d.put("truthcen_rejected_count", r.truthPhotonRejectedCount);
  d.put("truthcen_rejected_invalid_kinematics_count", r.truthPhotonRejectedKinematicsCount);
  d.put("truthcen_rejected_nonfinite_isolation_count", r.truthPhotonRejectedIsolationCount);
  d.put("truthcen_duplicate_track_count", r.truthPhotonDuplicateTrackCount);
  d.put("truthcen_isolation_input_incomplete_count", r.truthPhotonIsolationIncompleteCount);
  d.put("truthcen_analysis_signal_count", r.truthPhotonAnalysisSignalCount);
  d.put("recopho_capture_version", r.recoPhotonCaptureVersion);
  d.put("recopho_capture_domain", r.recoPhotonCaptureDomain);
  d.put("recopho_inspected", r.recoPhotonInspectedCount);
  d.put("recopho_eligible", r.recoPhotonEligibleCount);
  d.put("recopho_written", r.recoPhotonWrittenCount);
  d.put("recopho_unclassifiable", r.recoPhotonUnclassifiableCount);
  d.put("truthvtx_primary_id", r.truthPrimaryVertexId);
  d.put("truthvtx_capture_state",m_config.dataKind==DataKind::Data?-1:
      r.truthVertexCaptureState==CaptureState::Complete?1:0);
  d.put("ue_flow_mode", r.jetBackgroundFlowMode);
  d.put("ue_v2", r.jetBackgroundV2);
  d.put("ue_psi2", r.jetBackgroundPsi2);
  d.all("ue_emcal", r.jetBackgroundUEEmcal);
  d.all("ue_ihcal", r.jetBackgroundUEIhcal);
  d.all("ue_ohcal", r.jetBackgroundUEOhcal);
  d.put("ue_n_strips", r.jetBackgroundNStrips);
  d.put("ue_n_towers", r.jetBackgroundNTowers);
  d.put("ue_flow_failure_state", static_cast<std::int32_t>(r.jetBackgroundFlowFailure));
  d.put("truthcen_capture_state",m_config.dataKind==DataKind::Data?-1:
      r.truthPhotonCaptureState==CaptureState::Complete?1:0);
  d.put("truthjetcap_state",m_config.dataKind==DataKind::Data?-1:
      r.truthJetCaptureState==CaptureState::NotApplicable?0:1);
  d.all("truthjetcap_container_valid",r.truthJetContainerValid);
  if(r.recoObjectVertexInDomain)for(bool unsub:{false,true})for(int rc:{2,3,4})
    for(std::size_t i=0;i<r.recoJetView.size();++i)
      if((r.recoJetView[i]==static_cast<int>(JetView::AuAuNoSub))==unsub&&r.recoJetRadiusCode[i]==rc*10)
        d.put("jetview_available",r.recoJetCaptureState[i]==static_cast<int>(CaptureState::Complete)?1:0);
  d.put("ue_contract_state",static_cast<std::int32_t>(r.jetBackgroundState));
  if(r.jetBackgroundUEEmcal.size()!=r.jetBackgroundUEIhcal.size()||
      r.jetBackgroundUEEmcal.size()!=r.jetBackgroundUEOhcal.size())
    throw std::runtime_error("background calorimeter eta-bin counts disagree");
  d.put("ue_eta_bins",r.jetBackgroundUEEmcal.size());
  // The snapshot is a cumulative hardware counter, never an event weight.
  d.put("gl1_bco",r.triggerBunchNumber);
  const TriggerScalerSnapshotRecord* scaler=nullptr;
  if(r.triggerScalerSnapshotId>=0) {
    if(static_cast<std::size_t>(r.triggerScalerSnapshotId)>=m_triggerScalerSnapshots.size())
      throw std::runtime_error("dangling trigger scaler snapshot");
    scaler=&m_triggerScalerSnapshots[r.triggerScalerSnapshotId];
  }
  d.put("gl1_state",m_config.dataKind==DataKind::Simulation?1:
      scaler?(r.triggerPacketStatus?4:0):r.triggerDecisionState);
  d.put("gl1_discontinuity",r.triggerScalerDiscontinuityFlags);
  const std::array<std::uint64_t,64> noScaler{};
  d.all("gl1_raw",scaler?scaler->raw:noScaler);
  d.all("gl1_live",scaler?scaler->live:noScaler);
  d.all("gl1_scaled",scaler?scaler->scaled:noScaler);
  const WeightComponentRecord* weights=nullptr;
  for(const auto& value:m_weightComponents)if(value.componentType=="final") {
    if(weights)throw std::runtime_error("duplicate final event weight");
    weights=&value;
  }
  if(!weights)throw std::runtime_error("missing final event weight");
  d.put("w_slice",weights->sliceWeight);d.put("w_xsec",weights->crossSectionWeight);
  d.put("w_vertex",weights->vertexWeight);d.put("w_sidi",weights->siDiWeight);
  d.put("w_period",weights->periodWeight);d.put("w_exposure",weights->exposureWeight);
  d.put("w_application_count",weights->applicationCount);
  d.put("ppg12_weight_active_state",m_config.weights.ppg12Nominal && weights->valid ? 1 : 0);
  d.put("ppg12_weight_lane_code",m_config.weights.ppg12Nominal ? m_config.weights.ppg12LaneCode : 0);
  d.put("recopho_capture_state",r.recoPhotonCaptureVersion==0?0:r.recoPhotonCaptureState==CaptureState::Complete?1:2);
  d.put("mbd_native_capture_version",1);
  d.put("mbd_native_t0_finite",std::isfinite(r.mbdT0Ns));
  d.put("mbd_native_south_time_finite",std::isfinite(r.mbdSouthTimeNs));
  d.put("mbd_native_north_time_finite",std::isfinite(r.mbdNorthTimeNs));
  // Direct observations use the same five-valued augmentation states. Their
  // provenance belongs to this execution; there is no historical donor join.
  d.put("evt_aug_donor_evidence",-1);
  for(const auto& value:std::vector<std::tuple<std::string,double,bool>>{
      {"mbd_t0",r.mbdT0Ns,r.mbdAvailable},{"mbd_time_south",r.mbdSouthTimeNs,r.mbdAvailable},
      {"mbd_time_north",r.mbdNorthTimeNs,r.mbdAvailable}}) {
    const auto& name=std::get<0>(value);const double measurement=std::get<1>(value);
    const bool present=std::get<2>(value),valid=present&&std::isfinite(measurement);
    if(name!="mbd_t0")d.put(name,measurement);
    d.put(name+"_state",valid?1:present?2:0);
    d.put(name+"_reason",d.code(name+"_reason",valid?"":present?
        "native MbdOut time is nonfinite":"native MbdOut is unavailable"));
  }
  d.put("mbd_total_charge_native",r.mbdTotalCharge);
  d.put("mbd_charge_sum",r.mbdTotalCharge);
  d.put("mbd_charge_sum_state",std::isfinite(r.mbdTotalCharge)?1:r.mbdAvailable?2:0);
  d.put("mbd_charge_sum_origin",std::isfinite(r.mbdTotalCharge)?1:0);
  d.put("calo_capture_version",1);d.put("calo_available_mask",r.caloAvailableMask);
  d.put("calo_valid_mask",(r.caloValidMask>>8)&7);d.put("calo_require_isgood",r.caloRequiredGoodTowers);
  const std::array<double,4> energy={r.cemcEnergy,r.ihcalEnergy,r.ohcalEnergy,r.totalCaloEnergy};
  const std::array<std::string,4> names={"calo_emcal","calo_ihcal","calo_ohcal","calo_total"};
  for(std::size_t i=0;i<energy.size();++i) {
    const unsigned mask=i==3?7U:1U<<i;
    const bool present=(r.caloAvailableMask&mask)==mask,valid=((r.caloValidMask>>8)&mask)==mask;
    const double value=present?static_cast<float>(energy[i]):missingValue;
    d.put(names[i],value);d.put(names[i]+"_state",!present?0:valid&&std::isfinite(value)?1:2);
    d.put(names[i]+"_reason",d.code(names[i]+"_reason",!present?"calibrated calorimeter input unavailable":
        valid&&std::isfinite(value)?"":"calibrated calorimeter input incomplete or nonfinite"));
  }
  if(r.centralityPmtCharge!=r.mbdPmtCharge && !r.centralityPmtCharge.empty()) {
    // NaNs compare by their stored bits; different populations need different
    // storage and cannot be silently dropped as redundant.
    if(r.centralityPmtCharge.size()!=r.mbdPmtCharge.size() ||
       std::memcmp(r.centralityPmtCharge.data(),r.mbdPmtCharge.data(),r.mbdPmtCharge.size()*sizeof(double)))
      throw std::runtime_error("centrality and collaborator PMT charge witnesses differ");
  }
  if(r.centralityPmtId!=r.mbdPmtId && !r.centralityPmtId.empty())
    throw std::runtime_error("centrality and collaborator PMT identities differ");
  // Historical tag count: candidates whose producer score exceeds its WP80
  // witness (recorded only inside the validated domain).
  {
    const auto& production=m_config.photon.productionScore;
    const auto& edges=production.workingPointEtEdges;
    int tags=0;
    for(const auto& photon:m_photons) {
      if(production.scoreName.empty()||!std::isfinite(photon.productionScore))continue;
      const double centrality=m_event.analysisCentrality;
      if(!(photon.et>=production.validatedEtMin&&photon.et<production.validatedEtMax&&
           (std::isnan(production.validatedCentralityMin)||
            (centrality>=production.validatedCentralityMin&&centrality<production.validatedCentralityMax))))continue;
      bool finiteInputs=true;
      for(const auto& name:m_config.output.modelFeatures)finiteInputs&=std::isfinite(feature(name,photon,m_event,nullptr));
      if(!finiteInputs)continue;
      for(std::size_t bin=0;bin+1<edges.size();++bin)
        if(photon.et>=edges[bin]&&photon.et<edges[bin+1]&&photon.productionScore>production.workingPoints[1][bin])++tags;
    }
    d.put("evt_n_tags",tags);
  }
  const auto pkeys=indices(m_photons,&PhotonRecord::photonId);
  const auto tkeys=indices(m_truthPhotons,&TruthPhotonRecord::truthPhotonId);

  // Photon records and ordered shower-definition witnesses.
  for(const auto& r:m_photons) {
    // Preserve exact optional NPB invocation witnesses. Version zero means no
    // invocation was observed; neither a score nor applicability is inferred.
    for(const auto& item:std::vector<std::pair<std::string,const photonjet::ModelInvocationRecord*>>{
        {"pho_npbcap_",&r.npbCapture},{"pho_aunpbcap_",&r.auauNpbCapture}}) {
      const auto& n=item.first;const auto& v=*item.second;
      d.put(n+"capture_version",v.captureVersion);d.put(n+"evaluation_state",v.evaluationState);
      d.put(n+"domain_mode",v.domainMode);d.put(n+"in_producer_domain",v.inProducerDomain);
      d.put(n+"model_available",v.modelAvailable);d.put(n+"input_order_is_npb25",v.inputOrderIsNpb25);
      d.put(n+"inputs_captured",v.inputsCaptured);d.put(n+"inputs_finite",v.inputsFinite);
      d.put(n+"min_et",v.minEt);d.put(n+"max_et",v.maxEt);d.put(n+"max_abs_eta",v.maxAbsEta);
      d.put(n+"raw_score",v.rawScore);d.put(n+"source_photon_id",v.sourcePhotonId);
      d.put(n+"source_cluster_et",v.sourceClusterEt);d.put(n+"source_eta",v.sourceEta);d.put(n+"source_phi",v.sourcePhi);
      d.put(n+"binding_policy",v.bindingPolicy);d.put(n+"source_identity_available",v.sourceIdentityAvailable);
      d.put(n+"source_native_cluster_key",v.sourceNativeClusterKey);d.put(n+"source_encounter_ordinal",v.sourceEncounterOrdinal);
      d.put(n+"kinematic_input_bits_available",v.kinematicInputBitsAvailable);
      d.put(n+"input_et_bits",v.inputEtBits);d.put(n+"input_eta_bits",v.inputEtaBits);
      // Frozen schema stores 25 slots per candidate. Absent/partial input
      // capture retains missing slots and zero finite flags, never fake data.
      for(std::size_t i=0;i<25;++i) {
        d.put(n+"inputs",i<v.orderedInputs.size()?v.orderedInputs[i]:std::numeric_limits<float>::quiet_NaN());
        d.put(n+"input_finite",i<v.orderedInputFinite.size()?v.orderedInputFinite[i]:0);
      }
      const std::string source=item.first=="pho_npbcap_"?"photons.npb_capture_":"photons.auau_npb_capture_";
      const json constants={{source+"source_cluster_node",v.sourceClusterNode},
        {source+"model_file",v.modelFile},{source+"feature_names",v.featureNames}};
      for(auto constant=constants.begin();constant!=constants.end();++constant) {
        if(d.meta.contains(constant.key())&&d.meta[constant.key()]!=constant.value())
          throw std::runtime_error("model invocation constant changed within source: "+constant.key());
        d.meta[constant.key()]=constant.value();
      }
    }
    d.put("pho_npb_score_native",r.npbCapture.rawScore);
    d.put("pho_auau_npb_score_native",r.auauNpbCapture.rawScore);
    // The collaborator NPB fields follow the NPB score the preselection uses:
    // the PPG12 split model's npb_score in both systems, as the historical
    // augmentation defined it. Each capture record keeps its own model's cut.
    const bool auauScore=m_config.photon.preselection.npbScoreName=="auau_npb_score";
    const auto& npb=auauScore?r.auauNpbCapture:r.npbCapture;
    const auto configuredCut=[&](const char* name) {
      double value=missingValue;
      for(const auto& model:m_config.photon.models)if(model.scoreName==name)value=model.configuredCut;
      return value;
    };
    const double cut=configuredCut(auauScore?"auau_npb_score":"npb_score");
    d.put("pho_aug_donor_evidence",-1);
    d.put("pho_npbcap_association",npb.sourceIdentityAvailable?npb.bindingPolicy:0);
    d.put("pho_npbcap_configured_cut",configuredCut("npb_score"));
    d.put("pho_aunpbcap_configured_cut",configuredCut("auau_npb_score"));
    d.put("pho_npb_cut",cut);d.put("pho_npb_applicability",0);
    d.put("pho_npb_evaluation",npb.evaluationState==4?1:npb.evaluationState==5?5:
        npb.evaluationState==3?3:npb.evaluationState?2:0);
    d.put("pho_npb_score",npb.rawScore);
    const int scoreState=npb.evaluationState==4?1:npb.evaluationState==5?2:npb.captureVersion?3:0;
    d.put("pho_npb_score_state",scoreState);
    d.put("pho_npb_score_reason",d.code("pho_npb_score_reason",scoreState==1?"":
        scoreState==2?"native Compute returned an invalid score":"no valid evaluated NPB score"));
    const int passState=scoreState==1?(std::isfinite(cut)?1:0):scoreState==2?2:3;
    d.put("pho_npb_pass",passState==1?int(npb.rawScore>cut):-1);d.put("pho_npb_pass_state",passState);
    d.put("pho_npb_pass_reason",d.code("pho_npb_pass_reason",passState==1?"":
        scoreState==1?"configured NPB cut was not captured/bound":"no valid evaluated NPB score"));
    d.put("pho_npb_scoring_policy",npb.bindingPolicy);
    d.put("pho_npb_scoring_node",d.code("pho_npb_scoring_node",npb.sourceClusterNode));
    d.put("pho_npb_scoring_key",npb.sourceNativeClusterKey);
    d.put("pho_npb_scoring_encounter",npb.sourceEncounterOrdinal);
    d.put("pho_npb_scoring_pt",npb.sourceClusterEt);d.put("pho_npb_scoring_eta",npb.sourceEta);
    d.put("pho_npb_scoring_phi",npb.sourcePhi);d.put("pho_npb_scoring_et_bits",npb.inputEtBits);
    d.put("pho_npb_scoring_eta_bits",npb.inputEtaBits);
    d.put("pho_npb_scoring_evidence",npb.captureVersion?
        d.code("pho_npb_scoring_evidence",m_config.provenance.producerSourceSha256):-1);
    for(std::size_t i=0;i<photonjet::npbFeatureNames().size();++i) {
      std::string name=photonjet::npbFeatureNames()[i];
      if(name=="cluster_Et")name="cluster_et";if(name=="cluster_Eta")name="cluster_eta";
      const std::string target="pho_npb_"+name;
      const bool captured=npb.inputOrderIsNpb25&&npb.inputsCaptured&&i<npb.orderedInputs.size();
      const double value=captured?npb.orderedInputs[i]:missingValue;
      d.put(target,value);d.put(target+"_state",captured?(std::isfinite(value)?1:2):0);
      d.put(target+"_reason_index",d.code(target+"_reason_index",!captured?"no exact ordered NPB input witness":
          std::isfinite(value)?"":"native NPB input is nonfinite"));
      d.put(target+"_origin_index",captured?d.code(target+"_origin_index","native evaluator input vector"):-1);
      d.put(target+"_evidence_index",captured?d.code(target+"_evidence_index",m_config.provenance.producerSourceSha256):-1);
    }
    // Reconstructed isolation witnesses, one per retained radius in order.
    if(r.isolationWitnesses.size()!=m_config.photon.isolationRadii.size())
      throw std::runtime_error("isolation witnesses do not cover the retained radii");
    for(const auto& w:r.isolationWitnesses) {
      d.put("isow_subtraction_method",w.subtractionMethod);d.put("isow_reconstructed_or_truth",w.reconstructedOrTruth);
      d.put("isow_cone_sum",w.coneSum);d.put("isow_threshold",w.threshold);
      d.put("isow_sideband_threshold",w.sidebandThreshold);d.put("isow_pass_state",w.passState);
      d.put("isow_constituent_count",w.constituentCount);
    }
    d.put("pho_time_raw",r.timing.meanTimeSamples);d.put("pho_time_ns",r.timing.timeNs);
    d.put("pho_time_ntowers",r.timing.contributingTowerCount);d.put("pho_time_den_e",r.timing.energyDenominator);
    d.put("pho_time_raw_state",r.timing.valid?1:2);
    d.put("pho_time_raw_reason",d.code("pho_time_raw_reason",r.timing.valid?"":"native cluster timing is invalid"));
    d.put("pho_time_den_e_state",std::isfinite(r.timing.energyDenominator)?1:2);
    d.put("pho_time_den_e_reason",d.code("pho_time_den_e_reason",std::isfinite(r.timing.energyDenominator)?"":"native cluster timing denominator is nonfinite"));
    d.put("pho_key", r.nativeClusterKey);
    d.put("pho_encounter", r.encounterOrdinal);
    d.put("pho_rank_keys",r.et);d.put("pho_rank_keys",-std::abs(r.eta));d.put("pho_rank_keys",r.encounterOrdinal);
    d.put("pho_et", r.et);
    d.put("pho_eta", r.eta);
    d.put("pho_phi", r.phi);
    d.put("pho_producer_vtx_z", r.producerVertexZ);
    d.put("pho_ntime_den", r.timing.energyDenominator);
    d.put("pho_ntime_num", r.timing.energyWeightedNumerator);
    d.put("pho_ntime_ntowers", r.timing.contributingTowerCount);
    d.put("pho_ntime_finite", r.timing.finite);
    d.put("pho_ntime_version",1);
    d.put("pho_dom_track_id", r.dominantTruth.trackId);
    d.put("pho_dom_pid", r.dominantTruth.pid);
    d.put("pho_dom_barcode", r.dominantTruth.barcode);
    d.put("pho_dom_embedding_id", r.dominantTruth.embeddingId);
    d.put("pho_dom_vertex_id", r.dominantTruth.vertexId);
    d.put("pho_dom_energy", r.dominantTruth.energyContribution);
    d.put("pho_dom_state", static_cast<std::int32_t>(r.dominantTruth.state));
    d.put("pho_dom_evaluator_mode", static_cast<std::int32_t>(r.dominantTruth.evaluator));
    d.put("pho_weta_cogx", r.primaryShapes[0]);
    d.put("pho_wphi_cogx", r.primaryShapes[1]);
    d.put("pho_weta33_cogx", r.primaryShapes[2]);
    d.put("pho_wphi33_cogx", r.primaryShapes[3]);
    d.put("pho_weta35_cogx", r.primaryShapes[4]);
    d.put("pho_wphi53_cogx", r.primaryShapes[5]);
    d.put("pho_et1", r.primaryShapes[6]);
    for(const auto& item:std::vector<std::pair<std::string,std::pair<int,int>>>{{"e11_over_e33",{10,11}},{"e32_over_e35",{12,13}}}) {
      const double a=r.primaryShapes[item.second.first],b=r.primaryShapes[item.second.second];
      d.put("pho_"+item.first,std::isfinite(a)&&std::isfinite(b)&&b>0?a/b:missingValue);
    }
    const double e11=r.primaryShapes[10],e33=r.primaryShapes[11],e32=r.primaryShapes[12],e35=r.primaryShapes[13];
    const double ratio11=std::isfinite(e11)&&std::isfinite(e33)&&e33>0?e11/e33:missingValue;
    const double ratio32=std::isfinite(e32)&&std::isfinite(e35)&&e35>0?e32/e35:missingValue;
    d.put("pho_presel_reference",ratio11<.98&&r.primaryShapes[6]>.60&&r.primaryShapes[6]<1.&&
        ratio32>.80&&ratio32<1.&&r.primaryShapes[0]<.60);
    // Active preselection: -1 when none is configured.
    {
      const auto& pre=m_config.photon.preselection;
      int active=-1;
      if(pre.variant=="newPPG12") {
        double cut=missingValue;
        for(const auto& model:m_config.photon.models)if(model.scoreName==pre.npbScoreName)cut=model.configuredCut;
        const double npbScore=(pre.npbScoreName=="auau_npb_score"?r.auauNpbCapture:r.npbCapture).rawScore;
        const double et1=r.primaryShapes[6];
        active=std::isfinite(npbScore)&&std::isfinite(et1)&&std::isfinite(ratio11)&&std::isfinite(ratio32)&&
            npbScore>cut&&et1>pre.et1Min&&et1<pre.et1Max&&ratio11<pre.e11OverE33Max&&
            ratio32>pre.e32OverE35Min&&ratio32<pre.e32OverE35Max;
      }
      d.put("pho_presel_active",active);
    }
    bool allFinite=true;
    for(const auto& name:m_config.output.modelFeatures) {
      const auto value=feature(name,r,m_event,nullptr);d.put("pho_prodbdt_inputs",value);allFinite&=std::isfinite(value);
    }
    // The producer score with its historical witnesses. Applicability:
    // 0 validated domain, 1 diagnostic extrapolation, 2 no producer score
    // configured, 3 a nonfinite model input. Working points are recorded
    // only inside the validated domain.
    {
      const auto& production=m_config.photon.productionScore;
      const double raw=r.productionScore;
      int applicability=2;
      if(!production.scoreName.empty()) {
        const double centrality=m_event.analysisCentrality;
        const bool validated=r.et>=production.validatedEtMin&&r.et<production.validatedEtMax&&
            (std::isnan(production.validatedCentralityMin)||
             (centrality>=production.validatedCentralityMin&&centrality<production.validatedCentralityMax));
        applicability=!allFinite?3:validated?0:1;
      }
      std::array<double,3> wp{missingValue,missingValue,missingValue};
      const auto& edges=production.workingPointEtEdges;
      if(applicability==0)for(std::size_t bin=0;bin+1<edges.size();++bin)
        if(r.et>=edges[bin]&&r.et<edges[bin+1])for(std::size_t k=0;k<wp.size();++k)wp[k]=production.workingPoints[k][bin];
      d.put("pho_prodbdt_score",raw);d.put("pho_prodbdt_finite",std::isfinite(raw));
      d.put("pho_prodbdt_applicability",applicability);
      const std::array<const char*,3> names{{"70","80","90"}};
      for(std::size_t k=0;k<wp.size();++k) {
        d.put(std::string("pho_prodbdt_wp")+names[k],wp[k]);
        d.put(std::string("pho_prodbdt_delta_wp")+names[k],raw-wp[k]);
      }
    }
    d.put("pho_finite_feature_state",allFinite);d.put("pho_below15_state",r.et<15.);
    d.put("pho_presel_bits",7ULL|(static_cast<std::uint64_t>(allFinite)<<3));
    for(const auto& name:m_config.showerDefinitions) {
      const auto it=std::find_if(r.showerShapes.begin(),r.showerShapes.end(),[&](const auto& v){return v.definitionName==name;});
      if(it==r.showerShapes.end())throw std::runtime_error("missing captured shower definition: "+name);
      for(const auto& featureName:m_config.output.modelFeatures)
        d.put("shv_features",feature(featureName,r,m_event,&*it));
      const auto& r=*it;
      d.put("shv_finite_feature_state",std::isfinite(r.wetaCogX)&&std::isfinite(r.wphiCogX)&&
        std::isfinite(r.weta33CogX)&&std::isfinite(r.wphi33CogX)&&std::isfinite(r.e11OverE33)&&std::isfinite(r.e32OverE35)&&
        std::isfinite(r.et1)&&std::isfinite(r.et2)&&std::isfinite(r.et3)&&std::isfinite(r.et4));
      d.put("shv_cog_eta", r.centerOfGravityEta);
      d.put("shv_cog_phi", r.centerOfGravityPhi);
      d.put("shv_raw_center_eta", r.rawCenterEta);
      d.put("shv_raw_center_phi", r.rawCenterPhi);
      d.put("shv_center_eta_index", r.centerEtaIndex);
      d.put("shv_center_phi_index", r.centerPhiIndex);
      d.put("shv_e11", r.e11);
      d.put("shv_e33", r.e33);
      d.put("shv_e32", r.e32);
      d.put("shv_e35", r.e35);
      d.put("shv_e11_over_e33", r.e11OverE33);
      d.put("shv_e32_over_e35", r.e32OverE35);
      d.put("shv_weta_cogx", r.wetaCogX);
      d.put("shv_wphi_cogx", r.wphiCogX);
      d.put("shv_weta33_cogx", r.weta33CogX);
      d.put("shv_wphi33_cogx", r.wphi33CogX);
      d.put("shv_moment_eta_numerator", r.momentSums[0]);
      d.put("shv_moment_phi_numerator", r.momentSums[1]);
      d.put("shv_moment_denominator", r.momentSums[2]);
      d.put("shv_moment33_eta_numerator", r.momentSums[3]);
      d.put("shv_moment33_phi_numerator", r.momentSums[4]);
      d.put("shv_moment33_denominator", r.momentSums[5]);
      d.put("shv_energy_source", r.energySource);
      d.put("shv_rectangular_membership", r.rectangularMembership);
      d.put("shv_moment_membership", r.momentMembership);
      d.put("shv_good_cell_count", r.goodCellCount);
      d.put("shv_owned_cell_count", r.ownedCellCount);
      d.put("shv_active_sum_cell_count", r.activeSumCellCount);
      d.put("shv_active_moment_cell_count", r.activeMomentCellCount);
      d.put("shv_et1", r.et1);
      d.put("shv_et2", r.et2);
      d.put("shv_et3", r.et3);
      d.put("shv_et4", r.et4);
      d.put("shv_exact_zero_count", r.zeroCellCount);
      d.put("shv_negative_count", r.negativeCellCount);
      d.put("shv_nonfinite_count", r.nonFiniteCellCount);
    }
  }

  // Jet references address the canonical collection order, regardless of the
  // order in which reconstruction registered the nodes.
  std::map<Key,int> jkeys,tjkeys;int ji=0,ti=0;
  for(bool unsub:{false,true})for(int rc:{2,3,4}) {
    const std::string prefix=(unsub?"jetunsub0":"jet0")+std::to_string(rc);
    for(const auto& r:m_jets)if(r.radius==rc/10. && (r.view==JetView::AuAuNoSub)==unsub) {
      if(!jkeys.emplace(key(r.jetId),ji++).second)throw std::runtime_error("duplicate jet identity");
      d.put(prefix+"_rawpt",r.rawPt);d.put(prefix+"_pt",r.correctedPt);
      d.put(prefix+"_eta",r.eta);d.put(prefix+"_phi",r.phi);d.put(prefix+"_area",r.area);
      if(!unsub) {
        d.put(prefix+"_jes_valid",r.calibrationValid);
        d.put(prefix+"_emcal_frac",r.calorimeterFractions[0]);
        d.put(prefix+"_ihcal_frac",r.calorimeterFractions[1]);
        d.put(prefix+"_ohcal_frac",r.calorimeterFractions[2]);
      }
    }
  }
  if(jkeys.size()!=m_jets.size())throw std::runtime_error("unsupported jet collection");
  for(int rc:{2,3,4})for(const auto& r:m_truthJets)if(r.radius==rc/10.) {
    if(!tjkeys.emplace(key(r.truthJetId),ti++).second)throw std::runtime_error("duplicate truth jet identity");
    const auto prefix="truthjet0"+std::to_string(rc);
    d.put(prefix+"_pt",r.pt);d.put(prefix+"_eta",r.eta);d.put(prefix+"_phi",r.phi);
    d.put(prefix+"_guard_state",r.pt<5?1:r.pt>=35?2:0);
  }
  if(tjkeys.size()!=m_truthJets.size())throw std::runtime_error("unsupported truth jet collection");
  for(const auto& r:m_photonJetPairs) {
    d.put("pair_pho",local(pkeys,r.photonId));d.put("pair_jet",local(jkeys,r.jetId));
    d.put("pair_delta_phi",r.deltaPhi);d.put("pair_xjgamma",r.xJGamma);
    d.put("pair_photon_rank",r.photonRank);d.put("pair_jet_rank",r.jetRank);
    d.put("pair_recoil_state",r.recoilWitness);
    // Retained pair rows did not run wrong-object classification; zero is the
    // explicitly defined legacy unclassified default, not a truth statement.
    d.put("pair_wrong_photon_class",0);d.put("pair_wrong_recoil_class",0);
  }
  for(const auto& r:m_truthPhotons) {
    d.put("truthpho_pt", r.pt);
    d.put("truthpho_eta", r.eta);
    d.put("truthpho_phi", r.phi);
    d.put("truthpho_native_track_id", r.trackId);
    d.put("truthpho_native_vertex_id", r.vertexId);
    d.put("truthpho_generator_barcode", r.barcode);
    d.put("truthpho_embedding_id", r.embeddingId);
    d.put("truthpho_g4_photon_valid", r.geantValid);
    d.put("truthpho_hepmc_association_valid", r.generatorAssociationValid);
    d.put("truthpho_truth_isolation_r03", r.isolationR03);
    d.put("truthpho_truth_isolation_r04", r.isolationR04);
    d.put("truthpho_truth_isolation_valid", r.isolationValid);
    d.put("truthpho_analysis_signal_r03", r.analysisSignal);
    d.put("truthpho_source_role",r.analysisSignal?1:0);
    d.put("truthpho_truth_isolation_witness",r.isolationR03);
    d.put("truthpho_reporting_guard_state",r.pt<15?1:r.pt>=35?2:0);
    d.put("truthpho_gen_embedding_id",r.embeddingId);
    d.put("truthpho_sample_source_role",m_config.simulationRole==SimulationRole::PhotonJet?1:
        m_config.simulationRole==SimulationRole::InclusiveJet?2:0);
    d.put("truthpho_prompt_class", static_cast<std::int32_t>(r.promptClass));
  }
  for(const auto& r:m_truthVertices) {
    d.put("truthvtx_id",r.vertexId);d.put("truthvtx_embedding_id",r.embeddingId);
    d.put("truthvtx_embedding_valid",r.embeddingValid);d.put("truthvtx_z",r.z);d.put("truthvtx_z_valid",r.valid);
  }
  for(const auto& r:m_isolationConstituents) {
    if(r.radius!=.4)throw std::runtime_error("constituents must be stored once at R=0.4");
    d.put("isoc_pho",local(pkeys,r.photonId));d.put("isoc_native_source",r.subsystem);
    d.put("isoc_phosub_residual",missingValue);
    d.put("isoc_delta_eta", r.deltaEta);
    d.put("isoc_delta_phi", r.deltaPhi);
    d.put("isoc_delta_r", r.deltaR);
    d.put("isoc_subsystem", r.subsystem);
    d.put("isoc_raw_energy", r.rawEnergy);
    d.put("isoc_calibrated_energy", r.transverseEnergy);
    d.put("isoc_sub1_energy", r.subtractedTransverseEnergy);
    d.put("isoc_quality_state", r.qualityState);
    d.put("isoc_mask_state", r.maskState);
    d.put("isoc_candidate_removal_state", r.candidateRemoved);
    d.put("isoc_native_key", r.nativeKey);
  }
  // Cell support follows the captured H70/H0 centers, including wrapped phi.
  for(const auto& r:m_photonCells) {
    const int parent=local(pkeys,r.photonId);const auto& photon=m_photons[parent];
    auto center=[&](const std::string& name) {
      const auto it=std::find_if(photon.showerShapes.begin(),photon.showerShapes.end(),
                               [&](const auto& v){return v.definitionName==name;});
      if(it==photon.showerShapes.end()||!std::isfinite(it->rawCenterEta)||!std::isfinite(it->rawCenterPhi))
        throw std::runtime_error("cell witness lacks a finite shower center");
      return std::pair<int,int>(static_cast<int>(std::floor(it->rawCenterEta)),
          (static_cast<int>(std::floor(it->rawCenterPhi))%256+256)%256);
    };
    auto wrap=[](int value){return (value+128+256)%256-128;};
    const auto c=center("H70");const int de=r.etaIndex-c.first,raw=r.phiIndex-c.second;
    int dp=wrap(raw);if(dp==-128&&raw>0)dp=128;
    int grid=0;for(const auto& item:std::vector<std::pair<std::string,int>>{{"H0",1},{"H70",2}}) {
      const auto c=center(item.first);
      if(std::abs(r.etaIndex-c.first)<=3&&std::abs(wrap(r.phiIndex-c.second))<=3)grid|=item.second;
    }
    const bool member=r.good&&std::isfinite(r.calibratedEnergy)&&r.calibratedEnergy>0.;
    d.put("cell_pho",parent);d.put("cell_ieta",r.etaIndex);d.put("cell_iphi",r.phiIndex);
    d.put("cell_local_ieta",de);d.put("cell_local_iphi",dp);d.put("cell_seed_state",de==0&&dp==0);
    d.put("cell_grid_bits",grid);d.put("cell_e",r.calibratedEnergy);d.put("cell_rawcluster_value",r.clusterMapEnergy);
    d.put("cell_is_good",r.good);d.put("cell_is_zero",r.zeroEnergy);d.put("cell_is_negative",r.negativeEnergy);
    d.put("cell_is_nonfinite",!std::isfinite(r.calibratedEnergy));d.put("cell_rawcluster_owned",r.ownedByCluster);
    d.put("cell_rawcluster_value_present",r.ownedByCluster);d.put("cell_time_samples",r.towerTimeSamples);
    d.put("cell_time_available",r.towerAvailable);d.put("cell_time_finite",std::isfinite(r.towerTimeSamples));
    d.put("cell_tower_status",r.towerStatus);d.put("cell_time_version",1);
    d.put("cell_denominator_member",member);d.put("cell_floor0_member",member);
    d.put("cell_floor70_member",r.good&&std::isfinite(r.calibratedEnergy)&&r.calibratedEnergy>.070);
  }
  // Keep raw inventory associations distinct from selected native response
  // links: a photon excluded from the native signal can still be a raw match.
  const bool sim=m_config.dataKind==DataKind::Simulation;
  std::vector<int> rawTruth(m_photons.size(),-1),assoc(m_photons.size(),sim?0:-1);
  std::vector<int> response(m_photons.size(),-128),truthResponse(m_truthPhotons.size(),-128);
  std::vector<double> distance(m_photons.size(),missingValue);
  std::set<int> seenReco,matchedRaw;bool unknown=false;
  for(const auto& r:m_photonTruthLinks) {
    const int pi=local(pkeys,r.photonId,true),ti=local(tkeys,r.truthPhotonId,true);
    const int raw=local(tkeys,r.associationTruthPhotonId,true),state=static_cast<int>(r.state);
    if(!sim||pi<0&&ti<0)throw std::runtime_error("invalid photon truth relation population");
    if(pi>=0) {
      if(!seenReco.insert(pi).second)throw std::runtime_error("duplicate photon truth relation");
      if((r.associationState==2)!=(raw>=0))throw std::runtime_error("raw association endpoint disagrees");
      assoc[pi]=r.associationState;rawTruth[pi]=raw;
      response[pi]=r.associationState==0?-1:r.associationState==2?0:1;
      unknown|=r.associationState==0;
      if(raw>=0) {
        matchedRaw.insert(raw);const auto& p=m_photons[pi];const auto& t=m_truthPhotons[raw];
        distance[pi]=std::hypot(p.eta-t.eta,std::remainder(p.phi-t.phi,2*std::acos(-1.)));
      }
    }
    if(state!=0&&state!=3&&state!=4&&state!=5)throw std::runtime_error("unsupported photon response state");
    d.put("plink_reco",pi);d.put("plink_truth",ti);d.put("plink_class",state==0?-1:state==3?0:state==4?1:2);
    d.put("plink_metric",r.deltaR);d.put("plink_policy",1);
  }
  if(sim&&seenReco.size()!=m_photons.size())throw std::runtime_error("incomplete photon relation population");
  const bool complete=m_event.recoPhotonCaptureVersion==1&&m_event.recoPhotonCaptureState==CaptureState::Complete&&
                      m_event.recoPhotonUnclassifiableCount==0;
  const bool outside=m_event.recoPhotonCaptureVersion==0&&!m_event.recoVertexValid;
  for(std::size_t i=0;i<truthResponse.size();++i)if(!matchedRaw.count(i))
    truthResponse[i]=outside?3:complete&&!unknown?2:-1;
  d.all("pho_assoc_truthpho",rawTruth);d.all("pho_assoc_state",assoc);d.all("pho_assoc_dr",distance);
  d.all("pho_resp_class",response);d.all("truthpho_resp_class",truthResponse);
  for(std::size_t i=0;i<m_photons.size();++i) {
    const auto& witness=m_photons[i].dominantTruth;
    const bool evaluator=sim&&witness.evaluator!=DominantTruthEvaluator::Unavailable&&
        witness.evaluator!=DominantTruthEvaluator::NotApplicable;
    const bool matched=rawTruth[i]>=0&&m_truthPhotons[rawTruth[i]].analysisSignal;
    d.put("pho_truth_match_state",evaluator?int(matched):-1);
    d.put("pho_truth_match_barcode",matched?m_truthPhotons[rawTruth[i]].barcode:-1);
  }
  for(const auto& r:m_jetTruthLinks) {
    const int ri=local(jkeys,r.jetId,true),ti=local(tjkeys,r.truthJetId,true),state=static_cast<int>(r.state);
    if(!sim||ri<0&&ti<0||state==1||state<0||state>5)throw std::runtime_error("invalid jet truth link");
    if(r.selectedMatch!=(state==3))throw std::runtime_error("jet selected-match flag differs");
    d.put("jlink_reco",ri);d.put("jlink_truth",ti);d.put("jlink_metric",r.deltaR);
    d.put("jlink_class",state==0?-1:state==2?5:state==3?0:state==4?1:2);
  }
  // Explicit schema contracts for currently unfilled optional collections.
  d.empty("jcon_constituent_ordinal");
  d.empty("jcon_energy");
  d.empty("jcon_eta");
  d.empty("jcon_native_key");
  d.empty("jcon_native_source");
  d.empty("jcon_phi");
  d.empty("jcon_quality_state");
  d.empty("jcon_subsystem");
  d.empty("snap_quota_class");
  d.empty("snap_selection_reason");
  d.empty("snap_serialized_payload_hash");
  d.finish();
  // Centrality-fit denominator precedes the analysis vertex cut and carries
  // unit event counts. Applying centrality weights here would bias their fit.
  if(d.centrality && r.recoVertexValid && std::isfinite(r.analysisCentrality) &&
      r.analysisCentrality>=0. && r.analysisCentrality<=100.)
    d.centrality->Fill(r.analysisCentrality);
  m_source.lastEntry=r.sourceEntry;
}

void PhotonJetTree::fillRejectedOutput(const UpstreamRejectedEventRecord& r) {
  auto& d=*m_directOutput;d.begin();
  if(r.sourceEntry!=d.nextEntry++)throw std::runtime_error("rejected source entry gap or duplicate");
  d.put("evt_record",1);d.put("evt_run",r.run);d.put("evt_phys_seq",r.physicalEventSequence);
  d.put("evt_phys_seq_valid",r.physicalEventSequence>=0);
  d.put("evt_src_entry",r.sourceEntry);d.put("evt_src_file_ordinal",m_source.sourceFileOrdinal);
  d.put("evt_src_global_entry",r.sourceEntry);
  d.finish(true);
}

void PhotonJetTree::writeCompletionMetadata() {
  auto& d=*m_directOutput;
  m_source.completed=!m_aborted;
  if(d.centrality)d.writer.object("hist/SIM/h_centrality",*d.centrality);
  if(!d.entries){d.begin();d.finish(false,true);}
  json counters=json::object();for(const auto& f:d.specs)if(f.length==-1&&f.type!="string")counters[f.name]=f.counter;
  const bool auau=m_config.system==CollisionSystem::AuAu;
  json collections={{"reco",json::array()},{"truth",json::array()}};
  for(bool unsub:{false,true})for(int rc:{2,3,4})if(!unsub||auau) {
    const std::string prefix=(unsub?"jetunsub0":"jet0")+std::to_string(rc);
    collections["reco"].push_back({{"name",prefix},{"radius_code",rc},{"radius",bits(rc/10.)},
      {"input_identity",auau?(unsub?"towerinfo_retower_unsubtracted":"towerinfo_sub1"):"towerinfo_calib"},
      {"algorithm","antikt"},{"subtraction_identity",auau&&!unsub?"SUB1":"none"}});
  }
  for(int rc:{2,3,4})collections["truth"].push_back({{"name","truthjet0"+std::to_string(rc)},
    {"radius_code",rc},{"radius",bits(rc/10.)}});
  // File constants shared by every row: shower views, isolation witness radii,
  // collection vocabularies, weight targets and the producer score identity.
  json constants=d.meta.is_null()?json::object():d.meta;
  {
    json names=json::array(),digests=json::array(),floors=json::array(),byFloor=json::array();
    for(const auto& name:m_config.showerDefinitions) {
      const auto& view=photonjet::showerViewDefinition(name);
      names.push_back(view.name);digests.push_back(view.semanticSha256);floors.push_back(bits(view.floorGeV));
    }
    // Candidate view list in the historical order: 70 MeV views, then 0 MeV.
    for(bool floor70:{true,false})for(const auto& name:m_config.showerDefinitions)
      if((photonjet::showerViewDefinition(name).floorGeV>0)==floor70)byFloor.push_back(name);
    constants["showerViews.definition_name"]=names;constants["showerViews.semantic_sha256"]=digests;
    constants["showerViews.floor_gev"]=floors;constants["photons.shower_definition_views"]=byFloor;
    json radii=json::array();for(const double radius:m_config.photon.isolationRadii)radii.push_back(bits(radius));
    constants["isolationWitnesses.radius"]=radii;
    json viewCodes=json::array(),viewInputs=json::array(),viewSubtraction=json::array();
    for(const auto& entry:collections["reco"]) {
      viewCodes.push_back(entry["radius_code"]);viewInputs.push_back(entry["input_identity"]);
      viewSubtraction.push_back(entry["subtraction_identity"]);
    }
    constants["view_radius_codes"]=viewCodes;constants["view_input_identity"]=viewInputs;
    constants["view_subtraction_identity"]=viewSubtraction;
    constants["weights.component_type"]="event";constants["weights.target_type"]=1;
    if(m_config.dataKind==DataKind::Simulation) {
      constants["truth_radius_codes"]={2,3,4};
      constants["truthJets.algorithm"]="antikt";constants["truthJets.ownership_state"]="source_owned";
    }
    const auto& production=m_config.photon.productionScore;
    if(!production.scoreName.empty()) {
      constants["models.model_sha256"]=production.modelSha256;
      constants["models.shower_definition_id"]=production.showerDefinition;
      constants["models.shower_semantic_sha256"]=photonjet::showerViewDefinition(production.showerDefinition).semanticSha256;
    }
  }
  json meta={{"schema",{{"schema","RJEventTreeV1"},{"version","1.0"},{"schema_sha256",m_config.output.schemaSha256},
        {"entries",d.entries},{"retained_events",m_source.retainedEvents},{"system",auau?"auau":"pp"},
        {"lane",m_config.source.lane},{"layout","native_direct"},{"encoding","counted"},
        {"counters",counters},{"absent_branches",d.missing},{"strides",{{"shv_features",m_config.output.modelFeatures.size()},{"pho_prodbdt_inputs",m_config.output.modelFeatures.size()},
          {"pho_rank_keys",3},{"pho_npbcap_inputs",25},{"pho_npbcap_input_finite",25},
          {"pho_aunpbcap_inputs",25},{"pho_aunpbcap_input_finite",25}}}, {"filled_reserved",json::array()}}},
    {"source",{{"lane",m_config.source.lane},{"dataset",m_source.dataset},{"sample",m_source.sample},{"period",m_source.period},
        {"si_di_role",m_source.siDiRole},{"run",m_source.run},{"segment",m_source.segment},
        {"source_file_ordinal",m_source.sourceFileOrdinal},{"first_entry",m_source.firstEntry},{"last_entry",m_source.lastEntry},
        {"input_uri_sha256",m_source.inputUriSha256},{"input_file_sha256",m_source.inputFileSha256},
        {"source_manifest_sha256",m_source.sourceManifestSha256},{"completed",m_source.completed},
        {"retained_events",m_source.retainedEvents},{"upstream_rejected_events",m_source.upstreamRejectedEvents}}},
    {"collections",collections},{"codes",d.codes},{"constants",constants},{"policies",d.policies},
    {"weighting",{{"schema","ProducerWeightV1"},
        {"composition",m_config.weights.ppg12Nominal?"slice*vertex*sidi*period":
            (m_config.weights.applyVertexReweight?"vertex":"unity")},
        {"vertex_enabled",m_config.weights.applyVertexReweight},
        {"vertex_coordinate",m_config.weights.applyVertexReweight?"truth_hard":"none"},
        {"vertex_lookup",m_config.weights.applyVertexReweight?"TH1::Interpolate":"none"},
        {"double_interaction",m_config.weights.doubleInteraction},
        {"vertex_file_sha256",m_config.weights.vertexReweightFileSha256},
        {"vertex_histogram",m_config.weights.vertexReweightHistogram},
        {"absolute_normalization_applied",false},{"centrality_applied",false},
        {"aliases",m_config.weights.ppg12Nominal?
            json{{"w_xsec","w_slice"},{"w_exposure","w_period"}}:json::object()}}},
    {"native_capture",{{"status",d.missing.empty()?"MAPPING_COMPLETE":"INCOMPLETE_MAPPING"},
        {"missing_branches",d.missing},{"coverage",d.covered},{"replacement_certified",false},
        {"producer_source_sha256",m_config.provenance.producerSourceSha256},
        {"producer_library_sha256",m_config.provenance.producerLibrarySha256},
        {"shared_nominal",sharedNominal},{"ready",d.missing==sharedNominal}}}};
  // DATA run trigger configuration: one entry per GL1 bit. Readers select a
  // stream by its exact name and enabled prescale, never by bit index.
  if(!m_triggerRunInfo.empty()) {
    if(m_triggerRunInfo.size()!=64)throw std::runtime_error("trigger run information does not hold 64 bits");
    json table={{"name",json::array()},{"state",json::array()},{"initial_prescale",json::array()},
                {"prescale",json::array()},{"run_raw",json::array()},{"run_live",json::array()},{"run_scaled",json::array()}};
    for(std::size_t bit=0;bit<m_triggerRunInfo.size();++bit) {
      const auto& row=m_triggerRunInfo[bit];
      if(row.bit!=static_cast<int>(bit))throw std::runtime_error("trigger run information is out of bit order");
      table["name"].push_back(row.name);table["state"].push_back(row.state);
      table["initial_prescale"].push_back(std::isfinite(row.hardwarePrescale)?static_cast<int>(row.hardwarePrescale):-1);
      table["prescale"].push_back(bits(row.runAveragePrescale));
      table["run_raw"].push_back(row.rawCount);table["run_live"].push_back(row.liveCount);
      table["run_scaled"].push_back(row.scaledCount);
    }
    meta["trigger_run_info"]=table;
  }
  for(auto it=meta.begin();it!=meta.end();++it)d.writer.text("meta/"+it.key(),it.value().dump());
  std::ostringstream text;
  appendLine(text,"completion_status",m_aborted?"aborted":"complete");
  appendLine(text,"encountered_events",m_source.encounteredEvents);
  appendLine(text,"retained_events",m_source.retainedEvents);
  appendLine(text,"upstream_rejected_events",m_source.upstreamRejectedEvents);
  appendLine(text,"upstream_observed_events",m_upstreamObservedEvents);
  appendLine(text,"first_entry",m_source.firstEntry);appendLine(text,"last_entry",m_source.lastEntry);
  d.writer.text("completion",text.str());
  photonjet::production::WriteProgress(m_config,m_aborted?"aborted":"complete",
      static_cast<std::int64_t>(m_upstreamObservedEvents),-1,
      m_aborted?"event loop aborted":"",true);
}

void PhotonJetTree::closeOutput() {
  if(!m_directOutput)return;
  m_directOutput->writer.close();m_directOutput.reset();
}
