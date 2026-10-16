/* Native source definitions from effective producer arguments and retained
   numeric model/lookup witnesses. Missing reconstruction ancestry is explicit. */
#ifndef PHOTONJET_SOURCE_SCIENTIFIC_DEFINITIONS_H
#define PHOTONJET_SOURCE_SCIENTIFIC_DEFINITIONS_H
#include "PrimaryShowerDefinition.h"
namespace photonjet_science {
inline Json bound(double value) { return std::isfinite(value)?Json(value):std::isinf(value)?Json(value>0?"positive_infinity":"negative_infinity"):absent("unset optional bound"); }
inline Json line(const photonjet::IsolationLine& x) {return {{"intercept",bound(x.interceptGeV)},{"slope",bound(x.slope)},{"sideband_gap_gev",bound(x.sidebandGapGeV)}};}
inline Json sourceDefinitions(const photonjet::Config& c,const std::set<std::string>& primary,
                              const Json& observed=Json::object()) {
  using namespace photonjet;
  const bool data=c.dataKind==DataKind::Data,auau=c.system==CollisionSystem::AuAu;
  const auto witness=[&](const std::string& key) {return observed.value(key,unknown("native reconstruction witness absent: "+key));};
  Json result={{"schema","NativeSourceScientificDefinitionsV1"},{"branch_definitions",primaryBranchDefinitions(primary)}};
  const auto& event=c.classification;
  result["event_acceptance"]={{"enabled",event.enabled},{"retain_every_encounter",c.retainEveryProducerEncounter},
    {"scaled_triggers",event.scaledTriggers},{"named_triggers",event.triggers},{"scaled_bit",event.scaledBit},
    {"minimum_bias",event.requireMinimumBias},{"vertex_cut_enabled",event.useVertexCut},{"vertex_abs_max_cm",bound(event.vertexCutCm)},
    {"filter_data_period",event.filterDataByPeriod},{"period_runs",{event.period.runMin,event.period.runMaxExclusive}},
    {"centrality_percent",{event.centralityMinPercent,event.centralityMaxPercent}},
    {"centrality_measurement",auau?witness("centrality"):absent("pp has no centrality selector")}};
  result["candidate_geometry"]={{"et_min_gev",c.photon.minEtGeV},{"et_max_gev",c.photon.maxEtGeV},{"abs_eta_max",c.photon.maxAbsEta},
    {"data_object_vertex_abs_max_cm",data?Json(c.photon.objectVertexAbsZMaxCm):absent("SIM object retention does not apply DATA vertex cut")},
    {"kinematics",witness("photon_kinematics")}};
  result["vertex"]=witness("vertex");
  result["jets"]={{"truth_match_delta_r",c.jets.truthMatchMaxDeltaR},{"energy_scale_uses_em_fraction",c.jets.energyScaleUsesEmFraction},
    {"energy_scale_function",witness("jet_energy_scale")},{"collections",Json::array()}};
  for(const auto& jet:c.jets.nodes)result["jets"]["collections"].push_back({{"radius",jet.radius},{"view",int(jet.view)},
    {"input",jet.inputIdentity},{"subtraction",jet.subtractionIdentity}});
  result["ownership"]=c.weights.ppg12Nominal?Json{{"lane",c.weights.ppg12LaneCode},{"low_gev",bound(c.weights.ownershipLowGeV)},
    {"high_gev",bound(c.weights.ownershipHighGeV)},{"cluster_et_upper_gev",bound(c.weights.clusterEtUpperGeV)},
    {"tower_mask",witness("tower_mask")}}:absent("producer does not apply PPG12 source ownership");
  result["event_weight"]=data?Json{{"rule","unit"}}:Json{{"ppg12_nominal",c.weights.ppg12Nominal},
    {"slice_factor",c.weights.sliceFactor},{"interaction_factor",c.weights.interactionFactor},{"period_factor",c.weights.periodFactor},
    {"vertex_enabled",c.weights.applyVertexReweight},{"vertex_lookup",c.weights.applyVertexReweight?witness("vertex_weight"):absent("vertex weighting disabled")}};
  result["models"]=observed.value("models",Json::object());
  if(!result["models"].contains("npb"))result["models"]["npb"]=c.photon.preselection.npbScoreName.empty()?absent("no NPB score configured"):witness("npb_model");
  const auto& id=c.photon.nominalId;
  result["models"]["pho_bdt_score"]=id.usesProducerScore?
    result["models"].value(c.photon.productionScore.scoreName,unknown("loaded producer score model absent")):
    witness("nominal_model");
  result["nominal_id"]={{"feature_order",id.featureOrder},{"view",id.view},{"et",{bound(id.etMin),bound(id.etMaxExclusive)}},
    {"centrality",{bound(id.centralityMin),bound(id.centralityMaxExclusive)}},{"uses_producer_score",id.usesProducerScore},
    {"wp70",line(id.wp70)},{"wp80",line(id.wp80)},{"model",witness("nominal_model")},
    {"regions",{{"tight_comparison",">WP70"},{"non_tight_low",.1},{"non_tight_comparison",">low and <WP80"}}}};
  const auto& iso=c.photon.nominalIsolation;
  Json witnessLines=Json::array();for(const auto& x:c.photon.isolationWitness.lines)witnessLines.push_back(line(x));
  result["isolation"]={{"radii",c.photon.isolationRadii},{"witness_lines",witnessLines},{"fallback",line(c.photon.isolationWitness.fallback)},
    {"centrality",{bound(c.photon.isolationWitness.centralityMin),bound(c.photon.isolationWitness.centralityMax)}},
    {"nominal_radius",bound(iso.radius)},{"witness_index",iso.witnessIndex},{"subtraction_method",iso.subtractionMethod},
    {"threshold",line(iso.threshold)},{"sideband_gap_gev",bound(iso.threshold.sidebandGapGeV)},{"stored_thresholds",iso.storedThresholds},
    {"cone_calculation",witness("isolation_calculation")}};
  result["cluster_timing"]={{"sample_ns",c.photon.timingSampleNs},{"calculation",witness("cluster_timing")}};
  result["truth_classification"]=data?absent("DATA has no truth"):Json{{"signal_isolation_radius",c.photon.truthIsolationSignalRadius},
    {"isolation_max_et_gev",c.photon.truthIsolationMaxEtGeV},{"core_radius",c.photon.truthIsolationCoreRadius},
    {"classification",witness("truth_classification")}};
  result["response_matching"]=data?absent("DATA has no response matching"):Json{{"jet_delta_r",c.jets.truthMatchMaxDeltaR},{"photon_matching",witness("photon_matching")}};
  return portable(result);
}
} // namespace photonjet_science
#endif
