// Production.cc
//
// What production job are we running?
// Owns configuration, profile resolution and plan validation.
// Does not own conditions binding or reconstruction registration.
// Called through Production.h before PhotonJetTree captures the event.
//
#include "Production.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace photonjet
{
namespace production
{
namespace detail
{
[[noreturn]] void fail(const std::string& message)
{
  throw std::runtime_error("PhotonJetTree production: " + message);
}

std::string trim(const std::string& value)
{
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
  {
    return {};
  }
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string lower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

}  // namespace detail

namespace
{
using photonjet::production::Conditions;
using photonjet::production::InputLayout;
using photonjet::production::InputStream;
using photonjet::production::Job;
using photonjet::production::Plan;
using photonjet::production::Reconstruction;
using photonjet::production::TruthJetMode;
using photonjet::production::detail::fail;
using photonjet::production::detail::lower;
// Frozen production constants
//
// These are the values the reference production accepted. They are checked
// against the configuration rather than read from it, so a configuration
// edit cannot silently move a nominal choice.
//
constexpr const char* kNominalCentralitySet = "sam_fit_20260913_v1";

constexpr const char* kNominalCentralityManifestSha256 =
    "4a11e0e10ae9b82062325c3275679935da2259f4c91f796950a6b5ddd2e1501c";

constexpr const char* kNominalJetEnergyScaleSha256 =
    "76a8788fdb4e0b3859d01361f884e295ea609bb8c1cc188105389214db2de27d";

constexpr int kQualifiedCentralityFallbackRun = 68144;

std::string readTextFile(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  if (!input)
  {
    fail("cannot read file: " + path);
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

// ${NAME} expansion for configured paths. An unset variable is an error, not
// an empty string, because an empty path would silently select nothing.
std::string expandEnvironment(const std::string& input)
{
  std::string output;
  for (std::size_t i = 0; i < input.size();)
  {
    if (input[i] == '$' && i + 1 < input.size() && input[i + 1] == '{')
    {
      const std::size_t close = input.find('}', i + 2);
      if (close == std::string::npos)
      {
        fail("unterminated environment variable in path: " + input);
      }
      const std::string key = input.substr(i + 2, close - i - 2);
      const char* value = std::getenv(key.c_str());
      if (!value)
      {
        fail("required environment variable is unset: " + key);
      }
      output += value;
      i = close + 1;
      continue;
    }
    output += input[i];
    ++i;
  }
  return output;
}
//
// YAML helpers
//
template <class T>
T required(const YAML::Node& node, const std::string& key, const std::string& context)
{
  if (!node || !node[key] || node[key].IsNull())
  {
    fail("missing required configuration key '" + key + "' in " + context);
  }
  return node[key].as<T>();
}

template <class T>
T optional(const YAML::Node& node, const std::string& key, const T& fallback)
{
  if (!node || !node[key] || node[key].IsNull())
  {
    return fallback;
  }
  return node[key].as<T>();
}

std::string optionalString(const YAML::Node& node, const std::string& key,
                           const std::string& fallback = {})
{
  if (!node || !node[key] || node[key].IsNull())
  {
    return fallback;
  }
  return node[key].as<std::string>();
}

std::vector<std::string> stringList(const YAML::Node& node)
{
  std::vector<std::string> result;
  if (node && node.IsSequence())
  {
    for (const YAML::Node& value : node)
    {
      result.push_back(value.as<std::string>());
    }
  }
  return result;
}

std::vector<double> doubleList(const YAML::Node& node)
{
  std::vector<double> result;
  if (node && node.IsSequence())
  {
    for (const YAML::Node& value : node)
    {
      result.push_back(value.as<double>());
    }
  }
  return result;
}

std::vector<int> intList(const YAML::Node& node)
{
  std::vector<int> result;
  if (node && node.IsSequence())
  {
    for (const YAML::Node& value : node)
    {
      result.push_back(value.as<int>());
    }
  }
  return result;
}

TruthJetMode parseTruthJetMode(const std::string& value)
{
  const std::string mode = lower(value);
  if (mode.empty() || mode == "none") return TruthJetMode::None;
  if (mode == "auto") return TruthJetMode::Auto;
  if (mode == "dst") return TruthJetMode::DST;
  if (mode == "build") return TruthJetMode::Build;
  fail("unknown truth-jet mode: " + value);
}

photonjet::JetView parseJetView(const std::string& value)
{
  const std::string view = lower(value);
  if (view == "pp") return photonjet::JetView::PP;
  if (view == "auau_sub1" || view == "sub1") return photonjet::JetView::AuAuSub1;
  if (view == "auau_nosub" || view == "nosub") return photonjet::JetView::AuAuNoSub;
  fail("unknown jet view: " + value);
}

}  // namespace

using detail::fail;
using detail::lower;
using detail::fileSha256;
using detail::validSha256;
using detail::loadInputs;
using detail::hasInputStream;
using detail::kArchivedDoubleInteractionRun;
using detail::kDataRunThreshold;

// LoadPlan
//
Plan LoadPlan(const std::string& yamlFile, const Job& job)
{
  if (job.profile.empty()) fail("production profile is empty");
  if (job.inputList.empty()) fail("input-list path is empty");
  if (job.outputFile.empty()) fail("output file name is empty");

  Plan plan;
  plan.job = job;

  plan.configurationText = readTextFile(yamlFile);
  plan.configurationSha256 = fileSha256(yamlFile);

  const YAML::Node root = YAML::Load(plan.configurationText);
  const YAML::Node profile = root["profiles"][job.profile];
  if (!profile)
  {
    fail("profile '" + job.profile + "' is absent from " + yamlFile);
  }
  const std::string context = "profiles." + job.profile;

  Reconstruction& reco = plan.reconstruction;
  Conditions& cond = plan.conditions;
  photonjet::Config& tree = plan.tree;

  // ---- mode ---------------------------------------------------------------

  const std::string system = lower(required<std::string>(profile, "system", context));
  const std::string dataKind = lower(required<std::string>(profile, "data_kind", context));
  const std::string role = lower(required<std::string>(profile, "simulation_role", context));
  const std::string inputMode = lower(required<std::string>(profile, "input_mode", context));

  reco.isPP = system == "pp";
  reco.isAuAu = system == "auau";
  reco.isData = dataKind == "data";
  reco.isSimulation = dataKind == "simulation";
  reco.isEmbedded = reco.isAuAu && reco.isSimulation;
  reco.isArchivedDoubleInteraction = inputMode == "archived_g4_only";

  if (reco.isPP == reco.isAuAu) fail("profile must name exactly one collision system");
  if (reco.isData == reco.isSimulation) fail("profile must name exactly one data kind");

  tree.system = reco.isPP ? CollisionSystem::PP : CollisionSystem::AuAu;
  tree.dataKind = reco.isData ? DataKind::Data : DataKind::Simulation;
  tree.inputMode = reco.isArchivedDoubleInteraction ? InputMode::ArchivedG4Only
                                                    : InputMode::ReconstructedDST;
  if (role == "photon_jet") tree.simulationRole = SimulationRole::PhotonJet;
  else if (role == "inclusive_jet") tree.simulationRole = SimulationRole::InclusiveJet;
  else if (role == "minimum_bias") tree.simulationRole = SimulationRole::MinimumBias;
  else tree.simulationRole = SimulationRole::None;

  tree.outputFile = job.outputFile;

  if (inputMode != "archived_g4_only" && inputMode != "reconstructed_dst") fail("unknown input mode: " + inputMode);
  if (role != "none" && role != "photon_jet" && role != "inclusive_jet" && role != "minimum_bias") fail("unknown simulation role: " + role);

  // ---- source binding -----------------------------------------------------

  tree.source.dataset = job.dataset.empty() ? optionalString(profile, "dataset_token") : job.dataset;
  tree.source.sample = job.sample;
  tree.source.period = job.period;
  tree.source.siDiRole = job.siDiRole;
  tree.source.inputUriSha256 = job.inputUriSha256;
  tree.source.inputFileSha256 = job.inputFileSha256;
  tree.source.sourceManifestSha256 = job.sourceManifestSha256;
  tree.source.run = job.run;
  tree.source.segment = job.segment;
  tree.source.sourceFileOrdinal = job.sourceFileOrdinal;
  tree.source.firstEntry = job.firstEntry;
  tree.source.maxEvents = job.nEvents;

  // ---- nodes --------------------------------------------------------------

  const YAML::Node nodes = root["nodes"];
  tree.nodes.eventHeader = optionalString(nodes, "event_header", "EventHeader");
  tree.nodes.gl1Packet = optionalString(nodes, "gl1_packet", "GL1Packet");
  tree.nodes.globalVertexMap = optionalString(nodes, "global_vertex_map", "GlobalVertexMap");
  tree.nodes.centralityInfo = optionalString(nodes, "centrality_info", "CentralityInfo");
  tree.nodes.minimumBiasInfo = optionalString(nodes, "minimum_bias_info", "MinimumBiasInfo");
  tree.nodes.mbdOut = optionalString(nodes, "mbd_out", "MbdOut");
  tree.nodes.mbdPmts = optionalString(nodes, "mbd_pmts", "MbdPmtContainer");
  tree.nodes.cemcTowers = optionalString(nodes, "cemc_towers", "TOWERINFO_CALIB_CEMC");
  tree.nodes.ihcalTowers = optionalString(nodes, "ihcal_towers", "TOWERINFO_CALIB_HCALIN");
  tree.nodes.ohcalTowers = optionalString(nodes, "ohcal_towers", "TOWERINFO_CALIB_HCALOUT");
  tree.nodes.cemcGeometry = optionalString(nodes, "cemc_geometry", "TOWERGEOM_CEMC");
  tree.nodes.truthInfo = optionalString(nodes, "truth_info", "G4TruthInfo");
  tree.nodes.hepmcEventMap = optionalString(nodes, "hepmc_event_map", "PHHepMCGenEventMap");

  // ---- calorimeter event state ---------------------------------------------

  tree.calorimeter.requireGoodTowers =
      optional<bool>(root["calorimeter"]["event_energy"], "require_good_towers", true);

  // ---- photons ------------------------------------------------------------

  const YAML::Node photons = root["photons"];
  const YAML::Node domain = photons["capture_domain"];
  tree.photon.minEtGeV = required<double>(domain, "min_et_gev", "photons.capture_domain");
  tree.photon.maxEtGeV = required<double>(domain, "max_et_gev", "photons.capture_domain");
  tree.photon.maxAbsEta = required<double>(domain, "max_abs_eta", "photons.capture_domain");
  tree.photon.objectVertexAbsZMaxCm =
      optional<double>(domain, "object_vertex_abs_z_max_cm", tree.photon.objectVertexAbsZMaxCm);
  tree.photon.rawClusterNode = optionalString(photons, "source_cluster_node", "CLUSTERINFO_CEMC");
  tree.photon.photonNode = optionalString(photons, "output_photon_node", "PHOTONCLUSTER_CEMC");
  tree.photon.topoclusterNode = optionalString(photons, "topocluster_node", "TOPOCLUSTER_ALLCALO");

  const std::vector<double> radii = doubleList(photons["isolation"]["radii"]);
  if (!radii.empty())
  {
    tree.photon.isolationRadii = radii;
  }
  tree.photon.truthIsolationSignalRadius =
      optional<double>(photons["isolation"]["truth"], "signal_radius", 0.30);
  tree.photon.truthIsolationMaxEtGeV =
      required<double>(photons["isolation"]["truth"], "max_et_gev", "photons.isolation.truth");
  tree.photon.truthIsolationCoreRadius =
      required<double>(photons["isolation"]["truth"], "self_removal_radius", "photons.isolation.truth");
  tree.photon.timingSampleNs = required<double>(photons["timing"], "sample_ns", "photons.timing");

  const std::vector<std::string> definitions = stringList(photons["shower_shapes"]["definitions"]);
  if (!definitions.empty())
  {
    tree.showerDefinitions = definitions;
  }

  // ---- reconstruction stages -----------------------------------------------

  reco.inputLayout = reco.isArchivedDoubleInteraction ? InputLayout::ArchivedDoubleInteraction
                     : reco.isData                    ? InputLayout::PairedData
                                                      : InputLayout::ReconstructedSimulation;

  reco.truthJetMode = parseTruthJetMode(
      optional<std::string>(profile, "truth_jets_mode", reco.isSimulation ? "auto" : "none"));

  reco.reconstructEventGeometry = !reco.isEmbedded;
  reco.reconstructZdc = reco.isData;
  reco.reconstructCalorimeter = reco.isData || reco.isEmbedded || reco.isArchivedDoubleInteraction;
  reco.runCalorimeterStatusSkimmer = reco.isData;
  reco.setTowerStatus = reco.reconstructCalorimeter &&
                        (!reco.isEmbedded || optional<bool>(profile, "force_embedded_tower_status", false));
  reco.reconstructCentrality = reco.isAuAu && reco.isData &&
                               optional<bool>(profile, "reconstruct_centrality", true);
  reco.runMinimumBiasClassifier = reco.reconstructCentrality &&
                                  optional<bool>(profile, "require_minimum_bias_classifier", true);
  reco.reconstructAuAuBackground = reco.isAuAu;
  reco.reconstructJets = true;
  reco.reconstructPhotons = true;

  const YAML::Node clusterBuilder = root["calorimeter"]["cluster_builder"];
  reco.clusterThresholdGeV =
      required<double>(clusterBuilder, "threshold_energy_gev", "calorimeter.cluster_builder");
  reco.clusterProfile = expandEnvironment(
      required<std::string>(clusterBuilder["profile"], "path", "calorimeter.cluster_builder.profile"));
  reco.simulationHotTowerMapKey =
      optionalString(root["calorimeter"]["tower_status"], "simulation_hot_map_key", "CEMC_hotTowers_status");

  const YAML::Node systemConfig = reco.isPP ? root["pp"] : root["auau"];

  reco.photonMinEtGeV = tree.photon.minEtGeV;
  reco.photonUseVertexCut = optional<bool>(photons["builder"], "use_vertex_cut", false);
  reco.photonVertexCutCm = optional<double>(photons["builder"], "vertex_cut_cm", 60.0);
  reco.showerShapeTowerMinEnergyGeV =
      optional<double>(systemConfig, "shower_shape_tower_min_energy_gev", reco.isPP ? 0.070 : 0.0);

  const YAML::Node isolation = systemConfig["isolation"];
  const std::string isolationMethod = lower(optionalString(isolation, "method", reco.isPP ? "topocluster" : "sub1"));
  reco.useTopoclusterIsolation = reco.isPP && isolationMethod == "topocluster";
  reco.useSubtractedIsolation = reco.isAuAu;
  reco.topoclusterNode = tree.photon.topoclusterNode;

  if (reco.isAuAu)
  {
    const YAML::Node background = root["auau"]["background_subtraction"];
    reco.retowerFractionCut = optional<double>(background["retower_cemc"], "fraction_cut", 0.5);
    reco.firstPassSeedJetD = optional<int>(background["iteration_1"], "seed_jet_D", 3);
    reco.secondPassSeedPtGeV = optional<double>(background["iteration_2"], "seed_pt_min_gev", 7.0);
    reco.flowModulation = required<int>(background, "flow_modulation", "auau.background_subtraction");
    tree.calorimeter.backgroundFlowMode = reco.flowModulation;
  }

  // ---- jets ---------------------------------------------------------------

  tree.jets.truthMatchMaxDeltaR = optional<double>(root["jets"]["matching"], "max_delta_r", 0.30);

  const YAML::Node calibration = root["jets"]["calibration"];
  cond.applyJetEnergyScale = optional<bool>(calibration, "enabled", true);
  cond.jetEnergyScaleCdbKey = optionalString(calibration, "cdb_key", "JES_Calib_Default");
  cond.jetEnergyScalePayloadSha256 =
      lower(required<std::string>(calibration, "expected_payload_sha256", "jets.calibration"));
  cond.jetEnergyScaleUsesEmFraction = optional<bool>(calibration, "use_em_fraction_calibration", false);
  tree.jets.energyScalePayloadSha256 = cond.jetEnergyScalePayloadSha256;
  tree.jets.energyScaleUsesEmFraction = cond.jetEnergyScaleUsesEmFraction;

  const YAML::Node views = systemConfig["jets"]["views"];
  if (!views || !views.IsSequence())
  {
    fail("missing jet view list for the selected collision system");
  }
  for (const YAML::Node& view : views)
  {
    photonjet::JetNodeConfig jet;
    jet.radius = required<double>(view, "radius", "jet view");
    jet.view = parseJetView(required<std::string>(view, "view", "jet view"));
    jet.rawNode = required<std::string>(view, "raw_node", "jet view");
    jet.correctedNode = required<std::string>(view, "corrected_node", "jet view");
    jet.truthNode = optionalString(view, "truth_node");
    jet.inputIdentity = required<std::string>(view, "view", "jet view");
    jet.subtractionIdentity = jet.inputIdentity;
    tree.jets.nodes.push_back(std::move(jet));
  }

  // ---- centrality -----------------------------------------------------------

  tree.centrality.edges = intList(root["auau"]["centrality"]["analysis_bin_edges_percent"]);
  if (tree.centrality.edges.empty())
  {
    tree.centrality.edges = {0, 20, 50, 80};
  }
  tree.centrality.reconstructDataCentrality = reco.reconstructCentrality;

  if (reco.reconstructCentrality)
  {
    const YAML::Node data = root["auau"]["centrality"]["data"];
    tree.centrality.requireOwnRunPayload = optional<bool>(data, "require_own_run_payload", true);
    tree.centrality.calibrationTag = required<std::string>(data, "calibration_set", "auau.centrality.data");
    cond.centralityCalibrationSet = tree.centrality.calibrationTag;
    cond.centralityManifestSha256 =
        lower(required<std::string>(data, "calibration_manifest_sha256", "auau.centrality.data"));
    cond.centralityLocalPayloads = lower(optionalString(data, "source", "cdb")) != "cdb";

    // A run-scoped override table pins the exact payload triplet per run.
    YAML::Node resolved = data;
    if (const YAML::Node runs = root["auau"]["centrality"]["runs"])
    {
      if (const YAML::Node override = runs[std::to_string(job.run)])
      {
        resolved = override;
      }
    }
    cond.centralityPayloadRun = optional<int>(resolved, "payload_run", job.run);
    cond.centralityDivisions = expandEnvironment(optionalString(resolved["divisions"], "path"));
    cond.centralityDivisionsSha256 = lower(optionalString(resolved["divisions"], "sha256"));
    cond.centralityRunScale = expandEnvironment(optionalString(resolved["run_scale"], "path"));
    cond.centralityRunScaleSha256 = lower(optionalString(resolved["run_scale"], "sha256"));
    cond.centralityVertexScale = expandEnvironment(optionalString(resolved["vertex_scale"], "path"));
    cond.centralityVertexScaleSha256 = lower(optionalString(resolved["vertex_scale"], "sha256"));
    cond.centralityFallbackEvidenceSha256 = lower(optionalString(resolved, "fallback_evidence_sha256"));

    tree.centrality.divisionPayloadPath = cond.centralityDivisions;
    tree.centrality.runScalePayloadPath = cond.centralityRunScale;
    tree.centrality.vertexScalePayloadPath = cond.centralityVertexScale;
  }

  // ---- CEMC status and zero suppression, Au+Au data -------------------------

  if (reco.isAuAu && reco.isData)
  {
    const YAML::Node status = root["auau"]["cemc_status"];
    cond.cemcStatusGuard = optional<bool>(status, "validate_paired_input_status", true);
    cond.cemcOriginalCalorimeterInput = expandEnvironment(optionalString(status, "original_calorimeter_input"));
    cond.cemcStatusRecoveryAuthorized = optional<bool>(status, "allow_missing_status_recovery", false);
    cond.cemcStatusRecoveryPayload = expandEnvironment(optionalString(status["recovery_payload"], "path"));
    cond.cemcStatusRecoveryPayloadSha256 = lower(optionalString(status["recovery_payload"], "sha256"));
  }
  cond.cemcZeroSuppressionCrossCalibration =
      reco.isData && optional<bool>(root["calorimeter"]["cemc"], "zero_suppression_cross_calibration", true);

  // ---- archived double interaction -------------------------------------------

  if (reco.isArchivedDoubleInteraction)
  {
    cond.cdbGlobalTag = required<std::string>(profile, "cdb_global_tag", context);
    cond.cdbTimestamp = required<std::uint64_t>(profile, "cdb_timestamp", context);
    const YAML::Node rng = profile["historical_rng_contract"];
    reco.expectedPedestalSequence = required<int>(rng, "expected_pedestal_sequence", context);
    reco.forbidRecoConstsRandomSeed = optional<bool>(rng, "forbid_reco_consts_random_seed", true);
    for (const YAML::Node& seed : rng["seed_sequence"])
    {
      reco.historicalSeeds.push_back(seed.as<unsigned int>());
    }
  }
  else
  {
    cond.cdbTimestamp = job.run > 0 ? static_cast<std::uint64_t>(job.run)
                                    : static_cast<std::uint64_t>(kArchivedDoubleInteractionRun);
    cond.cdbGlobalTag = optionalString(profile, "cdb_global_tag");
  }

  // ---- weights ---------------------------------------------------------------

  if (const YAML::Node reweight = profile["vertex_reweighting"])
  {
    cond.applyVertexReweight = optional<bool>(reweight, "enabled", false);
    cond.vertexReweightFile = expandEnvironment(optionalString(reweight, "file"));
    cond.vertexReweightFileSha256 = lower(optionalString(reweight, "sha256"));
    cond.vertexReweightHistogram = optionalString(reweight, "histogram");
  }
  tree.weights.applyVertexReweight = cond.applyVertexReweight;
  tree.weights.vertexReweightFile = cond.vertexReweightFile;
  tree.weights.vertexReweightFileSha256 = cond.vertexReweightFileSha256;
  tree.weights.vertexReweightHistogram = cond.vertexReweightHistogram;

  // ---- output tables ----------------------------------------------------------

  const YAML::Node output = root["output"];
  tree.output.writePhotonCells = optional<bool>(output, "write_photon_cell_table", true);
  tree.output.writeIsolationConstituents = optional<bool>(output, "write_isolation_constituent_table", true);
  tree.output.writePhotonJetPairs = optional<bool>(output, "write_photon_jet_pair_table", true);

  // ---- retention and provenance ----------------------------------------------

  tree.retainEveryProducerEncounter = optional<bool>(root["capture"], "retain_every_producer_encounter", true);
  tree.requireEmbeddedMinimumBias = optional<bool>(profile, "require_embedded_minimum_bias", false);

  tree.provenance.productionTag = required<std::string>(root["contract"], "name", "contract");
  tree.provenance.softwareRelease = optionalString(root["runtime"], "sphenix_release");
  tree.provenance.producerGitCommit = optionalString(root["runtime"], "producer_git_commit");
  tree.provenance.producerSourceSha256 = optionalString(root["runtime"], "producer_source_sha256");
  tree.provenance.producerLibrarySha256 = optionalString(root["runtime"], "producer_library_sha256");
  tree.provenance.macroSha256 = optionalString(root["runtime"], "macro_sha256");
  tree.provenance.configurationSha256 = plan.configurationSha256;
  tree.provenance.configurationText = plan.configurationText;
  tree.provenance.calibrationManifestSha256 = cond.centralityManifestSha256;

  // ---- physical inputs ----------------------------------------------------------

  plan.inputs = loadInputs(job, reco);

  if (reco.truthJetMode == TruthJetMode::Auto)
  {
    reco.truthJetMode = hasInputStream(plan, "TRUTH_JETS") ? TruthJetMode::DST : TruthJetMode::Build;
  }

  Validate(plan);
  return plan;
}

//
// Validate
//
void Validate(const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  const Conditions& cond = plan.conditions;

  if (plan.inputs.empty()) fail("no input streams were resolved");
  if (plan.job.firstEntry < 0 || plan.job.firstEntry > std::numeric_limits<int>::max() ||
      plan.job.nEvents < -1 || plan.job.nEvents == 0 || plan.job.nEvents > std::numeric_limits<int>::max())
    fail("invalid event range: use -1 for all events, or a positive bounded count");
  // One synchronized source bundle per invocation. Multiple list rows need a
  // source-boundary observer before they can safely share one SourceRecord.
  for (const auto& stream : plan.inputs)
    if (stream.files.size() != 1) fail("one source bundle per invocation is required: " + stream.name);
  if (cond.cdbGlobalTag.empty()) fail("profile must bind an explicit CDB global tag");
  if (!cond.applyJetEnergyScale) fail("canonical reconstructed jets require JES");
  if (cond.jetEnergyScaleCdbKey != "JES_Calib_Default") fail("JES key must match the key consumed by JetCalib");
  if (cond.applyVertexReweight && !validSha256(cond.vertexReweightFileSha256)) fail("vertex weight payload requires SHA256");
  for (const InputStream& stream : plan.inputs)
  {
    if (stream.required && stream.files.empty())
    {
      fail("required input stream is empty: " + stream.name);
    }
  }

  if (reco.isData && plan.tree.simulationRole != SimulationRole::None)
  {
    fail("a data profile cannot carry a simulation role");
  }
  if (reco.isSimulation && plan.tree.simulationRole == SimulationRole::None)
  {
    fail("a simulation profile must name its sample role");
  }

  if (reco.isArchivedDoubleInteraction)
  {
    if (!reco.isPP || !reco.isSimulation) fail("archived double interaction is p+p simulation only");
    // TODO(CODE): port the qualified G4-hit -> waveform/tower sequence and its
    // seed consumers before enabling this lane. CaloTowerCalib alone cannot
    // turn archived G4 hits into towers. Never manufacture an empty success.
    fail("archived DI waveform/tower registration is not yet bound");
    if (reco.truthJetMode != TruthJetMode::DST) fail("archived double interaction requires shipped truth jets");
    if (cond.cdbTimestamp != static_cast<std::uint64_t>(kArchivedDoubleInteractionRun))
    {
      fail("archived double interaction requires the run-28 conditions timestamp");
    }
  }

  if ((reco.isData || reco.isEmbedded) && plan.job.run <= kDataRunThreshold)
  {
    fail("data and embedded production require a real data run number");
  }

  if (reco.reconstructCentrality)
  {
    if (cond.centralityCalibrationSet != kNominalCentralitySet)
    {
      fail("Au+Au centrality is not bound to the nominal frozen calibration set");
    }
    if (cond.centralityManifestSha256 != kNominalCentralityManifestSha256)
    {
      fail("Au+Au centrality manifest SHA256 does not match the frozen nominal set");
    }
    if (cond.centralityLocalPayloads)
    {
      if (cond.centralityPayloadRun <= 0) fail("local centrality requires a positive payload run");
      if (cond.centralityDivisions.empty() || cond.centralityRunScale.empty() || cond.centralityVertexScale.empty())
      {
        fail("local centrality requires divisions, run-scale and vertex-scale payloads");
      }
      if (!validSha256(cond.centralityDivisionsSha256) || !validSha256(cond.centralityRunScaleSha256) ||
          !validSha256(cond.centralityVertexScaleSha256))
      {
        fail("local centrality requires an exact SHA256 for each of the three payloads");
      }
      if (cond.centralityPayloadRun != plan.job.run)
      {
        if (cond.centralityPayloadRun != kQualifiedCentralityFallbackRun)
        {
          fail("only the qualified run-68144 centrality fallback is permitted");
        }
        if (!validSha256(cond.centralityFallbackEvidenceSha256))
        {
          fail("the centrality fallback requires pinned qualification evidence");
        }
      }
      else if (!cond.centralityFallbackEvidenceSha256.empty())
      {
        fail("own-run centrality cannot carry a fallback approval");
      }
    }
  }

  if (cond.cemcStatusGuard)
  {
    if (cond.cemcOriginalCalorimeterInput.empty())
    {
      fail("the CEMC status guard requires the exact original calorimeter input file");
    }
    if (cond.cemcStatusRecoveryAuthorized &&
        (cond.cemcStatusRecoveryPayload.empty() || !validSha256(cond.cemcStatusRecoveryPayloadSha256)))
    {
      fail("CEMC status recovery requires a pinned run-specific payload and its SHA256");
    }
  }

  if (cond.applyJetEnergyScale)
  {
    if (cond.jetEnergyScaleUsesEmFraction)
    {
      fail("the nominal jet energy scale is the legacy method without the electromagnetic-fraction calibration");
    }
    if (cond.jetEnergyScalePayloadSha256 != kNominalJetEnergyScaleSha256)
    {
      fail("jet energy-scale payload SHA256 differs from the nominal common p+p v6 payload");
    }
  }

  if (cond.applyVertexReweight)
  {
    if (!reco.isSimulation) fail("vertex reweighting applies to simulation only");
    if (cond.vertexReweightFile.empty() || cond.vertexReweightHistogram.empty())
    {
      fail("vertex reweighting requires a file and a histogram name");
    }
  }

  if (reco.isAuAu && (reco.flowModulation < 0 || reco.flowModulation > 3))
  {
    fail("Au+Au flow modulation must be resolved to the accepted production value");
  }
}
// Print
//
void Print(const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  const Conditions& cond = plan.conditions;

  const char* truthMode = "none";
  switch (reco.truthJetMode)
  {
    case TruthJetMode::None: truthMode = "none"; break;
    case TruthJetMode::Auto: truthMode = "auto"; break;
    case TruthJetMode::DST: truthMode = "dst"; break;
    case TruthJetMode::Build: truthMode = "build"; break;
  }

  std::cout << "\n"
            << "============================================================\n"
            << " PhotonJetTree production plan\n"
            << "============================================================\n"
            << " profile              : " << plan.job.profile << "\n"
            << " system               : " << (reco.isPP ? "pp" : "AuAu") << "\n"
            << " data kind            : " << (reco.isData ? "DATA" : "SIM") << "\n"
            << " embedded             : " << (reco.isEmbedded ? "yes" : "no") << "\n"
            << " archived DI          : " << (reco.isArchivedDoubleInteraction ? "yes" : "no") << "\n"
            << " run / segment        : " << plan.job.run << " / " << plan.job.segment << "\n"
            << " first entry / events : " << plan.job.firstEntry << " / " << plan.job.nEvents << "\n"
            << " input list           : " << plan.job.inputList << "\n"
            << " output               : " << plan.job.outputFile << "\n"
            << " config sha256        : " << plan.configurationSha256 << "\n"
            << " truth jets           : " << truthMode << "\n"
            << " jet collections      : " << plan.tree.jets.nodes.size() << "\n"
            << " CDB timestamp / tag  : " << cond.cdbTimestamp << " / "
            << (cond.cdbGlobalTag.empty() ? "<preconfigured>" : cond.cdbGlobalTag) << "\n"
            << " calorimeter          : " << (reco.reconstructCalorimeter ? "reconstructed here" : "from input") << "\n"
            << " JES                  : " << (cond.applyJetEnergyScale ? "legacy, once, " + cond.jetEnergyScalePayload : "disabled") << "\n";

  if (reco.reconstructCentrality)
  {
    std::cout << " centrality           : " << cond.centralityCalibrationSet
              << (cond.centralityLocalPayloads ? " local payload run " + std::to_string(cond.centralityPayloadRun)
                                               : " from CDB")
              << "\n";
  }
  if (cond.cemcStatusGuard)
  {
    std::cout << " CEMC status guard    : " << cond.cemcOriginalCalorimeterInput
              << (cond.cemcStatusRecoveryAuthorized ? " (recovery authorised)" : "") << "\n";
  }
  if (cond.applyVertexReweight)
  {
    std::cout << " vertex reweight      : " << cond.vertexReweightFile << ":" << cond.vertexReweightHistogram << "\n";
  }

  std::cout << " input streams        : " << plan.inputs.size() << "\n";
  for (const InputStream& stream : plan.inputs)
  {
    std::cout << "   " << std::left << std::setw(12) << stream.name
              << " files=" << stream.files.size()
              << " sync=" << (stream.synchronized ? "yes" : "no") << "\n";
  }
  std::cout << "============================================================\n" << std::endl;
}

}  // namespace production
}  // namespace photonjet
