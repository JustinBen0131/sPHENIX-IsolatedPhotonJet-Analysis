/**
 * @file Production.cc
 * Resolve a source Job and reconstruction preset into a validated Plan.
 * Configuration values bind algorithms and payloads without duplicating them.
 * Missing conditions, incompatible modes and identity ambiguity fail closed.
 * Setup.cc performs the conditions and Fun4All registration steps.
 */

#include "Production.h"
#include "BuildIdentity.h"

#include <yaml-cpp/yaml.h>
#include <dlfcn.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <filesystem>
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

// ============================================================================
// Shared production utilities
// ============================================================================

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
  std::transform(
      value.begin(),
      value.end(),
      value.begin(),
      [](const unsigned char c)
      {
        return static_cast<char>(std::tolower(c));
      });

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


// ============================================================================
// Frozen nominal production identities
//
// These values are part of the accepted production contract rather than
// ordinary tunable YAML parameters. The configuration must reproduce them
// exactly before nominal production is admitted.
//
// Keeping the accepted identities in code provides an independent guard against
// a configuration edit silently moving the nominal centrality or JES basis.
// ============================================================================

constexpr const char* kNominalCentralitySet =
    "sam_fit_20260913_v1";

constexpr const char* kNominalCentralityManifestSha256 =
    "4a11e0e10ae9b82062325c3275679935da2259f4c91f796950a6b5ddd2e1501c";

constexpr const char* kNominalJetEnergyScaleSha256 =
    "76a8788fdb4e0b3859d01361f884e295ea609bb8c1cc188105389214db2de27d";

constexpr int kQualifiedCentralityFallbackRun = 68144;


// ============================================================================
// Configuration-file helpers
// ============================================================================

/**
 * Read the exact configuration bytes.
 *
 * The same bytes are retained in Plan and propagated into PhotonJetTree
 * provenance so an output file can later identify the configuration that
 * produced it without relying on a mutable external YAML file.
 */
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


/**
 * Expand ${NAME} expressions in configured filesystem paths.
 *
 * An unset variable is a hard configuration error. Replacing it with an empty
 * string could accidentally redirect a payload lookup or turn a required path
 * into a silently missing optional input.
 */
std::string expandEnvironment(const std::string& input)
{
  std::string output;

  for (std::size_t i = 0; i < input.size();)
  {
    if (input[i] == '$' &&
        i + 1 < input.size() &&
        input[i + 1] == '{')
    {
      const std::size_t close =
          input.find('}', i + 2);

      if (close == std::string::npos)
      {
        fail(
            "unterminated environment variable in path: " +
            input);
      }

      const std::string key =
          input.substr(
              i + 2,
              close - i - 2);

      const char* value =
          std::getenv(key.c_str());

      if (!value)
      {
        fail(
            "required environment variable is unset: " +
            key);
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


// ============================================================================
// YAML access helpers
//
// Required values fail immediately when absent or null. Optional values use an
// explicit caller-provided fallback. These helpers keep configuration semantics
// uniform throughout LoadPlan.
// ============================================================================

template <class T>
T required(const YAML::Node& node,
           const std::string& key,
           const std::string& context)
{
  if (!node ||
      !node[key] ||
      node[key].IsNull())
  {
    fail(
        "missing required configuration key '" +
        key +
        "' in " +
        context);
  }

  return node[key].as<T>();
}


template <class T>
T optional(const YAML::Node& node,
           const std::string& key,
           const T& fallback)
{
  if (!node ||
      !node[key] ||
      node[key].IsNull())
  {
    return fallback;
  }

  return node[key].as<T>();
}


std::string optionalString(
    const YAML::Node& node,
    const std::string& key,
    const std::string& fallback = {})
{
  if (!node ||
      !node[key] ||
      node[key].IsNull())
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
      result.push_back(
          value.as<std::string>());
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
      result.push_back(
          value.as<double>());
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
      result.push_back(
          value.as<int>());
    }
  }

  return result;
}


// ============================================================================
// Configuration vocabulary parsers
// ============================================================================

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

  if (view == "auau_sub1" ||
      view == "sub1")
  {
    return photonjet::JetView::AuAuSub1;
  }

  if (view == "auau_nosub" ||
      view == "nosub")
  {
    return photonjet::JetView::AuAuNoSub;
  }

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


// ============================================================================
// Plan construction
//
// LoadPlan is the single configuration entrypoint used by the steering macro.
//
// Inputs
// ------
//   yamlFile
//     shared reconstruction/storage configuration
//
//   job
//     source-specific values already resolved by the production frontend
//
// Output
// ------
//   Plan
//     typed source binding, PhotonJetTree configuration, physical input streams,
//     reconstruction steering, conditions steering, and exact configuration
//     provenance
//
// LoadPlan resolves values but does not bind CDB/local payload paths itself;
// ConfigureConditions performs that subsequent run-dependent step.
// ============================================================================

Plan LoadPlan(const std::string& yamlFile,
              const Job& job)
{
  if (job.profile.empty())
  {
    fail("production profile is empty");
  }

  if (job.inputList.empty())
  {
    fail("input-list path is empty");
  }

  if (job.outputFile.empty())
  {
    fail("output file name is empty");
  }

  Plan plan;

  plan.job = job;

  /*
   * Preserve both the exact configuration bytes and their digest before parsing.
   * PhotonJetTree receives these through its provenance configuration.
   */
  plan.configurationText =
      readTextFile(yamlFile);

  plan.configurationSha256 =
      fileSha256(yamlFile);

  const YAML::Node root =
      YAML::Load(plan.configurationText);

  const YAML::Node profile =
      root["profiles"][job.profile];

  if (!profile)
  {
    fail(
        "profile '" +
        job.profile +
        "' is absent from " +
        yamlFile);
  }

  const std::string context =
      "profiles." + job.profile;

  Reconstruction& reco =
      plan.reconstruction;

  Conditions& cond =
      plan.conditions;

  photonjet::Config& tree =
      plan.tree;


  // ==========================================================================
  // Production mode
  //
  // Resolve the four fundamental mode dimensions first. All later stage
  // registration is derived from these flags rather than from repeated profile
  // string tests.
  // ==========================================================================

  const std::string system =
      lower(
          required<std::string>(
              profile,
              "system",
              context));

  const std::string dataKind =
      lower(
          required<std::string>(
              profile,
              "data_kind",
              context));

  const std::string role =
      lower(
          required<std::string>(
              profile,
              "simulation_role",
              context));

  const std::string inputMode =
      lower(
          required<std::string>(
              profile,
              "input_mode",
              context));

  reco.isPP =
      system == "pp";

  reco.isAuAu =
      system == "auau";

  reco.isData =
      dataKind == "data";

  reco.isSimulation =
      dataKind == "simulation";

  reco.isEmbedded =
      reco.isAuAu &&
      reco.isSimulation;

  reco.isArchivedDoubleInteraction =
      inputMode == "archived_g4_only";

  if (reco.isPP == reco.isAuAu)
  {
    fail(
        "profile must name exactly one collision system");
  }

  if (reco.isData == reco.isSimulation)
  {
    fail(
        "profile must name exactly one data kind");
  }

  tree.system =
      reco.isPP
          ? CollisionSystem::PP
          : CollisionSystem::AuAu;

  tree.dataKind =
      reco.isData
          ? DataKind::Data
          : DataKind::Simulation;

  tree.inputMode =
      reco.isArchivedDoubleInteraction
          ? InputMode::ArchivedG4Only
          : InputMode::ReconstructedDST;

  if (role == "photon_jet")
  {
    tree.simulationRole =
        SimulationRole::PhotonJet;
  }
  else if (role == "inclusive_jet")
  {
    tree.simulationRole =
        SimulationRole::InclusiveJet;
  }
  else if (role == "minimum_bias")
  {
    tree.simulationRole =
        SimulationRole::MinimumBias;
  }
  else
  {
    tree.simulationRole =
        SimulationRole::None;
  }

  tree.outputFile =
      job.outputFile;

  if (inputMode != "archived_g4_only" &&
      inputMode != "reconstructed_dst")
  {
    fail(
        "unknown input mode: " +
        inputMode);
  }

  if (role != "none" &&
      role != "photon_jet" &&
      role != "inclusive_jet" &&
      role != "minimum_bias")
  {
    fail(
        "unknown simulation role: " +
        role);
  }


  // ==========================================================================
  // Source binding
  //
  // Per-job source coordinates come from produce.py and remain distinct from
  // shared reconstruction configuration. The profile supplies only a default
  // dataset token when the Job did not explicitly provide one.
  // ==========================================================================

  tree.source.lane = job.profile;

  tree.source.dataset =
      job.dataset.empty()
          ? optionalString(
                profile,
                "dataset_token")
          : job.dataset;

  tree.source.sample =
      job.sample;

  tree.source.period =
      job.period;

  tree.source.siDiRole =
      job.siDiRole;

  tree.source.inputUriSha256 =
      job.inputUriSha256;

  tree.source.inputFileSha256 =
      job.inputFileSha256;

  tree.source.sourceManifestSha256 =
      job.sourceManifestSha256;

  tree.source.run =
      job.run;

  tree.source.segment =
      job.segment;

  tree.source.sourceFileOrdinal =
      job.sourceFileOrdinal;

  tree.source.firstEntry =
      job.firstEntry;

  tree.source.maxEvents =
      job.nEvents;


  // ==========================================================================
  // Canonical input-node vocabulary
  //
  // These names tell PhotonJetTree where to read the reconstructed event.
  // Setup.cc is responsible for creating the required nodes
  // when reconstruction occurs in this pass.
  // ==========================================================================

  const YAML::Node nodes =
      root["nodes"];

  tree.nodes.eventHeader =
      optionalString(
          nodes,
          "event_header",
          "EventHeader");

  tree.nodes.gl1Packet =
      optionalString(
          nodes,
          "gl1_packet",
          "GL1Packet");

  tree.nodes.globalVertexMap =
      optionalString(
          nodes,
          "global_vertex_map",
          "GlobalVertexMap");

  tree.nodes.centralityInfo =
      optionalString(
          nodes,
          "centrality_info",
          "CentralityInfo");

  tree.nodes.minimumBiasInfo =
      optionalString(
          nodes,
          "minimum_bias_info",
          "MinimumBiasInfo");

  tree.nodes.mbdOut =
      optionalString(
          nodes,
          "mbd_out",
          "MbdOut");

  tree.nodes.mbdPmts =
      optionalString(
          nodes,
          "mbd_pmts",
          "MbdPmtContainer");

  tree.nodes.mbdPmtsFallback = optionalString(nodes, "mbd_pmts_fallback", "MbdPmtContainer_data");

  tree.nodes.cemcTowers =
      optionalString(
          nodes,
          "cemc_towers",
          "TOWERINFO_CALIB_CEMC");

  tree.nodes.ihcalTowers =
      optionalString(
          nodes,
          "ihcal_towers",
          "TOWERINFO_CALIB_HCALIN");

  tree.nodes.ohcalTowers =
      optionalString(
          nodes,
          "ohcal_towers",
          "TOWERINFO_CALIB_HCALOUT");

  tree.nodes.cemcGeometry =
      optionalString(
          nodes,
          "cemc_geometry",
          "TOWERGEOM_CEMC");

  tree.nodes.truthInfo =
      optionalString(
          nodes,
          "truth_info",
          "G4TruthInfo");

  tree.nodes.hepmcEventMap =
      optionalString(
          nodes,
          "hepmc_event_map",
          "PHHepMCGenEventMap");


  // ==========================================================================
  // Calorimeter event-state capture
  // ==========================================================================

  tree.calorimeter.requireGoodTowers =
      optional<bool>(
          root["calorimeter"]["event_energy"],
          "require_good_towers",
          true);


  // ==========================================================================
  // Photon capture configuration
  //
  // These values define the broad model-independent storage domain and the
  // reconstruction observables retained for later PhotonID/analysis stages.
  // They are not the final measurement selection.
  // ==========================================================================

  const YAML::Node photons =
      root["photons"];

  const YAML::Node domain =
      photons["capture_domain"];

  tree.photon.minEtGeV =
      required<double>(
          domain,
          "min_et_gev",
          "photons.capture_domain");

  tree.photon.maxEtGeV =
      required<double>(
          domain,
          "max_et_gev",
          "photons.capture_domain");

  tree.photon.maxAbsEta =
      required<double>(
          domain,
          "max_abs_eta",
          "photons.capture_domain");

  tree.photon.objectVertexAbsZMaxCm =
      optional<double>(
          domain,
          "object_vertex_abs_z_max_cm",
          tree.photon.objectVertexAbsZMaxCm);

  tree.photon.rawClusterNode =
      optionalString(
          photons,
          "source_cluster_node",
          "CLUSTERINFO_CEMC");

  tree.photon.photonNode =
      optionalString(
          photons,
          "output_photon_node",
          "PHOTONCLUSTER_CEMC");

  tree.photon.topoclusterNode =
      optionalString(
          photons,
          "topocluster_node",
          "TOPOCLUSTER_ALLCALO");

  const std::vector<double> radii =
      doubleList(
          photons["isolation"]["radii"]);

  if (!radii.empty())
  {
    tree.photon.isolationRadii =
        radii;
  }

  tree.photon.truthIsolationSignalRadius =
      optional<double>(
          photons["isolation"]["truth"],
          "signal_radius",
          0.30);

  tree.photon.truthIsolationMaxEtGeV =
      required<double>(
          photons["isolation"]["truth"],
          "max_et_gev",
          "photons.isolation.truth");

  tree.photon.truthIsolationCoreRadius =
      required<double>(
          photons["isolation"]["truth"],
          "self_removal_radius",
          "photons.isolation.truth");

  tree.photon.timingSampleNs =
      required<double>(
          photons["timing"],
          "sample_ns",
          "photons.timing");

  const std::vector<std::string> definitions =
      stringList(
          photons["shower_shapes"]["definitions"]);

  if (!definitions.empty())
  {
    tree.showerDefinitions =
        definitions;
  }

  // Optional exact named models. No implicit file, feature order or physics
  // domain is inferred from a score name; generic production defaults stay off.
  if(photons["models"])for(const auto& node:photons["models"]) {
    PhotonModelConfig model;
    model.scoreName=required<std::string>(node,"score_name","photons.models");
    model.file=required<std::string>(node,"file","photons.models");
    model.sha256=required<std::string>(node,"sha256","photons.models");
    model.features=stringList(node["features"]);
    if(model.features.empty())fail("photons.models requires an explicit feature order");
    // Each invocation bound is optional; an absent bound leaves that side open.
    model.minEt=optional<float>(node,"min_et",NAN);
    model.maxEt=optional<float>(node,"max_et",NAN);
    model.maxAbsEta=optional<float>(node,"max_abs_eta",NAN);
    model.unevaluatedScore=optional<float>(node,"unevaluated_score",NAN);
    model.configuredCut=optional<double>(node,"configured_cut",NAN);
    tree.photon.models.push_back(model);
  }

  // Active photon preselection (pho_presel_active); only the PPG12 newPPG12
  // definition is supported, on a captured NPB score with a configured cut.
  if (const YAML::Node pre = photons["preselection"])
  {
    const std::string where = "photons.preselection";
    PreselectionConfig& preselection = tree.photon.preselection;
    preselection.variant = required<std::string>(pre, "variant", where);
    if (preselection.variant != "newPPG12") fail("unsupported " + where + ".variant: " + preselection.variant);
    preselection.npbScoreName = required<std::string>(pre, "npb_score", where);
    if (preselection.npbScoreName != "npb_score" && preselection.npbScoreName != "auau_npb_score")
      fail(where + ".npb_score must name a captured NPB score");
    const auto window = [&](const char* key, double& low, double& high)
    {
      const YAML::Node range = pre[key];
      if (!range || !range.IsSequence() || range.size() != 2)
        fail(where + "." + key + " must be [low, high]");
      low = range[0].as<double>();
      high = range[1].as<double>();
      if (!(low < high)) fail("empty " + where + "." + key);
    };
    window("et1", preselection.et1Min, preselection.et1Max);
    window("e32_over_e35", preselection.e32OverE35Min, preselection.e32OverE35Max);
    preselection.e11OverE33Max = required<double>(pre, "e11_over_e33_max", where);
    const auto model = std::find_if(tree.photon.models.begin(), tree.photon.models.end(),
        [&](const PhotonModelConfig& value) { return value.scoreName == preselection.npbScoreName; });
    if (model == tree.photon.models.end() || !std::isfinite(model->configuredCut))
      fail(where + " needs the registered " + preselection.npbScoreName + " model with its configured_cut");
  }


  // ==========================================================================
  // Reconstruction-stage selection
  //
  // Translate the high-level profile into the stages Setup.cc
  // must register. The same PhotonJetTree implementation is used for every mode;
  // only the upstream reconstruction graph changes.
  // ==========================================================================

  reco.inputLayout =
      reco.isArchivedDoubleInteraction
          ? InputLayout::ArchivedDoubleInteraction
          : reco.isData
                ? InputLayout::PairedData
                : InputLayout::ReconstructedSimulation;

  reco.truthJetMode =
      parseTruthJetMode(
          optional<std::string>(
              profile,
              "truth_jets_mode",
              reco.isSimulation
                  ? "auto"
                  : "none"));

  reco.reconstructEventGeometry =
      !reco.isEmbedded;

  reco.reconstructZdc =
      reco.isData;

  reco.reconstructCalorimeter =
      reco.isData ||
      reco.isEmbedded ||
      reco.isArchivedDoubleInteraction;

  reco.runCalorimeterStatusSkimmer =
      reco.isData;

  reco.setTowerStatus =
      reco.reconstructCalorimeter &&
      (!reco.isEmbedded ||
       optional<bool>(
           profile,
           "force_embedded_tower_status",
           false));

  reco.reconstructCentrality =
      reco.isAuAu &&
      reco.isData &&
      optional<bool>(
          profile,
          "reconstruct_centrality",
          true);

  reco.runMinimumBiasClassifier =
      reco.reconstructCentrality &&
      optional<bool>(
          profile,
          "require_minimum_bias_classifier",
          true);

  reco.reconstructAuAuBackground =
      reco.isAuAu;

  reco.reconstructJets =
      true;

  reco.reconstructPhotons =
      true;


  // ==========================================================================
  // Calorimeter and photon reconstruction parameters
  // ==========================================================================

  const YAML::Node clusterBuilder =
      root["calorimeter"]["cluster_builder"];

  reco.clusterThresholdGeV =
      required<double>(
          clusterBuilder,
          "threshold_energy_gev",
          "calorimeter.cluster_builder");

  reco.clusterProfile =
      expandEnvironment(
          required<std::string>(
              clusterBuilder["profile"],
              "path",
              "calorimeter.cluster_builder.profile"));

  reco.simulationHotTowerMapKey =
      optionalString(
          root["calorimeter"]["tower_status"],
          "simulation_hot_map_key",
          "CEMC_hotTowers_status");

  const YAML::Node systemConfig =
      reco.isPP
          ? root["pp"]
          : root["auau"];

  reco.photonMinEtGeV =
      tree.photon.minEtGeV;

  reco.photonUseVertexCut =
      optional<bool>(
          photons["builder"],
          "use_vertex_cut",
          false);

  reco.photonVertexCutCm =
      optional<double>(
          photons["builder"],
          "vertex_cut_cm",
          60.0);

  reco.showerShapeTowerMinEnergyGeV =
      optional<double>(
          systemConfig,
          "shower_shape_tower_min_energy_gev",
          reco.isPP
              ? 0.070
              : 0.0);

  const YAML::Node isolation =
      systemConfig["isolation"];

  const std::string isolationMethod =
      lower(
          optionalString(
              isolation,
              "method",
              reco.isPP
                  ? "topocluster"
                  : "sub1"));

  reco.useTopoclusterIsolation =
      reco.isPP &&
      isolationMethod == "topocluster";

  reco.useSubtractedIsolation =
      reco.isAuAu;

  reco.topoclusterNode =
      tree.photon.topoclusterNode;

  // --------------------------------------------------------------------------
  // Isolation witnesses
  //
  // The historical working point recorded with each retained cone. pp lines
  // are in photon ET; Au+Au lines are in the truncated integer centrality
  // inside centrality_range, with the photon-ET fallback outside it.
  // --------------------------------------------------------------------------
  {
    const std::string where = std::string(reco.isPP ? "pp" : "auau") + ".isolation.witness";
    const YAML::Node witness = isolation["witness"];
    if (!witness) fail("missing required configuration block " + where);
    auto line = [&](const std::string& key)
    {
      const std::string context = where + "." + key;
      IsolationLine value;
      value.interceptGeV = required<double>(witness[key], "intercept_gev", context);
      value.slope = required<double>(witness[key], "slope", context);
      value.sidebandGapGeV = required<double>(witness[key], "sideband_gap_gev", context);
      if (!std::isfinite(value.interceptGeV) || !std::isfinite(value.slope) ||
          !std::isfinite(value.sidebandGapGeV) || value.sidebandGapGeV < 0.0)
      {
        fail("invalid isolation witness line " + context);
      }
      return value;
    };
    if (tree.photon.isolationRadii != std::vector<double>{0.30, 0.40})
      fail("isolation witnesses are defined for the retained radii 0.30 and 0.40");
    IsolationWitnessConfig& config = tree.photon.isolationWitness;
    config.lines = {line("r03"), line("r04")};
    if (reco.isAuAu)
    {
      config.fallback = line("fallback");
      const YAML::Node range = witness["centrality_range"];
      if (!range || !range.IsSequence() || range.size() != 2)
        fail(where + ".centrality_range must be [minimum, maximum)");
      config.centralityMin = range[0].as<double>();
      config.centralityMax = range[1].as<double>();
      if (!(config.centralityMin < config.centralityMax))
        fail("empty " + where + ".centrality_range");
    }
  }

  // --------------------------------------------------------------------------
  // Producer photon-identification score
  //
  // Bound to a canonical PhotonID model package: its file, digest and input
  // order come from the package manifest, never from this preset. The preset
  // supplies only the historical invocation (score name, minimum ET,
  // unevaluated sentinel) and the recorded working-point witnesses.
  // --------------------------------------------------------------------------
  if (const YAML::Node score = systemConfig["photon_id"]["production_score"])
  {
    const std::string where = std::string(reco.isPP ? "pp" : "auau") + ".photon_id.production_score";
    const auto folder = std::filesystem::path(yamlFile).parent_path().parent_path().parent_path() /
                        "PhotonID" / "models" / (reco.isPP ? "PP" : "AuAu") /
                        required<std::string>(score, "model", where);
    const YAML::Node manifest = YAML::LoadFile((folder / "manifest.json").string());
    const YAML::Node definition = manifest["definition"];
    const ShowerViewDefinition& primary = showerViewDefinition("H70");

    PhotonModelConfig model;
    model.scoreName = required<std::string>(score, "score_name", where);
    model.file = (folder / "model.tmva.root").string();
    model.sha256 = required<std::string>(manifest["artifacts"], "model.tmva.root", where + " manifest");
    model.features = stringList(definition["features"]["order"]);
    model.minEt = required<float>(score, "min_et", where);
    model.unevaluatedScore = required<float>(score, "unevaluated_score", where);
    if (!validSha256(model.sha256) ||
        model.sha256 != required<std::string>(definition["model"], "sha256", where + " manifest") ||
        fileSha256(model.file) != model.sha256)
    {
      fail(where + ": model file differs from its PhotonID manifest");
    }
    if (model.features.empty() ||
        required<std::string>(definition["features"], "view", where + " manifest") != primary.name ||
        required<std::string>(definition["features"], "view_semantic_sha256", where + " manifest") != primary.semanticSha256)
    {
      fail(where + ": the builder evaluates the primary H70 view; the package names another input");
    }
    if (definition["domain"]["et_min"] &&
        definition["domain"]["et_min"].as<double>() != static_cast<double>(model.minEt))
    {
      fail(where + ": min_et differs from the package domain");
    }
    for (const PhotonModelConfig& other : tree.photon.models)
      if (other.scoreName == model.scoreName) fail(where + ": score name registered twice");

    ProductionScoreConfig& production = tree.photon.productionScore;
    production.scoreName = model.scoreName;
    production.modelSha256 = model.sha256;
    production.showerDefinition = primary.name;
    const YAML::Node validated = score["validated_et"];
    if (!validated || !validated.IsSequence() || validated.size() != 2)
      fail(where + ".validated_et must be [minimum, maximum)");
    production.validatedEtMin = validated[0].as<double>();
    production.validatedEtMax = validated[1].as<double>();
    if (reco.isAuAu)
    {
      production.validatedCentralityMin = tree.photon.isolationWitness.centralityMin;
      production.validatedCentralityMax = tree.photon.isolationWitness.centralityMax;
    }
    if (const YAML::Node points = score["working_points"])
    {
      for (const YAML::Node& edge : points["et_edges"]) production.workingPointEtEdges.push_back(edge.as<double>());
      const std::array<const char*, 3> names{{"wp70", "wp80", "wp90"}};
      for (std::size_t index = 0; index < names.size(); ++index)
        for (const YAML::Node& value : points[names[index]]) production.workingPoints[index].push_back(value.as<double>());
      const auto& edges = production.workingPointEtEdges;
      bool ordered = edges.size() >= 2;
      for (std::size_t index = 1; index < edges.size(); ++index) ordered &= edges[index - 1] < edges[index];
      for (const auto& values : production.workingPoints)
      {
        ordered &= values.size() + 1 == edges.size();
        for (const double value : values) ordered &= std::isfinite(value);
      }
      if (!ordered) fail(where + ".working_points must give ascending ET edges and one finite value per bin");
    }
    tree.photon.models.push_back(model);
  }


  // ==========================================================================
  // Au+Au underlying-event reconstruction
  // ==========================================================================

  if (reco.isAuAu)
  {
    const YAML::Node background =
        root["auau"]["background_subtraction"];

    reco.retowerFractionCut =
        optional<double>(
            background["retower_cemc"],
            "fraction_cut",
            0.5);

    reco.firstPassSeedJetD =
        optional<int>(
            background["iteration_1"],
            "seed_jet_D",
            3);

    reco.secondPassSeedPtGeV =
        optional<double>(
            background["iteration_2"],
            "seed_pt_min_gev",
            7.0);

    reco.flowModulation =
        required<int>(
            background,
            "flow_modulation",
            "auau.background_subtraction");

    tree.calorimeter.backgroundFlowMode =
        reco.flowModulation;
  }


  // ==========================================================================
  // Reconstructed jets and jet-energy-scale contract
  // ==========================================================================

  tree.jets.truthMatchMaxDeltaR =
      optional<double>(
          root["jets"]["matching"],
          "max_delta_r",
          0.30);

  const YAML::Node calibration =
      root["jets"]["calibration"];

  cond.applyJetEnergyScale =
      optional<bool>(
          calibration,
          "enabled",
          true);

  cond.jetEnergyScaleCdbKey =
      optionalString(
          calibration,
          "cdb_key",
          "JES_Calib_Default");

  cond.jetEnergyScalePayloadSha256 =
      lower(
          required<std::string>(
              calibration,
              "expected_payload_sha256",
              "jets.calibration"));

  cond.jetEnergyScaleUsesEmFraction =
      optional<bool>(
          calibration,
          "use_em_fraction_calibration",
          false);

  tree.jets.energyScalePayloadSha256 =
      cond.jetEnergyScalePayloadSha256;

  tree.jets.energyScaleUsesEmFraction =
      cond.jetEnergyScaleUsesEmFraction;

  const YAML::Node views =
      systemConfig["jets"]["views"];

  if (!views || !views.IsSequence())
  {
    fail(
        "missing jet view list for the selected collision system");
  }

  for (const YAML::Node& view : views)
  {
    photonjet::JetNodeConfig jet;

    jet.radius =
        required<double>(
            view,
            "radius",
            "jet view");

    jet.view =
        parseJetView(
            required<std::string>(
                view,
                "view",
                "jet view"));

    jet.rawNode =
        required<std::string>(
            view,
            "raw_node",
            "jet view");

    jet.correctedNode =
        required<std::string>(
            view,
            "corrected_node",
            "jet view");

    jet.truthNode =
        optionalString(
            view,
            "truth_node");

    jet.inputIdentity =
        required<std::string>(
            view,
            "view",
            "jet view");

    jet.subtractionIdentity =
        jet.inputIdentity;

    tree.jets.nodes.push_back(
        std::move(jet));
  }


  // ==========================================================================
  // Au+Au centrality configuration
  //
  // Coarse analysis edges are retained alongside the native centile. DATA
  // centrality reconstruction can resolve either CDB-managed payloads or an
  // explicitly pinned local per-run triplet.
  // ==========================================================================

  tree.centrality.edges =
      intList(
          root["auau"]["centrality"]
              ["analysis_bin_edges_percent"]);

  if (tree.centrality.edges.empty())
  {
    tree.centrality.edges =
        {0, 20, 50, 80};
  }

  tree.centrality.reconstructDataCentrality =
      reco.reconstructCentrality;

  if (reco.reconstructCentrality)
  {
    const YAML::Node data =
        root["auau"]["centrality"]["data"];

    tree.centrality.requireOwnRunPayload =
        optional<bool>(
            data,
            "require_own_run_payload",
            true);

    tree.centrality.calibrationTag =
        required<std::string>(
            data,
            "calibration_set",
            "auau.centrality.data");

    cond.centralityCalibrationSet =
        tree.centrality.calibrationTag;

    cond.centralityManifestSha256 =
        lower(
            required<std::string>(
                data,
                "calibration_manifest_sha256",
                "auau.centrality.data"));

    cond.centralityLocalPayloads =
        lower(
            optionalString(
                data,
                "source",
                "cdb")) != "cdb";

    /*
     * Run-scoped entries override the shared defaults as one coherent payload
     * triplet. This keeps divisions, run scale, and vertex scale tied to the
     * same reviewed run binding.
     */
    YAML::Node resolved =
        data;

    if (const YAML::Node runs =
            root["auau"]["centrality"]["runs"])
    {
      if (const YAML::Node override =
              runs[std::to_string(job.run)])
      {
        resolved =
            override;
      }
    }

    cond.centralityPayloadRun =
        optional<int>(
            resolved,
            "payload_run",
            job.run);

    cond.centralityDivisions =
        expandEnvironment(
            optionalString(
                resolved["divisions"],
                "path"));

    cond.centralityDivisionsSha256 =
        lower(
            optionalString(
                resolved["divisions"],
                "sha256"));

    cond.centralityRunScale =
        expandEnvironment(
            optionalString(
                resolved["run_scale"],
                "path"));

    cond.centralityRunScaleSha256 =
        lower(
            optionalString(
                resolved["run_scale"],
                "sha256"));

    cond.centralityVertexScale =
        expandEnvironment(
            optionalString(
                resolved["vertex_scale"],
                "path"));

    cond.centralityVertexScaleSha256 =
        lower(
            optionalString(
                resolved["vertex_scale"],
                "sha256"));

    cond.centralityFallbackEvidenceSha256 =
        lower(
            optionalString(
                resolved,
                "fallback_evidence_sha256"));

    tree.centrality.divisionPayloadPath =
        cond.centralityDivisions;

    tree.centrality.runScalePayloadPath =
        cond.centralityRunScale;

    tree.centrality.vertexScalePayloadPath =
        cond.centralityVertexScale;
  }


  // ==========================================================================
  // Au+Au DATA CEMC status and zero-suppression calibration
  //
  // The status guard ties the current reconstruction to the original
  // calorimeter input. Missing historical status can be restored only when an
  // explicitly authorized, hash-pinned recovery payload is supplied.
  // ==========================================================================

  if (reco.isAuAu && reco.isData)
  {
    const YAML::Node status =
        root["auau"]["cemc_status"];

    cond.cemcStatusGuard =
        optional<bool>(
            status,
            "validate_paired_input_status",
            true);

    cond.cemcOriginalCalorimeterInput =
        expandEnvironment(
            optionalString(
                status,
                "original_calorimeter_input"));

    cond.cemcStatusRecoveryAuthorized =
        optional<bool>(
            status,
            "allow_missing_status_recovery",
            false);

    cond.cemcStatusRecoveryPayload =
        expandEnvironment(
            optionalString(
                status["recovery_payload"],
                "path"));

    cond.cemcStatusRecoveryPayloadSha256 =
        lower(
            optionalString(
                status["recovery_payload"],
                "sha256"));
  }

  cond.cemcZeroSuppressionCrossCalibration =
      reco.isData &&
      optional<bool>(
          root["calorimeter"]["cemc"],
          "zero_suppression_cross_calibration",
          true);


  // ==========================================================================
  // Archived p+p double-interaction conditions and RNG replay
  // ==========================================================================

  if (reco.isArchivedDoubleInteraction)
  {
    cond.cdbGlobalTag =
        required<std::string>(
            profile,
            "cdb_global_tag",
            context);

    cond.cdbTimestamp =
        required<std::uint64_t>(
            profile,
            "cdb_timestamp",
            context);

    const YAML::Node rng =
        profile["historical_rng_contract"];

    reco.expectedPedestalSequence =
        required<int>(
            rng,
            "expected_pedestal_sequence",
            context);

    reco.forbidRecoConstsRandomSeed =
        optional<bool>(
            rng,
            "forbid_reco_consts_random_seed",
            true);

    for (const YAML::Node& seed :
         rng["seed_sequence"])
    {
      reco.historicalSeeds.push_back(
          seed.as<unsigned int>());
    }
  }
  else
  {
    cond.cdbTimestamp =
        job.run > 0
            ? static_cast<std::uint64_t>(job.run)
            : static_cast<std::uint64_t>(
                  kArchivedDoubleInteractionRun);

    cond.cdbGlobalTag =
        optionalString(
            profile,
            "cdb_global_tag");
  }


  // ==========================================================================
  // Producer-owned simulation weights
  //
  // Retain PPG12 relative event weights; absolute sample normalization and
  // source ownership remain downstream. The period selects its own vertex fit.
  // ==========================================================================

  YAML::Node reweight;
  if (reco.isSimulation && !reco.isAuAu && (role == "photon_jet" || role == "inclusive_jet"))
  {
    const YAML::Node weights = root["ppg12_weights"];
    const YAML::Node period = weights["periods"][tree.source.period];
    const YAML::Node family = weights[role];
    if (!period || !family || !family["cross_sections_pb"][tree.source.sample])
      fail("PPG12 weights require an exact maintained sample and period binding");
    const std::string interaction = lower(tree.source.siDiRole);
    if (interaction != "si" && interaction != "di")
      fail("PPG12 weights require si_di_role=SI or DI");
    auto& w = tree.weights;
    w.ppg12Nominal = true;
    w.doubleInteraction = interaction == "di";
    if (w.doubleInteraction != (inputMode == "archived_g4_only"))
      fail("PPG12 interaction role differs from the selected reconstruction profile");
    w.ppg12LaneCode = role == "photon_jet" ? 1 : 2;
    const double sigma = family["cross_sections_pb"][tree.source.sample].as<double>();
    const double reference = required<double>(family, "reference_cross_section_pb", "ppg12_weights");
    const double luminosity = required<double>(period, "luminosity_pb_inv", "ppg12_weights.period");
    const double total = required<double>(weights, "total_luminosity_pb_inv", "ppg12_weights");
    const double fraction = required<double>(period, "double_fraction", "ppg12_weights.period");
    if (!std::isfinite(sigma) || sigma <= 0 || !std::isfinite(reference) || reference <= 0 ||
        !std::isfinite(luminosity) || luminosity <= 0 || !std::isfinite(total) || total <= 0 ||
        !std::isfinite(fraction) || fraction < 0 || fraction > 1)
      fail("invalid PPG12 weight constants");
    w.sliceFactor = sigma / reference;
    w.interactionFactor = w.doubleInteraction ? fraction : 1.0 - fraction;
    w.periodFactor = luminosity / total;
    reweight = period["vertex_reweighting"];
    if (!reweight || !required<bool>(reweight, "enabled", "ppg12_weights.vertex_reweighting"))
      fail("nominal PPG12 weights require the period truth-vertex fit");
  }
  else reweight = profile["vertex_reweighting"];
  if (reweight)
  {
    cond.applyVertexReweight =
        optional<bool>(
            reweight,
            "enabled",
            false);

    cond.vertexReweightFile =
        expandEnvironment(
            optionalString(
                reweight,
                "file"));

    cond.vertexReweightFileSha256 =
        lower(
            optionalString(
                reweight,
                "sha256"));

    cond.vertexReweightHistogram =
        optionalString(
            reweight,
            "histogram");
  }

  tree.weights.applyVertexReweight =
      cond.applyVertexReweight;

  tree.weights.vertexReweightFile =
      cond.vertexReweightFile;

  tree.weights.vertexReweightFileSha256 =
      cond.vertexReweightFileSha256;

  tree.weights.vertexReweightHistogram =
      cond.vertexReweightHistogram;


  // ==========================================================================
  // Optional high-volume output tables
  // ==========================================================================

  const YAML::Node output =
      root["output"];

  // Expand the schema shared with Python readers. Storage widths never
  // depend on the first event's collection sizes.
  const auto schemaPath = std::filesystem::path(yamlFile).parent_path() / "eventTree_v1.yaml";
  const auto schema = YAML::LoadFile(schemaPath.string());
  tree.output.schemaSha256 = fileSha256(schemaPath.string());
  for(const auto& section:{"branches","reserved"})for(const auto& field:schema[section]) {
    const auto pattern=field["name"].as<std::string>();
    std::vector<std::string> names;
    const auto pos=pattern.find("<RR>");
    if(pos==std::string::npos)names.push_back(pattern);
    else for(const auto& radius:{"02","03","04"}) { auto name=pattern;name.replace(pos,4,radius);names.push_back(name); }
    for(const auto& name:names) {
      OutputField f;f.name=name;f.type=field["type"].as<std::string>();
      f.title=field["title"].as<std::string>();
      const auto unit=field["unit"]?field["unit"].as<std::string>():"";
      if(!unit.empty()&&f.title.find('[')==std::string::npos)f.title+=" ["+unit+"]";
      if(f.type.rfind("vector<",0)==0){f.type=f.type.substr(7,f.type.size()-8);f.length=-1;}
      else if(const auto bracket=f.type.find('[');bracket!=std::string::npos){
        f.length=std::stoi(f.type.substr(bracket+1));f.type.resize(bracket);
      }
      f.counter="n_"+name;
      if(name.rfind("jet",0)==0||name.rfind("truthjet",0)==0)
        if(name.size()>3&&name.substr(name.size()-3)=="_pt")f.counter="n_"+name.substr(0,name.size()-3);
      tree.output.fields.push_back(std::move(f));
    }
  }

  const auto modelPath=schemaPath.parent_path().parent_path().parent_path()/"PhotonID/models/nominal.json";
  const auto nominal=YAML::LoadFile(modelPath.string());
  const auto featureSystem=tree.system==CollisionSystem::AuAu?"auau":"pp";
  tree.output.modelFeatures=nominal[featureSystem]["features"]["order"].as<std::vector<std::string>>();

  tree.output.writePhotonCells =
      optional<bool>(
          output,
          "write_photon_cell_table",
          true);

  tree.output.writeIsolationConstituents =
      optional<bool>(
          output,
          "write_isolation_constituent_table",
          true);

  tree.output.writePhotonJetPairs =
      optional<bool>(
          output,
          "write_photon_jet_pair_table",
          true);


  // ==========================================================================
  // Event-retention policy and provenance
  // ==========================================================================

  tree.retainEveryProducerEncounter =
      optional<bool>(
          root["capture"],
          "retain_every_producer_encounter",
          true);

  tree.requireEmbeddedMinimumBias =
      optional<bool>(
          profile,
          "require_embedded_minimum_bias",
          false);

  tree.provenance.productionTag =
      required<std::string>(
          root["contract"],
          "name",
          "contract");

  tree.provenance.softwareRelease =
      optionalString(
          root["runtime"],
          "sphenix_release");

  tree.provenance.producerGitCommit =
      optionalString(
          root["runtime"],
          "producer_git_commit");

  // Bind the loaded library and compiled sources, not a manually copied label.
  Dl_info library{};
  if (!dladdr(reinterpret_cast<const void*>(&LoadPlan), &library) || !library.dli_fname)
    fail("cannot identify the loaded PhotonJetTree library");
  auto pinnedRuntime = [&](const std::string& key, const std::string& actual) {
    if (!validSha256(actual)) fail("missing runtime digest: " + key);
    const std::string expected = optionalString(root["runtime"], key);
    if (!expected.empty() && expected != actual) fail("runtime digest differs: " + key);
    return actual;
  };
  tree.provenance.producerSourceSha256 = pinnedRuntime("producer_source_sha256", PHOTONJET_SOURCE_SHA256);
  tree.provenance.producerLibrarySha256 = pinnedRuntime("producer_library_sha256", detail::fileSha256(library.dli_fname));
  tree.provenance.macroSha256 = pinnedRuntime("macro_sha256", job.macroSha256);

  tree.provenance.configurationSha256 =
      plan.configurationSha256;

  tree.provenance.configurationText =
      plan.configurationText;

  tree.provenance.calibrationManifestSha256 =
      cond.centralityManifestSha256;


  // ==========================================================================
  // Physical inputs
  //
  // Input-list interpretation belongs to Setup.cc. At this point the
  // profile is sufficiently resolved to determine which physical streams are
  // required and how they must be synchronized.
  // ==========================================================================

  plan.inputs =
      loadInputs(
          job,
          reco);

  /*
   * Auto truth-jet mode is resolved only after the physical input inventory is
   * known. If shipped truth jets are present they are consumed; otherwise the
   * reconstruction layer must build them.
   */
  if (reco.truthJetMode == TruthJetMode::Auto)
  {
    reco.truthJetMode =
        hasInputStream(
            plan,
            "TRUTH_JETS")
            ? TruthJetMode::DST
            : TruthJetMode::Build;
  }

  Validate(plan);

  return plan;
}


// ============================================================================
// Plan validation
//
// Validate establishes the invariants that must hold before ConfigureConditions
// opens or resolves production payloads and before any Fun4All subsystem is
// registered.
//
// These checks intentionally mix structural consistency with frozen nominal
// production requirements. A Plan that cannot prove its intended semantics is
// rejected before the event loop.
// ============================================================================

void Validate(const Plan& plan)
{
  const Reconstruction& reco =
      plan.reconstruction;

  const Conditions& cond =
      plan.conditions;


  // ==========================================================================
  // Input and event-range integrity
  // ==========================================================================

  if (plan.inputs.empty())
  {
    fail(
        "no input streams were resolved");
  }

  if (plan.job.firstEntry < 0 ||
      plan.job.firstEntry > std::numeric_limits<int>::max() ||
      plan.job.nEvents < -1 ||
      plan.job.nEvents == 0 ||
      plan.job.nEvents > std::numeric_limits<int>::max())
  {
    fail(
        "invalid event range: use -1 for all events, or a positive bounded count");
  }

  /*
   * One synchronized physical source bundle is admitted per invocation.
   * Supporting several list rows in one SourceRecord would require explicit
   * source-boundary observation and accounting first.
   */
  for (const auto& stream : plan.inputs)
  {
    if (stream.files.size() != 1)
    {
      fail(
          "one source bundle per invocation is required: " +
          stream.name);
    }
  }

  if (cond.cdbGlobalTag.empty())
  {
    fail(
        "profile must bind an explicit CDB global tag");
  }

  if (!cond.applyJetEnergyScale)
  {
    fail(
        "canonical reconstructed jets require JES");
  }

  if (cond.jetEnergyScaleCdbKey != "JES_Calib_Default")
  {
    fail(
        "JES key must match the key consumed by JetCalib");
  }

  if (cond.applyVertexReweight &&
      !validSha256(
          cond.vertexReweightFileSha256))
  {
    fail(
        "vertex weight payload requires SHA256");
  }

  for (const InputStream& stream :
       plan.inputs)
  {
    if (stream.required &&
        stream.files.empty())
    {
      fail(
          "required input stream is empty: " +
          stream.name);
    }
  }


  // ==========================================================================
  // DATA / simulation role consistency
  // ==========================================================================

  if (reco.isData &&
      plan.tree.simulationRole != SimulationRole::None)
  {
    fail(
        "a data profile cannot carry a simulation role");
  }

  if (reco.isSimulation &&
      plan.tree.simulationRole == SimulationRole::None)
  {
    fail(
        "a simulation profile must name its sample role");
  }


  // ==========================================================================
  // Archived p+p double-interaction production
  //
  // This lane remains fail-closed until the historically qualified waveform /
  // tower reconstruction and its exact random-seed consumption are registered.
  // ==========================================================================

  if (reco.isArchivedDoubleInteraction)
  {
    if (!reco.isPP ||
        !reco.isSimulation)
    {
      fail(
          "archived double interaction is p+p simulation only");
    }

    /*
     * TODO(CODE): port the qualified G4-hit -> waveform/tower sequence and its
     * seed consumers before enabling this lane. CaloTowerCalib alone cannot
     * turn archived G4 hits into towers. Never manufacture an empty success.
     */
    fail(
        "archived DI waveform/tower registration is not yet bound");

    if (reco.truthJetMode != TruthJetMode::DST)
    {
      fail(
          "archived double interaction requires shipped truth jets");
    }

    if (cond.cdbTimestamp !=
        static_cast<std::uint64_t>(
            kArchivedDoubleInteractionRun))
    {
      fail(
          "archived double interaction requires the run-28 conditions timestamp");
    }
  }


  // ==========================================================================
  // Real-run requirement for DATA and embedded production
  // ==========================================================================

  if ((reco.isData || reco.isEmbedded) &&
      plan.job.run <= kDataRunThreshold)
  {
    fail(
        "data and embedded production require a real data run number");
  }


  // ==========================================================================
  // Au+Au DATA centrality contract
  //
  // Nominal centrality must use the frozen calibration set and manifest.
  // Locally supplied payloads must form one fully pinned run-specific triplet.
  //
  // Own-run payloads are the default. The only permitted cross-run fallback is
  // the explicitly qualified run-68144 package, accompanied by pinned evidence.
  // ==========================================================================

  if (reco.reconstructCentrality)
  {
    if (cond.centralityCalibrationSet !=
        kNominalCentralitySet)
    {
      fail(
          "Au+Au centrality is not bound to the nominal frozen calibration set");
    }

    if (cond.centralityManifestSha256 !=
        kNominalCentralityManifestSha256)
    {
      fail(
          "Au+Au centrality manifest SHA256 does not match the frozen nominal set");
    }

    if (cond.centralityLocalPayloads)
    {
      if (cond.centralityPayloadRun <= 0)
      {
        fail(
            "local centrality requires a positive payload run");
      }

      if (cond.centralityDivisions.empty() ||
          cond.centralityRunScale.empty() ||
          cond.centralityVertexScale.empty())
      {
        fail(
            "local centrality requires divisions, run-scale and vertex-scale payloads");
      }

      if (!validSha256(
              cond.centralityDivisionsSha256) ||
          !validSha256(
              cond.centralityRunScaleSha256) ||
          !validSha256(
              cond.centralityVertexScaleSha256))
      {
        fail(
            "local centrality requires an exact SHA256 for each of the three payloads");
      }

      if (cond.centralityPayloadRun !=
          plan.job.run)
      {
        if (cond.centralityPayloadRun !=
            kQualifiedCentralityFallbackRun)
        {
          fail(
              "only the qualified run-68144 centrality fallback is permitted");
        }

        if (!validSha256(
                cond.centralityFallbackEvidenceSha256))
        {
          fail(
              "the centrality fallback requires pinned qualification evidence");
        }
      }
      else if (!cond.centralityFallbackEvidenceSha256.empty())
      {
        fail(
            "own-run centrality cannot carry a fallback approval");
      }
    }
  }


  // ==========================================================================
  // Au+Au DATA CEMC tower-status contract
  // ==========================================================================

  if (cond.cemcStatusGuard)
  {
    if (cond.cemcOriginalCalorimeterInput.empty())
    {
      fail(
          "the CEMC status guard requires the exact original calorimeter input file");
    }

    if (cond.cemcStatusRecoveryAuthorized &&
        (cond.cemcStatusRecoveryPayload.empty() ||
         !validSha256(
             cond.cemcStatusRecoveryPayloadSha256)))
    {
      fail(
          "CEMC status recovery requires a pinned run-specific payload and its SHA256");
    }
  }


  // ==========================================================================
  // Jet-energy-scale contract
  //
  // The canonical analysis uses the legacy calibration method exactly once,
  // without electromagnetic-fraction calibration, and with the frozen common
  // p+p v6 payload identity.
  // ==========================================================================

  if (cond.applyJetEnergyScale)
  {
    if (cond.jetEnergyScaleUsesEmFraction)
    {
      fail(
          "the nominal jet energy scale is the legacy method without the electromagnetic-fraction calibration");
    }

    if (cond.jetEnergyScalePayloadSha256 !=
        kNominalJetEnergyScaleSha256)
    {
      fail(
          "jet energy-scale payload SHA256 differs from the nominal common p+p v6 payload");
    }
  }


  // ==========================================================================
  // Producer-side vertex reweighting
  // ==========================================================================

  if (cond.applyVertexReweight)
  {
    if (!reco.isSimulation)
    {
      fail(
          "vertex reweighting applies to simulation only");
    }

    if (cond.vertexReweightFile.empty() ||
        cond.vertexReweightHistogram.empty())
    {
      fail(
          "vertex reweighting requires a file and a histogram name");
    }
  }


  // ==========================================================================
  // Au+Au flow-modulation configuration
  // ==========================================================================

  if (reco.isAuAu &&
      (reco.flowModulation < 0 ||
       reco.flowModulation > 3))
  {
    fail(
        "Au+Au flow modulation must be resolved to the accepted production value");
  }
}


// ============================================================================
// Human-readable resolved-plan report
//
// Print is intentionally concise. Its purpose is to make the major resolved
// production decisions visible at job startup without requiring a user to
// inspect YAML, manifests, or the complete embedded configuration.
//
// The full machine-readable provenance remains in the output product.
// ============================================================================

void Print(const Plan& plan)
{
  const Reconstruction& reco =
      plan.reconstruction;

  const Conditions& cond =
      plan.conditions;

  const char* truthMode =
      "none";

  switch (reco.truthJetMode)
  {
    case TruthJetMode::None:
      truthMode = "none";
      break;

    case TruthJetMode::Auto:
      truthMode = "auto";
      break;

    case TruthJetMode::DST:
      truthMode = "dst";
      break;

    case TruthJetMode::Build:
      truthMode = "build";
      break;
  }

  std::cout
      << "\n"
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
    std::cout
        << " centrality           : "
        << cond.centralityCalibrationSet
        << (cond.centralityLocalPayloads
                ? " local payload run " +
                      std::to_string(
                          cond.centralityPayloadRun)
                : " from CDB")
        << "\n";
  }

  if (cond.cemcStatusGuard)
  {
    std::cout
        << " CEMC status guard    : "
        << cond.cemcOriginalCalorimeterInput
        << (cond.cemcStatusRecoveryAuthorized
                ? " (recovery authorised)"
                : "")
        << "\n";
  }

  if (cond.applyVertexReweight)
  {
    std::cout
        << " vertex reweight      : "
        << cond.vertexReweightFile
        << ":"
        << cond.vertexReweightHistogram
        << "\n";
  }

  std::cout
      << " input streams        : "
      << plan.inputs.size()
      << "\n";

  for (const InputStream& stream :
       plan.inputs)
  {
    std::cout
        << "   "
        << std::left
        << std::setw(12)
        << stream.name
        << " files="
        << stream.files.size()
        << " sync="
        << (stream.synchronized ? "yes" : "no")
        << "\n";
  }

  std::cout
      << "============================================================\n"
      << std::endl;
}

}  // namespace production
}  // namespace photonjet
