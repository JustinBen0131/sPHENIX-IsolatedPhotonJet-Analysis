/**
 * @file Production.cc
 * Resolve a source Job and reconstruction preset into a validated Plan.
 * Configuration values bind algorithms and payloads without duplicating them.
 * Missing conditions, incompatible modes and identity ambiguity fail closed.
 * Setup.cc performs the conditions and Fun4All registration steps.
 */

#include "Production.h"
#include "BuildIdentity.h"

#include "../EventTree.h"

#include <TMVA/RBDT.hxx>
#include <nlohmann/json.hpp>
#include "Configuration.h"
#include <dlfcn.h>

#include <algorithm>
#include <cctype>
#include <cmath>
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

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error("PhotonJetTree production: " + message); }


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
using photonjet::production::detail::ConfigurationNode;


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

constexpr const char* kNominalCentralitySet = "sam_fit_20260913_v1";

constexpr const char* kNominalCentralityManifestSha256 =
    "4a11e0e10ae9b82062325c3275679935da2259f4c91f796950a6b5ddd2e1501c";

constexpr const char* kNominalJetEnergyScaleSha256 = "76a8788fdb4e0b3859d01361f884e295ea609bb8c1cc188105389214db2de27d";

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


// ============================================================================
// Parsed-configuration access helpers
//
// Required values fail immediately when absent or null. Optional values use an
// explicit caller-provided fallback. These helpers keep configuration semantics
// uniform throughout LoadPlan.
// ============================================================================

template <class T>
T required(const ConfigurationNode& node, const std::string& key, const std::string& context)
{
  if (!node || !node[key] || node[key].IsNull())
  {
    fail("missing required configuration key '" + key + "' in " + context);
  }

  return node[key].as<T>();
}


template <class T>
T optional(const ConfigurationNode& node, const std::string& key, const T& fallback)
{
  if (!node || !node[key] || node[key].IsNull())
  {
    return fallback;
  }

  return node[key].as<T>();
}


std::string optionalString(const ConfigurationNode& node, const std::string& key, const std::string& fallback = {})
{
  if (!node || !node[key] || node[key].IsNull())
  {
    return fallback;
  }

  return node[key].as<std::string>();
}


std::vector<std::string> stringList(const ConfigurationNode& node)
{
  std::vector<std::string> result;

  if (node && node.IsSequence())
  {
    for (const ConfigurationNode& value : node)
    {
      result.push_back(value.as<std::string>());
    }
  }

  return result;
}


std::vector<double> doubleList(const ConfigurationNode& node)
{
  std::vector<double> result;

  if (node && node.IsSequence())
  {
    for (const ConfigurationNode& value : node)
    {
      result.push_back(value.as<double>());
    }
  }

  return result;
}


std::vector<int> intList(const ConfigurationNode& node)
{
  std::vector<int> result;

  if (node && node.IsSequence())
  {
    for (const ConfigurationNode& value : node)
    {
      result.push_back(value.as<int>());
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

  if (view == "auau_sub1" || view == "sub1")
  {
    return photonjet::JetView::AuAuSub1;
  }

  if (view == "auau_nosub" || view == "nosub")
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
// Nominal photon labels
//
// The frozen PhotonID package pins both systems' models, feature order,
// domains and working points; its nominal manifests and model digests are
// verified here once per job, exactly as the package defines them. The Au+Au
// isolation class comes from the freeze pinned in the preset; pp uses the
// PPG12 R = 0.4 definition on the stored witness.
// ============================================================================
namespace
{
nlohmann::json readJson(const std::filesystem::path& path)
{
  try
  {
    return nlohmann::json::parse(readTextFile(path.string()));
  }
  catch (const nlohmann::json::exception& error)
  {
    fail("invalid JSON in " + path.string() + ": " + error.what());
  }
}

double numberOrOpen(const nlohmann::json& value)
{
  return value.is_null() ? std::numeric_limits<double>::quiet_NaN() : value.get<double>();
}

// packageFolder is PhotonID/models; the Au+Au freeze and its digest come from
// the preset (auau.isolation.nominal) and are unused for pp.
void loadNominalDefinitions(Config& tree, const bool auau, const std::filesystem::path& packageFolder,
                            const std::string& freezeFile, const std::string& freezeSha256)
{
  NominalPhotonIdConfig& id = tree.photon.nominalId;
  id.system = auau ? "auau" : "pp";
  id.packageFile = (packageFolder / "nominal.json").string();
  id.packageSha256 = fileSha256(id.packageFile);
  const nlohmann::json package = readJson(id.packageFile);
  if (package.value("schema", "") != "CanonicalPhotonIDFreezeV1")
    fail("not a canonical photon-ID freeze: " + id.packageFile);
  for (const auto& [key, folder] : std::vector<std::pair<std::string, std::string>>{{"auau", "AuAu"}, {"pp", "PP"}})
  {
    const auto nominal = packageFolder / folder / "nominal";
    const nlohmann::json manifest = readJson(nominal / "manifest.json");
    const nlohmann::json& definition = package.at(key);
    const std::string modelSha256 = definition.at("model").at("sha256").get<std::string>();
    if (manifest.value("schema", "") != "PhotonIDFrozenNominalV1" ||
        manifest.value("freeze_sha256", "") != id.packageSha256 ||
        manifest.value("definition", nlohmann::json()) != definition || manifest.value("system", "") != folder ||
        manifest.value("artifacts", nlohmann::json::object()).value("model.tmva.root", "") != modelSha256)
      fail(key + " nominal manifest differs from the frozen definition");
    const std::string model = (nominal / "model.tmva.root").string();
    if (fileSha256(model) != modelSha256) fail(key + " model file differs from the freeze");
    if (key == id.system)
    {
      id.modelFile = model;
      id.modelSha256 = modelSha256;
    }
  }
  const nlohmann::json& definition = package.at(id.system);
  const nlohmann::json& features = definition.at("features");
  id.view = features.at("view").get<std::string>();
  const ShowerViewDefinition& view = showerViewDefinition(id.view);
  if (features.at("view_semantic_sha256").get<std::string>() != view.semanticSha256 ||
      features.at("floor_gev").get<double>() != view.floorGeV)
    fail("the nominal model view differs from the " + id.view + " shower-view definition");
  id.featureOrder = features.at("order").get<std::vector<std::string>>();
  if (id.featureOrder.empty() ||
      (features.contains("count") && features.at("count").get<std::size_t>() != id.featureOrder.size()))
    fail("the nominal feature order is empty or its count differs");
  const nlohmann::json& domain = definition.at("domain");
  id.etMin = domain.at("et_min").get<double>();
  id.etMaxExclusive = numberOrOpen(domain.at("et_max_exclusive"));
  id.usesProducerScore = !auau;
  if (auau)
  {
    id.centralityMin = domain.at("centrality_min").get<double>();
    id.centralityMaxExclusive = domain.at("centrality_max_exclusive").get<double>();
    for (const auto& [name, line] :
         std::vector<std::pair<std::string, IsolationLine*>>{{"WP70", &id.wp70}, {"WP80", &id.wp80}})
    {
      line->interceptGeV = definition.at("working_points").at(name).at("intercept").get<double>();
      line->slope = definition.at("working_points").at(name).at("slope").get<double>();
    }
  }
  else
  {
    // The pp nominal score is the producer score of this same model.
    const ProductionScoreConfig& production = tree.photon.productionScore;
    if (production.scoreName.empty() || production.modelSha256 != id.modelSha256 ||
        production.showerDefinition != id.view)
      fail("the pp nominal score requires the producer score of the nominal PPG12 model on the " + id.view + " view");
    for (const PhotonModelConfig& model : tree.photon.models)
      if (model.scoreName == production.scoreName && static_cast<double>(model.minEt) != id.etMin)
        fail("the pp producer score minimum differs from the nominal domain");
  }

  NominalIsolationConfig& iso = tree.photon.nominalIsolation;
  iso.system = id.system;
  if (auau)
  {
    iso.freezeFile = freezeFile;
    iso.freezeSha256 = freezeSha256;
    if (!validSha256(iso.freezeSha256) || fileSha256(iso.freezeFile) != iso.freezeSha256)
      fail("the Au+Au isolation freeze differs from its pinned digest: " + freezeFile);
    const nlohmann::json freeze = readJson(iso.freezeFile);
    if (freeze.value("schema", "") != "CanonicalIsolationFreezeV1" || freeze.value("system", "") != "auau" ||
        freeze.at("witness").value("reconstructed", false) != true)
      fail("not a reconstructed-witness Au+Au isolation freeze: " + freezeFile);
    iso.witnessIndex = freeze.at("witness").at("stride_index").get<int>();
    iso.radius = freeze.at("witness").at("radius").get<double>();
    iso.subtractionMethod = freeze.at("witness").at("subtraction_method").get<std::int32_t>();
    iso.threshold.interceptGeV = freeze.at("threshold").at("intercept_gev").get<double>();
    iso.threshold.slope = freeze.at("threshold").at("slope_gev_per_percent").get<double>();
    iso.threshold.sidebandGapGeV = freeze.at("sideband_gap_gev").get<double>();
    iso.centralityMin = freeze.at("domain").at("centrality_min").get<double>();
    iso.centralityMaxExclusive = freeze.at("domain").at("centrality_max_exclusive").get<double>();
    iso.derivationSha256 = freeze.at("derivation").at("package_sha256").get<std::string>();
  }
  else
  {
    const nlohmann::json ppg12 = {{"schema", "PPG12IsolationDefinitionV1"},
                                  {"system", "pp"},
                                  {"radius", 0.4},
                                  {"stride_index", 1},
                                  {"subtraction_method", 1},
                                  {"threshold", "stored reconstructed witness"},
                                  {"isolated", "cone < threshold"},
                                  {"nonisolated", "cone > sideband_threshold"},
                                  {"boundaries", "gap"},
                                  {"sim_correction", "already in stored witness; do not repeat"}};
    iso.definitionJson = ppg12.dump(-1, ' ', true);
    iso.freezeSha256 = eventtree::textSha256(iso.definitionJson);
    iso.derivationSha256 = iso.freezeSha256;
    iso.witnessIndex = 1;
    iso.radius = 0.4;
    iso.subtractionMethod = 1;
    iso.storedThresholds = true;
  }
  const auto& radii = tree.photon.isolationRadii;
  if (iso.witnessIndex < 0 || static_cast<std::size_t>(iso.witnessIndex) >= radii.size() ||
      radii[iso.witnessIndex] != iso.radius)
    fail("the nominal isolation witness radius differs from the retained isolation radii");
}
}  // namespace

PinnedModel::PinnedModel(const std::string& file, const std::string& sha256)
{
  if (fileSha256(file) != sha256) fail("model file changed after the plan was resolved: " + file);
  m_model = std::make_unique<TMVA::Experimental::RBDT>("myBDT", file);
}

PinnedModel::~PinnedModel() = default;

float PinnedModel::score(const std::vector<float>& features) const
{
  const auto prediction = m_model->Compute(features);
  if (prediction.size() != 1) fail("a model returned " + std::to_string(prediction.size()) + " outputs");
  return prediction[0];
}


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
// LoadPlan resolves values and the identities of local payloads;
// ConfigureConditions verifies those payloads and binds the run's conditions.
// ============================================================================

namespace
{
// An archived double-interaction source's random-number record,
// "s1,s2,s3,s4,s5:p": the five unsigned 32-bit seeds its production consumed
// and the pedestal file index in [0, 3260) they chose.
void parseHistoricalRng(const std::string& record, std::vector<unsigned int>& seeds, int& pedestal)
{
  const auto number = [&](const std::string& text, const unsigned long long maximum)
  {
    if (text.empty() || text.size() > 10 ||
        !std::all_of(text.begin(), text.end(), [](unsigned char c) { return std::isdigit(c) != 0; }))
      fail("malformed historical_rng record: " + record);
    const unsigned long long value = std::stoull(text);
    if (value > maximum) fail("historical_rng value out of range: " + record);
    return value;
  };
  const auto colon = record.find(':');
  if (colon == std::string::npos) fail("malformed historical_rng record: " + record);
  std::stringstream list(record.substr(0, colon));
  std::string seed;
  while (std::getline(list, seed, ',')) seeds.push_back(static_cast<unsigned int>(number(seed, 4294967295ULL)));
  pedestal = static_cast<int>(number(record.substr(colon + 1), 3259ULL));
  if (seeds.size() != 5) fail("historical_rng requires exactly five seeds: " + record);
}

// The data run an embedded DST was embedded in: <...>-data-<run>-<...>; 0
// when its name carries none.
std::uint64_t embeddedDataRun(const std::string& file)
{
  const std::string name = std::filesystem::path(file).filename().string();
  const std::size_t tag = name.find("-data-");
  if (tag == std::string::npos) return 0;
  const std::size_t begin = tag + 6, end = name.find_first_of("-.", begin);
  const std::string digits = name.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
  if (digits.empty() || digits.size() > 9 ||
      !std::all_of(digits.begin(), digits.end(), [](unsigned char c) { return std::isdigit(c) != 0; }))
    return 0;
  return std::stoull(digits);
}

// The run's centrality payload triplet from the frozen calibration set. The
// set's MANIFEST.json, pinned by SHA256, lists every payload with its digest.
// A run without its own triplet may use only the qualified run-68144 triplet,
// and only when the configuration pins that fallback's evidence.
void resolveCentralityPayloads(Conditions& cond, const ConfigurationNode& data, const int run)
{
  const std::filesystem::path root =
      expandEnvironment(required<std::string>(data, "calibration_root", "auau.centrality.data"));
  const std::filesystem::path manifestFile = root / "MANIFEST.json";
  if (fileSha256(manifestFile.string()) != cond.centralityManifestSha256)
    fail("the centrality calibration manifest differs from its pinned digest: " + manifestFile.string());
  const nlohmann::json manifest = readJson(manifestFile);
  if (manifest.value("schema", "") != "FrozenCentralityCalibrationV1")
    fail("not a frozen centrality calibration manifest: " + manifestFile.string());
  const nlohmann::json& files = manifest.at("files");

  const auto triplet = [](const int payloadRun)
  {
    const std::string r = std::to_string(payloadRun);
    return std::vector<std::string>{"divs/cdb_centrality_" + r + ".root", "scales/cdb_centrality_scale_" + r + ".root",
                                    "vertexscales/cdb_centrality_vertex_scale_" + r + ".root"};
  };
  const auto listed = [&](const int payloadRun)
  {
    const auto names = triplet(payloadRun);
    return std::all_of(names.begin(), names.end(), [&](const std::string& name) { return files.contains(name); });
  };

  cond.centralityFallbackEvidenceSha256 = lower(optionalString(data, "fallback_evidence_sha256"));
  cond.centralityPayloadRun = listed(run) ? run : kQualifiedCentralityFallbackRun;
  if (cond.centralityPayloadRun != run && cond.centralityFallbackEvidenceSha256.empty())
    fail("run " + std::to_string(run) + " has no centrality payloads in " + cond.centralityCalibrationSet);
  if (!listed(cond.centralityPayloadRun))
    fail("the frozen centrality set lacks the run-" + std::to_string(cond.centralityPayloadRun) + " triplet");

  const auto names = triplet(cond.centralityPayloadRun);
  const std::vector<std::pair<std::string*, std::string*>> targets = {
      {&cond.centralityDivisions, &cond.centralityDivisionsSha256},
      {&cond.centralityRunScale, &cond.centralityRunScaleSha256},
      {&cond.centralityVertexScale, &cond.centralityVertexScaleSha256}};
  for (std::size_t i = 0; i < names.size(); ++i)
  {
    *targets[i].first = (root / names[i]).string();
    *targets[i].second = lower(files.at(names[i]).get<std::string>());
  }
}
}  // namespace

Plan LoadPlan(const std::string& yamlFile, const Job& job)
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
  plan.configurationText = readTextFile(yamlFile);

  plan.configurationSha256 = fileSha256(yamlFile);

  nlohmann::json parsedConfiguration;
  if (job.resolvedConfigurationFile.empty())
  {
    // Direct native callers may provide a JSON preset. The ordinary launcher
    // supplies the envelope for the distributed YAML preset.
    parsedConfiguration = readJson(yamlFile);
  }
  else
  {
    const nlohmann::json envelope = readJson(job.resolvedConfigurationFile);
    if (envelope.value("schema", "") != "TreeProductionConfigurationV1" ||
        envelope.value("source_sha256", "") != plan.configurationSha256 || !envelope.contains("document") ||
        !envelope.at("document").is_object())
      fail("parsed configuration is not bound to the original preset");
    parsedConfiguration = envelope.at("document");
    plan.tree.provenance.resolvedConfigurationSha256 = fileSha256(job.resolvedConfigurationFile);
  }
  const ConfigurationNode root(parsedConfiguration);

  const ConfigurationNode profile = root["profiles"][job.profile];

  if (!profile)
  {
    fail("profile '" + job.profile + "' is absent from " + yamlFile);
  }

  const std::string context = "profiles." + job.profile;

  Reconstruction& reco = plan.reconstruction;

  Conditions& cond = plan.conditions;

  photonjet::Config& tree = plan.tree;


  // ==========================================================================
  // Production mode
  //
  // Resolve the four fundamental mode dimensions first. All later stage
  // registration is derived from these flags rather than from repeated profile
  // string tests.
  // ==========================================================================

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

  if (reco.isPP == reco.isAuAu)
  {
    fail("profile must name exactly one collision system");
  }

  if (reco.isData == reco.isSimulation)
  {
    fail("profile must name exactly one data kind");
  }

  tree.system = reco.isPP ? CollisionSystem::PP : CollisionSystem::AuAu;

  tree.dataKind = reco.isData ? DataKind::Data : DataKind::Simulation;

  tree.inputMode = reco.isArchivedDoubleInteraction ? InputMode::ArchivedG4Only : InputMode::ReconstructedDST;

  if (role == "photon_jet")
  {
    tree.simulationRole = SimulationRole::PhotonJet;
  }
  else if (role == "inclusive_jet")
  {
    tree.simulationRole = SimulationRole::InclusiveJet;
  }
  else if (role == "minimum_bias")
  {
    tree.simulationRole = SimulationRole::MinimumBias;
  }
  else
  {
    tree.simulationRole = SimulationRole::None;
  }

  tree.outputFile = job.outputFile;

  if (inputMode != "archived_g4_only" && inputMode != "reconstructed_dst")
  {
    fail("unknown input mode: " + inputMode);
  }

  if (role != "none" && role != "photon_jet" && role != "inclusive_jet" && role != "minimum_bias")
  {
    fail("unknown simulation role: " + role);
  }


  // ==========================================================================
  // Source binding
  //
  // Per-job source coordinates come from produce.py and remain distinct from
  // shared reconstruction configuration. The profile supplies only a default
  // dataset token when the Job did not explicitly provide one.
  // ==========================================================================

  tree.source.lane = job.profile;

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

  tree.source.sourceGlobalEntryBegin = job.sourceGlobalEntryBegin;

  tree.source.firstEntry = job.firstEntry;

  tree.source.maxEvents = job.nEvents;


  // ==========================================================================
  // Canonical input-node vocabulary
  //
  // These names tell PhotonJetTree where to read the reconstructed event.
  // Setup.cc is responsible for creating the required nodes
  // when reconstruction occurs in this pass.
  // ==========================================================================

  const ConfigurationNode nodes = root["nodes"];

  tree.nodes.eventHeader = optionalString(nodes, "event_header", "EventHeader");

  tree.nodes.gl1Packet = optionalString(nodes, "gl1_packet", "GL1Packet");

  tree.nodes.globalVertexMap = optionalString(nodes, "global_vertex_map", "GlobalVertexMap");

  tree.nodes.centralityInfo = optionalString(nodes, "centrality_info", "CentralityInfo");

  tree.nodes.minimumBiasInfo = optionalString(nodes, "minimum_bias_info", "MinimumBiasInfo");

  tree.nodes.mbdOut = optionalString(nodes, "mbd_out", "MbdOut");

  tree.nodes.mbdPmts = optionalString(nodes, "mbd_pmts", "MbdPmtContainer");

  tree.nodes.mbdPmtsFallback = optionalString(nodes, "mbd_pmts_fallback", "MbdPmtContainer_data");

  tree.nodes.cemcTowers = optionalString(nodes, "cemc_towers", "TOWERINFO_CALIB_CEMC");

  tree.nodes.ihcalTowers = optionalString(nodes, "ihcal_towers", "TOWERINFO_CALIB_HCALIN");

  tree.nodes.ohcalTowers = optionalString(nodes, "ohcal_towers", "TOWERINFO_CALIB_HCALOUT");

  tree.nodes.cemcGeometry = optionalString(nodes, "cemc_geometry", "TOWERGEOM_CEMC");

  tree.nodes.truthInfo = optionalString(nodes, "truth_info", "G4TruthInfo");

  tree.nodes.hepmcEventMap = optionalString(nodes, "hepmc_event_map", "PHHepMCGenEventMap");


  // ==========================================================================
  // Calorimeter event-state capture
  // ==========================================================================

  tree.calorimeter.requireGoodTowers = optional<bool>(root["calorimeter"]["event_energy"], "require_good_towers", true);


  // ==========================================================================
  // Photon capture configuration
  //
  // These values define the broad model-independent storage domain and the
  // reconstruction observables retained for later PhotonID/analysis stages.
  // They are not the final measurement selection.
  // ==========================================================================

  const ConfigurationNode photons = root["photons"];

  const ConfigurationNode domain = photons["capture_domain"];

  tree.photon.minEtGeV = required<double>(domain, "min_et_gev", "photons.capture_domain");

  tree.photon.maxEtGeV = required<double>(domain, "max_et_gev", "photons.capture_domain");

  tree.photon.maxAbsEta = required<double>(domain, "max_abs_eta", "photons.capture_domain");

  // A lane may raise the lower edge: Au+Au DATA captures photons from 12 GeV.
  tree.photon.minEtGeV = optional<double>(profile, "photon_min_et_gev", tree.photon.minEtGeV);

  if (!(tree.photon.minEtGeV >= 0.0 && tree.photon.minEtGeV < tree.photon.maxEtGeV))
  {
    fail(context + ".photon_min_et_gev must lie inside the photon capture domain");
  }

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

  tree.photon.truthIsolationSignalRadius = optional<double>(photons["isolation"]["truth"], "signal_radius", 0.30);

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

  // Optional exact named models. No implicit file, feature order or physics
  // domain is inferred from a score name; generic production defaults stay off.
  if (photons["models"])
    for (const auto& node : photons["models"])
    {
      PhotonModelConfig model;
      model.scoreName = required<std::string>(node, "score_name", "photons.models");
      model.file = required<std::string>(node, "file", "photons.models");
      model.sha256 = required<std::string>(node, "sha256", "photons.models");
      model.features = stringList(node["features"]);
      if (model.features.empty()) fail("photons.models requires an explicit feature order");
      // Each invocation bound is optional; an absent bound leaves that side open.
      model.minEt = optional<float>(node, "min_et", NAN);
      model.maxEt = optional<float>(node, "max_et", NAN);
      model.maxAbsEta = optional<float>(node, "max_abs_eta", NAN);
      model.unevaluatedScore = optional<float>(node, "unevaluated_score", NAN);
      model.configuredCut = optional<double>(node, "configured_cut", NAN);
      tree.photon.models.push_back(model);
    }

  // Active photon preselection (pho_presel_active): the profile's established
  // variant, reference or PPG12's newPPG12 on a captured NPB score with a
  // configured cut.
  if (const ConfigurationNode pre = photons["preselection"])
  {
    const std::string where = "photons.preselection";
    PreselectionConfig& preselection = tree.photon.preselection;
    preselection.variant = required<std::string>(profile, "photon_preselection", "profiles." + job.profile);
    if (preselection.variant != "reference" && preselection.variant != "newPPG12")
      fail("unsupported profiles." + job.profile + ".photon_preselection: " + preselection.variant);
    preselection.npbScoreName = required<std::string>(pre, "npb_score", where);
    if (preselection.npbScoreName != "npb_score" && preselection.npbScoreName != "auau_npb_score")
      fail(where + ".npb_score must name a captured NPB score");
    const auto window = [&](const char* key, double& low, double& high)
    {
      const ConfigurationNode range = pre[key];
      if (!range || !range.IsSequence() || range.size() != 2) fail(where + "." + key + " must be [low, high]");
      low = range[0].as<double>();
      high = range[1].as<double>();
      if (!(low < high)) fail("empty " + where + "." + key);
    };
    window("et1", preselection.et1Min, preselection.et1Max);
    window("e32_over_e35", preselection.e32OverE35Min, preselection.e32OverE35Max);
    preselection.e11OverE33Max = required<double>(pre, "e11_over_e33_max", where);
    const auto model =
        std::find_if(tree.photon.models.begin(), tree.photon.models.end(),
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

  reco.inputLayout = reco.isArchivedDoubleInteraction ? InputLayout::ArchivedDoubleInteraction
                     : reco.isData                    ? InputLayout::PairedData
                                                      : InputLayout::ReconstructedSimulation;

  reco.truthJetMode =
      parseTruthJetMode(optional<std::string>(profile, "truth_jets_mode", reco.isSimulation ? "auto" : "none"));

  reco.reconstructEventGeometry = !reco.isEmbedded;

  reco.reconstructZdc = reco.isData;

  reco.reconstructCalorimeter = reco.isData || reco.isEmbedded || reco.isArchivedDoubleInteraction;

  reco.runCalorimeterStatusSkimmer = reco.isData;

  // Simulation sets tower status from the conditions database. Paired DATA
  // keeps the status its DSTs were produced with (the Au+Au CEMC status guard
  // validates it), and embedded simulation keeps its producer status.
  reco.setTowerStatus = reco.reconstructCalorimeter && !reco.isData &&
                        (!reco.isEmbedded || optional<bool>(profile, "force_embedded_tower_status", false));

  reco.reconstructCentrality = reco.isAuAu && reco.isData && optional<bool>(profile, "reconstruct_centrality", true);

  reco.runMinimumBiasClassifier =
      reco.reconstructCentrality && optional<bool>(profile, "require_minimum_bias_classifier", true);

  reco.reconstructAuAuBackground = reco.isAuAu;

  reco.reconstructJets = true;

  reco.reconstructPhotons = true;


  // ==========================================================================
  // Calorimeter and photon reconstruction parameters
  // ==========================================================================

  const ConfigurationNode clusterBuilder = root["calorimeter"]["cluster_builder"];

  reco.clusterThresholdGeV = required<double>(clusterBuilder, "threshold_energy_gev", "calorimeter.cluster_builder");

  reco.clusterProfile = expandEnvironment(
      required<std::string>(clusterBuilder["profile"], "path", "calorimeter.cluster_builder.profile"));

  reco.simulationHotTowerMapKey =
      optionalString(root["calorimeter"]["tower_status"], "simulation_hot_map_key", "CEMC_hotTowers_status");

  const ConfigurationNode systemConfig = reco.isPP ? root["pp"] : root["auau"];

  reco.photonMinEtGeV = tree.photon.minEtGeV;

  reco.photonUseVertexCut = optional<bool>(photons["builder"], "use_vertex_cut", false);

  reco.photonVertexCutCm = optional<double>(photons["builder"], "vertex_cut_cm", 60.0);

  reco.showerShapeTowerMinEnergyGeV =
      optional<double>(systemConfig, "shower_shape_tower_min_energy_gev", reco.isPP ? 0.070 : 0.0);

  const ConfigurationNode isolation = systemConfig["isolation"];

  const std::string isolationMethod = lower(optionalString(isolation, "method", reco.isPP ? "topocluster" : "sub1"));

  reco.useTopoclusterIsolation = reco.isPP && isolationMethod == "topocluster";

  reco.useSubtractedIsolation = reco.isAuAu;

  reco.topoclusterNode = tree.photon.topoclusterNode;

  // --------------------------------------------------------------------------
  // Isolation witnesses
  //
  // The historical working point recorded with each retained cone. pp lines
  // are in photon ET; Au+Au lines are in the truncated integer centrality
  // inside centrality_range, with the photon-ET fallback outside it.
  // --------------------------------------------------------------------------
  {
    const std::string where = std::string(reco.isPP ? "pp" : "auau") + ".isolation.witness";
    const ConfigurationNode witness = isolation["witness"];
    if (!witness) fail("missing required configuration block " + where);
    auto line = [&](const std::string& key)
    {
      const std::string context = where + "." + key;
      IsolationLine value;
      value.interceptGeV = required<double>(witness[key], "intercept_gev", context);
      value.slope = required<double>(witness[key], "slope", context);
      value.sidebandGapGeV = required<double>(witness[key], "sideband_gap_gev", context);
      if (!std::isfinite(value.interceptGeV) || !std::isfinite(value.slope) || !std::isfinite(value.sidebandGapGeV) ||
          value.sidebandGapGeV < 0.0)
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
      const ConfigurationNode range = witness["centrality_range"];
      if (!range || !range.IsSequence() || range.size() != 2)
        fail(where + ".centrality_range must be [minimum, maximum)");
      config.centralityMin = range[0].as<double>();
      config.centralityMax = range[1].as<double>();
      if (!(config.centralityMin < config.centralityMax)) fail("empty " + where + ".centrality_range");
    }
  }

  // --------------------------------------------------------------------------
  // PPG12 photon-yield source arbitration (pp only)
  //
  // PPG12's tower mask, pinned by digest and verified before reconstruction.
  // Objects.cc records each candidate's mask state and the event's leading
  // unmasked source cluster.
  // --------------------------------------------------------------------------
  if (reco.isPP)
  {
    const std::string where = "pp.ppg12_tower_mask";
    const ConfigurationNode mask = systemConfig["ppg12_tower_mask"];
    if (!mask) fail("missing required configuration block " + where);
    tree.photon.ppg12TowerMaskFile = expandEnvironment(required<std::string>(mask, "file", where));
    tree.photon.ppg12TowerMaskSha256 = lower(required<std::string>(mask, "sha256", where));
    tree.photon.ppg12TowerMaskHistogram = required<std::string>(mask, "histogram", where);
    if (!validSha256(tree.photon.ppg12TowerMaskSha256)) fail(where + " requires a SHA256");
  }
  else if (systemConfig["ppg12_tower_mask"])
  {
    fail("the PPG12 tower mask applies to pp only");
  }

  // --------------------------------------------------------------------------
  // Producer photon-identification score
  //
  // Bound to a canonical PhotonID model package: its file, digest and input
  // order come from the package manifest, never from this preset. The preset
  // supplies only the historical invocation (score name, minimum ET,
  // unevaluated sentinel) and the recorded working-point witnesses.
  // --------------------------------------------------------------------------
  const ConfigurationNode productionScore = systemConfig["photon_id"]["production_score"];
  if (productionScore && reco.isAuAu)
  {
    // The established Au+Au production recorded its tight model's identity,
    // validated domain and centrality-linear working points with every
    // candidate, without evaluating the score (pho_prodbdt_score is NaN).
    const std::string where = "auau.photon_id.production_score";
    ProductionScoreConfig& production = tree.photon.productionScore;
    production.scoreName = required<std::string>(productionScore, "score_name", where);
    production.modelSha256 = required<std::string>(productionScore, "model_sha256", where);
    production.showerDefinition =
        showerViewDefinition(required<std::string>(productionScore, "shower_definition", where)).name;
    const ConfigurationNode validated = productionScore["validated_et"];
    if (!validSha256(production.modelSha256) || !validated || !validated.IsSequence() || validated.size() != 2)
      fail(where + " needs a model SHA256 and validated_et [minimum, maximum)");
    production.validatedEtMin = validated[0].as<double>();
    production.validatedEtMax = validated[1].as<double>();
    production.validatedCentralityMin = tree.photon.isolationWitness.centralityMin;
    production.validatedCentralityMax = tree.photon.isolationWitness.centralityMax;
    production.centralityLinearWorkingPoints = true;
    const std::array<const char*, 3> names{{"wp70", "wp80", "wp90"}};
    for (std::size_t index = 0; index < names.size(); ++index)
    {
      const ConfigurationNode line = productionScore["working_points"][names[index]];
      production.workingPointIntercepts[index] = required<double>(line, "intercept", where + ".working_points");
      production.workingPointSlopes[index] = required<double>(line, "slope", where + ".working_points");
      if (!std::isfinite(production.workingPointIntercepts[index]) ||
          !std::isfinite(production.workingPointSlopes[index]))
        fail(where + ".working_points must be finite lines");
    }
  }
  else if (const ConfigurationNode score = productionScore)
  {
    const std::string where = "pp.photon_id.production_score";
    const auto folder = std::filesystem::path(yamlFile).parent_path().parent_path().parent_path() / "PhotonID" /
                        "models" / (reco.isPP ? "PP" : "AuAu") / required<std::string>(score, "model", where);
    const ConfigurationNode manifest(readJson(folder / "manifest.json"));
    const ConfigurationNode definition = manifest["definition"];
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
        required<std::string>(definition["features"], "view_semantic_sha256", where + " manifest") !=
            primary.semanticSha256)
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
    const ConfigurationNode validated = score["validated_et"];
    if (!validated || !validated.IsSequence() || validated.size() != 2)
      fail(where + ".validated_et must be [minimum, maximum)");
    production.validatedEtMin = validated[0].as<double>();
    production.validatedEtMax = validated[1].as<double>();
    if (const ConfigurationNode points = score["working_points"])
    {
      for (const ConfigurationNode& edge : points["et_edges"])
        production.workingPointEtEdges.push_back(edge.as<double>());
      const std::array<const char*, 3> names{{"wp70", "wp80", "wp90"}};
      for (std::size_t index = 0; index < names.size(); ++index)
        for (const ConfigurationNode& value : points[names[index]])
          production.workingPoints[index].push_back(value.as<double>());
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

  {
    const auto configDirectory = std::filesystem::path(yamlFile).parent_path();
    std::string freezeFile, freezeSha256;
    if (reco.isAuAu)
    {
      const std::string where = "auau.isolation.nominal";
      const ConfigurationNode nominal = isolation["nominal"];
      if (!nominal) fail("missing required configuration block " + where);
      freezeFile = (configDirectory / required<std::string>(nominal, "freeze", where)).lexically_normal().string();
      freezeSha256 = required<std::string>(nominal, "freeze_sha256", where);
    }
    loadNominalDefinitions(tree, reco.isAuAu, configDirectory.parent_path().parent_path() / "PhotonID" / "models",
                           freezeFile, freezeSha256);
  }


  // ==========================================================================
  // Au+Au underlying-event reconstruction
  // ==========================================================================

  if (reco.isAuAu)
  {
    const ConfigurationNode background = root["auau"]["background_subtraction"];

    reco.retowerFractionCut = optional<double>(background["retower_cemc"], "fraction_cut", 0.5);

    reco.firstPassSeedJetD = optional<int>(background["iteration_1"], "seed_jet_D", 3);

    reco.secondPassSeedPtGeV = optional<double>(background["iteration_2"], "seed_pt_min_gev", 7.0);

    reco.flowModulation = required<int>(background, "flow_modulation", "auau.background_subtraction");

    tree.calorimeter.backgroundFlowMode = reco.flowModulation;
  }


  // ==========================================================================
  // Reconstructed jets and jet-energy-scale contract
  // ==========================================================================

  tree.jets.truthMatchMaxDeltaR = optional<double>(root["jets"]["matching"], "max_delta_r", 0.30);

  const ConfigurationNode calibration = root["jets"]["calibration"];

  cond.applyJetEnergyScale = optional<bool>(calibration, "enabled", true);

  cond.jetEnergyScaleCdbKey = optionalString(calibration, "cdb_key", "JES_Calib_Default");

  cond.jetEnergyScalePayloadSha256 =
      lower(required<std::string>(calibration, "expected_payload_sha256", "jets.calibration"));

  cond.jetEnergyScaleUsesEmFraction = optional<bool>(calibration, "use_em_fraction_calibration", false);

  tree.jets.energyScalePayloadSha256 = cond.jetEnergyScalePayloadSha256;

  tree.jets.energyScaleUsesEmFraction = cond.jetEnergyScaleUsesEmFraction;

  const ConfigurationNode views = systemConfig["jets"]["views"];

  if (!views || !views.IsSequence())
  {
    fail("missing jet view list for the selected collision system");
  }

  for (const ConfigurationNode& view : views)
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


  // ==========================================================================
  // Au+Au centrality configuration
  //
  // Coarse analysis edges are retained alongside the native centile. DATA
  // centrality reconstruction can resolve either CDB-managed payloads or an
  // explicitly pinned local per-run triplet.
  // ==========================================================================

  tree.centrality.edges = intList(root["auau"]["centrality"]["analysis_bin_edges_percent"]);

  if (tree.centrality.edges.empty())
  {
    tree.centrality.edges = {0, 20, 50, 80};
  }

  tree.centrality.reconstructDataCentrality = reco.reconstructCentrality;

  if (reco.reconstructCentrality)
  {
    const ConfigurationNode data = root["auau"]["centrality"]["data"];

    tree.centrality.calibrationTag = required<std::string>(data, "calibration_set", "auau.centrality.data");

    cond.centralityCalibrationSet = tree.centrality.calibrationTag;

    cond.centralityManifestSha256 =
        lower(required<std::string>(data, "calibration_manifest_sha256", "auau.centrality.data"));

    cond.centralityLocalPayloads = lower(optionalString(data, "source", "cdb")) != "cdb";

    if (cond.centralityLocalPayloads)
    {
      resolveCentralityPayloads(cond, data, job.run);
    }
  }


  // ==========================================================================
  // Au+Au DATA CEMC status and zero-suppression calibration
  //
  // The status guard ties the current reconstruction to the original
  // calorimeter input. Missing historical status can be restored only when the
  // configuration allows it and the source binds its run's own map, pinned by
  // SHA256 (the catalogue's cemc_recovery_map).
  // ==========================================================================

  if (reco.isAuAu && reco.isData)
  {
    const ConfigurationNode status = root["auau"]["cemc_status"];

    cond.cemcStatusGuard = optional<bool>(status, "validate_paired_input_status", true);

    cond.cemcStatusRecoveryAuthorized =
        optional<bool>(status, "allow_missing_status_recovery", false) && !job.cemcRecoveryMap.empty();

    cond.cemcStatusRecoveryPayload = job.cemcRecoveryMap;
    cond.cemcStatusRecoveryPayloadSha256 = lower(job.cemcRecoveryMapSha256);
  }
  else if (!job.cemcRecoveryMap.empty())
  {
    fail("cemc_recovery_map applies to Au+Au DATA sources only");
  }

  cond.cemcZeroSuppressionCrossCalibration =
      reco.isData && optional<bool>(root["calorimeter"]["cemc"], "zero_suppression_cross_calibration", true);


  // ==========================================================================
  // Archived p+p double-interaction conditions and RNG replay
  // ==========================================================================

  if (reco.isArchivedDoubleInteraction)
  {
    cond.cdbGlobalTag = required<std::string>(profile, "cdb_global_tag", context);

    cond.cdbTimestamp = required<std::uint64_t>(profile, "cdb_timestamp", context);

    const ConfigurationNode rng = profile["historical_rng_contract"];

    reco.forbidRecoConstsRandomSeed = optional<bool>(rng, "forbid_reco_consts_random_seed", true);

    // Each archived source replays its own production's seeds, recorded in
    // its catalogue entry (historical_rng).
    if (job.historicalRng.empty())
    {
      fail("archived double interaction requires the source's historical_rng record");
    }
    parseHistoricalRng(job.historicalRng, reco.historicalSeeds, reco.expectedPedestalSequence);
  }
  else
  {
    if (!job.historicalRng.empty())
    {
      fail("historical_rng applies to archived double-interaction sources only");
    }

    cond.cdbTimestamp =
        job.run > 0 ? static_cast<std::uint64_t>(job.run) : static_cast<std::uint64_t>(kArchivedDoubleInteractionRun);

    cond.cdbGlobalTag = optionalString(profile, "cdb_global_tag");
  }


  // ==========================================================================
  // Producer-owned simulation weights
  //
  // Retain PPG12 relative event weights; absolute sample normalization and
  // source ownership remain downstream. The period selects its own vertex fit.
  // ==========================================================================

  ConfigurationNode reweight;
  if (reco.isSimulation && !reco.isAuAu && (role == "photon_jet" || role == "inclusive_jet"))
  {
    const ConfigurationNode weights = root["ppg12_weights"];
    const ConfigurationNode period = weights["periods"][tree.source.period];
    const ConfigurationNode family = weights[role];
    const ConfigurationNode sliceName = family ? family["samples"][tree.source.sample] : ConfigurationNode();
    const ConfigurationNode slice = sliceName ? family["slices"][sliceName.as<std::string>()] : ConfigurationNode();
    if (!period || !family || !slice) fail("PPG12 weights require an exact maintained sample and period binding");
    const std::string interaction = lower(tree.source.siDiRole);
    if (interaction != "si" && interaction != "di") fail("PPG12 weights require si_di_role=SI or DI");
    auto& w = tree.weights;
    w.ppg12Nominal = true;
    w.doubleInteraction = interaction == "di";
    if (w.doubleInteraction != (inputMode == "archived_g4_only"))
      fail("PPG12 interaction role differs from the selected reconstruction profile");
    w.ppg12LaneCode = role == "photon_jet" ? 1 : 2;
    w.ppg12Slice = sliceName.as<std::string>();
    w.ppg12SliceCode = required<int>(slice, "code", "ppg12_weights.slices");
    const ConfigurationNode window = slice["window_gev"];
    if (!window || !window.IsSequence() || window.size() != 2 || !window[0])
      fail("a PPG12 slice requires its ownership window [low, high]");
    w.ownershipLowGeV = window[0].as<double>();
    w.ownershipHighGeV = window[1] ? window[1].as<double>() : std::numeric_limits<double>::infinity();
    if (!(w.ownershipLowGeV >= 0.0) || !(w.ownershipHighGeV > w.ownershipLowGeV))
      fail("invalid PPG12 slice ownership window");
    if (w.ppg12LaneCode == 2)
    {
      w.clusterEtUpperGeV = required<double>(slice, "cluster_et_max_gev", "ppg12_weights.slices");
      if (!(w.clusterEtUpperGeV > 0.0)) fail("invalid inclusive-jet slice cluster ET limit");
    }
    const double sigma = required<double>(slice, "cross_section_pb", "ppg12_weights.slices");
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
  else
    reweight = profile["vertex_reweighting"];
  if (reweight)
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


  // ==========================================================================
  // Historical event classification (Au+Au)
  //
  // The established Au+Au producer's gates, in its order: DATA scaled
  // triggers, the minimum-bias classifier and |vz|; embedded simulation |vz|
  // and its generator-slice stitching; then the analysis centrality range.
  // Recorded as evt_terminal_status; nothing is removed.
  // ==========================================================================

  if (const ConfigurationNode classification = systemConfig["event_classification"]; reco.isAuAu && classification)
  {
    const std::string where = "auau.event_classification";
    EventClassificationConfig& c = tree.classification;
    c.enabled = true;

    for (const ConfigurationNode& trigger : classification["scaled_triggers"])
    {
      const int bit = required<int>(trigger, "bit", where + ".scaled_triggers");
      if (bit < 0 || bit > 63) fail(where + ": a scaled trigger needs a GL1 bit in [0, 63]");
      c.scaledTriggers.emplace_back(bit, required<std::string>(trigger, "stream", where + ".scaled_triggers"));
    }
    c.requireMinimumBias = required<bool>(classification, "require_minimum_bias", where);
    c.useVertexCut = true;
    c.vertexCutCm = required<double>(classification, "vertex_cut_cm", where);

    const ConfigurationNode range = classification["centrality_range_percent"];
    if (!range || range.size() != 2) fail(where + ".centrality_range_percent needs [low, high)");
    c.centralityMinPercent = range[0].as<int>();
    c.centralityMaxPercent = range[1].as<int>();
    if (c.centralityMinPercent < 0 || c.centralityMaxPercent <= c.centralityMinPercent)
      fail(where + ".centrality_range_percent needs 0 <= low < high");

    const ConfigurationNode windows = classification["embedded_jet_windows_gev"];
    for (const std::string& key : windows.keys())
    {
      const ConfigurationNode window = windows[key];
      if (!window || window.size() != 2 || window[0].IsNull())
        fail(where + ".embedded_jet_windows_gev." + key + " needs [low, high or null]");
      const double high = window[1].IsNull() ? std::numeric_limits<double>::infinity() : window[1].as<double>();
      c.embeddedJetWindows[std::stoi(key)] = {window[0].as<double>(), high};
    }
  }


  // ==========================================================================
  // Au+Au embedded truth-photon misses
  //
  // The established producer's leading isolated photon passing its cut-based
  // tight identification; recorded as truthmiss_* rows.
  // ==========================================================================

  if (const ConfigurationNode misses = systemConfig["truth_photon_misses"]; reco.isAuAu && reco.isSimulation && misses)
  {
    const std::string where = "auau.truth_photon_misses";
    TruthMissConfig& m = tree.truthMisses;
    m.enabled = true;
    const ConfigurationNode et = misses["et_range_gev"];
    if (!et || et.size() != 2) fail(where + ".et_range_gev needs [low, high)");
    m.etMinGeV = et[0].as<double>();
    m.etMaxGeV = et[1].as<double>();
    m.absEtaMax = required<double>(misses, "abs_eta_max", where);
    m.isolationRadius = required<double>(misses, "isolation_radius", where);
    const auto& radii = tree.photon.isolationRadii;
    if (std::find(radii.begin(), radii.end(), m.isolationRadius) == radii.end())
      fail(where + ".isolation_radius is not a retained isolation radius");
    const ConfigurationNode tight = misses["tight_id"];
    const std::string at = where + ".tight_id";
    m.widthMin = required<double>(tight, "width_min", at);
    m.widthMaxIntercept = required<double>(tight["width_max"], "intercept", at + ".width_max");
    m.widthMaxSlopePerGeV = required<double>(tight["width_max"], "slope_per_gev", at + ".width_max");
    for (const auto& [key, window] : {std::pair<const char*, std::array<double, 2>*>{"e11_over_e33", &m.e11OverE33},
                                      {"et1", &m.et1},
                                      {"e32_over_e35", &m.e32OverE35}})
    {
      const ConfigurationNode range = tight[key];
      if (!range || range.size() != 2 || !(range[0].as<double>() < range[1].as<double>()))
        fail(at + "." + key + " needs an open window [low, high]");
      *window = {range[0].as<double>(), range[1].as<double>()};
    }
  }


  // ==========================================================================
  // Historical event classification (pp)
  //
  // Recorded as evt_terminal_status and the PPG12 audit rows; nothing is
  // removed. The period contract is the source's own period; DATA run labels
  // use every period's run range.
  // ==========================================================================

  if (const ConfigurationNode classification = systemConfig["event_classification"]; reco.isPP && classification)
  {
    EventClassificationConfig& c = tree.classification;
    c.enabled = true;

    for (const ConfigurationNode& trigger : classification["triggers"])
      c.triggers.emplace_back(required<std::string>(trigger, "gl1_name", "event_classification.triggers"),
                              required<std::string>(trigger, "stream", "event_classification.triggers"));

    if (const ConfigurationNode audit = classification["scaled_bit_audit"])
    {
      c.scaledBitStream = required<std::string>(audit, "stream", "event_classification.scaled_bit_audit");
      c.scaledBit = required<int>(audit, "bit", "event_classification.scaled_bit_audit");
      if (c.scaledBit < 0 || c.scaledBit > 63) fail("the scaled-bit audit stream needs a GL1 bit in [0, 63]");
    }

    c.useVertexCut = true;
    c.vertexCutCm = required<double>(classification, "vertex_cut_cm", "event_classification");
    c.useLuminosityWeight = optional<bool>(classification, "use_luminosity_weight", true);

    const ConfigurationNode weights = root["ppg12_weights"];
    const double target = required<double>(weights, "total_luminosity_pb_inv", "ppg12_weights");
    for (const std::string& key : weights["periods"].keys())
    {
      const ConfigurationNode node = weights["periods"][key];
      const ConfigurationNode runs = node["runs"];
      if (!runs || runs.size() != 2) fail("PPG12 period " + key + " needs runs: [first, last)");

      PPG12PeriodConfig period;
      period.key = key;
      period.runMin = runs[0].as<int>();
      period.runMaxExclusive = runs[1].as<int>();
      period.luminosityPbInv = required<double>(node, "luminosity_pb_inv", "ppg12_weights.period");
      period.targetLuminosityPbInv = target;
      period.doubleFraction = required<double>(node, "double_fraction", "ppg12_weights.period");
      period.closureZCutCm = required<double>(node, "closure_z_cut_cm", "ppg12_weights.period");
      c.periods.push_back(period);
      if (key == tree.source.period) c.period = period;
    }
    if (!tree.source.period.empty() && c.period.key.empty()) fail("unknown PPG12 period: " + tree.source.period);

    const bool contract = !c.period.key.empty();
    c.filterDataByPeriod =
        contract && !reco.isSimulation && optional<bool>(classification, "filter_data_by_period", true);

    // SIM takes its SI/DI mixture weight and truth-vertex fit from the period.
    if (reco.isSimulation && contract)
    {
      c.doubleInteraction = tree.weights.doubleInteraction;
      c.mixWeight = c.doubleInteraction ? c.period.doubleFraction : 1.0 - c.period.doubleFraction;
      c.mixWeightFromPeriod = true;
      c.vertexFileFromPeriod = true;
    }

    // The PPG12 bookkeeping of classified events.
    const ConfigurationNode book = systemConfig["ppg12_bookkeeping"];
    const auto edges = [&](const char* key, std::vector<double>& target)
    {
      for (const ConfigurationNode& edge : book[key]) target.push_back(edge.as<double>());
      if (target.size() < 2 || !std::is_sorted(target.begin(), target.end()) ||
          std::adjacent_find(target.begin(), target.end()) != target.end())
        fail(std::string("ppg12_bookkeeping.") + key + " needs at least two increasing edges");
    };
    edges("photon_et_bins_gev", tree.ppg12.photonEtBins);
    edges("unfold_reco_photon_et_bins_gev", tree.ppg12.unfoldRecoPhotonEtBins);
    edges("unfold_truth_photon_et_bins_gev", tree.ppg12.unfoldTruthPhotonEtBins);
    edges("unfold_xj_bins", tree.ppg12.unfoldXJBins);
    tree.ppg12.recoilJetPtMinGeV = required<double>(book, "recoil_jet_pt_min_gev", "ppg12_bookkeeping");
    // The producer's cut: the configured fraction times pi.
    tree.ppg12.recoilDeltaPhiMin = required<double>(book, "recoil_delta_phi_min_pi", "ppg12_bookkeeping") * M_PI;
    tree.ppg12.jetMatchMaxDeltaR = required<double>(book, "jet_match_max_delta_r", "ppg12_bookkeeping");
    tree.ppg12.responseSeedRun =
        reco.isArchivedDoubleInteraction
            ? cond.cdbTimestamp
            : required<std::uint64_t>(book, "established_simulation_timestamp", "ppg12_bookkeeping");
    if (!(tree.ppg12.recoilJetPtMinGeV > 0.0) || !(tree.ppg12.recoilDeltaPhiMin > 0.0) ||
        !(tree.ppg12.jetMatchMaxDeltaR > 0.0))
      fail("invalid ppg12_bookkeeping recoil definition");
  }


  // ==========================================================================
  // Optional high-volume output tables
  // ==========================================================================

  const ConfigurationNode output = root["output"];

  // The branch table built into libEventTree from this schema file (a changed
  // file is refused). Storage widths never depend on collection sizes.
  const auto schemaPath = std::filesystem::path(yamlFile).parent_path() / "eventTree_v1.yaml";
  const eventtree::Schema& schema = eventtree::Schema::load(schemaPath.string());
  tree.output.schemaSha256 = schema.sha256();
  for (const eventtree::Branch& branch : schema.branches())
    tree.output.fields.push_back(
        {branch.name, branch.type, branch.title, branch.counter, branch.length, branch.reserved});
  tree.output.modelFeatures = tree.photon.nominalId.featureOrder;

  tree.output.writePhotonCells = optional<bool>(output, "write_photon_cell_table", true);

  tree.output.writeIsolationConstituents = optional<bool>(output, "write_isolation_constituent_table", true);

  tree.output.writePhotonJetPairs = optional<bool>(output, "write_photon_jet_pair_table", true);


  // ==========================================================================
  // Event-retention policy and provenance
  // ==========================================================================

  tree.retainEveryProducerEncounter = optional<bool>(root["capture"], "retain_every_producer_encounter", true);

  tree.provenance.productionTag = required<std::string>(root["contract"], "name", "contract");

  tree.provenance.softwareRelease = optionalString(root["runtime"], "sphenix_release");

  tree.provenance.producerGitCommit = optionalString(root["runtime"], "producer_git_commit");

  // Bind the loaded library and compiled sources, not a manually copied label.
  Dl_info library{};
  if (!dladdr(reinterpret_cast<const void*>(&LoadPlan), &library) || !library.dli_fname)
    fail("cannot identify the loaded PhotonJetTree library");
  auto pinnedRuntime = [&](const std::string& key, const std::string& actual)
  {
    if (!validSha256(actual)) fail("missing runtime digest: " + key);
    const std::string expected = optionalString(root["runtime"], key);
    if (!expected.empty() && expected != actual) fail("runtime digest differs: " + key);
    return actual;
  };
  tree.provenance.producerSourceSha256 = pinnedRuntime("producer_source_sha256", PHOTONJET_SOURCE_SHA256);
  tree.provenance.producerLibrarySha256 =
      pinnedRuntime("producer_library_sha256", detail::fileSha256(library.dli_fname));
  tree.provenance.macroSha256 = pinnedRuntime("macro_sha256", job.macroSha256);

  tree.provenance.configurationSha256 = plan.configurationSha256;

  tree.provenance.configurationText = plan.configurationText;

  tree.provenance.calibrationManifestSha256 = cond.centralityManifestSha256;


  // ==========================================================================
  // Physical inputs
  //
  // Input-list interpretation belongs to Setup.cc. At this point the
  // profile is sufficiently resolved to determine which physical streams are
  // required and how they must be synchronized.
  // ==========================================================================

  plan.inputs = loadInputs(job, reco);

  // The CEMC status guard reads the RUN metadata of this job's own
  // calorimeter file: the paired DATA input's primary stream. Embedded
  // simulation uses the conditions of the data run it is embedded in, which
  // its calorimeter DST name carries as "-data-<run>-" (the established
  // producer's CDB run).
  for (const InputStream& stream : plan.inputs)
  {
    if (stream.files.size() != 1) continue;
    if (cond.cemcStatusGuard && stream.name == "PRIMARY") cond.cemcOriginalCalorimeterInput = stream.files.front();
    if (reco.isEmbedded && stream.name == "CALO") cond.cdbTimestamp = embeddedDataRun(stream.files.front());
  }

  /*
   * Auto truth-jet mode is resolved only after the physical input inventory is
   * known. If shipped truth jets are present they are consumed; otherwise the
   * reconstruction layer must build them.
   */
  if (reco.truthJetMode == TruthJetMode::Auto)
  {
    reco.truthJetMode = hasInputStream(plan, "TRUTH_JETS") ? TruthJetMode::DST : TruthJetMode::Build;
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
  const Reconstruction& reco = plan.reconstruction;

  const Conditions& cond = plan.conditions;


  // ==========================================================================
  // Input and event-range integrity
  // ==========================================================================

  if (plan.inputs.empty())
  {
    fail("no input streams were resolved");
  }

  if (plan.job.firstEntry < 0 || plan.job.firstEntry > std::numeric_limits<int>::max() || plan.job.nEvents < -1 ||
      plan.job.nEvents == 0 || plan.job.nEvents > std::numeric_limits<int>::max())
  {
    fail("invalid event range: use -1 for all events, or a positive bounded count");
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
      fail("one source bundle per invocation is required: " + stream.name);
    }
  }

  if (cond.cdbGlobalTag.empty())
  {
    fail("profile must bind an explicit CDB global tag");
  }

  if (!cond.applyJetEnergyScale)
  {
    fail("canonical reconstructed jets require JES");
  }

  if (cond.jetEnergyScaleCdbKey != "JES_Calib_Default")
  {
    fail("JES key must match the key consumed by JetCalib");
  }

  if (cond.applyVertexReweight && !validSha256(cond.vertexReweightFileSha256))
  {
    fail("vertex weight payload requires SHA256");
  }

  for (const InputStream& stream : plan.inputs)
  {
    if (stream.required && stream.files.empty())
    {
      fail("required input stream is empty: " + stream.name);
    }
  }


  // ==========================================================================
  // DATA / simulation role consistency
  // ==========================================================================

  if (reco.isData && plan.tree.simulationRole != SimulationRole::None)
  {
    fail("a data profile cannot carry a simulation role");
  }

  if (reco.isSimulation && plan.tree.simulationRole == SimulationRole::None)
  {
    fail("a simulation profile must name its sample role");
  }


  // ==========================================================================
  // Archived p+p double-interaction production
  //
  // Geant4 hits with shipped truth jets, reconstructed in the run-28
  // conditions with the source's own replayed random sequence (Setup.cc).
  // ==========================================================================

  if (reco.isArchivedDoubleInteraction)
  {
    if (!reco.isPP || !reco.isSimulation)
    {
      fail("archived double interaction is p+p simulation only");
    }

    if (reco.truthJetMode != TruthJetMode::DST)
    {
      fail("archived double interaction requires shipped truth jets");
    }

    if (cond.cdbTimestamp != static_cast<std::uint64_t>(kArchivedDoubleInteractionRun))
    {
      fail("archived double interaction requires the run-28 conditions timestamp");
    }

    // The replayed random sequence starts with the source's first event.
    if (plan.job.firstEntry != 0)
    {
      fail("an archived double-interaction source is replayed from its first event");
    }
  }


  // ==========================================================================
  // Real-run requirement for DATA and embedded production
  // ==========================================================================

  if ((reco.isData && plan.job.run <= kDataRunThreshold) ||
      (reco.isEmbedded && plan.conditions.cdbTimestamp <= kDataRunThreshold))
  {
    fail("data production requires a real run number, embedded simulation its embedded data run");
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
      if (cond.centralityPayloadRun <= 0)
      {
        fail("local centrality requires a positive payload run");
      }

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


  // ==========================================================================
  // Au+Au DATA CEMC tower-status contract
  // ==========================================================================

  if (cond.cemcStatusGuard)
  {
    if (cond.cemcOriginalCalorimeterInput.empty())
    {
      fail("the CEMC status guard requires the exact original calorimeter input file");
    }

    // The recovery map is the run's own CDB payload, <...>_<run>cdb.root.
    const std::string& map = cond.cemcStatusRecoveryPayload;
    const std::string ownRun = "_" + std::to_string(plan.job.run) + "cdb.root";
    if (cond.cemcStatusRecoveryAuthorized &&
        (!validSha256(cond.cemcStatusRecoveryPayloadSha256) || map.size() < ownRun.size() ||
         map.compare(map.size() - ownRun.size(), ownRun.size(), ownRun) != 0))
    {
      fail("CEMC status recovery requires the run's own bad-tower map and its SHA256");
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
      fail("the nominal jet energy scale is the legacy method without the electromagnetic-fraction calibration");
    }

    if (cond.jetEnergyScalePayloadSha256 != kNominalJetEnergyScaleSha256)
    {
      fail("jet energy-scale payload SHA256 differs from the nominal common p+p v6 payload");
    }
  }


  // ==========================================================================
  // Producer-side vertex reweighting
  // ==========================================================================

  if (cond.applyVertexReweight)
  {
    if (!reco.isSimulation)
    {
      fail("vertex reweighting applies to simulation only");
    }

    if (cond.vertexReweightFile.empty() || cond.vertexReweightHistogram.empty())
    {
      fail("vertex reweighting requires a file and a histogram name");
    }
  }


  // ==========================================================================
  // Au+Au flow-modulation configuration
  // ==========================================================================

  if (reco.isAuAu && (reco.flowModulation < 0 || reco.flowModulation > 3))
  {
    fail("Au+Au flow modulation must be resolved to the accepted production value");
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
  const Reconstruction& reco = plan.reconstruction;

  const Conditions& cond = plan.conditions;

  const char* truthMode = "none";

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
            << " JES                  : "
            << (cond.applyJetEnergyScale ? "legacy, once, " + cond.jetEnergyScalePayload : "disabled") << "\n";

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
    std::cout << "   " << std::left << std::setw(12) << stream.name << " files=" << stream.files.size()
              << " sync=" << (stream.synchronized ? "yes" : "no") << "\n";
  }

  std::cout << "============================================================\n" << std::endl;
}

}  // namespace production
}  // namespace photonjet
