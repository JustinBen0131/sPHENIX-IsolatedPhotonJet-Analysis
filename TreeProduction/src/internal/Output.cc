/**
 * @file Output.cc
 * Write one eventTree per source, directly from the in-memory event records.
 *
 * DirectEventWriter owns the ROOT storage of one file. The producer fills it in
 * three steps: writeFileMetadata() when the file opens (source binding and
 * configuration), fillEventOutput() or fillRejectedOutput() once per event in
 * source order, and writeCompletionMetadata() when the job ends (accounting,
 * code dictionaries, constants, policies, the run trigger table and the
 * photon-ID and isolation receipts).
 *
 * Every schema branch has exactly one filling rule here. filledBranches() at
 * the end of this file lists them; a schema branch without one is reported in
 * meta/native_capture and refuses production before the first event is read.
 * A value is written only from what reconstruction measured: storage defaults
 * are never evidence of a capture.
 */
#include "../PhotonJetTree.h"
#include "Production.h"
#include "../EventTree.h"

#include <nlohmann/json.hpp>

#include <TH1D.h>
#include <TROOT.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
// The schema branches this file fills (defined at the end of the file).
std::set<std::string> filledBranches();
}  // namespace

namespace photonjet
{
/**
 * The storage of one eventTree file: the schema-bound entry, the ROOT writer,
 * and what the completion metadata gathers while events are written (code
 * dictionaries, file constants, policies, counters and, for Au+Au simulation,
 * the centrality denominator).
 */
class DirectEventWriter
{
 public:
  std::vector<eventtree::FieldSpec> specs;  // the schema branches, in schema order
  eventtree::Entry entry;                   // the current event's values
  eventtree::Writer writer;
  std::set<std::string> covered, missing;              // schema branches with and without a filling rule
  nlohmann::json meta;                                 // file constants (meta/constants)
  nlohmann::json codes = nlohmann::json::object();     // one dictionary per coded branch (meta/codes)
  nlohmann::json policies = nlohmann::json::object();  // file-level rules (meta/policies)
  std::uint64_t entries = 0;
  std::unique_ptr<TH1D> centrality;  // Au+Au simulation: hist/SIM/h_centrality
  std::int64_t nextEntry;            // the source entry the next event must have
  std::set<std::string> filledReserved;
  std::set<std::string> truthMissPaths;  // the decision paths of the file's truth misses
  std::map<std::string, std::int64_t> invalidUnsubtractedJes;  // Au+Au: per unsubtracted view
  // What the file's events held. A file constant that describes rows exists only
  // when such rows were written, as the established output records it.
  std::uint64_t photonRows = 0, truthJetRows = 0;
  std::set<std::string> jetCollectionsWithJets;
  std::size_t jetViewsRecorded = 0;     // the most jet views one event recorded
  std::size_t truthRadiiRecorded = 0;   // the most truth-jet containers one event recorded
  static std::vector<eventtree::FieldSpec> fields(const Config& c)
  {
    std::vector<eventtree::FieldSpec> result;
    for (const auto& f : c.output.fields) result.push_back({f.name, f.type, f.title, f.counter, f.length});
    if (result.empty()) throw std::runtime_error("eventTree schema is not bound to the production plan");
    return result;
  }
  explicit DirectEventWriter(const Config& c)
      : specs(fields(c)),
        entry(specs),
        writer(c.outputFile, "eventTree", "RJEventTreeV1", 404, 0),
        nextEntry(c.source.firstEntry)
  {
    const auto rules = filledBranches();
    for (const auto& field : specs) (rules.count(field.name) ? covered : missing).insert(field.name);
    // Reserved branches this producer fills natively are declared filled.
    for (const auto& field : c.output.fields)
      if (field.reserved && covered.count(field.name)) filledReserved.insert(field.name);
    if (c.system == CollisionSystem::AuAu && c.dataKind == DataKind::Simulation)
    {
      centrality = std::make_unique<TH1D>("h_centrality", ";Centrality [%];Events", 100, 0., 100.);
      centrality->SetDirectory(nullptr);
      centrality->Sumw2();
    }
  }
  template <class T>
  void put(const std::string& n, T v)
  {
    entry.put(n, v);
  }
  template <class T>
  void all(const std::string& n, const T& v)
  {
    entry.putAll(n, v);
  }
  void empty(const std::string& n) { entry.empty(n); }
  void begin() { entry.reset(); }
  // One dictionary per coded branch (meta/codes/<branch>), the shared
  // consolidation and check contract; -1 encodes "no value" (a valid state).
  int code(const std::string& dictionary, const std::string& value)
  {
    if (value.empty()) return -1;
    auto& values = codes[dictionary];
    if (values.is_null()) values = nlohmann::json::array();
    for (std::size_t i = 0; i < values.size(); ++i)
      if (values[i] == value) return static_cast<int>(i);
    values.push_back(value);
    return static_cast<int>(values.size() - 1);
  }
  void finish(bool rejected = false, bool bookOnly = false)
  {
    // Empty collections still have a filling rule. Missing rules remain
    // distinct from an observed empty collection in the terminal metadata.
    for (const auto& f : specs)
      if (f.length == -1 && covered.count(f.name)) entry.empty(f.name);
    for (const auto& name : entry.missing())
    {
      const auto& f = *std::find_if(specs.begin(), specs.end(), [&](const auto& s) { return s.name == name; });
      if (!rejected && !bookOnly && covered.count(name) && f.length != -1)
        throw std::runtime_error("declared native scalar rule did not fill field: " + name);
      if (f.length == -1)
        entry.empty(name);
      else
        for (int i = 0; i < std::max(1, f.length); ++i)
        {
          if (f.type == "double" || f.type == "float")
            entry.put(name, std::numeric_limits<double>::quiet_NaN());
          else
            entry.put(name, 0);
        }
    }
    entry.write(writer, bookOnly);
    if (!bookOnly) ++entries;
  }
};
}

namespace
{
using json = nlohmann::json;
constexpr double missingValue = std::numeric_limits<double>::quiet_NaN();
using Key = std::pair<std::uint64_t, std::uint64_t>;
Key key(const photonjet::Identity& i) { return {i.hi, i.lo}; }
template <class Row, class Member>
std::map<Key, int> indices(const std::vector<Row>& rows, Member member)
{
  std::map<Key, int> out;
  for (std::size_t i = 0; i < rows.size(); ++i)
    if (!out.emplace(key(rows[i].*member), static_cast<int>(i)).second)
      throw std::runtime_error("duplicate event-local object identity");
  return out;
}
int local(const std::map<Key, int>& objects, const photonjet::Identity& id, bool optional = false)
{
  if (!id.valid() && optional) return -1;
  const auto it = objects.find(key(id));
  if (it == objects.end()) throw std::runtime_error("dangling event-local object reference");
  return it->second;
}
json bits(double value)
{
  std::uint64_t raw;
  std::memcpy(&raw, &value, sizeof(raw));
  return {{"f64", raw}};
}
template <class T>
void appendLine(std::ostringstream& out, const std::string& name, const T& value)
{
  out << name << '=' << value << '\n';
}

// Per-file receipts of the nominal labels: the frozen definitions they bind,
// the code dictionaries and the source totals of every value.
json tally(const std::map<int, std::uint64_t>& counts, std::initializer_list<int> codes)
{
  json out = json::object();
  for (const int code : codes)
  {
    const auto it = counts.find(code);
    out[std::to_string(code)] = it == counts.end() ? 0 : it->second;
  }
  return out;
}
json nominalPhotonIdReceipt(const photonjet::NominalPhotonIdConfig& id, const photonjet::NominalReceiptRecord& r)
{
  json receipt = {{"schema", "RJEventTreePhotonIDReceiptV1"},
                  {"freeze_sha256", id.packageSha256},
                  {"system", id.system},
                  {"model_sha256", id.modelSha256},
                  {"photons", r.photons},
                  {"codes",
                   {{"pho_bdt_state",
                     {{"1", "scored"},
                      {"0", "outside the model domain"},
                      {"2", "inside the domain with an invalid input feature"}}},
                    {"pho_id_class",
                     {{"3", "tight"},
                      {"2", "between non-tight and tight (excluded)"},
                      {"1", "non-tight"},
                      {"0", "below non-tight"},
                      {"-1", "not scored"}}}}},
                  {"states", tally(r.states, {1, 0, 2})},
                  {"classes", tally(r.idClasses, {3, 2, 1, 0, -1})}};
  if (r.photons)
  {
    receipt["root_version"] = gROOT->GetVersion();
    if (id.usesProducerScore)
    {
      receipt["pp_recompute_max_float32_ulp"] = r.maxRecomputeUlp;
      receipt["pp_recompute_region_flips"] = r.recomputeRegionFlips;
    }
  }
  return receipt;
}
json nominalIsolationReceipt(const photonjet::NominalIsolationConfig& c, const photonjet::NominalReceiptRecord& r)
{
  json receipt = {{"schema", "RJEventTreeIsolationReceiptV1"},
                  {"freeze_sha256", c.freezeSha256},
                  {"system", c.system},
                  {"photons", r.photons},
                  {"derivation_sha256", c.derivationSha256},
                  {"codes",
                   {{"pho_iso_class",
                     {{"2", "isolated"},
                      {"1", "between isolated and non-isolated (excluded)"},
                      {"0", "non-isolated"},
                      {"-1", "invalid input"}}}}},
                  {"classes", tally(r.isolationClasses, {2, 1, 0, -1})}};
  if (!c.definitionJson.empty()) receipt["definition"] = json::parse(c.definitionJson);
  return receipt;
}
}

std::string photonjet::production::OutputMappingStatusJson(const Config& config)
{
  const auto rules = filledBranches();
  nlohmann::json implemented = nlohmann::json::array(), missing = nlohmann::json::array();
  for (const auto& field : config.output.fields)
    (rules.count(field.name) ? implemented : missing).push_back(field.name);
  return nlohmann::json{{"schema", "NativeOutputMappingV1"},
                        {"schema_sha256", config.output.schemaSha256},
                        {"implemented", implemented},
                        {"missing_branches", missing},
                        {"ready", !config.output.fields.empty() && missing.empty()}}
      .dump();
}

void PhotonJetTree::initializeOutput()
{
  m_directOutput = std::make_shared<photonjet::DirectEventWriter>(m_config);
  writeFileMetadata();
  m_outputReady = true;
}

void photonjet::production::WriteProgress(const Config& c, const std::string& phase, std::int64_t processed,
                                          std::int64_t skipped, const std::string& error, bool force)
{
  struct State
  {
    std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now(), written = started;
    std::int64_t processed = 0, skipped = 0;
    std::string phase;
  };
  static std::map<std::string, State> states;
  auto& state = states[c.outputFile];
  if (c.source.maxEvents < 0) throw std::runtime_error("progress requires an exact finite event budget");
  const auto skipTotal = c.source.maxEvents == 0 ? 0 : c.source.firstEntry;
  if (processed < 0) processed = state.processed;
  if (skipped < 0) skipped = state.skipped;
  if (processed < state.processed || processed > c.source.maxEvents || skipped < state.skipped || skipped > skipTotal)
    throw std::runtime_error("invalid or regressing native progress count");
  state.processed = processed;
  state.skipped = skipped;
  const auto now = std::chrono::steady_clock::now();
  if (!force && state.phase == phase && processed % 128 != 0 &&
      std::chrono::duration<double>(now - state.written).count() < 5.)
    return;
  const double seconds = std::chrono::duration<double>(now - state.started).count();
  const double rate = seconds > 0 ? (processed + skipped) / seconds : 0.;
  const auto remaining = c.source.maxEvents + skipTotal - processed - skipped;
  const auto stamp = std::time(nullptr);
  std::ostringstream utc;
  utc << std::put_time(std::gmtime(&stamp), "%Y-%m-%dT%H:%M:%SZ");
  json record = {{"schema", "NativeTreeProductionProgressV1"},
                 {"phase", phase},
                 {"dataset", c.source.dataset},
                 {"sample", c.source.sample},
                 {"run", c.source.run},
                 {"segment", c.source.segment},
                 {"source_file_ordinal", c.source.sourceFileOrdinal},
                 {"source_manifest_sha256", c.source.sourceManifestSha256},
                 {"output", c.outputFile},
                 {"processed", processed + skipped},
                 {"total", c.source.maxEvents + skipTotal},
                 {"remaining", remaining},
                 {"row_processed", processed},
                 {"row_total", c.source.maxEvents},
                 {"row_remaining", c.source.maxEvents - processed},
                 {"deterministic_skip_processed", skipped},
                 {"deterministic_skip_total", skipTotal},
                 {"deterministic_skip_remaining", skipTotal - skipped},
                 {"rate", rate},
                 {"eta_seconds", rate > 0 ? json(remaining / rate) : json(nullptr)},
                 {"updated_at", utc.str()},
                 {"latest_error", error}};
  const auto target = c.outputFile + ".progress.json", temporary = target + ".tmp";
  {
    std::ofstream out(temporary, std::ios::trunc);
    out << record.dump() << '\n';
    out.flush();
    if (!out) throw std::runtime_error("cannot write native progress: " + temporary);
  }
  if (std::rename(temporary.c_str(), target.c_str()))
    throw std::runtime_error("cannot publish native progress: " + target);
  state.written = now;
  state.phase = phase;
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
  appendLine(text, "source_global_entry_begin", m_config.source.sourceGlobalEntryBegin);
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
  appendLine(text, "resolved_configuration_sha256", m_config.provenance.resolvedConfigurationSha256);
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
  appendLine(text, "ppg12_tower_mask_sha256", m_config.photon.ppg12TowerMaskSha256);
  appendLine(text, "ppg12_tower_mask_histogram", m_config.photon.ppg12TowerMaskHistogram);

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
    key << "jet_collection_view" << static_cast<int>(collection.view) << "_r"
        << static_cast<int>(photonjet::radiusCode(collection.radius));
    appendLine(text, key.str(),
               collection.rawNode + "|" + collection.correctedNode + "|" + collection.truthNode + "|" +
                   collection.inputIdentity + "|" + collection.subtractionIdentity);
  }

  m_directOutput->writer.text("metadata", text.str());
  m_directOutput->writer.text("configuration", m_config.provenance.configurationText);
}

namespace
{
// The photon's NPB model invocations, its NPB measurement in the shared
// five-state contract and the NPB input features.
void fillNpbBranches(photonjet::DirectEventWriter& d, const photonjet::Config& config, const photonjet::PhotonRecord& r)
{
  // Preserve exact optional NPB invocation witnesses. Version zero means no
  // invocation was observed; neither a score nor applicability is inferred.
  for (const auto& item : std::vector<std::pair<std::string, const photonjet::ModelInvocationRecord*>>{
           {"pho_npbcap_", &r.npbCapture}, {"pho_aunpbcap_", &r.auauNpbCapture}})
  {
    const auto& n = item.first;
    const auto& v = *item.second;
    d.put(n + "capture_version", v.captureVersion);
    d.put(n + "evaluation_state", v.evaluationState);
    d.put(n + "domain_mode", v.domainMode);
    d.put(n + "in_producer_domain", v.inProducerDomain);
    d.put(n + "model_available", v.modelAvailable);
    d.put(n + "input_order_is_npb25", v.inputOrderIsNpb25);
    d.put(n + "inputs_captured", v.inputsCaptured);
    d.put(n + "inputs_finite", v.inputsFinite);
    d.put(n + "min_et", v.minEt);
    d.put(n + "max_et", v.maxEt);
    d.put(n + "max_abs_eta", v.maxAbsEta);
    d.put(n + "raw_score", v.rawScore);
    d.put(n + "source_photon_id", v.sourcePhotonId);
    d.put(n + "source_cluster_et", v.sourceClusterEt);
    d.put(n + "source_eta", v.sourceEta);
    d.put(n + "source_phi", v.sourcePhi);
    d.put(n + "binding_policy", v.bindingPolicy);
    d.put(n + "source_identity_available", v.sourceIdentityAvailable);
    d.put(n + "source_native_cluster_key", v.sourceNativeClusterKey);
    d.put(n + "source_encounter_ordinal", v.sourceEncounterOrdinal);
    d.put(n + "kinematic_input_bits_available", v.kinematicInputBitsAvailable);
    d.put(n + "input_et_bits", v.inputEtBits);
    d.put(n + "input_eta_bits", v.inputEtaBits);
    // Frozen schema stores 25 slots per candidate. Absent/partial input
    // capture retains missing slots and zero finite flags, never fake data.
    for (std::size_t i = 0; i < 25; ++i)
    {
      d.put(n + "inputs", i < v.orderedInputs.size() ? v.orderedInputs[i] : std::numeric_limits<float>::quiet_NaN());
      d.put(n + "input_finite", i < v.orderedInputFinite.size() ? v.orderedInputFinite[i] : 0);
    }
    const std::string source = item.first == "pho_npbcap_" ? "photons.npb_capture_" : "photons.auau_npb_capture_";
    const json constants = {{source + "source_cluster_node", v.sourceClusterNode},
                            {source + "model_file", v.modelFile},
                            {source + "feature_names", v.featureNames}};
    for (auto constant = constants.begin(); constant != constants.end(); ++constant)
    {
      if (d.meta.contains(constant.key()) && d.meta[constant.key()] != constant.value())
        throw std::runtime_error("model invocation constant changed within source: " + constant.key());
      d.meta[constant.key()] = constant.value();
    }
  }
  // The producer computed its own NPB score only where its preselection or
  // PPG12 table QA used it: p+p. Au+Au keeps the builder's score (pho_npb_score).
  d.put("pho_npb_score_native", config.system == photonjet::CollisionSystem::PP
                                    ? r.npbCapture.rawScore
                                    : std::numeric_limits<float>::quiet_NaN());
  d.put("pho_auau_npb_score_native", r.auauNpbCapture.rawScore);
  // The collaborator NPB fields follow the NPB score the preselection uses:
  // the PPG12 split model's npb_score in both systems, as the historical
  // augmentation defined it. Each capture record keeps its own model's cut.
  const bool auauScore = config.photon.preselection.npbScoreName == "auau_npb_score";
  const auto& npb = auauScore ? r.auauNpbCapture : r.npbCapture;
  const auto configuredCut = [&](const char* name)
  {
    double value = missingValue;
    for (const auto& model : config.photon.models)
      if (model.scoreName == name) value = model.configuredCut;
    return value;
  };
  const double cut = configuredCut(auauScore ? "auau_npb_score" : "npb_score");
  d.put("pho_aug_donor_evidence", d.code("pho_aug_donor_evidence", config.provenance.producerSourceSha256));
  d.put("pho_npbcap_association", npb.sourceIdentityAvailable ? npb.bindingPolicy : 0);
  d.put("pho_npbcap_configured_cut", configuredCut("npb_score"));
  d.put("pho_aunpbcap_configured_cut", configuredCut("auau_npb_score"));
  d.put("pho_npb_cut", cut);
  // The NPB measurement in the shared contract. An invocation without a bound
  // scoring object or with nonfinite accessor kinematics is missing data.
  constexpr std::int32_t Missing = 0, Valid = 1, Invalid = 2, NotEvaluated = 3, NotApplicable = 4;
  auto accessor = [](std::uint32_t bits)
  {
    float value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
  };
  const char* missingReason = nullptr;
  if (npb.captureVersion != 1)
    missingReason = "native NPB capture absent";
  else if (!npb.sourceIdentityAvailable || !npb.kinematicInputBitsAvailable ||
           (npb.bindingPolicy != 1 && npb.bindingPolicy != 2))
    missingReason = "NPB scoring-object identity or actual kinematic accessor bits unavailable";
  else if (!std::isfinite(accessor(npb.inputEtBits)) || !std::isfinite(accessor(npb.inputEtaBits)))
    missingReason = "nonfinite native NPB kinematic inputs; raw capture retained, valid scoring binding unavailable";
  const bool measured = missingReason == nullptr;
  const bool inside = npb.inProducerDomain == 1, inputsCaptured = npb.inputsCaptured == 1,
             inputsFinite = npb.inputsFinite == 1;
  const double score = npb.rawScore;
  // evaluation: 0 missing, 1 evaluated, 2 not evaluated, 3 inputs unavailable,
  // 4 not applicable, 5 evaluated invalid; applicability: 0 unverified, 1, 2.
  std::int32_t evaluation = 0, scoreState = Missing;
  const char* scoreReason = "";
  if (!measured)
    scoreReason = missingReason;
  else if (npb.evaluationState == 4 && inputsCaptured && inputsFinite && score >= 0. && score <= 1.)
  {
    evaluation = 1;
    scoreState = Valid;
  }
  else if (npb.evaluationState == 4 || npb.evaluationState == 5)
  {
    evaluation = 5;
    scoreState = Invalid;
    scoreReason = "native Compute ran with invalid inputs or invalid output";
  }
  else if (!inside)
  {
    evaluation = 4;
    scoreState = NotApplicable;
    scoreReason = "outside the bound native NPB domain";
  }
  else
  {
    evaluation = inputsFinite ? 2 : 3;
    scoreState = NotEvaluated;
    scoreReason = "native Compute did not run";
  }
  d.put("pho_npb_applicability", measured ? (inside ? 1 : 2) : 0);
  d.put("pho_npb_evaluation", evaluation);
  d.put("pho_npb_score", scoreState == Missing || !std::isfinite(score) ? missingValue : score);
  d.put("pho_npb_score_state", scoreState);
  d.put("pho_npb_score_reason", d.code("pho_npb_score_reason", scoreReason));
  std::int32_t passState = Valid;
  const char* passReason = "";
  if (!measured)
  {
    passState = Missing;
    passReason = "no NPB donor";
  }
  else if (evaluation == 4)
  {
    passState = NotApplicable;
    passReason = "outside the explicitly bound NPB domain";
  }
  else if (evaluation == 5)
  {
    passState = Invalid;
    passReason = "native Compute ran but its inputs or output were invalid";
  }
  else if (evaluation != 1)
  {
    passState = NotEvaluated;
    passReason = "no valid evaluated NPB score";
  }
  else if (!std::isfinite(cut))
  {
    passState = Missing;
    passReason = "configured NPB cut was not captured/bound";
  }
  d.put("pho_npb_pass", passState == Valid ? int(score > cut) : -1);
  d.put("pho_npb_pass_state", passState);
  d.put("pho_npb_pass_reason", d.code("pho_npb_pass_reason", passReason));
  // Scoring-object witness: the same object (1) or the first native
  // kinematic match (2), whose encounter ordinal is then recorded.
  d.put("pho_npb_scoring_policy", measured ? npb.bindingPolicy : 0);
  d.put("pho_npb_scoring_node", measured ? d.code("pho_npb_scoring_node", npb.sourceClusterNode) : -1);
  d.put("pho_npb_scoring_key", measured ? npb.sourceNativeClusterKey : 0U);
  d.put("pho_npb_scoring_encounter",
        measured && npb.bindingPolicy == 2 ? npb.sourceEncounterOrdinal : std::int64_t{-1});
  d.put("pho_npb_scoring_pt", measured ? static_cast<double>(npb.sourceClusterEt) : missingValue);
  d.put("pho_npb_scoring_eta", measured ? static_cast<double>(npb.sourceEta) : missingValue);
  d.put("pho_npb_scoring_phi", measured ? static_cast<double>(npb.sourcePhi) : missingValue);
  d.put("pho_npb_scoring_et_bits", measured ? npb.inputEtBits : 0U);
  d.put("pho_npb_scoring_eta_bits", measured ? npb.inputEtaBits : 0U);
  d.put("pho_npb_scoring_evidence",
        measured ? d.code("pho_npb_scoring_evidence", config.provenance.producerSourceSha256) : -1);
  for (std::size_t i = 0; i < photonjet::npbFeatureNames().size(); ++i)
  {
    std::string name = photonjet::npbFeatureNames()[i];
    if (name == "cluster_Et") name = "cluster_et";
    if (name == "cluster_Eta") name = "cluster_eta";
    const std::string target = "pho_npb_" + name;
    const bool captured = measured && inputsCaptured && i < npb.orderedInputs.size();
    const double value = captured ? npb.orderedInputs[i] : missingValue;
    d.put(target, std::isfinite(value) ? value : missingValue);
    d.put(target + "_state", captured ? (std::isfinite(value) ? Valid : Invalid) : Missing);
    d.put(target + "_reason_index",
          d.code(target + "_reason_index", !measured              ? "no feature donor"
                                           : !captured            ? "native input not captured"
                                           : std::isfinite(value) ? ""
                                                                  : "nonfinite native NPB input"));
    d.put(target + "_origin_index",
          measured ? d.code(target + "_origin_index", captured ? "NATIVE_PRODUCER" : "UNAVAILABLE") : -1);
    d.put(target + "_evidence_index",
          captured ? d.code(target + "_evidence_index", config.provenance.producerSourceSha256) : -1);
  }
}

// The seven shower views of one photon (stride 7), in the configured order.
void fillShowerViews(photonjet::DirectEventWriter& d, const photonjet::Config& config,
                     const photonjet::EventRecord& event, const photonjet::PhotonRecord& r)
{
  for (const auto& name : config.showerDefinitions)
  {
    const auto it = std::find_if(r.showerShapes.begin(), r.showerShapes.end(),
                                 [&](const auto& v) { return v.definitionName == name; });
    if (it == r.showerShapes.end()) throw std::runtime_error("missing captured shower definition: " + name);
    for (const auto& featureName : config.output.modelFeatures)
      d.put("shv_features", photonjet::production::ModelFeature(featureName, r, event, &*it));
    const auto& view = *it;
    d.put("shv_finite_feature_state", photonjet::production::ViewShapesFinite(view));
    d.put("shv_cog_eta", view.centerOfGravityEta);
    d.put("shv_cog_phi", view.centerOfGravityPhi);
    d.put("shv_raw_center_eta", view.rawCenterEta);
    d.put("shv_raw_center_phi", view.rawCenterPhi);
    d.put("shv_center_eta_index", view.centerEtaIndex);
    d.put("shv_center_phi_index", view.centerPhiIndex);
    d.put("shv_e11", view.e11);
    d.put("shv_e33", view.e33);
    d.put("shv_e32", view.e32);
    d.put("shv_e35", view.e35);
    d.put("shv_e11_over_e33", view.e11OverE33);
    d.put("shv_e32_over_e35", view.e32OverE35);
    d.put("shv_weta_cogx", view.wetaCogX);
    d.put("shv_wphi_cogx", view.wphiCogX);
    d.put("shv_weta33_cogx", view.weta33CogX);
    d.put("shv_wphi33_cogx", view.wphi33CogX);
    d.put("shv_moment_eta_numerator", view.momentSums[0]);
    d.put("shv_moment_phi_numerator", view.momentSums[1]);
    d.put("shv_moment_denominator", view.momentSums[2]);
    d.put("shv_moment33_eta_numerator", view.momentSums[3]);
    d.put("shv_moment33_phi_numerator", view.momentSums[4]);
    d.put("shv_moment33_denominator", view.momentSums[5]);
    d.put("shv_energy_source", view.energySource);
    d.put("shv_rectangular_membership", view.rectangularMembership);
    d.put("shv_moment_membership", view.momentMembership);
    d.put("shv_good_cell_count", view.goodCellCount);
    d.put("shv_owned_cell_count", view.ownedCellCount);
    d.put("shv_active_sum_cell_count", view.activeSumCellCount);
    d.put("shv_active_moment_cell_count", view.activeMomentCellCount);
    d.put("shv_et1", view.et1);
    d.put("shv_et2", view.et2);
    d.put("shv_et3", view.et3);
    d.put("shv_et4", view.et4);
    d.put("shv_exact_zero_count", view.zeroCellCount);
    d.put("shv_negative_count", view.negativeCellCount);
    d.put("shv_nonfinite_count", view.nonFiniteCellCount);
  }
}

// One photon candidate: PPG12 eligibility, NPB, isolation witnesses, timing,
// kinematics and nominal labels, dominant truth contributor, shower shapes,
// preselection, the producer score and the shower views.
void fillPhoton(photonjet::DirectEventWriter& d, const photonjet::Config& config, const photonjet::EventRecord& event,
                const photonjet::PhotonRecord& r)
{
  fillNpbBranches(d, config, r);
  d.put("pho_ppg12_tower_mask_state", r.ppg12TowerMaskState);
  d.put("pho_ppg12_source_eligible", r.ppg12TowerMaskState == 0 ? 1 : 0);
  // Reconstructed isolation witnesses, one per retained radius in order.
  if (r.isolationWitnesses.size() != config.photon.isolationRadii.size())
    throw std::runtime_error("isolation witnesses do not cover the retained radii");
  for (const auto& w : r.isolationWitnesses)
  {
    d.put("isow_subtraction_method", w.subtractionMethod);
    d.put("isow_reconstructed_or_truth", w.reconstructedOrTruth);
    d.put("isow_cone_sum", w.coneSum);
    d.put("isow_threshold", w.threshold);
    d.put("isow_sideband_threshold", w.sidebandThreshold);
    d.put("isow_pass_state", w.passState);
    d.put("isow_constituent_count", w.constituentCount);
  }
  d.put("pho_time_raw", r.timing.meanTimeSamples);
  d.put("pho_time_ns", r.timing.timeNs);
  d.put("pho_time_ntowers", r.timing.contributingTowerCount);
  d.put("pho_time_den_e", r.timing.energyDenominator);
  d.put("pho_time_raw_state", r.timing.valid ? 1 : 2);
  d.put("pho_time_raw_reason", d.code("pho_time_raw_reason", r.timing.valid ? "" : "native cluster timing is invalid"));
  d.put("pho_time_den_e_state", std::isfinite(r.timing.energyDenominator) ? 1 : 2);
  d.put("pho_time_den_e_reason",
        d.code("pho_time_den_e_reason",
               std::isfinite(r.timing.energyDenominator) ? "" : "native cluster timing denominator is nonfinite"));
  d.put("pho_key", r.nativeClusterKey);
  d.put("pho_encounter", r.encounterOrdinal);
  d.put("pho_rank_keys", r.et);
  d.put("pho_rank_keys", -std::abs(r.eta));
  d.put("pho_rank_keys", r.encounterOrdinal);
  d.put("pho_et", r.et);
  d.put("pho_eta", r.eta);
  d.put("pho_bdt_score", r.nominal.score);
  d.put("pho_bdt_state", r.nominal.state);
  d.put("pho_id_class", r.nominal.idClass);
  d.put("pho_iso_class", r.nominal.isolationClass);
  d.put("pho_phi", r.phi);
  d.put("pho_producer_vtx_z", r.producerVertexZ);
  d.put("pho_ntime_den", r.timing.energyDenominator);
  d.put("pho_ntime_num", r.timing.energyWeightedNumerator);
  d.put("pho_ntime_ntowers", r.timing.contributingTowerCount);
  d.put("pho_ntime_finite", r.timing.finite);
  d.put("pho_ntime_version", 1);
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
  for (const auto& item :
       std::vector<std::pair<std::string, std::pair<int, int>>>{{"e11_over_e33", {10, 11}}, {"e32_over_e35", {12, 13}}})
  {
    const double a = r.primaryShapes[item.second.first], b = r.primaryShapes[item.second.second];
    d.put("pho_" + item.first, std::isfinite(a) && std::isfinite(b) && b > 0 ? a / b : missingValue);
  }
  d.put("pho_presel_reference", photonjet::production::ReferencePreselection(r));
  d.put("pho_presel_active", photonjet::production::ActivePreselection(config.photon, r));
  bool allFinite = true;
  for (const auto& name : config.output.modelFeatures)
  {
    const auto value = photonjet::production::ModelFeature(name, r, event, nullptr);
    d.put("pho_prodbdt_inputs", value);
    allFinite &= std::isfinite(value);
  }
  // The producer score with its historical witnesses. Applicability:
  // 0 validated domain, 1 diagnostic extrapolation, 2 no producer score
  // configured, 3 a nonfinite model input. Working points are recorded
  // only inside the validated domain.
  {
    const auto& production = config.photon.productionScore;
    const double raw = r.productionScore;
    const double centrality = event.analysisCentrality;
    int applicability = 2;
    if (!production.scoreName.empty())
    {
      const bool validated =
          r.et >= production.validatedEtMin && r.et < production.validatedEtMax &&
          (std::isnan(production.validatedCentralityMin) ||
           (centrality >= production.validatedCentralityMin && centrality < production.validatedCentralityMax));
      applicability = !allFinite ? 3 : validated ? 0 : 1;
    }
    std::array<double, 3> wp{missingValue, missingValue, missingValue};
    const auto& edges = production.workingPointEtEdges;
    if (applicability == 0 && production.centralityLinearWorkingPoints)
      for (std::size_t k = 0; k < wp.size(); ++k)
        wp[k] = production.workingPointIntercepts[k] + production.workingPointSlopes[k] * centrality;
    else if (applicability == 0)
      for (std::size_t bin = 0; bin + 1 < edges.size(); ++bin)
        if (r.et >= edges[bin] && r.et < edges[bin + 1])
          for (std::size_t k = 0; k < wp.size(); ++k) wp[k] = production.workingPoints[k][bin];
    d.put("pho_prodbdt_score", raw);
    d.put("pho_prodbdt_finite", std::isfinite(raw));
    d.put("pho_prodbdt_applicability", applicability);
    const std::array<const char*, 3> names{{"70", "80", "90"}};
    for (std::size_t k = 0; k < wp.size(); ++k)
    {
      d.put(std::string("pho_prodbdt_wp") + names[k], wp[k]);
      d.put(std::string("pho_prodbdt_delta_wp") + names[k], raw - wp[k]);
    }
  }
  d.put("pho_finite_feature_state", allFinite);
  d.put("pho_below15_state", r.et < 15.);
  d.put("pho_presel_bits", 7ULL | (static_cast<std::uint64_t>(allFinite) << 3));
  fillShowerViews(d, config, event, r);
}
}  // namespace

// Map one event's records to schema fields after reset. The source-entry check below
// rejects gaps and duplicates; collection counts are derived from the retained records.
void PhotonJetTree::fillEventOutput()
{
  auto& d = *m_directOutput;
  d.begin();
  const auto& r = m_event;
  if (r.sourceEntry != d.nextEntry++) throw std::runtime_error("source entry gap or duplicate");
  // The established producer recorded its own event energy and MBD charge sums
  // only in Au+Au, and only for an event its classification accepted (its
  // rejections returned earlier). In pp, and for every other event, the same
  // sums are carried by calo_* and mbd_charge_sum below.
  const bool nativeEventSums = m_config.system == CollisionSystem::AuAu && r.nodesReady && r.terminalStatus == 0;

  // --------------------------------------------------------------------------
  // Event identity and bookkeeping (evt_)
  // --------------------------------------------------------------------------
  d.put("evt_record", 0);
  d.put("evt_run", r.run);
  // The producer's own 1-based event count.
  d.put("evt_seq", r.producerEventOrdinal + 1);
  d.put("evt_phys_seq", r.physicalEventSequence);
  d.put("evt_phys_seq_valid", r.physicalEventSequenceValid);
  d.put("evt_src_file_ordinal", r.sourceFileOrdinal);
  d.put("evt_src_entry", r.sourceEntry);
  d.put("evt_src_global_entry", r.sourceGlobalEntry);
  d.put("evt_collab_iface_version", 2);
  d.put("evt_n_candidates", m_photons.size());
  d.put("evt_n_recoils", std::count_if(m_photonJetPairs.begin(), m_photonJetPairs.end(),
                                       [](const auto& pair) { return pair.recoilWitness; }));
  // The historical node check failed exactly when the reconstructed vertex
  // the analysis requires was unavailable.
  d.put("evt_nodes_ready", r.nodesReady ? 1 : 0);
  // The established event classification (EventClassificationConfig). It is
  // recorded only: every encounter keeps its row and objects, and readers
  // select triggers by exact name and events downstream.
  d.put("evt_terminal_status", r.terminalStatus);
  // The established productions never applied the embedded minimum-bias
  // classifier (-1); the MinimumBiasInfo decision itself is cent_mb_decision.
  d.put("evt_embedded_mb_decision", -1);
  // Every event value is this execution's own capture (evidence entry 0).
  const int captureEvidence = d.code("evt_aug_donor_evidence", m_config.provenance.producerSourceSha256);
  d.put("evt_aug_donor_evidence", captureEvidence);
  // Historical tag count: candidates whose producer score exceeds its WP80
  // witness (recorded only inside the validated domain).
  {
    const auto& production = m_config.photon.productionScore;
    const auto& edges = production.workingPointEtEdges;
    int tags = 0;
    for (const auto& photon : m_photons)
    {
      if (production.scoreName.empty() || !std::isfinite(photon.productionScore)) continue;
      const double centrality = m_event.analysisCentrality;
      if (!(photon.et >= production.validatedEtMin && photon.et < production.validatedEtMax &&
            (std::isnan(production.validatedCentralityMin) ||
             (centrality >= production.validatedCentralityMin && centrality < production.validatedCentralityMax))))
        continue;
      bool finiteInputs = true;
      for (const auto& name : m_config.output.modelFeatures)
        finiteInputs &= std::isfinite(photonjet::production::ModelFeature(name, photon, m_event, nullptr));
      if (!finiteInputs) continue;
      for (std::size_t bin = 0; bin + 1 < edges.size(); ++bin)
        if (photon.et >= edges[bin] && photon.et < edges[bin + 1] &&
            photon.productionScore > production.workingPoints[1][bin])
          ++tags;
    }
    d.put("evt_n_tags", tags);
  }

  // --------------------------------------------------------------------------
  // Trigger and GL1 scalers (trig_, gl1_)
  // --------------------------------------------------------------------------
  d.put("trig_live", r.triggerLiveBits);
  d.put("trig_scaled", r.triggerScaledBits);
  d.put("trig_packet_status", r.triggerPacketStatus);
  d.put("trig_packet_version", r.triggerPacketVersion);
  d.put("trig_decision_available", r.triggerDecisionsAvailable);
  d.put("trig_capture_version", 1);
  d.put("trig_decision_state", r.triggerDecisionState);
  // The snapshot is a cumulative hardware counter, never an event weight.
  // GL1 exists only in DATA; simulation keeps the schema's not-populated value.
  const bool gl1 = m_config.dataKind == DataKind::Data;
  d.put("gl1_bco", gl1 ? r.triggerBunchNumber : 0);
  const TriggerScalerSnapshotRecord* scaler = nullptr;
  if (r.triggerScalerSnapshotId >= 0)
  {
    if (static_cast<std::size_t>(r.triggerScalerSnapshotId) >= m_triggerScalerSnapshots.size())
      throw std::runtime_error("dangling trigger scaler snapshot");
    scaler = &m_triggerScalerSnapshots[r.triggerScalerSnapshotId];
  }
  d.put("gl1_state", !gl1 ? 0 : scaler ? (r.triggerPacketStatus ? 4 : 0) : r.triggerDecisionState);
  d.put("gl1_discontinuity", gl1 ? r.triggerScalerDiscontinuityFlags : 0);
  const std::array<std::uint64_t, 64> noScaler{};
  d.all("gl1_raw", scaler ? scaler->raw : noScaler);
  d.all("gl1_live", scaler ? scaler->live : noScaler);
  d.all("gl1_scaled", scaler ? scaler->scaled : noScaler);

  // --------------------------------------------------------------------------
  // Vertex, MBD and centrality (vtx_, mbd_, cent_)
  // --------------------------------------------------------------------------
  d.put("vtx_z", r.recoVertexZ);
  d.put("vtx_valid", r.recoVertexValid);
  d.put("mbd_pmt_available", static_cast<std::int32_t>(r.mbdPmtAvailable));
  d.all("mbd_pmt_q", r.mbdPmtCharge);
  // PMT arms and validity are not stored per channel. The file records the
  // one rule that derives them; a channel that follows none is refused.
  if (!r.mbdPmtId.empty())
  {
    if (r.mbdPmtArm.size() != r.mbdPmtId.size() || r.mbdPmtValid.size() != r.mbdPmtId.size())
      throw std::runtime_error("MBD PMT witnesses are not aligned");
    bool divide = true, unknown = true, valid = true;
    for (std::size_t i = 0; i < r.mbdPmtId.size(); ++i)
    {
      divide &= r.mbdPmtArm[i] == r.mbdPmtId[i] / 64;
      unknown &= r.mbdPmtArm[i] == -1;
      valid &= r.mbdPmtValid[i] == 1;
    }
    if (!divide && !unknown) throw std::runtime_error("MBD PMT arms follow no recordable rule");
    const std::string arm = divide ? "id_div_64" : "minus_one";
    if (d.policies.contains("mbd_pmt_arm") && d.policies["mbd_pmt_arm"] != arm)
      throw std::runtime_error("MBD PMT arm rule changed within the source");
    d.policies["mbd_pmt_arm"] = arm;
    d.policies["mbd_pmt_charge_valid"] = valid && d.policies.value("mbd_pmt_charge_valid", 1) == 1 ? 1 : 0;
  }
  d.put("mbd_native_available", r.mbdAvailable);
  d.put("mbd_native_valid", r.mbdNativeValid);
  d.put("mbd_npmt_south", r.mbdSouthNpmt);
  d.put("mbd_npmt_north", r.mbdNorthNpmt);
  d.put("mbd_native_north_npmt", r.mbdNorthNpmt);
  d.put("mbd_t0", r.mbdT0Ns);
  d.put("mbd_native_south_time_ns", r.mbdSouthTimeNs);
  d.put("mbd_native_north_time_ns", r.mbdNorthTimeNs);
  d.put("mbd_native_capture_version", 1);
  d.put("mbd_native_t0_finite", std::isfinite(r.mbdT0Ns));
  d.put("mbd_native_south_time_finite", std::isfinite(r.mbdSouthTimeNs));
  d.put("mbd_native_north_time_finite", std::isfinite(r.mbdNorthTimeNs));
  // Direct observations use the five-valued measurement states of the shared
  // contract: 0 missing, 1 valid, 2 invalid, 3 not evaluated, 4 not applicable.
  {
    // MbdOut: t0 is valid when MbdOut says so; an arm time needs hits and a
    // finite time above the -9999 sentinel. A missing node is missing data.
    const std::array<std::string, 3> names{{"mbd_t0", "mbd_time_south", "mbd_time_north"}};
    const std::array<double, 3> values{{r.mbdT0Ns, r.mbdSouthTimeNs, r.mbdNorthTimeNs}};
    const std::array<int, 2> hits{{r.mbdSouthNpmt, r.mbdNorthNpmt}};
    for (std::size_t i = 0; i < names.size(); ++i)
    {
      const double value = values[i];
      const bool usable = std::isfinite(value) && value > -9999.;
      const bool valid = i == 0 ? r.mbdNativeValid == 1 && usable : hits[i - 1] > 0 && usable;
      if (i) d.put(names[i], r.mbdAvailable && std::isfinite(value) ? value : missingValue);
      d.put(names[i] + "_state", !r.mbdAvailable ? 0 : valid ? 1 : 2);
      d.put(names[i] + "_reason", d.code(names[i] + "_reason", !r.mbdAvailable ? "native MbdOut node unavailable"
                                                               : valid         ? ""
                                                               : i == 0        ? "native MbdOut invalid"
                                                                        : "empty MBD arm or invalid native arm time"));
    }
  }
  // Event charge: the producer's own total where it recorded one, otherwise
  // the exact finite-positive sum over a complete, valid 128-PMT array.
  const double nativeCharge = nativeEventSums ? r.mbdTotalCharge : missingValue;
  const int chargeOrigin = std::isfinite(nativeCharge) ? 1 : r.mbdPmtArrayComplete ? 2 : 0;
  d.put("mbd_total_charge_native", nativeCharge);
  d.put("mbd_charge_sum", chargeOrigin == 1 ? nativeCharge : chargeOrigin == 2 ? r.mbdTotalCharge : missingValue);
  d.put("mbd_charge_sum_state", chargeOrigin ? 1 : 0);
  d.put("mbd_charge_sum_origin", chargeOrigin);
  if (r.centralityPmtCharge != r.mbdPmtCharge && !r.centralityPmtCharge.empty())
  {
    // NaNs compare by their stored bits; different populations need different
    // storage and cannot be silently dropped as redundant.
    if (r.centralityPmtCharge.size() != r.mbdPmtCharge.size() ||
        std::memcmp(r.centralityPmtCharge.data(), r.mbdPmtCharge.data(), r.mbdPmtCharge.size() * sizeof(double)))
      throw std::runtime_error("centrality and collaborator PMT charge witnesses differ");
  }
  if (r.centralityPmtId != r.mbdPmtId && !r.centralityPmtId.empty())
    throw std::runtime_error("centrality and collaborator PMT identities differ");
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
  // pp has no native centrality decision: -1, not "invalid".
  d.put("cent_native_valid",
        m_config.system == CollisionSystem::AuAu ? static_cast<std::int32_t>(r.centralityNativeValid) : -1);
  d.put("cent_selected_charge", r.centralitySelectedCharge);
  d.put("cent_mb_decision", static_cast<std::int32_t>(r.minimumBiasDecision));
  // Au+Au without a centrality for this event is undefined (NaN); p+p has none (-1).
  d.put("cent_value", m_config.system == CollisionSystem::AuAu && r.analysisCentrality < 0
                          ? missingValue
                          : static_cast<double>(r.analysisCentrality));

  // --------------------------------------------------------------------------
  // Calorimeter sums (calo_)
  // --------------------------------------------------------------------------
  d.put("calo_emcal_native", nativeEventSums ? r.cemcEnergy : missingValue);
  d.put("calo_ihcal_native", nativeEventSums ? r.ihcalEnergy : missingValue);
  d.put("calo_ohcal_native", nativeEventSums ? r.ohcalEnergy : missingValue);
  d.put("calo_total_native", nativeEventSums ? r.totalCaloEnergy : missingValue);
  d.put("calo_require_isgood_native", nativeEventSums ? static_cast<std::int32_t>(r.caloRequiredGoodTowers) : -1);
  d.put("calo_max_cluster_e", r.maximumRawClusterEnergy);
  // Bits 0-2: the EMCal, inner and outer HCal sums; bit 3: their total.
  const std::uint32_t totalBit = (r.caloAvailableMask & 7U) == 7U ? 8U : 0U;
  d.put("calo_native_mask", nativeEventSums ? r.caloAvailableMask | totalBit : 0U);
  d.put("calo_capture_version", 1);
  d.put("calo_available_mask", r.caloAvailableMask);
  d.put("calo_valid_mask", (r.caloValidMask >> 8) & 7);
  d.put("calo_require_isgood", r.caloRequiredGoodTowers);
  const std::array<double, 4> energy = {r.cemcEnergy, r.ihcalEnergy, r.ohcalEnergy, r.totalCaloEnergy};
  const std::array<std::string, 4> names = {"calo_emcal", "calo_ihcal", "calo_ohcal", "calo_total"};
  for (std::size_t i = 0; i < energy.size(); ++i)
  {
    const unsigned mask = i == 3 ? 7U : 1U << i;
    const bool present = (r.caloAvailableMask & mask) == mask, valid = ((r.caloValidMask >> 8) & mask) == mask;
    const double value = present ? static_cast<float>(energy[i]) : missingValue;
    d.put(names[i], value);
    d.put(names[i] + "_state", !present ? 0 : valid && std::isfinite(value) ? 1 : 2);
    d.put(names[i] + "_reason",
          d.code(names[i] + "_reason", !present ? "calibrated calorimeter input unavailable"
                                       : valid && std::isfinite(value)
                                           ? ""
                                           : "calibrated calorimeter input incomplete or nonfinite"));
  }

  // --------------------------------------------------------------------------
  // Underlying event and jet views (ue_, jetview_)
  // --------------------------------------------------------------------------
  d.put("ue_flow_mode", r.jetBackgroundFlowMode);
  d.put("ue_v2", r.jetBackgroundV2);
  d.put("ue_psi2", r.jetBackgroundPsi2);
  d.all("ue_emcal", r.jetBackgroundUEEmcal);
  d.all("ue_ihcal", r.jetBackgroundUEIhcal);
  d.all("ue_ohcal", r.jetBackgroundUEOhcal);
  d.put("ue_n_strips", r.jetBackgroundNStrips);
  d.put("ue_n_towers", r.jetBackgroundNTowers);
  d.put("ue_flow_failure_state", static_cast<std::int32_t>(r.jetBackgroundFlowFailure));
  // The jet views this event recorded, in collection order.
  std::size_t jetViews = 0;
  if (r.recoObjectVertexInDomain)
    for (bool unsub : {false, true})
      for (int rc : {2, 3, 4})
        for (std::size_t i = 0; i < r.recoJetView.size(); ++i)
          if ((r.recoJetView[i] == static_cast<int>(JetView::AuAuNoSub)) == unsub && r.recoJetRadiusCode[i] == rc * 10)
          {
            d.put("jetview_available", r.recoJetCaptureState[i] == static_cast<int>(CaptureState::Complete) ? 1 : 0);
            ++jetViews;
          }
  d.jetViewsRecorded = std::max(d.jetViewsRecorded, jetViews);
  d.put("ue_contract_state", static_cast<std::int32_t>(r.jetBackgroundState));
  if (r.jetBackgroundUEEmcal.size() != r.jetBackgroundUEIhcal.size() ||
      r.jetBackgroundUEEmcal.size() != r.jetBackgroundUEOhcal.size())
    throw std::runtime_error("background calorimeter eta-bin counts disagree");
  d.put("ue_eta_bins", r.jetBackgroundUEEmcal.size());

  // --------------------------------------------------------------------------
  // PPG12 event quantities and weights (ppg12_, w_)
  // --------------------------------------------------------------------------
  // The lane component of an applied PPG12 weight: photon+jet SI 1, DI 2,
  // inclusive-jet SI 3, DI 4; 0 when no PPG12 weight applies.
  const auto& w = m_config.weights;
  d.put("ppg12_weight_active_state", r.ppg12WeightActive ? 1 : 0);
  d.put("ppg12_weight_lane_code", r.ppg12WeightActive ? photonjet::production::PPG12WeightLaneCode(w) : 0);
  d.put("ppg12_arb_state", r.ppg12ArbitrationState);
  d.put("ppg12_lead_key", r.ppg12LeadKey);
  d.put("ppg12_lead_encounter", r.ppg12LeadEncounter);
  d.put("ppg12_lead_et", r.ppg12LeadEt);
  d.put("ppg12_lead_eta", r.ppg12LeadEta);
  d.put("ppg12_lead_phi", r.ppg12LeadPhi);
  d.put("w_event", r.eventWeight);
  d.put("w_legacy", r.eventWeight);
  // Saved collaborator interface convention: SIM sample/vertex normalization
  // remains unverified independently of a finite composed producer weight.
  d.put("w_sample_valid", m_config.dataKind == DataKind::Data ? 1 : -1);
  d.put("w_vertex_valid", m_config.dataKind == DataKind::Data ? 1 : -1);
  const WeightComponentRecord* weights = nullptr;
  for (const auto& value : m_weightComponents)
    if (value.componentType == "final")
    {
      if (weights) throw std::runtime_error("duplicate final event weight");
      weights = &value;
    }
  if (!weights) throw std::runtime_error("missing final event weight");
  d.put("w_slice", weights->sliceWeight);
  d.put("w_xsec", weights->crossSectionWeight);
  d.put("w_vertex", weights->vertexWeight);
  d.put("w_sidi", weights->siDiWeight);
  d.put("w_period", weights->periodWeight);
  d.put("w_exposure", weights->exposureWeight);
  d.put("w_application_count", weights->applicationCount);

  // --------------------------------------------------------------------------
  // Capture bookkeeping and the truth event (recopho_, truth_, truthvtx_, truthcen_, truthjetcap_)
  // --------------------------------------------------------------------------
  d.put("recopho_capture_version", r.recoPhotonCaptureVersion);
  d.put("recopho_capture_domain", r.recoPhotonCaptureDomain);
  d.put("recopho_inspected", r.recoPhotonInspectedCount);
  d.put("recopho_eligible", r.recoPhotonEligibleCount);
  d.put("recopho_written", r.recoPhotonWrittenCount);
  d.put("recopho_unclassifiable", r.recoPhotonUnclassifiableCount);
  d.put("recopho_capture_state", r.recoPhotonCaptureVersion == 0                      ? 0
                                 : r.recoPhotonCaptureState == CaptureState::Complete ? 1
                                                                                      : 2);
  d.put("truth_vtx_z", r.truthHardVertexZ);
  d.put("truth_mbvtx_z", r.truthMinimumBiasVertexZ);
  d.put("truth_denominator_complete",
        m_config.dataKind == DataKind::Data ? 1 : static_cast<std::int32_t>(r.truthDenominatorComplete));
  d.put("truthcen_native_embedded_primary_pid22_count", r.truthPhotonEmbeddedPrimaryCount);
  d.put("truthcen_serialized_raw_count", r.truthPhotonWrittenCount);
  d.put("truthcen_rejected_count", r.truthPhotonRejectedCount);
  d.put("truthcen_rejected_invalid_kinematics_count", r.truthPhotonRejectedKinematicsCount);
  d.put("truthcen_rejected_nonfinite_isolation_count", r.truthPhotonRejectedIsolationCount);
  d.put("truthcen_duplicate_track_count", r.truthPhotonDuplicateTrackCount);
  d.put("truthcen_isolation_input_incomplete_count", r.truthPhotonIsolationIncompleteCount);
  d.put("truthcen_analysis_signal_count", r.truthPhotonAnalysisSignalCount);
  d.put("truthvtx_primary_id", r.truthPrimaryVertexId);
  d.put("truthvtx_capture_state", m_config.dataKind == DataKind::Data                   ? -1
                                  : r.truthVertexCaptureState == CaptureState::Complete ? 1
                                                                                        : 0);
  d.put("truthcen_capture_state", m_config.dataKind == DataKind::Data                   ? -1
                                  : r.truthPhotonCaptureState == CaptureState::Complete ? 1
                                                                                        : 0);
  d.put("truthjetcap_state", m_config.dataKind == DataKind::Data                     ? -1
                             : r.truthJetCaptureState == CaptureState::NotApplicable ? 0
                                                                                     : 1);
  d.all("truthjetcap_container_valid", r.truthJetContainerValid);
  d.truthRadiiRecorded = std::max(d.truthRadiiRecorded, r.truthJetContainerValid.size());

  // --------------------------------------------------------------------------
  // Au+Au embedded diagnostics (embdiag_, truthmiss_)
  // --------------------------------------------------------------------------
  // Au+Au embedded diagnostics in emission order: the stitch decision (the
  // sample ownership witness TreeToHists reads) and the truth-isolation audit.
  for (const auto& row : r.embeddedDiagnostics)
  {
    d.put("embdiag_kind", d.code("embdiag_kind", row.kind));
    d.put("embdiag_audit_bin", row.auditBin);
    d.put("embdiag_sample_code", row.sampleCode);
    d.put("embdiag_centrality_bin_index", row.centralityBinIndex);
    d.put("embdiag_photon_class", row.photonClass);
    d.put("embdiag_truth_pt", row.truthPt);
    d.put("embdiag_truth_isolation_witness", row.truthIsolationWitness);
    d.put("embdiag_stitch_decision", row.stitchDecision);
    d.put("embdiag_stitch_photon_pt", row.stitchPhotonPt);
    d.meta["embdiag.trigger"] = "SIM";
  }
  // Au+Au embedded truth-photon misses: one trigger, the canonical isolation
  // view, and the decision paths of the file's rows, '|'-joined in sorted
  // order (a single path, as the established producer kept, when only one occurs).
  for (const auto& miss : r.truthMisses)
  {
    d.put("truthmiss_centrality_bin", miss.centralityBin);
    d.put("truthmiss_truth_pt", miss.truthPt);
    d.meta["truthPhotonMisses.trigger"] = "SIM";
    d.meta["truthPhotonMisses.iso_view"] = "canonical";
    d.truthMissPaths.insert(miss.decisionPath);
    std::string paths;
    for (const std::string& path : d.truthMissPaths) paths += (paths.empty() ? "" : "|") + path;
    d.meta["truthPhotonMisses.decision_path"] = paths;
  }

  // References between objects are positions in their event-local collections.
  const auto pkeys = indices(m_photons, &PhotonRecord::photonId);
  const auto tkeys = indices(m_truthPhotons, &TruthPhotonRecord::truthPhotonId);

  // --------------------------------------------------------------------------
  // Photons (pho_, shv_, isow_)
  // --------------------------------------------------------------------------
  // Photon records and ordered shower-definition witnesses.
  d.photonRows += m_photons.size();
  for (const auto& photon : m_photons) fillPhoton(d, m_config, m_event, photon);

  // --------------------------------------------------------------------------
  // Jets (jet<RR>_, jetunsub<RR>_) and truth jets (truthjet<RR>_)
  // --------------------------------------------------------------------------
  // Jet references address the canonical collection order, regardless of the
  // order in which reconstruction registered the nodes.
  std::map<Key, int> jkeys, tjkeys;
  int ji = 0, ti = 0;
  for (const JetRecord* jet : photonjet::production::CanonicalJetOrder(m_jets))
  {
    const JetRecord& r = *jet;
    const bool unsub = r.view == JetView::AuAuNoSub;
    const std::string prefix = (unsub ? "jetunsub0" : "jet0") + std::to_string(std::lround(r.radius * 10.));
    if (!jkeys.emplace(key(r.jetId), ji++).second) throw std::runtime_error("duplicate jet identity");
    if (unsub && !r.calibrationValid) ++d.invalidUnsubtractedJes[prefix];
    d.jetCollectionsWithJets.insert(prefix);
    d.put(prefix + "_rawpt", r.rawPt);
    d.put(prefix + "_pt", r.correctedPt);
    d.put(prefix + "_eta", r.eta);
    d.put(prefix + "_phi", r.phi);
    d.put(prefix + "_area", r.area);
    if (!unsub)
    {
      d.put(prefix + "_jes_valid", r.calibrationValid);
      d.put(prefix + "_emcal_frac", r.calorimeterFractions[0]);
      d.put(prefix + "_ihcal_frac", r.calorimeterFractions[1]);
      d.put(prefix + "_ohcal_frac", r.calorimeterFractions[2]);
    }
  }
  if (jkeys.size() != m_jets.size()) throw std::runtime_error("unsupported jet collection");
  for (int rc : {2, 3, 4})
    for (const auto& r : m_truthJets)
      if (r.radius == rc / 10.)
      {
        if (!tjkeys.emplace(key(r.truthJetId), ti++).second) throw std::runtime_error("duplicate truth jet identity");
        const auto prefix = "truthjet0" + std::to_string(rc);
        d.put(prefix + "_pt", r.pt);
        d.put(prefix + "_eta", r.eta);
        d.put(prefix + "_phi", r.phi);
        d.put(prefix + "_guard_state", r.pt < 5 ? 1 : r.pt >= 35 ? 2 : 0);
      }
  if (tjkeys.size() != m_truthJets.size()) throw std::runtime_error("unsupported truth jet collection");
  d.truthJetRows += m_truthJets.size();

  // --------------------------------------------------------------------------
  // Photon-jet pairs (pair_)
  // --------------------------------------------------------------------------
  for (const auto& r : m_photonJetPairs)
  {
    d.put("pair_pho", local(pkeys, r.photonId));
    d.put("pair_jet", local(jkeys, r.jetId));
    d.put("pair_delta_phi", r.deltaPhi);
    d.put("pair_xjgamma", r.xJGamma);
    d.put("pair_photon_rank", r.photonRank);
    d.put("pair_jet_rank", r.jetRank);
    d.put("pair_recoil_state", r.recoilWitness);
    // Retained pair rows did not run wrong-object classification; zero is the
    // explicitly defined legacy unclassified default, not a truth statement.
    d.put("pair_wrong_photon_class", 0);
    d.put("pair_wrong_recoil_class", 0);
  }

  // --------------------------------------------------------------------------
  // Truth photons and vertices (truthpho_, truthvtx_)
  // --------------------------------------------------------------------------
  const int truthSampleRole = m_config.system != CollisionSystem::AuAu                                       ? 0
                              : photonjet::production::EmbeddedPhotonSampleCode(m_config.source.sample) != 0 ? 1
                              : photonjet::production::EmbeddedJetSampleCode(m_config.source.sample) != 0    ? 2
                                                                                                             : 0;
  for (const auto& r : m_truthPhotons)
  {
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
    d.put("truthpho_source_role", r.analysisSignal ? 1 : 0);
    d.put("truthpho_truth_isolation_witness", r.isolationR03);
    d.put("truthpho_reporting_guard_state", r.pt < 15 ? 1 : r.pt >= 35 ? 2 : 0);
    d.put("truthpho_gen_embedding_id", r.embeddingId);
    // The sample family of an Au+Au embedded source: 1 photon+jet, 2 inclusive
    // jet. The p+p productions declared no photon-ID training role (0).
    d.put("truthpho_sample_source_role", truthSampleRole);
    d.put("truthpho_prompt_class", static_cast<std::int32_t>(r.promptClass));
  }
  for (const auto& r : m_truthVertices)
  {
    d.put("truthvtx_id", r.vertexId);
    d.put("truthvtx_embedding_id", r.embeddingId);
    d.put("truthvtx_embedding_valid", r.embeddingValid);
    d.put("truthvtx_z", r.z);
    d.put("truthvtx_z_valid", r.valid);
  }

  // --------------------------------------------------------------------------
  // Isolation constituents (isoc_) and shower cells (cell_)
  // --------------------------------------------------------------------------
  for (const auto& r : m_isolationConstituents)
  {
    if (r.radius != .4) throw std::runtime_error("constituents must be stored once at R=0.4");
    d.put("isoc_pho", local(pkeys, r.photonId));
    d.put("isoc_native_source", r.subsystem);
    d.put("isoc_phosub_residual", missingValue);
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
  for (const auto& r : m_photonCells)
  {
    const int parent = local(pkeys, r.photonId);
    const auto& photon = m_photons[parent];
    auto center = [&](const std::string& name)
    {
      const auto it = std::find_if(photon.showerShapes.begin(), photon.showerShapes.end(),
                                   [&](const auto& v) { return v.definitionName == name; });
      if (it == photon.showerShapes.end() || !std::isfinite(it->rawCenterEta) || !std::isfinite(it->rawCenterPhi))
        throw std::runtime_error("cell witness lacks a finite shower center");
      return std::pair<int, int>(static_cast<int>(std::floor(it->rawCenterEta)),
                                 (static_cast<int>(std::floor(it->rawCenterPhi)) % 256 + 256) % 256);
    };
    auto wrap = [](int value) { return (value + 128 + 256) % 256 - 128; };
    const auto c = center("H70");
    const int de = r.etaIndex - c.first, raw = r.phiIndex - c.second;
    int dp = wrap(raw);
    if (dp == -128 && raw > 0) dp = 128;
    int grid = 0;
    for (const auto& item : std::vector<std::pair<std::string, int>>{{"H0", 1}, {"H70", 2}})
    {
      const auto c = center(item.first);
      if (std::abs(r.etaIndex - c.first) <= 3 && std::abs(wrap(r.phiIndex - c.second)) <= 3) grid |= item.second;
    }
    const bool member = r.good && std::isfinite(r.calibratedEnergy) && r.calibratedEnergy > 0.;
    d.put("cell_pho", parent);
    d.put("cell_ieta", r.etaIndex);
    d.put("cell_iphi", r.phiIndex);
    d.put("cell_local_ieta", de);
    d.put("cell_local_iphi", dp);
    d.put("cell_seed_state", de == 0 && dp == 0);
    d.put("cell_grid_bits", grid);
    d.put("cell_e", r.calibratedEnergy);
    d.put("cell_rawcluster_value", r.clusterMapEnergy);
    d.put("cell_is_good", r.good);
    d.put("cell_is_zero", r.zeroEnergy);
    d.put("cell_is_negative", r.negativeEnergy);
    d.put("cell_is_nonfinite", !std::isfinite(r.calibratedEnergy));
    d.put("cell_rawcluster_owned", r.ownedByCluster);
    d.put("cell_rawcluster_value_present", r.ownedByCluster);
    d.put("cell_time_samples", r.towerTimeSamples);
    d.put("cell_time_available", r.towerAvailable);
    d.put("cell_time_finite", std::isfinite(r.towerTimeSamples));
    d.put("cell_tower_status", r.towerStatus);
    d.put("cell_time_version", 1);
    d.put("cell_denominator_member", member);
    d.put("cell_floor0_member", member);
    d.put("cell_floor70_member", r.good && std::isfinite(r.calibratedEnergy) && r.calibratedEnergy > .070);
  }

  // --------------------------------------------------------------------------
  // Reconstruction-to-truth associations (pho_assoc_, pho_resp_, plink_, jlink_)
  // --------------------------------------------------------------------------
  // The identity association gives the raw association and the response
  // classes; a photon excluded from the native signal can still be a raw
  // match. The established producer's distance links are the plink_ rows.
  const bool sim = m_config.dataKind == DataKind::Simulation;
  std::vector<int> rawTruth(m_photons.size(), -1), assoc(m_photons.size(), sim ? 0 : -1);
  std::vector<int> response(m_photons.size(), -128), truthResponse(m_truthPhotons.size(), -128);
  std::vector<double> distance(m_photons.size(), missingValue);
  std::set<int> seenReco, matchedRaw;
  bool unknown = false;
  for (const auto& r : m_photonTruthLinks)
  {
    const int pi = local(pkeys, r.photonId, true), ti = local(tkeys, r.truthPhotonId, true);
    const int raw = local(tkeys, r.associationTruthPhotonId, true), state = static_cast<int>(r.state);
    if (!sim || pi < 0 && ti < 0) throw std::runtime_error("invalid photon truth relation population");
    if (pi >= 0)
    {
      if (!seenReco.insert(pi).second) throw std::runtime_error("duplicate photon truth relation");
      if ((r.associationState == 2) != (raw >= 0)) throw std::runtime_error("raw association endpoint disagrees");
      assoc[pi] = r.associationState;
      rawTruth[pi] = raw;
      response[pi] = r.associationState == 0 ? -1 : r.associationState == 2 ? 0 : 1;
      unknown |= r.associationState == 0;
      if (raw >= 0)
      {
        matchedRaw.insert(raw);
        distance[pi] = r.associationDeltaR;
      }
    }
    if (state != 0 && state != 3 && state != 4 && state != 5)
      throw std::runtime_error("unsupported photon response state");
  }
  if (sim && seenReco.size() != m_photons.size()) throw std::runtime_error("incomplete photon relation population");
  for (const auto& r : m_photonDistanceLinks)
  {
    const int pi = local(pkeys, r.photonId, true), ti = local(tkeys, r.truthPhotonId, true),
              state = static_cast<int>(r.state);
    if (!sim || pi < 0 && ti < 0 || (state != 3 && state != 4 && state != 5))
      throw std::runtime_error("invalid photon distance link");
    d.put("plink_reco", pi);
    d.put("plink_truth", ti);
    d.put("plink_class", state == 3 ? 0 : state == 4 ? 1 : 2);
    d.put("plink_metric", r.deltaR);
    d.put("plink_policy", 1);
  }
  const bool complete = m_event.recoPhotonCaptureVersion == 1 &&
                        m_event.recoPhotonCaptureState == CaptureState::Complete &&
                        m_event.recoPhotonUnclassifiableCount == 0;
  const bool outside = m_event.recoPhotonCaptureVersion == 0 && !m_event.recoVertexValid;
  for (std::size_t i = 0; i < truthResponse.size(); ++i)
    if (!matchedRaw.count(i)) truthResponse[i] = outside ? 3 : complete && !unknown ? 2 : -1;
  d.all("pho_assoc_truthpho", rawTruth);
  d.all("pho_assoc_state", assoc);
  d.all("pho_assoc_dr", distance);
  d.all("pho_resp_class", response);
  d.all("truthpho_resp_class", truthResponse);
  for (std::size_t i = 0; i < m_photons.size(); ++i)
  {
    const auto& photon = m_photons[i];
    const auto& witness = photon.dominantTruth;
    const bool evaluator = sim && witness.evaluator != DominantTruthEvaluator::Unavailable &&
                           witness.evaluator != DominantTruthEvaluator::NotApplicable;
    const bool matched = rawTruth[i] >= 0 && m_truthPhotons[rawTruth[i]].analysisSignal;
    if (m_config.system == CollisionSystem::AuAu)
    {
      // The Au+Au producer labelled only the photons its analysis loop examined
      // (an accepted event and the truth-miss photon domain), each with its
      // dominant primary's barcode, matched or not.
      const auto& c = m_config.truthMisses;
      const bool examined = sim && m_event.terminalStatus == 0 && std::fabs(photon.eta) < c.absEtaMax &&
                            photon.et >= c.etMinGeV && photon.et < c.etMaxGeV;
      const bool labelled = witness.state != DominantTruthState::EvaluatorUnavailable &&
                            witness.state != DominantTruthState::NotSimulation;
      const bool primary =
          witness.state == DominantTruthState::ValidPrimary || witness.state == DominantTruthState::InvalidPrimary;
      d.put("pho_truth_match_state", examined && labelled ? int(matched) : -1);
      d.put("pho_truth_match_barcode", examined && primary ? witness.barcode : -1);
    }
    else
    {
      d.put("pho_truth_match_state", evaluator ? int(matched) : -1);
      d.put("pho_truth_match_barcode", matched ? m_truthPhotons[rawTruth[i]].barcode : -1);
    }
  }
  for (const auto& r : m_jetTruthLinks)
  {
    const int ri = local(jkeys, r.jetId, true), ti = local(tjkeys, r.truthJetId, true),
              state = static_cast<int>(r.state);
    if (!sim || ri < 0 && ti < 0 || state == 1 || state < 0 || state > 5)
      throw std::runtime_error("invalid jet truth link");
    if (r.selectedMatch != (state == 3)) throw std::runtime_error("jet selected-match flag differs");
    d.put("jlink_reco", ri);
    d.put("jlink_truth", ti);
    d.put("jlink_metric", r.deltaR);
    d.put("jlink_class", state == 0 ? -1 : state == 2 ? 5 : state == 3 ? 0 : state == 4 ? 1 : 2);
  }

  // --------------------------------------------------------------------------
  // The established producer's PPG12 diagnostics (pdiag_, ppg12_wfill, ppg12_rfill)
  // --------------------------------------------------------------------------
  // PPG12 fill-audit occurrences, in the order the producer logged them.
  for (const auto& row : r.ppg12Diagnostics)
  {
    d.put("pdiag_trigger", d.code("pdiag_trigger", row.trigger));
    d.put("pdiag_kind", d.code("pdiag_kind", row.kind));
    d.put("pdiag_stage_code", row.stageCode);
    d.put("pdiag_source_code", row.sourceCode);
    d.put("pdiag_sample_code", row.sampleCode);
    d.put("pdiag_decision", row.decision);
    d.put("pdiag_cluster_ordinal", row.clusterOrdinal);
    d.put("pdiag_truth_track_id", row.truthTrackId);
    d.put("pdiag_generator_barcode", row.generatorBarcode);
    d.put("pdiag_source_photon_pt", row.sourcePhotonPt);
    d.put("pdiag_flow_value", row.flowValue);
    d.put("pdiag_reco_photon_pt", row.recoPhotonPt);
    d.put("pdiag_response_photon_pt", row.responsePhotonPt);
    d.put("pdiag_truth_photon_pt", row.truthPhotonPt);
    d.put("pdiag_truth_prior_weight", row.truthPriorWeight);
    d.put("pdiag_reco_vertex_z", row.recoVertexZ);
    d.put("pdiag_hard_truth_vertex_z", row.hardTruthVertexZ);
    d.put("pdiag_mb_truth_vertex_z", row.minimumBiasTruthVertexZ);
    d.put("pdiag_vertex_weight", row.vertexWeight);
    d.put("pdiag_period_event_weight", row.periodEventWeight);
    d.put("pdiag_cluster_energy", row.clusterEnergy);
    d.put("pdiag_cluster_eta", row.clusterEta);
    d.put("pdiag_cluster_et", row.clusterEt);
    d.put("pdiag_truth_energy", row.truthEnergy);
    d.put("pdiag_truth_eta", row.truthEta);
    d.put("pdiag_truth_et", row.truthEt);
    d.put("pdiag_energy_contribution", row.energyContribution);
    d.put("pdiag_occurrence_weight", row.occurrenceWeight);
    d.put("pdiag_selection_bits", row.selectionBits);
  }
  // Fill counters of the PPG12 simulation histogram flow, by weighting lane:
  // DATA and Au+Au have no PPG12 weighting lane, so their counts are zero.
  std::int64_t weighted = 0, raw = 0;
  for (int lane = 1; lane <= 4; ++lane)
  {
    weighted += r.ppg12WeightedFills[lane];
    raw += r.ppg12RawFills[lane];
    d.put("ppg12_wfill_code" + std::to_string(lane), r.ppg12WeightedFills[lane]);
    d.put("ppg12_rfill_code" + std::to_string(lane), r.ppg12RawFills[lane]);
  }
  d.put("ppg12_wfill", weighted);
  d.put("ppg12_rfill", raw);

  // --------------------------------------------------------------------------
  // Collections the current production leaves empty, and the event
  // --------------------------------------------------------------------------
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
  // The centrality-fit denominator counts the accepted events (vertex cut,
  // generator-slice stitch and centrality range passed), as the established
  // producer did, with unit counts: centrality weights would bias their fit.
  if (d.centrality && r.terminalStatus == 0 && r.recoVertexValid && std::isfinite(r.analysisCentrality) &&
      r.analysisCentrality >= 0. && r.analysisCentrality <= 100.)
    d.centrality->Fill(r.analysisCentrality);
  m_source.lastEntry = r.sourceEntry;
}

void PhotonJetTree::fillRejectedOutput(const UpstreamRejectedEventRecord& r)
{
  auto& d = *m_directOutput;
  d.begin();
  if (r.sourceEntry != d.nextEntry++) throw std::runtime_error("rejected source entry gap or duplicate");
  d.put("evt_record", 1);
  d.put("evt_run", r.run);
  d.put("evt_phys_seq", r.physicalEventSequence);
  d.put("evt_phys_seq_valid", r.physicalEventSequence >= 0);
  d.put("evt_src_entry", r.sourceEntry);
  d.put("evt_src_file_ordinal", m_source.sourceFileOrdinal);
  d.put("evt_src_global_entry", m_config.source.sourceGlobalEntryBegin >= 0
                                    ? m_config.source.sourceGlobalEntryBegin + r.sourceEntry
                                    : std::int64_t{-1});
  d.finish(true);
}

void PhotonJetTree::writeCompletionMetadata()
{
  auto& d = *m_directOutput;
  m_source.completed = !m_aborted;
  if (d.centrality) d.writer.object("hist/SIM/h_centrality", *d.centrality);
  if (!d.entries)
  {
    d.begin();
    d.finish(false, true);
  }
  json counters = json::object();
  for (const auto& f : d.specs)
    if (f.length == -1 && f.type != "string") counters[f.name] = f.counter;
  const bool auau = m_config.system == CollisionSystem::AuAu;
  // The jet collections in collection order. A collection records its algorithm
  // and subtraction with its jets: one without any keeps only its definition.
  // The view constants list the jet views in the same order, as many as the
  // most any event recorded.
  json collections = {{"reco", json::array()}, {"truth", json::array()}};
  json viewCodes = json::array(), viewInputs = json::array(), viewSubtraction = json::array();
  for (bool unsub : {false, true})
    for (int rc : {2, 3, 4})
      if (!unsub || auau)
      {
        const std::string prefix = (unsub ? "jetunsub0" : "jet0") + std::to_string(rc);
        const std::string input =
            auau ? (unsub ? "towerinfo_retower_unsubtracted" : "towerinfo_sub1") : "towerinfo_calib";
        const std::string subtraction = auau && !unsub ? "SUB1" : "none";
        json entry = {{"name", prefix}, {"radius_code", rc}, {"radius", bits(rc / 10.)}, {"input_identity", input}};
        if (d.jetCollectionsWithJets.count(prefix))
        {
          entry["algorithm"] = "antikt";
          entry["subtraction_identity"] = subtraction;
        }
        collections["reco"].push_back(entry);
        if (viewCodes.size() < d.jetViewsRecorded)
        {
          viewCodes.push_back(rc);
          viewInputs.push_back(input);
          viewSubtraction.push_back(subtraction);
        }
      }
  for (int rc : {2, 3, 4})
    collections["truth"].push_back(
        {{"name", "truthjet0" + std::to_string(rc)}, {"radius_code", rc}, {"radius", bits(rc / 10.)}});
  // File constants shared by every row: shower views, isolation witness radii,
  // collection vocabularies, weight targets and the producer score identity.
  // A jet whose energy-scale correction is undefined keeps its original pT
  // (jet<RR>_jes_valid = 0), and its photon pairs carry xJ = NaN.
  if (!m_config.jets.energyScalePayloadSha256.empty())
    d.policies["invalid_jes_policy"] = "PRESERVE_ORIGINAL_JET_CORRECTED_PT_SET_PAIR_XJGAMMA_NAN";
  // The unsubtracted Au+Au views carry no jes_valid branch: their invalid-JES
  // jets are counted per view instead (zero included).
  if (auau)
  {
    json invalid = json::object();
    for (const char* view : {"jetunsub02", "jetunsub03", "jetunsub04"})
      invalid[view] = d.invalidUnsubtractedJes.count(view) ? d.invalidUnsubtractedJes.at(view) : 0;
    d.policies["invalid_jes_unsubtracted_jets"] = invalid;
  }
  // File constants describe the rows the file holds: a file without photons has
  // no shower-view, isolation-witness or model constants, and one without truth
  // jets no truth-jet constants, as the established output records it.
  json constants = d.meta.is_null() ? json::object() : d.meta;
  {
    if (d.photonRows)
    {
      json names = json::array(), digests = json::array(), floors = json::array(), byFloor = json::array();
      for (const auto& name : m_config.showerDefinitions)
      {
        const auto& view = photonjet::showerViewDefinition(name);
        names.push_back(view.name);
        digests.push_back(view.semanticSha256);
        floors.push_back(bits(view.floorGeV));
      }
      // Candidate view list in the historical order: 70 MeV views, then 0 MeV.
      for (bool floor70 : {true, false})
        for (const auto& name : m_config.showerDefinitions)
          if ((photonjet::showerViewDefinition(name).floorGeV > 0) == floor70) byFloor.push_back(name);
      constants["showerViews.definition_name"] = names;
      constants["showerViews.semantic_sha256"] = digests;
      constants["showerViews.floor_gev"] = floors;
      constants["photons.shower_definition_views"] = byFloor;
      json radii = json::array();
      for (const double radius : m_config.photon.isolationRadii) radii.push_back(bits(radius));
      constants["isolationWitnesses.radius"] = radii;
    }
    constants["view_radius_codes"] = viewCodes;
    constants["view_input_identity"] = viewInputs;
    constants["view_subtraction_identity"] = viewSubtraction;
    if (m_source.retainedEvents)
    {
      constants["weights.component_type"] = "event";
      constants["weights.target_type"] = 1;
    }
    // Radius codes of the truth-jet collections, as many as the most any event
    // recorded (DATA records none).
    json truthRadii = json::array();
    for (int rc : {2, 3, 4})
      if (truthRadii.size() < d.truthRadiiRecorded) truthRadii.push_back(rc);
    constants["truth_radius_codes"] = truthRadii;
    if (d.truthJetRows)
    {
      constants["truthJets.algorithm"] = "antikt";
      constants["truthJets.ownership_state"] = "source_owned";
    }
    const auto& production = m_config.photon.productionScore;
    if (!production.scoreName.empty() && d.photonRows)
    {
      constants["models.model_sha256"] = production.modelSha256;
      constants["models.shower_definition_id"] = production.showerDefinition;
      constants["models.shower_semantic_sha256"] =
          photonjet::showerViewDefinition(production.showerDefinition).semanticSha256;
    }
  }
  // A stride is the observed length of a fixed-length per-photon vector; a file
  // without photons observes none (0), as the eventTree contract records it.
  const std::uint64_t observed = d.photonRows ? 1 : 0, features = m_config.output.modelFeatures.size();
  const json strides = {{"shv_features", observed * features},
                        {"pho_prodbdt_inputs", observed * features},
                        {"pho_rank_keys", observed * 3},
                        {"pho_npbcap_inputs", observed * 25},
                        {"pho_npbcap_input_finite", observed * 25},
                        {"pho_aunpbcap_inputs", observed * 25},
                        {"pho_aunpbcap_input_finite", observed * 25}};
  json meta = {
      {"schema",
       {{"schema", "RJEventTreeV1"},
        {"version", "1.0"},
        {"schema_sha256", m_config.output.schemaSha256},
        {"entries", d.entries},
        {"retained_events", m_source.retainedEvents},
        {"system", auau ? "auau" : "pp"},
        {"lane", m_config.source.lane},
        {"layout", "native_direct"},
        {"encoding", "counted"},
        {"counters", counters},
        {"absent_branches", d.missing},
        {"strides", strides},
        {"filled_reserved", d.filledReserved}}},
      {"source",
       {{"lane", m_config.source.lane},
        {"dataset", m_source.dataset},
        {"sample", m_source.sample},
        {"period", m_source.period},
        {"si_di_role", m_source.siDiRole},
        {"run", m_source.run},
        {"segment", m_source.segment},
        {"source_file_ordinal", m_source.sourceFileOrdinal},
        {"first_entry", m_source.firstEntry},
        {"last_entry", m_source.lastEntry},
        {"input_uri_sha256", m_source.inputUriSha256},
        {"input_file_sha256", m_source.inputFileSha256},
        {"source_manifest_sha256", m_source.sourceManifestSha256},
        {"completed", m_source.completed},
        {"retained_events", m_source.retainedEvents},
        {"upstream_rejected_events", m_source.upstreamRejectedEvents}}},
      {"collections", collections},
      {"codes", d.codes},
      {"constants", constants},
      {"policies", d.policies},
      {"weighting",
       {{"schema", "ProducerWeightV1"},
        {"composition", m_config.weights.ppg12Nominal ? "slice*vertex*sidi*period"
                                                      : (m_config.weights.applyVertexReweight ? "vertex" : "unity")},
        {"vertex_enabled", m_config.weights.applyVertexReweight},
        {"vertex_coordinate", m_config.weights.applyVertexReweight ? "truth_hard" : "none"},
        {"vertex_lookup", m_config.weights.applyVertexReweight ? "TH1::Interpolate" : "none"},
        {"double_interaction", m_config.weights.doubleInteraction},
        {"vertex_file_sha256", m_config.weights.vertexReweightFileSha256},
        {"vertex_histogram", m_config.weights.vertexReweightHistogram},
        {"absolute_normalization_applied", false},
        {"centrality_applied", false},
        {"aliases",
         m_config.weights.ppg12Nominal ? json{{"w_xsec", "w_slice"}, {"w_exposure", "w_period"}} : json::object()}}},
      {"native_capture",
       {{"status", d.missing.empty() ? "MAPPING_COMPLETE" : "INCOMPLETE_MAPPING"},
        {"missing_branches", d.missing},
        {"producer_source_sha256", m_config.provenance.producerSourceSha256},
        {"producer_library_sha256", m_config.provenance.producerLibrarySha256}}},
      {"photon_id", nominalPhotonIdReceipt(m_config.photon.nominalId, m_nominalReceipt)},
      {"isolation", nominalIsolationReceipt(m_config.photon.nominalIsolation, m_nominalReceipt)}};
  // DATA run trigger configuration: one entry per GL1 bit. Readers select a
  // stream by its exact name and enabled prescale, never by bit index.
  if (!m_triggerRunInfo.empty())
  {
    if (m_triggerRunInfo.size() != 64) throw std::runtime_error("trigger run information does not hold 64 bits");
    json table = {{"name", json::array()},      {"state", json::array()},   {"initial_prescale", json::array()},
                  {"prescale", json::array()},  {"run_raw", json::array()}, {"run_live", json::array()},
                  {"run_scaled", json::array()}};
    for (std::size_t bit = 0; bit < m_triggerRunInfo.size(); ++bit)
    {
      const auto& row = m_triggerRunInfo[bit];
      if (row.bit != static_cast<int>(bit)) throw std::runtime_error("trigger run information is out of bit order");
      table["name"].push_back(row.name);
      table["state"].push_back(row.state);
      table["initial_prescale"].push_back(std::isfinite(row.hardwarePrescale) ? static_cast<int>(row.hardwarePrescale)
                                                                              : -1);
      table["prescale"].push_back(bits(row.runAveragePrescale));
      table["run_raw"].push_back(row.rawCount);
      table["run_live"].push_back(row.liveCount);
      table["run_scaled"].push_back(row.scaledCount);
    }
    meta["trigger_run_info"] = table;
  }
  for (auto it = meta.begin(); it != meta.end(); ++it) d.writer.text("meta/" + it.key(), it.value().dump());
  std::ostringstream text;
  appendLine(text, "completion_status", m_aborted ? "aborted" : "complete");
  appendLine(text, "encountered_events", m_source.encounteredEvents);
  appendLine(text, "retained_events", m_source.retainedEvents);
  appendLine(text, "upstream_rejected_events", m_source.upstreamRejectedEvents);
  appendLine(text, "upstream_observed_events", m_upstreamObservedEvents);
  appendLine(text, "first_entry", m_source.firstEntry);
  appendLine(text, "last_entry", m_source.lastEntry);
  d.writer.text("completion", text.str());
  photonjet::production::WriteProgress(m_config, m_aborted ? "aborted" : "complete",
                                       static_cast<std::int64_t>(m_upstreamObservedEvents), -1,
                                       m_aborted ? "event loop aborted" : "", true);
}

void PhotonJetTree::closeOutput()
{
  if (!m_directOutput) return;
  m_directOutput->writer.close();
  m_directOutput.reset();
}

// ============================================================================
// The filled branches
// ============================================================================

namespace
{
// Every schema branch with a filling rule above, in name order. The list is
// fixed, independent of what an event holds: an empty collection still has
// its rule. A branch added to config/eventTree_v1.yaml must be filled above
// and listed here.
std::set<std::string> filledBranches()
{
  return {"calo_available_mask",
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
          "embdiag_audit_bin",
          "embdiag_centrality_bin_index",
          "embdiag_kind",
          "embdiag_photon_class",
          "embdiag_sample_code",
          "embdiag_stitch_decision",
          "embdiag_stitch_photon_pt",
          "embdiag_truth_isolation_witness",
          "embdiag_truth_pt",
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
          "pdiag_cluster_energy",
          "pdiag_cluster_et",
          "pdiag_cluster_eta",
          "pdiag_cluster_ordinal",
          "pdiag_decision",
          "pdiag_energy_contribution",
          "pdiag_flow_value",
          "pdiag_generator_barcode",
          "pdiag_hard_truth_vertex_z",
          "pdiag_kind",
          "pdiag_mb_truth_vertex_z",
          "pdiag_occurrence_weight",
          "pdiag_period_event_weight",
          "pdiag_reco_photon_pt",
          "pdiag_reco_vertex_z",
          "pdiag_response_photon_pt",
          "pdiag_sample_code",
          "pdiag_selection_bits",
          "pdiag_source_code",
          "pdiag_source_photon_pt",
          "pdiag_stage_code",
          "pdiag_trigger",
          "pdiag_truth_energy",
          "pdiag_truth_et",
          "pdiag_truth_eta",
          "pdiag_truth_photon_pt",
          "pdiag_truth_prior_weight",
          "pdiag_truth_track_id",
          "pdiag_vertex_weight",
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
          "pho_bdt_score",
          "pho_bdt_state",
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
          "pho_id_class",
          "pho_iso_class",
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
          "pho_ppg12_source_eligible",
          "pho_ppg12_tower_mask_state",
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
          "ppg12_arb_state",
          "ppg12_lead_encounter",
          "ppg12_lead_et",
          "ppg12_lead_eta",
          "ppg12_lead_key",
          "ppg12_lead_phi",
          "ppg12_rfill",
          "ppg12_rfill_code1",
          "ppg12_rfill_code2",
          "ppg12_rfill_code3",
          "ppg12_rfill_code4",
          "ppg12_weight_active_state",
          "ppg12_weight_lane_code",
          "ppg12_wfill",
          "ppg12_wfill_code1",
          "ppg12_wfill_code2",
          "ppg12_wfill_code3",
          "ppg12_wfill_code4",
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
          "truthmiss_centrality_bin",
          "truthmiss_truth_pt",
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
          "w_xsec"};
}
}  // namespace
