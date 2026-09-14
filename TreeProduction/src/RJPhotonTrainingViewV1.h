#ifndef RJ_PHOTON_TRAINING_VIEW_V1_H
#define RJ_PHOTON_TRAINING_VIEW_V1_H

// Opt-in, normalized photon-ID training rows for the seven shower
// definitions.  This is deliberately a separate ROOT artifact: it does not
// add a tree to ReplayFoundationV1 and does not alter the legacy
// AuAuPhotonIDTrainingTree contract.

#include "RJReplayFoundationV1.h"
#include "RJShowerFactorialV1.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace RJPhotonTrainingViewV1
{
using RJReplayFoundationV1::Identity128;
using RJReplayFoundationV1::ModelApplicability;
using RJReplayFoundationV1::PhotonCandidateRow;
using RJReplayFoundationV1::SerializedIdentity128;
using RJReplayFoundationV1::SerializedInt64;
using RJReplayFoundationV1::SerializedUInt64;
using RJReplayFoundationV1::ShowerFeatureViewRow;
using RJReplayFoundationV1::SourceOccurrenceRow;
using RJReplayFoundationV1::makeIdentity;
using RJReplayFoundationV1::sha256Hex;

constexpr const char* kSchemaName = "RJ_PHOTON_TRAINING_VIEW_V1";
constexpr int kSchemaVersion = 1;
constexpr const char* kTreeName = "RJPhotonTrainingViewV1";
constexpr int kSystemPP = 1;
constexpr int kSystemAuAu = 2;

inline const std::vector<std::string>& featureNames(int systemCode)
{
  static const std::vector<std::string> pp{
      "cluster_Et", "cluster_weta_cogx", "cluster_wphi_cogx", "vertexz",
      "cluster_Eta", "e11_over_e33", "cluster_et1", "cluster_et2",
      "cluster_et3", "cluster_et4", "e32_over_e35"};
  static const std::vector<std::string> auau{
      "cluster_Et", "cluster_weta_cogx", "cluster_wphi_cogx",
      "cluster_weta33_cogx", "cluster_wphi33_cogx", "vertexz",
      "cluster_Eta", "e11_over_e33", "cluster_et1", "cluster_et2",
      "cluster_et3", "cluster_et4", "e32_over_e35", "centrality"};
  if (systemCode == kSystemPP) return pp;
  if (systemCode == kSystemAuAu) return auau;
  throw std::invalid_argument("unknown photon-training system code");
}

inline std::string featureContractText(int systemCode)
{
  std::ostringstream out;
  out << "RJ_PHOTON_TRAINING_FEATURE_ORDER_V1|system="
      << (systemCode == kSystemPP ? "pp" : "auau") << '|';
  const auto& names = featureNames(systemCode);
  for (std::size_t index = 0; index < names.size(); ++index)
  {
    if (index != 0) out << ',';
    out << names[index];
  }
  return out.str();
}

inline std::string featureContractSha256(int systemCode)
{
  return sha256Hex(featureContractText(systemCode));
}

inline Identity128 sourceOccurrenceIdentity(const SourceOccurrenceRow& source)
{
  return makeIdentity(
      source.lane + "|" + source.dataset + "|" + source.sample + "|" +
      source.period + "|" + std::to_string(source.run) + "|" +
      std::to_string(source.segment) + "|" + source.input_uri_hash + "|" +
      source.input_file_sha256 + "|" + source.source_manifest_sha256);
}

inline Identity128 eventIdentity(const SourceOccurrenceRow& source,
                                 int run,
                                 std::int64_t eventSequence)
{
  return makeIdentity(
      source.lane + "|" + source.sample + "|" + std::to_string(run) + "|" +
      std::to_string(source.segment) + "|" + std::to_string(eventSequence));
}

struct Label
{
  int cluster_index = -1;
  int training_label = -1;
  int is_signal = -1;
  std::string label_authority;

  int truth_match_found = -1;
  int truth_hepmc_association_valid = -1;
  int truth_photon_class = -1;
  int truth_is_prompt = -1;
  double truth_iso_et = std::numeric_limits<double>::quiet_NaN();
  double truth_iso_et_r03 = std::numeric_limits<double>::quiet_NaN();
  double truth_iso_et_r04 = std::numeric_limits<double>::quiet_NaN();
  int truth_iso_valid = -1;
  int truth_iso_pass = -1;
  int cluster_truth_track_id = -1;
  int cluster_truth_pid = 0;
  int cluster_truth_barcode = -1;
  double truth_energy_contribution = std::numeric_limits<double>::quiet_NaN();

  int source_role = 0;
  int source_sample_code = 0;
  int ppg12_source_role_label = -1;
  int npb_label = -1;
  int is_npb = -1;
  int minimum_bias_classifier_decision = -1;

  int ppg12_analysis_window_pass = -1;
  int ppg12_response_window_pass = -1;
  int ppg12_truth_window_pass_r04 = -1;
  int ppg12_sample_bin = 0;
  double ppg12_xsec_pb = std::numeric_limits<double>::quiet_NaN();
  double ppg12_xsec_weight = 1.0;
  double ppg12_window_low = std::numeric_limits<double>::quiet_NaN();
  double ppg12_window_high = std::numeric_limits<double>::quiet_NaN();
  double max_truth_jet_pt_r04 = std::numeric_limits<double>::quiet_NaN();

  double weight_slice = 1.0;
  double weight_cross_section = 1.0;
  double weight_vertex = 1.0;
  double weight_si_di = 1.0;
  double weight_period = 1.0;
  double weight_exposure = 1.0;
  double weight_final = 1.0;
  int weight_application_count = 1;
};

struct CandidateContext
{
  Identity128 event_id;
  Identity128 candidate_id;
  // Compact native IDs are shard-local. The producer must bind the context
  // to its initialized source occurrence; a compact SOURCE ID is insufficient.
  // Legacy callers may omit this field, retaining their existing contract.
  Identity128 source_occurrence_id;
  int run = 0;
  std::int64_t event_sequence = 0;
  int encounter_ordinal = -1;
  std::uint64_t cluster_map_key = 0;
  double cluster_et = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double vertex_z = std::numeric_limits<double>::quiet_NaN();
  double centrality = -1.0;
  double event_weight = 1.0;
};

struct Row
{
  Identity128 id, source_id, event_id, candidate_id, definition_id;
  std::string definition_name, shower_semantic_sha256, feature_contract_sha256;
  std::vector<float> ordered_features;
  int feature_count = 0;
  int finite_feature_state = 0;
  int system_code = 0;

  std::string source_lane, source_dataset, source_sample, source_period;
  std::string source_si_di_role, source_ownership_state;
  std::string source_manifest_sha256, input_uri_sha256, input_file_sha256;
  int source_run = 0, source_segment = 0;

  int run = 0;
  std::int64_t event_sequence = 0;
  int encounter_ordinal = -1;
  std::uint64_t cluster_map_key = 0;
  int cluster_index = -1;
  double cluster_et = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double vertex_z = std::numeric_limits<double>::quiet_NaN();
  double centrality = -1.0;
  double event_weight = 1.0;
  int model_domain_state = static_cast<int>(ModelApplicability::MODEL_NOT_APPLICABLE);
  int below15_retention_state = 0;
  int nominal_training_eligible = 0;
  int working_point_state = -1;
  int tag_state = -1;

  int training_label = -1;
  int is_signal = -1;
  std::string label_authority;
  int truth_match_found = -1;
  int truth_hepmc_association_valid = -1;
  int truth_photon_class = -1;
  int truth_is_prompt = -1;
  double truth_iso_et = std::numeric_limits<double>::quiet_NaN();
  double truth_iso_et_r03 = std::numeric_limits<double>::quiet_NaN();
  double truth_iso_et_r04 = std::numeric_limits<double>::quiet_NaN();
  int truth_iso_valid = -1;
  int truth_iso_pass = -1;
  int cluster_truth_track_id = -1;
  int cluster_truth_pid = 0;
  int cluster_truth_barcode = -1;
  double truth_energy_contribution = std::numeric_limits<double>::quiet_NaN();
  int source_role = 0;
  int source_sample_code = 0;
  int ppg12_source_role_label = -1;
  int npb_label = -1;
  int is_npb = -1;
  int minimum_bias_classifier_decision = -1;
  int ppg12_analysis_window_pass = -1;
  int ppg12_response_window_pass = -1;
  int ppg12_truth_window_pass_r04 = -1;
  int ppg12_sample_bin = 0;
  double ppg12_xsec_pb = std::numeric_limits<double>::quiet_NaN();
  double ppg12_xsec_weight = 1.0;
  double ppg12_window_low = std::numeric_limits<double>::quiet_NaN();
  double ppg12_window_high = std::numeric_limits<double>::quiet_NaN();
  double max_truth_jet_pt_r04 = std::numeric_limits<double>::quiet_NaN();
  double weight_slice = 1.0;
  double weight_cross_section = 1.0;
  double weight_vertex = 1.0;
  double weight_si_di = 1.0;
  double weight_period = 1.0;
  double weight_exposure = 1.0;
  double weight_final = 1.0;
  int weight_application_count = 1;
};

inline std::vector<float> canonicalFeatures(const ShowerFeatureViewRow& view,
                                            const CandidateContext& candidate,
                                            int systemCode)
{
  return RJShowerFactorialV1::modelFeatures(
      view,candidate.cluster_et,candidate.vertex_z,candidate.eta,
      candidate.centrality,systemCode==kSystemAuAu);
}

class Runtime
{
 public:
  Runtime() = default;
  Runtime(const Runtime&) = delete;
  Runtime& operator=(const Runtime&) = delete;
  ~Runtime()
  {
    if (m_file && m_file->IsOpen()) m_file->Close();
  }

  bool initialize(const std::string& outputPath,
                  int systemCode,
                  const SourceOccurrenceRow& source,
                  const std::string& configSha256,
                  const std::string& codeSha256,
                  std::string* error = nullptr)
  {
    if (m_initialized || m_failed)
      return fail(error, "training-view writer cannot be reinitialized");
    if (outputPath.empty()) return fail(error, "training-view output path is empty");
    if (systemCode != kSystemPP && systemCode != kSystemAuAu)
      return fail(error, "training-view system code must be pp or AuAu");
    if (source.lane.empty() || source.dataset.empty() || source.sample.empty() ||
        source.period.empty() || source.si_di_role.empty() ||
        source.ownership_state.empty() || source.run < 0 || source.segment < 0)
      return fail(error, "training-view source contract is incomplete");
    const Identity128 stableSourceId=sourceOccurrenceIdentity(source);
    const Identity128 compactSourceId=RJReplayFoundationV1::makeScopedIdentity(
        RJReplayFoundationV1::CompactIdentityDomain::SOURCE,0,0);
    if (source.id.isNull() ||
        (source.id != stableSourceId && source.id != compactSourceId) ||
        !validHex64(source.source_manifest_sha256) ||
        !validHex64(source.input_uri_hash) || !validHex64(source.input_file_sha256))
      return fail(error, "training-view source stable identity or input/manifest hash is invalid");
    if (!validHex64(configSha256) || !validHex64(codeSha256))
      return fail(error, "training-view config/code hashes are invalid");

    TDirectory* savedDirectory = gDirectory;
    m_file.reset(TFile::Open(outputPath.c_str(), "RECREATE"));
    if (!m_file || !m_file->IsOpen() || m_file->IsZombie())
    {
      if (savedDirectory) savedDirectory->cd();
      return fail(error, "could not create training-view ROOT file");
    }
    m_file->SetCompressionAlgorithm(static_cast<int>(ROOT::RCompressionSetting::EAlgorithm::kZSTD));
    m_file->SetCompressionLevel(5);
    m_systemCode = systemCode;
    m_source = source;
    m_compactInputs = source.id == compactSourceId;
    // Keep the persisted sidecar contract unchanged. Compact IDs are input
    // references only; persisted joins remain qualified by the stable source.
    m_source.id = stableSourceId;
    m_configSha256 = configSha256;
    m_codeSha256 = codeSha256;
    bookTree();
    if (!m_tree)
    {
      if (savedDirectory) savedDirectory->cd();
      return fail(error, "could not book training-view tree");
    }

    TNamed schema("rj_photon_training_schema", kSchemaName);
    schema.Write("rj_photon_training_schema", TObject::kOverwrite);
    TNamed version("rj_photon_training_schema_version", "1");
    version.Write("rj_photon_training_schema_version", TObject::kOverwrite);
    writeMeta("schema_sha256", schemaSha256());
    writeMeta("pp_feature_contract_sha256", featureContractSha256(kSystemPP));
    writeMeta("auau_feature_contract_sha256", featureContractSha256(kSystemAuAu));
    writeMeta("source_manifest_sha256", source.source_manifest_sha256);
    writeMeta("config_sha256", configSha256);
    writeMeta("code_sha256", codeSha256);
    writeMeta("input_identity_mode",m_compactInputs?"scoped_compact":"legacy_hash");
    writeMeta("persisted_identity_mode","legacy_hash_source_qualified");
    m_initialized = true;
    if (savedDirectory) savedDirectory->cd();
    return true;
  }

  bool recordLabel(std::int64_t eventSequence,
                   int encounterOrdinal,
                   const Label& label,
                   std::string* error = nullptr)
  {
    if (!ready(error)) return false;
    if (eventSequence < 0 || encounterOrdinal < 0 || label.cluster_index < 0)
      return fail(error, "invalid training-label event or candidate ordinal");
    if (!bindEvent(eventSequence,error)) return false;
    if (label.cluster_index != encounterOrdinal)
      return fail(error, "training-label cluster index does not match encounter ordinal");
    if (label.training_label != 0 && label.training_label != 1 &&
        !(m_systemCode == kSystemAuAu && label.training_label == -1))
      return fail(error, "training label is outside the declared source-role states");
    if (label.is_signal != 0 && label.is_signal != 1)
      return fail(error, "legacy signal label is not binary");
    if (label.label_authority.empty() || label.weight_application_count != 1)
      return fail(error, "training label authority or weight application count is invalid");
    const int expectedSampleCode = sourceSampleCode();
    if (expectedSampleCode > 0 && label.source_sample_code != expectedSampleCode)
      return fail(error, "training label sample code does not match the frozen source identity");
    const double weightComponents[] = {
        label.weight_slice, label.weight_cross_section, label.weight_vertex,
        label.weight_si_di, label.weight_period, label.weight_exposure,
        label.weight_final};
    if (!std::all_of(std::begin(weightComponents), std::end(weightComponents),
                     [](double value) { return std::isfinite(value); }))
      return fail(error, "training label weight ledger contains a nonfinite value");
    // The accepted p+p replay contract retains cross_section as the
    // provenance alias of slice and exposure as the provenance alias of
    // period.  Each physical factor is nevertheless applied exactly once.
    // Au+Au has independent vertex and exposure (centrality) components.
    double componentProduct = 1.0;
    if (m_systemCode == kSystemPP)
    {
      const double aliasScale = std::max(
          {1.0, std::fabs(label.weight_slice),
           std::fabs(label.weight_cross_section),
           std::fabs(label.weight_period), std::fabs(label.weight_exposure)});
      if (std::fabs(label.weight_slice - label.weight_cross_section) >
              1.0e-12 * aliasScale ||
          std::fabs(label.weight_period - label.weight_exposure) >
              1.0e-12 * aliasScale)
        return fail(error, "p+p training weight provenance aliases do not close");
      componentProduct = label.weight_slice * label.weight_vertex *
          label.weight_si_di * label.weight_period;
    }
    else
    {
      componentProduct = label.weight_slice * label.weight_cross_section *
          label.weight_vertex * label.weight_si_di * label.weight_period *
          label.weight_exposure;
    }
    const double weightScale = std::max(
        {1.0, std::fabs(componentProduct), std::fabs(label.weight_final)});
    if (std::fabs(componentProduct - label.weight_final) > 1.0e-12 * weightScale)
      return fail(error, "training label weight components do not close to final weight");
    const Key key{eventSequence, encounterOrdinal};
    if (m_seenCandidates.count(key) != 0)
      return fail(error, "training label arrived after its candidate");
    if (!m_labels.emplace(key, label).second)
      return fail(error, "duplicate training label for event/candidate ordinal");
    return true;
  }

  bool appendCandidate(const CandidateContext& candidate,
                       const std::vector<ShowerFeatureViewRow>& views,
                       std::string* error = nullptr)
  {
    if (!ready(error)) return false;
    if (candidate.event_sequence < 0 || candidate.encounter_ordinal < 0)
      return fail(error, "invalid training candidate event or encounter ordinal");
    if (!bindEvent(candidate.event_sequence,error)) return false;
    const Key key{candidate.event_sequence, candidate.encounter_ordinal};
    const auto labelIt = m_labels.find(key);
    if (m_seenCandidates.count(key) != 0 ||
        m_seenClusterKeys.count(candidate.cluster_map_key) != 0)
      return fail(error, "duplicate training candidate encounter or cluster key");
    if ((m_compactInputs && candidate.source_occurrence_id != m_source.id) ||
        (!candidate.source_occurrence_id.isNull() &&
         candidate.source_occurrence_id != m_source.id))
      return fail(error, "training candidate source occurrence mismatch");
    if (candidate.event_id.isNull() || candidate.candidate_id.isNull())
      return fail(error, "training candidate event/candidate identity is null");
    const Identity128 stableEventId=eventIdentity(
        m_source,candidate.run,candidate.event_sequence);
    const Identity128 expectedEventId=m_compactInputs
        ? RJReplayFoundationV1::makeScopedIdentity(
            RJReplayFoundationV1::CompactIdentityDomain::EVENT,0,
            static_cast<std::uint64_t>(candidate.event_sequence))
        : stableEventId;
    if (candidate.run != m_source.run ||
        candidate.event_id != expectedEventId)
      return fail(error, "training candidate event stable identity input mismatch");
    const Identity128 stableCandidateId = makeIdentity(
        stableEventId.hex() + "|candidate|" +
        std::to_string(candidate.encounter_ordinal) + "|" +
        std::to_string(candidate.cluster_map_key));
    if (m_compactInputs && static_cast<std::uint64_t>(candidate.event_sequence) >= (1ULL<<56U))
      return fail(error, "training candidate compact event scope is out of range");
    const Identity128 expectedCandidateId=m_compactInputs
        ? RJReplayFoundationV1::makeScopedIdentity(
            RJReplayFoundationV1::CompactIdentityDomain::CANDIDATE,
            static_cast<std::uint64_t>(candidate.event_sequence),m_nextCompactCandidate)
        : stableCandidateId;
    if (candidate.candidate_id != expectedCandidateId)
      return fail(error, "training candidate stable identity input mismatch");
    // Validate identity even for unlabeled replay candidates. Those candidates
    // advance the compact ordinal just as they do in the foundation producer.
    m_seenCandidates.insert(key);
    m_seenClusterKeys.insert(candidate.cluster_map_key);
    ++m_nextCompactCandidate;
    if (labelIt == m_labels.end()) return true;
    if (!std::isfinite(candidate.event_weight))
      return fail(error, "training candidate event weight is nonfinite");
    const double eventWeightScale = std::max(
        {1.0, std::fabs(candidate.event_weight),
         std::fabs(labelIt->second.weight_final)});
    if (std::fabs(candidate.event_weight - labelIt->second.weight_final) >
        1.0e-12 * eventWeightScale)
      return fail(error, "training candidate event weight does not match final weight ledger");
    if (views.size() != 7)
      return fail(error, "labeled training candidate does not have seven shower views");

    std::set<std::string> observed;
    for (const auto& view : views)
    {
      if (!observed.insert(view.definition_name).second)
        return fail(error, "duplicate shower definition for labeled training candidate");
    }
    static const std::set<std::string> expected{"H70", "H0", "G70", "G0", "O70", "O0", "R70"};
    if (observed != expected)
      return fail(error, "labeled training candidate shower-definition set is incomplete");

    for (const auto& view : views)
    {
      const auto& definition=RJShowerFactorialV1::definition(view.definition_name);
      const auto expectedDefinitionId=makeIdentity(
          std::string("shower-definition|")+definition.name+"|"+
          RJShowerFactorialV1::semanticText(definition));
      if(view.candidate_id!=candidate.candidate_id||
         view.definition_id!=expectedDefinitionId||
         view.semantic_sha256!=RJShowerFactorialV1::semanticSha256(definition))
        return fail(error,"training-view candidate/definition semantic identity mismatch");
      Row row;
      row.source_id = m_source.id;
      row.event_id = stableEventId;
      row.candidate_id = stableCandidateId;
      row.definition_id = view.definition_id;
      row.definition_name = view.definition_name;
      row.shower_semantic_sha256 = view.semantic_sha256;
      row.feature_contract_sha256 = featureContractSha256(m_systemCode);
      row.id = makeIdentity(stableCandidateId.hex() + "|training-view|" + view.definition_id.hex());
      row.ordered_features = canonicalFeatures(view, candidate, m_systemCode);
      row.feature_count = static_cast<int>(row.ordered_features.size());
      if (row.ordered_features.size() != featureNames(m_systemCode).size())
        return fail(error, "canonical feature helper disagrees with the feature contract");
      // Model-input validity is system-specific.  In particular, p+p owns an
      // 11-feature contract and must not become INPUT_INVALID merely because
      // Au+Au-only 3x3 width diagnostics are absent.  Au+Au includes those
      // widths in its canonical vector, so the same check remains strict.
      const bool finite = std::all_of(
          row.ordered_features.begin(), row.ordered_features.end(),
          [](float value) { return std::isfinite(value); });
      row.finite_feature_state = finite ? 1 : 0;
      row.system_code = m_systemCode;

      row.source_lane = m_source.lane;
      row.source_dataset = m_source.dataset;
      row.source_sample = m_source.sample;
      row.source_period = m_source.period;
      row.source_si_di_role = m_source.si_di_role;
      row.source_ownership_state = m_source.ownership_state;
      row.source_manifest_sha256 = m_source.source_manifest_sha256;
      row.input_uri_sha256 = m_source.input_uri_hash;
      row.input_file_sha256 = m_source.input_file_sha256;
      row.source_run = m_source.run;
      row.source_segment = m_source.segment;

      row.run = candidate.run;
      row.event_sequence = candidate.event_sequence;
      row.encounter_ordinal = candidate.encounter_ordinal;
      row.cluster_map_key = candidate.cluster_map_key;
      row.cluster_index = labelIt->second.cluster_index;
      row.cluster_et = candidate.cluster_et;
      row.eta = candidate.eta;
      row.phi = candidate.phi;
      row.vertex_z = candidate.vertex_z;
      row.centrality = candidate.centrality;
      row.event_weight = candidate.event_weight;
      row.below15_retention_state = candidate.cluster_et < 15.0 ? 1 : 0;
      const bool centralityValid = m_systemCode == kSystemPP ||
          (candidate.centrality >= 0.0 && candidate.centrality < 80.0);
      row.model_domain_state = !finite
          ? static_cast<int>(ModelApplicability::INPUT_INVALID)
          : (candidate.cluster_et >= 15.0 && candidate.cluster_et < 35.0 && centralityValid
             ? static_cast<int>(ModelApplicability::VALIDATED_DOMAIN)
             : static_cast<int>(ModelApplicability::DIAGNOSTIC_EXTRAPOLATION));
      // A training-view artifact never owns a working-point or tag decision.
      // In particular, below-15 rows remain diagnostic with explicit null state.
      row.working_point_state = -1;
      row.tag_state = -1;
      assignLabel(row, labelIt->second);
      row.nominal_training_eligible =
          row.model_domain_state == static_cast<int>(ModelApplicability::VALIDATED_DOMAIN) &&
          (row.training_label == 0 || row.training_label == 1) ? 1 : 0;
      m_pendingRows.push_back(std::move(row));
    }
    m_consumed.insert(key);
    return true;
  }

  bool finishEvent(std::int64_t eventSequence, std::string* error = nullptr)
  {
    if (!ready(error)) return false;
    if (!bindEvent(eventSequence,error)) return false;
    for (const auto& item : m_labels)
      if (item.first.first != eventSequence || m_consumed.count(item.first) == 0)
        return fail(error, "accepted legacy training label has no replay candidate/view rows");
    for (const auto& row : m_pendingRows)
    {
      if (row.event_sequence != eventSequence)
        return fail(error, "training-view event buffer contains a cross-event row");
    }
    // Validate the complete event before any ROOT Fill. I/O failures still
    // poison the artifact and can never publish a completion certificate.
    for (const auto& row : m_pendingRows)
    {
      m_row = row;
      syncSerializedRow();
      if (m_tree->Fill() <= 0) return fail(error, "training-view tree fill failed");
      ++m_entries;
    }
    eraseEvent(eventSequence);
    m_pendingRows.clear();
    m_seenCandidates.clear();
    m_seenClusterKeys.clear();
    m_nextCompactCandidate=0;
    m_lastFinishedEvent=eventSequence;
    m_activeEvent=-1;
    return true;
  }

  bool finish(std::string* error = nullptr)
  {
    if (!ready(error)) return false;
    if (m_activeEvent >= 0 || !m_labels.empty() || !m_pendingRows.empty())
      return fail(error, "training-view writer has unfinished event state");
    TDirectory* savedDirectory = gDirectory;
    m_file->cd();
    if (!m_tree || m_tree->Write("", TObject::kOverwrite) <= 0)
    {
      if (savedDirectory) savedDirectory->cd();
      return fail(error, "training-view tree write failed");
    }
    TNamed complete("rj_photon_training_complete", "1");
    complete.Write("rj_photon_training_complete", TObject::kOverwrite);
    const std::string entryCount = std::to_string(m_entries);
    TNamed entries("rj_photon_training_entries", entryCount.c_str());
    entries.Write("rj_photon_training_entries", TObject::kOverwrite);
    m_file->Write("", TObject::kOverwrite);
    m_file->Close();
    m_finished = true;
    if (savedDirectory) savedDirectory->cd();
    return true;
  }

  std::int64_t entries() const { return m_entries; }
  Identity128 sourceOccurrenceId() const { return m_source.id; }

  int sourceSampleCode() const
  {
    const std::string& sample = m_source.sample;
    std::size_t begin = sample.size();
    while (begin > 0 && sample[begin - 1] >= '0' && sample[begin - 1] <= '9')
      --begin;
    if (begin == sample.size()) return 0;
    int code = 0;
    for (std::size_t index = begin; index < sample.size(); ++index)
      code = code * 10 + static_cast<int>(sample[index] - '0');
    return code;
  }

  static std::string schemaSha256()
  {
    return sha256Hex(
        "RJ_PHOTON_TRAINING_VIEW_V1|tree=RJPhotonTrainingViewV1|"
        "identity=source,event,candidate,definition|"
        "source_stable_inputs=lane,dataset,sample,period,run,segment,input_uri,input_file,manifest|"
        "event_stable_inputs=lane,sample,run,segment,event_sequence|"
        "candidate_stable_inputs=event,encounter_ordinal,cluster_map_key|"
        "features=ordered+contract|"
        "labels=truth+source+npb|weights=component-ledger|domain=15to35");
  }

 private:
  using Key = std::pair<std::int64_t, int>;

  static bool validHex64(const std::string& value)
  {
    if (value.size() != 64) return false;
    return std::all_of(value.begin(), value.end(), [](char c)
    {
      return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    });
  }

  bool fail(std::string* error, const std::string& message)
  {
    m_failed = true;
    m_pendingRows.clear();
    if (error) *error = message;
    return false;
  }

  bool ready(std::string* error)
  {
    if (!m_initialized || m_finished || m_failed)
    {
      if (error) *error = m_failed ? "training-view writer is failed" : "training-view writer is not active";
      return false;
    }
    return true;
  }

  bool bindEvent(std::int64_t eventSequence,std::string* error)
  {
    if (eventSequence < 0 || eventSequence <= m_lastFinishedEvent)
      return fail(error,"training-view event is invalid or already committed");
    if (m_activeEvent >= 0 && m_activeEvent != eventSequence)
      return fail(error,"training-view cross-event transaction");
    m_activeEvent=eventSequence;
    return true;
  }

  void eraseEvent(std::int64_t eventSequence)
  {
    for (auto it = m_labels.begin(); it != m_labels.end();)
      it = it->first.first == eventSequence ? m_labels.erase(it) : std::next(it);
    for (auto it = m_consumed.begin(); it != m_consumed.end();)
      it = it->first == eventSequence ? m_consumed.erase(it) : std::next(it);
  }

  static void assignLabel(Row& row, const Label& label)
  {
    row.training_label = label.training_label;
    row.is_signal = label.is_signal;
    row.label_authority = label.label_authority;
    row.truth_match_found = label.truth_match_found;
    row.truth_hepmc_association_valid = label.truth_hepmc_association_valid;
    row.truth_photon_class = label.truth_photon_class;
    row.truth_is_prompt = label.truth_is_prompt;
    row.truth_iso_et = label.truth_iso_et;
    row.truth_iso_et_r03 = label.truth_iso_et_r03;
    row.truth_iso_et_r04 = label.truth_iso_et_r04;
    row.truth_iso_valid = label.truth_iso_valid;
    row.truth_iso_pass = label.truth_iso_pass;
    row.cluster_truth_track_id = label.cluster_truth_track_id;
    row.cluster_truth_pid = label.cluster_truth_pid;
    row.cluster_truth_barcode = label.cluster_truth_barcode;
    row.truth_energy_contribution = label.truth_energy_contribution;
    row.source_role = label.source_role;
    row.source_sample_code = label.source_sample_code;
    row.ppg12_source_role_label = label.ppg12_source_role_label;
    row.npb_label = label.npb_label;
    row.is_npb = label.is_npb;
    row.minimum_bias_classifier_decision = label.minimum_bias_classifier_decision;
    row.ppg12_analysis_window_pass = label.ppg12_analysis_window_pass;
    row.ppg12_response_window_pass = label.ppg12_response_window_pass;
    row.ppg12_truth_window_pass_r04 = label.ppg12_truth_window_pass_r04;
    row.ppg12_sample_bin = label.ppg12_sample_bin;
    row.ppg12_xsec_pb = label.ppg12_xsec_pb;
    row.ppg12_xsec_weight = label.ppg12_xsec_weight;
    row.ppg12_window_low = label.ppg12_window_low;
    row.ppg12_window_high = label.ppg12_window_high;
    row.max_truth_jet_pt_r04 = label.max_truth_jet_pt_r04;
    row.weight_slice = label.weight_slice;
    row.weight_cross_section = label.weight_cross_section;
    row.weight_vertex = label.weight_vertex;
    row.weight_si_di = label.weight_si_di;
    row.weight_period = label.weight_period;
    row.weight_exposure = label.weight_exposure;
    row.weight_final = label.weight_final;
    row.weight_application_count = label.weight_application_count;
  }

  void writeMeta(const char* name, const std::string& value)
  {
    TNamed metadata(name, value.c_str());
    metadata.Write(name, TObject::kOverwrite);
  }

  static SerializedUInt64 serialize(std::uint64_t value)
  {
    return static_cast<SerializedUInt64>(value);
  }

  static SerializedInt64 serialize(std::int64_t value)
  {
    return static_cast<SerializedInt64>(value);
  }

  static void serialize(SerializedIdentity128& output, const Identity128& input)
  {
    output.hi = serialize(input.hi);
    output.lo = serialize(input.lo);
  }

  void syncSerializedRow()
  {
    serialize(m_serializedId, m_row.id);
    serialize(m_serializedSourceId, m_row.source_id);
    serialize(m_serializedEventId, m_row.event_id);
    serialize(m_serializedCandidateId, m_row.candidate_id);
    serialize(m_serializedDefinitionId, m_row.definition_id);
    m_serializedEventSequence = serialize(m_row.event_sequence);
    m_serializedClusterMapKey = serialize(m_row.cluster_map_key);
  }

  static void unsigned64(TTree* tree, const char* name, SerializedUInt64* address)
  {
    const std::string leaf = std::string(name) + "/l";
    tree->Branch(name, address, leaf.c_str());
  }

  static void signed64(TTree* tree, const char* name, SerializedInt64* address)
  {
    const std::string leaf = std::string(name) + "/L";
    tree->Branch(name, address, leaf.c_str());
  }

  static void idBranches(TTree* tree, const char* prefix, SerializedIdentity128* id)
  {
    const std::string hi = std::string(prefix) + "_hi";
    const std::string lo = std::string(prefix) + "_lo";
    unsigned64(tree, hi.c_str(), &id->hi);
    unsigned64(tree, lo.c_str(), &id->lo);
  }

  template <class T>
  static void scalar(TTree* tree, const char* name, T* address, const char* type)
  {
    const std::string leaf = std::string(name) + "/" + type;
    tree->Branch(name, address, leaf.c_str());
  }

  static void text(TTree* tree, const char* name, std::string* address)
  {
    tree->Branch(name, address);
  }

  void bookTree()
  {
    m_file->cd();
    m_tree = new TTree(kTreeName, "Labeled photon-ID rows by shower-definition view");
    m_tree->SetDirectory(m_file.get());
    m_tree->SetAutoFlush(-5000000);
    idBranches(m_tree, "training_view_id", &m_serializedId);
    idBranches(m_tree, "source_occurrence_id", &m_serializedSourceId);
    idBranches(m_tree, "event_id", &m_serializedEventId);
    idBranches(m_tree, "candidate_id", &m_serializedCandidateId);
    idBranches(m_tree, "definition_id", &m_serializedDefinitionId);
    text(m_tree, "definition_name", &m_row.definition_name);
    text(m_tree, "shower_semantic_sha256", &m_row.shower_semantic_sha256);
    text(m_tree, "feature_contract_sha256", &m_row.feature_contract_sha256);
    m_tree->Branch("ordered_features", &m_row.ordered_features);
    scalar(m_tree, "feature_count", &m_row.feature_count, "I");
    scalar(m_tree, "finite_feature_state", &m_row.finite_feature_state, "I");
    scalar(m_tree, "system_code", &m_row.system_code, "I");
    text(m_tree, "source_lane", &m_row.source_lane);
    text(m_tree, "source_dataset", &m_row.source_dataset);
    text(m_tree, "source_sample", &m_row.source_sample);
    text(m_tree, "source_period", &m_row.source_period);
    text(m_tree, "source_si_di_role", &m_row.source_si_di_role);
    text(m_tree, "source_ownership_state", &m_row.source_ownership_state);
    text(m_tree, "source_manifest_sha256", &m_row.source_manifest_sha256);
    text(m_tree, "input_uri_sha256", &m_row.input_uri_sha256);
    text(m_tree, "input_file_sha256", &m_row.input_file_sha256);
    scalar(m_tree, "source_run", &m_row.source_run, "I");
    scalar(m_tree, "source_segment", &m_row.source_segment, "I");
    scalar(m_tree, "run", &m_row.run, "I");
    signed64(m_tree, "event_sequence", &m_serializedEventSequence);
    scalar(m_tree, "encounter_ordinal", &m_row.encounter_ordinal, "I");
    unsigned64(m_tree, "cluster_map_key", &m_serializedClusterMapKey);
    scalar(m_tree, "cluster_index", &m_row.cluster_index, "I");
    scalar(m_tree, "cluster_Et", &m_row.cluster_et, "D");
    scalar(m_tree, "cluster_Eta", &m_row.eta, "D");
    scalar(m_tree, "cluster_Phi", &m_row.phi, "D");
    scalar(m_tree, "vertexz", &m_row.vertex_z, "D");
    scalar(m_tree, "centrality", &m_row.centrality, "D");
    scalar(m_tree, "event_weight", &m_row.event_weight, "D");
    scalar(m_tree, "model_domain_state", &m_row.model_domain_state, "I");
    scalar(m_tree, "below15_retention_state", &m_row.below15_retention_state, "I");
    scalar(m_tree, "nominal_training_eligible", &m_row.nominal_training_eligible, "I");
    scalar(m_tree, "working_point_state", &m_row.working_point_state, "I");
    scalar(m_tree, "tag_state", &m_row.tag_state, "I");
    scalar(m_tree, "training_label", &m_row.training_label, "I");
    scalar(m_tree, "is_signal", &m_row.is_signal, "I");
    text(m_tree, "label_authority", &m_row.label_authority);
    scalar(m_tree, "truth_match_found", &m_row.truth_match_found, "I");
    scalar(m_tree, "truth_hepmc_association_valid", &m_row.truth_hepmc_association_valid, "I");
    scalar(m_tree, "truth_photon_class", &m_row.truth_photon_class, "I");
    scalar(m_tree, "truth_is_prompt", &m_row.truth_is_prompt, "I");
    scalar(m_tree, "truth_iso_et", &m_row.truth_iso_et, "D");
    scalar(m_tree, "truth_iso_et_r03", &m_row.truth_iso_et_r03, "D");
    scalar(m_tree, "truth_iso_et_r04", &m_row.truth_iso_et_r04, "D");
    scalar(m_tree, "truth_iso_valid", &m_row.truth_iso_valid, "I");
    scalar(m_tree, "truth_iso_pass", &m_row.truth_iso_pass, "I");
    scalar(m_tree, "cluster_truth_track_id", &m_row.cluster_truth_track_id, "I");
    scalar(m_tree, "cluster_truth_pid", &m_row.cluster_truth_pid, "I");
    scalar(m_tree, "cluster_truth_barcode", &m_row.cluster_truth_barcode, "I");
    scalar(m_tree, "truth_energy_contribution", &m_row.truth_energy_contribution, "D");
    scalar(m_tree, "source_role", &m_row.source_role, "I");
    scalar(m_tree, "source_sample_code", &m_row.source_sample_code, "I");
    scalar(m_tree, "ppg12_source_role_label", &m_row.ppg12_source_role_label, "I");
    scalar(m_tree, "npb_label", &m_row.npb_label, "I");
    scalar(m_tree, "is_npb", &m_row.is_npb, "I");
    scalar(m_tree, "minimum_bias_classifier_decision", &m_row.minimum_bias_classifier_decision, "I");
    scalar(m_tree, "ppg12_analysis_window_pass", &m_row.ppg12_analysis_window_pass, "I");
    scalar(m_tree, "ppg12_response_window_pass", &m_row.ppg12_response_window_pass, "I");
    scalar(m_tree, "ppg12_truth_window_pass_r04", &m_row.ppg12_truth_window_pass_r04, "I");
    scalar(m_tree, "ppg12_sample_bin", &m_row.ppg12_sample_bin, "I");
    scalar(m_tree, "ppg12_xsec_pb", &m_row.ppg12_xsec_pb, "D");
    scalar(m_tree, "ppg12_xsec_weight", &m_row.ppg12_xsec_weight, "D");
    scalar(m_tree, "ppg12_window_low", &m_row.ppg12_window_low, "D");
    scalar(m_tree, "ppg12_window_high", &m_row.ppg12_window_high, "D");
    scalar(m_tree, "max_truth_jet_pt_r04", &m_row.max_truth_jet_pt_r04, "D");
    scalar(m_tree, "weight_slice", &m_row.weight_slice, "D");
    scalar(m_tree, "weight_cross_section", &m_row.weight_cross_section, "D");
    scalar(m_tree, "weight_vertex", &m_row.weight_vertex, "D");
    scalar(m_tree, "weight_si_di", &m_row.weight_si_di, "D");
    scalar(m_tree, "weight_period", &m_row.weight_period, "D");
    scalar(m_tree, "weight_exposure", &m_row.weight_exposure, "D");
    scalar(m_tree, "weight_final", &m_row.weight_final, "D");
    scalar(m_tree, "weight_application_count", &m_row.weight_application_count, "I");
  }

  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  Row m_row;
  SerializedIdentity128 m_serializedId;
  SerializedIdentity128 m_serializedSourceId;
  SerializedIdentity128 m_serializedEventId;
  SerializedIdentity128 m_serializedCandidateId;
  SerializedIdentity128 m_serializedDefinitionId;
  SerializedInt64 m_serializedEventSequence = 0;
  SerializedUInt64 m_serializedClusterMapKey = 0;
  int m_systemCode = 0;
  SourceOccurrenceRow m_source;
  std::string m_configSha256, m_codeSha256;
  std::map<Key, Label> m_labels;
  std::set<Key> m_consumed;
  std::set<Key> m_seenCandidates;
  std::set<std::uint64_t> m_seenClusterKeys;
  std::uint64_t m_nextCompactCandidate=0;
  std::int64_t m_activeEvent=-1,m_lastFinishedEvent=-1;
  bool m_compactInputs=false;
  std::vector<Row> m_pendingRows;
  std::int64_t m_entries = 0;
  bool m_initialized = false, m_finished = false, m_failed = false;
};
}  // namespace RJPhotonTrainingViewV1

#endif
