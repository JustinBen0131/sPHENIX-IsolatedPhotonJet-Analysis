#ifndef RJ_REPLAY_RUNTIME_V1_H
#define RJ_REPLAY_RUNTIME_V1_H

#include "RJReplayFoundationV1.h"
#include "RJSchema10DataRetentionContractV1.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace RJReplayRuntimeV1
{
using namespace RJReplayFoundationV1;

inline std::string env(const char* key)
{
  const char* value = std::getenv(key);
  return value ? std::string(value) : std::string{};
}

inline bool envEnabled(const char* key)
{
  const std::string value = env(key);
  return value == "1" || value == "true" || value == "TRUE" || value == "yes" || value == "on";
}

inline int envInt(const char* key, int fallback)
{
  const std::string value = env(key);
  if (value.empty()) return fallback;
  try { return std::stoi(value); }
  catch (...) { return fallback; }
}

struct LegacyOrderCoordinates
{
  bool configured=false;
  std::int64_t source_file_ordinal=-1;
  std::int64_t source_entry_begin=-1;
  std::int64_t source_global_entry_begin=-1;
};

inline bool parseNonnegativeInt64(const std::string& raw,
                                  std::int64_t& value)
{
  if (raw.empty()) return false;
  std::int64_t parsed=0;
  for (const char character : raw)
  {
    if (character<'0'||character>'9') return false;
    const std::int64_t digit=static_cast<std::int64_t>(character-'0');
    if (parsed>(std::numeric_limits<std::int64_t>::max()-digit)/10)
      return false;
    parsed=parsed*10+digit;
  }
  value=parsed;
  return true;
}

inline bool readLegacyOrderCoordinates(LegacyOrderCoordinates& coordinates,
                                       std::string* error=nullptr)
{
  coordinates=LegacyOrderCoordinates{};
  const std::string file=env("RJ_REPLAY_SOURCE_FILE_ORDINAL");
  const std::string entry=env("RJ_REPLAY_SOURCE_ENTRY_BEGIN");
  const std::string global=env("RJ_REPLAY_SOURCE_GLOBAL_ENTRY_BEGIN");
  const unsigned int present=static_cast<unsigned int>(!file.empty())+
      static_cast<unsigned int>(!entry.empty())+
      static_cast<unsigned int>(!global.empty());
  if (present==0) return true;
  if (present!=3)
  {
    if (error) *error="partial PPG12 source-order environment contract";
    return false;
  }
  if (!parseNonnegativeInt64(file,coordinates.source_file_ordinal)||
      !parseNonnegativeInt64(entry,coordinates.source_entry_begin)||
      !parseNonnegativeInt64(global,coordinates.source_global_entry_begin))
  {
    if (error) *error="invalid PPG12 source-order environment coordinate";
    coordinates=LegacyOrderCoordinates{};
    return false;
  }
  coordinates.configured=true;
  return true;
}

inline bool applyLegacyOrderCoordinates(const LegacyOrderCoordinates& coordinates,
                                        std::int64_t local_entry_offset,
                                        EventRow& event,
                                        std::string* error=nullptr)
{
  if (!coordinates.configured) return true;
  if (local_entry_offset<0||
      coordinates.source_entry_begin>
          std::numeric_limits<std::int64_t>::max()-local_entry_offset||
      coordinates.source_global_entry_begin>
          std::numeric_limits<std::int64_t>::max()-local_entry_offset)
  {
    if (error) *error="PPG12 source-order event coordinate overflow";
    return false;
  }
  event.source_file_ordinal=coordinates.source_file_ordinal;
  event.source_entry_ordinal=coordinates.source_entry_begin+local_entry_offset;
  event.source_global_entry_ordinal=
      coordinates.source_global_entry_begin+local_entry_offset;
  return true;
}

struct EventBundle
{
  Schema10DataRetentionContractV1::EventClass retention_class =
      Schema10DataRetentionContractV1::EventClass::FULL_UNFILTERED;
  bool is_data = false;
  EventRow event;
  std::vector<PhotonCandidateRow> candidates;
  std::vector<ModelEvaluationRow> models;
  std::vector<ShowerCellRow> shower_cells;
  std::vector<ShowerFeatureViewRow> shower_feature_views;
  std::vector<IsolationConstituentRow> isolation_constituents;
  std::vector<IsolationWitnessRow> isolation_witnesses;
  std::vector<JetRow> jets;
  std::vector<JetConstituentRow> jet_constituents;
  std::vector<PhotonJetPairRow> pairs;
  std::vector<TruthPhotonRow> truth_photons;
  std::vector<TruthPhotonMissOccurrenceRow> truth_photon_miss_occurrences;
  std::vector<EmbeddedPhotonDiagnosticOccurrenceRow> embedded_photon_diagnostic_occurrences;
  std::vector<PPG12DiagnosticOccurrenceRow> ppg12_diagnostic_occurrences;
  std::vector<TruthJetRow> truth_jets;
  std::vector<RecoTruthLinkRow> links;
  std::vector<WeightComponentRow> weights;
  std::vector<EventDisplaySnapshotRow> snapshots;
};

class Runtime
{
 public:
  bool initialize(TFile* file, const SourceOccurrenceRow& source, std::string* error = nullptr)
  {
    return initialize(file,source,WriterMode::SERIALIZE,error);
  }

  bool initialize(TFile* file,
                  const SourceOccurrenceRow& source,
                  WriterMode mode,
                  std::string* error = nullptr)
  {
    Metadata metadata;
    metadata.schema_sha256 = env("RJ_REPLAY_SCHEMA_SHA256");
    metadata.semantic_sha256 = env("RJ_REPLAY_SEMANTIC_SHA256");
    metadata.source_sha256 = env("RJ_REPLAY_SOURCE_SHA256");
    metadata.model_sha256 = env("RJ_REPLAY_MODEL_SHA256");
    metadata.config_sha256 = env("RJ_REPLAY_CONFIG_SHA256");
    metadata.code_sha256 = env("RJ_REPLAY_CODE_SHA256");
    metadata.provenance_manifest_sha256 = env("RJ_REPLAY_PROVENANCE_MANIFEST_SHA256");
    metadata.analysis_jet_node_suffix = env("RJ_ANALYSIS_JET_NODE_SUFFIX");
    m_data_retention_profile = env("RJ_SCHEMA10_DATA_RETENTION_PROFILE");
    // Broad trigger-independent DATA capture is the default for new trees.
    // Explicit historical sparse profiles remain reproducible by name/digest.
    if (m_data_retention_profile.empty())
      m_data_retention_profile = Schema10DataRetentionContractV1::kWorkfestProfile;
    if (!Schema10DataRetentionContractV1::validProfile(m_data_retention_profile))
    {
      if (error) *error = "unknown RJ_SCHEMA10_DATA_RETENTION_PROFILE";
      return false;
    }
    if (!m_data_retention_profile.empty())
    {
      metadata.data_retention_profile = m_data_retention_profile;
      metadata.data_retention_contract_sha256 =
          Schema10DataRetentionContractV1::contractSha256(m_data_retention_profile);
      if (metadata.data_retention_contract_sha256.empty())
      {
        if (error) *error = "no contract digest for RJ_SCHEMA10_DATA_RETENTION_PROFILE";
        return false;
      }
    }
    if (!readLegacyOrderCoordinates(m_legacy_order_coordinates,error)) return false;
    if (!m_writer.initialize(file, metadata, mode, error)) return false;
    if (!m_writer.fill(source, error)) return false;
    m_source_id = source.id;
    return true;
  }

  Schema10DataRetentionContractV1::Decision classifyAndRecord(
      bool isSimulation,
      bool isAuAu,
      double vertexZ,
      double centrality,
      const std::vector<PhotonCandidateRow>& candidates)
  {
    std::vector<Schema10DataRetentionContractV1::CandidateKinematics> view;
    view.reserve(candidates.size());
    for (const auto& candidate : candidates)
      view.push_back({candidate.cluster_et,candidate.eta,
                      candidate.active_preselection_state != 0});
    const auto decision = Schema10DataRetentionContractV1::classify(
        m_data_retention_profile,isSimulation,isAuAu,vertexZ,centrality,view);
    ++m_processed_events;
    m_primary_candidates += decision.primary_candidates;
    m_extension_candidates += decision.extension_candidates;
    using EventClass = Schema10DataRetentionContractV1::EventClass;
    using OmitReason = Schema10DataRetentionContractV1::OmitReason;
    using OmittedBand = Schema10DataRetentionContractV1::OmittedCentralityBand;
    if (decision.event_class == EventClass::OMIT_DATA)
    {
      ++m_omitted_events;
      if (decision.omit_reason == OmitReason::VERTEX) ++m_omitted_vertex;
      else if (decision.omit_reason == OmitReason::CENTRALITY)
      {
        ++m_omitted_centrality;
        switch (decision.omitted_centrality_band)
        {
          case OmittedBand::NONFINITE: ++m_omitted_centrality_nonfinite; break;
          case OmittedBand::NEGATIVE: ++m_omitted_centrality_negative; break;
          case OmittedBand::NEAR_BOUNDARY: ++m_omitted_centrality_80_90; break;
          default: ++m_omitted_centrality_ge90; break;
        }
      }
      else ++m_omitted_no_candidate;
    }
    else
    {
      ++m_retained_events;
      if (decision.event_class == EventClass::PRIMARY_FULL)
        ++m_primary_events;
      else if (decision.event_class == EventClass::EXTENSION_ONLY)
        ++m_extension_events;
      else if (decision.event_class == EventClass::EVENT_ONLY)
        ++m_event_only_events;
      else
        ++m_full_unfiltered_events;
    }
    return decision;
  }

  bool write(EventBundle& bundle, std::string* error = nullptr)
  {
    using EventClass = Schema10DataRetentionContractV1::EventClass;
    if (bundle.retention_class == EventClass::OMIT_DATA)
    {
      if (error) *error = "an omitted data event cannot be serialized";
      return false;
    }
    if (bundle.is_data)
    {
      if (!bundle.truth_photons.empty() || !bundle.truth_jets.empty() ||
          !bundle.links.empty() ||
          !bundle.truth_photon_miss_occurrences.empty() ||
          !bundle.embedded_photon_diagnostic_occurrences.empty())
      {
        if (error) *error = "data truth/link tables must remain empty";
        return false;
      }
      if (bundle.retention_class == EventClass::EXTENSION_ONLY &&
          (!bundle.jets.empty() || !bundle.jet_constituents.empty() ||
           !bundle.pairs.empty()))
      {
        if (error) *error = "extension-only data event contains jet payload";
        return false;
      }
      if (bundle.retention_class == EventClass::EVENT_ONLY &&
          (!bundle.candidates.empty() || !bundle.models.empty() ||
           !bundle.shower_cells.empty() || !bundle.shower_feature_views.empty() ||
           !bundle.isolation_constituents.empty() || !bundle.isolation_witnesses.empty() ||
           !bundle.jets.empty() || !bundle.jet_constituents.empty() || !bundle.pairs.empty()))
      {
        if (error) *error = "event-only data event contains object payload";
        return false;
      }
    }
    const bool trace = envEnabled("RJ_REPLAY_TRACE");
    auto mark = [&](const char* stage)
    {
      if (trace)
      {
        std::cerr << "RJ_REPLAY_TRACE event=" << bundle.event.event_sequence
                  << " stage=" << stage
                  << " candidates=" << bundle.candidates.size()
                  << " models=" << bundle.models.size()
                  << " shower_cells=" << bundle.shower_cells.size()
                  << " shower_feature_views=" << bundle.shower_feature_views.size()
                  << " iso_constituents=" << bundle.isolation_constituents.size()
                  << " iso_witnesses=" << bundle.isolation_witnesses.size()
                  << " jets=" << bundle.jets.size()
                  << " jet_constituents=" << bundle.jet_constituents.size()
                  << " pairs=" << bundle.pairs.size()
                  << " truth_photons=" << bundle.truth_photons.size()
                  << " truth_photon_miss_occurrences=" << bundle.truth_photon_miss_occurrences.size()
                  << " embedded_photon_diagnostic_occurrences=" << bundle.embedded_photon_diagnostic_occurrences.size()
                  << " ppg12_diagnostic_occurrences=" << bundle.ppg12_diagnostic_occurrences.size()
                  << " truth_jets=" << bundle.truth_jets.size()
                  << " links=" << bundle.links.size()
                  << std::endl;
      }
    };

    if (bundle.event.source_id.isNull()) bundle.event.source_id = m_source_id;
    mark("before_event_fill");
    if (!m_writer.fill(bundle.event, error)) return false;
    mark("after_event_fill");
    for (const auto& row : bundle.candidates) if (!m_writer.fill(row, error)) return false;
    mark("after_candidate_fill");
    for (const auto& row : bundle.models) if (!m_writer.fill(row, error)) return false;
    mark("after_model_fill");
    for (const auto& row : bundle.shower_cells) if (!m_writer.fill(row, error)) return false;
    mark("after_shower_fill");
    for (const auto& row : bundle.shower_feature_views) if (!m_writer.fill(row, error)) return false;
    mark("after_shower_feature_view_fill");
    for (const auto& row : bundle.isolation_constituents) if (!m_writer.fill(row, error)) return false;
    mark("after_iso_constituent_fill");
    for (const auto& row : bundle.isolation_witnesses) if (!m_writer.fill(row, error)) return false;
    mark("after_iso_witness_fill");
    for (const auto& row : bundle.jets) if (!m_writer.fill(row, error)) return false;
    mark("after_jet_fill");
    for (const auto& row : bundle.jet_constituents) if (!m_writer.fill(row, error)) return false;
    mark("after_jet_constituent_fill");
    for (const auto& row : bundle.pairs) if (!m_writer.fill(row, error)) return false;
    mark("after_pair_fill");
    for (const auto& row : bundle.truth_photons) if (!m_writer.fill(row, error)) return false;
    mark("after_truth_photon_fill");
    for (const auto& row : bundle.truth_photon_miss_occurrences) if (!m_writer.fill(row, error)) return false;
    mark("after_truth_photon_miss_occurrence_fill");
    for (const auto& row : bundle.embedded_photon_diagnostic_occurrences) if (!m_writer.fill(row, error)) return false;
    mark("after_embedded_photon_diagnostic_occurrence_fill");
    for (const auto& row : bundle.ppg12_diagnostic_occurrences) if (!m_writer.fill(row, error)) return false;
    mark("after_ppg12_diagnostic_occurrence_fill");
    for (const auto& row : bundle.truth_jets) if (!m_writer.fill(row, error)) return false;
    mark("after_truth_jet_fill");
    for (const auto& row : bundle.links) if (!m_writer.fill(row, error)) return false;
    mark("after_link_fill");
    for (const auto& row : bundle.weights) if (!m_writer.fill(row, error)) return false;
    mark("after_weight_fill");
    for (const auto& row : bundle.snapshots) if (!m_writer.fill(row, error)) return false;
    mark("complete");
    ++m_table_counts[0];
    m_table_counts[1] += bundle.candidates.size();
    m_table_counts[2] += bundle.models.size();
    m_table_counts[3] += bundle.shower_cells.size();
    m_table_counts[4] += bundle.shower_feature_views.size();
    m_table_counts[5] += bundle.isolation_constituents.size();
    m_table_counts[6] += bundle.isolation_witnesses.size();
    m_table_counts[7] += bundle.jets.size();
    m_table_counts[8] += bundle.jet_constituents.size();
    m_table_counts[9] += bundle.pairs.size();
    m_table_counts[10] += bundle.truth_photons.size();
    m_table_counts[11] += bundle.truth_photon_miss_occurrences.size();
    m_table_counts[12] += bundle.embedded_photon_diagnostic_occurrences.size();
    m_table_counts[13] += bundle.ppg12_diagnostic_occurrences.size();
    m_table_counts[14] += bundle.truth_jets.size();
    m_table_counts[15] += bundle.links.size();
    m_table_counts[16] += bundle.weights.size();
    m_table_counts[17] += bundle.snapshots.size();
    return true;
  }

  bool finish(std::string* error = nullptr)
  {
    static const std::array<const char*,18> tableNames{{
      "RJEventV1","RJPhotonCandidateV1","RJModelEvaluationV1",
      "RJShowerCellV1","RJShowerFeatureViewV1","RJIsolationConstituentV1",
      "RJIsolationWitnessV1","RJJetV1","RJJetConstituentV1",
      "RJPhotonJetPairV1","RJTruthPhotonV1",
      "RJTruthPhotonMissOccurrenceV1","RJEmbeddedPhotonDiagnosticOccurrenceV1",
      "RJPPG12DiagnosticOccurrenceV1","RJTruthJetV1","RJRecoTruthLinkV1",
      "RJWeightComponentV1","RJEventDisplaySnapshotV1"}};
    CompletionMetadata completion;
    completion.counters = {
      {"processed_events",m_processed_events},
      {"retained_events",m_retained_events},
      {"omitted_events",m_omitted_events},
      {"omitted_vertex_events",m_omitted_vertex},
      {"omitted_centrality_events",m_omitted_centrality},
      {"omitted_centrality_nonfinite_events",m_omitted_centrality_nonfinite},
      {"omitted_centrality_negative_events",m_omitted_centrality_negative},
      {"omitted_centrality_80_to_90_events",m_omitted_centrality_80_90},
      {"omitted_centrality_ge90_events",m_omitted_centrality_ge90},
      {"omitted_no_candidate_events",m_omitted_no_candidate},
      {"retained_primary_events",m_primary_events},
      {"retained_extension_events",m_extension_events},
      {"retained_event_only_events",m_event_only_events},
      {"full_unfiltered_events",m_full_unfiltered_events},
      {"retained_primary_candidates",m_primary_candidates},
      {"retained_extension_candidates",m_extension_candidates},
      {"RJSourceOccurrenceV1_entries",1U}};
    for (std::size_t index=0; index<tableNames.size(); ++index)
      completion.counters.emplace_back(
          std::string(tableNames[index])+"_entries",m_table_counts[index]);
    return m_writer.finish(completion,error);
  }

  const Metadata& metadata() const { return m_writer.metadata(); }

  bool captureObjectDomain(bool isSimulation, double vertexZ) const
  {
    return Schema10DataRetentionContractV1::captureObjectDomain(
        m_data_retention_profile, isSimulation, vertexZ);
  }

  bool applyLegacyOrderCoordinates(EventRow& event,
                                   std::int64_t local_entry_offset,
                                   std::string* error=nullptr) const
  {
    return RJReplayRuntimeV1::applyLegacyOrderCoordinates(
        m_legacy_order_coordinates,local_entry_offset,event,error);
  }

 private:
  Writer m_writer;
  Identity128 m_source_id;
  LegacyOrderCoordinates m_legacy_order_coordinates;
  std::string m_data_retention_profile;
  std::uint64_t m_processed_events=0,m_retained_events=0,m_omitted_events=0;
  std::uint64_t m_omitted_vertex=0,m_omitted_centrality=0,m_omitted_no_candidate=0;
  std::uint64_t m_omitted_centrality_nonfinite=0,m_omitted_centrality_negative=0;
  std::uint64_t m_omitted_centrality_80_90=0,m_omitted_centrality_ge90=0;
  std::uint64_t m_primary_events=0,m_extension_events=0,m_full_unfiltered_events=0;
  std::uint64_t m_event_only_events=0;
  std::uint64_t m_primary_candidates=0,m_extension_candidates=0;
  std::array<std::uint64_t,18> m_table_counts{};
};

template <class Function>
class ScopeExit
{
 public:
  explicit ScopeExit(Function&& function) : m_function(std::move(function)) {}
  ScopeExit(const ScopeExit&) = delete;
  ScopeExit& operator=(const ScopeExit&) = delete;
  ~ScopeExit() { m_function(); }

 private:
  Function m_function;
};

template <class Function>
ScopeExit<Function> onScopeExit(Function&& function)
{
  return ScopeExit<Function>(std::forward<Function>(function));
}
}  // namespace RJReplayRuntimeV1

#endif
