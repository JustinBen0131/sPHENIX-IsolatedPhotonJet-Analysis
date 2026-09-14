#ifndef RECOILJETS_THE106OBSERVATIONDISABLED_H
#define RECOILJETS_THE106OBSERVATIONDISABLED_H

#include <cstddef>
#include <cstdint>
#include <string>

// Compile-time hard-disable adapter for the optional THE-106 observer API.
// Every function is intentionally inline and inert.  It introduces no runtime
// component and preserves the THE-117 direct-path independence certificate.
namespace the106
{
namespace c0h2
{
struct PairToken
{
  std::uint64_t generation = 0;
  std::uint64_t pair_ordinal = 0;
};

struct CandidateContext
{
  PairToken pair = {};
  std::uint64_t delivered_event_ordinal = 0;
  bool run_valid = false;
  std::int64_t run_number = 0;
  bool event_valid = false;
  std::int64_t event_number = 0;
  std::uint64_t container_key = 0;
  std::uint64_t cluster_id = 0;
  std::uint64_t producer_encounter_ordinal = 0;
  const char* module_name = nullptr;
  bool canonical_view = false;
  const char* view_key = nullptr;
};

enum class ScoreStatus : std::uint8_t
{
  valid,
  feature_nonfinite,
  output_empty,
  output_nonfinite
};

class ScopedCandidateContext
{
 public:
  explicit ScopedCandidateContext(const CandidateContext&) noexcept {}
};

inline bool scoreObservationEnabled() noexcept { return false; }
inline PairToken currentFrameworkPairToken() noexcept { return {}; }
inline void emitScoreObservation(ScoreStatus, const float*, std::size_t,
                                 const std::string*, const char*, const char*,
                                 bool, float, std::size_t) noexcept {}
}  // namespace c0h2

namespace c0rh
{
enum class EventProcessOutcome : std::uint8_t
{
  entered,
  mandatory_nodes_missing,
  first_event_gate_rejected,
  embedded_inclusive_stitch_rejected,
  embedded_photon_stitch_rejected,
  centrality_invalid,
  vertex_rejected,
  embedded_minbias_rejected,
  audit_abort_run,
  candidate_processing_exception,
  completed,
  unclassified_exit
};

enum class RawQATightTag : std::uint8_t
{
  not_evaluated,
  preselection_fail,
  tight,
  non_tight,
  neither
};

struct RawQACandidateObservation
{
  c0h2::CandidateContext context = {};
  const double* values = nullptr;
  const std::uint64_t* value_bits = nullptr;
  std::size_t value_count = 0;
  std::uint16_t valid_mask = 0;
  const std::string* variable_names = nullptr;
  const std::string* active_triggers = nullptr;
  std::size_t active_trigger_count = 0;
  int photon_pt_slice = 0;
  int centrality_slice = 0;
  bool canonical_view = false;
  const char* view_suffix = nullptr;
  RawQATightTag tight_tag = RawQATightTag::not_evaluated;
  bool preselection_pass = false;
  bool direct_candidate_admitted = false;
};

struct RawQAFillObservation
{
  c0h2::CandidateContext context = {};
  const char* variable_name = nullptr;
  const char* trigger_name = nullptr;
  const char* tag_name = nullptr;
  const char* view_suffix = nullptr;
  const char* directory_name = nullptr;
  const char* object_name = nullptr;
  const char* object_class = nullptr;
  const char* object_contract_key = nullptr;
  int photon_pt_slice = 0;
  int centrality_slice = 0;
  bool canonical_view = false;
  RawQATightTag tight_tag = RawQATightTag::not_evaluated;
  double value = 0.0;
  std::uint64_t value_bits = 0;
  double weight = 0.0;
  std::uint64_t weight_bits = 0;
  std::uint64_t fill_ordinal = 0;
  bool value_valid = false;
  bool filled = false;
  bool sumw2_enabled = false;
};

class ScopedEventObservation
{
 public:
  ScopedEventObservation(std::uint64_t, const char*) noexcept {}
  bool active() const noexcept { return false; }
  void bindIdentity(bool, std::int64_t, bool, std::int64_t) noexcept {}
  void setCentrality(bool, double) noexcept {}
  void setActiveTriggers(const std::string*, std::size_t) noexcept {}
  void setOutcome(EventProcessOutcome, std::int64_t, const char*) noexcept {}
};

inline bool observationEnabled() noexcept { return false; }
inline bool candidateObservationEnabled() noexcept { return false; }
inline bool fillObservationEnabled() noexcept { return false; }
inline bool emitRawQACandidateObservation(
    const RawQACandidateObservation&) noexcept { return false; }
inline void emitRawQAFillObservation(const RawQAFillObservation&) noexcept {}
inline void noteEncounteredCandidate() noexcept {}
inline void noteObservationFailure() noexcept {}
inline void noteCandidateProcessingException(const char*) noexcept {}
inline std::uint64_t nextFillOrdinal() noexcept { return 0; }
}  // namespace c0rh
}  // namespace the106

#endif
