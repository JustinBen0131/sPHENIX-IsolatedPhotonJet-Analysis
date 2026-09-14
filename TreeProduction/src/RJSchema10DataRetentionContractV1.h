#ifndef RJ_SCHEMA10_DATA_RETENTION_CONTRACT_V1_H
#define RJ_SCHEMA10_DATA_RETENTION_CONTRACT_V1_H

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace Schema10DataRetentionContractV1
{
// Two named data retention profiles exist and they retain different event
// populations.  `sparse_photon_analysis_v1` is the profile every already
// Produced schema-10 data row was written under: an event is retained
// when any stored photon candidate is inside the kinematic window.
// `sparse_photon_analysis_v2` additionally requires that candidate to pass the
// active loose-photon preselection.  v2 therefore drops every event whose
// candidates all fail preselection, which removes the data-side denominator
// for a preselection efficiency measured on those events.  The profile name is
// part of the product identity; never silently move rows between the two.
constexpr const char* kProfileName = "sparse_photon_analysis_v1";
constexpr const char* kProfileNameV2 = "sparse_photon_analysis_v2";
constexpr const char* kWorkfestProfile = "workfest_replay_capture_v1";

// These texts are the single content-addressed authority for the schema-10
// data retention policy.  They are deliberately independent of the physical
// ROOT schema: changing one changes the configuration/provenance identity, not
// the schema version or any simulation table.
inline const std::string& canonicalContractTextV1()
{
  static const std::string value =
      "Schema10DataRetentionContractV1|"
      "data_only=1|photon_et=[5,40)|photon_abs_eta<0.7|abs_vertex_z<60|"
      "auau_centrality=[0,80)|primary_et=[15,35)|"
      "extension_et=[5,15)+[35,40)|empty=histograms_only|"
      "primary=event+photon+model+shower+isolation+ueflow+jets+areas+pairs+"
      "jet_constituents_corrected_pt_ge_5|"
      "extension=event+photon+model+shower+isolation+ueflow|"
      "data_truth_links=empty|simulation=unchanged";
  return value;
}

inline const std::string& canonicalContractTextV2()
{
  static const std::string value =
      "Schema10DataRetentionContractV1|"
      "data_only=1|photon_et=[5,40)|photon_abs_eta<0.7|abs_vertex_z<60|"
      "loose_photon_active_preselection=required|"
      "auau_centrality=[0,80)|primary_et=[15,35)|"
      "extension_et=[5,15)+[35,40)|empty=histograms_only|"
      "primary=event+photon+model+shower+isolation+ueflow+jets+areas+pairs+"
      "jet_constituents_corrected_pt_ge_5|"
      "extension=event+photon+model+shower+isolation+ueflow|"
      "data_truth_links=empty|simulation=unchanged";
  return value;
}

// Retained for source compatibility with the single-profile header.
inline const std::string& canonicalContractText() { return canonicalContractTextV1(); }

// Separate workfest profile: old sparse products keep their exact semantics.
// Event-only means an event that reached this producer, not a claim that every
// original DST event survived the upstream reconstruction/trigger chain.
inline const std::string& canonicalWorkfestContractText()
{
  static const std::string value =
      "Schema10DataRetentionContractV1|data_only=1|"
      "event=every_producer_encounter+identity+raw_mbd+trigger+quality+weights|"
      "photon_et=[5,40)|photon_abs_eta<0.7|object_abs_vertex_z<60|"
      "auau_centrality=no_retention_cut|loose_photon_active_preselection=not_required|"
      "bdt=not_a_retention_gate|empty=event_only|"
      "photon_event=event+photon+model+shower+isolation+ueflow+jets+areas+pairs+"
      "jet_constituents_corrected_pt_ge_5|data_truth_links=empty|"
      "simulation=unchanged|upstream_population=separately_accounted";
  return value;
}

inline bool captureObjectDomain(const std::string& profile,
                                bool isSimulation, double vertexZ)
{
  return profile != kWorkfestProfile || isSimulation ||
         (std::isfinite(vertexZ) && std::fabs(vertexZ) < 60.0);
}

inline bool validProfile(const std::string& profile)
{
  return profile.empty() || profile == kProfileName || profile == kProfileNameV2 ||
         profile == kWorkfestProfile;
}

inline bool requiresLoosePhoton(const std::string& profile)
{
  return profile == kProfileNameV2;
}

// SHA-256 of the matching canonical text.  Standalone constants keep the
// policy unit-testable without ROOT while the Python contract test recomputes
// both digests from the same strings.
inline std::string contractSha256(const std::string& profile)
{
  if (profile == kWorkfestProfile)
    return "025304c37c51bcaae16b87ec7094f2a2ff03c7917f069bb51ae58b3fa55073b9";
  if (profile == kProfileNameV2)
    return "227167f8e98fdee6e44de0b57c843ff0842638f0d7c747fa68256a1a4b546d74";
  if (profile == kProfileName)
    return "1388f6e3c89703d3f7a6d9b7f66371403b2aa4d861376af18cc0d21bd3b54df4";
  return std::string{};
}

enum class EventClass
{
  OMIT_DATA = 0,
  EXTENSION_ONLY = 1,
  PRIMARY_FULL = 2,
  FULL_UNFILTERED = 3,
  EVENT_ONLY = 4
};

// First failing gate, in evaluation order.  Diagnostics only: the omitted
// population is a counted, reason-resolved quantity instead of one scalar.
enum class OmitReason
{
  NONE = 0,
  VERTEX = 1,
  CENTRALITY = 2,
  NO_QUALIFYING_CANDIDATE = 3
};

// Coarse Au+Au centrality band of an event the centrality gate omitted.  A
// later centrality recalibration can only move events across the 80% boundary
// from a band that was actually written, so the size of NEAR_BOUNDARY is the
// honest bound on what a tree-only centrality update cannot recover.
enum class OmittedCentralityBand
{
  NOT_APPLICABLE = 0,
  NONFINITE = 1,
  NEGATIVE = 2,
  NEAR_BOUNDARY = 3,  // [80,90)
  FAR_PERIPHERAL = 4  // >=90
};

struct CandidateKinematics
{
  double et = 0.0;
  double eta = 0.0;
  bool loose_photon = false;
};

struct Decision
{
  EventClass event_class = EventClass::FULL_UNFILTERED;
  std::size_t primary_candidates = 0;
  std::size_t extension_candidates = 0;
  OmitReason omit_reason = OmitReason::NONE;
  OmittedCentralityBand omitted_centrality_band =
      OmittedCentralityBand::NOT_APPLICABLE;

  bool retainEvent() const { return event_class != EventClass::OMIT_DATA; }
  bool retainJetPayload() const
  {
    return event_class == EventClass::PRIMARY_FULL ||
           event_class == EventClass::FULL_UNFILTERED;
  }
};

inline OmittedCentralityBand centralityBand(double centrality)
{
  if (!std::isfinite(centrality)) return OmittedCentralityBand::NONFINITE;
  if (centrality < 0.0) return OmittedCentralityBand::NEGATIVE;
  if (centrality < 90.0) return OmittedCentralityBand::NEAR_BOUNDARY;
  return OmittedCentralityBand::FAR_PERIPHERAL;
}

inline Decision classify(const std::string& profile,
                         bool isSimulation,
                         bool isAuAu,
                         double vertexZ,
                         double centrality,
                         const std::vector<CandidateKinematics>& candidates)
{
  // An absent profile preserves the historical full writer.  Simulation is
  // full even when a sparse data profile is present in a shared runtime.
  if (profile.empty() || isSimulation)
    return {EventClass::FULL_UNFILTERED, 0U, 0U, OmitReason::NONE,
            OmittedCentralityBand::NOT_APPLICABLE};

  if (profile == kWorkfestProfile)
  {
    Decision result;
    result.event_class = EventClass::EVENT_ONLY;
    if (!captureObjectDomain(profile, isSimulation, vertexZ)) return result;
    for (const auto& candidate : candidates)
    {
      if (!std::isfinite(candidate.et) || !std::isfinite(candidate.eta) ||
          candidate.et < 5.0 || candidate.et >= 40.0 ||
          std::fabs(candidate.eta) >= 0.7) continue;
      // All captured photon energies get the same full recoil support.
      ++result.primary_candidates;
    }
    if (result.primary_candidates) result.event_class = EventClass::PRIMARY_FULL;
    return result;
  }

  const bool requireLoosePhoton = requiresLoosePhoton(profile);
  Decision result;
  result.event_class = EventClass::OMIT_DATA;
  if (!std::isfinite(vertexZ) || std::fabs(vertexZ) >= 60.0)
  {
    result.omit_reason = OmitReason::VERTEX;
    return result;
  }
  if (isAuAu &&
      (!std::isfinite(centrality) || centrality < 0.0 || centrality >= 80.0))
  {
    result.omit_reason = OmitReason::CENTRALITY;
    result.omitted_centrality_band = centralityBand(centrality);
    return result;
  }

  for (const auto& candidate : candidates)
  {
    if ((requireLoosePhoton && !candidate.loose_photon) ||
        !std::isfinite(candidate.et) || !std::isfinite(candidate.eta) ||
        candidate.et < 5.0 || candidate.et >= 40.0 ||
        std::fabs(candidate.eta) >= 0.7)
      continue;
    if (candidate.et >= 15.0 && candidate.et < 35.0)
      ++result.primary_candidates;
    else
      ++result.extension_candidates;
  }

  if (result.primary_candidates != 0U)
    result.event_class = EventClass::PRIMARY_FULL;
  else if (result.extension_candidates != 0U)
    result.event_class = EventClass::EXTENSION_ONLY;
  else
    result.omit_reason = OmitReason::NO_QUALIFYING_CANDIDATE;
  return result;
}
}  // namespace Schema10DataRetentionContractV1

#endif
