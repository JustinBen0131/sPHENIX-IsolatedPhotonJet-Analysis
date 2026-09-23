/**
 * @file Output.cc
 * @brief ROOT serialization, provenance, and completion of PhotonJetTree output.
 *
 * This file is the persistence boundary of TreeProduction. It opens the ROOT file, books the canonical tables, copies completed records into ROOT branch
 * buffers, writes provenance, and certifies successful completion.
 *
 * Output.cc does not define physics. Every scientific value written here has already been measured or derived by the capture code that owns it. This layer
 * only preserves those records and their relationships on disk.
 *
 * Data model
 * ----------
 * The output is relational:
 *
 *   - one TTree per record type;
 *   - one row per record;
 *   - flat scalar/vector leaves;
 *   - relations carried by stable 128-bit identities stored as <name>_hi and
 *     <name>_lo uint64 leaves;
 *   - persisted enums written using their fixed integer representation;
 *   - physical observables stored as double;
 *   - missing measurements preserved through their declared NaN/state
 *     semantics rather than replaced with plausible defaults.
 *
 * Canonical tables
 * ----------------
 *   Sources
 *       One row describing the physical input source and final exposure counts.
 *
 *   UpstreamRejectedEvents
 *       Input events rejected by reconstruction before PhotonJetTree ran.
 *
 *   Events
 *       One row for every producer encounter retained by TreeProduction.
 *
 *   Photons
 *       Reconstructed photon candidates and candidate-level witnesses.
 *
 *   PhotonShowerViews
 *       One row per photon and stored shower-shape definition.
 *
 *   PhotonCells
 *       Calorimeter-cell witnesses retained around photon candidates.
 *
 *   Isolation
 *       One row per photon, cone radius, and reconstructed isolation method.
 *
 *   IsolationConstituents
 *       Constituents retained for reconstruction-level isolation auditing.
 *
 *   Jets
 *       Reconstructed jets across the configured view/radius collections.
 *
 *   PhotonJetPairs
 *       Stored relationships between reconstructed photons and jets.
 *
 *   TruthVertices
 *       Geant4 truth-vertex census.
 *
 *   TruthPhotons
 *       Retained primary truth photons and their truth-level properties.
 *
 *   TruthJets
 *       Retained truth-jet collections.
 *
 *   PhotonTruthLinks
 *       Reconstructed-photon to truth-photon associations.
 *
 *   JetTruthLinks
 *       Reconstructed-jet to truth-jet associations.
 *
 *   WeightComponents
 *       Producer-owned event-weight components and their validity.
 *
 *   TriggerScalers
 *       Cumulative trigger-scaler snapshots over source-entry intervals.
 *
 *   TriggerRunInfo
 *       Run-level trigger names, prescales, and scaler totals.
 *
 * File metadata
 * -------------
 * "metadata"
 *     Written during initialization. Contains the production mode, source
 *     identity, software/configuration hashes, calibration bindings, storage
 *     policy, and other file-level provenance as key=value records.
 *
 * "configuration"
 *     The exact resolved TreeProduction configuration bytes used for the job.
 *
 * "completion"
 *     Written only during finalization. Records completion_status and terminal
 *     source/accounting counts. A file without this object was never finalized
 *     and must not be treated as a completed production product.
 */
#include "../PhotonJetTree.h"

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
// ============================================================================
// ROOT branch helpers
//
// These helpers are the only place where in-memory C++ types are translated
// into ROOT branch representations. Fixed-width scalars and persisted enums are
// given explicit leaf types so the file schema does not depend on compiler
// choices. Strings and vectors use ROOT's object-branch interface directly.
//
// Identity values are always serialized as two uint64 leaves:
//
//   <name>_hi
//   <name>_lo
//
// Together they preserve the 128-bit relational key used throughout the file.
// ============================================================================

void branchIdentity(TTree* tree, const std::string& name, photonjet::Identity& id)
{
  tree->Branch((name + "_hi").c_str(), &id.hi, (name + "_hi/l").c_str());
  tree->Branch((name + "_lo").c_str(), &id.lo, (name + "_lo/l").c_str());
}

template <typename Enum>
void branchEnum32(TTree* tree, const std::string& name, Enum& value)
{
  static_assert(sizeof(Enum) == sizeof(std::int32_t), "32-bit enumerator expected");
  tree->Branch(name.c_str(), reinterpret_cast<std::int32_t*>(&value), (name + "/I").c_str());
}

template <typename Enum>
void branchEnum8(TTree* tree, const std::string& name, Enum& value)
{
  static_assert(sizeof(Enum) == sizeof(std::uint8_t), "8-bit enumerator expected");
  tree->Branch(name.c_str(), reinterpret_cast<std::uint8_t*>(&value), (name + "/b").c_str());
}

void branchBool(TTree* tree, const std::string& name, bool& value)
{
  tree->Branch(name.c_str(), &value, (name + "/O").c_str());
}

void branchDouble(TTree* tree, const std::string& name, double& value)
{
  tree->Branch(name.c_str(), &value, (name + "/D").c_str());
}

void branchInt32(TTree* tree, const std::string& name, std::int32_t& value)
{
  tree->Branch(name.c_str(), &value, (name + "/I").c_str());
}

void branchInt(TTree* tree, const std::string& name, int& value)
{
  tree->Branch(name.c_str(), &value, (name + "/I").c_str());
}

void branchUInt32(TTree* tree, const std::string& name, std::uint32_t& value)
{
  tree->Branch(name.c_str(), &value, (name + "/i").c_str());
}

void branchInt64(TTree* tree, const std::string& name, std::int64_t& value)
{
  tree->Branch(name.c_str(), &value, (name + "/L").c_str());
}

void branchUInt64(TTree* tree, const std::string& name, std::uint64_t& value)
{
  tree->Branch(name.c_str(), &value, (name + "/l").c_str());
}

template <std::size_t N>
void branchUInt64Array(TTree* tree, const std::string& name, std::array<std::uint64_t, N>& value)
{
  tree->Branch(name.c_str(), value.data(), (name + "[" + std::to_string(N) + "]/l").c_str());
}

void branchString(TTree* tree, const std::string& name, std::string& value)
{
  tree->Branch(name.c_str(), &value);
}

template <typename T>
void branchVector(TTree* tree, const std::string& name, std::vector<T>& value)
{
  tree->Branch(name.c_str(), &value);
}


/**
 * Create one TTree owned by the output file.
 *
 * All canonical tables are booked through this helper so file ownership is
 * explicit and uniform.
 */
TTree* bookTree(TFile* file, const char* name, const char* title)
{
  file->cd();
  auto* tree = new TTree(name, title);
  tree->SetDirectory(file);
  return tree;
}


/**
 * Fill one already-prepared row and promote ROOT write failure to a producer
 * error. Serialization failure must never look like an ordinary missing row.
 */
void fillChecked(TTree* tree)
{
  if (!tree || tree->Fill() < 0) throw std::runtime_error("ROOT TTree fill failed");
}


// Metadata objects use a deliberately simple, human-readable key=value format.
void appendLine(std::ostringstream& out, const std::string& key, const std::string& value)
{
  out << key << '=' << value << '\n';
}

template <typename T>
void appendLine(std::ostringstream& out, const std::string& key, const T& value)
{
  out << key << '=' << value << '\n';
}

}  // namespace


// ============================================================================
// Output initialization
// ============================================================================

/**
 * Create the output product and establish its schema before event processing.
 *
 * CREATE mode deliberately refuses to overwrite an existing ROOT file. Once
 * the file exists, all canonical TTrees are booked and immutable file-level
 * provenance is written before the producer reports itself ready.
 */
void PhotonJetTree::initializeOutput()
{
  if (m_config.outputFile.empty())
  {
    throw std::runtime_error("no output file configured");
  }

  m_outputFile = TFile::Open(m_config.outputFile.c_str(), "CREATE");

  if (!m_outputFile || m_outputFile->IsZombie())
  {
    throw std::runtime_error("cannot create output file '" + m_config.outputFile + "'");
  }

  bookTrees();
  writeFileMetadata();
  m_outputReady = true;
}


// ============================================================================
// Canonical ROOT schema
//
// Each block below binds one record type to one relational TTree. The branch
// definitions are intentionally explicit: this function is the readable
// physical schema of the produced ROOT file.
//
// Object tables carry stable foreign-key identities rather than relying on
// entry ordering between TTrees.
// ============================================================================

void PhotonJetTree::bookTrees()
{
  // ==========================================================================
  // Sources
  //
  // One terminal row for the physical source processed by this invocation.
  // Stores source provenance, source range, and final exposure accounting.
  // ==========================================================================

  m_sourceTree = bookTree(m_outputFile, "Sources", "input sources and exposure accounting");
  {
    TTree* t = m_sourceTree;
    SourceRecord& r = m_sourceRow;
    branchIdentity(t, "source", r.sourceId);
    branchString(t, "dataset", r.dataset);
    branchString(t, "sample", r.sample);
    branchString(t, "period", r.period);
    branchString(t, "si_di_role", r.siDiRole);
    branchString(t, "input_uri_sha256", r.inputUriSha256);
    branchString(t, "input_file_sha256", r.inputFileSha256);
    branchString(t, "source_manifest_sha256", r.sourceManifestSha256);
    branchInt(t, "run", r.run);
    branchInt(t, "segment", r.segment);
    branchInt64(t, "source_file_ordinal", r.sourceFileOrdinal);
    branchInt64(t, "first_entry", r.firstEntry);
    branchInt64(t, "last_entry", r.lastEntry);
    branchUInt64(t, "encountered_events", r.encounteredEvents);
    branchUInt64(t, "retained_events", r.retainedEvents);
    branchUInt64(t, "upstream_rejected_events", r.upstreamRejectedEvents);
    branchBool(t, "completed", r.completed);
  }

  // ==========================================================================
  // UpstreamRejectedEvents
  //
  // Physical input events observed before reconstruction but prevented from
  // reaching PhotonJetTree::process_event(). These rows preserve the difference
  // between input exposure and producer encounters.
  // ==========================================================================

  m_upstreamRejectedTree = bookTree(m_outputFile, "UpstreamRejectedEvents",
                                    "events stopped by reconstruction before this module");
  {
    TTree* t = m_upstreamRejectedTree;
    UpstreamRejectedEventRecord& r = m_upstreamRejectedRow;
    branchIdentity(t, "source", r.sourceId);
    branchInt(t, "run", r.run);
    branchInt64(t, "physical_event_sequence", r.physicalEventSequence);
    branchInt64(t, "source_entry", r.sourceEntry);
  }

  // ==========================================================================
  // Events
  //
  // One row per producer encounter. This table is the event-level join point
  // for reconstructed objects, truth objects, weights, and event-quality
  // witnesses. Zero-photon and zero-jet events remain ordinary event rows.
  // ==========================================================================

  m_eventTree = bookTree(m_outputFile, "Events", "retained events");
  {
    TTree* t = m_eventTree;
    EventRecord& r = m_eventRow;
    branchIdentity(t, "source", r.sourceId);
    branchIdentity(t, "event", r.eventId);
    branchInt(t, "run", r.run);
    branchInt(t, "segment", r.segment);
    branchInt64(t, "source_file_ordinal", r.sourceFileOrdinal);
    branchInt64(t, "source_entry", r.sourceEntry);
    branchInt64(t, "source_global_entry", r.sourceGlobalEntry);
    branchUInt64(t, "producer_event_ordinal", r.producerEventOrdinal);
    branchInt64(t, "physical_event_sequence", r.physicalEventSequence);
    branchBool(t, "physical_event_sequence_valid", r.physicalEventSequenceValid);
    branchUInt64(t, "trigger_input_bits", r.triggerInputBits);
    branchUInt64(t, "trigger_live_bits", r.triggerLiveBits);
    branchUInt64(t, "trigger_scaled_bits", r.triggerScaledBits);
    branchUInt32(t, "trigger_packet_number", r.triggerPacketNumber);
    branchUInt64(t, "trigger_bunch_number", r.triggerBunchNumber);
    branchUInt64(t, "trigger_packet_status", r.triggerPacketStatus);
    branchInt32(t, "trigger_packet_version", r.triggerPacketVersion);
    branchUInt32(t, "trigger_decisions_available", r.triggerDecisionsAvailable);
    branchBool(t, "trigger_packet_valid", r.triggerPacketValid);
    branchInt64(t, "trigger_scaler_snapshot_id", r.triggerScalerSnapshotId);
    branchDouble(t, "reco_vertex_x", r.recoVertexX);
    branchDouble(t, "reco_vertex_y", r.recoVertexY);
    branchDouble(t, "reco_vertex_z", r.recoVertexZ);
    branchBool(t, "reco_vertex_valid", r.recoVertexValid);
    branchInt32(t, "reco_vertex_source", r.recoVertexSource);
    branchBool(t, "reco_object_vertex_in_domain", r.recoObjectVertexInDomain);
    branchDouble(t, "mbd_t0_ns", r.mbdT0Ns);
    branchDouble(t, "mbd_south_time_ns", r.mbdSouthTimeNs);
    branchDouble(t, "mbd_north_time_ns", r.mbdNorthTimeNs);
    branchDouble(t, "mbd_south_charge", r.mbdSouthCharge);
    branchDouble(t, "mbd_north_charge", r.mbdNorthCharge);
    branchDouble(t, "mbd_total_charge", r.mbdTotalCharge);
    branchBool(t, "mbd_valid", r.mbdValid);
    branchEnum32(t, "mbd_pmt_available", r.mbdPmtAvailable);
    branchVector(t, "mbd_pmt_id", r.mbdPmtId);
    branchVector(t, "mbd_pmt_arm", r.mbdPmtArm);
    branchVector(t, "mbd_pmt_charge", r.mbdPmtCharge);
    branchVector(t, "mbd_pmt_time_ns", r.mbdPmtTimeNs);
    branchVector(t, "mbd_pmt_valid", r.mbdPmtValid);
    branchBool(t, "centrality_applicable", r.centralityApplicable);
    branchEnum32(t, "minimum_bias_decision", r.minimumBiasDecision);
    branchDouble(t, "centrality_selected_charge", r.centralitySelectedCharge);
    branchInt32(t, "centrality_native_bin", r.centralityNativeBin);
    branchDouble(t, "centrality_native_percent", r.centralityNativePercent);
    branchBool(t, "centrality_native_valid", r.centralityNativeValid);
    branchInt32(t, "centrality_bin", r.centralityBin);
    branchDouble(t, "centrality_percent", r.centralityPercent);
    branchBool(t, "centrality_valid", r.centralityValid);
    branchDouble(t, "cemc_energy", r.cemcEnergy);
    branchDouble(t, "ihcal_energy", r.ihcalEnergy);
    branchDouble(t, "ohcal_energy", r.ohcalEnergy);
    branchDouble(t, "total_calo_energy", r.totalCaloEnergy);
    branchUInt32(t, "calo_available_mask", r.caloAvailableMask);
    branchUInt32(t, "calo_valid_mask", r.caloValidMask);
    branchBool(t, "calo_required_good_towers", r.caloRequiredGoodTowers);
    branchEnum32(t, "jet_background_state", r.jetBackgroundState);
    branchInt32(t, "jet_background_flow_mode", r.jetBackgroundFlowMode);
    branchDouble(t, "jet_background_v2", r.jetBackgroundV2);
    branchDouble(t, "jet_background_psi2", r.jetBackgroundPsi2);
    branchVector(t, "jet_background_ue_emcal", r.jetBackgroundUEEmcal);
    branchVector(t, "jet_background_ue_ihcal", r.jetBackgroundUEIhcal);
    branchVector(t, "jet_background_ue_ohcal", r.jetBackgroundUEOhcal);
    branchInt32(t, "jet_background_n_strips", r.jetBackgroundNStrips);
    branchInt32(t, "jet_background_n_towers", r.jetBackgroundNTowers);
    branchEnum32(t, "jet_background_flow_failure", r.jetBackgroundFlowFailure);
    branchVector(t, "reco_jet_view", r.recoJetView);
    branchVector(t, "reco_jet_radius_code", r.recoJetRadiusCode);
    branchVector(t, "reco_jet_capture_state", r.recoJetCaptureState);
    branchVector(t, "reco_jet_count", r.recoJetCount);
    branchUInt64Array(t, "calo_tower_count", r.caloTowerCount);
    branchUInt64Array(t, "calo_accepted_tower_count", r.caloAcceptedTowerCount);
    branchUInt64Array(t, "calo_bad_quality_tower_count", r.caloBadQualityTowerCount);
    branchUInt64Array(t, "calo_null_tower_count", r.caloNullTowerCount);
    branchUInt64Array(t, "calo_non_finite_tower_count", r.caloNonFiniteTowerCount);
    branchEnum32(t, "truth_vertex_capture_state", r.truthVertexCaptureState);
    branchInt32(t, "truth_primary_vertex_id", r.truthPrimaryVertexId);
    branchDouble(t, "truth_hard_vertex_z", r.truthHardVertexZ);
    branchBool(t, "truth_hard_vertex_valid", r.truthHardVertexValid);
    branchDouble(t, "truth_minimum_bias_vertex_z", r.truthMinimumBiasVertexZ);
    branchBool(t, "truth_minimum_bias_vertex_valid", r.truthMinimumBiasVertexValid);
    branchEnum32(t, "embedded_minimum_bias", r.embeddedMinimumBias);
    branchEnum32(t, "truth_photon_capture_state", r.truthPhotonCaptureState);
    branchUInt32(t, "truth_photon_embedded_primary_count", r.truthPhotonEmbeddedPrimaryCount);
    branchUInt32(t, "truth_photon_written_count", r.truthPhotonWrittenCount);
    branchUInt32(t, "truth_photon_rejected_count", r.truthPhotonRejectedCount);
    branchUInt32(t, "truth_photon_rejected_kinematics_count", r.truthPhotonRejectedKinematicsCount);
    branchUInt32(t, "truth_photon_rejected_isolation_count", r.truthPhotonRejectedIsolationCount);
    branchUInt32(t, "truth_photon_duplicate_track_count", r.truthPhotonDuplicateTrackCount);
    branchUInt32(t, "truth_photon_isolation_incomplete_count", r.truthPhotonIsolationIncompleteCount);
    branchUInt32(t, "truth_photon_analysis_signal_count", r.truthPhotonAnalysisSignalCount);
    branchBool(t, "truth_denominator_complete", r.truthDenominatorComplete);
    branchEnum32(t, "truth_jet_capture_state", r.truthJetCaptureState);
    branchVector(t, "truth_jet_radius_code", r.truthJetRadiusCode);
    branchVector(t, "truth_jet_container_valid", r.truthJetContainerValid);
    branchEnum32(t, "reco_photon_capture_state", r.recoPhotonCaptureState);
    branchUInt32(t, "reco_photon_unclassifiable_count", r.recoPhotonUnclassifiableCount);
    branchDouble(t, "event_weight", r.eventWeight);
    branchBool(t, "event_weight_valid", r.eventWeightValid);
    branchUInt32(t, "photon_count", r.photonCount);
    branchUInt32(t, "jet_count", r.jetCount);
    branchUInt32(t, "photon_jet_pair_count", r.photonJetPairCount);
    branchUInt32(t, "truth_photon_count", r.truthPhotonCount);
    branchUInt32(t, "truth_jet_count", r.truthJetCount);
    branchInt32(t, "terminal_status", r.terminalStatus);
  }

  // ==========================================================================
  // Photons
  //
  // One row per retained reconstructed photon candidate. Detailed shower views,
  // cells, and isolation are normalized into their own tables rather than
  // duplicating variable-sized payloads here.
  // ==========================================================================

  m_photonTree = bookTree(m_outputFile, "Photons", "reconstructed photon candidates");
  {
    TTree* t = m_photonTree;
    PhotonRecord& r = m_photonRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "photon", r.photonId);
    branchUInt64(t, "native_cluster_key", r.nativeClusterKey);
    branchUInt32(t, "encounter_ordinal", r.encounterOrdinal);
    branchDouble(t, "energy", r.energy);
    branchDouble(t, "et", r.et);
    branchDouble(t, "eta", r.eta);
    branchDouble(t, "phi", r.phi);
    branchBool(t, "kinematics_finite", r.kinematicsFinite);
    branchDouble(t, "producer_vertex_z", r.producerVertexZ);
    branchDouble(t, "timing_mean_time_samples", r.timing.meanTimeSamples);
    branchDouble(t, "timing_energy_weighted_numerator", r.timing.energyWeightedNumerator);
    branchDouble(t, "timing_energy_denominator", r.timing.energyDenominator);
    branchUInt32(t, "timing_contributing_tower_count", r.timing.contributingTowerCount);
    branchBool(t, "timing_finite", r.timing.finite);
    branchBool(t, "timing_valid", r.timing.valid);
    branchDouble(t, "timing_ns", r.timing.timeNs);
    branchEnum32(t, "dominant_truth_state", r.dominantTruth.state);
    branchEnum32(t, "dominant_truth_evaluator", r.dominantTruth.evaluator);
    branchInt32(t, "dominant_truth_track_id", r.dominantTruth.trackId);
    branchInt32(t, "dominant_truth_pid", r.dominantTruth.pid);
    branchInt32(t, "dominant_truth_barcode", r.dominantTruth.barcode);
    branchInt32(t, "dominant_truth_embedding_id", r.dominantTruth.embeddingId);
    branchInt32(t, "dominant_truth_vertex_id", r.dominantTruth.vertexId);
    branchDouble(t, "dominant_truth_energy_contribution", r.dominantTruth.energyContribution);
  }

  // ==========================================================================
  // PhotonShowerViews
  //
  // One row per photon and shower-definition variant. Event and photon
  // identities are carried alongside each view so consumers join explicitly
  // rather than relying on row ordering.
  // ==========================================================================

  m_showerShapeTree = bookTree(m_outputFile, "PhotonShowerViews", "shower-shape definitions per candidate");
  {
    TTree* t = m_showerShapeTree;

    // The row carries the owning identities beside the record so the table
    // joins without a parallel index.
    branchIdentity(t, "event", m_photonRow.eventId);
    branchIdentity(t, "photon", m_photonRow.photonId);

    ShowerShapeRecord& r = m_showerShapeRow;
    branchString(t, "definition", r.definitionName);
    branchDouble(t, "energy_floor_gev", r.energyFloorGeV);
    branchDouble(t, "e11", r.e11); branchDouble(t, "e13", r.e13); branchDouble(t, "e15", r.e15);
    branchDouble(t, "e17", r.e17); branchDouble(t, "e22", r.e22); branchDouble(t, "e31", r.e31);
    branchDouble(t, "e32", r.e32); branchDouble(t, "e33", r.e33); branchDouble(t, "e35", r.e35);
    branchDouble(t, "e37", r.e37); branchDouble(t, "e51", r.e51); branchDouble(t, "e52", r.e52);
    branchDouble(t, "e53", r.e53); branchDouble(t, "e55", r.e55); branchDouble(t, "e57", r.e57);
    branchDouble(t, "e71", r.e71); branchDouble(t, "e72", r.e72); branchDouble(t, "e73", r.e73);
    branchDouble(t, "e75", r.e75); branchDouble(t, "e77", r.e77);
    branchDouble(t, "weta", r.weta); branchDouble(t, "wphi", r.wphi);
    branchDouble(t, "weta_cog", r.wetaCog); branchDouble(t, "wphi_cog", r.wphiCog);
    branchDouble(t, "weta_cogx", r.wetaCogX); branchDouble(t, "wphi_cogx", r.wphiCogX);
    branchDouble(t, "weta33_cogx", r.weta33CogX); branchDouble(t, "wphi33_cogx", r.wphi33CogX);
    branchDouble(t, "w32", r.w32); branchDouble(t, "w52", r.w52); branchDouble(t, "w72", r.w72);
    branchDouble(t, "et1", r.et1); branchDouble(t, "et2", r.et2);
    branchDouble(t, "et3", r.et3); branchDouble(t, "et4", r.et4);
    branchInt32(t, "center_eta_index", r.centerEtaIndex);
    branchInt32(t, "center_phi_index", r.centerPhiIndex);
    branchDouble(t, "center_of_gravity_eta", r.centerOfGravityEta);
    branchDouble(t, "center_of_gravity_phi", r.centerOfGravityPhi);
    branchInt32(t, "delta_eta_max", r.deltaEtaMax);
    branchInt32(t, "delta_phi_max", r.deltaPhiMax);
    branchDouble(t, "delta_eta_cog", r.deltaEtaCog);
    branchDouble(t, "delta_phi_cog", r.deltaPhiCog);
    branchDouble(t, "radial_spread", r.radialSpread);
    branchUInt32(t, "saturated_tower_count", r.saturatedTowerCount);
    branchUInt32(t, "owned_cell_count", r.ownedCellCount);
    branchUInt32(t, "good_cell_count", r.goodCellCount);
    branchUInt32(t, "zero_cell_count", r.zeroCellCount);
    branchUInt32(t, "negative_cell_count", r.negativeCellCount);
    branchUInt32(t, "non_finite_cell_count", r.nonFiniteCellCount);
    branchBool(t, "valid", r.valid);
  }

  // ==========================================================================
  // PhotonCells
  //
  // Calorimeter-cell witnesses around retained photon candidates. These rows
  // preserve the low-level inputs needed to audit or reproduce shower-derived
  // quantities without embedding them repeatedly in the photon table.
  // ==========================================================================

  m_photonCellTree = bookTree(m_outputFile, "PhotonCells", "calorimeter cells around each candidate");
  {
    TTree* t = m_photonCellTree;
    PhotonCellRecord& r = m_photonCellRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "photon", r.photonId);
    branchUInt64(t, "tower_key", r.towerKey);
    branchInt32(t, "subsystem", r.subsystem);
    branchInt32(t, "eta_index", r.etaIndex);
    branchInt32(t, "phi_index", r.phiIndex);
    branchDouble(t, "calibrated_energy", r.calibratedEnergy);
    branchDouble(t, "cluster_map_energy", r.clusterMapEnergy);
    branchDouble(t, "tower_time_samples", r.towerTimeSamples);
    branchUInt32(t, "tower_status", r.towerStatus);
    branchBool(t, "owned_by_cluster", r.ownedByCluster);
    branchBool(t, "in_shower_support", r.inShowerSupport);
    branchBool(t, "good", r.good);
    branchBool(t, "zero_energy", r.zeroEnergy);
    branchBool(t, "negative_energy", r.negativeEnergy);
    branchBool(t, "finite", r.finite);
  }

  // ==========================================================================
  // Isolation
  //
  // Reconstructed isolation measurements keyed to the photon, cone radius, and
  // method. These are observables only; isolated/non-isolated decisions are
  // deliberately downstream.
  // ==========================================================================

  m_isolationTree = bookTree(m_outputFile, "Isolation", "isolation cones per candidate");
  {
    TTree* t = m_isolationTree;
    IsolationRecord& r = m_isolationRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "photon", r.photonId);
    branchDouble(t, "axis_eta", r.axisEta);
    branchDouble(t, "axis_phi", r.axisPhi);
    branchDouble(t, "radius", r.radius);
    branchEnum32(t, "method", r.method);
    branchDouble(t, "cone_sum", r.coneSum);
    branchDouble(t, "electromagnetic_component", r.electromagneticComponent);
    branchDouble(t, "inner_hadronic_component", r.innerHadronicComponent);
    branchDouble(t, "outer_hadronic_component", r.outerHadronicComponent);
    branchBool(t, "candidate_removed", r.candidateRemoved);
    branchBool(t, "valid", r.valid);
  }

  // ==========================================================================
  // IsolationConstituents
  //
  // Optional constituent-level witnesses underlying reconstructed isolation.
  // They preserve source, geometry, energy/subtraction state, masking, and
  // candidate-removal information for later auditing or recomputation.
  // ==========================================================================

  m_isolationConstituentTree = bookTree(m_outputFile, "IsolationConstituents", "objects entering isolation cones");
  {
    TTree* t = m_isolationConstituentTree;
    IsolationConstituentRecord& r = m_isolationConstituentRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "photon", r.photonId);
    branchDouble(t, "radius", r.radius);
    branchString(t, "source", r.source);
    branchUInt64(t, "native_key", r.nativeKey);
    branchInt32(t, "subsystem", r.subsystem);
    branchDouble(t, "delta_eta", r.deltaEta);
    branchDouble(t, "delta_phi", r.deltaPhi);
    branchDouble(t, "delta_r", r.deltaR);
    branchDouble(t, "raw_energy", r.rawEnergy);
    branchDouble(t, "transverse_energy", r.transverseEnergy);
    branchDouble(t, "subtracted_energy", r.subtractedEnergy);
    branchDouble(t, "subtracted_transverse_energy", r.subtractedTransverseEnergy);
    branchInt32(t, "quality_state", r.qualityState);
    branchBool(t, "masked", r.masked);
    branchUInt32(t, "mask_state", r.maskState);
    branchBool(t, "candidate_removed", r.candidateRemoved);
  }

  // ==========================================================================
  // Jets
  //
  // Reconstructed jets from every configured view/radius. Raw and corrected
  // momenta, collection identities, area, ordering, and calibration witnesses
  // are stored together in each row.
  // ==========================================================================

  m_jetTree = bookTree(m_outputFile, "Jets", "reconstructed jets");
  {
    TTree* t = m_jetTree;
    JetRecord& r = m_jetRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "jet", r.jetId);
    branchEnum8(t, "view", r.view);
    branchDouble(t, "radius", r.radius);
    branchString(t, "input_identity", r.inputIdentity);
    branchString(t, "subtraction_identity", r.subtractionIdentity);
    branchUInt64(t, "native_key", r.nativeKey);
    branchUInt64(t, "native_raw_key", r.nativeRawKey);
    branchUInt32(t, "encounter_ordinal", r.encounterOrdinal);
    branchUInt32(t, "deterministic_order", r.deterministicOrder);
    branchDouble(t, "raw_pt", r.rawPt);
    branchDouble(t, "corrected_pt", r.correctedPt);
    branchDouble(t, "eta", r.eta);
    branchDouble(t, "phi", r.phi);
    branchDouble(t, "area", r.area);
    branchUInt64(t, "quality_bitmask", r.qualityBitmask);
    branchBool(t, "calibration_valid", r.calibrationValid);
  }

  // ==========================================================================
  // PhotonJetPairs
  //
  // Explicit reconstructed photon-jet relations. The table stores the
  // relationship observables and recoil witness without imposing the final
  // recoil selection on the underlying object inventory.
  // ==========================================================================

  m_photonJetPairTree = bookTree(m_outputFile, "PhotonJetPairs", "photon and jet pairs");
  {
    TTree* t = m_photonJetPairTree;
    PhotonJetPairRecord& r = m_photonJetPairRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "pair", r.pairId);
    branchIdentity(t, "photon", r.photonId);
    branchIdentity(t, "jet", r.jetId);
    branchEnum8(t, "jet_view", r.jetView);
    branchDouble(t, "jet_radius", r.jetRadius);
    branchDouble(t, "delta_phi", r.deltaPhi);
    branchDouble(t, "xj_gamma", r.xJGamma);
    branchBool(t, "recoil_witness", r.recoilWitness);
    branchUInt32(t, "photon_rank", r.photonRank);
    branchUInt32(t, "jet_rank", r.jetRank);
  }

  // ==========================================================================
  // TruthVertices
  //
  // Geant4 truth-vertex census, including embedding identity and validity.
  // ==========================================================================

  m_truthVertexTree = bookTree(m_outputFile, "TruthVertices", "Geant4 vertex census");
  {
    TTree* t = m_truthVertexTree;
    TruthVertexRecord& r = m_truthVertexRow;
    branchIdentity(t, "event", r.eventId);
    branchInt32(t, "vertex_id", r.vertexId);
    branchInt32(t, "embedding_id", r.embeddingId);
    branchBool(t, "embedding_valid", r.embeddingValid);
    branchDouble(t, "z", r.z);
    branchBool(t, "valid", r.valid);
  }

  // ==========================================================================
  // TruthPhotons
  //
  // Retained primary truth photons with generator/embedding provenance,
  // classification, truth isolation, and analysis-signal witness.
  // ==========================================================================

  m_truthPhotonTree = bookTree(m_outputFile, "TruthPhotons", "embedded primary photons");
  {
    TTree* t = m_truthPhotonTree;
    TruthPhotonRecord& r = m_truthPhotonRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "truth_photon", r.truthPhotonId);
    branchInt32(t, "track_id", r.trackId);
    branchInt32(t, "vertex_id", r.vertexId);
    branchInt32(t, "barcode", r.barcode);
    branchInt32(t, "embedding_id", r.embeddingId);
    branchInt32(t, "pid", r.pid);
    branchDouble(t, "pt", r.pt);
    branchDouble(t, "eta", r.eta);
    branchDouble(t, "phi", r.phi);
    branchBool(t, "geant_valid", r.geantValid);
    branchBool(t, "generator_association_valid", r.generatorAssociationValid);
    branchEnum32(t, "prompt_class", r.promptClass);
    branchDouble(t, "isolation_r03", r.isolationR03);
    branchDouble(t, "isolation_r04", r.isolationR04);
    branchBool(t, "isolation_valid", r.isolationValid);
    branchBool(t, "analysis_signal", r.analysisSignal);
    branchEnum8(t, "sample_role", r.sampleRole);
    branchBool(t, "signal_source_role", r.signalSourceRole);
  }

  // ==========================================================================
  // TruthJets
  //
  // Truth jets remain independent of reconstructed JES. Radius and native key
  // identify the originating truth collection/object.
  // ==========================================================================

  m_truthJetTree = bookTree(m_outputFile, "TruthJets", "truth jets");
  {
    TTree* t = m_truthJetTree;
    TruthJetRecord& r = m_truthJetRow;
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "truth_jet", r.truthJetId);
    branchDouble(t, "radius", r.radius);
    branchUInt64(t, "native_key", r.nativeKey);
    branchDouble(t, "pt", r.pt);
    branchDouble(t, "eta", r.eta);
    branchDouble(t, "phi", r.phi);
    branchBool(t, "container_valid", r.containerValid);
  }

  // ==========================================================================
  // PhotonTruthLinks
  //
  // Explicit reco-photon <-> truth-photon relations. Relation state and
  // association witnesses remain separate from the two object tables.
  // ==========================================================================

  m_photonTruthLinkTree = bookTree(m_outputFile, "PhotonTruthLinks", "reconstructed photon to truth photon relations");
  {
    TTree* t = m_photonTruthLinkTree;
    PhotonTruthLinkRecord& r = m_photonTruthLinkRow;
    branchIdentity(t, "link", r.linkId);
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "photon", r.photonId);
    branchIdentity(t, "truth_photon", r.truthPhotonId);
    branchEnum32(t, "state", r.state);
    branchDouble(t, "energy_contribution", r.energyContribution);
    branchDouble(t, "delta_r", r.deltaR);
  }

  // ==========================================================================
  // JetTruthLinks
  //
  // Explicit reco-jet <-> truth-jet relations, retaining jet view/radius and
  // whether a candidate relation was selected by the matching algorithm.
  // ==========================================================================

  m_jetTruthLinkTree = bookTree(m_outputFile, "JetTruthLinks", "reconstructed jet to truth jet relations");
  {
    TTree* t = m_jetTruthLinkTree;
    JetTruthLinkRecord& r = m_jetTruthLinkRow;
    branchIdentity(t, "link", r.linkId);
    branchIdentity(t, "event", r.eventId);
    branchIdentity(t, "jet", r.jetId);
    branchIdentity(t, "truth_jet", r.truthJetId);
    branchEnum8(t, "jet_view", r.jetView);
    branchDouble(t, "jet_radius", r.jetRadius);
    branchEnum32(t, "state", r.state);
    branchDouble(t, "delta_r", r.deltaR);
    branchBool(t, "selected_match", r.selectedMatch);
  }

  // ==========================================================================
  // WeightComponents
  //
  // Producer-owned event-weight factors and their validity/application count.
  // Downstream sample stitching and final normalization are not invented here.
  // ==========================================================================

  m_weightTree = bookTree(m_outputFile, "WeightComponents", "producer weight factors");
  {
    TTree* t = m_weightTree;
    WeightComponentRecord& r = m_weightRow;
    branchIdentity(t, "event", r.eventId);
    branchString(t, "component_type", r.componentType);
    branchDouble(t, "slice_weight", r.sliceWeight);
    branchDouble(t, "cross_section_weight", r.crossSectionWeight);
    branchDouble(t, "vertex_weight", r.vertexWeight);
    branchDouble(t, "si_di_weight", r.siDiWeight);
    branchDouble(t, "period_weight", r.periodWeight);
    branchDouble(t, "exposure_weight", r.exposureWeight);
    branchDouble(t, "final_weight", r.finalWeight);
    branchBool(t, "valid", r.valid);
    branchInt32(t, "application_count", r.applicationCount);
  }

  // ==========================================================================
  // TriggerScalers
  //
  // Run-length-like snapshots of cumulative raw/live/scaled GL1 counters.
  // Source-entry and bunch ranges identify the interval represented by a row.
  // ==========================================================================

  m_triggerScalerTree = bookTree(m_outputFile, "TriggerScalers", "cumulative trigger scaler snapshots");
  {
    TTree* t = m_triggerScalerTree;
    TriggerScalerSnapshotRecord& r = m_triggerScalerRow;
    branchIdentity(t, "source", r.sourceId);
    branchInt64(t, "snapshot_id", r.snapshotId);
    branchUInt64(t, "bunch_start", r.bunchStart);
    branchUInt64(t, "bunch_end", r.bunchEnd);
    branchInt64(t, "source_entry_start", r.sourceEntryStart);
    branchInt64(t, "source_entry_end", r.sourceEntryEnd);
    branchUInt64(t, "observed_events", r.observedEvents);
    branchUInt64Array(t, "raw", r.raw);
    branchUInt64Array(t, "live", r.live);
    branchUInt64Array(t, "scaled", r.scaled);
    branchBool(t, "valid", r.valid);
    branchBool(t, "discontinuity", r.discontinuity);
  }

  // ==========================================================================
  // TriggerRunInfo
  //
  // One row per trigger bit containing run-level trigger names, prescales, and
  // cumulative counts. Missing/unsupported configuration is represented by row
  // validity rather than by silently omitting the table.
  // ==========================================================================

  m_triggerRunInfoTree = bookTree(m_outputFile, "TriggerRunInfo", "run-level trigger configuration");
  {
    TTree* t = m_triggerRunInfoTree;
    TriggerRunInfoRecord& r = m_triggerRunInfoRow;
    branchIdentity(t, "source", r.sourceId);
    branchInt(t, "run", r.run);
    branchInt32(t, "bit", r.bit);
    branchString(t, "name", r.name);
    branchDouble(t, "hardware_prescale", r.hardwarePrescale);
    branchDouble(t, "run_average_prescale", r.runAveragePrescale);
    branchUInt64(t, "raw_count", r.rawCount);
    branchUInt64(t, "live_count", r.liveCount);
    branchUInt64(t, "scaled_count", r.scaledCount);
    branchBool(t, "valid", r.valid);
  }
}


// ============================================================================
// File-level provenance
// ============================================================================

/**
 * Persist everything needed to identify how this ROOT product was produced.
 *
 * "metadata" is intentionally compact and human-readable. It binds production
 * mode, source identity, software/configuration hashes, calibration payloads,
 * storage policy, and collection definitions.
 *
 * The resolved configuration itself is written separately as "configuration"
 * so the output remains self-describing without reconstructing settings from
 * the metadata summary.
 */
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

  m_outputFile->cd();

  TObjString metadata(text.str().c_str());
  if (metadata.Write("metadata", TObject::kOverwrite) <= 0)
    throw std::runtime_error("ROOT metadata write failed");

  // Persist the exact resolved configuration bytes beside the summary metadata.
  TObjString configuration(m_config.provenance.configurationText.c_str());
  if (configuration.Write("configuration", TObject::kOverwrite) <= 0)
    throw std::runtime_error("ROOT configuration write failed");
}


// ============================================================================
// Per-event serialization
// ============================================================================

/**
 * Finalize event-level counts and write the single Events row.
 *
 * Object counts are taken from the completed in-memory collections immediately
 * before serialization so the event row describes exactly what will be written
 * to the corresponding object tables.
 */
void PhotonJetTree::fillEventOutput()
{
  m_event.sourceId = m_source.sourceId;
  m_event.photonCount = static_cast<std::uint32_t>(m_photons.size());
  m_event.jetCount = static_cast<std::uint32_t>(m_jets.size());
  m_event.photonJetPairCount = static_cast<std::uint32_t>(m_photonJetPairs.size());
  m_event.truthPhotonCount = static_cast<std::uint32_t>(m_truthPhotons.size());
  m_event.truthJetCount = static_cast<std::uint32_t>(m_truthJets.size());
  m_eventRow = m_event;
  fillChecked(m_eventTree);
}


/**
 * Serialize every object/relation record accumulated for the current event.
 *
 * Each record is copied into the persistent row buffer whose address is bound
 * to ROOT, then TTree::Fill() snapshots that row. No physics calculation is
 * introduced by this function.
 */
void PhotonJetTree::fillObjectOutput()
{
  // Reconstructed photons and their shower-definition views.
  for (const PhotonRecord& photon : m_photons)
  {
    m_photonRow = photon;
    fillChecked(m_photonTree);

    for (const ShowerShapeRecord& shape : photon.showerShapes)
    {
      m_showerShapeRow = shape;
      fillChecked(m_showerShapeTree);
    }
  }

  // Photon-associated calorimeter and isolation witnesses.
  for (const PhotonCellRecord& cell : m_photonCells)
  {
    m_photonCellRow = cell;
    fillChecked(m_photonCellTree);
  }

  for (const IsolationRecord& cone : m_isolationRecords)
  {
    m_isolationRow = cone;
    fillChecked(m_isolationTree);
  }

  for (const IsolationConstituentRecord& constituent : m_isolationConstituents)
  {
    m_isolationConstituentRow = constituent;
    fillChecked(m_isolationConstituentTree);
  }

  // Reconstructed jets and photon-jet relationships.
  for (const JetRecord& jet : m_jets)
  {
    m_jetRow = jet;
    fillChecked(m_jetTree);
  }

  for (const PhotonJetPairRecord& pair : m_photonJetPairs)
  {
    m_photonJetPairRow = pair;
    fillChecked(m_photonJetPairTree);
  }

  // Simulation truth inventory.
  for (const TruthVertexRecord& vertex : m_truthVertices)
  {
    m_truthVertexRow = vertex;
    fillChecked(m_truthVertexTree);
  }

  for (const TruthPhotonRecord& photon : m_truthPhotons)
  {
    m_truthPhotonRow = photon;
    fillChecked(m_truthPhotonTree);
  }

  for (const TruthJetRecord& jet : m_truthJets)
  {
    m_truthJetRow = jet;
    fillChecked(m_truthJetTree);
  }

  // Reconstruction-to-truth relations.
  for (const PhotonTruthLinkRecord& link : m_photonTruthLinks)
  {
    m_photonTruthLinkRow = link;
    fillChecked(m_photonTruthLinkTree);
  }

  for (const JetTruthLinkRecord& link : m_jetTruthLinks)
  {
    m_jetTruthLinkRow = link;
    fillChecked(m_jetTruthLinkTree);
  }

  // Producer-owned event-weight decomposition.
  for (const WeightComponentRecord& weight : m_weightComponents)
  {
    m_weightRow = weight;
    fillChecked(m_weightTree);
  }
}


// ============================================================================
// Final source-level serialization and completion certificate
// ============================================================================

/**
 * Flush source/run-lifetime tables and write the terminal completion object.
 *
 * Event/object rows are written during event processing. Source accounting,
 * upstream rejections, scaler snapshots, trigger-run configuration, and the
 * final source record become authoritative only after the input has finished,
 * so they are serialized here.
 *
 * The ROOT file is flushed before "completion" is written. Presence of that
 * object therefore acts as the terminal product marker; its status still
 * distinguishes a successful production from an explicitly aborted one.
 */
void PhotonJetTree::writeCompletionMetadata()
{
  if (!m_outputFile)
  {
    return;
  }

  // --------------------------------------------------------------------------
  // Source/run-lifetime tables
  // --------------------------------------------------------------------------

  for (const UpstreamRejectedEventRecord& rejected : m_upstreamRejectedEvents)
  {
    m_upstreamRejectedRow = rejected;
    fillChecked(m_upstreamRejectedTree);
  }

  for (const TriggerScalerSnapshotRecord& snapshot : m_triggerScalerSnapshots)
  {
    m_triggerScalerRow = snapshot;
    fillChecked(m_triggerScalerTree);
  }

  for (const TriggerRunInfoRecord& info : m_triggerRunInfo)
  {
    m_triggerRunInfoRow = info;
    fillChecked(m_triggerRunInfoTree);
  }

  /*
   * Finalize the source record only after all exposure accounting is known.
   * completed describes producer status, not merely whether a ROOT file exists.
   */
  m_source.upstreamRejectedEvents = static_cast<std::uint64_t>(m_upstreamRejectedEvents.size());
  m_source.completed = !m_aborted;
  m_sourceRow = m_source;
  fillChecked(m_sourceTree);

  /*
   * Flush every TTree before declaring terminal completion. A ROOT write error
   * here prevents creation of a completion certificate.
   */
  if (m_outputFile->Write(nullptr, TObject::kOverwrite) <= 0 || m_outputFile->TestBit(TFile::kWriteError))
    throw std::runtime_error("ROOT flush failed before completion");

  // --------------------------------------------------------------------------
  // Terminal accounting summary
  // --------------------------------------------------------------------------

  std::ostringstream text;
  appendLine(text, "completion_status", std::string(m_aborted ? "aborted" : "complete"));
  appendLine(text, "encountered_events", m_source.encounteredEvents);
  appendLine(text, "retained_events", m_source.retainedEvents);
  appendLine(text, "upstream_rejected_events", m_source.upstreamRejectedEvents);
  appendLine(text, "upstream_observed_events", m_upstreamObservedEvents);
  appendLine(text, "first_entry", m_source.firstEntry);
  appendLine(text, "last_entry", m_source.lastEntry);
  appendLine(text, "trigger_scaler_snapshots", m_triggerScalerSnapshots.size());
  appendLine(text, "trigger_run_info_rows", m_triggerRunInfo.size());

  m_outputFile->cd();

  TObjString completion(text.str().c_str());
  if (completion.Write("completion", TObject::kOverwrite) <= 0)
    throw std::runtime_error("ROOT completion write failed");
}


/**
 * Close the ROOT product and invalidate all file-owned handles.
 *
 * The TTrees have already been flushed before completion metadata was written,
 * so this function does not intentionally rewrite the schema after the product
 * has been declared complete. ROOT's write-error state is still checked after
 * closing and propagated as a producer failure.
 */
void PhotonJetTree::closeOutput()
{
  if (!m_outputFile)
  {
    return;
  }

  m_outputFile->cd();

  // Trees were flushed before the completion object; do not rewrite them
  // after declaring completion. Close still checks ROOT's write-error bit.
  m_outputFile->Close();

  const bool writeError = m_outputFile->TestBit(TFile::kWriteError);

  delete m_outputFile;
  m_outputFile = nullptr;

  /*
   * All TTrees were owned by the TFile and were released with it. Clearing the
   * non-owning pointers prevents accidental use after close.
   */
  m_sourceTree = m_upstreamRejectedTree = m_eventTree = nullptr;
  m_photonTree = m_showerShapeTree = m_photonCellTree = nullptr;
  m_isolationTree = m_isolationConstituentTree = nullptr;
  m_jetTree = m_photonJetPairTree = nullptr;
  m_truthVertexTree = m_truthPhotonTree = m_truthJetTree = nullptr;
  m_photonTruthLinkTree = m_jetTruthLinkTree = nullptr;
  m_weightTree = m_triggerScalerTree = m_triggerRunInfoTree = nullptr;

  if (writeError) throw std::runtime_error("ROOT close failed");
}
