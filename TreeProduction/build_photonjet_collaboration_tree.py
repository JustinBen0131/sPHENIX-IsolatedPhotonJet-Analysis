#!/usr/bin/env python3
"""Build a collaborator-facing TTree package from the internal ROOT files.

The input files remain authoritative and unchanged. The normalized_v1 layout
stores each event/object payload once and preserves unscored candidates. The
legacy expanded view remains available for existing consumers; unlike raw
capture, that historical model view can omit unscorable guard-range candidates.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
from pathlib import Path
from typing import Any, Iterable

import awkward as ak
import numpy as np
import uproot

import centrality_replay
import trigger_scaler_interface
import gl1_trigger_interface


PREFIX = "ReplayFoundationV1/"
UINT64_ZERO = np.uint64(0)
MODEL_ET_MIN_GEV = 15.0
MODEL_ET_MAX_GEV = 35.0
ROW_BATCH_SIZE = 2_000
EVENT_BATCH_SIZE = 512
LAYOUTS = ("expanded", "normalized_v1")
JET_MATCH_DELTA_R_MAX = 0.3


EVENT_SCHEMA = {
    "source_file_index": "int32",
    "source_entry": "int64",
    "event_id_hi": "uint64",
    "event_id_lo": "uint64",
    "run": "int32",
    "event_sequence": "int64",
    "physical_event_sequence": "int64",
    "trigger_bits": "uint64",
    "live_trigger_bits": "uint64",
    "scaled_trigger_bits": "uint64",
    "scaled_bit30": "int32",
    "vertex_z": "float64",
    "centrality": "float64",
    "event_weight": "float64",
    "total_calo_energy": "float64",
    "terminal_status": "int32",
    "trigger_scaler_snapshot_id": "uint64",
}

# Additive collaborator interface; the upstream schema number is independent.
LEGACY_EVENT_SCHEMA = dict(EVENT_SCHEMA)
EVENT_SCHEMA.update(gl1_trigger_interface.EVENT_SCHEMA)
EVENT_SCHEMA["physical_event_sequence_valid"] = "int32"
PPG12_EVENT_WITNESS_SCHEMA = {
    "source_file_ordinal": "int64", "source_entry_ordinal": "int64",
    "source_global_entry_ordinal": "int64", "ppg12_source_arbitration_state": "int32",
    "ppg12_leading_source_cluster_key": "int64",
    "ppg12_leading_source_encounter_ordinal": "int32",
    "ppg12_leading_source_et": "float64", "ppg12_leading_source_eta": "float64",
    "ppg12_leading_source_phi": "float64",
}
EVENT_SCHEMA.update(PPG12_EVENT_WITNESS_SCHEMA)
TRUTH_INVENTORY_SCHEMA = {"truth_photon_"+name: "int32" for name in (
    "capture_state", "native_embedded_primary_pid22_count", "serialized_raw_count",
    "rejected_count", "rejected_invalid_kinematics_count",
    "rejected_nonfinite_isolation_count", "duplicate_track_count",
    "isolation_input_incomplete_count", "analysis_signal_count")}
WEIGHT_FIELDS = (
    "slice_weight", "cross_section_weight", "vertex_weight", "si_di_weight",
    "period_weight", "exposure_weight", "final_weight",
)
INTERFACE_SCHEMA = {
    "collaborator_interface_version": "int32",
    "source_occurrence_id_hi": "uint64", "source_occurrence_id_lo": "uint64",
    "source_occurrence_id_available": "int32",
    "truth_vertex_z": "float64", "truth_mb_vertex_z": "float64",
    "truth_vertex_z_finite": "int32", "truth_mb_vertex_z_finite": "int32",
    "truth_vertex_capture_state": "int32", "truth_primary_vertex_id": "int32",
    "truth_jet_capture_state": "int32",
    "reco_vertex_z": "float64", "reco_vertex_valid": "int32",
    "reco_vertex_abs_lt_30": "int32", "reco_vertex_abs_lt_60": "int32",
    "mbd_total_charge": "float64", "mbd_total_charge_finite": "int32",
    "mbd_pmt_available": "int32", "mbd_pmt_sums_valid": "int32",
    "mbd_pmt_arm_from_channel_map": "int32",
    "mbd_pmt_sum_charge": "float64", "mbd_pmt_south_charge": "float64",
    "mbd_pmt_north_charge": "float64",
    "nodes_ready": "int32", "ppg12_weight_active_state": "int32",
    "ppg12_weight_lane_component_code": "int32",
    "sample_weight": "float64", "sample_weight_valid": "int32",
    "vertex_weight_valid": "int32", "truth_denominator_complete": "int32",
    "event_weight_without_vertex": "float64",
    "legacy_event_weight": "float64",
    "weight_components_recorded": "int32", "weight_identity_valid": "int32",
    "weight_application_count": "int32",
    **{("weight_final" if name == "final_weight" else name): "float64"
       for name in WEIGHT_FIELDS},
}
INTERFACE_SCHEMA.update(centrality_replay.SCALARS)
INTERFACE_SCHEMA.update(TRUTH_INVENTORY_SCHEMA)
EVENT_SCHEMA.update(INTERFACE_SCHEMA)
MBD_ARRAY_SCHEMA = {
    "mbd_pmt_id": "var * int32", "mbd_pmt_arm": "var * int32",
    "mbd_pmt_charge": "var * float64", "mbd_pmt_charge_valid": "var * int32",
}
TRUTH_VERTEX_ARRAY_SCHEMA = {
    "truth_vertex_id": "var * int32", "truth_vertex_embedding_id": "var * int32",
    "truth_vertex_embedding_valid": "var * int32", "truth_vertex_z_valid": "var * int32",
    "truth_vertex_z_values": "var * float64",
}
TRUTH_JET_AVAILABILITY_SCHEMA = {
    "truth_jet_radius_code": "var * int32", "truth_jet_container_valid": "var * int32",
}
RECO_JET_VIEW_AVAILABILITY_INPUT_BRANCHES = (
    "reco_jet_view_radius_code", "reco_jet_view_input_identity",
    "reco_jet_view_subtraction_identity", "reco_jet_view_available",
)
JET_VIEW_AVAILABILITY_SCHEMA = {
    "source_file_index": "int32", "source_entry": "int64",
    "event_id_hi": "uint64", "event_id_lo": "uint64",
    "reco_jet_view_radius_code": "int32",
    "reco_jet_view_input_identity": "string",
    "reco_jet_view_subtraction_identity": "string",
    "reco_jet_view_available": "int32",
}

PHOTON_FIELDS = {
    "candidate_id_hi": "uint64",
    "candidate_id_lo": "uint64",
    "photon_encounter_ordinal": "int32",
    "native_cluster_key": "int64",
    "ppg12_tower_mask_state": "int32",
    "ppg12_source_eligible_state": "int32",
    "reference_preselection_state": "int32",
    "active_preselection_state": "int32",
    "photon_et": "float64",
    "photon_eta": "float64",
    "photon_phi": "float64",
    "bdt_score": "float64",
    # 0: unscored, 1: evaluated with the supplied model, 2: source-recorded.
    "bdt_evaluation_state": "int32",
    "bdt_tight_threshold": "float64",
    "bdt_nontight_low_threshold": "float64",
    "bdt_nontight_high_threshold": "float64",
    "bdt_is_tight": "int32",
    "bdt_is_nontight": "int32",
    "bdt_is_not_tight": "int32",
    "bdt_input_count": "int32",
    "iso_r03": "float64",
    "iso_r03_threshold": "float64",
    "iso_r03_nonisolated_threshold": "float64",
    "iso_r03_pass": "int32",
    "iso_r04": "float64",
    "iso_r04_threshold": "float64",
    "iso_r04_nonisolated_threshold": "float64",
    "iso_r04_pass": "int32",
    "truth_matched": "int32",
    "truth_barcode": "int32",
    "dominant_truth_state": "int32",
    "dominant_truth_evaluator_mode": "int32",
    "dominant_truth_track_id": "int32",
    "dominant_truth_pid": "int32",
    "dominant_truth_barcode": "int32",
    "dominant_truth_embedding_id": "int32",
    "dominant_truth_energy_contribution": "float64",
    "dominant_truth_vertex_id": "int32",
    "native_weta_cogx": "float64",
    "native_wphi_cogx": "float64",
    "native_weta33_cogx": "float64",
    "native_wphi33_cogx": "float64",
    "native_e11_over_e33": "float64",
    "native_e32_over_e35": "float64",
    "native_et1": "float64",
}
for _index in range(14):
    PHOTON_FIELDS[f"bdt_input_{_index:02d}"] = "float64"

JET_FIELDS = {
    "jet_id_hi": "uint64",
    "jet_id_lo": "uint64",
    "jet_radius": "float64",
    "jet_raw_pt": "float64",
    "jet_pt": "float64",
    "jet_eta": "float64",
    "jet_phi": "float64",
    "jet_mass": "float64",
    "jet_area": "float64",
    "jet_quality_bitmask": "uint64",
    "jet_order": "int32",
    "jet_input_identity": "string",
    "jet_subtraction_identity": "string",
    # 1: exact producer labels, 0: legacy source without view evidence.
    "jet_view_provenance_state": "int32",
}

PAIR_FIELDS = {
    "pair_id_hi": "uint64",
    "pair_id_lo": "uint64",
    "photon_index": "int32",
    "jet_index": "int32",
    "delta_phi": "float64",
    "xjgamma": "float64",
    "recoil_state": "int32",
    "photon_rank": "int32",
    "jet_rank": "int32",
    "wrong_photon_class": "int32",
    "wrong_recoil_class": "int32",
}


# A familiar event-oriented view.  The flat object and pair TTrees remain
# available for direct loops; this tree collects the same rows into arrays.
EVENT_ARRAY_SCHEMA = {
    "nphotons": "int32",
    "njets": "int32",
    "npairs": "int32",
    "photon_candidate_id_hi": "var * uint64",
    "photon_candidate_id_lo": "var * uint64",
    "photon_encounter_ordinal": "var * int32",
    "photon_et": "var * float64",
    "photon_eta": "var * float64",
    "photon_phi": "var * float64",
    "photon_bdt_score": "var * float64",
    "photon_bdt_tight_threshold": "var * float64",
    "photon_bdt_nontight_low_threshold": "var * float64",
    "photon_bdt_nontight_high_threshold": "var * float64",
    "photon_bdt_is_tight": "var * int32",
    "photon_bdt_is_nontight": "var * int32",
    "photon_bdt_is_not_tight": "var * int32",
    "photon_iso_r03": "var * float64",
    "photon_iso_r03_threshold": "var * float64",
    "photon_iso_r03_nonisolated_threshold": "var * float64",
    "photon_iso_r04": "var * float64",
    "photon_iso_r04_threshold": "var * float64",
    "photon_iso_r04_nonisolated_threshold": "var * float64",
    "jet_id_hi": "var * uint64",
    "jet_id_lo": "var * uint64",
    "jet_pt": "var * float64",
    "jet_raw_pt": "var * float64",
    "jet_eta": "var * float64",
    "jet_phi": "var * float64",
    "jet_radius": "var * float64",
    "jet_mass": "var * float64",
    "jet_area": "var * float64",
    "pair_photon_index": "var * int32",
    "pair_jet_index": "var * int32",
    "pair_delta_phi": "var * float64",
    "pair_xjgamma": "var * float64",
    "pair_recoil_state": "var * int32",
}
for _region in "ABCD":
    EVENT_ARRAY_SCHEMA[f"leader_{_region}_r04_index"] = "int32"
for _region in "CD":
    EVENT_ARRAY_SCHEMA[f"leader_{_region}_r04_complement_index"] = "int32"

TRUTH_PHOTON_SCHEMA = {
    "source_file_index": "int32",
    "event_id_hi": "uint64",
    "event_id_lo": "uint64",
    "truth_photon_id_hi": "uint64",
    "truth_photon_id_lo": "uint64",
    "truth_photon_pt": "float64",
    "truth_photon_eta": "float64",
    "truth_photon_phi": "float64",
    "prompt_class": "int32",
    "source_role": "int32",
    "sample_source_role": "int32",
    "generator_occurrence_embedding_id": "int32",
    "generator_barcode": "int32",
    "truth_isolation": "float64",
    "truth_isolation_r03": "float64",
    "truth_isolation_r04": "float64",
    "truth_isolation_valid": "int32",
    "g4_photon_valid": "int32",
    "hepmc_association_valid": "int32",
    "analysis_signal_r03": "int32",
    "native_track_id": "int32",
    "native_vertex_id": "int32",
    "embedding_id": "int32",
}

TRUTH_JET_SCHEMA = {
    "source_file_index": "int32",
    "event_id_hi": "uint64",
    "event_id_lo": "uint64",
    "truth_jet_id_hi": "uint64",
    "truth_jet_id_lo": "uint64",
    "truth_jet_radius": "float64",
    "truth_jet_pt": "float64",
    "truth_jet_eta": "float64",
    "truth_jet_phi": "float64",
}

LINK_SCHEMA = {
    "source_file_index": "int32",
    "event_id_hi": "uint64",
    "event_id_lo": "uint64",
    "link_id_hi": "uint64",
    "link_id_lo": "uint64",
    "reco_type": "int32",
    "reco_id_hi": "uint64",
    "reco_id_lo": "uint64",
    "reco_index": "int32",
    "truth_type": "int32",
    "truth_id_hi": "uint64",
    "truth_id_lo": "uint64",
    "truth_index": "int32",
    "match_metric": "float64",
    "link_class": "int32",
}

# Authoritative per-view jet response. RJRecoTruthLinkV1 remains a lossless raw
# transport table because its historical final rows may mix reconstruction
# views. Candidate edges retain their source identities; rebuilt final rows use
# fresh deterministic identities above the source high-water mark.
JET_RESPONSE_SCHEMA = {
    **LINK_SCHEMA,
    "jet_input_identity": "string",
    "jet_subtraction_identity": "string",
    # 0: unsupported legacy view evidence, 1: rebuilt final, 2: source candidate.
    "jet_response_state": "int32",
}


EVENT_JOIN_SCHEMA = {key: EVENT_SCHEMA[key] for key in (
    "source_file_index", "source_entry", "event_id_hi", "event_id_lo")}
PAIR_OBJECT_IDS = {key: "uint64" for key in (
    "candidate_id_hi", "candidate_id_lo", "jet_id_hi", "jet_id_lo")}


def validate_io_paths(input_paths, output_path):
    """Reject duplicate physical files and output aliases before any write.

    Distinct copied files still require the campaign's source/range manifest;
    a file-list position is not proof of independent physical exposure.
    """
    paths, seen_paths, seen_files = [], set(), set()
    for supplied in input_paths:
        path = Path(supplied).resolve(strict=True)
        if not path.is_file():
            raise ValueError(f"input is not a regular file: {path}")
        stat = path.stat()
        identity = (stat.st_dev, stat.st_ino)
        if path in seen_paths or identity in seen_files:
            raise ValueError(f"duplicate input file or filesystem alias: {path}")
        paths.append(path)
        seen_paths.add(path)
        seen_files.add(identity)
    if not paths:
        raise ValueError("at least one input ROOT is required")
    output = Path(output_path).resolve()
    if output in seen_paths or (output.exists() and
            (output.stat().st_dev, output.stat().st_ino) in seen_files):
        raise ValueError("output aliases an input file; refusing before output creation")
    seen_uuids = set()
    for path in paths:
        with uproot.open(path) as root:
            uuid = str(root.file.uuid)
            if uuid in seen_uuids:
                raise ValueError(f"duplicate ROOT file UUID (copied input): {path}")
            seen_uuids.add(uuid)
    return paths, output


def source_file_record(root, path, source_index):
    """Small once-per-source locator and preserved upstream provenance.

    This is not a substitute for the accepted source/range campaign manifest.
    In particular, compact SOURCE(0,0) identities are not globally unique.
    """
    metadata = {}
    for key, classname in root[PREFIX.rstrip("/")].classnames(recursive=False, cycle=False).items():
        if classname == "TNamed":
            metadata[key] = root[PREFIX + key].member("fTitle")
        elif classname == "TObjString":
            metadata[key] = str(root[PREFIX + key])
    source_tree = PREFIX + "RJSourceOccurrenceV1"
    occurrences = (ak.to_list(root[source_tree].arrays(library="ak"))
                   if source_tree in root else [])
    return {"source_file_index": source_index, "input_path": str(path.resolve()),
            "root_uuid": str(root.file.uuid), "upstream_metadata": metadata,
            "source_occurrences": occurrences}


def output_schemas(layout="expanded"):
    if layout not in LAYOUTS:
        raise ValueError(f"unknown storage layout: {layout}")
    compact = layout == "normalized_v1"
    join = EVENT_JOIN_SCHEMA if compact else EVENT_SCHEMA
    schemas = {
        "events": {**EVENT_SCHEMA, **MBD_ARRAY_SCHEMA, **centrality_replay.ARRAYS,
                   **TRUTH_VERTEX_ARRAY_SCHEMA, **TRUTH_JET_AVAILABILITY_SCHEMA}
                  if compact else EVENT_SCHEMA,
        "photons": {**join, **PHOTON_FIELDS},
        "jets": {**join, **JET_FIELDS},
        "photonJets": {**join, **(PAIR_OBJECT_IDS if compact else
                                 {**PHOTON_FIELDS, **JET_FIELDS}), **PAIR_FIELDS},
        "truthPhotons": TRUTH_PHOTON_SCHEMA,
        "truthJets": TRUTH_JET_SCHEMA,
        "recoTruthLinks": LINK_SCHEMA,
        "jetResponseLinks": JET_RESPONSE_SCHEMA,
        "jetViewAvailability": JET_VIEW_AVAILABILITY_SCHEMA,
        "triggerScalers": trigger_scaler_interface.SCHEMA,
        "triggerRunInfo": gl1_trigger_interface.RUN_SCHEMA,
    }
    if not compact:
        schemas["eventTree"] = {**EVENT_SCHEMA, **EVENT_ARRAY_SCHEMA,
                                **MBD_ARRAY_SCHEMA, **centrality_replay.ARRAYS,
                                **TRUTH_VERTEX_ARRAY_SCHEMA, **TRUTH_JET_AVAILABILITY_SCHEMA}
    return schemas


def _tree(root: uproot.ReadOnlyDirectory, name: str) -> uproot.behaviors.TTree.TTree:
    key = PREFIX + name
    if key not in root:
        raise ValueError(f"missing required tree: {key}")
    return root[key]


def _arrays(
    tree: uproot.behaviors.TTree.TTree,
    branches: Iterable[str] | None = None,
) -> dict[str, np.ndarray]:
    """Read only the columns consumed by this adapter.

    Some normalized input trees contain large audit-only payloads.  Reading
    every branch multiplied the resident set without changing the output.
    """
    if branches is None:
        return dict(tree.arrays(library="np"))
    available = set(tree.keys())
    selected = [name for name in branches if name in available]
    return dict(tree.arrays(selected, library="np"))


EVENT_INPUT_BRANCHES = tuple(
    name for name in EVENT_SCHEMA if name not in {"source_file_index", "source_entry"}
) + tuple(MBD_ARRAY_SCHEMA) + tuple(centrality_replay.ARRAYS) + tuple(TRUTH_VERTEX_ARRAY_SCHEMA) + tuple(TRUTH_JET_AVAILABILITY_SCHEMA) + RECO_JET_VIEW_AVAILABILITY_INPUT_BRANCHES
WEIGHT_INPUT_BRANCHES = (
    "target_id_hi", "target_id_lo", "component_type", *WEIGHT_FIELDS,
    "application_count",
)
CANDIDATE_INPUT_BRANCHES = (
    "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
    "encounter_ordinal", "cluster_et", "eta", "phi",
    "native_cluster_key", "ppg12_tower_mask_state", "ppg12_source_eligible_state",
    "reference_preselection_state", "active_preselection_state",
    "truth_signal_match_state", "truth_signal_match_barcode",
    "dominant_truth_state", "dominant_truth_evaluator_mode",
    "dominant_truth_track_id", "dominant_truth_pid", "dominant_truth_barcode",
    "dominant_truth_embedding_id", "dominant_truth_energy_contribution",
    "dominant_truth_vertex_id",
    "native_weta_cogx", "native_wphi_cogx", "native_weta33_cogx",
    "native_wphi33_cogx", "native_e11_over_e33", "native_e32_over_e35", "native_et1",
)
MODEL_INPUT_BRANCHES = (
    "candidate_id_hi", "candidate_id_lo", "raw_score",
    "shower_definition_id", "ordered_input_witnesses",
)
ISOLATION_INPUT_BRANCHES = (
    "candidate_id_hi", "candidate_id_lo", "radius", "cone_sum",
    "threshold", "sideband_threshold", "pass_state",
)
SHOWER_VIEW_INPUT_BRANCHES = (
    "candidate_id_hi", "candidate_id_lo", "definition_name",
    "ordered_features", "finite_feature_state",
)
JET_INPUT_BRANCHES = (
    "event_id_hi", "event_id_lo", "jet_id_hi", "jet_id_lo", "radius",
    "raw_pt", "corrected_pt", "eta", "phi", "mass", "area",
    "quality_bitmask", "deterministic_order", "input_identity",
    "subtraction_identity",
)
PAIR_INPUT_BRANCHES = (
    "event_id_hi", "event_id_lo", "pair_id_hi", "pair_id_lo",
    "candidate_id_hi", "candidate_id_lo", "jet_id_hi", "jet_id_lo",
    "delta_phi", "xjgamma", "recoil_state", "photon_rank", "jet_rank",
    "wrong_photon_class", "wrong_recoil_class",
)
TRUTH_PHOTON_INPUT_BRANCHES = (
    "event_id_hi", "event_id_lo", "truth_photon_id_hi", "truth_photon_id_lo",
    "pt", "eta", "phi", "prompt_class", "source_role", "generator_barcode",
    "truth_isolation_witness",
    "truth_isolation_r03", "truth_isolation_r04", "truth_isolation_valid",
    "g4_photon_valid", "hepmc_association_valid", "analysis_signal_r03",
    "native_track_id", "native_vertex_id", "embedding_id",
    "sample_source_role", "generator_occurrence_embedding_id",
)
TRUTH_JET_INPUT_BRANCHES = (
    "event_id_hi", "event_id_lo", "truth_jet_id_hi", "truth_jet_id_lo",
    "radius", "pt", "eta", "phi",
)
LINK_INPUT_BRANCHES = (
    "link_id_hi", "link_id_lo", "reco_type", "reco_id_hi", "reco_id_lo",
    "truth_type", "truth_id_hi", "truth_id_lo", "match_metric", "link_class",
)


def _column(arrays: dict[str, np.ndarray], name: str, length: int, dtype: str, default: Any) -> np.ndarray:
    if name in arrays:
        return np.asarray(arrays[name], dtype=dtype)
    return np.full(length, default, dtype=dtype)


def _key(hi: Any, lo: Any) -> tuple[int, int]:
    return int(hi), int(lo)


def _insert_unique(mapping: dict[Any, Any], key: Any, value: Any, label: str) -> None:
    if key in mapping:
        raise ValueError(f"duplicate {label}: {key}")
    mapping[key] = value


def _event_columns(events: dict[str, np.ndarray], source_index: int) -> dict[str, np.ndarray]:
    n = len(events["event_id_hi"])
    decisions = gl1_trigger_interface.event_columns(events)
    scaled = _column(events, "scaled_trigger_bits", n, "uint64", UINT64_ZERO)
    recorded = _column(events, "scaled_bit30", n, "int32", -1)
    reconstructed = ((scaled >> np.uint64(30)) & np.uint64(1)).astype("int32")
    has_scaled_word = "scaled_trigger_bits" in events
    has_scaled_scalar = "scaled_bit30" in events
    if has_scaled_word and has_scaled_scalar and not np.array_equal(recorded, reconstructed):
        bad = int(np.flatnonzero(recorded != reconstructed)[0])
        raise ValueError(f"scaled-bit-30 mismatch at event row {bad}")
    bit30 = recorded if has_scaled_scalar else (reconstructed if has_scaled_word else np.full(n, -1, dtype="int32"))
    modern = decisions["trigger_capture_version"] == 1
    bit30 = np.where(modern & (decisions["trigger_decision_state"] != 0), -1, bit30).astype("int32")
    return {
        **decisions,
        "source_file_index": np.full(n, source_index, dtype="int32"),
        "source_entry": np.arange(n, dtype="int64"),
        "event_id_hi": _column(events, "event_id_hi", n, "uint64", UINT64_ZERO),
        "event_id_lo": _column(events, "event_id_lo", n, "uint64", UINT64_ZERO),
        "run": _column(events, "run", n, "int32", -1),
        "event_sequence": _column(events, "event_sequence", n, "int64", -1),
        "physical_event_sequence": _column(events, "physical_event_sequence", n, "int64", -1),
        "physical_event_sequence_valid": _column(events, "physical_event_sequence_valid", n, "int32", -1),
        **{name: _column(events, name, n, dtype, math.nan if dtype == "float64" else -1)
           for name, dtype in PPG12_EVENT_WITNESS_SCHEMA.items()},
        "trigger_bits": _column(events, "trigger_bits", n, "uint64", UINT64_ZERO),
        "live_trigger_bits": _column(events, "live_trigger_bits", n, "uint64", UINT64_ZERO),
        "scaled_trigger_bits": scaled,
        "scaled_bit30": bit30,
        "vertex_z": _column(events, "vertex_z", n, "float64", math.nan),
        "centrality": _column(events, "centrality", n, "float64", math.nan),
        "event_weight": _column(events, "event_weight", n, "float64", 1.0),
        "total_calo_energy": _column(events, "total_calo_energy", n, "float64", math.nan),
        "terminal_status": _column(events, "terminal_status", n, "int32", -1),
        "trigger_scaler_snapshot_id": _column(events, "trigger_scaler_snapshot_id", n, "uint64", UINT64_ZERO),
    }


def _truth_inventory_columns(events):
    """Retain the native census, distinguishing absent from certified-empty."""
    n = len(events["event_id_hi"])
    present = set(TRUTH_INVENTORY_SCHEMA) & set(events)
    if present and present != set(TRUTH_INVENTORY_SCHEMA):
        raise ValueError("partial truth-photon inventory")
    result = {name: _column(events, name, n, "int32", -2 if name.endswith("capture_state") else -1)
              for name in TRUTH_INVENTORY_SCHEMA}
    if not present:
        return result
    def get(name): return result["truth_photon_"+name]
    state = get("capture_state")
    if not np.all(np.isin(state, [-2, -1, 0, 1])):
        raise ValueError("invalid truth-photon inventory state")
    captured = state == 1
    counts = [value for name, value in result.items() if not name.endswith("capture_state")]
    if any(np.any(captured & (value < 0)) for value in counts):
        raise ValueError("negative captured truth-photon inventory")
    closes = ((get("native_embedded_primary_pid22_count") == get("serialized_raw_count") + get("rejected_count"))
              & (get("rejected_count") == get("rejected_invalid_kinematics_count")
                 + get("rejected_nonfinite_isolation_count") + get("duplicate_track_count"))
              & (get("analysis_signal_count") <= get("serialized_raw_count")))
    certified = _column(events, "truth_denominator_complete", n, "int32", -1) == 1
    complete = (closes & (get("isolation_input_incomplete_count") == 0)
                & (get("duplicate_track_count") == 0)
                & (get("rejected_nonfinite_isolation_count") == 0))
    # DATA is explicitly not applicable. An incomplete SIM census is retained
    # for diagnosis but must never carry a successful denominator flag.
    if np.any(certified & (state != -1) & (~captured | ~complete)):
        raise ValueError("truth denominator certification contradicts inventory")
    return result


def _mbd_columns(events):
    """Keep channel order and observed charge; sums require 128 valid PMTs."""
    n = len(events["event_id_hi"])
    available = _column(events, "mbd_pmt_available", n, "int32", -1)
    if not np.all(np.isin(available, [-1, 0, 1])):
        raise ValueError("invalid mbd_pmt_available")
    present = set(MBD_ARRAY_SCHEMA) & set(events)
    if present and (present != set(MBD_ARRAY_SCHEMA) or "mbd_pmt_available" not in events):
        raise ValueError("incomplete MBD per-PMT source fields")
    arrays = {key: [] for key in MBD_ARRAY_SCHEMA}
    scalars = {key: np.full(n, np.nan) for key in
               ("mbd_pmt_sum_charge", "mbd_pmt_south_charge", "mbd_pmt_north_charge")}
    scalars["mbd_pmt_available"] = available
    scalars["mbd_pmt_sums_valid"] = np.zeros(n, dtype="int32")
    scalars["mbd_pmt_arm_from_channel_map"] = np.full(n, -1, dtype="int32")
    for i in range(n):
        values = {key: np.asarray(events[key][i], dtype=kind.split()[-1]) if present
                  else np.asarray([], dtype=kind.split()[-1])
                  for key, kind in MBD_ARRAY_SCHEMA.items()}
        lengths = {len(value) for value in values.values()}
        if len(lengths) != 1:
            raise ValueError("MBD per-PMT array lengths differ")
        ids, arms, charge, valid = (values[key] for key in MBD_ARRAY_SCHEMA)
        if available[i] != 1 and len(ids):
            raise ValueError("MBD values require available=1")
        if available[i] == 1 and not present:
            raise ValueError("MBD marked available without per-PMT arrays")
        if len(set(ids.tolist())) != len(ids) or np.any((ids < 0) | (ids >= 128)):
            raise ValueError("duplicate or out-of-range MBD PMT identity")
        if not np.all(np.isin(arms, [-1, 0, 1])) or not np.all(np.isin(valid, [0, 1])):
            raise ValueError("invalid MBD arm or charge validity")
        if np.any((valid == 1) & ~np.isfinite(charge)):
            raise ValueError("valid MBD PMT requires finite calibrated charge")
        if len(ids):
            missing_arm = arms == -1
            # MbdGeom::get_arm(pmtch) is pmtch / 64 in the pinned release.
            # Embedded DSTs can retain calibrated PMTs without an MbdGeom node.
            # Derive only absent arm labels; keep observed labels and charges.
            scalars["mbd_pmt_arm_from_channel_map"][i] = int(np.any(missing_arm))
            arms = np.where(missing_arm, ids // 64, arms).astype("int32")
            values["mbd_pmt_arm"] = arms
        if len(ids) == 128 and np.all(valid == 1) and np.all(arms == ids // 64):
            scalars["mbd_pmt_sums_valid"][i] = 1
            scalars["mbd_pmt_south_charge"][i] = np.sum(charge[arms == 0])
            scalars["mbd_pmt_north_charge"][i] = np.sum(charge[arms == 1])
            scalars["mbd_pmt_sum_charge"][i] = np.sum(charge)
        for key in arrays:
            arrays[key].append(values[key].tolist())
    return scalars, {key: _jagged(values, MBD_ARRAY_SCHEMA[key].split()[-1])
                     for key, values in arrays.items()}


def _truth_vertex_columns(events):
    """Preserve native vertex identities without guessing hard/MB roles."""
    n = len(events["event_id_hi"])
    fields = {"truth_vertex_capture_state", "truth_primary_vertex_id", *TRUTH_VERTEX_ARRAY_SCHEMA}
    present = fields & set(events)
    if present and present != fields:
        raise ValueError("partial truth-vertex capture group")
    scalar = {
        "truth_vertex_capture_state": _column(events, "truth_vertex_capture_state", n, "int32", -2),
        "truth_primary_vertex_id": _column(events, "truth_primary_vertex_id", n, "int32", -1),
    }
    arrays = {key: _jagged([np.asarray(events[key][i], dtype=dtype.split()[-1]).tolist()
                            if present else []
                            for i in range(n)], dtype.split()[-1])
              for key, dtype in TRUTH_VERTEX_ARRAY_SCHEMA.items()}
    for i, state in enumerate(scalar["truth_vertex_capture_state"]):
        if state not in (-2, -1, 0, 1):
            raise ValueError("invalid truth-vertex capture state")
        values = {key: np.asarray(array[i]) for key, array in arrays.items()}
        lengths = {len(v) for v in values.values()}
        if len(lengths) != 1 or (state != 1 and lengths != {0}):
            raise ValueError("malformed truth-vertex capture arrays")
        ids = values["truth_vertex_id"]
        if len(set(ids)) != len(ids) or np.any(ids <= 0):
            raise ValueError("duplicate or nonprimary truth-vertex identity")
        for flag in ("truth_vertex_embedding_valid", "truth_vertex_z_valid"):
            if not np.all(np.isin(values[flag], [0, 1])):
                raise ValueError("invalid truth-vertex validity flag")
        if np.any((values["truth_vertex_z_valid"] == 1) & ~np.isfinite(values["truth_vertex_z_values"])):
            raise ValueError("valid truth vertex has nonfinite z")
    return scalar, arrays


def _truth_jet_availability(events):
    n = len(events["event_id_hi"])
    fields = {"truth_jet_capture_state", *TRUTH_JET_AVAILABILITY_SCHEMA}
    present = fields & set(events)
    if present and present != fields:
        raise ValueError("partial truth-jet availability group")
    state = _column(events, "truth_jet_capture_state", n, "int32", -2)
    arrays = {key: _jagged([np.asarray(events[key][i], dtype="int32").tolist()
                           if present else [] for i in range(n)], "int32")
              for key in TRUTH_JET_AVAILABILITY_SCHEMA}
    for i, value in enumerate(state):
        codes = np.asarray(arrays["truth_jet_radius_code"][i])
        valid = np.asarray(arrays["truth_jet_container_valid"][i])
        if value not in (-2, -1, 0, 1) or len(codes) != len(valid):
            raise ValueError("invalid truth-jet availability")
        if (value != 1 and len(codes)) or len(set(codes)) != len(codes) or np.any(codes <= 0):
            raise ValueError("invalid truth-jet radius availability")
        if not np.all(np.isin(valid, [0, 1])):
            raise ValueError("invalid truth-jet container validity")
    return {"truth_jet_capture_state": state}, arrays


def _reco_jet_view_availability(events):
    """Expand exact per-event producer witnesses; None means unsupported legacy."""
    present = set(RECO_JET_VIEW_AVAILABILITY_INPUT_BRANCHES) & set(events)
    if present and present != set(RECO_JET_VIEW_AVAILABILITY_INPUT_BRANCHES):
        raise ValueError("partial reconstructed-jet view availability group")
    if not present:
        return None, []
    by_event, rows = {}, []
    for i in range(len(events["event_id_hi"])):
        columns = [list(events[name][i]) for name in RECO_JET_VIEW_AVAILABILITY_INPUT_BRANCHES]
        if len({len(values) for values in columns}) != 1:
            raise ValueError("reconstructed-jet view availability lengths differ")
        event = _key(events["event_id_hi"][i], events["event_id_lo"][i])
        views, seen = [], set()
        for radius_code, input_name, subtraction_name, available in zip(*columns):
            input_name = input_name.decode() if isinstance(input_name, bytes) else str(input_name)
            subtraction_name = (subtraction_name.decode() if isinstance(subtraction_name, bytes)
                                else str(subtraction_name))
            key = (input_name, subtraction_name, int(radius_code) / 10.0)
            if (int(radius_code) not in (2, 3, 4) or int(available) != 1 or
                    not input_name or input_name.strip() != input_name or
                    not subtraction_name or subtraction_name.strip() != subtraction_name or
                    key in seen):
                raise ValueError("invalid reconstructed-jet view availability")
            seen.add(key)
            views.append(key)
            rows.append({
                "event_index": i,
                "reco_jet_view_radius_code": int(radius_code),
                "reco_jet_view_input_identity": input_name,
                "reco_jet_view_subtraction_identity": subtraction_name,
                "reco_jet_view_available": 1,
            })
        by_event[event] = views
    return by_event, rows


def _interface_columns(
    events: dict[str, np.ndarray], weights: dict[str, np.ndarray] | None,
    system: str, *, require_complete: bool = False,
) -> dict[str, np.ndarray]:
    """Expose source witnesses without treating finite reset values as validity."""
    n = len(events["event_id_hi"])
    result = {
        name: np.full(n, 0 if dtype == "uint64" else
                      (-1 if dtype == "int32" else math.nan), dtype=dtype)
        for name, dtype in INTERFACE_SCHEMA.items()
    }
    result["collaborator_interface_version"][:] = 2
    result["legacy_event_weight"] = _column(
        events, "legacy_event_weight", n, "float64", math.nan)
    if "legacy_event_weight" not in events:
        result["legacy_event_weight"] = _column(events, "event_weight", n, "float64", math.nan)
    for name in ("source_occurrence_id_hi", "source_occurrence_id_lo",
                 "truth_vertex_z", "truth_mb_vertex_z", "mbd_total_charge",
                 "nodes_ready", "ppg12_weight_active_state",
                 "ppg12_weight_lane_component_code"):
        if name in events:
            result[name] = np.asarray(events[name], dtype=INTERFACE_SCHEMA[name])
    ids_present = all(name in events for name in
                      ("source_occurrence_id_hi", "source_occurrence_id_lo"))
    result["source_occurrence_id_available"] = (
        ((result["source_occurrence_id_hi"] != 0) |
         (result["source_occurrence_id_lo"] != 0)) & ids_present
    ).astype("int32")
    for name in ("truth_vertex_z", "truth_mb_vertex_z", "mbd_total_charge"):
        result[name + "_finite"] = np.isfinite(result[name]).astype("int32")
    result.update(_mbd_columns(events)[0])
    result.update(centrality_replay.capture_columns(events)[0])
    result.update(_truth_vertex_columns(events)[0])
    result.update(_truth_jet_availability(events)[0])
    result.update(_truth_inventory_columns(events))
    for name in ("reco_vertex_valid", "sample_weight_valid",
                 "vertex_weight_valid", "truth_denominator_complete"):
        if name in events:
            values = np.asarray(events[name])
            if not np.all(np.isin(values, [-1, 0, 1])):
                raise ValueError(f"invalid tri-state source witness: {name}")
            result[name] = values.astype("int32")
    reco = _column(events, "vertex_z", n, "float64", math.nan)
    valid_reco = result["reco_vertex_valid"] == 1
    if np.any(valid_reco & ~np.isfinite(reco)):
        raise ValueError("reco_vertex_valid=1 requires finite vertex_z")
    result["reco_vertex_z"][valid_reco] = reco[valid_reco]
    for cut in (30, 60):
        flags = np.full(n, -1, dtype="int32")
        flags[result["reco_vertex_valid"] == 0] = 0
        flags[valid_reco] = (np.abs(reco[valid_reco]) < cut).astype("int32")
        result[f"reco_vertex_abs_lt_{cut}"] = flags

    event_lookup: dict[tuple[int, int], int] = {}
    for i, (hi, lo) in enumerate(zip(events["event_id_hi"], events["event_id_lo"])):
        _insert_unique(event_lookup, _key(hi, lo), i, "event identity")
    result["weight_components_recorded"][:] = 0
    if weights is not None:
        missing = set(WEIGHT_INPUT_BRANCHES) - set(weights)
        if missing:
            raise ValueError(f"incomplete RJWeightComponentV1 fields: {sorted(missing)}")
        seen: set[tuple[int, int]] = set()
        for j, component in enumerate(weights["component_type"]):
            if isinstance(component, bytes):
                component = component.decode("utf-8")
            if component != "event":
                continue
            key = _key(weights["target_id_hi"][j], weights["target_id_lo"][j])
            if key in seen:
                raise ValueError(f"duplicate event weight target: {key}")
            seen.add(key)
            if key not in event_lookup:
                raise ValueError(f"orphan event weight target: {key}")
            i = event_lookup[key]
            for name in WEIGHT_FIELDS:
                result["weight_final" if name == "final_weight" else name][i] = weights[name][j]
            result["weight_application_count"][i] = weights["application_count"][j]
            result["weight_components_recorded"][i] = 1
        if len(seen) != n:
            raise ValueError("RJWeightComponentV1 lacks exactly one event component per event")

    recorded = result["weight_components_recorded"] == 1
    result["sample_weight"] = result["slice_weight"].copy()
    if system == "pp":
        aliases = (("slice_weight", "cross_section_weight"),
                   ("period_weight", "exposure_weight"))
        factor_names = ("sample_weight", "si_di_weight", "period_weight")
    elif system == "auau":
        aliases = (("slice_weight", "cross_section_weight"),)
        factor_names = ("sample_weight", "exposure_weight")
        for name in ("slice_weight", "cross_section_weight", "si_di_weight", "period_weight"):
            finite = recorded & np.isfinite(result[name])
            if np.any(finite & (result[name] != 1.0)):
                raise ValueError(f"unsupported AuAu producer factor: {name}")
    else:
        raise ValueError("system must be pp or auau")
    for left, right in aliases:
        if np.any(recorded & ~np.isclose(result[left], result[right],
                                        rtol=1e-12, atol=0, equal_nan=True)):
            raise ValueError(f"producer weight alias mismatch: {left}/{right}")
    with np.errstate(invalid="ignore", over="ignore"):
        without_vertex = np.prod([result[name] for name in factor_names], axis=0)
        reconstructed = without_vertex * result["vertex_weight"]
    result["event_weight_without_vertex"] = without_vertex
    event_weight = _column(events, "event_weight", n, "float64", math.nan)
    finite = recorded & np.isfinite(reconstructed) & np.isfinite(event_weight)
    finite &= np.isfinite(result["weight_final"])
    matches = (np.isclose(reconstructed, event_weight, rtol=1e-10, atol=1e-14) &
               np.isclose(result["weight_final"], event_weight, rtol=1e-10, atol=1e-14))
    if np.any(finite & ~matches):
        raise ValueError("event weight multiplication identity mismatch")
    result["weight_identity_valid"][recorded] = 0
    result["weight_identity_valid"][finite & matches] = 1
    for witness, factor in (("sample_weight_valid", "sample_weight"),
                            ("vertex_weight_valid", "vertex_weight")):
        if np.any((result[witness] == 1) & ~np.isfinite(result[factor])):
            raise ValueError(f"{witness}=1 requires finite recorded {factor}")
    if require_complete:
        for witness in ("sample_weight_valid", "vertex_weight_valid",
                        "truth_denominator_complete", "weight_components_recorded",
                        "weight_identity_valid"):
            if np.any(result[witness] != 1):
                raise ValueError(f"complete interface unavailable: {witness} is not certified by producer")
    return result


def _row(columns: dict[str, np.ndarray], index: int) -> dict[str, Any]:
    return {name: values[index] for name, values in columns.items()}


def _nan_photon() -> dict[str, Any]:
    row: dict[str, Any] = {}
    for name, dtype in PHOTON_FIELDS.items():
        if dtype == "uint64":
            row[name] = UINT64_ZERO
        elif dtype.startswith("int"):
            row[name] = -1
        else:
            row[name] = math.nan
    return row


def _nan_jet() -> dict[str, Any]:
    row: dict[str, Any] = {}
    for name, dtype in JET_FIELDS.items():
        if dtype == "string":
            row[name] = ""
        elif dtype == "uint64":
            row[name] = UINT64_ZERO
        elif dtype.startswith("int"):
            row[name] = -1
        else:
            row[name] = math.nan
    return row


def _to_columns(rows: list[dict[str, Any]], schema: dict[str, str]) -> dict[str, np.ndarray]:
    return {
        name: np.asarray([row[name] for row in rows], dtype=object if dtype == "string" else dtype)
        for name, dtype in schema.items()
    }


def _indices_by_event(arrays: dict[str, np.ndarray]) -> dict[tuple[int, int], list[int]]:
    grouped: dict[tuple[int, int], list[int]] = {}
    for index in range(len(arrays.get("event_id_hi", []))):
        grouped.setdefault(_key(arrays["event_id_hi"][index], arrays["event_id_lo"][index]), []).append(index)
    return grouped


def _filter_by_identity(
    arrays: dict[str, np.ndarray],
    id_hi_name: str,
    id_lo_name: str,
    retained: set[tuple[int, int]],
) -> dict[str, np.ndarray]:
    mask = np.asarray(
        [
            _key(high, low) in retained
            for high, low in zip(arrays.get(id_hi_name, []), arrays.get(id_lo_name, []))
        ],
        dtype="bool",
    )
    if bool(np.all(mask)):
        return arrays
    return {name: values[mask] for name, values in arrays.items()}


def _contiguous_event_slices(
    arrays: dict[str, np.ndarray],
    label: str,
) -> dict[tuple[int, int], tuple[int, int]]:
    """Map each event to one contiguous source slice without per-row Python objects."""
    high = np.asarray(arrays.get("event_id_hi", []), dtype="uint64")
    low = np.asarray(arrays.get("event_id_lo", []), dtype="uint64")
    if len(high) != len(low):
        raise ValueError(f"{label} event identity columns differ in length")
    if len(high) == 0:
        return {}
    starts_mask = np.empty(len(high), dtype="bool")
    starts_mask[0] = True
    starts_mask[1:] = (high[1:] != high[:-1]) | (low[1:] != low[:-1])
    starts = np.flatnonzero(starts_mask)
    stops = np.concatenate((starts[1:], np.asarray([len(high)], dtype=starts.dtype)))
    slices: dict[tuple[int, int], tuple[int, int]] = {}
    for start, stop in zip(starts, stops):
        event_key = _key(high[start], low[start])
        _insert_unique(slices, event_key, (int(start), int(stop)), f"{label} event segment")
    return slices


def _local_indices_by_contiguous_event(
    arrays: dict[str, np.ndarray],
    label: str,
) -> np.ndarray:
    """Return each row's zero-based index within its contiguous event segment."""
    high = np.asarray(arrays.get("event_id_hi", []), dtype="uint64")
    low = np.asarray(arrays.get("event_id_lo", []), dtype="uint64")
    if len(high) != len(low):
        raise ValueError(f"{label} event identity columns differ in length")
    if len(high) == 0:
        return np.asarray([], dtype="int32")
    starts_mask = np.empty(len(high), dtype="bool")
    starts_mask[0] = True
    starts_mask[1:] = (high[1:] != high[:-1]) | (low[1:] != low[:-1])
    starts = np.flatnonzero(starts_mask)
    stops = np.concatenate((starts[1:], np.asarray([len(high)], dtype=starts.dtype)))
    lengths = stops - starts
    local = np.arange(len(high), dtype="int64") - np.repeat(starts, lengths)
    if local.size and int(np.max(local)) > np.iinfo(np.int32).max:
        raise ValueError(f"{label} has more than int32 rows in one event")
    return local.astype("int32", copy=False)


def _resolve_identity_source_indices(
    source_arrays: dict[str, np.ndarray],
    source_hi_name: str,
    source_lo_name: str,
    query_arrays: dict[str, np.ndarray],
    query_hi_name: str,
    query_lo_name: str,
    label: str,
) -> np.ndarray:
    """Vectorized exact two-word identity join with bounded array memory."""
    identity_dtype = np.dtype([("hi", "<u8"), ("lo", "<u8")])
    source_high = np.asarray(source_arrays[source_hi_name], dtype="uint64")
    source_low = np.asarray(source_arrays[source_lo_name], dtype="uint64")
    query_high = np.asarray(query_arrays[query_hi_name], dtype="uint64")
    query_low = np.asarray(query_arrays[query_lo_name], dtype="uint64")
    if len(source_high) != len(source_low) or len(query_high) != len(query_low):
        raise ValueError(f"{label} identity columns differ in length")

    source_ids = np.empty(len(source_high), dtype=identity_dtype)
    source_ids["hi"] = source_high
    source_ids["lo"] = source_low
    order = np.argsort(source_ids, order=("hi", "lo"), kind="stable")
    sorted_ids = source_ids[order]
    if len(sorted_ids) > 1 and bool(np.any(sorted_ids[1:] == sorted_ids[:-1])):
        raise ValueError(f"duplicate {label} source identity")

    query_ids = np.empty(len(query_high), dtype=identity_dtype)
    query_ids["hi"] = query_high
    query_ids["lo"] = query_low
    positions = np.searchsorted(sorted_ids, query_ids)
    missing = positions == len(sorted_ids)
    if bool(np.any(~missing)):
        present_positions = positions[~missing]
        missing[~missing] = sorted_ids[present_positions] != query_ids[~missing]
    if bool(np.any(missing)):
        first = int(np.flatnonzero(missing)[0])
        raise ValueError(
            f"{label} query identity has no source join at row {first}: "
            f"{_key(query_high[first], query_low[first])}"
        )
    return np.asarray(order[positions], dtype="int64")


def _local_identity_lookup(
    arrays: dict[str, np.ndarray],
    id_hi_name: str,
    id_lo_name: str,
    label: str,
) -> dict[tuple[int, int], tuple[tuple[int, int], int]]:
    lookup: dict[tuple[int, int], tuple[tuple[int, int], int]] = {}
    for event_key, indices in _indices_by_event(arrays).items():
        for local_index, source_index in enumerate(indices):
            _insert_unique(
                lookup,
                _key(arrays[id_hi_name][source_index], arrays[id_lo_name][source_index]),
                (event_key, local_index),
                label,
            )
    return lookup


def _leader_index(rows: list[dict[str, Any]], region: str, *, complement: bool = False) -> int:
    eligible: list[int] = []
    for index, row in enumerate(rows):
        score = float(row["bdt_score"])
        tight_threshold = float(row["bdt_tight_threshold"])
        nontight_low = float(row["bdt_nontight_low_threshold"])
        nontight_high = float(row["bdt_nontight_high_threshold"])
        isolation = float(row["iso_r04"])
        isolated_threshold = float(row["iso_r04_threshold"])
        nonisolated_threshold = float(row["iso_r04_nonisolated_threshold"])
        if not all(
            math.isfinite(value)
            for value in (
                score, tight_threshold, nontight_low, nontight_high,
                isolation, isolated_threshold, nonisolated_threshold,
            )
        ):
            continue
        tight = score > tight_threshold
        nontight = (not tight) if complement else nontight_low < score < nontight_high
        isolated = isolation < isolated_threshold
        nonisolated = isolation > nonisolated_threshold
        accepted = {
            "A": tight and isolated,
            "B": tight and nonisolated,
            "C": nontight and isolated,
            "D": nontight and nonisolated,
        }[region]
        if accepted:
            eligible.append(index)
    if not eligible:
        return -1
    return min(
        eligible,
        key=lambda index: (
            -float(rows[index]["photon_et"]),
            int(rows[index]["photon_encounter_ordinal"]),
            index,
        ),
    )


def _jagged(rows: list[list[Any]], dtype: str) -> ak.Array:
    counts = np.asarray([len(values) for values in rows], dtype="int64")
    nonempty = [np.asarray(values, dtype=dtype) for values in rows if values]
    flat = np.concatenate(nonempty) if nonempty else np.asarray([], dtype=dtype)
    return ak.unflatten(flat, counts)


def _event_tree_context(
    candidate_arrays: dict[str, np.ndarray],
    jet_arrays: dict[str, np.ndarray],
    pair_arrays: dict[str, np.ndarray],
) -> tuple[
    dict[tuple[int, int], tuple[int, int]],
    dict[tuple[int, int], tuple[int, int]],
    dict[tuple[int, int], tuple[int, int]],
]:
    """Build one compact slice per event; output rows remain batched."""
    return (
        _contiguous_event_slices(candidate_arrays, "photon candidate"),
        _contiguous_event_slices(jet_arrays, "jet"),
        _contiguous_event_slices(pair_arrays, "photon-jet pair"),
    )


def _pair_join_indices(
    candidate_arrays: dict[str, np.ndarray],
    jet_arrays: dict[str, np.ndarray],
    pair_arrays: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Resolve pair joins and local event indices without source-wide Python maps."""
    photon_source = _resolve_identity_source_indices(
        candidate_arrays,
        "candidate_id_hi",
        "candidate_id_lo",
        pair_arrays,
        "candidate_id_hi",
        "candidate_id_lo",
        "photon candidate",
    )
    photon_local_by_source = _local_indices_by_contiguous_event(candidate_arrays, "photon candidate")
    photon_local = photon_local_by_source[photon_source]
    del photon_local_by_source

    jet_source = _resolve_identity_source_indices(
        jet_arrays,
        "jet_id_hi",
        "jet_id_lo",
        pair_arrays,
        "jet_id_hi",
        "jet_id_lo",
        "jet",
    )
    jet_local_by_source = _local_indices_by_contiguous_event(jet_arrays, "jet")
    jet_local = jet_local_by_source[jet_source]
    del jet_local_by_source

    pair_event_hi = np.asarray(pair_arrays["event_id_hi"], dtype="uint64")
    pair_event_lo = np.asarray(pair_arrays["event_id_lo"], dtype="uint64")
    photon_event_ok = (
        np.asarray(candidate_arrays["event_id_hi"], dtype="uint64")[photon_source] == pair_event_hi
    ) & (
        np.asarray(candidate_arrays["event_id_lo"], dtype="uint64")[photon_source] == pair_event_lo
    )
    jet_event_ok = (
        np.asarray(jet_arrays["event_id_hi"], dtype="uint64")[jet_source] == pair_event_hi
    ) & (
        np.asarray(jet_arrays["event_id_lo"], dtype="uint64")[jet_source] == pair_event_lo
    )
    event_ok = photon_event_ok & jet_event_ok
    if not bool(np.all(event_ok)):
        first = int(np.flatnonzero(~event_ok)[0])
        raise ValueError(f"photon-jet pair crosses event boundaries at row {first}")
    return photon_source, jet_source, photon_local, jet_local


def _event_tree_columns(
    event_columns: dict[str, np.ndarray],
    candidate_arrays: dict[str, np.ndarray],
    photon_rows: list[dict[str, Any]],
    jet_arrays: dict[str, np.ndarray],
    jet_rows: dict[str, np.ndarray] | list[dict[str, Any]],
    pair_arrays: dict[str, np.ndarray],
    *,
    event_start: int = 0,
    event_stop: int | None = None,
    context: tuple[
        dict[tuple[int, int], tuple[int, int]],
        dict[tuple[int, int], tuple[int, int]],
        dict[tuple[int, int], tuple[int, int]],
    ] | None = None,
    pair_photon_local: np.ndarray | None = None,
    pair_jet_local: np.ndarray | None = None,
) -> dict[str, Any]:
    if context is None:
        context = _event_tree_context(candidate_arrays, jet_arrays, pair_arrays)
    photon_groups, jet_groups, pair_groups = context
    if pair_photon_local is None or pair_jet_local is None:
        _, _, pair_photon_local, pair_jet_local = _pair_join_indices(
            candidate_arrays, jet_arrays, pair_arrays
        )
    total_events = len(event_columns["event_id_hi"])
    if event_stop is None:
        event_stop = total_events
    if not 0 <= event_start <= event_stop <= total_events:
        raise ValueError("event batch bounds differ")

    branch_rows: dict[str, list[list[Any]]] = {
        name: []
        for name, branch_type in EVENT_ARRAY_SCHEMA.items()
        if branch_type.startswith("var *")
    }
    scalar_rows: dict[str, list[int]] = {
        name: []
        for name, branch_type in EVENT_ARRAY_SCHEMA.items()
        if not branch_type.startswith("var *")
    }

    for event_index in range(event_start, event_stop):
        event_key = _key(event_columns["event_id_hi"][event_index], event_columns["event_id_lo"][event_index])
        p_slice = photon_groups.get(event_key, (0, 0))
        j_slice = jet_groups.get(event_key, (0, 0))
        q_slice = pair_groups.get(event_key, (0, 0))
        p_indices = range(*p_slice)
        j_indices = range(*j_slice)
        q_indices = range(*q_slice)
        event_photons = [photon_rows[index] for index in p_indices]

        scalar_rows["nphotons"].append(p_slice[1] - p_slice[0])
        scalar_rows["njets"].append(j_slice[1] - j_slice[0])
        scalar_rows["npairs"].append(q_slice[1] - q_slice[0])
        for region in "ABCD":
            scalar_rows[f"leader_{region}_r04_index"].append(
                _leader_index(event_photons, region)
            )
        for region in "CD":
            scalar_rows[f"leader_{region}_r04_complement_index"].append(
                _leader_index(event_photons, region, complement=True)
            )

        photon_value_map = {
            "photon_candidate_id_hi": "candidate_id_hi",
            "photon_candidate_id_lo": "candidate_id_lo",
            "photon_encounter_ordinal": "photon_encounter_ordinal",
            "photon_et": "photon_et",
            "photon_eta": "photon_eta",
            "photon_phi": "photon_phi",
            "photon_bdt_score": "bdt_score",
            "photon_bdt_tight_threshold": "bdt_tight_threshold",
            "photon_bdt_nontight_low_threshold": "bdt_nontight_low_threshold",
            "photon_bdt_nontight_high_threshold": "bdt_nontight_high_threshold",
            "photon_bdt_is_tight": "bdt_is_tight",
            "photon_bdt_is_nontight": "bdt_is_nontight",
            "photon_bdt_is_not_tight": "bdt_is_not_tight",
            "photon_iso_r03": "iso_r03",
            "photon_iso_r03_threshold": "iso_r03_threshold",
            "photon_iso_r03_nonisolated_threshold": "iso_r03_nonisolated_threshold",
            "photon_iso_r04": "iso_r04",
            "photon_iso_r04_threshold": "iso_r04_threshold",
            "photon_iso_r04_nonisolated_threshold": "iso_r04_nonisolated_threshold",
        }
        for output_name, row_name in photon_value_map.items():
            branch_rows[output_name].append([photon_rows[index][row_name] for index in p_indices])

        jet_value_map = {
            "jet_id_hi": "jet_id_hi",
            "jet_id_lo": "jet_id_lo",
            "jet_pt": "jet_pt",
            "jet_raw_pt": "jet_raw_pt",
            "jet_eta": "jet_eta",
            "jet_phi": "jet_phi",
            "jet_radius": "jet_radius",
            "jet_mass": "jet_mass",
            "jet_area": "jet_area",
        }
        for output_name, row_name in jet_value_map.items():
            if isinstance(jet_rows, dict):
                branch_rows[output_name].append(
                    np.asarray(jet_rows[row_name][j_slice[0]:j_slice[1]]).tolist()
                )
            else:
                branch_rows[output_name].append([jet_rows[index][row_name] for index in j_indices])

        branch_rows["pair_photon_index"].append(
            np.asarray(pair_photon_local[q_slice[0]:q_slice[1]], dtype="int32").tolist()
        )
        branch_rows["pair_jet_index"].append(
            np.asarray(pair_jet_local[q_slice[0]:q_slice[1]], dtype="int32").tolist()
        )
        branch_rows["pair_delta_phi"].append([float(pair_arrays["delta_phi"][index]) for index in q_indices])
        branch_rows["pair_xjgamma"].append([float(pair_arrays["xjgamma"][index]) for index in q_indices])
        branch_rows["pair_recoil_state"].append([int(pair_arrays["recoil_state"][index]) for index in q_indices])

    columns: dict[str, Any] = {
        name: values[event_start:event_stop] for name, values in event_columns.items()
    }
    for name, values in scalar_rows.items():
        columns[name] = np.asarray(values, dtype=EVENT_ARRAY_SCHEMA[name])
    for name, values in branch_rows.items():
        dtype = EVENT_ARRAY_SCHEMA[name].split("*", 1)[1].strip()
        columns[name] = _jagged(values, dtype)
    return columns


def _make_rbdt_evaluator(model_path: Path | None, model_input_count: int | None):
    if model_path is None:
        if model_input_count is not None:
            raise ValueError("model_input_count was supplied without a model")
        return None
    if model_input_count not in (11, 14):
        raise ValueError("a model requires an exact --model-input-count of 11 (p+p) or 14 (Au+Au)")
    model_path = Path(model_path).resolve()
    if not model_path.is_file():
        raise FileNotFoundError(model_path)
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError("PyROOT is required when --model is supplied") from exc
    model = ROOT.TMVA.Experimental.RBDT("myBDT", str(model_path))

    def evaluate(values: np.ndarray) -> float:
        if len(values) != model_input_count:
            raise ValueError(f"BDT input dimension mismatch: {len(values)} != {model_input_count}")
        inputs = ROOT.std.vector("float")()
        for value in values:
            inputs.push_back(float(value))
        scores = model.Compute(inputs)
        if len(scores) != 1:
            raise ValueError(f"BDT returned {len(scores)} scores instead of one")
        return float(scores[0])

    return evaluate


def _selection_thresholds(system: str, photon_et: float, centrality: float) -> tuple[float, float, float]:
    """Return the established tight and bounded non-tight score surfaces."""

    if system == "pp":
        return (
            0.815625 - 0.0015625 * photon_et,
            0.7333333333333333 - 0.01333333333333333 * photon_et,
            0.684375 + 0.0015625 * photon_et,
        )
    if system == "auau":
        if not math.isfinite(centrality):
            raise ValueError("Au+Au photon has no finite event centrality")
        return (
            0.6927385742113304 + 0.0010476808935335573 * centrality,
            0.4634577166843198 + 0.0015466253667456042 * centrality,
            0.6071724560222583 + 0.0013019695025559648 * centrality,
        )
    raise ValueError(f"unsupported collision system: {system!r}")


def _photon_rows(
    candidates: dict[str, np.ndarray],
    models: dict[str, np.ndarray],
    isolations: dict[str, np.ndarray],
    shower_views: dict[str, np.ndarray],
    centrality_by_event: dict[tuple[int, int], float],
    system: str,
    evaluate_bdt=None,
    model_input_count: int | None = None,
    *, preserve_capture: bool = False,
) -> tuple[list[dict[str, Any]], dict[tuple[int, int], dict[str, Any]]]:
    model_by_candidate: dict[tuple[int, int], int] = {}
    for i in range(len(models.get("candidate_id_hi", []))):
        candidate_key = _key(models["candidate_id_hi"][i], models["candidate_id_lo"][i])
        shower = str(models.get("shower_definition_id", np.asarray([""] * len(models["candidate_id_hi"]), dtype=object))[i])
        if shower.upper() == "H70":
            _insert_unique(model_by_candidate, candidate_key, i, "H70 model evaluation")

    h70_by_candidate: dict[tuple[int, int], tuple[np.ndarray, int]] = {}
    shower_view_count = len(shower_views.get("candidate_id_hi", []))
    finite_feature_states = shower_views.get(
        "finite_feature_state", np.full(shower_view_count, -1, dtype="int32")
    )
    for i in range(len(shower_views.get("candidate_id_hi", []))):
        definition = str(shower_views["definition_name"][i]).upper()
        if definition != "H70":
            continue
        candidate_key = _key(shower_views["candidate_id_hi"][i], shower_views["candidate_id_lo"][i])
        if candidate_key in h70_by_candidate:
            raise ValueError(f"duplicate H70 shower-feature view for candidate {candidate_key}")
        h70_by_candidate[candidate_key] = (
            np.asarray(shower_views["ordered_features"][i], dtype="float64"),
            int(finite_feature_states[i]),
        )

    iso_by_candidate: dict[tuple[int, int], dict[str, int]] = {}
    for i in range(len(isolations.get("candidate_id_hi", []))):
        candidate_key = _key(isolations["candidate_id_hi"][i], isolations["candidate_id_lo"][i])
        radius = float(isolations["radius"][i])
        if abs(radius - 0.3) < 1e-6:
            _insert_unique(iso_by_candidate.setdefault(candidate_key, {}), "r03", i, f"R=0.3 isolation witness for {candidate_key}")
        elif abs(radius - 0.4) < 1e-6:
            _insert_unique(iso_by_candidate.setdefault(candidate_key, {}), "r04", i, f"R=0.4 isolation witness for {candidate_key}")

    rows: list[dict[str, Any]] = []
    lookup: dict[tuple[int, int], dict[str, Any]] = {}
    n = len(candidates["candidate_id_hi"])
    for i in range(n):
        candidate_key = _key(candidates["candidate_id_hi"][i], candidates["candidate_id_lo"][i])
        event_key = _key(candidates["event_id_hi"][i], candidates["event_id_lo"][i])
        if event_key not in centrality_by_event:
            raise ValueError(f"photon candidate {candidate_key} has no event centrality join")
        row = _nan_photon()
        row["bdt_evaluation_state"] = 0
        row.update(
            candidate_id_hi=candidates["candidate_id_hi"][i],
            candidate_id_lo=candidates["candidate_id_lo"][i],
            photon_encounter_ordinal=int(candidates.get("encounter_ordinal", np.full(n, -1))[i]),
            photon_et=float(candidates["cluster_et"][i]),
            photon_eta=float(candidates["eta"][i]),
            photon_phi=float(candidates["phi"][i]),
            truth_matched=int(candidates.get("truth_signal_match_state", np.full(n, -1))[i]),
            truth_barcode=int(candidates.get("truth_signal_match_barcode", np.full(n, -1))[i]),
            # -2 means legacy input did not record this witness, not background.
            dominant_truth_state=int(candidates.get("dominant_truth_state", np.full(n, -2))[i]),
            dominant_truth_evaluator_mode=int(candidates.get("dominant_truth_evaluator_mode", np.full(n, 0))[i]),
            dominant_truth_track_id=int(candidates.get("dominant_truth_track_id", np.full(n, -1))[i]),
            dominant_truth_pid=int(candidates.get("dominant_truth_pid", np.full(n, 0))[i]),
            dominant_truth_barcode=int(candidates.get("dominant_truth_barcode", np.full(n, -1))[i]),
            dominant_truth_embedding_id=int(candidates.get("dominant_truth_embedding_id", np.full(n, 0))[i]),
            dominant_truth_energy_contribution=float(candidates.get("dominant_truth_energy_contribution", np.full(n, math.nan))[i]),
            dominant_truth_vertex_id=int(candidates.get("dominant_truth_vertex_id", np.full(n, -1))[i]),
            native_weta_cogx=float(candidates.get("native_weta_cogx", np.full(n, math.nan))[i]),
            native_wphi_cogx=float(candidates.get("native_wphi_cogx", np.full(n, math.nan))[i]),
            native_weta33_cogx=float(candidates.get("native_weta33_cogx", np.full(n, math.nan))[i]),
            native_wphi33_cogx=float(candidates.get("native_wphi33_cogx", np.full(n, math.nan))[i]),
            native_e11_over_e33=float(candidates.get("native_e11_over_e33", np.full(n, math.nan))[i]),
            native_e32_over_e35=float(candidates.get("native_e32_over_e35", np.full(n, math.nan))[i]),
            native_et1=float(candidates.get("native_et1", np.full(n, math.nan))[i]),
        )
        for name in ("native_cluster_key", "ppg12_tower_mask_state", "ppg12_source_eligible_state",
                     "reference_preselection_state", "active_preselection_state"):
            row[name] = int(candidates[name][i]) if name in candidates else -1
        model_index = model_by_candidate.get(candidate_key)
        h70_view = h70_by_candidate.get(candidate_key)
        witnesses = None if h70_view is None else h70_view[0]
        finite_feature_state = -1 if h70_view is None else h70_view[1]
        if evaluate_bdt is None and witnesses is None and model_index is not None:
            model_definition = str(models.get("shower_definition_id", np.asarray([""] * len(models["candidate_id_hi"]), dtype=object))[model_index]).upper()
            if model_definition == "H70":
                witnesses = np.asarray(models.get("ordered_input_witnesses", np.asarray([], dtype=object))[model_index], dtype="float64")
        if model_index is not None:
            row["bdt_score"] = float(models["raw_score"][model_index])
            if math.isfinite(row["bdt_score"]):
                row["bdt_evaluation_state"] = 2
        centrality = centrality_by_event[event_key]
        model_domain = (MODEL_ET_MIN_GEV <= float(row["photon_et"]) < MODEL_ET_MAX_GEV
                        and (system == "pp" or
                             (math.isfinite(centrality) and 0 <= centrality < 80)))
        # The collaborator branch has one meaning everywhere: the included
        # model evaluated on the exact H70 shower view.  A retained source
        # candidate outside the accepted model window may lack a valid H70
        # view; omit that candidate instead of imputing, falling back to a
        # producer witness, or inventing a finite sentinel score.  The
        # accepted 15--35 GeV population remains fail-closed.
        expected_dim = model_input_count or (11 if system == "pp" else 14)
        exact_inputs_valid = (witnesses is not None and finite_feature_state == 1
                              and len(witnesses) == expected_dim
                              and bool(np.all(np.isfinite(witnesses))))
        if preserve_capture and (not model_domain or not exact_inputs_valid):
            # This is a raw transport row, not a failed photon/background
            # selection. Preserve its features and associations without scoring.
            row["bdt_score"] = math.nan
            row["bdt_evaluation_state"] = 0
        elif evaluate_bdt is not None:
            if not exact_inputs_valid:
                if MODEL_ET_MIN_GEV <= float(row["photon_et"]) < MODEL_ET_MAX_GEV:
                    raise ValueError(
                        f"candidate {candidate_key} lacks valid exact H70 inputs "
                        "inside the accepted 15--35 GeV model window"
                    )
                continue
            row["bdt_score"] = evaluate_bdt(witnesses)
            row["bdt_evaluation_state"] = 1
            if not math.isfinite(row["bdt_score"]):
                raise ValueError(f"candidate {candidate_key} produced a non-finite BDT score")
        if witnesses is not None:
            row["bdt_input_count"] = min(len(witnesses), 14)
            for j, value in enumerate(witnesses[:14]):
                row[f"bdt_input_{j:02d}"] = float(value)
        tight_threshold, nontight_low, nontight_high = (
            (math.nan, math.nan, math.nan) if preserve_capture and not model_domain
            else _selection_thresholds(system, float(row["photon_et"]), centrality))
        if (not preserve_capture or model_domain) and not nontight_low < nontight_high:
            raise ValueError(f"photon candidate {candidate_key} has unordered non-tight thresholds")
        row.update(
            bdt_tight_threshold=tight_threshold,
            bdt_nontight_low_threshold=nontight_low,
            bdt_nontight_high_threshold=nontight_high,
        )
        if math.isfinite(float(row["bdt_score"])) and math.isfinite(tight_threshold):
            row["bdt_is_tight"] = int(float(row["bdt_score"]) > tight_threshold)
            row["bdt_is_nontight"] = int(nontight_low < float(row["bdt_score"]) < nontight_high)
            row["bdt_is_not_tight"] = int(not bool(row["bdt_is_tight"]))
        for radius_name in ("r03", "r04"):
            iso_index = iso_by_candidate.get(candidate_key, {}).get(radius_name)
            if iso_index is None:
                continue
            row.update(
                {
                    f"iso_{radius_name}": float(isolations["cone_sum"][iso_index]),
                    f"iso_{radius_name}_threshold": float(isolations["threshold"][iso_index]),
                    f"iso_{radius_name}_nonisolated_threshold": float(isolations["sideband_threshold"][iso_index]),
                    f"iso_{radius_name}_pass": int(isolations["pass_state"][iso_index]),
                }
            )
        rows.append(row)
        _insert_unique(lookup, candidate_key, row, "photon candidate identity")
    return rows, lookup


def _jet_rows(jets: dict[str, np.ndarray]) -> tuple[list[dict[str, Any]], dict[tuple[int, int], dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    lookup: dict[tuple[int, int], dict[str, Any]] = {}
    n = len(jets["jet_id_hi"])
    for i in range(n):
        row = _nan_jet()
        row.update(
            jet_id_hi=jets["jet_id_hi"][i],
            jet_id_lo=jets["jet_id_lo"][i],
            jet_radius=float(jets["radius"][i]),
            jet_raw_pt=float(jets["raw_pt"][i]),
            jet_pt=float(jets["corrected_pt"][i]),
            jet_eta=float(jets["eta"][i]),
            jet_phi=float(jets["phi"][i]),
            jet_mass=float(jets["mass"][i]),
            jet_area=float(jets["area"][i]),
            jet_quality_bitmask=jets["quality_bitmask"][i],
            jet_order=int(jets["deterministic_order"][i]),
        )
        rows.append(row)
        _insert_unique(lookup, _key(jets["jet_id_hi"][i], jets["jet_id_lo"][i]), row, "jet identity")
    return rows, lookup


def _jet_view_columns(jets: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return exact producer view labels, or an explicit unsupported state."""
    n = len(jets["jet_id_hi"])
    has_input = "input_identity" in jets
    has_subtraction = "subtraction_identity" in jets
    if has_input != has_subtraction:
        raise ValueError("RJJetV1 has partial reconstruction-view provenance")
    if not has_input:
        return (
            np.full(n, "", dtype=object),
            np.full(n, "", dtype=object),
            np.zeros(n, dtype="int32"),
        )
    inputs = np.asarray(jets["input_identity"], dtype=object)
    subtractions = np.asarray(jets["subtraction_identity"], dtype=object)
    for i, (input_name, subtraction_name) in enumerate(zip(inputs, subtractions)):
        if not isinstance(input_name, str) or not isinstance(subtraction_name, str):
            raise ValueError(f"RJJetV1 reconstruction-view provenance is not text at row {i}")
        if not input_name or not subtraction_name or input_name.strip() != input_name or subtraction_name.strip() != subtraction_name:
            raise ValueError(f"RJJetV1 has invalid reconstruction-view provenance at row {i}")
    return inputs, subtractions, np.ones(n, dtype="int32")


def _jet_columns(jets: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Project jet fields as compact typed arrays instead of per-jet dictionaries."""
    n = len(jets["jet_id_hi"])
    inputs, subtractions, provenance = _jet_view_columns(jets)
    return {
        "jet_id_hi": _column(jets, "jet_id_hi", n, "uint64", UINT64_ZERO),
        "jet_id_lo": _column(jets, "jet_id_lo", n, "uint64", UINT64_ZERO),
        "jet_radius": _column(jets, "radius", n, "float64", math.nan),
        "jet_raw_pt": _column(jets, "raw_pt", n, "float64", math.nan),
        "jet_pt": _column(jets, "corrected_pt", n, "float64", math.nan),
        "jet_eta": _column(jets, "eta", n, "float64", math.nan),
        "jet_phi": _column(jets, "phi", n, "float64", math.nan),
        "jet_mass": _column(jets, "mass", n, "float64", math.nan),
        "jet_area": _column(jets, "area", n, "float64", math.nan),
        "jet_quality_bitmask": _column(jets, "quality_bitmask", n, "uint64", UINT64_ZERO),
        "jet_order": _column(jets, "deterministic_order", n, "int32", -1),
        "jet_input_identity": inputs,
        "jet_subtraction_identity": subtractions,
        "jet_view_provenance_state": provenance,
    }


def _extend_if_nonempty(tree: uproot.WritableTree, columns: dict[str, np.ndarray]) -> None:
    if columns and len(next(iter(columns.values()))) > 0:
        tree.extend(columns)


def _extend_row_batches(
    tree: uproot.WritableTree,
    rows: Iterable[dict[str, Any]],
    schema: dict[str, str],
    *,
    batch_size: int = ROW_BATCH_SIZE,
) -> int:
    """Serialize flat rows without retaining a source-wide denormalized copy."""
    if batch_size < 1:
        raise ValueError("row batch size must be positive")
    batch: list[dict[str, Any]] = []
    total = 0
    for row in rows:
        batch.append(row)
        if len(batch) >= batch_size:
            _extend_if_nonempty(tree, _to_columns(batch, schema))
            total += len(batch)
            batch.clear()
    if batch:
        _extend_if_nonempty(tree, _to_columns(batch, schema))
        total += len(batch)
    return total


def _jet_response_rows(
    source_index: int,
    jets: dict[str, np.ndarray],
    truth_jets: dict[str, np.ndarray],
    links: dict[str, np.ndarray],
    reco_lookup: dict[tuple[int, int], tuple[tuple[int, int], int]],
    truth_lookup: dict[tuple[int, int], tuple[tuple[int, int], int]],
    available_views: dict[tuple[int, int], list[tuple[str, str, float]]] | None = None,
) -> list[dict[str, Any]]:
    """Recover the final jet partition per captured event/view/radius.

    This is the Python transcription of buildDeterministicJetLinks. It changes
    only the partition scope: each explicit reconstruction view is evaluated
    independently. The frozen candidate gate and tie breakers are unchanged.
    """
    inputs, subtractions, provenance = _jet_view_columns(jets)

    def source_row(i: int, input_name: str, subtraction_name: str, state: int) -> dict[str, Any]:
        reco_type = int(links["reco_type"][i])
        truth_type = int(links["truth_type"][i])
        reco_target = (None if reco_type == 0 else
                       reco_lookup.get(_key(links["reco_id_hi"][i], links["reco_id_lo"][i])))
        truth_target = (None if truth_type == 0 else
                        truth_lookup.get(_key(links["truth_id_hi"][i], links["truth_id_lo"][i])))
        if reco_type not in (0, 2) or truth_type not in (0, 2):
            raise ValueError(f"jet response row {i} has non-jet target type")
        if reco_type and reco_target is None:
            raise ValueError(f"jet response row {i} has no reconstructed jet target")
        if truth_type and truth_target is None:
            raise ValueError(f"jet response row {i} has no truth jet target")
        event = reco_target[0] if reco_target is not None else truth_target[0]
        if reco_target is not None and truth_target is not None and reco_target[0] != truth_target[0]:
            raise ValueError(f"jet response row {i} crosses event boundaries")
        return {
            "source_file_index": source_index,
            "event_id_hi": np.uint64(event[0]), "event_id_lo": np.uint64(event[1]),
            "link_id_hi": links["link_id_hi"][i], "link_id_lo": links["link_id_lo"][i],
            "reco_type": reco_type,
            "reco_id_hi": links["reco_id_hi"][i], "reco_id_lo": links["reco_id_lo"][i],
            "reco_index": -1 if reco_target is None else reco_target[1],
            "truth_type": truth_type,
            "truth_id_hi": links["truth_id_hi"][i], "truth_id_lo": links["truth_id_lo"][i],
            "truth_index": -1 if truth_target is None else truth_target[1],
            "match_metric": links["match_metric"][i], "link_class": int(links["link_class"][i]),
            "jet_input_identity": input_name,
            "jet_subtraction_identity": subtraction_name,
            "jet_response_state": state,
        }

    jet_link_indices = [i for i in range(len(links["link_id_hi"]))
                        if int(links["reco_type"][i]) == 2 or int(links["truth_type"][i]) == 2]
    if len(provenance) and np.any(provenance != provenance[0]):
        raise ValueError("RJJetV1 mixes supported and unsupported view provenance")
    if available_views is None or (len(provenance) and provenance[0] == 0):
        rows = [source_row(i, "", "", 0) for i in jet_link_indices]
        represented_truth = {
            _key(links["truth_id_hi"][i], links["truth_id_lo"][i])
            for i in jet_link_indices if int(links["truth_type"][i]) == 2
        }
        occupied = {_key(links["link_id_hi"][i], links["link_id_lo"][i])
                    for i in range(len(links["link_id_hi"]))}
        next_low = max((identity[1] for identity in occupied), default=0) + 1
        for truth_index in range(len(truth_jets["truth_jet_id_hi"])):
            truth_id = _key(truth_jets["truth_jet_id_hi"][truth_index],
                            truth_jets["truth_jet_id_lo"][truth_index])
            if truth_id in represented_truth:
                continue
            while (14, next_low) in occupied:
                next_low += 1
            event = _key(truth_jets["event_id_hi"][truth_index],
                         truth_jets["event_id_lo"][truth_index])
            rows.append({
                "source_file_index": source_index,
                "event_id_hi": np.uint64(event[0]), "event_id_lo": np.uint64(event[1]),
                "link_id_hi": np.uint64(14), "link_id_lo": np.uint64(next_low),
                "reco_type": 0, "reco_id_hi": UINT64_ZERO, "reco_id_lo": UINT64_ZERO,
                "reco_index": -1, "truth_type": 2,
                "truth_id_hi": np.uint64(truth_id[0]), "truth_id_lo": np.uint64(truth_id[1]),
                "truth_index": truth_lookup[truth_id][1], "match_metric": math.nan,
                "link_class": 2, "jet_input_identity": "",
                "jet_subtraction_identity": "", "jet_response_state": 0,
            })
            occupied.add((14, next_low))
            next_low += 1
        return rows

    jet_by_id = {_key(jets["jet_id_hi"][i], jets["jet_id_lo"][i]): i
                 for i in range(len(jets["jet_id_hi"]))}
    truth_by_id = {_key(truth_jets["truth_jet_id_hi"][i],
                        truth_jets["truth_jet_id_lo"][i]): i
                   for i in range(len(truth_jets["truth_jet_id_hi"]))}
    groups: dict[tuple[tuple[int, int], str, str, float], list[int]] = {
        (event, input_name, subtraction_name, radius): []
        for event, views in available_views.items()
        for input_name, subtraction_name, radius in views
    }
    for i in range(len(jets["jet_id_hi"])):
        event = _key(jets["event_id_hi"][i], jets["event_id_lo"][i])
        key = (event, str(inputs[i]), str(subtractions[i]), float(jets["radius"][i]))
        if key not in groups:
            raise ValueError("captured jet lacks exact per-event view availability witness")
        groups[key].append(i)
    truth_by_event: dict[tuple[int, int], list[int]] = {}
    for i in range(len(truth_jets["truth_jet_id_hi"])):
        event = _key(truth_jets["event_id_hi"][i], truth_jets["event_id_lo"][i])
        truth_by_event.setdefault(event, []).append(i)

    rows: list[dict[str, Any]] = []
    # Candidate edges are the frozen raw diagnostic and keep exact identities.
    for i in jet_link_indices:
        if int(links["link_class"][i]) != 5:
            continue
        if int(links["reco_type"][i]) != 2 or int(links["truth_type"][i]) != 2:
            raise ValueError(f"jet candidate row {i} has invalid typed targets")
        jet = jet_by_id.get(_key(links["reco_id_hi"][i], links["reco_id_lo"][i]))
        if jet is None:
            raise ValueError(f"jet candidate row {i} has no reconstructed jet")
        rows.append(source_row(i, str(inputs[jet]), str(subtractions[jet]), 2))

    occupied = {_key(links["link_id_hi"][i], links["link_id_lo"][i])
                for i in range(len(links["link_id_hi"]))}
    next_low = max((identity[1] for identity in occupied), default=0) + 1

    def generated_row(event, view, radius, reco_index, truth_index, metric, link_class, state=1):
        nonlocal next_low
        identity = (14, next_low)
        while identity in occupied:
            next_low += 1
            identity = (14, next_low)
        occupied.add(identity)
        next_low += 1
        reco_id = ((0, 0) if reco_index is None else
                   _key(jets["jet_id_hi"][reco_index], jets["jet_id_lo"][reco_index]))
        truth_id = ((0, 0) if truth_index is None else
                    _key(truth_jets["truth_jet_id_hi"][truth_index], truth_jets["truth_jet_id_lo"][truth_index]))
        return {
            "source_file_index": source_index,
            "event_id_hi": np.uint64(event[0]), "event_id_lo": np.uint64(event[1]),
            "link_id_hi": np.uint64(identity[0]), "link_id_lo": np.uint64(identity[1]),
            "reco_type": 0 if reco_index is None else 2,
            "reco_id_hi": np.uint64(reco_id[0]), "reco_id_lo": np.uint64(reco_id[1]),
            "reco_index": -1 if reco_index is None else reco_lookup[reco_id][1],
            "truth_type": 0 if truth_index is None else 2,
            "truth_id_hi": np.uint64(truth_id[0]), "truth_id_lo": np.uint64(truth_id[1]),
            "truth_index": -1 if truth_index is None else truth_lookup[truth_id][1],
            "match_metric": metric, "link_class": link_class,
            "jet_input_identity": view[0], "jet_subtraction_identity": view[1],
            "jet_response_state": state,
        }

    # A truth-only event/radius does not prove which reconstruction view ran.
    # Preserve a raw miss when present, otherwise materialize one explicit
    # unsupported row per truth object. A future RunMeta/per-event availability
    # witness can seed an authoritative empty view; cross-event view inventory
    # alone must never do so.
    covered_event_radii = {(key[0], key[3]) for key in groups}
    preserved_unscoped_truth: set[tuple[int, int]] = set()
    for i in jet_link_indices:
        if int(links["link_class"][i]) == 5:
            continue
        reco_id = _key(links["reco_id_hi"][i], links["reco_id_lo"][i])
        truth_id = _key(links["truth_id_hi"][i], links["truth_id_lo"][i])
        reco_index = jet_by_id.get(reco_id) if int(links["reco_type"][i]) == 2 else None
        truth_index = truth_by_id.get(truth_id) if int(links["truth_type"][i]) == 2 else None
        if reco_index is not None:
            event = _key(jets["event_id_hi"][reco_index], jets["event_id_lo"][reco_index])
            radius = float(jets["radius"][reco_index])
        elif truth_index is not None:
            event = _key(truth_jets["event_id_hi"][truth_index], truth_jets["event_id_lo"][truth_index])
            radius = float(truth_jets["radius"][truth_index])
        else:
            continue
        if (event, radius) not in covered_event_radii:
            rows.append(source_row(i, "", "", 0))
            if truth_index is not None:
                preserved_unscoped_truth.add(truth_id)
    for truth_index in range(len(truth_jets["truth_jet_id_hi"])):
        event = _key(truth_jets["event_id_hi"][truth_index], truth_jets["event_id_lo"][truth_index])
        radius = float(truth_jets["radius"][truth_index])
        truth_id = _key(truth_jets["truth_jet_id_hi"][truth_index],
                        truth_jets["truth_jet_id_lo"][truth_index])
        if ((event, radius) not in covered_event_radii and
                truth_id not in preserved_unscoped_truth):
            rows.append(generated_row(event, ("", ""), radius, None, truth_index,
                                      math.nan, 2, state=0))

    for (event, input_name, subtraction_name, radius), reco_indices in groups.items():
        truth_indices = [i for i in truth_by_event.get(event, [])
                         if abs(float(truth_jets["radius"][i]) - radius) <= 1.0e-6]
        nearest_reco = {i: math.nan for i in reco_indices}
        nearest_truth = {i: math.nan for i in truth_indices}
        candidates = []
        for reco_index in reco_indices:
            for truth_index in truth_indices:
                dphi = math.atan2(
                    math.sin(float(jets["phi"][reco_index]) - float(truth_jets["phi"][truth_index])),
                    math.cos(float(jets["phi"][reco_index]) - float(truth_jets["phi"][truth_index])),
                )
                delta_r = math.hypot(
                    float(jets["eta"][reco_index]) - float(truth_jets["eta"][truth_index]), dphi
                )
                if not math.isfinite(delta_r):
                    continue
                if not math.isfinite(nearest_reco[reco_index]) or delta_r < nearest_reco[reco_index]:
                    nearest_reco[reco_index] = delta_r
                if not math.isfinite(nearest_truth[truth_index]) or delta_r < nearest_truth[truth_index]:
                    nearest_truth[truth_index] = delta_r
                if delta_r < JET_MATCH_DELTA_R_MAX:
                    candidates.append((delta_r, reco_index, truth_index))
        candidates.sort(key=lambda edge: (
            edge[0], -float(jets["corrected_pt"][edge[1]]),
            _key(jets["jet_id_hi"][edge[1]], jets["jet_id_lo"][edge[1]]),
            _key(truth_jets["truth_jet_id_hi"][edge[2]], truth_jets["truth_jet_id_lo"][edge[2]]),
        ))
        matched_reco, matched_truth = set(), set()
        view = (input_name, subtraction_name)
        for delta_r, reco_index, truth_index in candidates:
            if reco_index in matched_reco or truth_index in matched_truth:
                continue
            matched_reco.add(reco_index)
            matched_truth.add(truth_index)
            rows.append(generated_row(event, view, radius, reco_index, truth_index, delta_r, 0))
        for reco_index in reco_indices:
            if reco_index not in matched_reco:
                rows.append(generated_row(event, view, radius, reco_index, None,
                                          nearest_reco[reco_index], 1))
        for truth_index in truth_indices:
            if truth_index not in matched_truth:
                rows.append(generated_row(event, view, radius, None, truth_index,
                                          nearest_truth[truth_index], 2))
    return rows


def build(
    input_paths: Iterable[Path],
    output_path: Path,
    system: str,
    require_scaled_bit30: bool = False,
    model_path: Path | None = None,
    model_input_count: int | None = None,
    require_complete_interface: bool = False,
    require_mbd_pmt: bool = False,
    require_centrality_replay: bool = False,
    require_trigger_scalers: bool = False,
    layout: str = "expanded",
) -> dict[str, int]:
    paths, output_path = validate_io_paths(input_paths, output_path)
    schemas = output_schemas(layout)
    compact = layout == "normalized_v1"
    if system not in {"pp", "auau"}:
        raise ValueError("system must be pp or auau")
    expected_model_input_count = 11 if system == "pp" else 14
    if (
        model_path is not None
        and model_input_count is not None
        and model_input_count != expected_model_input_count
    ):
        raise ValueError(
            f"BDT input dimension mismatch for {system}: "
            f"{model_input_count} != {expected_model_input_count}"
        )
    evaluate_bdt = _make_rbdt_evaluator(model_path, model_input_count)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = {name: 0 for name in schemas}
    photon_schema, jet_schema, pair_schema = (
        schemas[key] for key in ("photons", "jets", "photonJets"))

    with uproot.recreate(output_path, compression=uproot.ZSTD(5) if compact else uproot.ZLIB(1)) as out:
        out["storage_layout"] = layout
        source_records = []
        trees = {name: out.mktree(name, schema, title=name)
                 for name, schema in schemas.items()}

        for source_index, path in enumerate(paths):
            with uproot.open(path) as root:
                source_records.append(source_file_record(root, path, source_index))
                event_arrays = _arrays(_tree(root, "RJEventV1"), EVENT_INPUT_BRANCHES)
                if require_scaled_bit30:
                    if "scaled_trigger_bits" not in event_arrays or "scaled_bit30" not in event_arrays:
                        raise ValueError(f"canonical pp input lacks scaled-trigger branches: {path}")
                    reconstructed = ((event_arrays["scaled_trigger_bits"] >> np.uint64(30)) & np.uint64(1)).astype("int32")
                    if not np.all(reconstructed == 1) or not np.all(np.asarray(event_arrays["scaled_bit30"], dtype="int32") == 1):
                        raise ValueError(f"canonical pp input contains an event without scaled bit 30: {path}")
                event_columns = _event_columns(event_arrays, source_index)
                counts["triggerScalers"] += trigger_scaler_interface.copy_scalers(
                    root, trees["triggerScalers"], event_arrays, source_index,
                    require=require_trigger_scalers)
                counts["triggerRunInfo"] += gl1_trigger_interface.copy_run_info(
                    root, trees["triggerRunInfo"], event_arrays, source_index)
                weights = (_arrays(root[PREFIX + "RJWeightComponentV1"], WEIGHT_INPUT_BRANCHES)
                           if PREFIX + "RJWeightComponentV1" in root else None)
                event_columns.update(_interface_columns(
                    event_arrays, weights, system,
                    require_complete=require_complete_interface,
                ))
                if require_centrality_replay and np.any(event_columns["centrality_replay_version"] != 1):
                    raise ValueError("complete centrality replay capture is required")
                mbd_arrays = _mbd_columns(event_arrays)[1]
                mbd_arrays.update(centrality_replay.capture_columns(event_arrays)[1])
                mbd_arrays.update(_truth_vertex_columns(event_arrays)[1])
                mbd_arrays.update(_truth_jet_availability(event_arrays)[1])
                if require_mbd_pmt and np.any(event_columns["mbd_pmt_available"] != 1):
                    raise ValueError("per-PMT MBD unavailable in source; a total cannot supply PMT charges")
                event_lookup: dict[tuple[int, int], int] = {}
                for i in range(len(event_columns["event_id_hi"])):
                    _insert_unique(
                        event_lookup,
                        _key(event_columns["event_id_hi"][i], event_columns["event_id_lo"][i]),
                        i,
                        "event identity",
                    )
                centrality_by_event = {
                    event_key: float(event_columns["centrality"][event_index])
                    for event_key, event_index in event_lookup.items()
                }
                _extend_if_nonempty(trees["events"],
                                    {**event_columns, **mbd_arrays} if compact else event_columns)
                counts["events"] += len(event_columns["event_id_hi"])
                object_events = ({key: event_columns[key] for key in EVENT_JOIN_SCHEMA}
                                 if compact else event_columns)
                available_views, availability_rows = _reco_jet_view_availability(event_arrays)

                def availability_output_rows() -> Iterable[dict[str, Any]]:
                    for row in availability_rows:
                        event_index = row["event_index"]
                        yield {
                            **_row(event_columns, event_index),
                            **{key: value for key, value in row.items() if key != "event_index"},
                        }

                counts["jetViewAvailability"] += _extend_row_batches(
                    trees["jetViewAvailability"], availability_output_rows(),
                    JET_VIEW_AVAILABILITY_SCHEMA,
                )

                candidate_arrays = _arrays(_tree(root, "RJPhotonCandidateV1"), CANDIDATE_INPUT_BRANCHES)
                model_arrays = _arrays(_tree(root, "RJModelEvaluationV1"), MODEL_INPUT_BRANCHES)
                isolation_arrays = _arrays(_tree(root, "RJIsolationWitnessV1"), ISOLATION_INPUT_BRANCHES)
                shower_view_arrays = _arrays(_tree(root, "RJShowerFeatureViewV1"), SHOWER_VIEW_INPUT_BRANCHES)
                photon_rows, photon_lookup = _photon_rows(
                    candidate_arrays,
                    model_arrays,
                    isolation_arrays,
                    shower_view_arrays,
                    centrality_by_event,
                    system,
                    evaluate_bdt,
                    model_input_count,
                    preserve_capture=compact,
                )
                retained_candidate_keys = set(photon_lookup)
                candidate_arrays = _filter_by_identity(
                    candidate_arrays,
                    "candidate_id_hi",
                    "candidate_id_lo",
                    retained_candidate_keys,
                )
                del model_arrays, isolation_arrays, shower_view_arrays
                gc.collect()

                def photon_output_rows() -> Iterable[dict[str, Any]]:
                    for i, row in enumerate(photon_rows):
                        event_key = _key(candidate_arrays["event_id_hi"][i], candidate_arrays["event_id_lo"][i])
                        if event_key not in event_lookup:
                            raise ValueError(f"photon row {i} has no event join in {path}")
                        yield {**_row(object_events, event_lookup[event_key]), **row}

                counts["photons"] += _extend_row_batches(
                    trees["photons"], photon_output_rows(), photon_schema
                )

                jet_arrays = _arrays(_tree(root, "RJJetV1"), JET_INPUT_BRANCHES)
                jet_columns = _jet_columns(jet_arrays)

                def jet_output_rows() -> Iterable[dict[str, Any]]:
                    for i in range(len(jet_arrays["jet_id_hi"])):
                        event_key = _key(jet_arrays["event_id_hi"][i], jet_arrays["event_id_lo"][i])
                        if event_key not in event_lookup:
                            raise ValueError(f"jet row {i} has no event join in {path}")
                        yield {
                            **_row(object_events, event_lookup[event_key]),
                            **_row(jet_columns, i),
                        }

                counts["jets"] += _extend_row_batches(trees["jets"], jet_output_rows(), jet_schema)

                pair_arrays = _filter_by_identity(
                    _arrays(_tree(root, "RJPhotonJetPairV1"), PAIR_INPUT_BRANCHES),
                    "candidate_id_hi",
                    "candidate_id_lo",
                    retained_candidate_keys,
                )
                pair_n = len(pair_arrays["pair_id_hi"])
                wrong_photon = pair_arrays.get("wrong_photon_class", np.full(pair_n, -1, dtype="int32"))
                wrong_recoil = pair_arrays.get("wrong_recoil_class", np.full(pair_n, -1, dtype="int32"))
                (
                    pair_photon_source,
                    pair_jet_source,
                    pair_photon_local,
                    pair_jet_local,
                ) = _pair_join_indices(candidate_arrays, jet_arrays, pair_arrays)

                def pair_output_rows() -> Iterable[dict[str, Any]]:
                    for i in range(pair_n):
                        event_key = _key(pair_arrays["event_id_hi"][i], pair_arrays["event_id_lo"][i])
                        if event_key not in event_lookup:
                            raise ValueError(f"pair row {i} has no event join in {path}")
                        photon_source_index = int(pair_photon_source[i])
                        jet_source_index = int(pair_jet_source[i])
                        pair_row = {
                            "pair_id_hi": pair_arrays["pair_id_hi"][i],
                            "pair_id_lo": pair_arrays["pair_id_lo"][i],
                            "photon_index": pair_photon_local[i],
                            "jet_index": pair_jet_local[i],
                            "delta_phi": float(pair_arrays["delta_phi"][i]),
                            "xjgamma": float(pair_arrays["xjgamma"][i]),
                            "recoil_state": int(pair_arrays["recoil_state"][i]),
                            "photon_rank": int(pair_arrays["photon_rank"][i]),
                            "jet_rank": int(pair_arrays["jet_rank"][i]),
                            "wrong_photon_class": int(wrong_photon[i]),
                            "wrong_recoil_class": int(wrong_recoil[i]),
                        }
                        yield {
                            **_row(object_events, event_lookup[event_key]),
                            **({key: pair_arrays[key][i] for key in PAIR_OBJECT_IDS}
                               if compact else {**photon_rows[photon_source_index],
                                                **_row(jet_columns, jet_source_index)}),
                            **pair_row,
                        }

                counts["photonJets"] += _extend_row_batches(
                    trees["photonJets"], pair_output_rows(), pair_schema
                )

                event_context = (None if compact else
                                 _event_tree_context(candidate_arrays, jet_arrays, pair_arrays))
                event_n = len(event_columns["event_id_hi"])
                for event_start in range(0, 0 if compact else event_n, EVENT_BATCH_SIZE):
                    event_stop = min(event_start + EVENT_BATCH_SIZE, event_n)
                    event_tree_columns = _event_tree_columns(
                        event_columns,
                        candidate_arrays,
                        photon_rows,
                        jet_arrays,
                        jet_columns,
                        pair_arrays,
                        event_start=event_start,
                        event_stop=event_stop,
                        context=event_context,
                        pair_photon_local=pair_photon_local,
                        pair_jet_local=pair_jet_local,
                    )
                    event_tree_columns.update({name: values[event_start:event_stop]
                                               for name, values in mbd_arrays.items()})
                    _extend_if_nonempty(trees["eventTree"], event_tree_columns)
                    counts["eventTree"] += event_stop - event_start
                    del event_tree_columns

                del event_context, pair_arrays, wrong_photon, wrong_recoil
                del pair_photon_source, pair_jet_source, pair_photon_local, pair_jet_local
                del photon_rows, photon_lookup, jet_columns
                gc.collect()

                truth_photons = _arrays(_tree(root, "RJTruthPhotonV1"), TRUTH_PHOTON_INPUT_BRANCHES)
                truth_witnesses = {
                    field: _column(truth_photons, field, len(truth_photons["truth_photon_id_hi"]),
                        dtype, math.nan if dtype == "float64" else -1)
                    for field, dtype in TRUTH_PHOTON_SCHEMA.items()
                    if field in ("truth_isolation_r03", "truth_isolation_r04", "truth_isolation_valid",
                        "g4_photon_valid", "hepmc_association_valid", "analysis_signal_r03",
                        "native_track_id", "native_vertex_id", "embedding_id",
                        "sample_source_role", "generator_occurrence_embedding_id")
                }

                def truth_photon_rows() -> Iterable[dict[str, Any]]:
                    for i in range(len(truth_photons["truth_photon_id_hi"])):
                        yield {
                        "source_file_index": source_index,
                        "event_id_hi": truth_photons["event_id_hi"][i],
                        "event_id_lo": truth_photons["event_id_lo"][i],
                        "truth_photon_id_hi": truth_photons["truth_photon_id_hi"][i],
                        "truth_photon_id_lo": truth_photons["truth_photon_id_lo"][i],
                        "truth_photon_pt": truth_photons["pt"][i],
                        "truth_photon_eta": truth_photons["eta"][i],
                        "truth_photon_phi": truth_photons["phi"][i],
                        "prompt_class": truth_photons["prompt_class"][i],
                        "source_role": truth_photons["source_role"][i],
                        "generator_barcode": truth_photons["generator_barcode"][i],
                        "truth_isolation": truth_photons["truth_isolation_witness"][i],
                        **{field: values[i] for field, values in truth_witnesses.items()},
                        }

                counts["truthPhotons"] += _extend_row_batches(
                    trees["truthPhotons"], truth_photon_rows(), TRUTH_PHOTON_SCHEMA
                )

                truth_jets = _arrays(_tree(root, "RJTruthJetV1"), TRUTH_JET_INPUT_BRANCHES)

                def truth_jet_rows() -> Iterable[dict[str, Any]]:
                    for i in range(len(truth_jets["truth_jet_id_hi"])):
                        yield {
                        "source_file_index": source_index,
                        "event_id_hi": truth_jets["event_id_hi"][i],
                        "event_id_lo": truth_jets["event_id_lo"][i],
                        "truth_jet_id_hi": truth_jets["truth_jet_id_hi"][i],
                        "truth_jet_id_lo": truth_jets["truth_jet_id_lo"][i],
                        "truth_jet_radius": truth_jets["radius"][i],
                        "truth_jet_pt": truth_jets["pt"][i],
                        "truth_jet_eta": truth_jets["eta"][i],
                        "truth_jet_phi": truth_jets["phi"][i],
                        }

                counts["truthJets"] += _extend_row_batches(
                    trees["truthJets"], truth_jet_rows(), TRUTH_JET_SCHEMA
                )

                links = _arrays(_tree(root, "RJRecoTruthLinkV1"), LINK_INPUT_BRANCHES)
                if len(links["link_id_hi"]):
                    reco_lookup_by_type = {
                        1: _local_identity_lookup(
                            candidate_arrays,
                            "candidate_id_hi",
                            "candidate_id_lo",
                            "reco photon identity",
                        ),
                        2: _local_identity_lookup(
                            jet_arrays, "jet_id_hi", "jet_id_lo", "reco jet identity"
                        ),
                    }
                    truth_lookup_by_type = {
                        1: _local_identity_lookup(
                            truth_photons,
                            "truth_photon_id_hi",
                            "truth_photon_id_lo",
                            "truth photon identity",
                        ),
                        2: _local_identity_lookup(
                            truth_jets,
                            "truth_jet_id_hi",
                            "truth_jet_id_lo",
                            "truth jet identity",
                        ),
                    }
                else:
                    reco_lookup_by_type = {}
                    truth_lookup_by_type = {}
                def link_rows() -> Iterable[dict[str, Any]]:
                    for i in range(len(links["link_id_hi"])):
                        reco_type = int(links["reco_type"][i])
                        truth_type = int(links["truth_type"][i])
                        if reco_type == 1 and _key(
                            links["reco_id_hi"][i], links["reco_id_lo"][i]
                        ) not in retained_candidate_keys:
                            continue
                        reco_target = None
                        truth_target = None
                        if reco_type != 0:
                            reco_target = reco_lookup_by_type[reco_type].get(
                                _key(links["reco_id_hi"][i], links["reco_id_lo"][i])
                            )
                            if reco_target is None:
                                raise ValueError(f"truth link row {i} has no reconstructed target")
                        if truth_type != 0:
                            truth_target = truth_lookup_by_type[truth_type].get(
                                _key(links["truth_id_hi"][i], links["truth_id_lo"][i])
                            )
                            if truth_target is None:
                                raise ValueError(f"truth link row {i} has no truth target")
                        event_key = reco_target[0] if reco_target is not None else truth_target[0]
                        if reco_target is not None and truth_target is not None and reco_target[0] != truth_target[0]:
                            raise ValueError(f"truth link row {i} crosses event boundaries")
                        yield {
                            "source_file_index": source_index,
                            "event_id_hi": np.uint64(event_key[0]),
                            "event_id_lo": np.uint64(event_key[1]),
                            "link_id_hi": links["link_id_hi"][i],
                            "link_id_lo": links["link_id_lo"][i],
                            "reco_type": reco_type,
                            "reco_id_hi": links["reco_id_hi"][i],
                            "reco_id_lo": links["reco_id_lo"][i],
                            "reco_index": -1 if reco_target is None else reco_target[1],
                            "truth_type": truth_type,
                            "truth_id_hi": links["truth_id_hi"][i],
                            "truth_id_lo": links["truth_id_lo"][i],
                            "truth_index": -1 if truth_target is None else truth_target[1],
                            "match_metric": links["match_metric"][i],
                            "link_class": links["link_class"][i],
                        }

                counts["recoTruthLinks"] += _extend_row_batches(
                    trees["recoTruthLinks"], link_rows(), LINK_SCHEMA
                )
                jet_response_rows = _jet_response_rows(
                    source_index,
                    jet_arrays,
                    truth_jets,
                    links,
                    _local_identity_lookup(
                        jet_arrays, "jet_id_hi", "jet_id_lo", "reco jet identity"
                    ),
                    _local_identity_lookup(
                        truth_jets, "truth_jet_id_hi", "truth_jet_id_lo", "truth jet identity"
                    ),
                    available_views,
                )
                counts["jetResponseLinks"] += _extend_row_batches(
                    trees["jetResponseLinks"], jet_response_rows, JET_RESPONSE_SCHEMA
                )
                del links, truth_photons, truth_jets
                del candidate_arrays, jet_arrays, event_arrays, event_columns, event_lookup
                del reco_lookup_by_type, truth_lookup_by_type, retained_candidate_keys
                del jet_response_rows
                del available_views, availability_rows
                gc.collect()

        out["source_files"] = json.dumps(source_records, sort_keys=True, separators=(",", ":"))
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", default=[], type=Path, help="producer ROOT output file; repeatable")
    parser.add_argument("--input-list", action="append", default=[], type=Path, help="text file containing one input ROOT path per line")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--model", type=Path, help="TMVA RBDT used to calculate every collaborator-facing bdt_score")
    parser.add_argument("--model-input-count", type=int, choices=(11, 14), help="exact ordered input dimension for --model")
    parser.add_argument("--require-scaled-bit30", action="store_true", help="require canonical pp scaled-trigger branches and membership")
    parser.add_argument("--require-complete-interface", action="store_true",
                        help="require producer-certified sample factors and complete truth denominators")
    parser.add_argument("--require-mbd-pmt", action="store_true")
    parser.add_argument("--require-centrality-replay", action="store_true")
    parser.add_argument("--require-trigger-scalers", action="store_true",
                        help="require valid lossless source-bound GL1 scaler capture for DATA")
    parser.add_argument("--layout", choices=LAYOUTS, default="expanded",
                        help="normalized_v1 stores each payload once and retains unscored capture; expanded is the legacy convenience view")
    args = parser.parse_args()
    inputs = list(args.input)
    for list_path in args.input_list:
        inputs.extend(
            Path(line.strip())
            for line in list_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    counts = build(
        inputs, args.output, args.system, args.require_scaled_bit30,
        args.model, args.model_input_count, args.require_complete_interface, args.require_mbd_pmt,
        args.require_centrality_replay, args.require_trigger_scalers,
        args.layout,
    )
    print("OUTPUT", args.output.resolve())
    for name, count in counts.items():
        print(f"{name}={count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
