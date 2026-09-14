"""Write small collaborator-style tree files with the branches the stages read.

The engineering fixtures predate the truth-isolation branches of the current
tree contract, so tests of the nominal truth-signal paths (photon response,
sample weights, training labels) build their own inputs here.  Only the
branches consumed by the code under test are written; these files are not
full contract files and are never validated against contracts/.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import uproot


def _columns(rows: list[dict[str, Any]], schema: dict[str, str]) -> dict[str, np.ndarray]:
    return {name: np.asarray([row.get(name, 0) for row in rows], dtype=dtype) for name, dtype in schema.items()}


EVENT_SCHEMA = {
    "source_file_index": "int32", "event_id_hi": "uint64", "event_id_lo": "uint64", "run": "int32",
    "event_weight": "float64", "centrality": "float64", "vertex_z": "float64", "terminal_status": "int32",
}
PHOTON_SCHEMA = {
    **EVENT_SCHEMA,
    "candidate_id_hi": "uint64", "candidate_id_lo": "uint64", "photon_encounter_ordinal": "int32",
    "photon_et": "float64", "photon_eta": "float64", "photon_phi": "float64",
    "bdt_score": "float64", "bdt_tight_threshold": "float64", "bdt_nontight_low_threshold": "float64",
    "bdt_nontight_high_threshold": "float64", "bdt_is_tight": "int32", "bdt_is_nontight": "int32",
    "bdt_is_not_tight": "int32", "bdt_input_count": "int32",
    "iso_r03": "float64", "iso_r03_threshold": "float64", "iso_r03_nonisolated_threshold": "float64", "iso_r03_pass": "int32",
    "iso_r04": "float64", "iso_r04_threshold": "float64", "iso_r04_nonisolated_threshold": "float64", "iso_r04_pass": "int32",
    "truth_matched": "int32", "truth_barcode": "int32",
    **{f"bdt_input_{i:02d}": "float64" for i in range(14)},
}
PAIR_SCHEMA = {
    **PHOTON_SCHEMA,
    "pair_id_hi": "uint64", "pair_id_lo": "uint64", "jet_id_hi": "uint64", "jet_id_lo": "uint64",
    "photon_index": "int32", "jet_index": "int32", "jet_pt": "float64", "jet_eta": "float64", "jet_phi": "float64",
    "jet_radius": "float64", "delta_phi": "float64", "xjgamma": "float64", "recoil_state": "int32",
}
TRUTH_PHOTON_SCHEMA = {
    "source_file_index": "int32", "event_id_hi": "uint64", "event_id_lo": "uint64",
    "truth_photon_id_hi": "uint64", "truth_photon_id_lo": "uint64",
    "truth_photon_pt": "float64", "truth_photon_eta": "float64", "truth_photon_phi": "float64",
    "prompt_class": "int32", "source_role": "int32", "generator_barcode": "int32", "truth_isolation": "float64",
    "truth_isolation_r03": "float64", "truth_isolation_r04": "float64", "truth_isolation_valid": "int32",
    "g4_photon_valid": "int32", "hepmc_association_valid": "int32", "analysis_signal_r03": "int32",
}
TRUTH_JET_SCHEMA = {
    "source_file_index": "int32", "event_id_hi": "uint64", "event_id_lo": "uint64",
    "truth_jet_id_hi": "uint64", "truth_jet_id_lo": "uint64",
    "truth_jet_radius": "float64", "truth_jet_pt": "float64", "truth_jet_eta": "float64", "truth_jet_phi": "float64",
}
LINK_SCHEMA = {
    "source_file_index": "int32", "event_id_hi": "uint64", "event_id_lo": "uint64",
    "link_id_hi": "uint64", "link_id_lo": "uint64", "reco_id_hi": "uint64", "reco_id_lo": "uint64",
    "truth_id_hi": "uint64", "truth_id_lo": "uint64", "reco_type": "int32", "truth_type": "int32",
    "reco_index": "int32", "truth_index": "int32", "link_class": "int32", "match_metric": "float64",
}


def write_tree_file(
    path: Path,
    *,
    events: list[dict[str, Any]],
    photons: list[dict[str, Any]],
    pairs: list[dict[str, Any]] = (),
    truth_photons: list[dict[str, Any]] = (),
    truth_jets: list[dict[str, Any]] = (),
    links: list[dict[str, Any]] = (),
) -> Path:
    """Photon rows inherit their event's columns; pair rows inherit their photon's."""

    by_event = {(r["source_file_index"], r["event_id_hi"], r["event_id_lo"]): r for r in events}
    photon_rows = [{**by_event[(p["source_file_index"], p["event_id_hi"], p["event_id_lo"])], **p} for p in photons]
    by_candidate = {(r["source_file_index"], r["candidate_id_hi"], r["candidate_id_lo"]): r for r in photon_rows}
    pair_rows = [{**by_candidate[(q["source_file_index"], q["candidate_id_hi"], q["candidate_id_lo"])], **q} for q in pairs]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(path) as root:
        root["events"] = _columns(events, EVENT_SCHEMA)
        root["photons"] = _columns(photon_rows, PHOTON_SCHEMA)
        root["photonJets"] = _columns(pair_rows, PAIR_SCHEMA)
        root["truthPhotons"] = _columns(list(truth_photons), TRUTH_PHOTON_SCHEMA)
        root["truthJets"] = _columns(list(truth_jets), TRUTH_JET_SCHEMA)
        root["recoTruthLinks"] = _columns(list(links), LINK_SCHEMA)
    return path


def signal_truth_photon(source: int, hi: int, lo: int, truth_hi: int, truth_lo: int, pt: float, **overrides: Any) -> dict[str, Any]:
    """A truth photon satisfying the nominal contract unless overridden."""

    row = {
        "source_file_index": source, "event_id_hi": hi, "event_id_lo": lo,
        "truth_photon_id_hi": truth_hi, "truth_photon_id_lo": truth_lo,
        "truth_photon_pt": pt, "truth_photon_eta": 0.1, "truth_photon_phi": 0.5,
        "prompt_class": 1, "source_role": 1, "generator_barcode": 7, "truth_isolation": 1.0,
        "truth_isolation_r03": 1.0, "truth_isolation_r04": 1.5, "truth_isolation_valid": 1,
        "g4_photon_valid": 1, "hepmc_association_valid": 1, "analysis_signal_r03": 1,
    }
    row.update(overrides)
    return row


def match_link(source: int, hi: int, lo: int, reco_hi: int, reco_lo: int, truth_hi: int, truth_lo: int, metric: float = 0.01) -> dict[str, Any]:
    return {
        "source_file_index": source, "event_id_hi": hi, "event_id_lo": lo,
        "link_id_hi": reco_hi, "link_id_lo": reco_lo + 1000, "reco_id_hi": reco_hi, "reco_id_lo": reco_lo,
        "truth_id_hi": truth_hi, "truth_id_lo": truth_lo, "reco_type": 1, "truth_type": 1,
        "reco_index": 0, "truth_index": 0, "link_class": 0, "match_metric": metric,
    }
