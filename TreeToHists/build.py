#!/usr/bin/env python3
"""Build additive histogram statistics for one logical sample family.

One event-batched reader joins exact-base score sidecars by candidate identity
and fans out ABCD/recoil, leakage, responses, fakes/misses and boundary counts.
Sample/period/SI-DI partitions remain additive inside the family package.
Measurement settings and normalization must be explicitly approved; diagnostic
packages are labelled and cannot enter FinalAnalysis. See README.md for the
retained covariance scope and unresolved response population contract.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import sqlite3
import tempfile
from contextlib import closing, contextmanager
from pathlib import Path
import sys
from typing import Any, Mapping

import numpy as np
import uproot
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from histogram_contract import (  # noqa: E402
    ASSIGNED_HARD_NONFIDUCIAL, COMBINATORIC, Grids, Package, REGION_ORDER, UNMATCHED_RECO, axis_bin,
    file_sha256 as sha256,
)

# ---- enumerators persisted by the producer (internal/Types.h) -----------------------
JET_VIEW = {"PP": 1, "AuAuSub1": 2, "AuAuNoSub": 3}
ISOLATION_METHOD = {"calorimeter_raw": 1, "calorimeter_sub1": 2, "topocluster": 3}
ASSOCIATION_MATCHED = 3
TRISTATE_TRUE = 1
CAPTURE_COMPLETE = 1
SCORE_EVALUATED_FINITE = 1

EventKey = tuple[int, int]
PhotonKey = tuple[int, int, int, int]


# ============================================================================
# Sample weights (the one implementation of w = event_weight x source x centrality)
# ============================================================================

SOURCE_FACTORS = ("family_reference", "cross_section_over_generated", "ownership_stitch")


class SampleManifest:
    def __init__(self, path: Path) -> None:
        raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        self.path = str(path)
        self.campaign = raw.get("campaign")
        self.sha256 = sha256(Path(path))
        self.status = raw.get("status", "unresolved")
        references = raw.get("reference_cross_section_pb", {}) or {}
        self.samples: dict[str, dict[str, Any]] = {}
        for name, spec in (raw.get("samples") or {}).items():
            entry = {"name": name, "system": spec["system"], "kind": spec["kind"], "family": spec["family"],
                     "source_factor": None, "factor": None, "window": None, "note": "",
                     "evidence": spec.get("evidence"), "source_manifest_sha256": spec.get("source_manifest_sha256"),
                     "mixture_factor": spec.get("mixture_factor"), "mixture_evidence": spec.get("mixture_evidence"),
                     "approved": spec.get("status") == "approved"}
            if spec["kind"] == "simulation":
                kind = spec.get("source_factor")
                if kind not in SOURCE_FACTORS:
                    raise ValueError(f"{name}: source_factor must be one of {SOURCE_FACTORS}")
                entry["source_factor"] = kind
                if kind == "family_reference":
                    if spec.get("generated_events") is None:
                        entry["note"] = "generated_events unresolved: histogram construction forbidden"
                    else:
                        if references.get(spec["reference"]) is not None:
                            entry["factor"] = float(references[spec["reference"]]) / float(spec["generated_events"])
                elif kind == "cross_section_over_generated":
                    if spec.get("cross_section_pb") is not None and spec.get("generated_events") is not None:
                        entry["factor"] = float(spec["cross_section_pb"]) / float(spec["generated_events"])
                else:
                    low, high = spec["ownership_window_gev"]
                    entry["window"] = (float(low), None if high is None else float(high))
                    if spec.get("ownership_effective_cross_section_pb") is not None and spec.get("owned_events") is not None:
                        entry["factor"] = float(spec["ownership_effective_cross_section_pb"]) / float(spec["owned_events"])
                if entry["factor"] is not None and not (math.isfinite(entry["factor"]) and entry["factor"] > 0):
                    raise ValueError(f"{name}: source factor must be finite and positive")
            self.samples[name] = entry
        block = raw.get("centrality_weights")
        self.centrality_edges = None
        self.centrality_maps: dict[str, np.ndarray] = {}
        if block and block.get("status") == "approved" and block.get("evidence"):
            if block.get("campaign") != self.campaign:
                raise ValueError("centrality map and sample campaign differ")
            self.centrality_edges = np.asarray(block["bin_edges_percent"], dtype=float)
            if np.any(~np.isfinite(self.centrality_edges)) or np.any(np.diff(self.centrality_edges) <= 0):
                raise ValueError("centrality edges must be finite and increasing")
            for family in ("photonjet", "inclusivejet"):
                if family in block:
                    table = np.asarray(block[family], dtype=float)
                    if table.shape != (len(self.centrality_edges) - 1,) or np.any(~np.isfinite(table)) or np.any(table <= 0):
                        raise ValueError(f"centrality map {family} is malformed")
                    self.centrality_maps[family] = table

    def centrality_factor(self, family: str, centrality: np.ndarray) -> np.ndarray:
        if self.centrality_edges is None or family not in self.centrality_maps:
            raise ValueError(f"no centrality map for family {family!r} in {self.path}")
        c = np.asarray(centrality, dtype=float)
        if np.any(~np.isfinite(c)):
            raise ValueError("Au+Au simulation event has a non-finite centrality")
        index = np.searchsorted(self.centrality_edges, c, side="right") - 1
        inside = (c >= self.centrality_edges[0]) & (c < self.centrality_edges[-1])
        out = np.full(c.shape, np.nan)
        out[inside] = self.centrality_maps[family][index[inside]]
        return out


def complete_event_weights(manifest: SampleManifest, sample_name: str, *, event_weight: np.ndarray,
                           event_weight_valid: np.ndarray, centrality: np.ndarray,
                           leading_truth_jet_pt: np.ndarray | None) -> dict[str, Any]:
    sample = manifest.samples.get(sample_name)
    if sample is None:
        raise KeyError(f"sample {sample_name!r} is not in {manifest.path}")
    if manifest.status != "approved" or not manifest.campaign or not sample["approved"] or not sample["evidence"]:
        raise ValueError(f"{sample_name}: sample normalization contract is not approved")
    w = np.asarray(event_weight, dtype=float).copy()
    valid = np.asarray(event_weight_valid, dtype=bool)
    if np.any(valid & (~np.isfinite(w) | (w < 0))):
        raise ValueError(f"{sample_name}: a valid producer event_weight is not finite and non-negative")
    if not valid.all():
        raise ValueError(f"{sample_name}: producer weights are incomplete")
    components = ["producer_event_weight"]
    dropped = {"invalid_producer_weight": int((~valid).sum()), "outside_ownership_window": 0, "outside_centrality_support": 0}
    if sample["kind"] == "simulation":
        if sample["factor"] is None:
            raise ValueError(f"{sample_name}: missing generated/owned count or cross-section binding")
        if sample["mixture_factor"] is None or not sample["mixture_evidence"]:
            raise ValueError(f"{sample_name}: SI/DI mixture ownership is unresolved")
        mixture = float(sample["mixture_factor"])
        if not math.isfinite(mixture) or mixture < 0:
            raise ValueError("invalid mixture factor")
        w *= mixture
        components.append(f"mixture:{sample_name}:{mixture}")
        if sample["source_factor"] == "ownership_stitch":
            if leading_truth_jet_pt is None:
                raise ValueError(f"{sample_name}: ownership stitching needs the leading truth-jet pT per event")
            pt = np.asarray(leading_truth_jet_pt, dtype=float)
            low, high = sample["window"]
            owned = np.isfinite(pt) & (pt >= low) & (True if high is None else (pt < high))
            dropped["outside_ownership_window"] = int((~owned & np.isfinite(w)).sum())
            w[~owned] = np.nan
            w[owned] *= sample["factor"]
            components.append(f"source_stitch:{sample_name}:{sample['factor']:.8g}pb")
        elif sample["factor"] is not None:
            w *= sample["factor"]
            components.append(f"source_factor:{sample_name}:{sample['source_factor']}={sample['factor']:.8g}pb")
        else:
            raise ValueError("source factor was not applied")
        if sample["system"] == "auau":
            factor = manifest.centrality_factor(sample["family"], centrality)
            alive = np.isfinite(w)
            outside = alive & ~np.isfinite(factor)
            dropped["outside_centrality_support"] = int(outside.sum())
            w[outside] = np.nan
            w[alive & ~outside] *= factor[alive & ~outside]
            components.append(f"centrality:{sample['family']}:{manifest.campaign}")
    return {"weight": w, "kept": np.isfinite(w), "components": components, "dropped": dropped, "family": sample["family"]}


# ============================================================================
# Working points and classes (strict boundaries)
# ============================================================================

def threshold_values(spec: Any, photon_et: np.ndarray, centrality: np.ndarray) -> np.ndarray:
    et = np.asarray(photon_et, dtype=float)
    cent = np.asarray(centrality, dtype=float)
    if spec is None:
        return np.full(et.shape, np.nan)
    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        return np.full(et.shape, float(spec))
    if not isinstance(spec, Mapping):
        raise ValueError(f"unsupported threshold specification: {spec!r}")
    if spec.get("mode") == "binned":
        pt_edges = np.asarray(spec["pt_edges"], dtype=float)
        values = np.asarray(spec["values"], dtype=float)
        if len(pt_edges) < 2 or np.any(~np.isfinite(pt_edges)) or np.any(np.diff(pt_edges) <= 0) or np.any(~np.isfinite(values)):
            raise ValueError("threshold grid/values are invalid")
        pt_bin = np.searchsorted(pt_edges, et, side="right") - 1
        inside = (pt_bin >= 0) & (pt_bin < len(pt_edges) - 1)
        out = np.full(et.shape, np.nan)
        if "cent_edges" in spec:
            cent_edges = np.asarray(spec["cent_edges"], dtype=float)
            if values.shape != (len(pt_edges) - 1, len(cent_edges) - 1):
                raise ValueError("binned threshold values do not match pt_edges x cent_edges")
            if np.any(~np.isfinite(cent_edges)) or np.any(np.diff(cent_edges) <= 0):
                raise ValueError("centrality threshold edges are invalid")
            cent_bin = np.searchsorted(cent_edges, cent, side="right") - 1
            inside &= (cent_bin >= 0) & (cent_bin < len(cent_edges) - 1)
            out[inside] = values[pt_bin[inside], cent_bin[inside]]
        else:
            if values.shape != (len(pt_edges) - 1,):
                raise ValueError("binned threshold values do not match pt_edges")
            out[inside] = values[pt_bin[inside]]
        return out
    if spec.get("mode", "linear") != "linear" or "intercept" not in spec:
        raise ValueError("threshold requires an explicit linear or binned definition")
    out = float(spec["intercept"]) + float(spec.get("slope_photon_et", 0.0)) * et
    if float(spec.get("slope_centrality", 0.0)) != 0.0:
        out = out + float(spec["slope_centrality"]) * cent
    return out


def working_point_thresholds(block: Mapping[str, Any], photon_et: np.ndarray, centrality: np.ndarray) -> dict[str, np.ndarray]:
    wp = block.get("working_points", {}) or {}
    return {
        "tight": threshold_values(wp.get("tight_threshold"), photon_et, centrality),
        "nontight_low": threshold_values(wp.get("nontight_low"), photon_et, centrality),
        "nontight_high": threshold_values(wp.get("nontight_high"), photon_et, centrality),
        "isolated_max": threshold_values(wp.get("isolated_max_gev"), photon_et, centrality),
        "nonisolated_min": threshold_values(wp.get("nonisolated_min_gev"), photon_et, centrality),
    }


def missing_working_points(block: Mapping[str, Any]) -> list[str]:
    wp = block.get("working_points", {}) or {}
    return [k for k in ("tight_threshold", "nontight_low", "nontight_high", "isolated_max_gev", "nonisolated_min_gev") if wp.get(k) is None]


def classify(score: np.ndarray, isolation: np.ndarray, thresholds: Mapping[str, np.ndarray], non_tight_definition: str) -> dict[str, np.ndarray]:
    """Region membership with strict boundaries; invalid inputs belong to no region."""

    s = np.asarray(score, dtype=float); iso = np.asarray(isolation, dtype=float)
    t = thresholds["tight"]; lo = thresholds["nontight_low"]; hi = thresholds["nontight_high"]
    iso_max = thresholds["isolated_max"]; noniso_min = thresholds["nonisolated_min"]
    id_valid = np.isfinite(s) & np.isfinite(t) & np.isfinite(lo) & np.isfinite(hi)
    iso_valid = np.isfinite(iso) & np.isfinite(iso_max) & np.isfinite(noniso_min)
    with np.errstate(invalid="ignore"):
        tight = id_valid & (s > t)
        bounded = id_valid & (s > lo) & (s <= hi)
        complement = id_valid & ~(s > t)
        isolated = iso_valid & (iso < iso_max)
        nonisolated = iso_valid & (iso > noniso_min)
    if non_tight_definition not in ("bounded", "complement"):
        raise ValueError("unknown non-tight definition")
    if np.any(id_valid & (lo >= hi)) or np.any(id_valid & (hi > t)) or np.any(iso_valid & (iso_max > noniso_min)):
        raise ValueError("ABCD regions overlap or have invalid threshold order")
    non_tight = bounded if non_tight_definition == "bounded" else complement
    return {"A": tight & isolated, "B": tight & nonisolated, "C": non_tight & isolated, "D": non_tight & nonisolated,
            "id_valid": id_valid, "iso_valid": iso_valid, "tight": tight}


def event_leading_index(photon_et: np.ndarray, encounter_ordinal: np.ndarray, eligible: np.ndarray) -> int:
    best = -1
    for index in np.flatnonzero(np.asarray(eligible, dtype=bool)):
        if best < 0 or photon_et[index] > photon_et[best] or (
            photon_et[index] == photon_et[best] and encounter_ordinal[index] < encounter_ordinal[best]):
            best = int(index)
    return best


# ============================================================================
# Reading one canonical file
# ============================================================================

def _pairs(arrays: Mapping[str, np.ndarray], prefix: str) -> list[tuple[int, int]]:
    return list(zip(np.asarray(arrays[f"{prefix}_hi"]).astype(np.uint64).tolist(),
                    np.asarray(arrays[f"{prefix}_lo"]).astype(np.uint64).tolist()))


# Columns are deliberate: heavy cells, shower witnesses and Cartesian pair
# tables are not needed by this loop. Pair kinematics come from the retained
# collections using the explicit producer formula, pending equivalence closure.
TABLE_COLUMNS = {
    "Photons": "event_hi event_lo photon_hi photon_lo encounter_ordinal et eta phi",
    "Isolation": "event_hi event_lo photon_hi photon_lo radius method cone_sum valid",
    "Jets": "event_hi event_lo jet_hi jet_lo view radius corrected_pt eta phi calibration_valid",
    "TruthPhotons": "event_hi event_lo truth_photon_hi truth_photon_lo pt eta phi geant_valid generator_association_valid prompt_class isolation_r03 isolation_r04 isolation_valid analysis_signal",
    "TruthJets": "event_hi event_lo truth_jet_hi truth_jet_lo radius pt eta phi",
    "PhotonTruthLinks": "event_hi event_lo photon_hi photon_lo truth_photon_hi truth_photon_lo state",
    "JetTruthLinks": "event_hi event_lo jet_hi jet_lo truth_jet_hi truth_jet_lo jet_view jet_radius state selected_match",
}
EVENT_COLUMNS = "photon_count jet_count truth_photon_count truth_jet_count trigger_live_bits trigger_scaled_bits trigger_packet_valid trigger_packet_status trigger_decisions_available source_hi source_lo event_hi event_lo terminal_status reco_vertex_z reco_vertex_valid reco_object_vertex_in_domain centrality_percent centrality_valid minimum_bias_decision event_weight event_weight_valid reco_photon_capture_state truth_denominator_complete reco_jet_view reco_jet_radius_code reco_jet_capture_state truth_jet_radius_code truth_jet_container_valid".split()


def base_metadata(root) -> dict[str, str]:
    """Read the producer's key=value objects, never infer DATA/SIM from rows."""
    meta = dict(line.split("=", 1) for line in str(root["metadata"]).splitlines() if "=" in line)
    complete = dict(line.split("=", 1) for line in str(root["completion"]).splitlines() if "=" in line)
    if complete.get("completion_status") != "complete" or meta.get("contract") != "PhotonJetTrees":
        raise ValueError("base file lacks a complete producer receipt")
    if root["Events"].num_entries != int(complete["retained_events"]):
        raise ValueError("event rows disagree with completion accounting")
    if root["Sources"].num_entries != 1 or not bool(root["Sources"]["completed"].array(library="np")[0]):
        raise ValueError("base source is incomplete or has unsupported source boundaries")
    return meta


def _candidate_key(event, photon) -> str:
    return "".join(f"{int(v):016x}" for v in (*event, *photon))


@contextmanager
def score_index(sidecar: Path | None, base_hash: str, model: str, bindings: Mapping[str, Any]):
    """Disk-backed identity join: arbitrary sidecar row order, bounded RAM."""
    if sidecar is None:
        yield None, None
        return
    with tempfile.TemporaryDirectory(prefix="photonjet-scores-") as tmp, \
            closing(sqlite3.connect(str(Path(tmp) / "scores.sqlite"))) as db, uproot.open(sidecar) as root:
        receipt = json.loads(str(root["metadata"]))
        if receipt.get("schema") != "PhotonScoresSidecarV1" or receipt.get("input", {}).get("sha256") != base_hash:
            raise ValueError(f"{sidecar}: sidecar is not bound to this exact base file")
        if receipt.get("model") != model:
            raise ValueError(f"{sidecar}: wrong model")
        for field in ("model_sha256", "feature_definition_sha256"):
            if bindings.get(field) is not None and receipt.get(field) != bindings[field]:
                raise ValueError(f"{sidecar}: {field} differs from measurement binding")
        db.execute("CREATE TABLE scores (candidate TEXT PRIMARY KEY, source TEXT, score REAL, state INTEGER, used INTEGER DEFAULT 0)")
        n = 0
        for table in root["PhotonScores"].iterate(
                ["source_hi", "source_lo", "event_hi", "event_lo", "photon_hi", "photon_lo", "score", "state", "model"],
                step_size=8192, library="np"):
            rows = []
            for e, p, source, value, state, name in zip(_pairs(table, "event"), _pairs(table, "photon"),
                    _pairs(table, "source"), table["score"], table["state"], table["model"]):
                state = int(state)
                if str(name) != model or state not in (1, 2, 3, 4, 5) or (state == 1 and not math.isfinite(float(value))):
                    raise ValueError(f"{sidecar}: malformed score row")
                rows.append((_candidate_key(e, p), _candidate_key(source, ()),
                             float(value) if state == 1 else None, state))
            db.executemany("INSERT INTO scores(candidate,source,score,state) VALUES (?,?,?,?)", rows)
            n += len(rows)
        if n != int(receipt["candidates"]):
            raise ValueError("score receipt count differs from rows")
        db.commit()
        yield db, receipt
        if db.execute("SELECT COUNT(*) FROM scores WHERE used=0").fetchone()[0]:
            raise ValueError("sidecar contains candidates absent from the base file")


def _object_batches(tree, columns: list[str], event_batches):
    """Read contiguous producer event groups once; reject orphan/reordered rows."""
    chunks = iter(tree.iterate(columns, step_size=8192, library="np"))
    chunk = next(chunks, None)
    offset = 0
    for events in event_batches:
        keys = set(events)
        pieces = {c: [] for c in columns}
        count = 0
        while chunk is not None:
            end = offset
            while end < len(chunk["event_hi"]) and (int(chunk["event_hi"][end]), int(chunk["event_lo"][end])) in keys:
                end += 1
            if end == offset:
                break
            count += end - offset
            if count > 1000000:
                raise ValueError("object batch exceeds one million rows; reduce event_chunk_size")
            for column in columns:
                pieces[column].append(chunk[column][offset:end])
            offset = end
            if offset == len(chunk["event_hi"]):
                chunk = next(chunks, None)
                offset = 0
        yield {c: np.concatenate(v) if v else tree[c].array(entry_start=0, entry_stop=0, library="np")
               for c, v in pieces.items()}
    if chunk is not None:
        raise ValueError(f"{tree.name}: orphan or noncontiguous event rows")


def read_file(path: Path, *, system: str, config: dict[str, Any], model: str, scores_dir: Path | None):
    """Yield bounded event batches. Only compact event identities span a file."""
    file_hash = sha256(path)
    sidecar = scores_dir / f"{path.stem}.{model}.scores.root" if scores_dir else None
    sidecar_hash = sha256(sidecar) if sidecar else None
    block = config["systems"][system]
    view = JET_VIEW[config["recoil"]["jet_view"][system]]
    radius = float(config["isolation"]["cone_radius"])
    method = ISOLATION_METHOD[config["isolation"]["method"][system]]
    chunk_size = int(config.get("event_chunk_size", 256))
    if not 1 <= chunk_size <= 4096:
        raise ValueError("event_chunk_size must be in [1, 4096]")
    with uproot.open(path) as root, score_index(sidecar, file_hash, model, block) as (scores, score_receipt):
        metadata = base_metadata(root)
        if metadata.get("collision_system") != ("1" if system == "pp" else "2") or metadata.get("data_kind") not in ("1", "2"):
            raise ValueError(f"{path}: sample system/kind metadata mismatch")
        simulation = metadata["data_kind"] == "2"
        identities = root["Events"].arrays(["event_hi", "event_lo"], library="np")
        event_keys = _pairs(identities, "event")
        if len(set(event_keys)) != len(event_keys):
            raise ValueError(f"{path}: duplicate event identity")
        # Empty complete files still yield an empty batch with its schema/receipt.
        batches = [event_keys[i:i+chunk_size] for i in range(0, len(event_keys), chunk_size)] or [[]]
        tables = {name: _object_batches(root[name], columns.split(), batches)
                  for name, columns in TABLE_COLUMNS.items() if simulation or not name.startswith(("Truth", "PhotonTruth", "JetTruth"))}
        for batch_number, keys in enumerate(batches):
            start = batch_number * chunk_size
            ev = root["Events"].arrays(EVENT_COLUMNS, entry_start=start, entry_stop=start+len(keys), library="np")
            data = {name: next(reader) for name, reader in tables.items()}
            for table_name, count_name in (("Photons", "photon_count"), ("Jets", "jet_count"),
                                           ("TruthPhotons", "truth_photon_count"), ("TruthJets", "truth_jet_count")):
                if table_name not in data:
                    continue
                observed = Counter(_pairs(data[table_name], "event"))
                if any(observed.get(key, 0) != int(count) for key, count in zip(keys, ev[count_name])):
                    raise ValueError(f"{table_name}: object rows disagree with event accounting")
            ph = data["Photons"]
            photon_keys = [(*e, *p) for e, p in zip(_pairs(ph, "event"), _pairs(ph, "photon"))]
            if len(set(photon_keys)) != len(photon_keys):
                raise ValueError("duplicate photon identity")
            index = {key: i for i, key in enumerate(photon_keys)}
            sources = dict(zip(keys, _pairs(ev, "source")))
            score = np.full(len(photon_keys), np.nan)
            state = np.zeros(len(photon_keys), dtype=np.int32)
            if scores is not None:
                for i, key in enumerate(photon_keys):
                    encoded = _candidate_key(key[:2], key[2:])
                    row = scores.execute("SELECT source,score,state,used FROM scores WHERE candidate=?", (encoded,)).fetchone()
                    if row is None or row[3] or row[0] != _candidate_key(sources[key[:2]], ()):
                        raise ValueError("missing, duplicate, or wrong-source score identity")
                    score[i] = row[1] if row[2] == 1 else np.nan
                    state[i] = row[2]
                    scores.execute("UPDATE scores SET used=1 WHERE candidate=?", (encoded,))
                scores.commit()
            iso = np.full(len(photon_keys), np.nan)
            seen = set()
            table = data["Isolation"]
            for e, p, r, m, value, valid in zip(_pairs(table, "event"), _pairs(table, "photon"), table["radius"],
                                               table["method"], table["cone_sum"], table["valid"]):
                if int(m) != method or abs(float(r) - radius) > 1e-6:
                    continue
                key = (*e, *p)
                if key not in index or key in seen:
                    raise ValueError("unknown or duplicate photon isolation identity")
                seen.add(key)
                iso[index[key]] = float(value) if int(valid) == 1 else np.nan
            yield {"path": path, "file_hash": file_hash, "sidecar": sidecar, "sidecar_hash": sidecar_hash,
                   "score_receipt": score_receipt, "simulation": simulation, "metadata": metadata, "entry_start": start,
                   "events": ev, "photons": ph, "photon_keys": photon_keys, "score": score, "score_state": state,
                   "isolation": iso, "jets": data["Jets"], "truth_photons": data.get("TruthPhotons"),
                   "truth_jets": data.get("TruthJets"), "photon_links": data.get("PhotonTruthLinks"), "jet_links": data.get("JetTruthLinks"),
                   "view": view, "jet_radius": float(config["recoil"]["jet_radius"])}
        for reader in tables.values():
            if next(reader, None) is not None:
                raise ValueError("unconsumed object batch")
    if sha256(path) != file_hash or (sidecar and sha256(sidecar) != sidecar_hash):
        raise ValueError("input changed during histogram construction")


# ============================================================================
# Truth-signal predicate (re-evaluated from primitives)
# ============================================================================

def truth_signal_mask(truth: Mapping[str, np.ndarray], contract: Mapping[str, Any]) -> np.ndarray:
    classes = np.asarray(contract["valid_prompt_classes"], dtype=np.int64)
    radius = float(contract["isolation_cone_radius"])
    if abs(radius - 0.3) < 1e-9:
        iso = np.asarray(truth["isolation_r03"], dtype=float)
    elif abs(radius - 0.4) < 1e-9:
        iso = np.asarray(truth["isolation_r04"], dtype=float)
    else:
        raise ValueError("truth isolation cone must be 0.3 or 0.4 (the trees store both)")
    mask = np.asarray(truth["geant_valid"], dtype=bool)
    if contract.get("require_generator_association", False):
        mask &= np.asarray(truth["generator_association_valid"], dtype=bool)
    mask &= np.isin(np.asarray(truth["prompt_class"], dtype=np.int64), classes)
    mask &= np.asarray(truth["isolation_valid"], dtype=bool) & np.isfinite(iso) & (iso < float(contract["isolation_max_gev"]))
    eta = np.asarray(truth["eta"], dtype=float)
    mask &= np.isfinite(eta) & (np.abs(eta) < float(contract["abs_eta_max"]))
    return mask


def leading_truth_jet_pt(truth_jets: Mapping[str, np.ndarray], radius: float = 0.4) -> dict[EventKey, float]:
    best: dict[EventKey, float] = {}
    for e, r, pt in zip(_pairs(truth_jets, "event"), truth_jets["radius"], truth_jets["pt"]):
        if abs(float(r) - radius) > 1e-9 or not math.isfinite(float(pt)) or float(pt) <= 0:
            continue
        if float(pt) > best.get(e, -math.inf):
            best[e] = float(pt)
    return best


# ============================================================================
# The build
# ============================================================================

def event_multiplicity(replicas: int, seed: int, key: tuple[int, ...]) -> np.ndarray:
    """Stable event draws independent of file order, chunk size and model sidecar order."""
    digest = hashlib.sha256((str(seed) + ":" + ":".join(str(v) for v in key)).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:16], "big")).poisson(1.0, size=replicas)


def build(entries: list[tuple[Path, str | None]], *, system: str, config: dict[str, Any], manifest: SampleManifest | None,
          scores_dir: Path | None, truth_signal_only: bool, allow_missing_working_points: bool) -> Package:
    block = config["systems"][system]
    diagnostic = allow_missing_working_points
    approved = config.get("contract", {}).get("status") == "approved" and config.get("contract", {}).get("evidence")
    if not approved and not diagnostic:
        raise ValueError("measurement contract is unresolved; only explicit diagnostic mode is permitted")
    if not diagnostic and (scores_dir is None or any(not block.get(k) for k in ("model_sha256", "feature_definition_sha256"))):
        raise ValueError("nominal histogram construction requires exact model/feature bindings and score sidecars")
    model = str(block["model"])
    grids = Grids.from_config(config)
    n_pt, n_xj, n_reco, n_truth = grids.n_pt, grids.n_xj, grids.n_reco, grids.n_truth
    for edges in (grids.reco_pt_edges, grids.truth_pt_edges):
        if any(not np.any(np.isclose(edges, edge, rtol=0, atol=1e-9)) for edge in grids.pt_edges):
            raise ValueError("measurement boundaries must exist in fine response support")
    event_cfg = block["event"]
    photon_cfg = config["photon"]; recoil_cfg = config["recoil"]
    et_min, et_max, abs_eta_max = float(photon_cfg["et_min_gev"]), float(photon_cfg["et_max_gev"]), float(photon_cfg["abs_eta_max"])
    jet_pt_min, jet_abs_eta_max = float(recoil_cfg["jet_pt_min_gev"]), float(recoil_cfg["jet_abs_eta_max"])
    dphi_min = float(recoil_cfg["delta_phi_min_over_pi"]) * math.pi
    non_tight = str(config["abcd"]["non_tight_definition"])
    replicas = int(config["bootstrap"]["replicas"])

    missing = missing_working_points(block)
    if missing and not allow_missing_working_points:
        raise SystemExit(f"working points for {system} are not derived: {missing}. Complete config/measurement.yaml "
                         "or pass --allow-missing-working-points for a mechanism run (regions will be empty).")

    package: Package | None = None
    partitions: dict[str, Package] = {}
    partition_identity = {}
    seen_events = set()
    seen_file_hashes = set()
    family: str | None = None
    inputs: list[dict[str, Any]] = []
    counts = defaultdict(int)
    truth_flag_disagreements = 0

    for path, sample in entries:
        for data in read_file(path, system=system, config=config, model=model, scores_dir=scores_dir):
            if data["entry_start"] == 0:
                if data["file_hash"] in seen_file_hashes:
                    raise ValueError("duplicate base product content")
                seen_file_hashes.add(data["file_hash"])
            simulation = data["simulation"]
            if package is not None and ("pair_response" in package.arrays) != simulation:
                raise SystemExit(f"{path}: cannot mix data and simulation in one package")
            identity = {"sample": sample or "data", "period": data["metadata"].get("period", ""),
                        "si_di_role": data["metadata"].get("si_di_role", "")}
            partition_name = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
            partition_identity[partition_name] = identity
            if "::" in partition_name:
                raise ValueError("sample names cannot contain ::")
            if partition_name not in partitions:
                partitions[partition_name] = Package.empty(grids, simulation=simulation, replicas=replicas)
            package = partitions[partition_name]
            memory = sum(v.nbytes for part in partitions.values() for v in part.arrays.values())
            if 2 * memory > 512 * 1024**2:
                raise ValueError("histogram/bootstrap partitions exceed 512 MiB; reduce fanout or partition the build")
            ev = data["events"]

            # ---- weights ---------------------------------------------------------------
            event_keys = _pairs(ev, "event")
            if seen_events.intersection(event_keys):
                raise ValueError("overlapping event identities across base products")
            seen_events.update(event_keys)
            if len(set(event_keys)) != len(event_keys):
                raise ValueError(f"{path}: event identities are not unique")
            if sample is None:
                if simulation:
                    raise SystemExit(f"{path}: simulation input needs a sample name for its weights")
                weights = np.asarray(ev["event_weight"], dtype=float)
                if np.any(np.asarray(ev["event_weight_valid"]) != 1) or np.any(~np.isfinite(weights)) or np.any(weights < 0):
                    raise ValueError("data producer weights are incomplete or invalid")
                receipt = {"sample": None, "family": "data", "components": ["producer_event_weight"], "dropped": {}}
                this_family = "data"
            else:
                if manifest is None:
                    raise SystemExit("a sample manifest is required when a sample name is given")
                spec = manifest.samples[sample]
                if spec["system"] != system or (spec["kind"] == "simulation") != simulation:
                    raise ValueError("sample manifest and base system/kind disagree")
                if spec["source_manifest_sha256"] != data["metadata"].get("source_manifest_sha256"):
                    raise ValueError("sample is not bound to this base source manifest")
                leading = None
                if manifest.samples[sample]["source_factor"] == "ownership_stitch":
                    by_event = leading_truth_jet_pt(data["truth_jets"])
                    leading = np.asarray([by_event.get(k, math.nan) for k in event_keys], dtype=float)
                result = complete_event_weights(manifest, sample, event_weight=ev["event_weight"], event_weight_valid=ev["event_weight_valid"],
                                                centrality=np.asarray(ev["centrality_percent"], dtype=float), leading_truth_jet_pt=leading)
                weights = result["weight"]
                receipt = {"sample": sample, "family": result["family"], "components": result["components"], "dropped": result["dropped"]}
                this_family = result["family"]
            if family is None:
                family = this_family
            elif family != this_family:
                raise SystemExit(f"{path}: sample family {this_family!r} differs from the package family {family!r}")

            # ---- event selection --------------------------------------------------------
            accepted = np.isfinite(weights)
            if not simulation:
                trigger = event_cfg.get("trigger_selection")
                if not trigger or not trigger.get("evidence"):
                    if not diagnostic:
                        raise ValueError("data trigger selection/inclusive ownership must be explicitly bound")
                elif trigger.get("mode") == "bits_any":
                    word = trigger.get("word")
                    if word not in ("live", "scaled"):
                        raise ValueError("trigger word must be live or scaled; legacy word is not a raw input word")
                    bits = trigger.get("bits", [])
                    if not bits or any(type(b) is not int or not 0 <= b < 64 for b in bits):
                        raise ValueError("invalid trigger bit list")
                    available_bit = 1 if word == "live" else 2
                    known = (np.asarray(ev["trigger_packet_valid"]) == 1) & (np.asarray(ev["trigger_packet_status"]) == 0)
                    known &= (np.asarray(ev["trigger_decisions_available"]) & available_bit) != 0
                    if np.any(accepted & ~known):
                        raise ValueError("selected trigger word is unavailable or invalid")
                    mask = np.uint64(sum(1 << b for b in set(bits)))
                    accepted &= (np.asarray(ev[f"trigger_{word}_bits"], dtype=np.uint64) & mask) != 0
                elif trigger.get("mode") != "inclusive":
                    raise ValueError("unsupported trigger selection mode")
            if event_cfg.get("terminal_status") is not None:
                accepted &= np.asarray(ev["terminal_status"]) == int(event_cfg["terminal_status"])
            vz = np.asarray(ev["reco_vertex_z"], dtype=float)
            accepted &= np.asarray(ev["reco_vertex_valid"], dtype=bool) & np.isfinite(vz) & (np.abs(vz) < float(event_cfg["abs_vertex_z_max_cm"]))
            if event_cfg.get("require_object_vertex_in_domain", True):
                accepted &= np.asarray(ev["reco_object_vertex_in_domain"], dtype=bool)
            cent = np.asarray(ev["centrality_percent"], dtype=float)
            if event_cfg.get("centrality_percent") is not None:
                low, high = (float(v) for v in event_cfg["centrality_percent"])
                accepted &= np.asarray(ev["centrality_valid"], dtype=bool) & (cent >= low) & (cent < high)
            if event_cfg.get("require_minimum_bias", False):
                accepted &= np.asarray(ev["minimum_bias_decision"]) == TRISTATE_TRUE
            event_row = {k: i for i, k in enumerate(event_keys)}

            # ---- photons ------------------------------------------------------------------
            ph = data["photons"]; keys = data["photon_keys"]
            et = np.asarray(ph["et"], dtype=float); eta = np.asarray(ph["eta"], dtype=float)
            ordinal = np.asarray(ph["encounter_ordinal"], dtype=np.int64)
            photon_event = np.asarray([event_row[k[:2]] for k in keys], dtype=np.int64)
            photon_cent = cent[photon_event]
            thresholds = working_point_thresholds(block, et, photon_cent)
            classes = classify(data["score"], data["isolation"], thresholds, non_tight)
            kinematic = np.isfinite(et) & np.isfinite(eta) & (et >= et_min) & (et < et_max) & (np.abs(eta) < abs_eta_max)
            kinematic &= accepted[photon_event]

            if np.any(accepted & (np.asarray(ev["reco_photon_capture_state"]) != CAPTURE_COMPLETE)):
                raise ValueError("selected event has incomplete photon capture")
            if simulation and np.any(accepted & (np.asarray(ev["truth_denominator_complete"]) != 1)):
                raise ValueError("selected simulation event has an incomplete truth census")
            for ei in np.flatnonzero(accepted):
                views = list(zip(ev["reco_jet_view"][ei], ev["reco_jet_radius_code"][ei], ev["reco_jet_capture_state"][ei]))
                states = [int(st) for v, r, st in views if int(v) == data["view"] and abs(float(r)/100.0-data["jet_radius"]) < 1e-6]
                if states != [CAPTURE_COMPLETE]:
                    raise ValueError("selected reconstructed jet view is missing or incomplete")
                if simulation:
                    states = [int(st) for r, st in zip(ev["truth_jet_radius_code"][ei], ev["truth_jet_container_valid"][ei])
                              if abs(float(r)/100.0-data["jet_radius"]) < 1e-6]
                    if states != [CAPTURE_COMPLETE]:
                        raise ValueError("selected truth jet collection is missing or incomplete")

            support = kinematic.copy()
            if simulation:
                support = accepted[photon_event] & np.isfinite(et) & np.isfinite(eta) & (np.abs(eta) < abs_eta_max)
                support &= (et >= grids.reco_pt_edges[0]) & (et < grids.reco_pt_edges[-1])
            if not diagnostic and np.any(support & (~classes["id_valid"] | ~classes["iso_valid"])):
                raise ValueError("classification support lacks finite scores/isolation/working points; resolve model application domain")

            # truth context
            signal_keys: set[PhotonKey] = set()
            truth_by_event: dict[EventKey, list[int]] = defaultdict(list)
            photon_to_truth: dict[PhotonKey, tuple[int, int]] = {}
            jet_to_truth: dict[tuple[int, int, int, int], tuple[int, int]] = {}
            truth_jets_by_event: dict[EventKey, list[int]] = defaultdict(list)
            if simulation:
                tp = data["truth_photons"]
                signal = truth_signal_mask(tp, config["truth_signal"])
                truth_flag_disagreements += int(np.sum(signal != np.asarray(tp["analysis_signal"], dtype=bool)))
                truth_keys = [(*e, *t) for e, t in zip(_pairs(tp, "event"), _pairs(tp, "truth_photon"))]
                signal_truth = {k for k, ok in zip(truth_keys, signal) if ok}
                for i, k in enumerate(truth_keys):
                    if signal[i]:
                        truth_by_event[k[:2]].append(i)
                pl = data["photon_links"]
                for e, p, t, state in zip(_pairs(pl, "event"), _pairs(pl, "photon"), _pairs(pl, "truth_photon"), pl["state"]):
                    if int(state) == ASSOCIATION_MATCHED and (*e, *t) in signal_truth:
                        photon_to_truth[(*e, *p)] = t
                signal_keys = set(photon_to_truth)
                jl = data["jet_links"]
                for e, j, t, jv, jr, state, sel in zip(_pairs(jl, "event"), _pairs(jl, "jet"), _pairs(jl, "truth_jet"), jl["jet_view"],
                                                        jl["jet_radius"], jl["state"], jl["selected_match"]):
                    if int(jv) == data["view"] and abs(float(jr) - data["jet_radius"]) < 1e-9 and int(state) == ASSOCIATION_MATCHED and bool(sel):
                        jet_to_truth[(*e, *j)] = t
                tj = data["truth_jets"]
                for i, (e, r) in enumerate(zip(_pairs(tj, "event"), tj["radius"])):
                    if abs(float(r) - data["jet_radius"]) < 1e-9:
                        truth_jets_by_event[e].append(i)
            if simulation:
                unknown_photons = {(*e, *p) for e, p, st in zip(_pairs(pl, "event"), _pairs(pl, "photon"), pl["state"])
                                   if int(st) not in (3, 4, 5)}
                if any(k in unknown_photons for k, keep in zip(keys, support) if keep):
                    raise ValueError("selected photon has unknown truth association")
                known_photons = {(*e, *p) for e, p, st in zip(_pairs(pl, "event"), _pairs(pl, "photon"), pl["state"])
                                 if int(st) in (3, 4)}
                if any(k not in known_photons for k, keep in zip(keys, support) if keep):
                    raise ValueError("classification support lacks complete photon-truth relations")
                known_jets = {(*e, *j) for e, j, v, r, st, sel in zip(_pairs(jl, "event"), _pairs(jl, "jet"),
                              jl["jet_view"], jl["jet_radius"], jl["state"], jl["selected_match"])
                              if int(v) == data["view"] and abs(float(r)-data["jet_radius"]) < 1e-6
                              and (int(st) == 4 or (int(st) == 3 and int(sel) == 1))}
                jt = data["jets"]
                if any((*e, *j) not in known_jets for e, j, v, r in zip(_pairs(jt, "event"), _pairs(jt, "jet"), jt["view"], jt["radius"])
                       if accepted[event_row[e]] and int(v) == data["view"] and abs(float(r)-data["jet_radius"]) < 1e-6):
                    raise ValueError("selected jet collection lacks complete truth relations")
                if any(accepted[event_row[e]] and int(st) not in (2, 3, 4, 5)
                       for e, v, r, st in zip(_pairs(jl, "event"), jl["jet_view"], jl["jet_radius"], jl["state"])
                       if int(v) == data["view"] and abs(float(r)-data["jet_radius"]) < 1e-6):
                    raise ValueError("selected jet view has unknown truth association")
            if truth_signal_only and not simulation:
                raise ValueError("truth-signal-only requires simulation")
            if truth_signal_only:
                kinematic &= np.asarray([k in signal_keys for k in keys], dtype=bool)

            by_event: dict[EventKey, list[int]] = defaultdict(list)
            for i, k in enumerate(keys):
                by_event[k[:2]].append(i)

            jets = data["jets"]
            jets_by_event = defaultdict(list)
            seen_jets = set()
            for j, (e, key) in enumerate(zip(_pairs(jets, "event"), _pairs(jets, "jet"))):
                if (*e, *key) in seen_jets:
                    raise ValueError("duplicate jet identity")
                seen_jets.add((*e, *key))
                if int(jets["view"][j]) == data["view"] and abs(float(jets["radius"][j]) - data["jet_radius"]) < 1e-6:
                    jets_by_event[e].append(j)

            def accepted_recoil(photon_index: int) -> list[tuple[int, int, float, float]]:
                """Derive only the selected photon's recoil pairs, without a Cartesian table."""
                out = []
                if not (math.isfinite(et[photon_index]) and et[photon_index] > 0):
                    return out
                for j in jets_by_event.get(keys[photon_index][:2], []):
                    pt, eta_j, phi = (float(jets[k][j]) for k in ("corrected_pt", "eta", "phi"))
                    if not all(math.isfinite(v) for v in (pt, eta_j, phi)) or int(jets["calibration_valid"][j]) != 1:
                        raise ValueError("selected jet view has invalid kinematics/calibration")
                    dphi = abs(math.atan2(math.sin(phi-float(ph["phi"][photon_index])), math.cos(phi-float(ph["phi"][photon_index]))))
                    if pt > jet_pt_min and abs(eta_j) < jet_abs_eta_max and dphi >= dphi_min:
                        out.append((j, j, pt, pt / et[photon_index]))
                return out

            # ---- per-event fills ----------------------------------------------------------
            for ekey in event_keys:
                rows = by_event.get(ekey, [])
                ei = event_row[ekey]
                if not accepted[ei]:
                    continue
                w = float(weights[ei]); w_rep = event_multiplicity(replicas, int(config["bootstrap"]["seed"]),
                                                   (int(ev["source_hi"][ei]), int(ev["source_lo"][ei]), *ekey))
                rows_a = np.asarray(rows, dtype=np.int64)
                counts["events_accepted"] += 1
                counts["events_with_candidates"] += bool(rows)
                package.fill("accepted_event_weight", 0, w)

                # kinematic denominator (no identification)
                for i in rows_a[kinematic[rows_a]]:
                    b = axis_bin(et[i], grids.pt_edges)
                    if b is not None:
                        package.fill("photon_et_denominator", b, w)

                for r_index, region in enumerate(REGION_ORDER):
                    eligible = kinematic[rows_a] & classes[region][rows_a]
                    local = event_leading_index(et[rows_a], ordinal[rows_a], eligible)
                    if local < 0:
                        continue
                    leader = int(rows_a[local])
                    b = axis_bin(et[leader], grids.pt_edges)
                    if b is None:
                        continue
                    name = "leakage_abcd_counts" if truth_signal_only else "abcd_counts"
                    package.fill(name, (r_index, b), w)
                    if not truth_signal_only:
                        package.arrays["abcd_events"][r_index, b] += 1
                        package.arrays["bootstrap_abcd_counts"][:, r_index, b] += w * w_rep
                    if region in ("A", "C"):
                        s_index = 0 if region == "A" else 1
                        for _, _, _, xj in accepted_recoil(leader):
                            xb = axis_bin(xj, grids.xj_edges)
                            if xb is None:
                                continue
                            package.fill("recoil_spectra", (s_index, b, xb), w)
                            package.arrays["bootstrap_recoil_spectra"][:, s_index, b, xb] += w * w_rep
                            counts["recoil_pairs"] += 1

                if not simulation:
                    continue

                # ---- responses ----------------------------------------------------------
                truth_rows = truth_by_event.get(ekey, [])
                if len(truth_rows) > 1:
                    raise ValueError(f"{path}: event {ekey} has {len(truth_rows)} truth-signal photons; the response needs one")
                tp = data["truth_photons"]
                truth_index = truth_rows[0] if truth_rows else None
                truth_pt = float(tp["pt"][truth_index]) if truth_index is not None else None
                truth_bin = axis_bin(truth_pt, grids.truth_pt_edges) if truth_index is not None else None
                truth_key = (int(tp["truth_photon_hi"][truth_index]), int(tp["truth_photon_lo"][truth_index])) if truth_index is not None else None

                # region-A leader over the classification grid (wider than the window)
                wide = np.isfinite(et[rows_a]) & np.isfinite(eta[rows_a]) & (et[rows_a] > 0) \
                    & (np.abs(eta[rows_a]) < abs_eta_max) & classes["A"][rows_a]
                if truth_signal_only:
                    wide &= np.asarray([keys[i] in signal_keys for i in rows_a], dtype=bool)
                local = event_leading_index(et[rows_a], ordinal[rows_a], wide)
                wide_leader = int(rows_a[local]) if local >= 0 else -1
                reco_bin = axis_bin(et[wide_leader], grids.reco_pt_edges) if wide_leader >= 0 else None
                if wide_leader >= 0 and keys[wide_leader] in signal_keys:
                    denominator_bin = axis_bin(et[wide_leader], grids.pt_edges)
                    if denominator_bin is not None:
                        package.fill("combinatoric_photon_denominator", denominator_bin, w)
                matched = wide_leader >= 0 and truth_key is not None and photon_to_truth.get(keys[wide_leader]) == truth_key

                if truth_bin is not None:
                    package.fill("truth_photons", truth_bin, w)
                if reco_bin is not None:
                    package.fill("photon_reco", reco_bin, w)
                if matched:
                    if truth_bin is not None and reco_bin is not None:
                        package.fill("photon_response", (truth_bin, reco_bin), w)
                    elif truth_bin is not None:
                        package.fill("photon_misses", truth_bin, w)
                        package.component("photon_boundary_misses", (n_truth,))
                        package.fill("photon_boundary_misses", truth_bin, w)
                    elif reco_bin is not None:
                        package.fill("photon_boundary_fakes", reco_bin, w)
                else:
                    if truth_bin is not None:
                        package.fill("photon_misses", truth_bin, w)
                    if reco_bin is not None:
                        package.fill("photon_fakes", reco_bin, w)

                # truth pairs: truth-signal photon x accepted truth jets
                truth_pairs: dict[tuple[int, int], tuple[int | None, float]] = {}
                if truth_index is not None:
                    tj = data["truth_jets"]; tphi = float(tp["phi"][truth_index])
                    for j in truth_jets_by_event.get(ekey, []):
                        pt_j, eta_j, phi_j = float(tj["pt"][j]), float(tj["eta"][j]), float(tj["phi"][j])
                        if not (math.isfinite(pt_j) and math.isfinite(eta_j) and math.isfinite(phi_j)):
                            continue
                        dphi = abs(math.atan2(math.sin(phi_j - tphi), math.cos(phi_j - tphi)))
                        if not (pt_j > jet_pt_min and abs(eta_j) < jet_abs_eta_max and dphi >= dphi_min):
                            continue
                        xj_t = pt_j / truth_pt
                        xb = axis_bin(xj_t, grids.xj_edges)
                        g = grids.flatten(truth_bin, xb) if (truth_bin is not None and xb is not None) else None
                        truth_pairs[(int(tj["truth_jet_hi"][j]), int(tj["truth_jet_lo"][j]))] = (g, xj_t)
                        if g is not None:
                            package.fill("pair_truth", g, w)

                matched_truth_jets: set[tuple[int, int]] = set()
                if wide_leader >= 0:
                    leader_key = keys[wide_leader]
                    linked_photon = photon_to_truth.get(leader_key)
                    photon_ok = matched  # the leader is the truth-signal photon
                    for row, j, _, xj in accepted_recoil(wide_leader):
                        xb = axis_bin(xj, grids.xj_edges)
                        g_reco = grids.flatten(reco_bin, xb) if reco_bin is not None and xb is not None else None
                        if g_reco is not None:
                            package.fill("pair_reco", g_reco, w)
                        jkey = (ekey[0], ekey[1], int(jets["jet_hi"][j]), int(jets["jet_lo"][j]))
                        linked_jet = jet_to_truth.get(jkey)
                        truth_pair = truth_pairs.get(linked_jet) if (photon_ok and linked_jet is not None) else None
                        if truth_pair is not None:
                            matched_truth_jets.add(linked_jet)
                            g_truth, _ = truth_pair
                            if g_truth is not None and g_reco is not None:
                                package.fill("pair_response", (g_truth, g_reco), w)
                            elif g_truth is not None:
                                package.component("pair_boundary_truth_RECO_OFF_GRID", (n_truth * n_xj,))
                                package.fill("pair_boundary_truth_RECO_OFF_GRID", g_truth, w)
                            elif g_reco is not None:
                                package.component("pair_boundary_reco_TRUTH_OFF_GRID", (n_reco * n_xj,))
                                package.fill("pair_boundary_reco_TRUTH_OFF_GRID", g_reco, w)
                            continue
                        if g_reco is None:
                            continue
                        if linked_photon is not None and linked_jet is not None:
                            package.fill(f"pair_boundary_reco_{ASSIGNED_HARD_NONFIDUCIAL}", g_reco, w)
                        elif linked_photon is None:
                            package.fill(f"pair_fakes_{UNMATCHED_RECO}", g_reco, w)
                        else:
                            package.fill(f"pair_fakes_{COMBINATORIC}", g_reco, w)
                for tkey, (g_truth, _) in truth_pairs.items():
                    if g_truth is not None and tkey not in matched_truth_jets:
                        package.fill("pair_misses", g_truth, w)

            inputs.append({"path": str(path), "sha256": data["file_hash"], "sample": sample, "entry_start": data["entry_start"],
                           "scores": None if data["sidecar"] is None else {"path": str(data["sidecar"]), "sha256": data["sidecar_hash"]},
                           "events": int(len(event_keys)), "events_accepted": int(accepted.sum()), "photons": int(len(keys)),
                           "score_states": {int(s): int(n) for s, n in zip(*np.unique(data["score_state"], return_counts=True))},
                           "weights": receipt})
    if package is None:
        raise SystemExit("no inputs")
    if truth_flag_disagreements:
        raise ValueError("configured truth signal disagrees with producer reference flags")
    combined = {}
    for partition in partitions.values():
        partition.check()
        for name, values in partition.arrays.items():
            combined[name] = combined.get(name, np.zeros_like(values)) + values
    package = Package(grids, combined, sample_arrays={name: part.arrays for name, part in partitions.items()})
    package.metadata = {
        "status": "diagnostic" if diagnostic else "analysis_ready",
        "partition_identity": partition_identity,
        "implementation_sha256": {"build.py": sha256(Path(__file__)),
                                  "histogram_contract.py": sha256(HERE / "histogram_contract.py")},
        "resolved_configuration": config,
        "measurement_sha256": hashlib.sha256(json.dumps(config, sort_keys=True, allow_nan=False).encode()).hexdigest(),
        "sample_manifest_sha256": manifest.sha256 if manifest else None,
        "model_sha256": block.get("model_sha256"), "feature_definition_sha256": block.get("feature_definition_sha256"),
        "exposure_support": "producer event weights/counts only; luminosity requires external approved exposure contract",
        "family": family, "system": system, "measurement_version": config.get("measurement_version"),
        "model": model, "truth_signal_only": truth_signal_only,
        "selection": {"photon_et_gev": [et_min, et_max], "photon_abs_eta_max": abs_eta_max, "event": dict(event_cfg),
                      "isolation": {"cone_radius": float(config["isolation"]["cone_radius"]), "method": config["isolation"]["method"][system]},
                      "non_tight_definition": non_tight, "jet_view": recoil_cfg["jet_view"][system], "jet_radius": float(recoil_cfg["jet_radius"]),
                      "jet_pt_min_gev": jet_pt_min, "jet_abs_eta_max": jet_abs_eta_max, "delta_phi_min_rad": dphi_min, "delta_phi_inclusive": True,
                      "working_points": block.get("working_points"), "working_points_missing": missing,
                      "truth_signal": config["truth_signal"]},
        "bootstrap": {"replicas": replicas, "seed": int(config["bootstrap"]["seed"]), "rule": "source-event-keyed Poisson(1); common draws for ABCD and recoil",
                      "scope": "data ABCD/recoil covariance only; no joint simulation response/leakage covariance"},
        "counts": dict(counts), "truth_flag_disagreements": truth_flag_disagreements, "inputs": inputs,
    }
    return package


def read_inputs(inputs: list[Path], lists: list[Path], sample: str | None) -> list[tuple[Path, str | None]]:
    entries = [(Path(p), sample) for p in inputs]
    for listing in lists:
        for line in Path(listing).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            p = Path(parts[0])
            entries.append((p if p.is_absolute() else Path(listing).parent / p, parts[1] if len(parts) > 1 else sample))
    if not entries:
        raise SystemExit("no input files given (use --input or --input-list)")
    if len({p.resolve() for p, _ in entries}) != len(entries):
        raise ValueError("duplicate input path")
    return entries


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", action="append", default=[], type=Path)
    parser.add_argument("--input-list", action="append", default=[], type=Path, help="lines: path [sample]")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path, help="config/measurement.yaml")
    parser.add_argument("--samples", type=Path, help="config/samples.yaml (default: next to --config)")
    parser.add_argument("--sample", help="sample name applied to every --input")
    parser.add_argument("--scores", type=Path, help="directory of PhotonID/augment.py sidecars")
    parser.add_argument("--truth-signal-only", action="store_true", help="simulation: truth-signal candidates only (leakage package)")
    parser.add_argument("--diagnostic", "--allow-missing-working-points", dest="allow_missing_working_points", action="store_true",
                        help="write a diagnostic package that FinalAnalysis will refuse")
    parser.add_argument("--output-stem", required=True, type=Path, help="writes <stem>.npz and <stem>.json")
    args = parser.parse_args(argv)

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    manifest_path = args.samples or (args.config.parent / "samples.yaml")
    manifest = SampleManifest(manifest_path) if manifest_path.is_file() else None
    entries = read_inputs(args.input, args.input_list, args.sample)
    if any(s is not None for _, s in entries) and manifest is None:
        raise SystemExit(f"sample manifest not found: {manifest_path}")
    package = build(entries, system=args.system, config=config, manifest=manifest, scores_dir=args.scores,
                    truth_signal_only=args.truth_signal_only, allow_missing_working_points=args.allow_missing_working_points)
    npz, meta = package.save(args.output_stem)
    a = package.arrays
    print(f"[build] {npz}: family={package.metadata['family']} A/B/C/D={a['abcd_counts'].sum(axis=1).tolist()} "
          f"recoil_pairs={package.metadata['counts'].get('recoil_pairs', 0)} check={package.check()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
