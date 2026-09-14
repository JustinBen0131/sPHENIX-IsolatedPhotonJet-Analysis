"""Nominal photon selection shared by scoring, working points and histograms.

This is the single implementation of

* the truth-signal contract (which simulated photons count as signal, both for
  training and for the selected truth photon),
* working-point threshold evaluation from the configuration, and
* tight / non-tight / isolated / non-isolated flags with explicit boundaries.

Boundaries are strict: tight means ``score > tight_threshold``; the bounded
non-tight band is ``nontight_low < score <= nontight_high``; isolated means
``iso < isolated_max``; non-isolated means ``iso > nonisolated_min``.  A band
may leave a gap to the tight threshold and a candidate between the isolated
and non-isolated thresholds belongs to neither region.  Candidates with a
non-finite score, isolation or threshold are invalid, not background.

Units: energies in GeV, centrality in percent.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml


FEATURE_COUNT = {"pp": 11, "auau": 14}
REGION_ORDER = ("A", "B", "C", "D")

# Fields the truth-signal contract needs from the ``truthPhotons`` tree.
TRUTH_SIGNAL_FIELDS = (
    "prompt_class",
    "g4_photon_valid",
    "hepmc_association_valid",
    "truth_isolation_valid",
    "truth_isolation_r03",
    "truth_photon_eta",
)


def load_config(path: Path | str) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict) or "systems" not in config:
        raise ValueError(f"{path}: not a nominal configuration file")
    return config


def system_config(config: Mapping[str, Any], system: str) -> dict[str, Any]:
    if system not in FEATURE_COUNT:
        raise ValueError(f"collision system must be pp or auau, not {system!r}")
    block = config["systems"].get(system)
    if block is None:
        raise ValueError(f"configuration has no '{system}' section")
    features = list(block.get("features", []))
    if len(features) != FEATURE_COUNT[system]:
        raise ValueError(
            f"{system} feature list has {len(features)} entries; the stored trees carry "
            f"{FEATURE_COUNT[system]} ordered model inputs"
        )
    return block


# --------------------------------------------------------------------------
# Truth-signal contract
# --------------------------------------------------------------------------
def truth_signal_mask(truth: Mapping[str, np.ndarray], contract: Mapping[str, Any]) -> np.ndarray:
    """Boolean mask over ``truthPhotons`` rows satisfying the nominal signal contract."""

    missing = [name for name in TRUTH_SIGNAL_FIELDS if name not in truth]
    if missing:
        raise KeyError(
            "truthPhotons is missing the fields needed for the nominal signal contract: "
            f"{missing}. Trees produced before the R=0.3 truth-isolation capture cannot "
            "define the nominal training signal; regenerate them with TreeProduction."
        )
    n = len(truth["prompt_class"])
    mask = np.ones(n, dtype=bool)
    if contract.get("require_geant_photon", True):
        mask &= np.asarray(truth["g4_photon_valid"]) == 1
    if contract.get("require_generator_association", True):
        mask &= np.asarray(truth["hepmc_association_valid"]) == 1
    classes = np.asarray(contract["valid_prompt_classes"], dtype=np.int64)
    mask &= np.isin(np.asarray(truth["prompt_class"], dtype=np.int64), classes)
    iso = np.asarray(truth["truth_isolation_r03"], dtype=float)
    mask &= np.asarray(truth["truth_isolation_valid"]) == 1
    mask &= np.isfinite(iso) & (iso < float(contract["truth_isolation_max_gev"]))
    eta = np.asarray(truth["truth_photon_eta"], dtype=float)
    mask &= np.isfinite(eta) & (np.abs(eta) < float(contract["abs_eta_max"]))
    return mask


# --------------------------------------------------------------------------
# Thresholds
# --------------------------------------------------------------------------
def threshold_values(spec: Any, photon_et: np.ndarray, centrality: np.ndarray) -> np.ndarray:
    """Evaluate one configured threshold for every candidate (NaN when not derived)."""

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
        pt_bin = np.searchsorted(pt_edges, et, side="right") - 1
        inside = (pt_bin >= 0) & (pt_bin < len(pt_edges) - 1)
        out = np.full(et.shape, np.nan)
        if "cent_edges" in spec:
            cent_edges = np.asarray(spec["cent_edges"], dtype=float)
            if values.shape != (len(pt_edges) - 1, len(cent_edges) - 1):
                raise ValueError("binned threshold values do not match pt_edges x cent_edges")
            cent_bin = np.searchsorted(cent_edges, cent, side="right") - 1
            inside &= (cent_bin >= 0) & (cent_bin < len(cent_edges) - 1)
            out[inside] = values[pt_bin[inside], cent_bin[inside]]
        else:
            if values.shape != (len(pt_edges) - 1,):
                raise ValueError("binned threshold values do not match pt_edges")
            out[inside] = values[pt_bin[inside]]
        return out
    intercept = float(spec.get("intercept", 0.0))
    slope_et = float(spec.get("slope_photon_et", 0.0))
    slope_cent = float(spec.get("slope_centrality", 0.0))
    out = intercept + slope_et * et
    if slope_cent != 0.0:
        out = out + slope_cent * cent
    return out


def working_point_thresholds(
    block: Mapping[str, Any], photon_et: np.ndarray, centrality: np.ndarray
) -> dict[str, np.ndarray]:
    """Evaluate the five working-point thresholds of one system for every candidate."""

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
    names = ("tight_threshold", "nontight_low", "nontight_high", "isolated_max_gev", "nonisolated_min_gev")
    return [name for name in names if wp.get(name) is None]


# --------------------------------------------------------------------------
# Flags
# --------------------------------------------------------------------------
def classify(
    score: np.ndarray,
    isolation: np.ndarray,
    thresholds: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Return int32 flag arrays plus validity masks.

    ``id_valid`` requires a finite score and finite tight/non-tight thresholds;
    ``iso_valid`` requires a finite isolation and finite isolation thresholds.
    Invalid candidates get every flag set to 0.
    """

    s = np.asarray(score, dtype=float)
    iso = np.asarray(isolation, dtype=float)
    t = np.asarray(thresholds["tight"], dtype=float)
    lo = np.asarray(thresholds["nontight_low"], dtype=float)
    hi = np.asarray(thresholds["nontight_high"], dtype=float)
    iso_max = np.asarray(thresholds["isolated_max"], dtype=float)
    noniso_min = np.asarray(thresholds["nonisolated_min"], dtype=float)

    id_valid = np.isfinite(s) & np.isfinite(t) & np.isfinite(lo) & np.isfinite(hi)
    iso_valid = np.isfinite(iso) & np.isfinite(iso_max) & np.isfinite(noniso_min)
    with np.errstate(invalid="ignore"):
        tight = id_valid & (s > t)
        nontight = id_valid & (s > lo) & (s <= hi)
        not_tight = id_valid & ~(s > t)
        isolated = iso_valid & (iso < iso_max)
        nonisolated = iso_valid & (iso > noniso_min)
    return {
        "bdt_is_tight": tight.astype(np.int32),
        "bdt_is_nontight": nontight.astype(np.int32),
        "bdt_is_not_tight": not_tight.astype(np.int32),
        "iso_pass": isolated.astype(np.int32),
        "iso_nonisolated": nonisolated.astype(np.int32),
        "id_valid": id_valid,
        "iso_valid": iso_valid,
    }


def region_masks(flags: Mapping[str, np.ndarray], non_tight_definition: str) -> dict[str, np.ndarray]:
    """ABCD region membership from flags. A = tight & isolated, B = tight & non-isolated,
    C = non-tight & isolated, D = non-tight & non-isolated."""

    if non_tight_definition not in {"bounded", "complement"}:
        raise ValueError("non_tight_definition must be bounded or complement")
    tight = flags["bdt_is_tight"] == 1
    non_tight = (
        flags["bdt_is_nontight"] == 1
        if non_tight_definition == "bounded"
        else flags["bdt_is_not_tight"] == 1
    )
    isolated = flags["iso_pass"] == 1
    nonisolated = flags["iso_nonisolated"] == 1
    return {
        "A": tight & isolated,
        "B": tight & nonisolated,
        "C": non_tight & isolated,
        "D": non_tight & nonisolated,
    }


def event_leading_index(
    photon_et: np.ndarray, encounter_ordinal: np.ndarray, eligible: np.ndarray
) -> int:
    """Index of the highest-ET eligible photon; ties resolve to the lowest ordinal."""

    best = -1
    for index in np.flatnonzero(np.asarray(eligible, dtype=bool)):
        if best < 0:
            best = int(index)
            continue
        if photon_et[index] > photon_et[best] or (
            photon_et[index] == photon_et[best]
            and encounter_ordinal[index] < encounter_ordinal[best]
        ):
            best = int(index)
    return best


def delta_phi_min(config: Mapping[str, Any]) -> float:
    return float(config["recoil"]["delta_phi_min_over_pi"]) * math.pi


__all__ = [
    "FEATURE_COUNT",
    "REGION_ORDER",
    "TRUTH_SIGNAL_FIELDS",
    "classify",
    "delta_phi_min",
    "event_leading_index",
    "load_config",
    "missing_working_points",
    "region_masks",
    "system_config",
    "threshold_values",
    "truth_signal_mask",
    "working_point_thresholds",
]
