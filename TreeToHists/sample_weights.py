"""Complete per-event analysis weights from the sample manifest.

This is the single implementation of

    w = event_weight x source factor x centrality factor

for every sample in ``config/samples.yaml``.  The histogram stage and the
response build both call :func:`complete_event_weights`; nothing else
multiplies weights.  Sumw2 is always accumulated from the complete weight
squared by the caller.

Rules
  * The producer's stored ``event_weight`` is read, never assumed to be 1.
  * Embedded inclusive-jet slices own an event only when its leading anti-kT
    R=0.4 truth-jet pT falls in the slice's half-open window; other events of
    that slice are dropped, never reweighted.
  * The centrality factor is applied after the source factor, once, and only
    for Au+Au simulation with centrality inside the map support [0, 80).
    Events outside the support are dropped.  Non-finite centrality or
    non-positive map values are errors, never clamped.
  * Every returned weight carries the ordered list of applied components so a
    receipt can show what was multiplied.

Units: cross sections in pb, centrality in percent, pT in GeV.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml


SOURCE_FACTORS = ("family_reference", "cross_section_over_generated", "ownership_stitch")
FAMILIES = ("data", "photonjet", "inclusivejet")
TRUTH_JET_RADIUS = 0.4


@dataclass(frozen=True)
class Sample:
    name: str
    system: str
    kind: str
    family: str
    source_factor: str | None = None
    factor_value: float | None = None          # pb per (generated or owned) event, or None
    ownership_window: tuple[float, float | None] | None = None
    note: str = ""

    @property
    def is_simulation(self) -> bool:
        return self.kind == "simulation"


@dataclass(frozen=True)
class CentralityMap:
    campaign: str
    edges: np.ndarray
    values: dict[str, np.ndarray]

    def factor(self, family: str, centrality: np.ndarray) -> np.ndarray:
        """Per-event factor; NaN outside the support, error on bad inputs."""

        c = np.asarray(centrality, dtype=float)
        if np.any(~np.isfinite(c)):
            raise ValueError("Au+Au simulation event has a non-finite centrality")
        if family not in self.values:
            raise ValueError(f"no centrality map for family {family!r}")
        table = self.values[family]
        index = np.searchsorted(self.edges, c, side="right") - 1
        inside = (c >= self.edges[0]) & (c < self.edges[-1])
        out = np.full(c.shape, np.nan)
        out[inside] = table[index[inside]]
        return out


@dataclass(frozen=True)
class SampleManifest:
    campaign: str
    samples: dict[str, Sample]
    centrality: CentralityMap | None
    path: str


def load_manifest(path: Path | str) -> SampleManifest:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    references = raw.get("reference_cross_section_pb", {}) or {}
    samples: dict[str, Sample] = {}
    for name, spec in (raw.get("samples") or {}).items():
        system, kind, family = spec["system"], spec["kind"], spec["family"]
        if family not in FAMILIES:
            raise ValueError(f"{name}: unknown family {family!r}")
        if kind == "data":
            samples[name] = Sample(name, system, kind, family)
            continue
        factor_kind = spec.get("source_factor")
        if factor_kind not in SOURCE_FACTORS:
            raise ValueError(f"{name}: source_factor must be one of {SOURCE_FACTORS}")
        value: float | None
        window = None
        note = ""
        if factor_kind == "family_reference":
            generated = spec.get("generated_events")
            reference = references[spec["reference"]]
            if generated is None:
                value, note = None, "generated_events not recorded: event_weight only (relative normalisation)"
            else:
                value = float(reference) / float(generated)
        elif factor_kind == "cross_section_over_generated":
            value = float(spec["cross_section_pb"]) / float(spec["generated_events"])
        else:
            low, high = spec["ownership_window_gev"]
            window = (float(low), None if high is None else float(high))
            value = float(spec["ownership_effective_cross_section_pb"]) / float(spec["owned_events"])
        if value is not None and not (math.isfinite(value) and value > 0):
            raise ValueError(f"{name}: source factor must be finite and positive")
        samples[name] = Sample(name, system, kind, family, factor_kind, value, window, note)

    centrality = None
    block = raw.get("centrality_weights")
    if block:
        edges = np.asarray(block["bin_edges_percent"], dtype=float)
        values = {}
        for family in ("photonjet", "inclusivejet"):
            if family in block:
                table = np.asarray(block[family], dtype=float)
                if table.shape != (len(edges) - 1,):
                    raise ValueError(f"centrality map {family} does not match the bin edges")
                if np.any(~np.isfinite(table)) or np.any(table <= 0):
                    raise ValueError(f"centrality map {family} has a non-positive or non-finite value")
                values[family] = table
        centrality = CentralityMap(str(block.get("campaign", "")), edges, values)
    return SampleManifest(str(raw.get("campaign", "")), samples, centrality, str(path))


def leading_truth_jet_pt(truth_jets: Mapping[str, np.ndarray]) -> dict[tuple[int, int, int], float]:
    """Maximum R=0.4 truth-jet pT per event from the truthJets tree columns."""

    best: dict[tuple[int, int, int], float] = {}
    radius = np.asarray(truth_jets["truth_jet_radius"], dtype=float)
    pt = np.asarray(truth_jets["truth_jet_pt"], dtype=float)
    for s, hi, lo, r, value in zip(truth_jets["source_file_index"], truth_jets["event_id_hi"], truth_jets["event_id_lo"], radius, pt):
        if not math.isclose(r, TRUTH_JET_RADIUS, abs_tol=1e-9) or not math.isfinite(value) or value <= 0.0:
            continue
        key = (int(s), int(hi), int(lo))
        if value > best.get(key, -math.inf):
            best[key] = float(value)
    return best


def complete_event_weights(
    manifest: SampleManifest,
    sample_name: str,
    *,
    event_weight: np.ndarray,
    centrality: np.ndarray,
    leading_truth_jet_pt_gev: np.ndarray | None = None,
) -> dict[str, Any]:
    """Return the complete weight per event (NaN where the event is dropped).

    ``leading_truth_jet_pt_gev`` is required for ownership-stitched samples and
    may contain NaN for events without an eligible truth jet.
    """

    sample = manifest.samples.get(sample_name)
    if sample is None:
        raise KeyError(f"sample {sample_name!r} is not in {manifest.path}")
    w = np.asarray(event_weight, dtype=float).copy()
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError(f"{sample_name}: producer event_weight must be finite and non-negative")
    components = ["producer_event_weight"]
    dropped = {"outside_ownership_window": 0, "outside_centrality_support": 0}

    if sample.is_simulation:
        if sample.source_factor == "ownership_stitch":
            if leading_truth_jet_pt_gev is None:
                raise ValueError(f"{sample_name}: ownership stitching needs the leading truth-jet pT per event")
            pt = np.asarray(leading_truth_jet_pt_gev, dtype=float)
            low, high = sample.ownership_window
            owned = np.isfinite(pt) & (pt >= low) & (True if high is None else (pt < high))
            dropped["outside_ownership_window"] = int((~owned).sum())
            w[~owned] = np.nan
            w[owned] *= sample.factor_value
            components.append(f"source_stitch:{sample_name}:sigma_eff/owned_events={sample.factor_value:.8g}pb")
        elif sample.factor_value is not None:
            w *= sample.factor_value
            components.append(f"source_factor:{sample_name}:{sample.source_factor}={sample.factor_value:.8g}pb")
        else:
            components.append(f"source_factor:{sample_name}:not_applied({sample.note})")

        if sample.system == "auau":
            if manifest.centrality is None:
                raise ValueError(f"{sample_name}: Au+Au simulation needs centrality_weights in {manifest.path}")
            factor = manifest.centrality.factor(sample.family, centrality)
            alive = np.isfinite(w)
            outside = alive & ~np.isfinite(factor)
            dropped["outside_centrality_support"] = int(outside.sum())
            w[outside] = np.nan
            w[alive & ~outside] *= factor[alive & ~outside]
            components.append(f"centrality:{sample.family}:{manifest.centrality.campaign}")

    return {
        "sample": sample_name,
        "family": sample.family,
        "kind": sample.kind,
        "weight": w,
        "kept": np.isfinite(w),
        "components": tuple(components),
        "dropped": dropped,
        "campaign": manifest.campaign,
        "manifest": manifest.path,
    }


__all__ = [
    "CentralityMap",
    "Sample",
    "SampleManifest",
    "complete_event_weights",
    "leading_truth_jet_pt",
    "load_manifest",
]
