"""Common, information-preserving response semantics for p+p and Au+Au.

The measured unfolding state is intentionally distinct from classification
support.  Current aggregate adapters may carry unresolved witnesses; future
event adapters can use ``classify_state`` without changing downstream APIs.
ROOT underflow and overflow are never physical states in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


MEASURED_PTGAMMA_EDGES = np.asarray([15.0, 20.0, 25.0, 35.0])
CLASSIFICATION_RECO_PTGAMMA_EDGES = np.asarray([10.0, 15.0, 20.0, 25.0, 35.0, 40.0])
CLASSIFICATION_TRUTH_PTGAMMA_EDGES = np.asarray([5.0, 10.0, 15.0, 20.0, 25.0, 35.0, 40.0])
COMMON_XJ_EDGES = np.asarray([round(0.1 * i, 10) for i in range(21)] + [3.0])
REPORTED_XJ_LOW = 0.3
REPORTED_XJ_HIGH = 1.8


class Category(str, Enum):
    MATCHED_CENTRAL = "MATCHED_CENTRAL"
    PHOTON_LOW_FEED_IN = "PHOTON_LOW_FEED_IN"
    PHOTON_HIGH_FEED_IN = "PHOTON_HIGH_FEED_IN"
    PHOTON_LOW_FEED_OUT = "PHOTON_LOW_FEED_OUT"
    PHOTON_HIGH_FEED_OUT = "PHOTON_HIGH_FEED_OUT"
    PHOTON_RECO_MISS = "PHOTON_RECO_MISS"
    UNMATCHED_RECO = "UNMATCHED_RECO"
    COMBINATORIC = "COMBINATORIC"
    ASSIGNED_HARD_NONFIDUCIAL = "ASSIGNED_HARD_NONFIDUCIAL"
    XJ_LOW_FEED_IN = "XJ_LOW_FEED_IN"
    XJ_HIGH_FEED_IN = "XJ_HIGH_FEED_IN"
    XJ_LOW_FEED_OUT = "XJ_LOW_FEED_OUT"
    XJ_HIGH_FEED_OUT = "XJ_HIGH_FEED_OUT"
    XJ_GE3_BOUNDARY = "XJ_GE3_BOUNDARY"
    OUTSIDE_OUTSIDE = "OUTSIDE_OUTSIDE"
    UNRESOLVED_FROM_INPUT = "UNRESOLVED_FROM_INPUT"


def _region(value: float | None, low: float, high: float) -> str:
    if value is None:
        return "NONE"
    if value < low:
        return "LOW"
    if value >= high:
        return "HIGH"
    return "CENTRAL"


def classify_state(
    *,
    truth_ptgamma: float | None,
    reco_ptgamma: float | None,
    truth_xj: float | None = None,
    reco_xj: float | None = None,
) -> tuple[Category, ...]:
    """Classify one future event/pair with identical pp/Au+Au semantics.

    Photon categories are exclusive.  An optional xJ category is an orthogonal
    label on a matched photon/jet pair, not a second response destination.
    """
    truth_region = _region(truth_ptgamma, 15.0, 35.0)
    reco_region = _region(reco_ptgamma, 15.0, 35.0)
    if truth_region == "NONE" and reco_region == "CENTRAL":
        photon = Category.UNMATCHED_RECO
    elif truth_region == "CENTRAL" and reco_region == "NONE":
        photon = Category.PHOTON_RECO_MISS
    elif truth_region == "LOW" and reco_region == "CENTRAL":
        photon = Category.PHOTON_LOW_FEED_IN
    elif truth_region == "HIGH" and reco_region == "CENTRAL":
        photon = Category.PHOTON_HIGH_FEED_IN
    elif truth_region == "CENTRAL" and reco_region == "LOW":
        photon = Category.PHOTON_LOW_FEED_OUT
    elif truth_region == "CENTRAL" and reco_region == "HIGH":
        photon = Category.PHOTON_HIGH_FEED_OUT
    elif truth_region == "CENTRAL" and reco_region == "CENTRAL":
        photon = Category.MATCHED_CENTRAL
    else:
        photon = Category.OUTSIDE_OUTSIDE

    labels = [photon]
    if truth_xj is None and reco_xj is None:
        return tuple(labels)
    if (truth_xj is not None and truth_xj >= 3.0) or (reco_xj is not None and reco_xj >= 3.0):
        labels.append(Category.XJ_GE3_BOUNDARY)
        return tuple(labels)
    truth_x = _region(truth_xj, REPORTED_XJ_LOW, REPORTED_XJ_HIGH)
    reco_x = _region(reco_xj, REPORTED_XJ_LOW, REPORTED_XJ_HIGH)
    mapping = {
        ("LOW", "CENTRAL"): Category.XJ_LOW_FEED_IN,
        ("HIGH", "CENTRAL"): Category.XJ_HIGH_FEED_IN,
        ("CENTRAL", "LOW"): Category.XJ_LOW_FEED_OUT,
        ("CENTRAL", "HIGH"): Category.XJ_HIGH_FEED_OUT,
    }
    if (truth_x, reco_x) in mapping:
        labels.append(mapping[(truth_x, reco_x)])
    return tuple(labels)


def classify_photon_marginals(
    *, truth_ptgamma: float | None, reco_ptgamma: float | None, matched: bool
) -> tuple[Category, ...]:
    """Return the exhaustive 1D photon marginal destinations for one event.

    A fake reconstructed photon and an independently present central truth
    photon occupy one reco marginal and one truth marginal, respectively.  The
    two labels are therefore not double counting within either marginal.
    """
    truth_region = _region(truth_ptgamma, 15.0, 35.0)
    reco_region = _region(reco_ptgamma, 15.0, 35.0)
    if matched and (truth_region == "NONE" or reco_region == "NONE"):
        # A missing reconstructed photon is represented by matched=False: the
        # retained link cannot name a nonexistent object.
        matched = False
    labels: list[Category] = []
    if matched:
        if reco_region == "CENTRAL":
            if truth_region == "CENTRAL":
                labels.append(Category.MATCHED_CENTRAL)
            elif truth_region == "LOW":
                labels.append(Category.PHOTON_LOW_FEED_IN)
            elif truth_region == "HIGH":
                labels.append(Category.PHOTON_HIGH_FEED_IN)
        if truth_region == "CENTRAL" and reco_region != "CENTRAL":
            labels.append(
                Category.PHOTON_LOW_FEED_OUT
                if reco_region == "LOW"
                else Category.PHOTON_HIGH_FEED_OUT
                if reco_region == "HIGH"
                else Category.PHOTON_RECO_MISS
            )
    else:
        if reco_region == "CENTRAL":
            labels.append(Category.UNMATCHED_RECO)
        if truth_region == "CENTRAL":
            labels.append(Category.PHOTON_RECO_MISS)
    if not labels and (truth_region != "NONE" or reco_region != "NONE"):
        labels.append(Category.OUTSIDE_OUTSIDE)
    return tuple(labels)


@dataclass(frozen=True)
class GlobalBinMap:
    ptgamma_edges: np.ndarray
    xj_edges: np.ndarray

    def __post_init__(self) -> None:
        p = np.asarray(self.ptgamma_edges, dtype=float)
        x = np.asarray(self.xj_edges, dtype=float)
        if len(p) < 2 or len(x) < 2 or np.any(np.diff(p) <= 0) or np.any(np.diff(x) <= 0):
            raise ValueError("strictly increasing axes required")
        object.__setattr__(self, "ptgamma_edges", p)
        object.__setattr__(self, "xj_edges", x)

    @property
    def n_ptgamma(self) -> int:
        return len(self.ptgamma_edges) - 1

    @property
    def n_xj(self) -> int:
        return len(self.xj_edges) - 1

    @property
    def n_global(self) -> int:
        return self.n_ptgamma * self.n_xj

    def flatten(self, i_ptgamma: int, i_xj: int) -> int:
        if not (0 <= i_ptgamma < self.n_ptgamma and 0 <= i_xj < self.n_xj):
            raise IndexError((i_ptgamma, i_xj))
        return i_ptgamma * self.n_xj + i_xj

    def invert(self, global_bin: int) -> tuple[int, int]:
        if not 0 <= global_bin < self.n_global:
            raise IndexError(global_bin)
        return divmod(global_bin, self.n_xj)


def _array(value: Any, shape: tuple[int, ...], name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.shape != shape or np.any(~np.isfinite(result)):
        raise ValueError(f"{name} shape/finite contract differs: {result.shape} != {shape}")
    return result


@dataclass
class ResponseBundle:
    system: str
    dimension: str
    truth_ptgamma_edges: np.ndarray
    reco_ptgamma_edges: np.ndarray
    xj_edges: np.ndarray | None
    matrix: np.ndarray
    matrix_sumw2: np.ndarray
    truth: np.ndarray
    truth_sumw2: np.ndarray
    reco: np.ndarray
    reco_sumw2: np.ndarray
    misses: np.ndarray
    misses_sumw2: np.ndarray
    fakes: np.ndarray
    fakes_sumw2: np.ndarray
    fake_causes: dict[str, np.ndarray] = field(default_factory=dict)
    fake_causes_sumw2: dict[str, np.ndarray] = field(default_factory=dict)
    boundary_reco: dict[str, np.ndarray] = field(default_factory=dict)
    boundary_reco_sumw2: dict[str, np.ndarray] = field(default_factory=dict)
    boundary_truth: dict[str, np.ndarray] = field(default_factory=dict)
    boundary_truth_sumw2: dict[str, np.ndarray] = field(default_factory=dict)
    unfolding_fakes_override: np.ndarray | None = None
    unresolved: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    status: str = "INPUT_CONTRACT_INCOMPLETE"

    def __post_init__(self) -> None:
        self.truth_ptgamma_edges = np.asarray(self.truth_ptgamma_edges, dtype=float)
        self.reco_ptgamma_edges = np.asarray(self.reco_ptgamma_edges, dtype=float)
        self.xj_edges = None if self.xj_edges is None else np.asarray(self.xj_edges, dtype=float)
        nt = len(self.truth_ptgamma_edges) - 1
        nr = len(self.reco_ptgamma_edges) - 1
        if self.dimension == "2D":
            if self.xj_edges is None:
                raise ValueError("2D response requires xJ edges")
            nx = len(self.xj_edges) - 1
            nt *= nx
            nr *= nx
        elif self.dimension != "1D":
            raise ValueError(self.dimension)
        self.matrix = _array(self.matrix, (nt, nr), "matrix")
        self.matrix_sumw2 = _array(self.matrix_sumw2, (nt, nr), "matrix_sumw2")
        self.truth = _array(self.truth, (nt,), "truth")
        self.truth_sumw2 = _array(self.truth_sumw2, (nt,), "truth_sumw2")
        self.reco = _array(self.reco, (nr,), "reco")
        self.reco_sumw2 = _array(self.reco_sumw2, (nr,), "reco_sumw2")
        self.misses = _array(self.misses, (nt,), "misses")
        self.misses_sumw2 = _array(self.misses_sumw2, (nt,), "misses_sumw2")
        self.fakes = _array(self.fakes, (nr,), "fakes")
        self.fakes_sumw2 = _array(self.fakes_sumw2, (nr,), "fakes_sumw2")
        self.fake_causes = {k: _array(v, (nr,), k) for k, v in self.fake_causes.items()}
        self.fake_causes_sumw2 = {
            k: _array(v, (nr,), k) for k, v in self.fake_causes_sumw2.items()
        }
        if self.fake_causes and not np.allclose(
            sum(self.fake_causes.values(), np.zeros(nr)), self.fakes, rtol=1e-9, atol=1e-9
        ):
            raise ValueError("fake-cause partition does not equal fakes")
        if self.fake_causes_sumw2 and not np.allclose(
            sum(self.fake_causes_sumw2.values(), np.zeros(nr)),
            self.fakes_sumw2,
            rtol=1e-9,
            atol=1e-9,
        ):
            raise ValueError("fake-cause sumw2 partition does not equal fakes_sumw2")
        self.boundary_reco = {k: _array(v, (nr,), k) for k, v in self.boundary_reco.items()}
        self.boundary_reco_sumw2 = {
            k: _array(v, (nr,), k) for k, v in self.boundary_reco_sumw2.items()
        }
        self.boundary_truth = {k: _array(v, (nt,), k) for k, v in self.boundary_truth.items()}
        self.boundary_truth_sumw2 = {
            k: _array(v, (nt,), k) for k, v in self.boundary_truth_sumw2.items()
        }
        if self.unfolding_fakes_override is not None:
            self.unfolding_fakes_override = _array(
                self.unfolding_fakes_override, (nr,), "unfolding_fakes_override"
            )
        for collection in (
            self.matrix, self.matrix_sumw2, self.truth, self.truth_sumw2,
            self.reco, self.reco_sumw2, self.misses, self.misses_sumw2,
            self.fakes, self.fakes_sumw2, *self.boundary_reco.values(),
            *self.boundary_truth.values(),
        ):
            if np.any(collection < -1.0e-12):
                raise ValueError("response accumulators must be nonnegative")

    @property
    def all_fakes(self) -> np.ndarray:
        result = self.fakes.copy()
        for value in self.boundary_reco.values():
            result += value
        return result

    @property
    def all_misses(self) -> np.ndarray:
        result = self.misses.copy()
        for value in self.boundary_truth.values():
            result += value
        return result

    @property
    def all_fakes_sumw2(self) -> np.ndarray:
        result = self.fakes_sumw2.copy()
        for value in self.boundary_reco_sumw2.values():
            result += value
        return result

    @property
    def all_misses_sumw2(self) -> np.ndarray:
        result = self.misses_sumw2.copy()
        for value in self.boundary_truth_sumw2.values():
            result += value
        return result

    @property
    def unfolding_fakes(self) -> np.ndarray:
        """Fake cause consumed by Bayes after system-specific DATA correction."""
        return self.all_fakes if self.unfolding_fakes_override is None else self.unfolding_fakes_override

    def conservation(self) -> dict[str, Any]:
        reco_partition = self.matrix.sum(axis=0) + self.all_fakes
        truth_partition = self.matrix.sum(axis=1) + self.all_misses
        reco_partition_sumw2 = self.matrix_sumw2.sum(axis=0) + self.all_fakes_sumw2
        truth_partition_sumw2 = self.matrix_sumw2.sum(axis=1) + self.all_misses_sumw2
        reco_delta = self.reco - reco_partition
        truth_delta = self.truth - truth_partition
        reco_sumw2_delta = self.reco_sumw2 - reco_partition_sumw2
        truth_sumw2_delta = self.truth_sumw2 - truth_partition_sumw2
        return {
            "reco_max_abs_residual": float(np.max(np.abs(reco_delta), initial=0.0)),
            "truth_max_abs_residual": float(np.max(np.abs(truth_delta), initial=0.0)),
            "reco_sum_residual": float(reco_delta.sum()),
            "truth_sum_residual": float(truth_delta.sum()),
            "reco_sumw2_max_abs_residual": float(
                np.max(np.abs(reco_sumw2_delta), initial=0.0)
            ),
            "truth_sumw2_max_abs_residual": float(
                np.max(np.abs(truth_sumw2_delta), initial=0.0)
            ),
            "reco_sumw2_sum_residual": float(reco_sumw2_delta.sum()),
            "truth_sumw2_sum_residual": float(truth_sumw2_delta.sum()),
            "closes_rtol_1e-9_atol_1e-9": bool(
                np.allclose(self.reco, reco_partition, rtol=1e-9, atol=1e-9)
                and np.allclose(self.truth, truth_partition, rtol=1e-9, atol=1e-9)
                and np.allclose(self.reco_sumw2, reco_partition_sumw2, rtol=1e-9, atol=1e-9)
                and np.allclose(self.truth_sumw2, truth_partition_sumw2, rtol=1e-9, atol=1e-9)
            ),
            "closes_rtol_1e-10_atol_1e-10": bool(
                np.allclose(self.reco, reco_partition, rtol=1e-10, atol=1e-10)
                and np.allclose(self.truth, truth_partition, rtol=1e-10, atol=1e-10)
                and np.allclose(self.reco_sumw2, reco_partition_sumw2, rtol=1e-10, atol=1e-10)
                and np.allclose(self.truth_sumw2, truth_partition_sumw2, rtol=1e-10, atol=1e-10)
            ),
        }

    def diagnostics(self) -> dict[str, Any]:
        total = self.matrix.sum()
        occupied = self.matrix > 0
        effective = np.divide(
            self.matrix * self.matrix,
            self.matrix_sumw2,
            out=np.zeros_like(self.matrix),
            where=self.matrix_sumw2 > 0,
        )
        totals = self.matrix.sum(axis=1) + self.all_misses
        conditional = np.divide(
            self.matrix, totals[:, None], out=np.zeros_like(self.matrix), where=totals[:, None] > 0
        )
        active_rows = np.any(conditional > 0, axis=1)
        active_cols = np.any(conditional > 0, axis=0)
        active = conditional[np.ix_(active_rows, active_cols)]
        singular = np.linalg.svd(active, compute_uv=False) if active.size else np.asarray([])
        nonzero_singular = singular[singular > max(active.shape, default=0) * np.finfo(float).eps * (singular[0] if singular.size else 1)]
        condition = (
            float(nonzero_singular[0] / nonzero_singular[-1])
            if nonzero_singular.size else None
        )
        efficiency = np.divide(
            self.matrix.sum(axis=1), totals, out=np.zeros_like(totals), where=totals > 0
        )
        return {
            "shape": list(self.matrix.shape),
            "sumw": float(total),
            "sumw2": float(self.matrix_sumw2.sum()),
            "occupied_cells": int(occupied.sum()),
            "empty_cells": int(occupied.size - occupied.sum()),
            "occupancy_fraction": float(occupied.mean()),
            "near_empty_effective_lt_5": int(np.sum(occupied & (effective < 5.0))),
            "effective_occupancy_sum": float(effective.sum()),
            "rank": int(np.linalg.matrix_rank(active)) if active.size else 0,
            "active_shape": list(active.shape),
            "condition_nonzero_singular": condition,
            "zero_efficiency_states": np.flatnonzero((totals > 0) & (efficiency == 0)).tolist(),
            "efficiency": efficiency.tolist(),
            "conservation": self.conservation(),
        }

    def save(self, stem: Path) -> tuple[Path, Path]:
        stem.parent.mkdir(parents=True, exist_ok=True)
        npz = stem.with_suffix(".npz")
        meta = stem.with_suffix(".json")
        arrays: dict[str, Any] = {
            "truth_ptgamma_edges": self.truth_ptgamma_edges,
            "reco_ptgamma_edges": self.reco_ptgamma_edges,
            "xj_edges": np.asarray([]) if self.xj_edges is None else self.xj_edges,
            "matrix": self.matrix,
            "matrix_sumw2": self.matrix_sumw2,
            "truth": self.truth,
            "truth_sumw2": self.truth_sumw2,
            "reco": self.reco,
            "reco_sumw2": self.reco_sumw2,
            "misses": self.misses,
            "misses_sumw2": self.misses_sumw2,
            "fakes": self.fakes,
            "fakes_sumw2": self.fakes_sumw2,
            "unfolding_fakes": self.unfolding_fakes,
        }
        for name, value in self.fake_causes.items():
            arrays[f"fake_causes__{name}"] = value
        for name, value in self.fake_causes_sumw2.items():
            arrays[f"fake_causes_sumw2__{name}"] = value
        for name, value in self.boundary_reco.items():
            arrays[f"boundary_reco__{name}"] = value
        for name, value in self.boundary_reco_sumw2.items():
            arrays[f"boundary_reco_sumw2__{name}"] = value
        for name, value in self.boundary_truth.items():
            arrays[f"boundary_truth__{name}"] = value
        for name, value in self.boundary_truth_sumw2.items():
            arrays[f"boundary_truth_sumw2__{name}"] = value
        np.savez_compressed(npz, **arrays)
        payload = {
            "schema": "PhotonJetResponseBundleV1",
            "status": self.status,
            "system": self.system,
            "dimension": self.dimension,
            "global_bin_index_formula": "g = i_ptgamma * N_xJ + i_xJ",
            "root_underflow_policy": "NOT_A_PHYSICAL_STATE",
            "root_overflow_policy": "NOT_A_PHYSICAL_STATE",
            "unresolved": self.unresolved,
            "provenance": self.provenance,
            "boundary_reco_categories": sorted(self.boundary_reco),
            "boundary_truth_categories": sorted(self.boundary_truth),
            "diagnostics": self.diagnostics(),
        }
        meta.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz, meta


def exact_edge_ranges(source: np.ndarray, target: np.ndarray) -> list[tuple[int, int]]:
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    indices: list[int] = []
    for edge in target:
        match = np.flatnonzero(np.isclose(source, edge, rtol=0.0, atol=1e-9))
        if len(match) != 1:
            raise ValueError(f"target edge {edge} is not exactly represented")
        indices.append(int(match[0]))
    return list(zip(indices[:-1], indices[1:]))


def rebin_2d(values: np.ndarray, source_x: np.ndarray, source_y: np.ndarray,
             target_x: np.ndarray, target_y: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    xr = exact_edge_ranges(source_x, target_x)
    yr = exact_edge_ranges(source_y, target_y)
    return np.asarray([[values[x0:x1, y0:y1].sum() for y0, y1 in yr] for x0, x1 in xr])


def rebin_global_matrix(
    matrix: np.ndarray,
    truth_pt_source: np.ndarray,
    reco_pt_source: np.ndarray,
    xj_source: np.ndarray,
    truth_pt_target: np.ndarray,
    reco_pt_target: np.ndarray,
    xj_target: np.ndarray,
) -> np.ndarray:
    nt, nr, nx = len(truth_pt_source) - 1, len(reco_pt_source) - 1, len(xj_source) - 1
    source = np.asarray(matrix, dtype=float).reshape(nt, nx, nr, nx)
    tr = exact_edge_ranges(truth_pt_source, truth_pt_target)
    rr = exact_edge_ranges(reco_pt_source, reco_pt_target)
    xr = exact_edge_ranges(xj_source, xj_target)
    out = np.zeros(((len(truth_pt_target) - 1) * (len(xj_target) - 1),
                    (len(reco_pt_target) - 1) * (len(xj_target) - 1)))
    mapping_t = GlobalBinMap(truth_pt_target, xj_target)
    mapping_r = GlobalBinMap(reco_pt_target, xj_target)
    for it, (t0, t1) in enumerate(tr):
        for ix, (x0, x1) in enumerate(xr):
            gt = mapping_t.flatten(it, ix)
            for ir, (r0, r1) in enumerate(rr):
                for iy, (y0, y1) in enumerate(xr):
                    gr = mapping_r.flatten(ir, iy)
                    out[gt, gr] = source[t0:t1, x0:x1, r0:r1, y0:y1].sum()
    return out
