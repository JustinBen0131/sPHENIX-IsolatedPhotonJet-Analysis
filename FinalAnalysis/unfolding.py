"""Iterative-Bayes unfolding kernel and the analysis-window restriction.

Ported from the validated collaboration implementation
(``photonjet/analysis/unfolding.py`` and the ``restrict`` step of
``photonjet/analysis/chain.py`` at main a7f15eb). The kernel is unchanged;
the restriction reads a ``TreeToHists`` package instead of a response
bundle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Sequence

import numpy as np


def iterative_bayes(
    measured: Sequence[float],
    response_truth_x_reco: np.ndarray,
    misses: Sequence[float],
    iterations: int,
    prior: Sequence[float] | None = None,
    *,
    fakes: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """RooUnfold-compatible D'Agostini update including one fake cause."""

    measured_array = np.asarray(measured, dtype=float)
    matrix = np.asarray(response_truth_x_reco, dtype=float)
    miss = np.asarray(misses, dtype=float)
    if matrix.ndim != 2 or matrix.shape != (miss.size, measured_array.size):
        raise ValueError("unfolding dimensions differ")
    if np.any(~np.isfinite(measured_array)) or np.any(~np.isfinite(matrix)) or np.any(~np.isfinite(miss)) or np.any(matrix < 0) or np.any(miss < 0):
        raise ValueError("non-finite measured values or invalid response counts")
    if iterations < 1:
        raise ValueError("iterations must be positive")
    totals = matrix.sum(axis=1) + miss
    conditional = np.divide(matrix, totals[:, None], out=np.zeros_like(matrix), where=totals[:, None] > 0)
    efficiency = conditional.sum(axis=1)
    state_counts = totals.copy() if prior is None else np.asarray(prior, dtype=float).copy()
    if np.any(~np.isfinite(state_counts)) or np.any(state_counts < 0):
        raise ValueError("prior must be finite and nonnegative")
    if state_counts.shape != totals.shape:
        raise ValueError("prior dimensions differ")
    ntruth = matrix.shape[0]
    fake = np.zeros(matrix.shape[1]) if fakes is None else np.asarray(fakes, dtype=float)
    if fake.shape != (matrix.shape[1],) or np.any(fake < -1e-12) or np.any(~np.isfinite(fake)):
        raise ValueError("fake cause differs")
    fake = np.clip(fake, 0.0, None)
    if fake.sum() > 0:
        conditional = np.vstack((conditional, fake / fake.sum()))
        efficiency = np.append(efficiency, 1.0)
        state_counts = np.append(state_counts, fake.sum())
    if state_counts.sum() <= 0:
        raise ValueError("empty prior")
    state = state_counts / state_counts.sum()
    unfolded = np.zeros(conditional.shape[0])
    for _ in range(iterations):
        denominator = conditional.T @ state
        posterior = np.divide(conditional * state[:, None], denominator[None, :],
                              out=np.zeros_like(conditional), where=denominator[None, :] > 0)
        unfolded = np.divide(posterior @ measured_array, efficiency, out=np.zeros(conditional.shape[0]), where=efficiency > 0)
        total = float(unfolded.sum())
        if not math.isfinite(total) or abs(total) <= 1e-15:
            raise ValueError("empty signed estimate")
        state = unfolded / total
    return unfolded[:ntruth], conditional.T @ unfolded


# ----------------------------------------------------------------------------
# Analysis-window restriction of the package responses
# ----------------------------------------------------------------------------

@dataclass
class AnalysisResponse:
    """Pair and photon responses cut to the measurement photon-pT window."""

    pt_edges: np.ndarray
    xj_edges: np.ndarray
    xj_response: np.ndarray            # (n_pt nx, n_pt nx) truth x reco
    xj_misses: np.ndarray              # (n_pt nx)
    xj_truth: np.ndarray               # (n_pt, nx)
    xj_reco: np.ndarray                # (n_pt, nx)
    xj_boundary_fakes: np.ndarray      # (n_pt, nx) matched pairs with truth outside the window
    xj_detector_fakes: np.ndarray      # (n_pt, nx) linked photon and jet, not the selected truth pair
    xj_fake_photon: np.ndarray         # (n_pt, nx) reco photon without a truth link (removed by ABCD)
    xj_combinatoric: np.ndarray        # (n_pt, nx) linked photon, jet without a truth jet (removed by K)
    xj_combinatoric_sumw2: np.ndarray
    photon_response: np.ndarray        # (n_pt, n_pt) truth x reco
    photon_misses: np.ndarray
    photon_truth: np.ndarray
    photon_reco: np.ndarray
    photon_boundary_fakes: np.ndarray
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def n_pt(self) -> int:
        return len(self.pt_edges) - 1

    @property
    def n_xj(self) -> int:
        return len(self.xj_edges) - 1

    @property
    def unfolding_fakes(self) -> np.ndarray:
        return (self.xj_detector_fakes + self.xj_boundary_fakes).reshape(-1)


def _window_projection(edges: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Add fine support bins; never interpolate a boundary lost at histogram time."""
    for edge in target:
        if not np.any(np.isclose(edges, edge, rtol=0, atol=1e-9)):
            raise ValueError(f"measurement boundary {edge} is absent from the retained response support")
    return np.asarray([(edges[:-1] >= low - 1e-9) & (edges[1:] <= high + 1e-9)
                       for low, high in zip(target[:-1], target[1:])], dtype=float)


def restrict(package) -> AnalysisResponse:
    """Aggregate fine response support into the measurement window.

    Reco migration out of the window joins misses; truth migration in joins
    boundary fakes. Counts and sumw2 are added before any nonlinear operation.
    """
    g, a = package.grids, package.arrays
    if "pair_response" not in a:
        raise ValueError("the package carries no simulation response")
    tp = _window_projection(g.truth_pt_edges, g.pt_edges)
    rp = _window_projection(g.reco_pt_edges, g.pt_edges)
    tpair, rpair = np.kron(tp, np.eye(g.n_xj)), np.kron(rp, np.eye(g.n_xj))
    t_out, r_out = 1.0 - tpair.sum(axis=0), 1.0 - rpair.sum(axis=0)
    matrix = a["pair_response"]
    n_pt, nx = g.n_pt, g.n_xj
    response = tpair @ matrix @ rpair.T
    misses = tpair @ (package.all_pair_misses() + matrix @ r_out)
    boundary = (rpair @ (t_out @ matrix)).reshape(n_pt, nx)
    detector = np.zeros((n_pt, nx))
    for key, values in package.boundary_reco().items():
        block = (rpair @ values).reshape(n_pt, nx)
        if key == "ASSIGNED_HARD_NONFIDUCIAL":
            detector += block
        else:
            boundary += block
    causes = {name: (rpair @ values).reshape(n_pt, nx) for name, values in package.fake_causes().items()}
    zeros = np.zeros((n_pt, nx))
    photon_matrix = a["photon_response"]
    photon_misses = tp @ (a["photon_misses"] + photon_matrix @ (1.0-rp.sum(axis=0)))
    photon_boundary = rp @ (a["photon_boundary_fakes"] + (1.0-tp.sum(axis=0)) @ photon_matrix)
    return AnalysisResponse(
        pt_edges=g.pt_edges, xj_edges=g.xj_edges,
        xj_response=response, xj_misses=misses,
        xj_truth=(tpair @ a["pair_truth"]).reshape(n_pt, nx),
        xj_reco=(rpair @ a["pair_reco"]).reshape(n_pt, nx),
        xj_boundary_fakes=boundary, xj_detector_fakes=detector,
        xj_fake_photon=causes.get("UNMATCHED_RECO", zeros),
        xj_combinatoric=causes.get("COMBINATORIC", zeros),
        xj_combinatoric_sumw2=(rpair @ a["pair_fakes_COMBINATORIC_sumw2"]).reshape(n_pt, nx),
        photon_response=tp @ photon_matrix @ rp.T, photon_misses=photon_misses,
        photon_truth=tp @ a["truth_photons"], photon_reco=rp @ a["photon_reco"],
        photon_boundary_fakes=photon_boundary, provenance={"package": package.metadata},
    )


def unfold_pairs(measured: np.ndarray, window: AnalysisResponse, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    return iterative_bayes(np.asarray(measured).reshape(-1), window.xj_response, window.xj_misses, iterations,
                           window.xj_truth.reshape(-1), fakes=window.unfolding_fakes)


def unfold_photons(photons: np.ndarray, window: AnalysisResponse, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    return iterative_bayes(np.asarray(photons), window.photon_response, window.photon_misses, iterations,
                           window.photon_truth, fakes=window.photon_boundary_fakes)


__all__ = ["AnalysisResponse", "iterative_bayes", "restrict", "unfold_pairs", "unfold_photons"]
