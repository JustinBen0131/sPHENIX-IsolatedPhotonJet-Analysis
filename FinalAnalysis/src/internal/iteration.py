"""Iteration-count scan with its frozen gate contract, and the photon-yield refold test.

Ported from the validated collaboration implementation: ``scan_problem`` and ``gate_contract`` of
``photonjet/analysis/unfolding.py`` (frozen 2026-09-20) and the measured-covariance refold test that the
2026-09-21 review of the one-dimensional gate adopted for the photon-yield problem. Definitions, gate values and
the selection rule are unchanged. Unfolding goes through ``unfolding.unfold_photons`` and
``unfolding.unfold_pairs``; chi2, the reported mask and the central-68 half width come from ``statistics``.

A scan unfolds one problem of an ``AnalysisResponse`` for k = 1..12:

  ``"1D"``  the purity-corrected photon yield per photon-pT bin (``unfold_photons``);
  ``"2D"``  the joint (photon pT, xJ) pair spectrum (``unfold_pairs``).

Every row records

  S_data      median, over reported bins with a positive central value, of the central-68 half width of the
              unfolded DATA toys divided by the central unfolded value;
  S_response  the same for the unfolded response toys;
  D_k         median symmetric relative change of the central unfolded spectrum from k-1 to k, over reported
              bins positive at both;
  Q_k         sqrt(S_data^2 + S_response^2 + D_k^2).

Only k = 2..4 can be selected: the candidate with the smallest Q_k among those passing every gate, ties to the
smaller k. k = 5..12 are diagnostics. When no candidate passes nothing is selected; the cause belongs upstream
and is never absorbed by a larger k.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

import numpy as np

from . import statistics as _statistics
from . import unfolding as _unfolding

unfold_photons = _unfolding.unfold_photons
unfold_pairs = _unfolding.unfold_pairs
central_68 = _statistics.central_68
chi2_ndf = _statistics.chi2_ndf
reported_mask = _statistics.reported_mask
covariance_chi2_ndf_from_matrix = _statistics.covariance_chi2_ndf_from_matrix

COMPUTED_ITERATIONS = tuple(range(1, 13))
CANDIDATE_ITERATIONS = tuple(range(2, 5))
TOY_SUCCESS_MIN = 0.95
REFOLD_CHI2_NDF_MAX = 1.5
MC_CHI2_NDF_MAX = 5.0
NEGATIVE_FRACTION_MAX = 0.50
CONDITION_MAX = 1.0e8
# Every scan flags the row of this k as ``reference_selected``; a comparison flag only, it never selects.
REFERENCE_ITERATIONS = 12
RESPONSE_STAT_STATUS = "SOURCE_OCCURRENCE_BOOTSTRAP"
REFOLD_COMPONENTS = "DATA_RESIDUAL_TOYS_PLUS_RESPONSE_LEAKAGE_TCOMB_RESIDUAL_TOYS"
REFOLD_COMPONENTS_INCOMPLETE = "INCOMPLETE_RESPONSE_RESIDUAL_COVARIANCE"
PROBLEMS = ("1D", "2D")


def problem_arrays(window, dimension: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Response (truth x reco), misses, truth prior and unfolding fake cause of one problem of a window."""

    if dimension == "1D":
        return window.photon_response, window.photon_misses, window.photon_truth, window.photon_boundary_fakes
    if dimension == "2D":
        return window.xj_response, window.xj_misses, window.xj_truth.reshape(-1), window.unfolding_fakes
    raise ValueError(f"unknown unfolding problem {dimension!r}; expected one of {PROBLEMS}")


def unfold(measured: np.ndarray, window, dimension: str, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    """(unfolded, refolded) of one problem through the public unfolding kernels."""

    if dimension == "1D":
        return unfold_photons(measured, window, iterations)
    if dimension == "2D":
        return unfold_pairs(measured, window, iterations)
    raise ValueError(f"unknown unfolding problem {dimension!r}; expected one of {PROBLEMS}")


def reported_bins(window, dimension: str, reported_xj: tuple[float, float]) -> np.ndarray:
    """Every photon-pT bin of the photon problem; the reported xJ bins of every photon-pT bin of the pair problem."""

    if dimension == "1D":
        return np.ones(window.n_pt, dtype=bool)
    if dimension == "2D":
        return np.tile(reported_mask(window.xj_edges, *reported_xj), window.n_pt)
    raise ValueError(f"unknown unfolding problem {dimension!r}; expected one of {PROBLEMS}")


def conditional_probabilities(matrix: np.ndarray, misses: np.ndarray) -> np.ndarray:
    """P(reco bin | truth bin) with misses in the normalisation; empty truth rows stay zero."""

    totals = matrix.sum(axis=1) + misses
    return np.divide(matrix, totals[:, None], out=np.zeros_like(matrix), where=totals[:, None] > 0)


def condition_number(matrix: np.ndarray, misses: np.ndarray) -> float | None:
    """Largest over smallest nonzero singular value of the conditional response on its occupied rows and columns."""

    conditional = conditional_probabilities(matrix, misses)
    active = conditional[np.ix_(np.any(conditional > 0, axis=1), np.any(conditional > 0, axis=0))]
    singular = np.linalg.svd(active, compute_uv=False) if active.size else np.asarray([])
    floor = max(active.shape, default=0) * np.finfo(float).eps * (singular[0] if singular.size else 1)
    nonzero = singular[singular > floor]
    return float(nonzero[0] / nonzero[-1]) if nonzero.size else None


def zero_efficiency_states(matrix: np.ndarray, misses: np.ndarray, mask: np.ndarray) -> int:
    """Reported truth bins that are populated but never reconstructed inside the window."""

    totals = matrix.sum(axis=1) + misses
    efficiency = np.divide(matrix.sum(axis=1), totals, out=np.zeros_like(totals), where=totals > 0)
    return int(np.sum(mask & (totals > 0) & (efficiency == 0)))


def iteration_deviation(current: np.ndarray, previous: np.ndarray, mask: np.ndarray) -> float:
    """Median symmetric relative change over reported bins positive at k and k-1."""

    change = current - previous
    scale = np.maximum(0.5 * (np.abs(current) + np.abs(previous)), 1e-12)
    valid = mask & np.isfinite(change) & np.isfinite(scale) & (current > 0) & (previous > 0)
    values = np.abs(change[valid]) / scale[valid]
    return float(np.median(values)) if values.size else math.inf


def relative_half_width(central: np.ndarray, toys: np.ndarray, mask: np.ndarray) -> float:
    """Median central-68 half width of unfolded toys relative to the central value, over positive reported bins."""

    toys = np.asarray(toys, dtype=float)
    if toys.ndim != 2 or not len(toys):
        return math.inf
    _, _, half = central_68(toys)
    valid = mask & np.isfinite(central) & np.isfinite(half) & (central > 0)
    values = half[valid] / np.maximum(np.abs(central[valid]), 1e-12)
    return float(np.median(values)) if values.size else math.inf


def covariance_status(toys: np.ndarray) -> tuple[bool, bool, float | None]:
    """Symmetry, positive semi-definiteness and min(0, smallest eigenvalue) of the toy covariance.

    PSD allows lambda_min >= -1e-10 max|lambda|; an exactly zero covariance has a zero floor.
    """

    if toys.ndim != 2 or len(toys) < 2 or np.any(~np.isfinite(toys)):
        return False, False, None
    covariance = np.atleast_2d(np.cov(toys, rowvar=False, ddof=1))
    symmetric = bool(np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12))
    eigenvalues = np.linalg.eigvalsh(0.5 * (covariance + covariance.T))
    scale = float(np.max(np.abs(eigenvalues), initial=0.0))
    floor = -1e-10 * scale if scale > 0.0 else 0.0
    minimum = float(eigenvalues.min(initial=0.0))
    return symmetric, bool(minimum >= floor), minimum


def _unfold_toys(toys, dimension: str, iterations: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """Unfolded toys and refold residuals of (window, measured) pairs.

    A rejected toy (None), a failed unfolding and a non-finite result are dropped.
    """

    unfolded: list[np.ndarray] = []
    residuals: list[np.ndarray] = []
    for toy in toys:
        if toy is None:
            continue
        window, measured = toy
        try:
            value, refolded = unfold(measured, window, dimension, iterations)
        except ValueError:
            continue
        if np.all(np.isfinite(value)):
            unfolded.append(value)
            residuals.append(np.asarray(measured, dtype=float) - refolded)
    return np.asarray(unfolded), residuals


@dataclass
class IterationScan:
    """Every row of one scan, the selected k (None when no candidate passes) and the per-k spectra."""

    system: str
    dimension: str
    measured: np.ndarray
    rows: list[dict[str, Any]]
    selected_k: int | None
    unfolded: dict[int, np.ndarray]
    refolded: dict[int, np.ndarray]
    data_toy_unfolded: dict[int, np.ndarray]

    def row(self, iterations: int) -> dict[str, Any]:
        return next(row for row in self.rows if row["k"] == iterations)


def scan_iterations(
    window,
    dimension: str,
    measured: Sequence[float],
    data_toys: np.ndarray,
    response_toys: Sequence[tuple[Any, np.ndarray] | None],
    *,
    reported_xj: tuple[float, float],
    data_toys_requested: int | None = None,
    system: str = "",
    reference_iterations: int = REFERENCE_ITERATIONS,
) -> IterationScan:
    """Scan k = 1..12 of one problem with the frozen gates and select k.

    ``data_toys`` are the accepted DATA toys of the measured vector (the central response unfolds them);
    ``data_toys_requested`` counts rejected draws as failures. ``response_toys`` holds one (window, measured)
    pair per requested response toy, or None where the toy was rejected; every entry counts in the denominator
    of the response toy success.
    """

    matrix, misses, prior, fakes = problem_arrays(window, dimension)
    measured = np.asarray(measured, dtype=float)
    data_toys = np.asarray(data_toys, dtype=float)
    if response_toys is None or len(response_toys) < 2:
        raise ValueError("iteration selection requires at least two response/leakage/Tcomb toys")
    mask = reported_bins(window, dimension, reported_xj)
    zero_efficiency = zero_efficiency_states(matrix, misses, mask)
    condition = condition_number(matrix, misses)
    pseudo_measured = matrix.sum(axis=0) + fakes
    pseudo_variance = np.clip(pseudo_measured, 1e-12, None)
    unfolded_by_k: dict[int, np.ndarray] = {}
    refolded_by_k: dict[int, np.ndarray] = {}
    toys_by_k: dict[int, np.ndarray] = {}
    rows: list[dict[str, Any]] = []
    for k in COMPUTED_ITERATIONS:
        central, refolded = unfold(measured, window, dimension, k)
        unfolded_by_k[k], refolded_by_k[k] = central, refolded
        toy_array, data_residuals = _unfold_toys(((window, toy) for toy in data_toys), dimension, k)
        toys_by_k[k] = toy_array
        success = len(toy_array) / max(1, data_toys_requested or len(data_toys))
        stat = relative_half_width(central, toy_array, mask)
        response_array, response_residuals = _unfold_toys(response_toys, dimension, k)
        response_success = len(response_array) / len(response_toys)
        response_stat = relative_half_width(central, response_array, mask)
        deviation = math.nan if k == 1 else iteration_deviation(central, unfolded_by_k[k - 1], mask)
        if k == 1:
            objective = math.nan
        else:
            objective = float(math.sqrt(stat * stat + response_stat * response_stat + deviation * deviation))
        symmetric, psd, min_eigen = covariance_status(toy_array[:, mask] if len(toy_array) else toy_array)

        # DATA and response residual-toy covariances are independent and add; each ensemble is centred by its
        # own covariance so a mean offset never counts as spread.
        ensembles = [np.asarray(data_residuals, dtype=float)]
        if response_residuals:
            ensembles.append(np.asarray(response_residuals, dtype=float))
        covariance = np.zeros((measured.size, measured.size), dtype=float)
        response_covariance = False
        for index, ensemble in enumerate(ensembles):
            if ensemble.ndim == 2 and len(ensemble) >= 2:
                component = np.atleast_2d(np.cov(ensemble, rowvar=False, ddof=1))
                covariance += component
                if index > 0 and float(np.max(np.abs(component), initial=0.0)) > 0.0:
                    response_covariance = True
        if response_covariance:
            refold, rank, null_norm = covariance_chi2_ndf_from_matrix(measured - refolded, covariance)
            components = REFOLD_COMPONENTS
        else:
            refold, rank, null_norm, components = math.inf, 0, math.inf, REFOLD_COMPONENTS_INCOMPLETE

        closure, closure_refolded = unfold(pseudo_measured, window, dimension, k)
        mc_closure = chi2_ndf(prior, closure, np.clip(prior + np.abs(closure), 1e-12, None))
        mc_refold = chi2_ndf(pseudo_measured, closure_refolded, pseudo_variance)
        negative_fraction = float(np.abs(central[mask & (central < 0)]).sum() / max(np.abs(central[mask]).sum(), 1e-15))
        gates = {
            "finite_unfolded": bool(np.all(np.isfinite(central))),
            "covariance_symmetric": symmetric,
            "covariance_psd": psd,
            "reported_zero_efficiency_states": zero_efficiency == 0,
            "data_toy_success": success >= TOY_SUCCESS_MIN,
            "response_toy_success": response_success >= TOY_SUCCESS_MIN,
            "response_residual_covariance": response_covariance,
            "data_refold": refold < REFOLD_CHI2_NDF_MAX,
            "mc_closure": mc_closure < MC_CHI2_NDF_MAX,
            "mc_refold": mc_refold < MC_CHI2_NDF_MAX,
            "negative_fraction": negative_fraction < NEGATIVE_FRACTION_MAX,
            "conditioning": condition is not None and condition < CONDITION_MAX,
        }
        failures = [name for name, passed in gates.items() if not passed]
        rows.append({
            "system": system,
            "dimension": dimension,
            "k": k,
            "data_stat_term": stat,
            "response_stat_term": response_stat,
            "response_stat_status": RESPONSE_STAT_STATUS,
            "response_toy_success": response_success,
            "iteration_deviation": None if k == 1 else deviation,
            "combined_objective": None if k == 1 else objective,
            "stat_total": float(math.sqrt(stat * stat + response_stat * response_stat)),
            "photon_data_refold": refold if dimension == "1D" else None,
            "joint_data_refold": refold if dimension == "2D" else None,
            "refold_covariance_rank": rank,
            "refold_null_residual_norm": null_norm,
            "refold_covariance_components": components,
            "mc_closure": mc_closure,
            "mc_refold": mc_refold,
            "negative_fraction": negative_fraction,
            "toy_success": success,
            "covariance_symmetric": symmetric,
            "covariance_psd": psd,
            "covariance_min_eigenvalue": min_eigen,
            "zero_efficiency_reported_states": zero_efficiency,
            "condition_nonzero_singular": condition,
            "passes_gate": bool(k in CANDIDATE_ITERATIONS and not failures),
            "passes_all_non_k_gates": not failures,
            "non_k_gate_failures": "|".join(failures),
            "iteration_selection_eligible": k in CANDIDATE_ITERATIONS,
            "reference_selected": k == reference_iterations,
            "candidate_selected": False,
        })
    eligible = [row for row in rows if row["passes_gate"] and row["combined_objective"] is not None]
    selected_k = None
    if eligible:
        chosen = min(eligible, key=lambda row: (float(row["combined_objective"]), int(row["k"])))
        chosen["candidate_selected"] = True
        selected_k = int(chosen["k"])
    return IterationScan(system, dimension, measured, rows, selected_k, unfolded_by_k, refolded_by_k, toys_by_k)


def gate_contract() -> dict[str, Any]:
    return {
        "candidate_k": list(CANDIDATE_ITERATIONS),
        "computed_k": list(COMPUTED_ITERATIONS),
        "finite_unfolded": True,
        "covariance_symmetric": True,
        "covariance_psd_tolerance": "lambda_min >= -1e-10 * max_abs_eigenvalue (zero floor for an exactly zero covariance)",
        "reported_zero_efficiency_states": 0,
        "toy_success_fraction_min": TOY_SUCCESS_MIN,
        "response_toys_required_for_iteration_selection": True,
        "response_toy_success_fraction_min": TOY_SUCCESS_MIN,
        "response_residual_covariance_required": True,
        "data_refold_chi2_ndf_max": REFOLD_CHI2_NDF_MAX,
        "data_refold_metric": "generalized chi2/rank of central measured-minus-refolded residual using DATA residual-toy "
                              "covariance plus response/leakage/Tcomb residual-toy covariance; eigenmodes <=1e-10 of "
                              "lambda_max excluded with fail-closed null-residual test",
        "mc_truth_closure_chi2_ndf_max": MC_CHI2_NDF_MAX,
        "mc_refold_chi2_ndf_max": MC_CHI2_NDF_MAX,
        "negative_unfolded_fraction_max": NEGATIVE_FRACTION_MAX,
        "condition_nonzero_singular_max": CONDITION_MAX,
        "selection": "argmin Q among passing candidates; deterministic lower-k tie break",
        "high_k_policy": "k=5..12 are diagnostic-only; if no k=2..4 passes, stop on upstream/model closure rather than "
                         "using high k to absorb it",
    }


# ----------------------------------------------------------------------------
# Photon-yield refold test against the measured covariance
# ----------------------------------------------------------------------------

def measured_covariance_refold_chi2(residual: np.ndarray, data_toys: np.ndarray,
                                    response_toys: np.ndarray | None = None) -> float:
    """Refold residual against the covariance of the measured photon yields (DATA toys, plus response toys).

    D'Agostini conserves the measured yield, so the residual of the two-bin photon problem has one free
    component and this chi2 one degree of freedom. It asks whether the refolded yields agree with the measured
    ones within the measurement uncertainty; the residual-toy form of the frozen scan instead compares with the
    spread of the toys' own distance from convergence, which carries no information for two bins.
    """

    covariance = np.cov(np.asarray(data_toys, dtype=float), rowvar=False, ddof=1)
    if response_toys is not None:
        covariance = covariance + np.cov(np.asarray(response_toys, dtype=float), rowvar=False, ddof=1)
    return float(residual @ np.linalg.solve(covariance, residual))


def select_with_refold_statistic(scan: IterationScan, statistic: Mapping[int, float]) -> dict[str, Any]:
    """The frozen selection with the frozen refold statistic replaced by ``statistic``.

    Candidates, threshold, every other gate and the Q_k ordering with its lower-k tie break stay those of the scan.
    """

    eligible = []
    for row in scan.rows:
        others = [name for name in row["non_k_gate_failures"].split("|") if name and name != "data_refold"]
        if row["k"] in CANDIDATE_ITERATIONS and not others and statistic[row["k"]] < REFOLD_CHI2_NDF_MAX:
            eligible.append(row)
    chosen = min(eligible, key=lambda row: (row["combined_objective"], row["k"])) if eligible else None
    return {"eligible_k": [row["k"] for row in eligible], "selected_k": chosen["k"] if chosen else None}


def photon_refold_selection(scan: IterationScan, data_toys: np.ndarray, response_toys: Sequence[np.ndarray]) -> dict[str, Any]:
    """Select the photon-yield k with the measured-covariance refold test (DATA plus response covariance).

    ``data_toys`` are the accepted DATA toys of the photon yields, ``response_toys`` the photon yields of every
    accepted response toy. Returns the per-k residual and statistics and the selection.
    """

    if scan.dimension != "1D":
        raise ValueError("the measured-covariance refold test is defined for the photon-yield problem")
    response = np.asarray(response_toys, dtype=float)
    residual = {k: scan.measured - refolded for k, refolded in scan.refolded.items()}
    data_only = {k: measured_covariance_refold_chi2(r, data_toys) for k, r in residual.items()}
    combined = {k: measured_covariance_refold_chi2(r, data_toys, response) for k, r in residual.items()}
    return {"refold_residual": residual, "chi2_ndf1_data": data_only, "chi2_ndf1_data_plus_response": combined,
            **select_with_refold_statistic(scan, combined)}


__all__ = [
    "CANDIDATE_ITERATIONS", "COMPUTED_ITERATIONS", "CONDITION_MAX", "IterationScan", "MC_CHI2_NDF_MAX",
    "NEGATIVE_FRACTION_MAX", "REFOLD_CHI2_NDF_MAX", "TOY_SUCCESS_MIN", "condition_number", "covariance_status",
    "gate_contract", "iteration_deviation", "measured_covariance_refold_chi2", "photon_refold_selection",
    "problem_arrays", "relative_half_width", "reported_bins", "scan_iterations", "select_with_refold_statistic",
    "unfold", "zero_efficiency_states",
]
