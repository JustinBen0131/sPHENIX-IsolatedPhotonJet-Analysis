"""Frozen iteration-selection procedure; numerical k is configuration-specific."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

import numpy as np

from .response import REPORTED_XJ_HIGH, REPORTED_XJ_LOW, ResponseBundle


K_COMPUTE_RANGE = tuple(range(1, 13))
K_CANDIDATE_RANGE = tuple(range(2, 5))
TOY_SUCCESS_GATE = 0.95
REFOLD_CHI2_NDF_GATE = 1.5
MC_CHI2_NDF_GATE = 5.0
NEGATIVE_FRACTION_GATE = 0.50
CONDITION_GATE = 1.0e8
HIGH_K_SCIENCE_REVIEW_THRESHOLD = 4


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
    if iterations < 1:
        raise ValueError("iterations must be positive")
    totals = matrix.sum(axis=1) + miss
    conditional = np.divide(matrix, totals[:, None], out=np.zeros_like(matrix), where=totals[:, None] > 0)
    efficiency = conditional.sum(axis=1)
    state_counts = totals.copy() if prior is None else np.asarray(prior, dtype=float).copy()
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
        posterior = np.divide(
            conditional * state[:, None], denominator[None, :],
            out=np.zeros_like(conditional), where=denominator[None, :] > 0,
        )
        unfolded = np.divide(
            posterior @ measured_array, efficiency,
            out=np.zeros(conditional.shape[0]), where=efficiency > 0,
        )
        total = float(unfolded.sum())
        if not math.isfinite(total) or abs(total) <= 1e-15:
            raise ValueError("empty signed estimate")
        state = unfolded / total
    return unfolded[:ntruth], conditional.T @ unfolded


def chi2_ndf(observed: np.ndarray, expected: np.ndarray, variance: np.ndarray) -> float:
    observed, expected, variance = map(lambda x: np.asarray(x, dtype=float), (observed, expected, variance))
    mask = np.isfinite(observed) & np.isfinite(expected) & (variance > 0)
    if not np.any(mask):
        return 0.0 if np.allclose(observed, expected, rtol=1e-9, atol=1e-9) else math.inf
    return float(np.sum((observed[mask] - expected[mask]) ** 2 / variance[mask]) / int(mask.sum()))


def _covariance_chi2_ndf_from_matrix(
    residual: np.ndarray,
    covariance: np.ndarray,
    *,
    eigen_rtol: float = 1.0e-10,
) -> tuple[float, int, float]:
    """Generalized residual chi2 using a supported covariance subspace."""
    vector = np.asarray(residual, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    covariance = 0.5 * (covariance + covariance.T)
    if covariance.shape != (vector.size, vector.size) or np.any(~np.isfinite(covariance)):
        return math.inf, 0, math.inf
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    maximum = float(np.max(eigenvalues, initial=0.0))
    if maximum <= 0.0:
        return (0.0, 0, 0.0) if np.allclose(vector, 0.0, rtol=1e-9, atol=1e-9) else (math.inf, 0, math.inf)
    keep = eigenvalues > eigen_rtol * maximum
    rank = int(np.sum(keep))
    if rank == 0:
        return math.inf, 0, math.inf
    coordinates = eigenvectors.T @ vector
    null_norm = float(np.linalg.norm(coordinates[~keep]))
    reference = max(float(np.linalg.norm(vector)), 1.0)
    if null_norm > 1.0e-9 * reference:
        return math.inf, rank, null_norm
    value = float(np.sum(coordinates[keep] ** 2 / eigenvalues[keep]) / rank)
    return value, rank, null_norm


def covariance_chi2_ndf(
    residual: np.ndarray,
    residual_toys: Sequence[np.ndarray],
    *,
    eigen_rtol: float = 1.0e-10,
) -> tuple[float, int, float]:
    """Generalized residual chi2 using the supported toy-covariance subspace.

    Corrected recoil bins are signed and can contain zero observed counts, so
    an observed-count diagonal variance incorrectly assigns zero uncertainty
    after ABCD/K subtraction.  The residual-toy covariance preserves DATA and
    response/leakage/Tcomb correlations.  Numerically null modes are excluded
    by a fixed relative eigenvalue threshold; a residual component outside the
    supported subspace fails closed.
    """
    vector = np.asarray(residual, dtype=float)
    toys = np.asarray(residual_toys, dtype=float)
    if toys.ndim != 2 or toys.shape[1:] != vector.shape or len(toys) < 2:
        return math.inf, 0, math.inf
    covariance = np.atleast_2d(np.cov(toys, rowvar=False, ddof=1))
    return _covariance_chi2_ndf_from_matrix(
        vector, covariance, eigen_rtol=eigen_rtol
    )


def reported_mask(bundle: ResponseBundle) -> np.ndarray:
    if bundle.dimension == "1D":
        return np.ones(bundle.truth.size, dtype=bool)
    assert bundle.xj_edges is not None
    xmask = (bundle.xj_edges[:-1] >= REPORTED_XJ_LOW - 1e-12) & (bundle.xj_edges[1:] <= REPORTED_XJ_HIGH + 1e-12)
    return np.tile(xmask, len(bundle.truth_ptgamma_edges) - 1)


def relative_bin_values(numerator: np.ndarray, denominator_a: np.ndarray, denominator_b: np.ndarray,
                        mask: np.ndarray) -> np.ndarray:
    scale = np.maximum(0.5 * (np.abs(denominator_a) + np.abs(denominator_b)), 1e-12)
    valid = mask & np.isfinite(numerator) & np.isfinite(scale) & (denominator_a > 0) & (denominator_b > 0)
    return np.abs(numerator[valid]) / scale[valid]


def iteration_deviation(current: np.ndarray, previous: np.ndarray, mask: np.ndarray) -> float:
    """Median symmetric relative movement over reported bins positive at k,k-1."""
    values = relative_bin_values(current - previous, current, previous, mask)
    return float(np.median(values)) if values.size else math.inf


def data_stat_term(central: np.ndarray, toy_unfolded: np.ndarray, mask: np.ndarray) -> float:
    """Median central-68% relative half-width over positive reported bins."""
    if toy_unfolded.ndim != 2 or not len(toy_unfolded):
        return math.inf
    low, high = np.quantile(toy_unfolded, [0.16, 0.84], axis=0)
    half = 0.5 * (high - low)
    valid = mask & np.isfinite(central) & np.isfinite(half) & (central > 0)
    values = half[valid] / np.maximum(np.abs(central[valid]), 1e-12)
    return float(np.median(values)) if values.size else math.inf


def covariance_status(toys: np.ndarray) -> tuple[bool, bool, float | None]:
    if toys.ndim != 2 or len(toys) < 2 or np.any(~np.isfinite(toys)):
        return False, False, None
    cov = np.cov(toys, rowvar=False, ddof=1)
    cov = np.atleast_2d(cov)
    symmetric = bool(np.allclose(cov, cov.T, rtol=1e-10, atol=1e-12))
    eigen = np.linalg.eigvalsh(0.5 * (cov + cov.T))
    scale = float(np.max(np.abs(eigen), initial=0.0))
    floor = -1e-10 * scale if scale > 0.0 else 0.0
    return symmetric, bool(float(eigen.min(initial=0.0)) >= floor), float(eigen.min(initial=0.0))


@dataclass
class ScanResult:
    rows: list[dict[str, Any]]
    selected_k: int
    selected_unfolded: np.ndarray
    selected_refolded: np.ndarray
    selected_toys: np.ndarray


class NoPassingIterationError(RuntimeError):
    """Frozen-gate failure carrying the complete diagnostic k scan."""

    def __init__(self, message: str, rows: list[dict[str, Any]]):
        super().__init__(message)
        self.rows = rows


def scan_problem(
    bundle: ResponseBundle,
    measured: np.ndarray,
    variance: np.ndarray,
    measured_toys: np.ndarray,
    *,
    reference_selected_k: int = 12,
    measured_toy_requested: int | None = None,
    response_toys: Sequence[tuple[ResponseBundle, np.ndarray]] | None = None,
) -> ScanResult:
    measured = np.asarray(measured, dtype=float)
    variance = np.asarray(variance, dtype=float)
    measured_toys = np.asarray(measured_toys, dtype=float)
    if response_toys is None or len(response_toys) < 2:
        raise ValueError(
            "iteration selection requires at least two response/leakage/Tcomb toys; "
            "the iterative_bayes kernel remains available for diagnostic scans"
        )
    mask = reported_mask(bundle)
    totals = bundle.matrix.sum(axis=1) + bundle.all_misses
    efficiency = np.divide(bundle.matrix.sum(axis=1), totals, out=np.zeros_like(totals), where=totals > 0)
    zero_eff_reported = int(np.sum(mask & (totals > 0) & (efficiency == 0)))
    condition = bundle.diagnostics()["condition_nonzero_singular"]
    pseudo_measured = bundle.matrix.sum(axis=0) + bundle.unfolding_fakes
    pseudo_variance = np.clip(pseudo_measured, 1e-12, None)
    central_by_k: dict[int, np.ndarray] = {}
    refold_by_k: dict[int, np.ndarray] = {}
    toys_by_k: dict[int, np.ndarray] = {}
    rows: list[dict[str, Any]] = []
    for k in K_COMPUTE_RANGE:
        central, refolded = iterative_bayes(
            measured, bundle.matrix, bundle.all_misses, k, bundle.truth, fakes=bundle.unfolding_fakes
        )
        central_by_k[k], refold_by_k[k] = central, refolded
        unfolded_toys: list[np.ndarray] = []
        data_refold_residual_toys: list[np.ndarray] = []
        for toy in measured_toys:
            try:
                value, toy_refolded = iterative_bayes(
                    toy, bundle.matrix, bundle.all_misses, k, bundle.truth, fakes=bundle.unfolding_fakes
                )
            except ValueError:
                continue
            if np.all(np.isfinite(value)):
                unfolded_toys.append(value)
                data_refold_residual_toys.append(np.asarray(toy, dtype=float) - toy_refolded)
        toy_array = np.asarray(unfolded_toys)
        toys_by_k[k] = toy_array
        success = len(toy_array) / max(1, measured_toy_requested or len(measured_toys))
        stat = data_stat_term(central, toy_array, mask)
        response_unfolded: list[np.ndarray] = []
        response_refold_residual_toys: list[np.ndarray] = []
        if response_toys is not None:
            for toy_bundle, toy_measured in response_toys:
                try:
                    value, toy_refolded = iterative_bayes(
                        toy_measured, toy_bundle.matrix, toy_bundle.all_misses, k,
                        toy_bundle.truth, fakes=toy_bundle.unfolding_fakes,
                    )
                except ValueError:
                    continue
                if np.all(np.isfinite(value)):
                    response_unfolded.append(value)
                    response_refold_residual_toys.append(
                        np.asarray(toy_measured, dtype=float) - toy_refolded
                    )
        response_array = np.asarray(response_unfolded)
        response_success = len(response_array) / len(response_toys)
        response_stat = data_stat_term(central, response_array, mask)
        deviation = math.nan if k == 1 else iteration_deviation(central, central_by_k[k - 1], mask)
        if k == 1:
            objective = math.nan
        else:
            objective = float(math.sqrt(
                stat * stat + response_stat * response_stat + deviation * deviation
            ))
        symmetric, psd, min_eigen = covariance_status(toy_array[:, mask] if len(toy_array) else toy_array)
        residual_ensembles = [np.asarray(data_refold_residual_toys, dtype=float)]
        if response_refold_residual_toys:
            # Independent DATA-statistical and source-bootstrap response terms
            # add in covariance.  Center each ensemble before concatenation so
            # their means cannot masquerade as statistical variance.
            residual_ensembles.append(np.asarray(response_refold_residual_toys, dtype=float))
        covariance = np.zeros((measured.size, measured.size), dtype=float)
        response_residual_covariance_available = False
        for ensemble_index, ensemble in enumerate(residual_ensembles):
            if ensemble.ndim == 2 and len(ensemble) >= 2:
                component = np.atleast_2d(np.cov(ensemble, rowvar=False, ddof=1))
                covariance += component
                if ensemble_index > 0 and float(np.max(np.abs(component), initial=0.0)) > 0.0:
                    response_residual_covariance_available = True
        if not response_residual_covariance_available:
            data_refold = math.inf
            refold_covariance_rank = 0
            refold_null_residual_norm = math.inf
            refold_covariance_components = "INCOMPLETE_RESPONSE_RESIDUAL_COVARIANCE"
        else:
            data_refold, refold_covariance_rank, refold_null_residual_norm = _covariance_chi2_ndf_from_matrix(
                measured - refolded, covariance,
            )
            refold_covariance_components = (
                "DATA_RESIDUAL_TOYS_PLUS_RESPONSE_LEAKAGE_TCOMB_RESIDUAL_TOYS"
            )
        closure, closure_refold = iterative_bayes(
            pseudo_measured, bundle.matrix, bundle.all_misses, k, bundle.truth, fakes=bundle.unfolding_fakes
        )
        mc_closure = chi2_ndf(bundle.truth, closure, np.clip(bundle.truth + np.abs(closure), 1e-12, None))
        mc_refold = chi2_ndf(pseudo_measured, closure_refold, pseudo_variance)
        negative_fraction = float(np.abs(central[mask & (central < 0)]).sum() / max(np.abs(central[mask]).sum(), 1e-15))
        non_k_gates = {
            "finite_unfolded": bool(np.all(np.isfinite(central))),
            "covariance_symmetric": symmetric,
            "covariance_psd": psd,
            "reported_zero_efficiency_states": zero_eff_reported == 0,
            "data_toy_success": success >= TOY_SUCCESS_GATE,
            "response_toy_success": response_success >= 0.95,
            "response_residual_covariance": response_residual_covariance_available,
            "data_refold": data_refold < REFOLD_CHI2_NDF_GATE,
            "mc_closure": mc_closure < MC_CHI2_NDF_GATE,
            "mc_refold": mc_refold < MC_CHI2_NDF_GATE,
            "negative_fraction": negative_fraction < NEGATIVE_FRACTION_GATE,
            "conditioning": condition is not None and condition < CONDITION_GATE,
        }
        non_k_failures = [name for name, passed in non_k_gates.items() if not passed]
        passes_non_k_gates = not non_k_failures
        passes = bool(k in K_CANDIDATE_RANGE and passes_non_k_gates)
        rows.append({
            "system": bundle.system,
            "dimension": bundle.dimension,
            "k": k,
            "data_stat_term": stat,
            "response_stat_term": response_stat,
            "response_stat_status": "SOURCE_OCCURRENCE_BOOTSTRAP",
            "response_toy_success": response_success,
            "iteration_deviation": None if k == 1 else deviation,
            "combined_objective": None if k == 1 else objective,
            "photon_data_refold": data_refold if bundle.dimension == "1D" else None,
            "joint_data_refold": data_refold if bundle.dimension == "2D" else None,
            "refold_covariance_rank": refold_covariance_rank,
            "refold_null_residual_norm": refold_null_residual_norm,
            "refold_covariance_components": refold_covariance_components,
            "mc_closure": mc_closure,
            "mc_refold": mc_refold,
            "negative_fraction": negative_fraction,
            "toy_success": success,
            "covariance_symmetric": symmetric,
            "covariance_psd": psd,
            "covariance_min_eigenvalue": min_eigen,
            "zero_efficiency_reported_states": zero_eff_reported,
            "condition_nonzero_singular": condition,
            "passes_gate": passes,
            "passes_all_non_k_gates": passes_non_k_gates,
            "non_k_gate_failures": "|".join(non_k_failures),
            "iteration_selection_eligible": k in K_CANDIDATE_RANGE,
            "reference_selected": k == reference_selected_k,
            "candidate_selected": False,
        })
    eligible = [row for row in rows if row["passes_gate"] and row["combined_objective"] is not None]
    if not eligible:
        raise NoPassingIterationError(
            f"no iteration passes frozen gates: {bundle.system} {bundle.dimension}",
            rows,
        )
    selected = min(eligible, key=lambda row: (float(row["combined_objective"]), int(row["k"])))
    selected["candidate_selected"] = True
    kstar = int(selected["k"])
    return ScanResult(rows, kstar, central_by_k[kstar], refold_by_k[kstar], toys_by_k[kstar])


def gate_contract() -> dict[str, Any]:
    return {
        "candidate_k": list(K_CANDIDATE_RANGE),
        "computed_k": list(K_COMPUTE_RANGE),
        "finite_unfolded": True,
        "covariance_symmetric": True,
        "covariance_psd_tolerance": "lambda_min >= -1e-10 * max_abs_eigenvalue (zero floor for an exactly zero covariance)",
        "reported_zero_efficiency_states": 0,
        "toy_success_fraction_min": TOY_SUCCESS_GATE,
        "response_toys_required_for_iteration_selection": True,
        "response_toy_success_fraction_min": 0.95,
        "response_residual_covariance_required": True,
        "data_refold_chi2_ndf_max": REFOLD_CHI2_NDF_GATE,
        "data_refold_metric": "generalized chi2/rank of central measured-minus-refolded residual using DATA residual-toy covariance plus response/leakage/Tcomb residual-toy covariance; eigenmodes <=1e-10 of lambda_max excluded with fail-closed null-residual test",
        "mc_truth_closure_chi2_ndf_max": MC_CHI2_NDF_GATE,
        "mc_refold_chi2_ndf_max": MC_CHI2_NDF_GATE,
        "negative_unfolded_fraction_max": NEGATIVE_FRACTION_GATE,
        "condition_nonzero_singular_max": CONDITION_GATE,
        "selection": "argmin Q among passing candidates; deterministic lower-k tie break",
        "high_k_policy": "k=5..12 are diagnostic-only; if no k=2..4 passes, stop on upstream/model closure rather than using high k to absorb it",
    }
