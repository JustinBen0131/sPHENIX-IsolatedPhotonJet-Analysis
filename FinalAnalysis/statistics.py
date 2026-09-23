"""Uncertainties, closure metrics and the iteration gate.

Ported from the validated collaboration implementation (``chi2_ndf``,
``covariance_chi2_ndf`` of ``photonjet/analysis/unfolding.py`` and
``chain_metrics`` / ``gate_and_score`` of ``photonjet/analysis/chain.py`` at
main a7f15eb). The definitions are unchanged.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np


def chi2_ndf(observed: np.ndarray, expected: np.ndarray, variance: np.ndarray) -> float:
    observed, expected, variance = (np.asarray(x, dtype=float) for x in (observed, expected, variance))
    mask = np.isfinite(observed) & np.isfinite(expected) & (variance > 0)
    if not np.any(mask):
        return 0.0 if np.allclose(observed, expected, rtol=1e-9, atol=1e-9) else math.inf
    return float(np.sum((observed[mask] - expected[mask]) ** 2 / variance[mask]) / int(mask.sum()))


def _covariance_chi2_ndf_from_matrix(residual: np.ndarray, covariance: np.ndarray, *, eigen_rtol: float = 1.0e-10) -> tuple[float, int, float]:
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


def covariance_chi2_ndf(residual: np.ndarray, residual_toys: Sequence[np.ndarray], *, eigen_rtol: float = 1.0e-10) -> tuple[float, int, float]:
    """Generalised residual chi2 on the supported toy-covariance subspace."""

    vector = np.asarray(residual, dtype=float)
    toys = np.asarray(residual_toys, dtype=float)
    if toys.ndim != 2 or toys.shape[1:] != vector.shape or len(toys) < 2:
        return math.inf, 0, math.inf
    covariance = np.atleast_2d(np.cov(toys, rowvar=False, ddof=1))
    return _covariance_chi2_ndf_from_matrix(vector, covariance, eigen_rtol=eigen_rtol)


def central_68(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lower = np.quantile(samples, 0.16, axis=0)
    upper = np.quantile(samples, 0.84, axis=0)
    return lower, upper, 0.5 * (upper - lower)


def toy_covariance(samples: np.ndarray, error: np.ndarray, lower: np.ndarray, upper: np.ndarray) -> np.ndarray:
    """Winsorised toy covariance rescaled to the central-68 errors."""

    n = samples.shape[1]
    winsorized = np.clip(samples, lower, upper)
    covariance = np.atleast_2d(np.cov(winsorized, rowvar=False, ddof=1))
    scale = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    correlation = np.divide(covariance, scale[:, None] * scale[None, :], out=np.eye(n), where=(scale[:, None] * scale[None, :]) > 0)
    return correlation * error[:, None] * error[None, :]


def reported_mask(xj_edges: np.ndarray, low: float, high: float) -> np.ndarray:
    edges = np.asarray(xj_edges, dtype=float)
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError("invalid reported support")
    if not all(np.any(np.isclose(edges, edge, rtol=0, atol=1e-9)) for edge in (low, high)):
        raise ValueError("reported boundaries must coincide with retained bin edges")
    return (edges[:-1] >= low - 1e-9) & (edges[1:] <= high + 1e-9)


def chain_metrics(result: Mapping[str, Any], reported_low: float, reported_high: float) -> dict[str, float | int | bool]:
    """Quality metrics of one chain result."""

    edges = np.asarray(result["xj_edges"], dtype=float)
    density = np.asarray(result["density"], dtype=float)
    error = np.asarray(result["density_error"], dtype=float)
    shown = reported_mask(edges, reported_low, reported_high)
    positive = shown & np.isfinite(density) & np.isfinite(error) & (density > 1.0e-10)
    relative = np.divide(error, density, out=np.full_like(error, np.inf), where=density > 0)
    selected = relative[positive]
    median_relative = float(np.median(selected)) if selected.size else math.inf
    max_relative = float(np.max(selected)) if selected.size else math.inf
    log_density = np.log(np.clip(density[positive], 1.0e-12, None))
    smoothness = float(np.mean(np.abs(np.diff(log_density, n=2)))) if log_density.size >= 3 else math.inf
    finite = bool(np.all(np.isfinite(density)) and np.all(density >= 0) and np.all(np.isfinite(error[shown])) and np.all(error[shown] >= 0))
    covariance = np.asarray(result.get("density_stat_covariance", np.diag(np.square(error))), dtype=float)
    shape_ok = covariance.shape == (density.size, density.size)
    symmetric = bool(shape_ok and np.all(np.isfinite(covariance)) and np.allclose(covariance, covariance.T, rtol=1e-8, atol=1e-12))
    minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(0.5 * (covariance + covariance.T)))) if symmetric else -math.inf
    covariance_scale = float(np.max(np.abs(np.diag(covariance)))) if shape_ok else 0.0
    psd = bool(minimum_eigenvalue >= -1.0e-9 * max(1.0, covariance_scale))
    purity_rows = result.get("purity", [])
    purity_values = np.asarray([row.get("purity", math.nan) for row in purity_rows], dtype=float)
    purity_signals = np.asarray([row.get("signal_a", math.nan) for row in purity_rows], dtype=float)
    residuals = np.asarray([row.get("factorization_residual", math.inf) for row in purity_rows], dtype=float)
    observability = result.get("response_observability", {})
    n_pt = len(result.get("pt_edges", [0, 0, 0, 0])) - 1
    return {
        "finite_nonnegative": finite,
        "stat_covariance_symmetric": symmetric,
        "stat_covariance_psd": psd,
        "stat_covariance_min_eigenvalue": minimum_eigenvalue,
        "purity_bin_count": int(len(purity_rows)),
        "purity_finite_physical": bool(len(purity_rows) == n_pt and np.all(np.isfinite(purity_values))
                                       and np.all((purity_values > 0) & (purity_values <= 1))
                                       and np.all(np.isfinite(purity_signals)) and np.all(purity_signals > 0)),
        "abcd_max_factorization_residual": float(np.max(residuals)) if residuals.size else math.inf,
        "purity_boundary_bins": int(np.sum(purity_values >= 0.9999)) if purity_values.size else -1,
        "xj_zero_efficiency_supported_bins": int(observability.get("xj_zero_efficiency_reported_bins", -1)),
        "photon_zero_efficiency_supported_bins": int(observability.get("photon_zero_efficiency_supported_bins", -1)),
        "shown_nonzero_bins": int(positive.sum()),
        "median_relative_stat": median_relative,
        "max_relative_stat": max_relative,
        "log_curvature": smoothness,
        "refold_chi2_ndf": float(result["refold_chi2_ndf"]),
        "photon_refold_chi2_ndf": float(result["photon_refold_chi2_ndf"]),
        "mc_closure_chi2_ndf": float(result["mc_closure_chi2_ndf"]),
        "mc_refold_chi2_ndf": float(result["mc_refold_chi2_ndf"]),
        "negative_input_fraction": float(result["negative_input_fraction"]),
        "toy_success_fraction": float(result["toy_successes"] / max(1, result["toy_requested"])),
    }


def gate_and_score(result: Mapping[str, Any], reported_low: float, reported_high: float) -> tuple[bool, float, dict[str, Any]]:
    """Acceptance gate and ranking score of one iteration count."""

    m = chain_metrics(result, reported_low, reported_high)
    passed = (
        bool(m["finite_nonnegative"]) and bool(m["stat_covariance_symmetric"]) and bool(m["stat_covariance_psd"])
        and bool(m["purity_finite_physical"]) and float(m["abcd_max_factorization_residual"]) < 0.02
        and int(m["purity_boundary_bins"]) == 0 and int(m["xj_zero_efficiency_supported_bins"]) == 0
        and int(m["photon_zero_efficiency_supported_bins"]) == 0 and int(m["shown_nonzero_bins"]) >= 3
        and float(m["toy_success_fraction"]) >= 0.90
        and float(m["refold_chi2_ndf"]) < 5.0 and float(m["photon_refold_chi2_ndf"]) < 5.0
        and float(m["mc_closure_chi2_ndf"]) < 5.0 and float(m["mc_refold_chi2_ndf"]) < 5.0
        and float(m["median_relative_stat"]) < 0.80 and float(m["max_relative_stat"]) < 2.00
        and float(m["negative_input_fraction"]) < 0.50
    )
    score = 0.0
    for key in ("refold_chi2_ndf", "photon_refold_chi2_ndf", "mc_closure_chi2_ndf", "mc_refold_chi2_ndf"):
        score += math.log1p(max(0.0, min(float(m[key]), 1.0e6)))
    score += 0.40 * min(float(m["median_relative_stat"]), 20.0)
    score += 0.15 * min(float(m["max_relative_stat"]), 50.0)
    score += 0.25 * min(float(m["negative_input_fraction"]), 10.0)
    score += 0.10 * min(float(m["log_curvature"]), 20.0)
    score += 3.0 * max(0.0, 0.90 - float(m["toy_success_fraction"]))
    if not passed:
        score += 1000.0
    return passed, float(score), m


__all__ = ["central_68", "chain_metrics", "chi2_ndf", "covariance_chi2_ndf", "gate_and_score", "reported_mask", "toy_covariance"]
