"""Iterative-Bayes unfolding kernel and chi2 helpers (iteration selection lives in chain.py)."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

import numpy as np

from .response import REPORTED_XJ_HIGH, REPORTED_XJ_LOW, ResponseBundle




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


__all__ = ["iterative_bayes", "chi2_ndf", "covariance_chi2_ndf"]
