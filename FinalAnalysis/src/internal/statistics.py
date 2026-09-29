"""The statistical summaries of the measurement: the central 68% interval, chi2/ndf with a diagonal or a full
covariance, and the mask of the reported xJ bins.

Ported from the validated collaboration implementation (``chi2_ndf``, ``covariance_chi2_ndf`` of
``photonjet/analysis/unfolding.py`` at main a7f15eb). The definitions are unchanged.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np


def chi2_ndf(observed: np.ndarray, expected: np.ndarray, variance: np.ndarray) -> float:
    observed, expected, variance = (np.asarray(x, dtype=float) for x in (observed, expected, variance))
    mask = np.isfinite(observed) & np.isfinite(expected) & (variance > 0)
    if not np.any(mask):
        return 0.0 if np.allclose(observed, expected, rtol=1e-9, atol=1e-9) else math.inf
    return float(np.sum((observed[mask] - expected[mask]) ** 2 / variance[mask]) / int(mask.sum()))


def covariance_chi2_ndf_from_matrix(residual: np.ndarray, covariance: np.ndarray, *, eigen_rtol: float = 1.0e-10) -> tuple[float, int, float]:
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
    return covariance_chi2_ndf_from_matrix(vector, covariance, eigen_rtol=eigen_rtol)


def central_68(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lower = np.quantile(samples, 0.16, axis=0)
    upper = np.quantile(samples, 0.84, axis=0)
    return lower, upper, 0.5 * (upper - lower)


def reported_mask(xj_edges: np.ndarray, low: float, high: float) -> np.ndarray:
    edges = np.asarray(xj_edges, dtype=float)
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError("invalid reported support")
    if not all(np.any(np.isclose(edges, edge, rtol=0, atol=1e-9)) for edge in (low, high)):
        raise ValueError("reported boundaries must coincide with retained bin edges")
    return (edges[:-1] >= low - 1e-9) & (edges[1:] <= high + 1e-9)


__all__ = ["central_68", "chi2_ndf", "covariance_chi2_ndf", "covariance_chi2_ndf_from_matrix", "reported_mask"]
