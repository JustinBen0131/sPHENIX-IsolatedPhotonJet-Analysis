"""Reusable weighting helpers with explicit normalization and caps."""

from __future__ import annotations

from typing import Iterable

import numpy as np


def normalize_mean_one(values: Iterable[float]) -> np.ndarray:
    weights = np.asarray(tuple(values), dtype=float)
    if weights.ndim != 1 or not len(weights) or np.any(~np.isfinite(weights)):
        raise ValueError("finite one-dimensional weights are required")
    mean = float(weights.mean())
    if mean <= 0:
        raise ValueError("weight mean must be positive")
    return weights / mean


def inverse_density_weights(
    values: Iterable[float],
    labels: Iterable[int],
    *,
    bins: Iterable[float],
    cap: float | None = None,
) -> np.ndarray:
    x = np.asarray(tuple(values), dtype=float)
    y = np.asarray(tuple(labels), dtype=int)
    edges = np.asarray(tuple(bins), dtype=float)
    if x.shape != y.shape or x.ndim != 1 or np.any(~np.isfinite(x)):
        raise ValueError("values and labels must be aligned and finite")
    if len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("bin edges must be strictly increasing")
    result = np.zeros(len(x), dtype=float)
    indices = np.searchsorted(edges, x, side="right") - 1
    valid = (indices >= 0) & (indices < len(edges) - 1)
    for label in np.unique(y):
        selected = valid & (y == label)
        counts = np.bincount(indices[selected], minlength=len(edges) - 1)
        local = np.where(selected, counts[np.clip(indices, 0, len(counts) - 1)], 0)
        result[selected] = np.divide(
            1.0,
            local[selected],
            out=np.zeros(int(selected.sum()), dtype=float),
            where=local[selected] > 0,
        )
    if cap is not None:
        if cap <= 0:
            raise ValueError("cap must be positive")
        result = np.minimum(result, float(cap))
    return normalize_mean_one(result)
