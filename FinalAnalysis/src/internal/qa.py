"""Numerical photon QA inputs prepared once, independently of figure drawing.

Purity uses the measurement's ABCD solver. Thresholds describe the saved signal
distribution; stitching and centrality fits are descriptive, never selections.
"""
import hashlib
import json

import numpy as np

from . import corrections


def identity(kind, *arrays):
    h = hashlib.sha256(kind.encode())
    for value in arrays:
        a = np.ascontiguousarray(value, dtype=float)
        h.update(json.dumps(a.shape).encode()); h.update(a.tobytes())
    return kind + ':' + h.hexdigest()


def purity_rows(counts, leakage):
    raw, corrected = [], []
    for k in range(counts.shape[1]):
        n = np.maximum(counts[:, k], 0); l = leakage[:, k]
        raw.append(corrections.solve_leakage_abcd(n, [1, 0, 0, 0]).purity if n[0] > 0 else np.nan)
        corrected.append(corrections.solve_leakage_abcd(n, l).purity if n[0] > 0 and l[0] > 0 else np.nan)
    return np.asarray([raw, corrected])


def quantiles(values, edges, efficiencies, above):
    total = values.sum()
    if total <= 0: return np.full(len(efficiencies), np.nan)
    cdf = np.r_[0, np.cumsum(values)] / total
    return np.interp(1 - np.asarray(efficiencies) if above else efficiencies, cdf, edges)


def stitching_fit(x, y, error):
    good = (y > 0) & np.isfinite(y) & (error > 0)
    if good.sum() < 5: return np.full(4, np.nan)
    lx = np.log(x[good]); weight = (y[good] / error[good]) ** 2
    design = np.column_stack([np.ones(good.sum()), -lx, -lx * lx, -x[good] * lx])
    return np.linalg.lstsq(design * np.sqrt(weight)[:, None], np.log(y[good]) * np.sqrt(weight), rcond=None)[0]


def centrality_fit(x, y):
    good = np.isfinite(y)
    return np.polyfit(x[good], y[good], 1) if good.sum() >= 2 else np.full(2, np.nan)
