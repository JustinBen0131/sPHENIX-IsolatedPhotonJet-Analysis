"""Data-driven combinatoric-jet subtraction from the photon-jet Delta-phi distributions (the sPHENIX Au+Au dijet method).

Jets that do not come from the photon's hard scattering (underlying-event fluctuations, jets of another scattering)
populate every Delta phi with a flow-modulated shape. Per region (A, C), photon-pT bin and the shape scope:

  1. shape: the Delta-phi distribution of the jets separated from the photon by |Delta eta| > 0.8, which removes the
     jets of the same hard scattering, is fitted over 0 < Delta phi < 2.5 with the flow form
         f(Delta phi) = a0 + a2 cos(2 Delta phi) + a3 cos(3 Delta phi),   v2,2 = a2 / (2 a0),  v3,3 = a3 / (2 a0)
     by weighted least squares on the bin integrals;
  2. scale: per xJ bin, the shape is normalised to the jets at every Delta eta in the sideband 0.8 < Delta phi < 2.5;
  3. subtract: its integral over the recoil window Delta phi > 7 pi / 8 is the combinatoric count of that bin.

The flat variant (the systematic) replaces the shape by a constant. The shape is fitted per photon-pT bin on the jets
summed over xJ (scope "pt"), or per photon-pT and xJ bin (scope "pt_xj"). Inputs are counts (and their variances) on
the Delta-phi edges FillPhotonJet.C writes, which contain the fit, sideband and window boundaries. The counts enter the
measured statistics as recoil_background, which corrections.correct_measured removes from the region A and C recoil
spectra before the purity correction. simulation_closure checks the estimate on simulation, where the combinatoric jets
(no truth jet nearby) are known.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

# pi/32 steps with 0.8 and 2.5 rad added (TreeToHists FillPhotonJet.C deltaPhiEdges).
DPHI_EDGES = np.unique(np.r_[np.arange(33) * math.pi / 32.0, 0.8, 2.5])
FIT_RANGE = (0.0, 2.5)
SIDEBAND = (0.8, 2.5)
WINDOW = (7.0 * math.pi / 8.0, math.pi)
SHAPES = ("flow", "flat")
SCOPES = ("pt", "pt_xj")


def bin_mask(edges: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Bins inside [lo, hi]; lo and hi must be bin edges."""
    edges = np.asarray(edges, dtype=float)
    for value in (lo, hi):
        if not np.any(np.isclose(edges, value, rtol=0, atol=1e-12)):
            raise ValueError(f"{value} is not a Delta-phi bin edge")
    return (edges[:-1] >= lo - 1e-12) & (edges[1:] <= hi + 1e-12)


def basis(edges: np.ndarray) -> np.ndarray:
    """(n_bins, 3): the integrals of 1, cos 2x and cos 3x over every bin."""
    lo, hi = np.asarray(edges[:-1], dtype=float), np.asarray(edges[1:], dtype=float)
    return np.column_stack([hi - lo, (np.sin(2 * hi) - np.sin(2 * lo)) / 2.0, (np.sin(3 * hi) - np.sin(3 * lo)) / 3.0])


def integral(parameters: np.ndarray, lo: float, hi: float) -> float:
    a0, a2, a3 = parameters
    return float(a0 * (hi - lo) + a2 * (math.sin(2 * hi) - math.sin(2 * lo)) / 2.0 + a3 * (math.sin(3 * hi) - math.sin(3 * lo)) / 3.0)


def fit_flow(counts: np.ndarray, variances: np.ndarray, edges: np.ndarray = DPHI_EDGES, fit_range=FIT_RANGE) -> dict[str, Any]:
    """Weighted least squares of the flow form on the bin integrals inside the fit range. A bin's weight is 1 / its
    variance, the variance taken as at least one count; empty ranges are refused."""
    counts, variances = np.asarray(counts, dtype=float), np.asarray(variances, dtype=float)
    inside = bin_mask(edges, *fit_range)
    if counts[inside].sum() <= 0:
        raise ValueError("no jets in the Delta-phi fit range")
    design = basis(edges)[inside]
    weight = 1.0 / np.maximum(variances[inside], 1.0)
    normal = design.T @ (weight[:, None] * design)
    covariance = np.linalg.inv(normal)
    parameters = covariance @ (design.T @ (weight * counts[inside]))
    residual = counts[inside] - design @ parameters
    return {"parameters": parameters, "covariance": covariance, "chi2": float(np.sum(weight * residual ** 2)),
            "ndf": int(inside.sum() - 3), "v2": float(parameters[1] / (2 * parameters[0])), "v3": float(parameters[2] / (2 * parameters[0]))}


def window_fraction(parameters: np.ndarray) -> float:
    """Integral of the shape over the recoil window divided by its integral over the sideband."""
    sideband = integral(parameters, *SIDEBAND)
    if not sideband > 0:
        raise ValueError("the fitted shape is not positive in the sideband")
    return integral(parameters, *WINDOW) / sideband


FLAT = np.array([1.0, 0.0, 0.0])


def estimate(dphi: np.ndarray, dphi_variance: np.ndarray, *, shape: str = "flow", scope: str = "pt",
             edges: np.ndarray = DPHI_EDGES, window_occupancy: np.ndarray | None = None) -> dict[str, Any]:
    """Combinatoric counts in the recoil window per region, photon-pT and xJ bin.

    dphi: (2 regions, n_pt, n_xj, n_dphi, 2 sides) counts, side 0 near (|Delta eta| <= 0.8), 1 far; dphi_variance the
    same shape. Returns counts (2, n_pt, n_xj), their variance (sideband statistics through the scale), and per shape
    fit its v2,2, v3,3, chi2 / ndf and window fraction. ``window_occupancy`` (n_pt,), when given, scales the window
    extrapolation of every region (window_occupancy()).
    """
    if shape not in SHAPES or scope not in SCOPES:
        raise ValueError(f"unknown shape {shape} or scope {scope}")
    dphi, dphi_variance = np.asarray(dphi, dtype=float), np.asarray(dphi_variance, dtype=float)
    n_region, n_pt, n_xj, n_dphi, n_side = dphi.shape
    if n_dphi != len(edges) - 1 or n_side != 2:
        raise ValueError("Delta-phi counts do not match the edges and the two sides")
    occupancy = np.ones(n_pt) if window_occupancy is None else np.asarray(window_occupancy, dtype=float)
    if occupancy.shape != (n_pt,) or not np.all(np.isfinite(occupancy) & (occupancy > 0)):
        raise ValueError("the window occupancy needs one finite positive factor per photon-pT bin")
    sideband = bin_mask(edges, *SIDEBAND)
    every = dphi.sum(axis=4)                     # jets at every Delta eta
    every_variance = dphi_variance.sum(axis=4)
    sideband_counts = every[..., sideband].sum(axis=3)           # (region, pt, xj)
    sideband_variance = every_variance[..., sideband].sum(axis=3)
    counts = np.zeros((n_region, n_pt, n_xj))
    variance = np.zeros_like(counts)
    fits = []
    for r in range(n_region):
        for p in range(n_pt):
            cells = [None] if scope == "pt" else range(n_xj)
            for x in cells:
                far = dphi[r, p, :, :, 1].sum(axis=0) if x is None else dphi[r, p, x, :, 1]
                far_variance = dphi_variance[r, p, :, :, 1].sum(axis=0) if x is None else dphi_variance[r, p, x, :, 1]
                columns = slice(None) if x is None else slice(x, x + 1)
                if shape == "flat":
                    parameters, record = FLAT, {"shape": "flat"}
                else:
                    fit = fit_flow(far, far_variance, edges)
                    parameters = fit["parameters"]
                    record = {"shape": "flow", "v2": fit["v2"], "v3": fit["v3"], "chi2": fit["chi2"], "ndf": fit["ndf"]}
                fraction = window_fraction(parameters) * occupancy[p]
                counts[r, p, columns] = fraction * sideband_counts[r, p, columns]
                variance[r, p, columns] = fraction ** 2 * sideband_variance[r, p, columns]
                fits.append({"region": "AC"[r], "pt_bin": p, "xj_bin": x, "window_over_sideband": fraction,
                             "window_occupancy": float(occupancy[p]), **record})
    return {"counts": counts, "variance": variance, "fits": fits, "shape": shape, "scope": scope,
            "window_occupancy": occupancy if window_occupancy is not None else None}


def window_occupancy(combinatoric: np.ndarray, combinatoric_variance: np.ndarray,
                     edges: np.ndarray = DPHI_EDGES) -> dict[str, np.ndarray]:
    """Per photon-pT bin, the density of the simulated combinatoric jets (no truth jet nearby) in the recoil window
    relative to the sideband, per unit Delta phi, from (n_pt, n_xj, n_dphi, 2) counts and their variances.

    The photon's own recoil jets take up part of the narrow window, so fewer separate combinatoric jets form there than in
    the sideband: an estimate that counts every recoil jet per photon scales its window extrapolation by this factor
    (estimate(window_occupancy=)). In simulation the combinatoric jets carry no correlation with the photon, so the
    ratio holds the occupancy alone. Returns the factor and its variance.
    """
    combinatoric, combinatoric_variance = np.asarray(combinatoric, dtype=float), np.asarray(combinatoric_variance, dtype=float)
    if combinatoric.ndim != 4 or combinatoric_variance.shape != combinatoric.shape:
        raise ValueError("the combinatoric Delta-phi arrays must be (n_pt, n_xj, n_dphi, 2)")
    window, sideband = bin_mask(edges, *WINDOW), bin_mask(edges, *SIDEBAND)
    inside, inside_variance = (a[:, :, window, :].sum(axis=(1, 2, 3)) for a in (combinatoric, combinatoric_variance))
    outside, outside_variance = (a[:, :, sideband, :].sum(axis=(1, 2, 3)) for a in (combinatoric, combinatoric_variance))
    if not np.all((inside > 0) & (outside > 0)):
        raise ValueError("the window occupancy needs simulated combinatoric jets in the window and the sideband of every photon-pT bin")
    widths = (WINDOW[1] - WINDOW[0]) / (SIDEBAND[1] - SIDEBAND[0])
    factor = inside / outside / widths
    return {"factor": factor, "variance": factor ** 2 * (inside_variance / inside ** 2 + outside_variance / outside ** 2)}


def simulation_closure(truth_jet: np.ndarray, combinatoric: np.ndarray, truth_jet_variance: np.ndarray,
                       combinatoric_variance: np.ndarray, *, shape: str = "flow", scope: str = "pt",
                       edges: np.ndarray = DPHI_EDGES, window_occupancy: np.ndarray | None = None) -> dict[str, Any]:
    """Closure of the estimate on simulation.

    Inputs are the jets around matched photons with and without a truth jet nearby (TreeToHists FillPhotonJetResponse
    dphi_truthjet, dphi_combinatoric), each (n_pt, n_xj, n_dphi, 2 sides) with its variance. The estimate is made from
    all of them, as on data (with ``window_occupancy`` when given); the truth is the combinatoric jets in the recoil
    window. Per photon-pT bin, and per xJ bin: estimate, truth, their variances and the pull. Also: the share of sideband
    jets with a truth jet nearby (jets of the hard scattering, which the method counts as combinatoric), and the flow fit
    of the far-side combinatoric jets alone.
    """
    arrays = [np.asarray(a, dtype=float) for a in (truth_jet, combinatoric, truth_jet_variance, combinatoric_variance)]
    truth_jet, combinatoric, truth_jet_variance, combinatoric_variance = arrays
    if truth_jet.ndim != 4 or any(a.shape != truth_jet.shape for a in arrays):
        raise ValueError("the simulation Delta-phi arrays must share one (n_pt, n_xj, n_dphi, 2) shape")
    result = estimate((truth_jet + combinatoric)[None], (truth_jet_variance + combinatoric_variance)[None], shape=shape,
                      scope=scope, edges=edges, window_occupancy=window_occupancy)
    window, sideband = bin_mask(edges, *WINDOW), bin_mask(edges, *SIDEBAND)
    counts, variance = result["counts"][0], result["variance"][0]
    truth = combinatoric[:, :, window, :].sum(axis=(2, 3))
    truth_variance = combinatoric_variance[:, :, window, :].sum(axis=(2, 3))

    def compare(e, ev, t, tv):
        spread = np.sqrt(ev + tv)
        return {"estimate": e, "estimate_variance": ev, "truth": t, "truth_variance": tv,
                "pull": np.divide(e - t, spread, out=np.zeros_like(spread), where=spread > 0)}

    sideband_all = (truth_jet + combinatoric)[:, :, sideband, :].sum(axis=(1, 2, 3))
    sideband_truth_jet = truth_jet[:, :, sideband, :].sum(axis=(1, 2, 3))
    combinatoric_flow = []
    for p in range(truth_jet.shape[0]):
        try:
            fit = fit_flow(combinatoric[p, :, :, 1].sum(axis=0), combinatoric_variance[p, :, :, 1].sum(axis=0), edges)
        except ValueError:
            combinatoric_flow.append(None)
        else:
            combinatoric_flow.append({"v2": fit["v2"], "v3": fit["v3"], "chi2": fit["chi2"], "ndf": fit["ndf"]})
    return {"shape": shape, "scope": scope, "fits": result["fits"],
            "per_pt_xj": compare(counts, variance, truth, truth_variance),
            "per_pt": compare(counts.sum(axis=1), variance.sum(axis=1), truth.sum(axis=1), truth_variance.sum(axis=1)),
            "sideband_share_with_truth_jet": np.divide(sideband_truth_jet, sideband_all, out=np.full_like(sideband_all, np.nan),
                                                       where=sideband_all > 0),
            "combinatoric_far_side_flow": combinatoric_flow}


__all__ = ["DPHI_EDGES", "FIT_RANGE", "SCOPES", "SHAPES", "SIDEBAND", "WINDOW", "basis", "bin_mask", "estimate", "fit_flow",
           "integral", "simulation_closure", "window_fraction", "window_occupancy"]
