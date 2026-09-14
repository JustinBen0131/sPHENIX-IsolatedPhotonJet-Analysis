"""Purity correction, unfolding, uncertainties and iteration selection.

One function, :func:`run_chain`, takes the four inputs of the final analysis

  data      TreeToHists payload of the data sample
  leakage   TreeToHists payload of the photon+jet simulation restricted to
            nominal truth-signal photons (prompt-photon leakage into B, C, D)
  window    the pair and photon responses restricted to the analysis window
            (:func:`restrict`)

and produces the unfolded (1/N_gamma) dN/dxJgamma with statistical
uncertainties and closure diagnostics, in this order:

  1. leakage-aware ABCD purity per photon-pT bin and the purity-corrected
     recoil spectrum (region A minus the background transfer times region C,
     with the region-C prompt leakage restored);
  2. optional combinatoric subtraction: the simulated recoil yield of matched
     prompt photons whose jet has no truth counterpart, scaled to the data
     photon count (Au+Au only in the maintained analysis);
  3. joint iterative-Bayes unfolding of the (photon pT, xJ) spectrum with the
     truth spectrum as prior and only detector and boundary fakes as the fake
     cause (fake photons are removed by 1, combinatoric recoil by 2);
  4. one-dimensional unfolding of the purity-corrected photon count with the
     photon response, prior the truth photon spectrum, fake cause the photons
     whose truth pT lies outside 15-35 GeV;
  5. density = unfolded pairs summed over pT / (unfolded photons x bin width);
  6. toys: Gaussian fluctuations of the data sufficient statistics and, when
     2 is on, of the combinatoric template; each toy reruns 1-5; the error is
     the central 68 percent half-width and the covariance the winsorized toy
     covariance rescaled to those errors;
  7. diagnostics: refolding chi2/ndf, photon refolding chi2/ndf, simulation
     closure (matched reco plus residual fakes unfolded against truth), zero
     efficiency bins, negative-input fraction, toy success fraction.

:func:`select_iterations` scans the iteration count 2..12, applies the
maintained gate and score to each candidate and picks the minimum score.
Response statistics are not propagated to the uncertainty here; the
maintained analysis does not either (leakage and matrix are held fixed).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Sequence

import numpy as np

from photonjet.analysis.background import (
    PT_EDGES,
    CorrectionInputs,
    correct_measured,
    draw_combinatoric_template_toy,
    draw_joint_data_toy,
)
from photonjet.analysis.photon_response import PhotonResponse
from photonjet.analysis.response import Category, ResponseBundle
from photonjet.analysis.unfolding import chi2_ndf, iterative_bayes

REPORTED_XJ_LOW = 0.30
REPORTED_XJ_HIGH = 1.80
ITERATION_CANDIDATES = tuple(range(2, 13))
FEED_IN_PREFIXES = ("PHOTON_LOW_FEED_IN", "PHOTON_HIGH_FEED_IN", "XJ_LOW_FEED_IN", "XJ_HIGH_FEED_IN")


@dataclass
class AnalysisResponse:
    """Pair and photon responses restricted to the 15/20/25/35 GeV photon-pT window."""

    xj_edges: np.ndarray
    xj_response: np.ndarray            # (3 nx, 3 nx) truth x reco
    xj_misses: np.ndarray              # (3 nx)
    xj_truth: np.ndarray               # (3, nx)
    xj_reco: np.ndarray                # (3, nx)
    xj_boundary_fakes: np.ndarray      # (3, nx) matched pairs with truth outside the window
    xj_detector_fakes: np.ndarray      # (3, nx) linked photon and jet, not the selected truth pair
    xj_fake_photon: np.ndarray         # (3, nx) reco photon without a truth link (removed by ABCD)
    xj_combinatoric: np.ndarray        # (3, nx) linked photon, jet without a truth jet (removed by K)
    xj_combinatoric_sumw2: np.ndarray
    photon_response: np.ndarray        # (3, 3) truth x reco
    photon_misses: np.ndarray          # (3)
    photon_truth: np.ndarray           # (3)
    photon_reco: np.ndarray            # (3)
    photon_boundary_fakes: np.ndarray  # (3)
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def n_xj(self) -> int:
        return len(self.xj_edges) - 1

    @property
    def unfolding_fakes(self) -> np.ndarray:
        return (self.xj_detector_fakes + self.xj_boundary_fakes).reshape(-1)


def _window_bins(edges: np.ndarray) -> list[int]:
    """Indices of the source bins equal to the three analysis pT bins."""

    edges = np.asarray(edges, dtype=float)
    indices = []
    for low, high in zip(PT_EDGES[:-1], PT_EDGES[1:]):
        match = [i for i in range(len(edges) - 1) if math.isclose(edges[i], low, abs_tol=1e-9) and math.isclose(edges[i + 1], high, abs_tol=1e-9)]
        if len(match) != 1:
            raise ValueError(f"analysis bin [{low}, {high}) is not a single bin of the response grid {edges.tolist()}")
        indices.append(match[0])
    return indices


def restrict(bundle: ResponseBundle, photon: PhotonResponse) -> AnalysisResponse:
    """Cut the analysis window out of the classification-grid responses.

    Truth rows inside the window that reconstruct outside it become misses;
    truth rows outside the window that reconstruct inside it become boundary
    fakes.  Feed-in categories recorded by the pair builder (truth off the
    classification support) join the boundary fakes.
    """

    if bundle.dimension != "2D" or bundle.xj_edges is None:
        raise ValueError("the pair response must be two-dimensional")
    nt, nr, nx = len(bundle.truth_ptgamma_edges) - 1, len(bundle.reco_ptgamma_edges) - 1, len(bundle.xj_edges) - 1
    t_in = _window_bins(bundle.truth_ptgamma_edges)
    r_in = _window_bins(bundle.reco_ptgamma_edges)
    t_out = [i for i in range(nt) if i not in t_in]
    r_out = [i for i in range(nr) if i not in r_in]
    matrix = bundle.matrix.reshape(nt, nx, nr, nx)

    xj_response = matrix[t_in][:, :, r_in, :].reshape(3 * nx, 3 * nx)
    feed_out = matrix[t_in][:, :, r_out, :].sum(axis=(2, 3)) if r_out else np.zeros((3, nx))
    xj_misses = (bundle.all_misses.reshape(nt, nx)[t_in] + feed_out).reshape(-1)
    feed_in = matrix[t_out][:, :, r_in, :].sum(axis=(0, 1)) if t_out else np.zeros((3, nx))
    boundary = feed_in.copy()
    detector = np.zeros((3, nx))
    for key, values in bundle.boundary_reco.items():
        block = values.reshape(nr, nx)[r_in]
        if key == Category.ASSIGNED_HARD_NONFIDUCIAL.value:
            detector += block
        else:
            boundary += block
    causes = {name: values.reshape(nr, nx)[r_in] for name, values in bundle.fake_causes.items()}
    causes_w2 = {name: values.reshape(nr, nx)[r_in] for name, values in bundle.fake_causes_sumw2.items()}
    zeros = np.zeros((3, nx))

    pt_in = _window_bins(photon.truth_edges)
    pr_in = _window_bins(photon.reco_edges)
    pt_out = [i for i in range(len(photon.truth_edges) - 1) if i not in pt_in]
    pr_out = [i for i in range(len(photon.reco_edges) - 1) if i not in pr_in]
    photon_response = photon.matrix[pt_in][:, pr_in]
    photon_misses = photon.misses[pt_in] + (photon.matrix[pt_in][:, pr_out].sum(axis=1) if pr_out else 0.0)
    photon_boundary = (photon.matrix[pt_out][:, pr_in].sum(axis=0) if pt_out else 0.0) + photon.boundary_fakes[pr_in]

    return AnalysisResponse(
        xj_edges=np.asarray(bundle.xj_edges, dtype=float),
        xj_response=xj_response, xj_misses=xj_misses,
        xj_truth=bundle.truth.reshape(nt, nx)[t_in], xj_reco=bundle.reco.reshape(nr, nx)[r_in],
        xj_boundary_fakes=boundary, xj_detector_fakes=detector,
        xj_fake_photon=causes.get(Category.UNMATCHED_RECO.value, zeros),
        xj_combinatoric=causes.get(Category.COMBINATORIC.value, zeros),
        xj_combinatoric_sumw2=causes_w2.get(Category.COMBINATORIC.value, zeros),
        photon_response=photon_response, photon_misses=photon_misses,
        photon_truth=photon.truth[pt_in], photon_reco=photon.reco[pr_in], photon_boundary_fakes=photon_boundary,
        provenance={"pair_response": bundle.provenance, "photon_response": photon.provenance},
    )


def correction_inputs(leakage: Mapping[str, Any], window: AnalysisResponse) -> CorrectionInputs:
    """Leakage from the truth-signal-only payload; combinatoric template from the window.

    The template normalisation is the Region-A prompt-photon count per reco
    pT bin (the same simulation photons that produced the template).
    """

    counts = np.asarray(leakage["abcd_event_leading_counts_ptgamma"], dtype=float)
    counts_w2 = np.asarray(leakage["abcd_event_leading_sumw2_ptgamma"], dtype=float)
    if counts.shape != (4, len(PT_EDGES) - 1):
        raise ValueError("leakage payload must carry 4 x 3 event-leading counts")
    return CorrectionInputs(
        leakage=counts, leakage_sumw2=counts_w2,
        combinatoric_reco=window.xj_combinatoric, combinatoric_reco_sumw2=window.xj_combinatoric_sumw2,
        combinatoric_normalization_reco=counts[0], combinatoric_normalization_reco_sumw2=counts_w2[0],
    )


def _measured_arrays(data: Mapping[str, Any]) -> dict[str, np.ndarray]:
    return {
        "abcd_event_leading_counts_ptgamma": np.asarray(data["abcd_event_leading_counts_ptgamma"], dtype=float),
        "abcd_event_leading_sumw2_ptgamma": np.asarray(data["abcd_event_leading_sumw2_ptgamma"], dtype=float),
        "inclusive_recoil_spectra_ptgamma_xj": np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float),
        "inclusive_recoil_sumw2_ptgamma_xj": np.asarray(data["inclusive_recoil_sumw2_ptgamma_xj"], dtype=float),
    }


def unfold_pairs(measured: np.ndarray, window: AnalysisResponse, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    return iterative_bayes(np.asarray(measured).reshape(-1), window.xj_response, window.xj_misses, iterations,
                           window.xj_truth.reshape(-1), fakes=window.unfolding_fakes)


def unfold_photons(photons: np.ndarray, window: AnalysisResponse, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    return iterative_bayes(np.asarray(photons), window.photon_response, window.photon_misses, iterations,
                           window.photon_truth, fakes=window.photon_boundary_fakes)


def densities(unfolded: np.ndarray, photons: np.ndarray, xj_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(1/N_gamma) dN/dxJ summed over pT and per pT bin (rows), bin-width divided."""

    n_pt, nx = len(PT_EDGES) - 1, len(xj_edges) - 1
    widths = np.diff(xj_edges)
    blocks = np.asarray(unfolded).reshape(n_pt, nx)
    total = float(np.sum(photons))
    summed = np.divide(blocks.sum(axis=0), total * widths, out=np.zeros(nx), where=(total * widths) > 0)
    per_pt = np.divide(blocks, photons[:, None] * widths[None, :], out=np.zeros((n_pt, nx)),
                       where=(photons[:, None] * widths[None, :]) > 0)
    return summed, per_pt


def reported_mask(xj_edges: np.ndarray) -> np.ndarray:
    edges = np.asarray(xj_edges, dtype=float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return (edges[:-1] >= REPORTED_XJ_LOW) & (centres <= REPORTED_XJ_HIGH)


def run_chain(
    data: Mapping[str, Any],
    leakage: Mapping[str, Any],
    window: AnalysisResponse,
    *,
    iterations: int,
    subtract_combinatoric: bool,
    purity_strategy: str = "per_pt",
    toys: int = 300,
    seed: int = 7,
) -> dict[str, Any]:
    xj_edges = window.xj_edges
    if not np.allclose(data["xj_edges"], xj_edges):
        raise ValueError("data xJ edges differ from the response xJ edges")
    nx = window.n_xj
    measured = _measured_arrays(data)
    inputs = correction_inputs(leakage, window)

    correction = correct_measured(measured, inputs, subtract_unmatched_recoil=subtract_combinatoric, purity_strategy=purity_strategy)
    corrected = np.asarray(correction["corrected"])
    variance = np.asarray(correction["variance"])
    photons = np.asarray(correction["photons"])
    unfolded, refolded = unfold_pairs(corrected, window, iterations)
    unfolded_photons, refolded_photons = unfold_photons(photons, window, iterations)
    density, density_per_pt = densities(unfolded, unfolded_photons, xj_edges)

    # --- toys ---------------------------------------------------------------
    rng = np.random.default_rng(seed)
    central_methods = tuple(row["method"] for row in correction["solutions"])
    toy_density: list[np.ndarray] = []
    toy_per_pt: list[np.ndarray] = []
    toy_measured: list[np.ndarray] = []
    for _ in range(toys):
        toy_data = draw_joint_data_toy(measured, rng)
        if subtract_combinatoric:
            template, normalisation = draw_combinatoric_template_toy(inputs, rng)
        else:
            template, normalisation = inputs.combinatoric_reco, inputs.combinatoric_normalization_reco
        try:
            toy = correct_measured(toy_data, inputs, subtract_unmatched_recoil=subtract_combinatoric,
                                   combinatoric_reco=template, combinatoric_normalization_reco=normalisation,
                                   purity_strategy=purity_strategy)
            if tuple(row["method"] for row in toy["solutions"]) != central_methods:
                continue                      # a different ABCD branch is another estimator, not a fluctuation
            toy_unfolded, _ = unfold_pairs(toy["corrected"], window, iterations)
            toy_photons, _ = unfold_photons(toy["photons"], window, iterations)
        except ValueError:
            continue
        if toy_photons.sum() <= 0:
            continue
        summed, per_pt = densities(toy_unfolded, toy_photons, xj_edges)
        toy_density.append(summed); toy_per_pt.append(per_pt)
        toy_measured.append(np.asarray(toy["corrected"]).reshape(-1))
    successes = len(toy_density)

    def central_68(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        lower = np.quantile(samples, 0.16, axis=0)
        upper = np.quantile(samples, 0.84, axis=0)
        return lower, upper, 0.5 * (upper - lower)

    if successes > 1:
        toy_array = np.asarray(toy_density)
        lower, upper, density_error = central_68(toy_array)
        winsorized = np.clip(toy_array, lower, upper)
        covariance = np.atleast_2d(np.cov(winsorized, rowvar=False, ddof=1))
        scale = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
        correlation = np.divide(covariance, scale[:, None] * scale[None, :], out=np.eye(nx), where=(scale[:, None] * scale[None, :]) > 0)
        density_covariance = correlation * density_error[:, None] * density_error[None, :]
        _, _, per_pt_error = central_68(np.asarray(toy_per_pt))
        measured_variance = np.asarray(toy_measured).var(axis=0, ddof=1)
    else:
        density_error = np.full(nx, np.nan); density_covariance = np.full((nx, nx), np.nan)
        per_pt_error = np.full_like(density_per_pt, np.nan)
        measured_variance = np.clip(variance.reshape(-1), 1e-12, None)

    # --- diagnostics --------------------------------------------------------
    measured_vector = corrected.reshape(-1)
    live = measured_variance > 0
    refold = chi2_ndf(refolded[live], measured_vector[live], measured_variance[live]) if live.any() else math.nan
    photon_refold = chi2_ndf(refolded_photons, photons, np.maximum(photons, 1.0))
    closure_reco = window.xj_response.sum(axis=0) + window.unfolding_fakes
    mc_unfolded, mc_refolded = iterative_bayes(closure_reco, window.xj_response, window.xj_misses, iterations,
                                               window.xj_truth.reshape(-1), fakes=window.unfolding_fakes)
    closure_truth = window.xj_truth.reshape(-1)
    supported = closure_truth > 0
    mc_closure = chi2_ndf(mc_unfolded[supported], closure_truth[supported], np.clip(closure_truth + mc_unfolded, 1e-12, None)[supported]) if supported.any() else math.nan
    reco_live = closure_reco > 0
    mc_refold = chi2_ndf(mc_refolded[reco_live], closure_reco[reco_live], np.maximum(closure_reco, 1.0)[reco_live]) if reco_live.any() else math.nan
    xj_totals = window.xj_response.sum(axis=1) + window.xj_misses
    xj_efficiency = np.divide(window.xj_response.sum(axis=1), xj_totals, out=np.zeros_like(xj_totals), where=xj_totals > 0)
    reported = np.tile(reported_mask(xj_edges), len(PT_EDGES) - 1)
    photon_totals = window.photon_response.sum(axis=1) + window.photon_misses
    photon_efficiency = np.divide(window.photon_response.sum(axis=1), photon_totals, out=np.zeros_like(photon_totals), where=photon_totals > 0)
    negative_fraction = float(np.mean(measured_vector[reported] < 0)) if reported.any() else 0.0

    return {
        "iterations": iterations,
        "subtract_combinatoric": subtract_combinatoric,
        "purity_strategy": purity_strategy,
        "xj_edges": xj_edges.tolist(),
        "pt_edges": PT_EDGES.tolist(),
        "purity": correction["solutions"],
        "corrected_ptgamma_xj": corrected.tolist(),
        "corrected_variance_ptgamma_xj": variance.tolist(),
        "purity_corrected_photons_ptgamma": photons.tolist(),
        "combinatoric_subtracted_ptgamma_xj": np.asarray(correction["combinatoric_subtracted"]).tolist(),
        "unfolded_ptgamma_xj": unfolded.reshape(len(PT_EDGES) - 1, nx).tolist(),
        "refolded_ptgamma_xj": refolded.reshape(len(PT_EDGES) - 1, nx).tolist(),
        "unfolded_photons_ptgamma": unfolded_photons.tolist(),
        "density": density.tolist(),
        "density_error": density_error.tolist(),
        "density_stat_covariance": density_covariance.tolist(),
        "density_per_ptgamma": density_per_pt.tolist(),
        "density_per_ptgamma_error": per_pt_error.tolist(),
        "reported_xj_mask": reported_mask(xj_edges).tolist(),
        "toy_requested": toys,
        "toy_successes": successes,
        "refold_chi2_ndf": float(refold),
        "photon_refold_chi2_ndf": float(photon_refold),
        "mc_closure_chi2_ndf": float(mc_closure),
        "mc_refold_chi2_ndf": float(mc_refold),
        "negative_input_fraction": negative_fraction,
        "response_observability": {
            "xj_zero_efficiency_supported_bins": int(np.sum(supported & (xj_efficiency == 0))),
            "xj_zero_efficiency_reported_bins": int(np.sum(supported & reported & (xj_efficiency == 0))),
            "photon_zero_efficiency_supported_bins": int(np.sum((window.photon_truth > 0) & (photon_efficiency == 0))),
        },
        "uncertainty_note": "data statistics through toys; response, leakage and template statistics held fixed",
    }


def chain_metrics(result: Mapping[str, Any]) -> dict[str, float | int | bool]:
    """Quality metrics of one chain result (maintained definitions)."""

    edges = np.asarray(result["xj_edges"], dtype=float)
    density = np.asarray(result["density"], dtype=float)
    error = np.asarray(result["density_error"], dtype=float)
    shown = reported_mask(edges)
    positive = shown & np.isfinite(density) & np.isfinite(error) & (density > 1.0e-10)
    relative = np.divide(error, density, out=np.full_like(error, np.inf), where=density > 0)
    selected = relative[positive]
    median_relative = float(np.median(selected)) if selected.size else math.inf
    max_relative = float(np.max(selected)) if selected.size else math.inf
    log_density = np.log(np.clip(density[positive], 1.0e-12, None))
    smoothness = float(np.mean(np.abs(np.diff(log_density, n=2)))) if log_density.size >= 3 else math.inf
    finite = bool(np.all(np.isfinite(density)) and np.all(density >= 0) and np.all(np.isfinite(error[positive])))
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
    return {
        "finite_nonnegative": finite,
        "stat_covariance_symmetric": symmetric,
        "stat_covariance_psd": psd,
        "stat_covariance_min_eigenvalue": minimum_eigenvalue,
        "purity_bin_count": int(len(purity_rows)),
        "purity_finite_physical": bool(len(purity_rows) == 3 and np.all(np.isfinite(purity_values))
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


def gate_and_score(result: Mapping[str, Any]) -> tuple[bool, float, dict[str, Any]]:
    """Maintained acceptance gate and ranking score of one iteration count."""

    m = chain_metrics(result)
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


def select_iterations(
    data: Mapping[str, Any],
    leakage: Mapping[str, Any],
    window: AnalysisResponse,
    *,
    candidates: Sequence[int] = ITERATION_CANDIDATES,
    scan_toys: int = 300,
    final_toys: int = 3000,
    seed: int = 7,
    **chain_options: Any,
) -> dict[str, Any]:
    """Scan the iteration count, gate and score each candidate, rerun the best with final toys."""

    table = []
    for k in candidates:
        result = run_chain(data, leakage, window, iterations=k, toys=scan_toys, seed=seed, **chain_options)
        passed, score, metrics = gate_and_score(result)
        table.append({"iterations": k, "passed": passed, "score": score, "metrics": metrics})
    best = min(table, key=lambda row: (row["score"], row["iterations"]))
    final = run_chain(data, leakage, window, iterations=best["iterations"], toys=final_toys, seed=seed + 1, **chain_options)
    final["selection"] = {"candidates": list(candidates), "scan_toys": scan_toys, "final_toys": final_toys,
                          "table": table, "selected_iterations": best["iterations"], "selected_passed_gate": best["passed"]}
    return final


def load_bundle(stem) -> ResponseBundle:
    """Read the pair response written by ``ResponseBundle.save`` (<stem>.npz and <stem>.json)."""

    import json
    from pathlib import Path

    stem = Path(stem)
    meta = json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
    if meta.get("schema") != "PhotonJetResponseBundleV1":
        raise ValueError(f"{stem}: not a PhotonJetResponseBundleV1 file")
    with np.load(stem.with_suffix(".npz")) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    groups: dict[str, dict[str, np.ndarray]] = {"fake_causes": {}, "fake_causes_sumw2": {},
                                                "boundary_reco": {}, "boundary_reco_sumw2": {},
                                                "boundary_truth": {}, "boundary_truth_sumw2": {}}
    plain = {}
    for name, values in arrays.items():
        prefix, _, key = name.partition("__")
        if key and prefix in groups:
            groups[prefix][key] = values
        else:
            plain[name] = values
    xj_edges = plain["xj_edges"] if plain["xj_edges"].size else None
    return ResponseBundle(
        system=meta["system"], dimension=meta["dimension"],
        truth_ptgamma_edges=plain["truth_ptgamma_edges"], reco_ptgamma_edges=plain["reco_ptgamma_edges"], xj_edges=xj_edges,
        matrix=plain["matrix"], matrix_sumw2=plain["matrix_sumw2"], truth=plain["truth"], truth_sumw2=plain["truth_sumw2"],
        reco=plain["reco"], reco_sumw2=plain["reco_sumw2"], misses=plain["misses"], misses_sumw2=plain["misses_sumw2"],
        fakes=plain["fakes"], fakes_sumw2=plain["fakes_sumw2"], **groups,
        unresolved=meta.get("unresolved", {}), provenance=meta.get("provenance", {}), status=meta.get("status", ""),
    )


__all__ = [
    "AnalysisResponse", "ITERATION_CANDIDATES", "chain_metrics", "correction_inputs", "densities",
    "gate_and_score", "load_bundle", "reported_mask", "restrict", "run_chain", "select_iterations", "unfold_pairs", "unfold_photons",
]
