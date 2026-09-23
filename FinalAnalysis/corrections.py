"""Leakage-aware ABCD purity and background corrections.

Ported from the validated collaboration implementation
(``FinalAnalysis/photonjet/analysis/background.py`` at main a7f15eb). The
arithmetic is unchanged; the inputs are the arrays of a
``TreeToHists`` histogram package instead of a JSON payload.

Inputs are selected, weighted sufficient statistics, not trees. ABCD counts
are ordered A/B/C/D; recoil spectra are ordered A/C, then photon-pT bin,
then xJ bin. The measured photon-pT grid comes from the package, never from
a plot label.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class CorrectionInputs:
    """The six response quantities used by the correction and toy kernels.

    Leakage arrays have shape (4, n_pt) in A/B/C/D order. Combinatoric arrays
    have shape (n_pt, n_xj); normalisations have shape (n_pt,). The
    normalisation counts truth-conditioned region-A photons, not every
    reconstructed photon. Arrays and their sumw2 come from one selection.
    """

    leakage: np.ndarray
    leakage_sumw2: np.ndarray
    combinatoric_reco: np.ndarray
    combinatoric_reco_sumw2: np.ndarray
    combinatoric_normalization_reco: np.ndarray
    combinatoric_normalization_reco_sumw2: np.ndarray


@dataclass(frozen=True)
class ABCDSolution:
    signal_a: float
    purity: float
    background_a_over_c: float
    residual: float
    method: str


def solve_leakage_abcd(counts: Sequence[float], leakage: Sequence[float]) -> ABCDSolution:
    """Solve the factorised ABCD background with signal leakage relative to A."""

    n = np.asarray(counts, dtype=float)
    l = np.asarray(leakage, dtype=float)
    if n.shape != (4,) or l.shape != (4,) or np.any(~np.isfinite(n)) or np.any(n < 0) or np.any(~np.isfinite(l)) or np.any(l < 0):
        raise ValueError("ABCD counts/leakage must be four finite nonnegative values")
    if l[0] <= 0:
        raise ValueError("ABCD signal leakage requires positive A normalization")
    l = l / l[0]
    a = l[3] - l[1] * l[2]
    b = -n[0] * l[3] - n[3] + n[1] * l[2] + n[2] * l[1]
    c = n[0] * n[3] - n[1] * n[2]
    upper = n[0]
    for value, ratio in zip(n[1:], l[1:]):
        if ratio > 0:
            upper = min(upper, value / ratio)
    naive = float(np.clip(n[0] - n[1] * n[2] / n[3], 0.0, upper)) if n[3] > 0 else upper
    roots: list[float] = []
    if abs(a) < 1.0e-14:
        if abs(b) > 1.0e-14:
            roots = [-c / b]
    else:
        discriminant = b * b - 4.0 * a * c
        if discriminant >= 0:
            root = math.sqrt(discriminant)
            roots = [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]
    physical = [float(value) for value in roots if -1.0e-9 <= value <= upper + 1.0e-9]
    method = "analytic_physical_root"
    if physical:
        signal = min(physical, key=lambda value: abs(value - naive))
    else:
        # When the leakage-aware equation has no physical solution, keep the
        # clipped raw ABCD estimate; a residual-minimising boundary point is
        # not an ABCD solution.
        signal = naive
        method = "raw_abcd_fallback_no_physical_leakage_root"
    backgrounds = n - l * signal
    bkg_c = backgrounds[2]
    alpha = backgrounds[0] / bkg_c if bkg_c > 0 else 0.0
    residual = abs(backgrounds[0] * backgrounds[3] - backgrounds[1] * backgrounds[2])
    scale = max(1.0, abs(backgrounds[0] * backgrounds[3]), abs(backgrounds[1] * backgrounds[2]))
    return ABCDSolution(
        signal_a=max(0.0, signal),
        purity=max(0.0, min(1.0, signal / n[0])) if n[0] > 0 else 0.0,
        background_a_over_c=max(0.0, alpha),
        residual=float(residual / scale),
        method=method,
    )


def correct_measured(
    data: Mapping[str, np.ndarray],
    response: CorrectionInputs,
    subtract_unmatched_recoil: bool,
    *,
    pt_edges: np.ndarray,
    combinatoric_reco: np.ndarray | None = None,
    combinatoric_normalization_reco: np.ndarray | None = None,
    purity_strategy: str = "per_pt",
) -> dict[str, Any]:
    """Purity-correct the region-A recoil spectrum in every photon-pT bin.

    ``data`` carries ``abcd_counts`` (4, n_pt), ``abcd_counts_sumw2``,
    ``recoil_spectra`` (2, n_pt, n_xj) and ``recoil_spectra_sumw2``.
    """

    pt_edges = np.asarray(pt_edges, dtype=float)
    n_pt = len(pt_edges) - 1
    counts = np.asarray(data["abcd_counts"], dtype=float)
    spectra = np.asarray(data["recoil_spectra"], dtype=float)
    spectra_sumw2 = np.asarray(data["recoil_spectra_sumw2"], dtype=float)
    if counts.shape != (4, n_pt) or spectra.ndim != 3 or spectra.shape[:2] != (2, n_pt) or spectra_sumw2.shape != spectra.shape:
        raise ValueError("measured sufficient-statistic dimensions differ")
    if any(np.any(~np.isfinite(v)) or np.any(v < 0) for v in (counts, spectra, spectra_sumw2)):
        raise ValueError("measured statistics must be finite and nonnegative")
    corrected = np.zeros_like(spectra[0])
    variance = np.zeros_like(spectra_sumw2[0])
    photons = np.zeros(n_pt)
    solutions: list[dict[str, Any]] = []
    comb_subtracted = np.zeros_like(corrected)
    comb_variance = np.zeros_like(corrected)
    comb_template = response.combinatoric_reco if combinatoric_reco is None else np.asarray(combinatoric_reco, dtype=float)
    comb_normalization = (response.combinatoric_normalization_reco if combinatoric_normalization_reco is None
                          else np.asarray(combinatoric_normalization_reco, dtype=float))
    if comb_template.shape != corrected.shape:
        raise ValueError("combinatoric template shape differs from measured recoil shape")
    if comb_normalization.shape != photons.shape:
        raise ValueError("combinatoric photon-normalization shape differs")
    if subtract_unmatched_recoil and (np.any(~np.isfinite(comb_normalization)) or np.any(comb_normalization <= 0)):
        raise ValueError("combinatoric subtraction requires positive finite denominators")
    if purity_strategy not in {"per_pt", "integrated_15_35_transfer"}:
        raise ValueError(f"unknown purity strategy: {purity_strategy}")

    integrated_solution: ABCDSolution | None = None
    if purity_strategy == "integrated_15_35_transfer":
        if not np.isclose(pt_edges[0], 15) or not np.isclose(pt_edges[-1], 35):
            raise ValueError("integrated_15_35_transfer is defined only over 15-35 GeV")
        # One leakage-aware background transfer from the integrated A/B/C/D
        # population, applied inside each reconstructed photon-pT bin so the
        # response still sees the measured pT spectrum. Region A and C
        # spectra stay bin-local.
        integrated_solution = solve_leakage_abcd(counts.sum(axis=1), response.leakage.sum(axis=1))
        if integrated_solution.residual >= 0.02:
            raise ValueError(f"integrated ABCD transfer does not factorize: residual={integrated_solution.residual}")
    for pt in range(n_pt):
        local_solution = solve_leakage_abcd(counts[:, pt], response.leakage[:, pt])
        solution = integrated_solution or local_solution
        transfer = solution.background_a_over_c
        leakage_a = response.leakage[0, pt]
        leakage_c = response.leakage[2, pt]
        if leakage_a <= 0:
            raise ValueError(f"missing Region-A prompt leakage normalization in photon-pT bin {pt}")
        relative_c_leakage = leakage_c / leakage_a
        leakage_denominator = 1.0 - transfer * relative_c_leakage
        if not math.isfinite(leakage_denominator) or leakage_denominator <= 1.0e-6:
            raise ValueError(f"nonphysical Region-C leakage denominator in photon-pT bin {pt}: 1-k*fC={leakage_denominator}")
        if purity_strategy == "integrated_15_35_transfer":
            photons[pt] = (counts[0, pt] - transfer * counts[2, pt]) / leakage_denominator
            if not math.isfinite(photons[pt]) or photons[pt] <= 0:
                raise ValueError(f"integrated ABCD transfer gives a nonpositive signal in photon-pT bin {pt}: {photons[pt]}")
        else:
            photons[pt] = local_solution.signal_a
        numerator = spectra[0, pt] - transfer * spectra[1, pt]
        corrected[pt] = numerator / leakage_denominator
        variance[pt] = (spectra_sumw2[0, pt] + transfer ** 2 * spectra_sumw2[1, pt]) / (leakage_denominator ** 2)
        if subtract_unmatched_recoil and comb_normalization[pt] > 0:
            # ABCD already removed fake-photon recoil; scale only the matched
            # prompt photon with an unmatched reconstructed recoil.
            comb_subtracted[pt] = comb_template[pt] * photons[pt] / comb_normalization[pt]
            corrected[pt] -= comb_subtracted[pt]
            template = comb_template[pt]
            template_var = response.combinatoric_reco_sumw2[pt]
            denominator = comb_normalization[pt]
            denominator_var = response.combinatoric_normalization_reco_sumw2[pt]
            scale = photons[pt] / denominator
            derivative_denominator = -photons[pt] * template / (denominator * denominator)
            comb_variance[pt] = np.clip(scale * scale * template_var + derivative_denominator ** 2 * denominator_var, 0.0, None)
            variance[pt] += comb_variance[pt]
        solutions.append({
            "pt_low": float(pt_edges[pt]), "pt_high": float(pt_edges[pt + 1]),
            "signal_a": float(photons[pt]),
            "purity": float(photons[pt] / counts[0, pt]) if counts[0, pt] > 0 else 0.0,
            "background_a_over_c": transfer,
            "relative_region_c_signal_leakage": float(relative_c_leakage),
            "region_c_signal_leakage_restoration": float(1.0 / leakage_denominator),
            "factorization_residual": solution.residual,
            "local_factorization_residual": local_solution.residual,
            "method": "integrated_15_35_leakage_transfer" if integrated_solution is not None else local_solution.method,
            "combinatoric_sim_yield_per_photon": (float(comb_template[pt].sum() / comb_normalization[pt])
                                                 if subtract_unmatched_recoil and comb_normalization[pt] > 0 else 0.0),
            "combinatoric_scaled_recoil": float(comb_subtracted[pt].sum()),
            "purity_corrected_recoil_before_combinatoric": float((corrected[pt] + comb_subtracted[pt]).sum()),
        })
    return {
        "corrected": corrected, "variance": variance, "photons": photons, "solutions": solutions,
        "combinatoric_subtracted": comb_subtracted, "combinatoric_variance": comb_variance,
        "negative_input_sum": float(np.abs(corrected[corrected < 0]).sum()), "purity_strategy": purity_strategy,
    }


def draw_joint_data_toy(data: Mapping[str, np.ndarray], rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Fluctuate the count and pair sufficient statistics independently.

    Several recoil jets may accompany one photon, so the xJ bins do not
    partition the event-leading counts; the two are fluctuated separately.
    """

    counts = np.asarray(data["abcd_counts"], dtype=float)
    counts_var = np.asarray(data["abcd_counts_sumw2"], dtype=float)
    spectra = np.asarray(data["recoil_spectra"], dtype=float)
    spectra_var = np.asarray(data["recoil_spectra_sumw2"], dtype=float)
    toy_counts = np.clip(rng.normal(counts, np.sqrt(np.clip(counts_var, 0.0, None))), 0.0, None)
    toy_spectra = np.clip(rng.normal(spectra, np.sqrt(np.clip(spectra_var, 0.0, None))), 0.0, None)
    return {**data, "abcd_counts": toy_counts, "recoil_spectra": toy_spectra}


def bootstrap_data_toy(data: Mapping[str, np.ndarray], replica: int) -> dict[str, np.ndarray]:
    """The event-bootstrap replica stored by TreeToHists, as a toy data set."""

    return {**data,
            "abcd_counts": np.asarray(data["bootstrap_abcd_counts"][replica], dtype=float),
            "recoil_spectra": np.asarray(data["bootstrap_recoil_spectra"][replica], dtype=float)}


def draw_combinatoric_template_toy(response: CorrectionInputs, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    template = np.asarray(response.combinatoric_reco, dtype=float)
    template_var = np.asarray(response.combinatoric_reco_sumw2, dtype=float)
    denominator = np.asarray(response.combinatoric_normalization_reco, dtype=float)
    denominator_var = np.asarray(response.combinatoric_normalization_reco_sumw2, dtype=float)
    toy_template = np.clip(rng.normal(template, np.sqrt(np.clip(template_var, 0.0, None))), 0.0, None)
    toy_denominator = np.clip(rng.normal(denominator, np.sqrt(np.clip(denominator_var, 0.0, None))), 0.0, None)
    return toy_template, toy_denominator


__all__ = ["ABCDSolution", "CorrectionInputs", "bootstrap_data_toy", "correct_measured",
           "draw_combinatoric_template_toy", "draw_joint_data_toy", "solve_leakage_abcd"]
