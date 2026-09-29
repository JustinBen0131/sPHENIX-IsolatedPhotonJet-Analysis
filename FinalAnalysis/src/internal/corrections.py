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
class PurityCorrections:
    """Optional corrections of the measured photons and pairs, one value per photon-pT bin.

    ``background_correlation`` is R_bkg = N_A^bkg N_D^bkg / (N_B^bkg N_C^bkg) of the non-signal candidates. One means the
    background factorises between the identification and isolation axes, the nominal assumption of ATLAS and PPG12; a
    measured value corrects for their correlation (ATLAS, Phys. Rev. C 93 (2016) 034914, Eq. 4). ``scale`` multiplies
    the leakage-corrected signal (a purity systematic such as the simulation nonclosure P_truth/P_ABCD); the background
    transfer is then re-derived from the scaled signal. Both need the per_pt purity strategy.
    ``efficiency_scale`` is the data/simulation photon identification and isolation efficiency per reconstructed
    photon-pT bin: after every background correction the bin's photons and pairs are divided by it together, so that the
    simulated response, which carries the simulated efficiency, unfolds the data's own. Within a bin it cancels in the
    per-photon ratio; it acts through the migrations and the weight of each bin in the pT-integrated distribution.
    None leaves the nominal correction unchanged, bit for bit.
    """

    background_correlation: tuple[float, ...] | None = None
    scale: tuple[float, ...] | None = None
    efficiency_scale: tuple[float, ...] | None = None


@dataclass(frozen=True)
class ABCDSolution:
    signal_a: float
    purity: float
    background_a_over_c: float
    residual: float
    method: str


def solve_leakage_abcd(counts: Sequence[float], leakage: Sequence[float], background_correlation: float = 1.0) -> ABCDSolution:
    """Solve the ABCD background with signal leakage relative to A: (A - s) = R (B - l_B s)(C - l_C s)/(D - l_D s), R the
    background correlation (1: factorised)."""

    n = np.asarray(counts, dtype=float)
    l = np.asarray(leakage, dtype=float)
    if n.shape != (4,) or l.shape != (4,) or np.any(~np.isfinite(n)) or np.any(n < 0) or np.any(~np.isfinite(l)) or np.any(l < 0):
        raise ValueError("ABCD counts/leakage must be four finite nonnegative values")
    if l[0] <= 0:
        raise ValueError("ABCD signal leakage requires positive A normalization")
    r = float(background_correlation)
    if not math.isfinite(r) or r <= 0:
        raise ValueError("the background correlation must be finite and positive")
    l = l / l[0]
    a = l[3] - r * l[1] * l[2]
    b = -n[0] * l[3] - n[3] + r * n[1] * l[2] + r * n[2] * l[1]
    c = n[0] * n[3] - r * n[1] * n[2]
    upper = n[0]
    for value, ratio in zip(n[1:], l[1:]):
        if ratio > 0:
            upper = min(upper, value / ratio)
    naive = float(np.clip(n[0] - r * n[1] * n[2] / n[3], 0.0, upper)) if n[3] > 0 else upper
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
    residual = abs(backgrounds[0] * backgrounds[3] - r * backgrounds[1] * backgrounds[2])
    scale = max(1.0, abs(backgrounds[0] * backgrounds[3]), abs(r * backgrounds[1] * backgrounds[2]))
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
    purity_corrections: PurityCorrections | None = None,
) -> dict[str, Any]:
    """Purity-correct the region-A recoil spectrum in every photon-pT bin.

    ``purity_corrections`` (PurityCorrections): the background correlation R_bkg enters the ABCD solution; a scale
    multiplies the leakage-corrected signal and the transfer k = (A - s)/(C - l_C s) is re-derived from it, so the
    corrected pairs still integrate to the scaled photons (both per_pt only); an efficiency scale divides each bin's
    photons and pairs as the last step (any strategy). None leaves the correction unchanged.

    ``data`` carries ``abcd_counts`` (4, n_pt), ``abcd_counts_sumw2``,
    ``recoil_spectra`` (2, n_pt, n_xj) and ``recoil_spectra_sumw2``, and
    optionally ``recoil_background`` with ``recoil_background_sumw2`` (the
    same shape): a data-driven combinatoric estimate (combinatoric.py),
    removed from the region A and C spectra before the purity correction.

    Besides the corrected spectrum it returns its intermediate steps, for display only and before an efficiency scale:
    the region A and C spectra it corrected (after a data-driven background), the transferred background
    ``scaled_region_c`` (k C), and the leakage-restored ``purity_corrected`` spectrum before a simulated-template
    subtraction, with its variance.
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
    if "recoil_background" in data:
        background = np.asarray(data["recoil_background"], dtype=float)
        background_sumw2 = np.asarray(data["recoil_background_sumw2"], dtype=float)
        if background.shape != spectra.shape or background_sumw2.shape != spectra.shape:
            raise ValueError("recoil background dimensions differ from the recoil spectra")
        if np.any(~np.isfinite(background)) or np.any(~np.isfinite(background_sumw2)) or np.any(background_sumw2 < 0):
            raise ValueError("recoil background must be finite")
        spectra, spectra_sumw2 = spectra - background, spectra_sumw2 + background_sumw2
    corrected = np.zeros_like(spectra[0])
    variance = np.zeros_like(spectra_sumw2[0])
    photons = np.zeros(n_pt)
    solutions: list[dict[str, Any]] = []
    comb_subtracted = np.zeros_like(corrected)
    comb_variance = np.zeros_like(corrected)
    scaled_region_c = np.zeros_like(corrected)
    purity_corrected = np.zeros_like(corrected)
    purity_corrected_variance = np.zeros_like(corrected)
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
    purity_scale = correlation = efficiency = None
    if purity_corrections is not None:
        for name in ("scale", "background_correlation", "efficiency_scale"):
            values = getattr(purity_corrections, name)
            if values is None:
                continue
            if name != "efficiency_scale" and purity_strategy != "per_pt":
                raise ValueError("the purity scale and background correlation apply to the per_pt purity strategy only")
            values = np.asarray(values, dtype=float)
            if values.shape != (n_pt,) or np.any(~np.isfinite(values)) or np.any(values <= 0):
                raise ValueError(f"the purity {name} needs one finite positive value per photon-pT bin")
            if name == "scale":
                purity_scale = values
            elif name == "background_correlation":
                correlation = values
            else:
                efficiency = values

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
        local_solution = (solve_leakage_abcd(counts[:, pt], response.leakage[:, pt]) if correlation is None
                          else solve_leakage_abcd(counts[:, pt], response.leakage[:, pt], correlation[pt]))
        solution = integrated_solution or local_solution
        transfer = solution.background_a_over_c
        leakage_a = response.leakage[0, pt]
        leakage_c = response.leakage[2, pt]
        if leakage_a <= 0:
            raise ValueError(f"missing Region-A prompt leakage normalization in photon-pT bin {pt}")
        relative_c_leakage = leakage_c / leakage_a
        scaled_signal = None
        if purity_scale is not None:
            scaled_signal = purity_scale[pt] * local_solution.signal_a
            background_c = counts[2, pt] - relative_c_leakage * scaled_signal
            if not (0.0 < scaled_signal < counts[0, pt]) or background_c <= 0:
                raise ValueError(f"the purity scale {purity_scale[pt]} leaves no physical ABCD solution in photon-pT bin {pt}")
            transfer = (counts[0, pt] - scaled_signal) / background_c
        leakage_denominator = 1.0 - transfer * relative_c_leakage
        if not math.isfinite(leakage_denominator) or leakage_denominator <= 1.0e-6:
            raise ValueError(f"nonphysical Region-C leakage denominator in photon-pT bin {pt}: 1-k*fC={leakage_denominator}")
        if purity_strategy == "integrated_15_35_transfer":
            photons[pt] = (counts[0, pt] - transfer * counts[2, pt]) / leakage_denominator
            if not math.isfinite(photons[pt]) or photons[pt] <= 0:
                raise ValueError(f"integrated ABCD transfer gives a nonpositive signal in photon-pT bin {pt}: {photons[pt]}")
        else:
            photons[pt] = local_solution.signal_a if scaled_signal is None else scaled_signal
        numerator = spectra[0, pt] - transfer * spectra[1, pt]
        corrected[pt] = numerator / leakage_denominator
        variance[pt] = (spectra_sumw2[0, pt] + transfer ** 2 * spectra_sumw2[1, pt]) / (leakage_denominator ** 2)
        scaled_region_c[pt] = transfer * spectra[1, pt]
        purity_corrected[pt] = corrected[pt]
        purity_corrected_variance[pt] = variance[pt]
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
        if purity_scale is not None:
            solutions[-1]["purity_scale"] = float(purity_scale[pt])
        if correlation is not None:
            solutions[-1]["background_correlation"] = float(correlation[pt])
    if efficiency is not None:
        # the last step: the bin's photons and pairs together (the ABCD solutions above keep their own purity)
        corrected, variance, photons = corrected / efficiency[:, None], variance / efficiency[:, None] ** 2, photons / efficiency
        for pt, row in enumerate(solutions):
            row["efficiency_scale"] = float(efficiency[pt])
    return {
        "corrected": corrected, "variance": variance, "photons": photons, "solutions": solutions,
        "combinatoric_subtracted": comb_subtracted, "combinatoric_variance": comb_variance,
        "region_spectra": spectra, "region_spectra_sumw2": spectra_sumw2, "scaled_region_c": scaled_region_c,
        "purity_corrected": purity_corrected, "purity_corrected_variance": purity_corrected_variance,
        "negative_input_sum": float(np.abs(corrected[corrected < 0]).sum()), "purity_strategy": purity_strategy,
    }


__all__ = ["ABCDSolution", "CorrectionInputs", "PurityCorrections", "correct_measured", "solve_leakage_abcd"]
