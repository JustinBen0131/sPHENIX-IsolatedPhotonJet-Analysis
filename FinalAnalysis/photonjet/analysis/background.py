"""Leakage-aware ABCD and inclusive-recoil background corrections.

The same numerical implementation serves p+p and Au+Au. Inputs are already
selected, weighted sufficient statistics, not raw trees. ABCD counts are ordered
A/B/C/D; recoil spectra are ordered A/C, then photon-ET bin, then xJgamma bin.
The maintained measured photon-ET grid is 15, 20, 25, 35 GeV. Neither this module
nor its caller should infer a different grid from plot labels.

These kernels preserve the established correction arithmetic. Their isolated
equivalence does not certify the source selection, physicality of a particular
purity solution, full covariance chain, or a nominal physics result.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

import numpy as np


PT_EDGES = np.asarray([15.0, 20.0, 25.0, 35.0])


@dataclass(frozen=True)
class CorrectionInputs:
    """Only the six response quantities used by the correction and toy kernels.

    Leakage arrays have shape (4, 3) in A/B/C/D order. Combinatoric arrays have
    shape (3, n_xj); normalizations have shape (3,). The denominator counts
    truth-conditioned Region-A photons, not all reconstructed photon tags.
    Arrays and their sumw2 must come from the same selection and source set.

    For response uncertainty, pass each correlated source-bootstrap replica's
    leakage, template and normalization together. The independent template
    fluctuation helper below is a separately named approximation, not a
    substitute for that correlated ensemble.
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
    """Solve factorized ABCD backgrounds with signal leakage relative to A."""
    n = np.asarray(counts, dtype=float)
    l = np.asarray(leakage, dtype=float)
    if n.shape != (4,) or l.shape != (4,) or np.any(~np.isfinite(n)) or np.any(n < 0):
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
        # Match the established ROOT correction contract exactly: when the
        # leakage-aware fixed-point equation has no physical solution, retain
        # the clipped raw ABCD estimate.  A residual-minimizing boundary point
        # is not an ABCD solution and previously injected a spurious zero
        # 25--35 GeV Au+Au photon bin into the joint unfolding.
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
    combinatoric_reco: np.ndarray | None = None,
    combinatoric_normalization_reco: np.ndarray | None = None,
    purity_strategy: str = "per_pt",
) -> dict[str, Any]:
    counts = data["abcd_event_leading_counts_ptgamma"]
    spectra = data["inclusive_recoil_spectra_ptgamma_xj"]
    spectra_sumw2 = data["inclusive_recoil_sumw2_ptgamma_xj"]
    counts_sumw2 = data["abcd_event_leading_sumw2_ptgamma"]
    corrected = np.zeros_like(spectra[0])
    variance = np.zeros_like(spectra_sumw2[0])
    photons = np.zeros(len(PT_EDGES) - 1)
    solutions: list[dict[str, Any]] = []
    comb_subtracted = np.zeros_like(corrected)
    comb_variance = np.zeros_like(corrected)
    comb_template = (
        response.combinatoric_reco
        if combinatoric_reco is None
        else np.asarray(combinatoric_reco, dtype=float)
    )
    comb_normalization = (
        response.combinatoric_normalization_reco
        if combinatoric_normalization_reco is None
        else np.asarray(combinatoric_normalization_reco, dtype=float)
    )
    if comb_template.shape != corrected.shape:
        raise ValueError("combinatoric template shape differs from measured recoil shape")
    if comb_normalization.shape != photons.shape:
        raise ValueError("combinatoric photon-normalization shape differs")
    if purity_strategy not in {"per_pt", "integrated_15_35_transfer"}:
        raise ValueError(f"unknown purity strategy: {purity_strategy}")

    integrated_solution: ABCDSolution | None = None
    if purity_strategy == "integrated_15_35_transfer":
        # The reported observable is integrated over 15--35 GeV.  When the
        # upper-pT sidebands are too sparse to support three independent ABCD
        # solves, determine one leakage-aware background transfer from the
        # exact integrated A/B/C/D population.  Apply that transfer inside
        # each reconstructed photon-pT bin so the response still sees the
        # measured pT spectrum and its boundary migrations.  This is not a
        # pooled recoil shape: Region A and C spectra remain bin-local.
        integrated_solution = solve_leakage_abcd(
            counts.sum(axis=1),
            response.leakage.sum(axis=1),
        )
        if integrated_solution.residual >= 0.02:
            raise ValueError(
                "integrated 15--35 GeV ABCD transfer does not factorize: "
                f"residual={integrated_solution.residual}"
            )
    for pt in range(len(PT_EDGES) - 1):
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
            raise ValueError(
                f"nonphysical Region-C leakage denominator in photon-pT bin {pt}: "
                f"1-k*fC={leakage_denominator}"
            )
        if purity_strategy == "integrated_15_35_transfer":
            photons[pt] = (
                counts[0, pt] - transfer * counts[2, pt]
            ) / leakage_denominator
            if not math.isfinite(photons[pt]) or photons[pt] <= 0:
                raise ValueError(
                    "integrated 15--35 GeV ABCD transfer gives a nonpositive "
                    f"signal in photon-pT bin {pt}: {photons[pt]}"
                )
        else:
            photons[pt] = local_solution.signal_a
        numerator = spectra[0, pt] - transfer * spectra[1, pt]
        corrected[pt] = numerator / leakage_denominator
        variance[pt] = (
            spectra_sumw2[0, pt]
            + transfer ** 2 * spectra_sumw2[1, pt]
        ) / (leakage_denominator ** 2)
        if subtract_unmatched_recoil and comb_normalization[pt] > 0:
            # ABCD has already removed fake-photon recoil.  Scale only the
            # response component with a matched prompt photon and an unmatched
            # reconstructed recoil, thereby avoiding double subtraction.
            comb_subtracted[pt] = (
                comb_template[pt]
                * photons[pt] / comb_normalization[pt]
            )
            corrected[pt] -= comb_subtracted[pt]
            # Propagate finite embedded-template and photon-normalization
            # statistics.  Inclusive recoil bins are jet-pair yields and are
            # not a partition of the photon denominator, so no Bernoulli
            # subset covariance is invented here.
            template = comb_template[pt]
            template_var = response.combinatoric_reco_sumw2[pt]
            denominator = comb_normalization[pt]
            denominator_var = response.combinatoric_normalization_reco_sumw2[pt]
            scale = photons[pt] / denominator
            derivative_denominator = -photons[pt] * template / (denominator * denominator)
            comb_variance[pt] = np.clip(
                scale * scale * template_var
                + derivative_denominator * derivative_denominator * denominator_var,
                0.0,
                None,
            )
            variance[pt] += comb_variance[pt]
        solutions.append({
            "pt_low": float(PT_EDGES[pt]), "pt_high": float(PT_EDGES[pt + 1]),
            "signal_a": float(photons[pt]),
            "purity": float(photons[pt] / counts[0, pt]) if counts[0, pt] > 0 else 0.0,
            "background_a_over_c": transfer,
            "relative_region_c_signal_leakage": float(relative_c_leakage),
            "region_c_signal_leakage_restoration": float(1.0 / leakage_denominator),
            "factorization_residual": solution.residual,
            "local_factorization_residual": local_solution.residual,
            "method": (
                "integrated_15_35_leakage_transfer"
                if integrated_solution is not None else local_solution.method
            ),
            "combinatoric_sim_yield_per_photon": (
                float(comb_template[pt].sum() / comb_normalization[pt])
                if subtract_unmatched_recoil and comb_normalization[pt] > 0 else 0.0
            ),
            "combinatoric_scaled_recoil": float(comb_subtracted[pt].sum()),
            "purity_corrected_recoil_before_combinatoric": float(
                (corrected[pt] + comb_subtracted[pt]).sum()
            ),
        })
    return {
        "corrected": corrected,
        "variance": variance,
        "photons": photons,
        "solutions": solutions,
        "combinatoric_subtracted": comb_subtracted,
        "combinatoric_variance": comb_variance,
        "negative_input_sum": float(np.abs(corrected[corrected < 0]).sum()),
        "purity_strategy": purity_strategy,
    }


def correct_integrated_measured(data: Mapping[str, np.ndarray], response: CorrectionInputs) -> dict[str, Any]:
    """Purity-correct the integrated 15--35 GeV recoil spectrum.

    The headline observable is integrated over the declared photon range.  A
    single leakage-aware ABCD solve avoids imposing an unstable three-bin
    purity decomposition on the sparse high-pT AuAu sidebands.  Region C is a
    per-tag fake-photon recoil template, so its normalization is the inferred
    Region-A background divided by the *observed* Region-C tag count.  The
    detector unmatched-recoil component is not subtracted here as an absolute
    MC rate; it is applied later as a reconstructed-bin matching fraction.
    """
    counts = np.asarray(data["abcd_event_leading_counts_ptgamma"], dtype=float).sum(axis=1)
    spectra = np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float).sum(axis=1)
    spectra_sumw2 = np.asarray(data["inclusive_recoil_sumw2_ptgamma_xj"], dtype=float).sum(axis=1)
    leakage = np.asarray(response.leakage, dtype=float).sum(axis=1)
    solution = solve_leakage_abcd(counts, leakage)
    background_a = max(0.0, float(counts[0] - solution.signal_a))
    region_c_tags = float(counts[2])
    alpha = background_a / region_c_tags if region_c_tags > 0 else 0.0
    corrected = spectra[0] - alpha * spectra[1]
    variance = spectra_sumw2[0] + alpha * alpha * spectra_sumw2[1]
    negative = float(np.abs(corrected[corrected < 0]).sum())
    return {
        "corrected": corrected,
        "variance": variance,
        "photons": float(solution.signal_a),
        "solutions": [{
            "pt_low": float(PT_EDGES[0]),
            "pt_high": float(PT_EDGES[-1]),
            "signal_a": solution.signal_a,
            "purity": solution.purity,
            "background_a_over_observed_c": alpha,
            "factorization_residual": solution.residual,
            "method": solution.method,
        }],
        "region_a_recoil": spectra[0],
        "region_c_recoil": spectra[1],
        "region_c_scale": alpha,
        "negative_input_sum": negative,
    }


def draw_joint_data_toy(data: Mapping[str, np.ndarray], rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Fluctuate inclusive-pair sufficient statistics without a false partition.

    Multiple recoil jets may accompany one photon, so the xJ bins do not
    partition the event-leading ABCD counts.  Until an event-bootstrap
    covariance product is supplied, fluctuate the recorded count and pair
    sumw2 terms independently.  This is conservative and explicit; it never
    forces the inclusive jet yield below the photon count.
    """
    counts = np.asarray(data["abcd_event_leading_counts_ptgamma"], dtype=float)
    counts_var = np.asarray(data["abcd_event_leading_sumw2_ptgamma"], dtype=float)
    spectra = np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float)
    spectra_var = np.asarray(data["inclusive_recoil_sumw2_ptgamma_xj"], dtype=float)
    toy_counts = np.clip(
        rng.normal(counts, np.sqrt(np.clip(counts_var, 0.0, None))),
        0.0,
        None,
    )
    toy_spectra = np.clip(
        rng.normal(spectra, np.sqrt(np.clip(spectra_var, 0.0, None))),
        0.0,
        None,
    )
    return {
        **data,
        "abcd_event_leading_counts_ptgamma": toy_counts,
        "inclusive_recoil_spectra_ptgamma_xj": toy_spectra,
    }


def draw_combinatoric_template_toy(
    response: CorrectionInputs,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Fluctuate inclusive combinatoric pairs and photon normalization.

    The template is a per-photon inclusive jet yield and can exceed one when
    integrated over xJ.  It is therefore never treated as a Bernoulli
    partition of the reconstructed-photon population.
    """
    template = np.asarray(response.combinatoric_reco, dtype=float)
    template_var = np.asarray(response.combinatoric_reco_sumw2, dtype=float)
    denominator = np.asarray(response.combinatoric_normalization_reco, dtype=float)
    denominator_var = np.asarray(response.combinatoric_normalization_reco_sumw2, dtype=float)
    toy_template = np.clip(
        rng.normal(template, np.sqrt(np.clip(template_var, 0.0, None))),
        0.0,
        None,
    )
    toy_denominator = np.clip(
        rng.normal(denominator, np.sqrt(np.clip(denominator_var, 0.0, None))),
        0.0,
        None,
    )
    return toy_template, toy_denominator

