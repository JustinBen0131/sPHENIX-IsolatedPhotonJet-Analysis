"""Simulation closure of the unfolding: full closure and independent half closure.

The response of one set of simulated source occurrences (train) unfolds the
pseudo-data of another set (test), and the result is compared with the test
truth in absolute units. Nothing is normalised to the test truth.

  * full closure: train and test are the same occurrences, in the central
    value and in every bootstrap toy;
  * half closure: two source-disjoint halves, both directions (A response on
    B pseudo-data, B response on A pseudo-data); train and test are resampled
    independently, and test pseudo-data and test truth move together.

Response, misses, fakes, leakage, the combinatoric template (Tcomb) and its
normalisation of one side always come from the same source occurrences,
centrally and in every toy. Every sample keeps its full-sample per-event
weight: a half is not rescaled to full exposure.

Problems (the iteration count is a declared diagnostic value):
  1D response_only   pseudo = matched reco photons + pT-boundary fakes.
  2D response_only   pseudo = matched reco pairs + detector and pT-boundary
                     fakes (checked equal to reco - fake-photon recoil - Tcomb).
  2D tcomb_chain     pseudo = test reco - test fake-photon recoil - train Tcomb
                     scaled by (test truth-matched Region-A photons / train
                     Tcomb normalisation). This is the data subtraction with
                     the corrected photon yield replaced by its simulated truth.
                     It is not a purity closure.

Statistics per problem: the bootstrap sigma of (unfolded - truth), pulls, the
generalized chi2/rank against the toy covariance on the reported bins, the
finite-toy inflation (n - 1) / (n - p - 2) of an inverse covariance estimated
from n toys in p dimensions (Hartlap et al. 2007; reported, never applied), and
the closure chi2 in the form of the iteration scan.

Source occurrences enter as per-sample accumulators with the source occurrence
as the leading axis, named as in the histogram package (``pair_response``,
``pair_misses``, ``truth_photons``, ...); the leakage counts sit on the
reconstructed classification grid as ``leakage_abcd_reco``. A response builder
maps per-sample sums to one side, ``(AnalysisResponse, CorrectionInputs)``;
``package_window`` is the histogram-package path through ``restrict``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from . import response as _contract

from . import corrections as _corrections
from . import resampling as _resampling
from . import statistics as _statistics
from . import unfolding as _unfolding

CorrectionInputs = _corrections.CorrectionInputs
AnalysisResponse = _unfolding.AnalysisResponse
restrict = _unfolding.restrict
unfold_pairs = _unfolding.unfold_pairs
unfold_photons = _unfolding.unfold_photons
chi2_ndf = _statistics.chi2_ndf
occurrence_sums = _resampling.occurrence_sums
occurrence_weights = _resampling.occurrence_weights
covariance_chi2_ndf = _statistics.covariance_chi2_ndf
reported_mask = _statistics.reported_mask
Package = _contract.Package

Side = tuple  # (AnalysisResponse, CorrectionInputs)
ResponseBuilder = Callable[[Mapping[str, Mapping[str, np.ndarray]]], Side]

PROBLEMS = (("1D", "response_only"), ("2D", "response_only"), ("2D", "tcomb_chain"))

RESPONSE_KEYS = (
    "photon_response", "photon_misses", "truth_photons", "photon_reco", "photon_fakes",
    "pair_response", "pair_misses", "pair_truth", "pair_reco",
    f"pair_fakes_{_contract.UNMATCHED_RECO}", f"pair_fakes_{_contract.COMBINATORIC}",
    f"pair_boundary_reco_{_contract.ASSIGNED_HARD_NONFIDUCIAL}",
)
SUMW2_KEYS = (f"pair_fakes_{_contract.COMBINATORIC}_sumw2",)


# ----------------------------------------------------------------------------
# Source occurrences
# ----------------------------------------------------------------------------

@dataclass
class SourceSample:
    """Accumulators of one production sample, one row per source occurrence, and its per-event weight."""

    arrays: dict[str, np.ndarray]
    weight: float

    @property
    def n_sources(self) -> int:
        return len(next(iter(self.arrays.values())))

    def take(self, rows: np.ndarray) -> "SourceSample":
        return SourceSample({key: values[rows] for key, values in self.arrays.items()}, self.weight)


def sample_arrays(samples: Mapping[str, SourceSample]) -> dict[str, Mapping[str, np.ndarray]]:
    """The per-occurrence accumulators of each sample, in the form resampling.occurrence_sums takes."""

    return {name: sample.arrays for name, sample in samples.items()}


def weighted_package(grids, sample_weights: Mapping[str, float], sums: Mapping[str, Mapping[str, np.ndarray]]):
    """The histogram package of weighted per-sample sums on the response grids, before ``restrict``: counts times the
    sample weight, sumw2 times its square. ``package_window`` restricts it; the IAN response figures draw it."""

    package = Package.empty(grids)
    arrays = package.arrays
    for name, sums_of_sample in sums.items():
        weight = sample_weights[name]
        for key in RESPONSE_KEYS:
            arrays[key] += np.asarray(sums_of_sample[key]).reshape(arrays[key].shape) * weight
        for key in SUMW2_KEYS:
            arrays[key] += np.asarray(sums_of_sample[key]).reshape(arrays[key].shape) * weight * weight
    return package


def package_window(grids, sample_weights: Mapping[str, float]) -> ResponseBuilder:
    """Response builder through the histogram package: weighted per-sample sums, then ``restrict``.

    The Tcomb normalisation is the truth-matched Region-A photon count (leakage A).
    """

    projection = _unfolding.window_projection(grids.reco_pt_edges, grids.pt_edges)

    def build(sums: Mapping[str, Mapping[str, np.ndarray]]) -> Side:
        package = weighted_package(grids, sample_weights, sums)
        leakage = np.zeros((4, grids.n_pt))
        leakage_sumw2 = np.zeros((4, grids.n_pt))
        for name, sums_of_sample in sums.items():
            weight = sample_weights[name]
            leakage += np.asarray(sums_of_sample["leakage_abcd_reco"]) @ projection.T * weight
            leakage_sumw2 += np.asarray(sums_of_sample["leakage_abcd_reco_sumw2"]) @ projection.T * weight * weight
        window = restrict(package)
        inputs = CorrectionInputs(leakage=leakage, leakage_sumw2=leakage_sumw2, combinatoric_reco=window.xj_combinatoric,
                                  combinatoric_reco_sumw2=window.xj_combinatoric_sumw2,
                                  combinatoric_normalization_reco=leakage[0].copy(),
                                  combinatoric_normalization_reco_sumw2=leakage_sumw2[0].copy())
        return window, inputs

    return build


# ----------------------------------------------------------------------------
# Closure
# ----------------------------------------------------------------------------

def pseudo_data(level: str, dimension: str, train: Side, test: Side) -> np.ndarray:
    """Measured vector of the test side for one closure problem."""

    (_, train_inputs), (test_window, test_inputs) = train, test
    if dimension == "1D":
        matched_plus_unfolding_fakes = test_window.photon_response.sum(axis=0) + test_window.photon_boundary_fakes
    else:
        matched_plus_unfolding_fakes = test_window.xj_response.sum(axis=0) + test_window.unfolding_fakes
    if level == "response_only":
        if dimension == "2D":  # the same vector written as the data-like subtraction
            subtracted = (test_window.xj_reco - test_window.xj_fake_photon - test_window.xj_combinatoric).ravel()
            scale = max(1.0, float(np.max(np.abs(subtracted))))
            if not np.allclose(subtracted, matched_plus_unfolding_fakes, rtol=0.0, atol=1e-8 * scale):
                raise ValueError("pseudo-data identity differs: reco - fake photon - combinatoric != matched + unfolding fakes")
        return matched_plus_unfolding_fakes
    if level == "tcomb_chain" and dimension == "2D":
        norm = train_inputs.combinatoric_normalization_reco
        if np.any(norm <= 0):
            raise ValueError("train Tcomb normalization is not positive")
        scaled = train_inputs.combinatoric_reco * (test_inputs.leakage[0] / norm)[:, None]
        return (test_window.xj_reco - test_window.xj_fake_photon - scaled).ravel()
    raise ValueError(f"level {level} is not defined for {dimension}")


def closure_residuals(train: Side, test: Side, iterations: Sequence[int]) -> dict[tuple[str, str, int], tuple[np.ndarray, np.ndarray]]:
    """(unfolded with the train response, test truth) for every problem and iteration count."""

    train_window, test_window = train[0], test[0]
    out = {}
    for dimension, level in PROBLEMS:
        measured = pseudo_data(level, dimension, train, test)
        truth = test_window.photon_truth if dimension == "1D" else test_window.xj_truth.reshape(-1)
        unfold = unfold_photons if dimension == "1D" else unfold_pairs
        for k in iterations:
            out[(dimension, level, k)] = (unfold(measured, train_window, k)[0], truth)
    return out


def closure_split(mode: str, sizes: Mapping[str, int], labels: Mapping[str, np.ndarray] | None = None) -> list:
    """(direction, train rows, test rows, shared resample) per direction. Half labels are 'A'/'B' per source row."""

    if mode == "full":
        everything = {name: np.arange(size) for name, size in sizes.items()}
        return [("full", everything, everything, True)]
    if mode == "half":
        half = {h: {name: np.flatnonzero(np.asarray(labels[name]) == h) for name in sizes} for h in "AB"}
        return [("A_to_B", half["A"], half["B"], False), ("B_to_A", half["B"], half["A"], False)]
    raise ValueError(f"unknown closure mode: {mode}")


def evaluate(samples: Mapping[str, SourceSample], build: ResponseBuilder, split: Sequence, iterations: Sequence[int],
             n_toys: int, rng: np.random.Generator, reported: tuple[float, float]) -> dict[str, dict[str, Any]]:
    """Closure results keyed ``direction|dimension|level|k<k>``; the source bootstrap draws in sample order."""

    results = {}
    for label, train_rows, test_rows, shared in split:
        train_sources = {name: sample.take(train_rows[name]) for name, sample in samples.items()}
        test_sources = {name: sample.take(test_rows[name]) for name, sample in samples.items()}
        test_central = build(occurrence_sums(sample_arrays(test_sources)))
        central = closure_residuals(build(occurrence_sums(sample_arrays(train_sources))), test_central, iterations)
        toys = {key: [] for key in central}
        failed = 0
        for _ in range(n_toys):
            w_train = occurrence_weights(rng, {name: len(train_rows[name]) for name in samples})
            w_test = w_train if shared else occurrence_weights(rng, {name: len(test_rows[name]) for name in samples})
            try:
                train = build(occurrence_sums(sample_arrays(train_sources), w_train))
                test = train if shared else build(occurrence_sums(sample_arrays(test_sources), w_test))
                toy = closure_residuals(train, test, iterations)
            except ValueError:
                failed += 1
                continue
            for key, (unfolded, truth) in toy.items():
                toys[key].append(unfolded - truth)
        window = test_central[0]
        pair_mask = np.tile(reported_mask(window.xj_edges, *reported), window.n_pt)
        for (dimension, level, k), (unfolded, truth) in central.items():
            d = unfolded - truth
            ensemble = np.asarray(toys[(dimension, level, k)], dtype=float)
            mask = np.ones(d.size, dtype=bool) if dimension == "1D" else pair_mask
            sigma = ensemble.std(axis=0, ddof=1) if len(ensemble) > 1 else np.full(d.size, np.nan)
            chi2, rank, null_norm = covariance_chi2_ndf(d[mask], ensemble[:, mask]) if len(ensemble) > 1 else (math.inf, 0, math.inf)
            closure_form = chi2_ndf(truth, unfolded, np.clip(truth + np.abs(unfolded), 1e-12, None))
            n, p = len(ensemble), rank
            inflation = (n - 1.0) / (n - p - 2.0) if n - p - 2 > 0 else math.inf
            ratio = np.divide(unfolded, truth, out=np.full(d.size, np.nan), where=truth > 0)
            pulls = np.divide(d, sigma, out=np.full(d.size, np.nan), where=sigma > 0)
            results[f"{label}|{dimension}|{level}|k{k}"] = {
                "direction": label, "dimension": dimension, "level": level, "k_diagnostic": k,
                "unfolded": unfolded.tolist(), "test_truth": truth.tolist(), "ratio_unfolded_over_truth": ratio.tolist(),
                "bootstrap_sigma_of_difference": sigma.tolist(), "pull": pulls.tolist(), "reported_mask": mask.tolist(),
                "chi2_ndf_bootstrap_covariance_reported_bins": chi2, "covariance_rank": rank, "null_residual_norm": null_norm,
                "finite_toy_chi2_inflation_expected": inflation,
                "chi2_ndf_divided_by_expected_inflation": chi2 / inflation if math.isfinite(inflation) else None,
                "chi2_ndf_frozen_mc_closure_form_all_bins": closure_form,
                "max_abs_pull_reported": float(np.nanmax(np.abs(pulls[mask]))) if np.any(np.isfinite(pulls[mask])) else None,
                "toys_requested": n_toys, "toys_accepted": int(len(ensemble)), "toys_failed": failed,
                "bootstrap": "train and test share one resample (same events)" if shared else
                             "train and test resampled independently; test pseudo-data and test truth move together"}
    return results


def exposure(samples: Mapping[str, SourceSample], labels: Mapping[str, np.ndarray]) -> dict[str, Any]:
    """Source count and truth-photon yield fraction of each half; halves keep the full-sample per-event weight."""

    out = {}
    for name, sample in samples.items():
        truth = sample.arrays["truth_photons"]
        halves = np.asarray(labels[name])
        total = float(truth.sum())
        out[name] = {h: {"sources": int(np.sum(halves == h)), "truth_photon_yield_fraction": float(truth[halves == h].sum() / total)}
                     for h in "AB"}
        out[name]["external_weight_per_event"] = sample.weight
    return out


__all__ = ["PROBLEMS", "ResponseBuilder", "SourceSample", "closure_residuals", "closure_split",
           "evaluate", "exposure", "package_window", "pseudo_data", "sample_arrays", "weighted_package"]
