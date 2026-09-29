"""Data-driven prior for the joint (photon pT, xJ) unfolding.

The simulated truth is reweighted at particle level by a smooth function of xJ
until the folded simulation describes the measured reconstructed spectrum, and
the unfolding then starts from that prior. P(reco | truth bin) is not touched:
one weight per truth xJ bin changes the prior and the pT-boundary migration
vector, nothing else.

    w(xJ) = exp(t0 + t1 (u - 0.8) + t2 (u - 0.8)^2),  u = truth xJ bin centre,

held constant outside the reported range, one w for every photon-pT bin
(``exp1``: t2 = 0; ``exp2``). t comes from a generalized least-squares fit
(Levenberg-Marquardt) of the folded simulation, scaled to the measured photon
yield of each reconstructed pT bin, to the measured corrected pair spectrum,
with the covariance of the data-toy ensemble. The estimator includes the fit:
every data toy and every response toy is refitted.

Prior bracket (ATLAS, arXiv:1809.07280 §6): the unfolding starts from w^power, with the fit itself unchanged; power 1
is the nominal, 0.5 and 1.5 the bracket (sqrt(w) and w^(3/2)). Every refitted toy applies the same power.

Fit range: ``reported`` uses the reported reconstructed bins of every pT bin;
``all`` every measured bin with a defined uncertainty, as the unfolding does.

Detector fakes under a reweighted truth: ``fixed`` keeps the simulated
detector fakes as an absolute floor; ``scaled`` keeps the simulated fake
fraction of every reconstructed bin (a multiplicative purity correction, for
fakes that follow the local reconstructed signal).

Uncertainty: the central-68 half width of the refitted toys, data and
response in quadrature. The covariance drops a toy when any reported bin lies
beyond five robust sigma (1.4826 MAD) of the ensemble median, a
result-independent trim for the few toys in which D'Agostini with signed inputs
diverges; for a Gaussian ensemble it removes about one toy in a million.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares

from . import closure as _closure
from . import corrections as _corrections
from . import resampling as _resampling
from . import statistics as _statistics
from . import unfolding as _unfolding

correct_measured = _corrections.correct_measured
iterative_bayes = _unfolding.iterative_bayes
unfold_pairs = _unfolding.unfold_pairs
unfold_photons = _unfolding.unfold_photons
central_68 = _statistics.central_68
reported_mask = _statistics.reported_mask
densities = _unfolding.densities
SourceSample = _closure.SourceSample
sample_arrays = _closure.sample_arrays
occurrence_sums = _resampling.occurrence_sums
occurrence_weights = _resampling.occurrence_weights

PRIOR_FORMS = {"exp1": 2, "exp2": 3}   # parameters of each form
PIVOT_XJ = 0.8                         # u - 0.8 in w(xJ)
MAD_TO_SIGMA = 1.4826
TRIM_ROBUST_SIGMA = 5.0
FAKE_MODELS = ("fixed", "scaled")
FIT_RANGES = ("reported", "all")
DIAGNOSTIC_ITERATIONS = (1, 2, 3, 4, 6, 8)


def scale_truth_xj(sums: Mapping[str, np.ndarray], weights: np.ndarray) -> dict[str, np.ndarray]:
    """Weight every truth xJ bin of one sample's sums (all truth pT bins, guard bands included); reco is recomputed."""

    out = dict(sums)
    response = sums["pair_response"]
    tiled = np.tile(weights, response.shape[0] // len(weights))
    scaled = response * tiled[:, None]
    out["pair_response"] = scaled
    out["pair_misses"] = sums["pair_misses"] * tiled
    out["pair_truth"] = sums["pair_truth"] * tiled
    out["pair_reco"] = sums["pair_reco"] + (scaled - response).sum(axis=0)
    return out


def per_photon(pairs: np.ndarray, photons: np.ndarray, xj_edges: np.ndarray) -> np.ndarray:
    """pT-integrated per-photon density: pairs summed over pT / photons summed over pT / dxJ."""

    photons = np.asarray(photons, dtype=float)
    return densities(pairs, photons, xj_edges, len(photons))[0]


def robust_covariance(rows: np.ndarray, reported: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, int]:
    """(trimmed covariance, central-68 half width, toys kept by the trim, finite toys) of a toy ensemble."""

    values = np.asarray(rows, dtype=float)
    values = values[np.all(np.isfinite(values), axis=1)]
    median = np.median(values, axis=0)
    mad = MAD_TO_SIGMA * np.median(np.abs(values - median), axis=0)
    keep = np.all(np.abs(values - median)[:, reported] <= TRIM_ROBUST_SIGMA * np.maximum(mad[reported], 1e-300), axis=1)
    half_width = central_68(values)[2]
    return np.atleast_2d(np.cov(values[keep], rowvar=False, ddof=1)), half_width, int(keep.sum()), int(len(values))


class PriorFit:
    """One response (a side of the response builder) with the linear pieces the truth reweighting needs."""

    def __init__(self, build: Callable, sums: Mapping[str, Mapping[str, np.ndarray]], *, reported: tuple[float, float],
                 fake_model: str = "fixed", fit_range: str = "reported", power: float = 1.0) -> None:
        if fake_model not in FAKE_MODELS or fit_range not in FIT_RANGES:
            raise ValueError(f"unknown fake model {fake_model} or fit range {fit_range}")
        if not (np.isfinite(power) and power > 0):
            raise ValueError(f"the prior power must be finite and positive, got {power}")
        self.fake_model, self.fit_range, self.power = fake_model, fit_range, float(power)
        self.window, self.inputs = build(sums)
        window = self.window
        self.n_pt, self.n_xj = window.n_pt, window.n_xj
        base = window.xj_boundary_fakes.ravel()
        self.boundary = np.zeros((self.n_xj, self.n_pt * self.n_xj))  # d(pT-boundary fakes) / d(truth xJ weight)
        for x in range(self.n_xj):
            weights = np.ones(self.n_xj)
            weights[x] = 2.0
            shifted, _ = build({name: scale_truth_xj(s, weights) for name, s in sums.items()})
            self.boundary[x] = shifted.xj_boundary_fakes.ravel() - base
        self.boundary_fixed = base - self.boundary.sum(axis=0)  # migration without a truth xJ bin in the support (expected ~0)
        if np.max(np.abs(self.boundary_fixed)) > 1e-6 * max(1.0, float(np.max(np.abs(base)))):
            raise RuntimeError("pT-boundary migration is not linear in the truth xJ weights")
        self.detector = window.xj_detector_fakes.ravel()
        nominal_signal = window.xj_response.sum(axis=0) + base  # matched plus pT-boundary migration, both carry a truth xJ
        self.fake_fraction = np.divide(self.detector, nominal_signal, out=np.zeros(self.n_pt * self.n_xj), where=nominal_signal > 0)
        self.detector_without_signal = np.where(nominal_signal > 0, 0.0, self.detector)
        self.sim_photons = window.photon_response.sum(axis=0) + window.photon_boundary_fakes
        centres = 0.5 * (window.xj_edges[:-1] + window.xj_edges[1:])
        self.reported = reported_mask(window.xj_edges, *reported)
        self.u = np.clip(centres, centres[self.reported][0], centres[self.reported][-1]) - PIVOT_XJ

    def weights(self, theta: Sequence[float]) -> np.ndarray:
        t = np.zeros(3)
        t[:len(theta)] = theta
        return np.exp(t[0] + t[1] * self.u + t[2] * self.u ** 2)

    def prior_weights(self, theta: Sequence[float]) -> np.ndarray:
        """The truth weights the unfolding starts from: the fitted w, raised to the prior power (1: w itself)."""
        w = self.weights(theta)
        return w if self.power == 1.0 else w ** self.power

    def signal(self, weights: np.ndarray) -> np.ndarray:
        return ((self.window.xj_response * np.tile(weights, self.n_pt)[:, None]).sum(axis=0) + weights @ self.boundary
                + self.boundary_fixed)

    def detector_fakes(self, weights: np.ndarray, model: str | None = None) -> np.ndarray:
        if (model or self.fake_model) == "fixed":
            return self.detector
        return self.fake_fraction * self.signal(weights) + self.detector_without_signal

    def fakes(self, weights: np.ndarray, model: str | None = None) -> np.ndarray:
        """Unfolding fake cause: detector fakes plus pT-boundary migration."""

        return self.detector_fakes(weights, model) + weights @ self.boundary + self.boundary_fixed

    def folded(self, weights: np.ndarray, model: str | None = None) -> np.ndarray:
        """Simulated reconstructed pairs in the measurement bins for truth weights (simulation units)."""

        return self.signal(weights) + self.detector_fakes(weights, model)

    def fit(self, measured: np.ndarray, photons: np.ndarray, covariance: np.ndarray, form: str) -> tuple[np.ndarray, float, int, bool]:
        """(theta, chi2, ndf, converged) of the generalized least-squares fit of the scaled folded simulation."""

        mask = np.tile(self.reported, self.n_pt)
        if self.fit_range == "all":
            mask = np.diag(covariance) > 0
        scale = np.repeat(np.asarray(photons, dtype=float) / self.sim_photons, self.n_xj)
        v = covariance[np.ix_(mask, mask)]
        v = v + np.eye(v.shape[0]) * 1e-12 * np.trace(v) / v.shape[0]
        chol = np.linalg.cholesky(v)
        residual = lambda theta: np.linalg.solve(chol, (measured - scale * self.folded(self.weights(theta)))[mask])
        solution = least_squares(residual, np.zeros(PRIOR_FORMS[form]), method="lm")
        return solution.x, float(np.sum(solution.fun ** 2)), int(mask.sum() - PRIOR_FORMS[form]), bool(solution.success)

    def unfold(self, measured: np.ndarray, iterations: int, weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Joint unfolding with the simulated prior (weights None) or the reweighted one."""

        if weights is None:
            return unfold_pairs(measured, self.window, iterations)
        window = self.window
        return iterative_bayes(measured, window.xj_response, window.xj_misses, iterations,
                               window.xj_truth.reshape(-1) * np.tile(weights, self.n_pt), fakes=self.fakes(weights))


def reweighted_result(problem: PriorFit, problem_for: Callable[[Mapping[str, np.ndarray]], PriorFit], sizes: Mapping[str, int],
                      data: Mapping[str, np.ndarray], data_toys: tuple[np.ndarray, np.ndarray], *, form: str,
                      pair_iterations: int, photon_iterations: int, pt_edges: np.ndarray, response_toys: int,
                      rng: np.random.Generator, forms: Sequence[str] = tuple(PRIOR_FORMS),
                      iterations: Sequence[int] = DIAGNOSTIC_ITERATIONS, subtract_combinatoric: bool = True,
                      purity_strategy: str = "per_pt", purity_corrections=None,
                      record: dict | None = None) -> dict[str, Any]:
    """The per-photon result with the data-driven prior of ``form``, and its data-toy and response-toy spread.

    ``problem`` is the central response; ``problem_for(multiplicities)`` the response of a source resample of
    the samples of ``sizes``. ``data_toys`` are (photon yields, measured pair vectors) of the accepted data toys;
    their covariance weights every fit. Response toys resample source occurrences with ``rng``, correct the
    central data with the toy response and refit; a toy fails on a ValueError, a changed ABCD branch or an
    unconverged fit. ``purity_strategy`` and ``purity_corrections`` are the purity correction of corrections.correct_measured,
    the same for the central value and every toy. ``record`` (a dict), when given, receives the central joint unfolding at
    the result's iterations with its refolded spectrum, the photon unfolding at the result's and at every diagnostic
    iteration, and the accepted data toys' corrected and unfolded vectors, for the figures; it changes nothing computed.
    """

    window = problem.window
    central = correct_measured(data, problem.inputs, subtract_combinatoric, pt_edges=pt_edges, purity_strategy=purity_strategy,
                               purity_corrections=purity_corrections)
    methods = tuple(row["method"] for row in central["solutions"])
    measured, photons = central["corrected"].ravel(), central["photons"]
    photon_toys, pair_toys = data_toys
    covariance = np.cov(pair_toys, rowvar=False, ddof=1)
    unfolded_photons = unfold_photons(photons, window, photon_iterations)[0]
    entry: dict[str, Any] = {"fits": {}, "per_photon_by_k": {}}
    for name in forms:
        theta, chi2, ndf, ok = problem.fit(measured, photons, covariance, name)
        entry["fits"][name] = {"theta": theta, "chi2": chi2, "ndf": ndf, "converged": ok, "truth_weight_by_xj_bin": problem.weights(theta),
                               "prior_power": problem.power, "prior_weight_by_xj_bin": problem.prior_weights(theta)}
    scale = np.repeat(photons / problem.sim_photons, problem.n_xj)
    mask = np.tile(problem.reported, problem.n_pt)
    v = covariance[np.ix_(mask, mask)]
    r0 = (measured - scale * problem.folded(np.ones(problem.n_xj)))[mask]
    entry["chi2_of_unweighted_sim_at_reconstructed_level"] = float(r0 @ np.linalg.solve(v, r0))
    for method in ("nominal_prior",) + tuple(forms):
        weights = None if method == "nominal_prior" else problem.prior_weights(entry["fits"][method]["theta"])
        entry["per_photon_by_k"][method] = {k: per_photon(problem.unfold(measured, k, weights)[0], unfolded_photons, window.xj_edges)
                                            for k in iterations}
    if record is not None:
        pairs, refolded = problem.unfold(measured, pair_iterations, problem.prior_weights(entry["fits"][form]["theta"]))
        record.update(unfolded_pairs=pairs, refolded_pairs=refolded, unfolded_photons=unfolded_photons,
                      unfolded_photons_by_k={k: unfold_photons(photons, window, k)[0] for k in iterations},
                      toy_corrected_pairs=[], toy_corrected_photons=[], toy_unfolded_pairs=[], toy_unfolded_photons=[])

    y_data, theta_data = [], []
    for toy_photons, toy_pairs in zip(photon_toys, pair_toys):
        theta, _, _, ok = problem.fit(toy_pairs, toy_photons, covariance, form)
        theta_data.append(theta)
        row = np.full(problem.n_xj, np.nan)
        if ok:
            pairs = problem.unfold(toy_pairs, pair_iterations, problem.prior_weights(theta))[0]
            toy_unfolded_photons = unfold_photons(toy_photons, window, photon_iterations)[0]
            try:
                row = per_photon(pairs, toy_unfolded_photons, window.xj_edges)
            except ValueError:
                pass  # a nonpositive unfolded photon yield has no per-photon density
            if record is not None:
                record["toy_corrected_pairs"].append(np.asarray(toy_pairs, dtype=float))
                record["toy_corrected_photons"].append(np.asarray(toy_photons, dtype=float))
                record["toy_unfolded_pairs"].append(pairs)
                record["toy_unfolded_photons"].append(toy_unfolded_photons)
        y_data.append(row)
    y_data = np.asarray(y_data)

    y_response, rejected = [], 0
    for _ in range(response_toys):
        multiplicities = occurrence_weights(rng, sizes)
        try:
            toy = problem_for(multiplicities)
            corrected = correct_measured(data, toy.inputs, subtract_combinatoric, pt_edges=pt_edges, purity_strategy=purity_strategy,
                               purity_corrections=purity_corrections)
            if tuple(row["method"] for row in corrected["solutions"]) != methods:
                raise ValueError("abcd branch change")
            theta, _, _, ok = toy.fit(corrected["corrected"].ravel(), corrected["photons"], covariance, form)
            if not ok:
                raise ValueError("fit")
            y_response.append(per_photon(toy.unfold(corrected["corrected"].ravel(), pair_iterations, toy.prior_weights(theta))[0],
                                         unfold_photons(corrected["photons"], toy.window, photon_iterations)[0], toy.window.xj_edges))
        except ValueError:
            rejected += 1
    y_response = np.asarray(y_response)

    result = entry["per_photon_by_k"][form][pair_iterations]
    cov_d, half_d, used_d, finite_d = robust_covariance(y_data, problem.reported)
    cov_r, half_r, used_r, finite_r = robust_covariance(y_response, problem.reported)
    truth_pairs = window.xj_truth.reshape(-1)
    entry.update(central=result, covariance_data=cov_d, covariance_response=cov_r, data_toys_used=used_d, data_toys_finite=finite_d,
                 response_toys_used=used_r, response_toys_finite=finite_r, half_width_68_data=half_d, half_width_68_response=half_r,
                 response_toys_rejected=rejected, sigma_total=np.sqrt(half_d ** 2 + half_r ** 2),
                 sigma_total_from_trimmed_covariance=np.sqrt(np.diag(cov_d + cov_r)),
                 theta_data_toy_mean=np.nanmean(theta_data, axis=0), theta_data_toy_sd=np.nanstd(theta_data, axis=0, ddof=1),
                 prior_per_photon=per_photon(truth_pairs, window.photon_truth, window.xj_edges),
                 reweighted_prior_per_photon=per_photon(truth_pairs * np.tile(problem.prior_weights(entry["fits"][form]["theta"]), problem.n_pt),
                                                        window.photon_truth, window.xj_edges))
    return {"entry": entry, "data_toy_rows": y_data, "response_toy_rows": y_response}


def prior_problem_factory(samples: Mapping[str, SourceSample], build: Callable, *, reported: tuple[float, float],
                          fake_model: str, fit_range: str, power: float = 1.0,
                          transform: Callable[[dict], dict] | None = None) -> Callable[[Mapping[str, np.ndarray] | None], PriorFit]:
    """PriorFit of the central sums (None) or of sums weighted by source multiplicities.

    ``transform`` optionally maps each sample's sums before the response is built (an alternative response
    hypothesis applied to the same source occurrences).
    """

    def problem_for(multiplicities: Mapping[str, np.ndarray] | None) -> PriorFit:
        sums = occurrence_sums(sample_arrays(samples), multiplicities)
        if transform is not None:
            sums = {name: transform(s) for name, s in sums.items()}
        return PriorFit(build, sums, reported=reported, fake_model=fake_model, fit_range=fit_range, power=power)

    return problem_for


__all__ = ["DIAGNOSTIC_ITERATIONS", "FAKE_MODELS", "FIT_RANGES", "PIVOT_XJ", "PRIOR_FORMS", "PriorFit", "per_photon",
           "prior_problem_factory", "reweighted_result", "robust_covariance", "scale_truth_xj"]
