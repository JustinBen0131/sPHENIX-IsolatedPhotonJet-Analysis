"""Whole-event bootstrap of the data from an event-contribution table.

A data event can put one photon into the ABCD counts and several jets into
several xJ bins, so the photon normalisation and the recoil bins are
correlated through the events that fill them. Per-fill histogram sumw2 cannot
carry that. The event-contribution table keeps, for every contributing event,
each of its photon and jet rows, so the data can be resampled event by event:

  * one multiplier is drawn per event and applied to every photon row and
    every recoil pair of that event together;
  * the correction (leakage-aware ABCD, Region-C restoration, combinatoric
    subtraction) is rerun on the resampled population with the central
    response, exactly as for the central value;
  * a toy is rejected when the correction raises, or when the ABCD solution
    takes a different branch in any photon-pT bin (another estimator, not a
    fluctuation).

Primary model: an independent Poisson(1) multiplier per event (independent
events, counting at fixed exposure). Cross-check: multinomial resampling of the
same events with replacement at a fixed event count.

Table layout (bin numbers are ROOT bin numbers, 0 is the underflow)::

    event__weight, event__cent_slice
    photon__event_index, photon__region_ABCD, photon__pt_bin
    jet__event_index, jet__region_A0_C1, jet__pt_bin, jet__xj_bin, jet__thresholds_passed
    pt_edges_gev, xj_edges, jet_thresholds_gev        (native grids)

A jet row enters jet threshold index t when ``jet__thresholds_passed > t``.
Tables that also carry the jets at every Delta phi add ``jet__window`` (1 for
a recoil jet, Delta phi >= 7 pi / 8), ``jet__dphi`` and ``jet__far`` (1 for
|Delta eta| > 0.8): the recoil spectra take the window jets, and the
Delta-phi resolved counts feed the data-driven combinatoric estimate
(combinatoric.py), which a table built with ``combinatoric`` settings adds to
its measured statistics as ``recoil_background``, recomputed for every
resampled population.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping, Sequence

import numpy as np

from . import combinatoric as _combinatoric
from . import corrections as _corrections
from . import statistics as _statistics
from . import unfolding as _unfolding

CorrectionInputs = _corrections.CorrectionInputs
correct_measured = _corrections.correct_measured
AnalysisResponse = _unfolding.AnalysisResponse
unfold_pairs = _unfolding.unfold_pairs
unfold_photons = _unfolding.unfold_photons
central_68 = _statistics.central_68
reported_mask = _statistics.reported_mask
densities = _unfolding.densities

MODELS = ("poisson", "multinomial")


def _measurement_bins(root_bins: np.ndarray, native_edges: np.ndarray, target_edges: np.ndarray) -> np.ndarray:
    """Measurement bin holding the lower edge of each native ROOT bin; -1 outside the measurement window."""

    inside = (root_bins >= 1) & (root_bins <= native_edges.size - 1)
    low = native_edges[np.clip(root_bins - 1, 0, native_edges.size - 2)]
    index = np.searchsorted(target_edges, low, side="right") - 1
    return np.where(inside & (index >= 0) & (index < target_edges.size - 1), index, -1)


class EventTable:
    """The selected events of a data event-contribution table on the measurement (photon pT, xJ) grid."""

    def __init__(self, table: Mapping[str, np.ndarray], pt_edges: Sequence[float], *, centrality_slices: Sequence[int],
                 jet_threshold_gev: float, xj_edges: Sequence[float] | None = None,
                 combinatoric: Mapping[str, str] | None = None) -> None:
        self.pt_edges = np.asarray(pt_edges, dtype=float)
        native_pt = np.asarray(table["pt_edges_gev"], dtype=float)
        self.xj_edges = np.asarray(table["xj_edges"], dtype=float)
        if xj_edges is not None and not np.allclose(np.asarray(xj_edges, dtype=float), self.xj_edges, rtol=0, atol=1e-12):
            raise ValueError("the table xJ grid differs from the measurement xJ grid")
        for edge in self.pt_edges:
            if not np.any(np.isclose(native_pt, edge, rtol=0, atol=1e-9)):
                raise ValueError(f"measurement photon-pT boundary {edge} is not on the table grid")
        thresholds = np.asarray(table["jet_thresholds_gev"], dtype=float)
        match = np.flatnonzero(np.isclose(thresholds, jet_threshold_gev, rtol=0, atol=1e-9))
        if match.size != 1:
            raise ValueError(f"jet threshold {jet_threshold_gev} GeV is not one of the table thresholds {thresholds.tolist()}")
        threshold_index = int(match[0])
        self.n_pt, self.n_xj = self.pt_edges.size - 1, self.xj_edges.size - 1

        selected = np.isin(table["event__cent_slice"], tuple(centrality_slices))
        self.events = np.flatnonzero(selected)
        renumber = np.full(selected.size, -1)
        renumber[self.events] = np.arange(self.events.size)
        self.weight = np.asarray(table["event__weight"])[self.events]

        event = np.asarray(table["photon__event_index"])
        region = np.asarray(table["photon__region_ABCD"]).astype(int)
        pt = _measurement_bins(np.asarray(table["photon__pt_bin"]).astype(int), native_pt, self.pt_edges)
        keep = selected[event] & (pt >= 0)
        self.photon_event = renumber[event[keep]]
        self.photon_flat = region[keep] * self.n_pt + pt[keep]

        event = np.asarray(table["jet__event_index"])
        region = np.asarray(table["jet__region_A0_C1"]).astype(int)
        pt = _measurement_bins(np.asarray(table["jet__pt_bin"]).astype(int), native_pt, self.pt_edges)
        xj = np.asarray(table["jet__xj_bin"]).astype(int)
        any_angle = (selected[event] & (pt >= 0) & (np.asarray(table["jet__thresholds_passed"]) > threshold_index)
                     & (xj >= 1) & (xj <= self.n_xj))
        window = np.asarray(table["jet__window"]) == 1 if "jet__window" in table else np.ones(event.size, dtype=bool)
        keep = any_angle & window
        self.jet_event = renumber[event[keep]]
        self.jet_flat = (region[keep] * self.n_pt + pt[keep]) * self.n_xj + xj[keep] - 1

        # every jet at any Delta phi, on the Delta-phi edges of the combinatoric estimate: (region, pT, xJ, dphi, side)
        self.combinatoric = dict(combinatoric) if combinatoric else None
        self.n_dphi = _combinatoric.DPHI_EDGES.size - 1
        self.dphi_event = None
        if "jet__dphi" in table:
            dphi_bin = np.searchsorted(_combinatoric.DPHI_EDGES, np.asarray(table["jet__dphi"], dtype=float), side="right") - 1
            rows = any_angle & (dphi_bin >= 0) & (dphi_bin < self.n_dphi)
            self.dphi_event = renumber[event[rows]]
            cell = (region[rows] * self.n_pt + pt[rows]) * self.n_xj + xj[rows] - 1
            self.dphi_flat = (cell * self.n_dphi + dphi_bin[rows]) * 2 + np.asarray(table["jet__far"]).astype(int)[rows]
        elif self.combinatoric:
            raise ValueError("a data-driven combinatoric estimate needs the jets at every Delta phi (jet__dphi, jet__far)")

    @property
    def n_events(self) -> int:
        return int(self.events.size)

    def measured(self, multipliers: np.ndarray | None = None) -> dict[str, np.ndarray]:
        """ABCD counts (4, n_pt) and Region A/C recoil spectra (2, n_pt, n_xj) of the (resampled) events.

        The sumw2 are set to the counts, as for unit-weight Poisson fills; the toys use only the counts.
        """

        weight = self.weight if multipliers is None else self.weight * multipliers
        counts = np.bincount(self.photon_flat, weights=weight[self.photon_event], minlength=4 * self.n_pt).reshape(4, self.n_pt)
        spectra = np.bincount(self.jet_flat, weights=weight[self.jet_event],
                              minlength=2 * self.n_pt * self.n_xj).reshape(2, self.n_pt, self.n_xj)
        out = {"abcd_counts": counts, "abcd_counts_sumw2": counts.copy(),
               "recoil_spectra": spectra, "recoil_spectra_sumw2": spectra.copy()}
        if self.dphi_event is not None:
            dphi = np.bincount(self.dphi_flat, weights=weight[self.dphi_event],
                               minlength=2 * self.n_pt * self.n_xj * self.n_dphi * 2).reshape(2, self.n_pt, self.n_xj, self.n_dphi, 2)
            out.update(dphi_spectra=dphi, dphi_spectra_sumw2=dphi.copy())
            if self.combinatoric:
                estimate = _combinatoric.estimate(dphi, dphi, **self.combinatoric)
                out.update(recoil_background=estimate["counts"], recoil_background_sumw2=estimate["variance"])
        return out

    def combinatoric_estimate(self, multipliers: np.ndarray | None = None) -> dict[str, Any]:
        """The data-driven combinatoric estimate of the (resampled) events with its fits, for the record."""

        measured = self.measured(multipliers)
        if "dphi_spectra" not in measured:
            raise ValueError("this table has no jets at every Delta phi")
        return _combinatoric.estimate(measured["dphi_spectra"], measured["dphi_spectra_sumw2"], **(self.combinatoric or {}))

    def summary(self) -> dict[str, Any]:
        pairs = np.bincount(self.jet_event, minlength=self.n_events)
        return {"events": self.n_events, "photon_rows": int(self.photon_event.size), "jet_rows": int(self.jet_event.size),
                "events_with_photon_and_pair": int(np.intersect1d(self.photon_event, self.jet_event).size),
                "pairs_per_event": {int(n): int(np.sum(pairs == n)) for n in np.unique(pairs)}}


def draw_event_multipliers(rng: np.random.Generator, n_events: int, model: str) -> np.ndarray:
    """One multiplier per event: Poisson(1), or the multiplicities of a resample of n_events with replacement."""

    if model == "poisson":
        return rng.poisson(1.0, n_events).astype(float)
    if model == "multinomial":
        return np.bincount(rng.integers(0, n_events, n_events), minlength=n_events).astype(float)
    raise ValueError(f"unknown event-bootstrap model: {model}")


def solution_methods(correction: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(row["method"] for row in correction["solutions"])


def event_bootstrap_toys(table: EventTable, inputs: CorrectionInputs, rng: np.random.Generator, n_toys: int, *,
                         model: str, central_methods: Sequence[str], subtract_combinatoric: bool = True,
                         purity_strategy: str = "per_pt", purity_corrections=None
                         ) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Corrected photon yields (n, n_pt) and measured pair vectors (n, n_pt n_xj) of the accepted toys, and the rejections."""

    photons, corrected = [], []
    rejected = {"exception": 0, "abcd_branch_change": 0}
    for _ in range(n_toys):
        multipliers = draw_event_multipliers(rng, table.n_events, model)
        try:
            toy = correct_measured(table.measured(multipliers), inputs, subtract_combinatoric, pt_edges=table.pt_edges,
                                   purity_strategy=purity_strategy, purity_corrections=purity_corrections)
        except ValueError:
            rejected["exception"] += 1
            continue
        if solution_methods(toy) != tuple(central_methods):
            rejected["abcd_branch_change"] += 1
            continue
        photons.append(toy["photons"].copy())
        corrected.append(toy["corrected"].ravel().copy())
    return np.asarray(photons), np.asarray(corrected), rejected


# ----------------------------------------------------------------------------
# Spread of a toy ensemble
# ----------------------------------------------------------------------------

def finite_row_covariance(rows: np.ndarray) -> tuple[np.ndarray, int]:
    """Sample covariance (ddof 1) of the rows with every entry finite, and how many rows entered."""

    values = np.asarray(rows, dtype=float)
    values = values[np.all(np.isfinite(values), axis=1)]
    return np.atleast_2d(np.cov(values, rowvar=False, ddof=1)), len(values)


def correlation_matrix(covariance: np.ndarray) -> np.ndarray:
    sigma = np.sqrt(np.clip(np.diag(covariance), 0, None))
    scale = sigma[:, None] * sigma[None, :]
    return np.divide(covariance, scale, out=np.zeros_like(covariance), where=scale > 0)


def unfold_toy_rows(toys: np.ndarray, unfold: Callable[[np.ndarray], np.ndarray], size: int) -> np.ndarray:
    """Unfold every toy; a failed toy is a NaN row so that paired ensembles stay index-aligned."""

    out = np.full((len(toys), size), np.nan)
    for index, toy in enumerate(toys):
        try:
            out[index] = unfold(toy)
        except ValueError:
            pass
    return out


def per_photon_rows(pair_rows: np.ndarray, photon_rows: np.ndarray, xj_edges: np.ndarray) -> np.ndarray:
    """pT-integrated per-photon density of every toy (sum of pairs over pT / sum of photons / dxJ); NaN when undefined."""

    xj_edges = np.asarray(xj_edges, dtype=float)
    out = np.full((len(pair_rows), xj_edges.size - 1), np.nan)
    for index, (pairs, photons) in enumerate(zip(pair_rows, photon_rows)):
        try:
            out[index] = densities(pairs, np.asarray(photons, dtype=float), xj_edges, len(photons))[0]
        except ValueError:
            pass
    return out


def per_photon_spread(pair_toys: np.ndarray, photon_toys: np.ndarray, window: AnalysisResponse, *, pair_iterations: int,
                      photon_iterations: int) -> dict[str, Any]:
    """Covariance, standard deviation and central-68 half width of the per-photon density over a data-toy ensemble."""

    pairs = unfold_toy_rows(pair_toys, lambda toy: unfold_pairs(toy, window, pair_iterations)[0], window.xj_truth.size)
    photons = unfold_toy_rows(photon_toys, lambda toy: unfold_photons(toy, window, photon_iterations)[0], window.photon_truth.size)
    rows = per_photon_rows(pairs, photons, window.xj_edges)
    covariance, used = finite_row_covariance(rows)
    good = rows[np.all(np.isfinite(rows), axis=1)]
    return {"covariance_data": covariance, "toys_used": used, "sigma_data": np.sqrt(np.diag(covariance)),
            "half_width_68": central_68(good)[2]}


def pair_toy_spread(corrected_toys: np.ndarray, photon_toys: np.ndarray, window: AnalysisResponse, *, iterations: int,
                    reported: tuple[float, float]) -> dict[str, Any]:
    """Spread of the measured vector, of the unfolded pairs at one iteration count and of the photon yield over a data-toy ensemble.

    The correlation summary is taken over the reported xJ bins: the mean off-diagonal correlation within each photon-pT
    block, and the mean correlation between blocks i and j (key "i-j").
    """

    n_pt, n_xj = window.n_pt, window.n_xj
    mask = np.tile(reported_mask(window.xj_edges, *reported), n_pt)
    measured_cov, _ = finite_row_covariance(corrected_toys)
    unfolded = unfold_toy_rows(corrected_toys, lambda toy: unfold_pairs(toy, window, iterations)[0], window.xj_truth.size)
    unfolded_cov, _ = finite_row_covariance(unfolded)
    reported_cov = unfolded_cov[np.ix_(mask, mask)]
    corr = correlation_matrix(unfolded_cov)[np.ix_(mask, mask)]
    block = corr.shape[0] // n_pt
    off_diagonal = ~np.eye(block, dtype=bool)
    within = [float(corr[i * block:(i + 1) * block, i * block:(i + 1) * block][off_diagonal].mean()) for i in range(n_pt)]
    between = {f"{i}-{j}": float(corr[i * block:(i + 1) * block, j * block:(j + 1) * block].mean())
               for i in range(n_pt) for j in range(i + 1, n_pt)}
    return {"measured_vector_sigma": np.sqrt(np.diag(measured_cov)).reshape(n_pt, n_xj),
            "unfolded_pairs_sigma": np.sqrt(np.diag(unfolded_cov)).reshape(n_pt, n_xj),
            "unfolded_covariance_eigenvalues_reported_bins": np.linalg.eigvalsh(reported_cov),
            "unfolded_covariance_rank_reported_bins": int(np.linalg.matrix_rank(reported_cov)),
            "mean_correlation_within_pt_bin": within, "mean_correlation_between_pt_bins": between,
            "photon_yield_sigma": np.sqrt(np.diag(finite_row_covariance(photon_toys)[0]))}


__all__ = ["EventTable", "MODELS", "correlation_matrix", "draw_event_multipliers", "event_bootstrap_toys",
           "finite_row_covariance", "pair_toy_spread", "per_photon_rows", "per_photon_spread", "solution_methods",
           "unfold_toy_rows"]
