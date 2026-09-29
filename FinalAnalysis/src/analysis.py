#!/usr/bin/env python3
"""The per-photon xJgamma distribution of one system, with its data and response uncertainties.

    python -m FinalAnalysis.src.analysis --config FinalAnalysis/config/final_analysis.yaml --system auau \
        --data data/*.hists.root \
        --simulation auau_photon12 sim/photon12/*.response.root \
        --simulation auau_photon20 sim/photon20/*.response.root \
        --output results/auau.json

At production scale the file lists do not fit a command line: any argument may be @FILE, one argument per line
(argparse fromfile), e.g. --data @data_files.txt.

Inputs are TreeToHists outputs: FillPhotonJet files of the data (their event tables) and FillPhotonJetResponse files of
the simulation, one file per source, grouped by sample. A source is the unit the response bootstrap resamples. The
configuration holds every choice: measurement grids, centrality slices, jet threshold, reported xJ range, sample
normalisations, Au+Au centrality weights, iteration counts (or the frozen rule that selects them), prior, toy counts
and seeds.

The chain:
  1. data: the event tables of every data file, the selected events on the measurement grid (event_bootstrap.EventTable);
  2. simulation: per sample, one row per source on the response axes (TreeToHists package.simulation_sums), weighted by
     the sample's pb per generated event; the response, fakes, leakage and combinatoric template of a population are
     built together (closure.package_window);
  3. combinatoric subtraction and purity correction of the data (corrections.correct_measured). The combinatoric
     jets come from the simulated template (sim_template, the reference analysis) or from the data themselves
     (data_dphi, combinatoric.py: the flow-shaped Delta-phi sideband estimate, refitted for every data toy, with its
     closure on simulation from the response files' Delta-phi fills, simulation_dphi). Every recoil jet of a photon is
     counted, and the photon's own recoil jets take up part of the window, so the option occupancy: simulation scales the
     window extrapolation by the simulation's window occupancy (window_occupancy);
  4. data toys: one Poisson(1) multiplier per event of the event table (event_bootstrap), or independent Poisson counts
     (poisson); a toy whose ABCD solution changes branch in any photon-pT bin is rejected;
  5. iteration counts: configured, or selected by the frozen rule (``scan``): the joint k from the frozen scan of the
     simulated-prior unfolding (iteration.scan_iterations), the photon-yield k from the measured-covariance refold test
     (iteration.photon_refold_selection); both on independent-Poisson data toys and the source-occurrence response
     bootstrap, over k = 2..4. When no k passes, nothing is selected and the measurement stops;
  6. the data-driven prior: the simulated truth reweighted in xJ to fit the corrected data (prior.PriorFit); joint
     (photon pT, xJ) unfolding and photon-yield unfolding with those iteration counts; per-photon
     normalisation. Every data toy and every source-resampled response toy is refitted and unfolded, and the
     uncertainty is the central-68 half width of each ensemble, in quadrature (prior.reweighted_result);
  7. simulation closure (when configured): the per-photon estimator of step 6, fitted and unfolded with the same
     iteration counts, applied to simulated pseudo-data, full (train = test = every source) and half (A -> B, B -> A,
     sources split by identity, even rank A); the error is the spread of unfolded minus truth over source resamples
     (estimator_closure). The closure of the unfolding alone at the same k (closure.evaluate) is run beside it.

Output: one JSON with the result, every intermediate the uncertainty uses, the settings and the sha256 of every input
and of this code, and a CSV of the points. Existing outputs are never overwritten.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import yaml

from .internal import response as _contract
from .internal import inputs as _package

from .internal import closure as _closure
from .internal import combinatoric as _combinatoric
from .internal import corrections as _corrections
from .internal import event_bootstrap as _event_bootstrap
from .internal import iteration as _iteration
from .internal import prior as _prior
from .internal import resampling as _resampling
from .internal import unfolding as _unfolding

HERE = Path(__file__).resolve().parent

SourceSample = _closure.SourceSample
EventTable = _event_bootstrap.EventTable
DATA_TOY_MODELS = ("event_bootstrap", "poisson")
CODE = ("analysis.py", "prior.py", "iteration.py", "closure.py", "combinatoric.py", "resampling.py", "event_bootstrap.py",
        "corrections.py", "unfolding.py", "statistics.py")


def purity_corrections(corrections: Mapping[str, Any]):
    """corrections.PurityCorrections from a systems entry's corrections block, one value per photon-pT bin each:
    background_correlation (R_bkg), purity_scale (a factor on the ABCD signal), efficiency_scale (data/simulation photon
    ID and isolation efficiency); None when it names none of them."""

    names = {"background_correlation": "background_correlation", "scale": "purity_scale", "efficiency_scale": "efficiency_scale"}
    values = {field: corrections.get(key) for field, key in names.items()}
    if all(v is None for v in values.values()):
        return None
    return _corrections.PurityCorrections(**{field: None if v is None else tuple(map(float, v)) for field, v in values.items()})


@dataclass(frozen=True)
class Settings:
    """Every choice of one system's result."""

    pt_edges: tuple[float, ...]
    xj_edges: tuple[float, ...]
    reco_pt_edges: tuple[float, ...]
    truth_pt_edges: tuple[float, ...]
    reported_xj: tuple[float, float]
    jet_pt_min_gev: float
    pair_iterations: int | None      # None: selected by the frozen rule
    photon_iterations: int | None
    prior_form: str
    fake_model: str
    fit_range: str
    data_toys: int
    data_toy_model: str
    data_seed: int
    response_toys: int
    response_seed: int
    subtract_combinatoric: bool
    purity_strategy: str
    combinatoric_method: str = "sim_template"   # or "data_dphi"
    combinatoric_shape: str = "flow"            # data_dphi: "flow" (nominal) or "flat" (systematic)
    combinatoric_scope: str = "pt"
    combinatoric_occupancy: str = "none"        # data_dphi: "none", or "simulation" (window_occupancy)
    scan_data_toys: int | None = None       # toys of the frozen iteration rule
    scan_data_seed: int | None = None
    scan_response_toys: int | None = None
    scan_response_seed: int | None = None
    closure_toys: int | None = None         # simulation closure of the estimator; None: not run
    closure_seed: int | None = None
    unfolding_closure_toys: int | None = None  # closure of the unfolding alone at the result's k; None: not run
    unfolding_closure_seed: int | None = None
    prior_power: float = 1.0                # the unfolding starts from w^power (1 nominal; 0.5, 1.5: the ATLAS bracket)
    purity_corrections: Any = None          # corrections.PurityCorrections: R_bkg, signal scale, efficiency scale per photon-pT bin

    @classmethod
    def from_config(cls, config: Mapping[str, Any], system: str) -> "Settings":
        m, r, p, u = config["measurement"], config["response"], config["prior"], config["uncertainty"]
        s = config["systems"][system]
        k = lambda value: None if value == "scan" else int(value)  # noqa: E731
        scan = config.get("iteration_scan") or {}
        closure, unfolding_closure = config.get("closure") or {}, config.get("unfolding_closure") or {}
        combinatoric = s["corrections"].get("combinatoric") or {}
        settings = cls(pt_edges=tuple(map(float, m["photon_pt_edges_gev"])), xj_edges=tuple(map(float, m["xj_edges"])),
                       reco_pt_edges=tuple(map(float, r["reco_ptgamma_edges"])), truth_pt_edges=tuple(map(float, r["truth_ptgamma_edges"])),
                       reported_xj=tuple(map(float, m["reported_xj"])), jet_pt_min_gev=float(m["jet_pt_min_gev"]),
                       pair_iterations=k(s["iterations"]["pairs"]), photon_iterations=k(s["iterations"]["photons"]),
                       prior_form=str(p["form"]), fake_model=str(p["fake_model"]), fit_range=str(p["fit_range"]),
                       prior_power=float(p.get("power", 1.0)),
                       data_toys=int(u["data_toys"]), data_toy_model=str(s["data_toys"]["model"]), data_seed=int(s["data_toys"]["seed"]),
                       response_toys=int(u["response_toys"]), response_seed=int(u["response_seed"]),
                       subtract_combinatoric=bool(s["corrections"]["subtract_combinatoric"]),
                       purity_strategy=str(s["corrections"]["purity_strategy"]),
                       purity_corrections=purity_corrections(s["corrections"]),
                       combinatoric_method=str(combinatoric.get("method", "sim_template")),
                       combinatoric_shape=str(combinatoric.get("shape", "flow")),
                       combinatoric_scope=str(combinatoric.get("scope", "pt")),
                       combinatoric_occupancy=str(combinatoric.get("occupancy", "none")),
                       scan_data_toys=scan.get("data_toys"), scan_data_seed=scan.get("data_seed"),
                       scan_response_toys=scan.get("response_toys"), scan_response_seed=scan.get("response_seed"),
                       closure_toys=closure.get("toys"), closure_seed=(closure.get("seeds") or {}).get(system),
                       unfolding_closure_toys=unfolding_closure.get("toys"),
                       unfolding_closure_seed=(unfolding_closure.get("seeds") or {}).get(system))
        return settings

    def __post_init__(self) -> None:
        self.check()

    def check(self) -> None:
        if self.prior_form not in _prior.PRIOR_FORMS or self.fake_model not in _prior.FAKE_MODELS or self.fit_range not in _prior.FIT_RANGES:
            raise ValueError("unknown prior form, fake model or fit range")
        if not (math.isfinite(self.prior_power) and self.prior_power > 0):
            raise ValueError(f"the prior power must be finite and positive, got {self.prior_power}")
        if self.data_toy_model not in DATA_TOY_MODELS:
            raise ValueError(f"unknown data toy model {self.data_toy_model}")
        if self.combinatoric_method not in ("sim_template", "data_dphi"):
            raise ValueError(f"unknown combinatoric method {self.combinatoric_method}")
        if self.combinatoric_method == "data_dphi":
            if self.combinatoric_shape not in _combinatoric.SHAPES or self.combinatoric_scope not in _combinatoric.SCOPES:
                raise ValueError("unknown data-driven combinatoric shape or scope")
            if self.combinatoric_occupancy not in ("none", "simulation"):
                raise ValueError(f"unknown window occupancy {self.combinatoric_occupancy}")
            if self.data_toy_model != "event_bootstrap" or not self.subtract_combinatoric:
                raise ValueError("the data-driven combinatoric estimate is refitted per event-bootstrap toy and is always subtracted")
        if any(k is not None and k < 1 for k in (self.pair_iterations, self.photon_iterations)) or min(self.data_toys, self.response_toys) < 2:
            raise ValueError("iterations must be positive and each toy ensemble needs at least two toys")
        if None in (self.pair_iterations, self.photon_iterations):
            scan = (self.scan_data_toys, self.scan_data_seed, self.scan_response_toys, self.scan_response_seed)
            if None in scan or min(self.scan_data_toys, self.scan_response_toys) < 2:
                raise ValueError("the frozen iteration rule needs iteration_scan data_toys, data_seed, response_toys, response_seed")
        for toys, seed, name in ((self.closure_toys, self.closure_seed, "closure"),
                                 (self.unfolding_closure_toys, self.unfolding_closure_seed, "unfolding_closure")):
            if toys is not None and (toys < 2 or seed is None):
                raise ValueError(f"{name} needs at least two toys and a seed for this system")

    @property
    def template_subtraction(self) -> bool:
        """Whether the corrections subtract the simulated combinatoric template."""
        return self.subtract_combinatoric and self.combinatoric_method == "sim_template"

    @property
    def data_driven(self) -> dict | None:
        """The EventTable combinatoric settings of the data-driven estimate, or None."""
        return {"shape": self.combinatoric_shape, "scope": self.combinatoric_scope} if self.combinatoric_method == "data_dphi" else None

    def grids(self):
        return _contract.Grids(pt_edges=self.pt_edges, xj_edges=self.xj_edges, reco_pt_edges=self.reco_pt_edges,
                               truth_pt_edges=self.truth_pt_edges)


DataToys = Callable[[Any, tuple[str, ...]], tuple[np.ndarray, np.ndarray, dict[str, int]]]


def event_bootstrap_data_toys(table: EventTable, settings: Settings) -> DataToys:
    """Whole-event Poisson(1) bootstrap of the event table: (photon toys, pair toys, rejections) for correction inputs."""

    def draw(inputs, central_methods):
        return _event_bootstrap.event_bootstrap_toys(table, inputs, np.random.default_rng(settings.data_seed), settings.data_toys,
                                                     model="poisson", central_methods=central_methods,
                                                     subtract_combinatoric=settings.template_subtraction,
                                                     purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
    return draw


def poisson_data_toys(data: Mapping[str, np.ndarray], settings: Settings) -> DataToys:
    """Independent-Poisson toys of the ABCD counts and recoil spectra: (photon toys, pair toys, rejections)."""

    def draw(inputs, central_methods):
        toys = _resampling.poisson_data_toys(data, inputs, pt_edges=np.asarray(settings.pt_edges), rng=np.random.default_rng(settings.data_seed),
                                             toys=settings.data_toys, subtract_combinatoric=settings.template_subtraction,
                                             purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
        return toys.photons, toys.pairs, toys.rejected
    return draw


def select_iterations(samples: Mapping[str, SourceSample], data: Mapping[str, np.ndarray], problem, build: Callable,
                      settings: Settings, table: EventTable | None = None) -> dict[str, Any]:
    """The frozen iteration rule on the simulated prior: the joint k from the frozen scan, the photon-yield k from the
    measured-covariance refold test. Raises when either finds no passing k."""

    pt_edges = np.asarray(settings.pt_edges)
    options = dict(pt_edges=pt_edges, subtract_combinatoric=settings.template_subtraction, purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
    if settings.data_driven:  # the data-driven background is refitted per toy: whole-event bootstrap on the scan's seed
        central_methods = tuple(row["method"] for row in _corrections.correct_measured(
            data, problem.inputs, settings.template_subtraction, pt_edges=pt_edges, purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)["solutions"])
        photon_toys, pair_toys, rejected = _event_bootstrap.event_bootstrap_toys(
            table, problem.inputs, np.random.default_rng(settings.scan_data_seed), settings.scan_data_toys, model="poisson",
            central_methods=central_methods, subtract_combinatoric=settings.template_subtraction, purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
        toys = _resampling.DataToys(photon_toys, pair_toys, settings.scan_data_toys, rejected)
    else:
        toys = _resampling.poisson_data_toys(data, problem.inputs, rng=np.random.default_rng(settings.scan_data_seed),
                                             toys=settings.scan_data_toys, **options)
    responses = _resampling.source_occurrence_toys(_closure.sample_arrays(samples), build, data, problem.inputs,
                                                   rng=np.random.default_rng(settings.scan_response_seed),
                                                   toys=settings.scan_response_toys, **options)
    central = _corrections.correct_measured(data, problem.inputs, settings.template_subtraction, pt_edges=pt_edges,
                                            purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
    measured = {"1D": central["photons"], "2D": central["corrected"].ravel()}
    ensembles = {"1D": toys.photons, "2D": toys.pairs}
    scans = {d: _iteration.scan_iterations(problem.window, d, measured[d], ensembles[d], responses.problem(d),
                                           reported_xj=settings.reported_xj, data_toys_requested=toys.requested)
             for d in ("1D", "2D")}
    photon_rule = _iteration.photon_refold_selection(scans["1D"], toys.photons, responses.accepted("1D"))
    selected = {"pairs": scans["2D"].selected_k, "photons": photon_rule["selected_k"]}
    record = {"selected": selected, "rows": {d: scan.rows for d, scan in scans.items()},
              "photon_refold": {key: photon_rule[key] for key in ("eligible_k", "selected_k", "chi2_ndf1_data_plus_response")},
              "data_toys": {"requested": toys.requested, "accepted": int(len(toys.pairs)), "rejected": dict(toys.rejected)},
              "response_toys": {"requested": responses.requested, "rejected": dict(responses.rejected)},
              "gate_contract": _iteration.gate_contract()}
    if None in selected.values():
        raise ValueError(f"FAILED_ANALYSIS: the frozen iteration rule selects no k: {selected}")
    return record


def per_photon_views(pairs, photons, xj_edges) -> np.ndarray:
    """(n_pt + 1, n_xj): pairs / photons / dxJ in each photon-pT bin, then the sum of pairs / the sum of photons / dxJ."""

    widths = np.diff(np.asarray(xj_edges, dtype=float))
    photons = np.asarray(photons, dtype=float)
    pairs = np.asarray(pairs, dtype=float).reshape(photons.size, widths.size)
    return np.vstack([pairs / photons[:, None] / widths, pairs.sum(axis=0) / photons.sum() / widths])


def half_split(identities: Mapping[str, Sequence[tuple[str, str]]]) -> dict[str, np.ndarray]:
    """'A'/'B' per source row: within each sample the sources sorted by (identity, path), even rank A, odd rank B."""

    labels = {}
    for name, rows in identities.items():
        order = sorted(range(len(rows)), key=lambda i: tuple(rows[i]))
        label = np.empty(len(rows), dtype="<U1")
        label[order] = np.where(np.arange(len(rows)) % 2 == 0, "A", "B")
        labels[name] = label
    return labels


def estimator_closure(samples: Mapping[str, SourceSample], build: Callable, labels: Mapping[str, np.ndarray],
                      data_photons: np.ndarray, data_pair_toys: np.ndarray, settings: Settings, pair_iterations: int,
                      photon_iterations: int) -> dict[str, Any]:
    """Simulation closure of the per-photon estimator.

    Modes: full (train = test = every source, one shared resample) and half (A -> B and B -> A: disjoint halves,
    resampled independently). The test side's pseudo-data (closure.pseudo_data: the Tcomb chain for the pairs, response
    only for the photon yield) is unfolded with the train side by the estimator of the data: the prior fitted with the
    data covariance in simulation units (the data-toy covariance divided by the data-to-simulation photon scale), then
    joint and photon-yield unfolding with the result's iteration counts, everything per photon (per_photon_views).
    Error: the spread of unfolded minus truth over source resamples. ``labels`` hold 'A'/'B' per source row; a toy
    draws one resample per labelled sample, in label order, and uses those of ``samples``.
    """

    problem_for = _prior.prior_problem_factory(samples, build, reported=settings.reported_xj, fake_model=settings.fake_model,
                                               fit_range=settings.fit_range, power=settings.prior_power)
    central_problem = problem_for(None)
    scale = np.repeat(np.asarray(data_photons, dtype=float) / central_problem.sim_photons, central_problem.n_xj)
    covariance = np.cov(np.asarray(data_pair_toys, dtype=float), rowvar=False, ddof=1) / np.outer(scale, scale)
    per_source = _closure.sample_arrays(samples)

    def one_pass(w_train, w_test):
        train = problem_for(w_train)
        test = build(_resampling.occurrence_sums(per_source, w_test))
        sides = ((train.window, train.inputs), test)
        pairs = _closure.pseudo_data("tcomb_chain" if settings.template_subtraction else "response_only", "2D", *sides)
        photons = _closure.pseudo_data("response_only", "1D", *sides)
        theta, _, _, converged = train.fit(pairs, photons, covariance, settings.prior_form)
        if not converged:
            raise ValueError("the prior fit did not converge")
        unfolded_pairs = train.unfold(pairs, pair_iterations, train.prior_weights(theta))[0]
        unfolded_photons = _unfolding.unfold_photons(photons, train.window, photon_iterations)[0]
        edges = settings.xj_edges
        return {"truth": per_photon_views(test[0].xj_truth, test[0].photon_truth, edges), "reco": per_photon_views(pairs, photons, edges),
                "unfolded": per_photon_views(unfolded_pairs, unfolded_photons, edges), "theta": theta}

    rng = np.random.default_rng(settings.closure_seed)
    labels = {name: np.asarray(label) for name, label in labels.items()}
    ones = {name: np.ones(len(label)) for name, label in labels.items()}
    half = {h: {name: (label == h).astype(float) for name, label in labels.items()} for h in "AB"}
    modes = {"full": (ones, ones, True), "A_to_B": (half["A"], half["B"], False), "B_to_A": (half["B"], half["A"], False)}

    def draw(base):
        return {name: np.bincount(rng.choice(np.flatnonzero(b > 0), size=int((b > 0).sum())), minlength=b.size).astype(float)
                for name, b in base.items()}

    reported = central_problem.reported
    results = {}
    for mode, (base_train, base_test, shared) in modes.items():
        central = one_pass(base_train, base_test)
        differences, failed = [], 0
        for _ in range(settings.closure_toys):
            w_train = draw(base_train)
            w_test = w_train if shared else draw(base_test)
            try:
                toy = one_pass(w_train, w_test)
            except (ValueError, RuntimeError):  # an unconverged fit, or a resample whose boundary migration is not linear
                failed += 1
                continue
            differences.append(toy["unfolded"] - toy["truth"])
        sigma = np.std(np.asarray(differences), axis=0, ddof=1)
        pull = np.divide(central["unfolded"] - central["truth"], sigma, out=np.zeros_like(sigma), where=sigma > 0)
        results[mode] = {**central, "sigma_of_difference": sigma, "toys": len(differences), "toys_failed": failed,
                         "max_abs_pull_reported": [float(np.max(np.abs(row[reported]))) for row in pull],
                         "chi2_ndf_diagonal_reported": [float(np.mean(row[reported] ** 2)) for row in pull]}
    return {"modes": results, "views": "photon-pT bins, then their integral", "reported_xj_mask": reported,
            "half_split": {name: {h: int(np.sum(label == h)) for h in "AB"} for name, label in labels.items()},
            "toys_requested": settings.closure_toys, "seed": settings.closure_seed,
            "iterations": {"pairs": pair_iterations, "photons": photon_iterations}}


def unfolding_closure(samples: Mapping[str, SourceSample], build: Callable, labels: Mapping[str, np.ndarray],
                      settings: Settings, iterations: Sequence[int]) -> dict[str, Any]:
    """Closure of the unfolding alone (closure.evaluate: photon yield, pairs, and the pairs through the Tcomb chain) at
    the declared iteration counts, full and half, with one generator for both modes."""

    rng = np.random.default_rng(settings.unfolding_closure_seed)
    sizes = {name: sample.n_sources for name, sample in samples.items()}
    rows = {name: np.asarray(labels[name]) for name in samples}
    out = {}
    for mode in ("full", "half"):
        out.update(_closure.evaluate(samples, build, _closure.closure_split(mode, sizes, rows), tuple(iterations),
                                     settings.unfolding_closure_toys, rng, settings.reported_xj))
    return {"results": out, "exposure": _closure.exposure(samples, rows), "iterations": list(iterations),
            "toys_requested": settings.unfolding_closure_toys, "seed": settings.unfolding_closure_seed}


def gates(window, data: Mapping[str, np.ndarray], central: Mapping[str, np.ndarray], reported: np.ndarray) -> dict[str, Any]:
    """Whole-chain data/simulation gates per photon-pT bin over the reported xJ range.

    Jets per photon: the data's region-A recoil pairs per region-A photon (raw) against the simulation's reco pairs per reco
    photon; and the corrected data (purity, combinatoric) against the simulation's signal-like pairs (reco pairs without the
    combinatoric jets and the jets of unmatched photons) per truth-linked reco photon (matched, and matched outside the
    window), which is what the purity correction estimates. Mean xJ of the same two pairs. The simulation's detector fakes
    per matched pair and combinatoric share of its reco pairs, which expose a broken reco-truth jet match.
    """

    shown = np.asarray(reported, dtype=bool)
    edges = np.asarray(window.xj_edges, dtype=float)
    centres = 0.5 * (edges[:-1] + edges[1:])[shown]

    def ratio(numerator, denominator):
        numerator, denominator = np.asarray(numerator, dtype=float), np.asarray(denominator, dtype=float)
        return np.divide(numerator, denominator, out=np.full(numerator.shape, np.nan), where=denominator != 0)

    def per_photon(pairs, photons):
        return ratio(np.asarray(pairs, dtype=float)[:, shown].sum(axis=1), photons)

    def mean_xj(pairs):
        pairs = np.asarray(pairs, dtype=float)[:, shown]
        return ratio((pairs * centres).sum(axis=1), pairs.sum(axis=1))

    n_pt, n_xj = np.asarray(window.xj_reco).shape
    matched = np.asarray(window.xj_response, dtype=float).sum(axis=0).reshape(n_pt, n_xj)
    signal = np.asarray(window.xj_reco) - np.asarray(window.xj_combinatoric) - np.asarray(window.xj_fake_photon)
    jets = {"data_raw": per_photon(data["recoil_spectra"][0], np.asarray(data["abcd_counts"])[0]),
            "simulation_reco": per_photon(window.xj_reco, window.photon_reco),
            "data_corrected": per_photon(central["corrected"], central["photons"]),
            "simulation_signal": per_photon(signal, np.asarray(window.photon_response, dtype=float).sum(axis=0)
                                            + np.asarray(window.photon_boundary_fakes, dtype=float))}
    xj = {"data_corrected": mean_xj(central["corrected"]), "simulation_signal": mean_xj(signal)}
    return {"jets_per_photon": {**jets, "raw_data_over_simulation": ratio(jets["data_raw"], jets["simulation_reco"]),
                                "corrected_data_over_simulation": ratio(jets["data_corrected"], jets["simulation_signal"])},
            "mean_xj": {**xj, "data_over_simulation": ratio(xj["data_corrected"], xj["simulation_signal"])},
            "simulation": {"detector_fakes_over_matched": ratio(np.asarray(window.xj_detector_fakes)[:, shown].sum(axis=1),
                                                                matched[:, shown].sum(axis=1)),
                           "combinatoric_share": ratio(np.asarray(window.xj_combinatoric)[:, shown].sum(axis=1),
                                                       np.asarray(window.xj_reco)[:, shown].sum(axis=1))}}


def measure(samples: Mapping[str, SourceSample], data: Mapping[str, np.ndarray], data_toys: DataToys, settings: Settings,
            build: Callable | None = None, half_labels: Mapping[str, np.ndarray] | None = None,
            table: EventTable | None = None, simulation_dphi: Mapping[str, np.ndarray] | None = None,
            simulation_dphi_halves: Mapping[str, Mapping[str, np.ndarray]] | None = None,
            record: dict | None = None) -> dict[str, Any]:
    """The per-photon result of one system from its simulated samples, its measured statistics and its data toys.

    ``build`` maps per-sample sums to (response window, correction inputs); by default closure.package_window with the
    samples' weights. ``half_labels`` ('A'/'B' per source row of each sample, half_split) are needed by the closures.
    ``simulation_dphi`` (simulation_dphi) adds, with the data-driven combinatoric estimate, its simulation closure
    (combinatoric.simulation_closure); a closure that cannot be fitted is recorded with its error. With the window
    occupancy, ``simulation_dphi_halves`` (simulation_dphi_halves) adds its half closure: the factor of one half of the
    simulated sources applied to the other. ``record`` (a dict), when given, receives the objects the result is computed
    from (response window, central corrections with their intermediate spectra, accepted data toys, the unfoldings at
    the result's iterations) for the IAN figures; it changes nothing computed.
    """

    pt_edges = np.asarray(settings.pt_edges)
    if settings.data_driven and ("recoil_background" not in data or table is None):
        raise ValueError("the data-driven combinatoric estimate needs the event table built with its settings")
    occupancy = (table.combinatoric or {}).get("window_occupancy") if table is not None else None
    if settings.data_driven and (occupancy is not None) != (settings.combinatoric_occupancy == "simulation"):
        raise ValueError("the event table's window occupancy does not match the settings")
    if (settings.closure_toys or settings.unfolding_closure_toys) and half_labels is None:
        raise ValueError("the closures need the half split of the simulated sources")
    if build is None:
        build = _closure.package_window(settings.grids(), {name: sample.weight for name, sample in samples.items()})
    problem_for = _prior.prior_problem_factory(samples, build, reported=settings.reported_xj, fake_model=settings.fake_model,
                                               fit_range=settings.fit_range, power=settings.prior_power)
    problem = problem_for(None)
    rule = None
    pair_iterations, photon_iterations = settings.pair_iterations, settings.photon_iterations
    if None in (pair_iterations, photon_iterations):
        rule = select_iterations(samples, data, problem, build, settings, table)
        pair_iterations = rule["selected"]["pairs"] if pair_iterations is None else pair_iterations
        photon_iterations = rule["selected"]["photons"] if photon_iterations is None else photon_iterations
    central = _corrections.correct_measured(data, problem.inputs, settings.template_subtraction, pt_edges=pt_edges,
                                            purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
    methods = tuple(row["method"] for row in central["solutions"])
    photon_toys, pair_toys, rejected = data_toys(problem.inputs, methods)
    if len(pair_toys) < 2:
        raise ValueError("fewer than two accepted data toys")
    if record is not None:
        record.update(window=problem.window, central=central, data=data, pair_iterations=pair_iterations,
                      photon_iterations=photon_iterations, settings=settings)
    result = _prior.reweighted_result(problem, problem_for, {name: sample.n_sources for name, sample in samples.items()}, data,
                                      (photon_toys, pair_toys), form=settings.prior_form, pair_iterations=pair_iterations,
                                      photon_iterations=photon_iterations, pt_edges=pt_edges,
                                      response_toys=settings.response_toys, rng=np.random.default_rng(settings.response_seed),
                                      subtract_combinatoric=settings.template_subtraction, purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections,
                                      record=record)
    entry = result["entry"]
    combinatoric = None
    if settings.data_driven:
        # the estimate with its fits, and the simulated template it replaces as the cross-check (same photons)
        estimate = table.combinatoric_estimate()
        without = {key: value for key, value in data.items() if not key.startswith("recoil_background")}
        plain = _corrections.correct_measured(without, problem.inputs, False, pt_edges=pt_edges, purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
        template = _corrections.correct_measured(without, problem.inputs, True, pt_edges=pt_edges, purity_strategy=settings.purity_strategy, purity_corrections=settings.purity_corrections)
        combinatoric = {"method": "data_dphi", "shape": settings.combinatoric_shape, "scope": settings.combinatoric_scope,
                        "window_occupancy": occupancy,
                        "counts_region_A_C": estimate["counts"], "variance_region_A_C": estimate["variance"], "fits": estimate["fits"],
                        "purity_corrected_background": plain["corrected"] - central["corrected"],
                        "sim_template_cross_check": template["combinatoric_subtracted"]}
        if simulation_dphi is not None:
            try:
                combinatoric["simulation_closure"] = _combinatoric.simulation_closure(
                    simulation_dphi["dphi_truthjet"], simulation_dphi["dphi_combinatoric"],
                    simulation_dphi["dphi_truthjet_sumw2"], simulation_dphi["dphi_combinatoric_sumw2"],
                    shape=settings.combinatoric_shape, scope=settings.combinatoric_scope, window_occupancy=occupancy)
            except ValueError as error:
                combinatoric["simulation_closure"] = {"error": str(error)}
            if occupancy is not None and simulation_dphi_halves:
                halves = {}
                for train, test in (("A", "B"), ("B", "A")):
                    try:
                        factor = _combinatoric.window_occupancy(simulation_dphi_halves[train]["dphi_combinatoric"],
                                                                simulation_dphi_halves[train]["dphi_combinatoric_sumw2"])["factor"]
                        side = simulation_dphi_halves[test]
                        halves[f"{train}_to_{test}"] = _combinatoric.simulation_closure(
                            side["dphi_truthjet"], side["dphi_combinatoric"], side["dphi_truthjet_sumw2"], side["dphi_combinatoric_sumw2"],
                            shape=settings.combinatoric_shape, scope=settings.combinatoric_scope, window_occupancy=factor)
                    except ValueError as error:
                        halves[f"{train}_to_{test}"] = {"error": str(error)}
                combinatoric["simulation_closure_occupancy_halves"] = halves
    closure = None
    if settings.closure_toys:
        closure = estimator_closure(samples, build, half_labels, central["photons"], pair_toys, settings, pair_iterations,
                                    photon_iterations)
    unfolding = None
    if settings.unfolding_closure_toys:
        unfolding = unfolding_closure(samples, build, half_labels, settings, sorted({pair_iterations, photon_iterations}))
    return {"per_photon": entry["central"], "sigma_total": entry["sigma_total"], "entry": entry,
            "closure": closure, "unfolding_closure": unfolding, "combinatoric": combinatoric,
            "gates": gates(problem.window, data, central, problem.reported),
            "iterations": {"pairs": pair_iterations, "photons": photon_iterations,
                           "source": "configured" if rule is None else "frozen rule"}, "iteration_rule": rule,
            "reported_xj_mask": problem.reported, "purity": central["solutions"], "corrected_ptgamma_xj": central["corrected"],
            "purity_corrected_photons_ptgamma": central["photons"], "data_toys_requested": settings.data_toys,
            "data_toys_accepted": int(len(pair_toys)), "data_toys_rejected": dict(rejected),
            "samples": {name: {"sources": sample.n_sources, "weight_pb_per_event": sample.weight} for name, sample in samples.items()}}


# ----------------------------------------------------------------------------
# Inputs from TreeToHists outputs
# ----------------------------------------------------------------------------

def _text(path: Path, key: str) -> str:
    import uproot
    with uproot.open(path) as f:
        return str(f[key])


NOMINAL_VARIATION = {"jet_scale": 1.0, "jet_smear": 0.0, "photon_scale": 1.0, "photon_smear": 0.0, "seed": 0}


def simulation_variation(files: Mapping[str, Sequence[Path]]) -> dict[str, float]:
    """The reco energy scale and resolution variation the simulation files were filled with (FillPhotonJetResponse
    "variation"; a file written before it existed is nominal). The files of one measurement must share it."""

    import uproot
    found = set()
    for paths in files.values():
        for path in paths:
            with uproot.open(path) as f:
                keys = {key.split(";")[0] for key in f.keys()}
                recorded = json.loads(str(f["variation"])) if "variation" in keys else NOMINAL_VARIATION
            found.add(tuple(sorted((key, int(value) if key == "seed" else float(value)) for key, value in recorded.items())))
    if len(found) != 1:
        raise ValueError(f"the simulation files carry different variations: {sorted(found)}")
    return dict(found.pop())


def physical_system(config: Mapping[str, Any], system: str) -> str:
    """The collision system of a `systems` entry: its `system` field, else the key itself (pp, auau). A centrality class of
    Au+Au is an entry of its own, e.g. auau_cent20_50: {system: auau, class: cent20_50, centrality_slices: [4, ..., 9], ...},
    with its seeds under the same key."""
    if system not in config["systems"]:
        raise ValueError(f"no systems entry {system!r}; entries: {sorted(config['systems'])}")
    physical = config["systems"][system].get("system", system)
    if physical not in ("pp", "auau"):
        raise ValueError(f"systems.{system}: unknown collision system {physical!r}")
    return physical


def simulation_samples(files: Mapping[str, Sequence[Path]], config: Mapping[str, Any], system: str,
                       settings: Settings) -> tuple[dict[str, SourceSample], dict[str, list[tuple[str, str]]]]:
    """One SourceSample per configured sample of the system, a row per FillPhotonJetResponse file on the response axes,
    after checking the file's profile, and in Au+Au its centrality weights, against the configuration, and that every
    file shares one variation (simulation_variation); and the identity (production input sha256, path) of each row's
    source. A file holds exactly one source."""

    block = config["systems"][system]
    simulation_variation(files)
    if set(files) != set(block["samples"]):
        raise ValueError(f"{system} needs simulation for exactly {sorted(block['samples'])}, got {sorted(files)}")
    weights = config["response"].get("auau_centrality_weights") if physical_system(config, system) == "auau" else []
    samples, identities = {}, {}
    for name in block["samples"]:
        rows, identities[name] = [], []
        for path in _package.source_inputs(files[name]):
            identity = _package.response_identity(path, profile=config['profile'], centrality_weights=weights,
                normalization=dict(sample=config['samples'][name], reference=config.get('reference_cross_section_pb') or {}))
            if identity in identities[name]: raise ValueError('duplicate response source identity: ' + str(identity))
            identities[name].append(identity)
            rows.append(_package.simulation_sums(path, block['class'], settings.jet_pt_min_gev, settings.reco_pt_edges,
                                                 settings.truth_pt_edges, settings.xj_edges))
        if not rows:
            raise ValueError(f"no simulation files for {name}")
        arrays = {key: np.stack([row[key] for row in rows]) for key in rows[0]}
        samples[name] = SourceSample(arrays, _package.sample_weight(config["samples"][name], config.get("reference_cross_section_pb") or {}))
    return samples, identities


def simulation_dphi(files: Mapping[str, Sequence[Path]], config: Mapping[str, Any], system: str,
                    settings: Settings, keep: Mapping[str, Sequence[bool]] | None = None) -> dict[str, np.ndarray] | None:
    """The Delta-phi fills of the simulation (jets of matched photons with and without a truth jet nearby), summed over
    every file (or those ``keep`` marks, per sample in file order) with its sample's weight and projected onto the
    measurement photon-pT bins: (pT, xJ, Delta phi, side), the input of the simulation closure of the data-driven
    combinatoric estimate. None when a file predates the fills or no file is kept."""

    block = config["systems"][system]
    projection = _unfolding.window_projection(np.asarray(settings.reco_pt_edges, float), np.asarray(settings.pt_edges, float))
    total = None
    for name in block["samples"]:
        weight = _package.sample_weight(config["samples"][name], config.get("reference_cross_section_pb") or {})
        for index, path in enumerate(_package.source_inputs(files[name])):
            if keep is not None and not keep[name][index]:
                continue
            sums = _package.simulation_dphi(path, block["class"], settings.jet_pt_min_gev, settings.reco_pt_edges,
                                            settings.xj_edges)
            if sums is None:
                return None
            scaled = {key: (weight ** 2 if key.endswith("_sumw2") else weight) * value for key, value in sums.items()}
            total = scaled if total is None else {key: total[key] + scaled[key] for key in total}
    if total is None:
        return None
    return {key: np.einsum("pr,rxds->pxds", projection, value) for key, value in total.items()}


def simulation_dphi_halves(files: Mapping[str, Sequence[Path]], config: Mapping[str, Any], system: str, settings: Settings,
                           labels: Mapping[str, np.ndarray]) -> dict[str, dict[str, np.ndarray]] | None:
    """simulation_dphi of each half of the sources (half_split labels): the window occupancy of one half is closed on the
    other. None unless both halves have the fills."""

    halves = {h: simulation_dphi(files, config, system, settings, {name: np.asarray(label) == h for name, label in labels.items()})
              for h in "AB"}
    return halves if all(value is not None for value in halves.values()) else None


def window_occupancy(simulation_dphi: Mapping[str, np.ndarray] | None, settings: Settings) -> np.ndarray | None:
    """The window-occupancy factors of the data-driven estimate (combinatoric.window_occupancy of the simulation's
    combinatoric jets) when the settings ask for them, else None."""

    if settings.combinatoric_method != "data_dphi" or settings.combinatoric_occupancy == "none":
        return None
    if simulation_dphi is None:
        raise ValueError("the window occupancy needs the Delta-phi fills of the simulation (FillPhotonJetResponse)")
    return _combinatoric.window_occupancy(simulation_dphi["dphi_combinatoric"], simulation_dphi["dphi_combinatoric_sumw2"])["factor"]


def data_table(files: Sequence[Path], config: Mapping[str, Any], system: str, settings: Settings,
               occupancy: np.ndarray | None = None) -> EventTable:
    """The selected data events on the measurement grid, from the event tables of FillPhotonJet files; ``occupancy``
    (window_occupancy) goes to the data-driven estimate."""

    block = config["systems"][system]
    for path in files:
        if _text(path, "profile") != config["profile"]:
            raise ValueError(f"{path}: profile {_text(path, 'profile')} is not {config['profile']}")
    table = _package.event_table_arrays(files, settings.xj_edges, all_jets=settings.data_driven is not None,
                                        data_trigger=block.get('data_trigger'), runs=block.get('runs'))
    combinatoric = settings.data_driven
    if combinatoric is not None and occupancy is not None:
        combinatoric = dict(combinatoric, window_occupancy=np.asarray(occupancy, dtype=float))
    return EventTable(table, settings.pt_edges, centrality_slices=tuple(block["centrality_slices"]),
                      jet_threshold_gev=settings.jet_pt_min_gev, xj_edges=settings.xj_edges, combinatoric=combinatoric)


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonable(value):
    if isinstance(value, Mapping):
        return {str(key): jsonable(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None  # diagnostic non-finites are null
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    return value


def write_points(result: Mapping[str, Any], xj_edges: Sequence[float], path: Path) -> None:
    entry = result["entry"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["xj_low", "xj_high", "one_over_ngamma_dn_dxj", "sigma_total", "sigma_data", "sigma_response", "reported"])
        for i in range(len(xj_edges) - 1):
            writer.writerow([xj_edges[i], xj_edges[i + 1], entry["central"][i], entry["sigma_total"][i], entry["half_width_68_data"][i],
                             entry["half_width_68_response"][i], int(result["reported_xj_mask"][i])])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     fromfile_prefix_chars="@")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--system", required=True,
                        help="a key of the configuration's systems: pp, auau, or a class entry such as auau_cent20_50")
    parser.add_argument("--data", required=True, type=Path, nargs="+", help="FillPhotonJet outputs of the data")
    parser.add_argument("--simulation", required=True, nargs="+", action="append", metavar=("SAMPLE", "FILE"),
                        help="a sample name and its FillPhotonJetResponse outputs, one file per source; repeat per sample")
    parser.add_argument("--output", required=True, type=Path, help="result JSON; the CSV of the points is written beside it")
    args = parser.parse_args(argv)
    points = args.output.with_suffix(".csv")
    if args.output.exists() or points.exists():
        raise FileExistsError(args.output)

    output, result, settings = compute(args.config, args.system, args.data,
                                       {group[0]: [Path(p) for p in group[1:]] for group in args.simulation})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(jsonable(output), indent=1, allow_nan=False) + "\n", encoding="utf-8")
    write_points(result, settings.xj_edges, points)
    print(f"[measurement] {args.output}: {args.system}, {result['data_toys_accepted']} data toys, "
          f"{result['entry']['response_toys_used']} response toys")
    return 0


def compute(config_path: Path, system: str, data_paths: Sequence[Path], files: Mapping[str, Sequence[Path]],
            record: dict | None = None) -> tuple[dict[str, Any], dict[str, Any], Settings]:
    """The result document of one systems entry from TreeToHists outputs (what main writes), the raw result and its
    settings; ``record`` is passed to measure (the IAN figure package reads it)."""
    config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    collision_system = physical_system(config, system)
    settings = Settings.from_config(config, system)
    files = {name: [Path(p) for p in paths] for name, paths in files.items()}
    data_paths = [Path(p) for p in data_paths]
    samples, identities = simulation_samples(files, config, system, settings)
    labels = half_split(identities)
    dphi = simulation_dphi(files, config, system, settings)
    table = data_table(data_paths, config, system, settings, occupancy=window_occupancy(dphi, settings))
    data = table.measured()
    toys = (event_bootstrap_data_toys(table, settings) if settings.data_toy_model == "event_bootstrap"
            else poisson_data_toys(data, settings))
    result = measure(samples, data, toys, settings, half_labels=labels, table=table, simulation_dphi=dphi,
                     simulation_dphi_halves=simulation_dphi_halves(files, config, system, settings, labels), record=record)
    output = {"schema": "PhotonJetPerPhotonResultV1", "system": system, "collision_system": collision_system,
              "analysis": config.get("analysis"),
              "settings": asdict(settings), "config": {"path": str(config_path), "sha256": sha256(config_path)},
              "simulation_variation": simulation_variation(files),
              "inputs": {"data": [{"path": str(p), "sha256": sha256(p)} for p in data_paths],
                         "simulation": {name: [{"path": str(p), "sha256": sha256(p)} for p in paths] for name, paths in files.items()}},
              "data_events": table.summary(), "code_sha256": {name: sha256(HERE / (name if name == "analysis.py" else "internal/" + name)) for name in CODE}, **result}
    return output, result, settings


if __name__ == "__main__":
    raise SystemExit(main())
