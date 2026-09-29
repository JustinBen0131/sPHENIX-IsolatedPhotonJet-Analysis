"""Toy ensembles of the iteration scan: independent-Poisson DATA toys and the source-occurrence response bootstrap.

Ported from the validated collaboration implementation (``data_toys`` and ``response_toys`` of the frozen
iteration-scan driver, 2026-09-20). The corrections are ``corrections.correct_measured``.

DATA toys fluctuate the measured sufficient statistics as independent Poisson counts: the region A/B/C/D
event-leading photon counts and the region A/C recoil spectra. Each toy is corrected with the central
correction inputs (leakage, combinatoric template and its denominator).

Response toys resample whole source occurrences with replacement within each production slice, so the response,
misses, fakes, leakage, combinatoric template and its denominator of one toy come from one resampled simulated
population. The central measured statistics are corrected again with every toy's inputs.

A toy whose correction fails, or whose ABCD solution takes another branch than the central correction in any
photon-pT bin, is rejected and counts against the toy success of the scan: another branch is another estimator,
not a fluctuation.

The draw order of the numpy generator is part of the definition. A DATA toy draws the A/B/C/D counts, then the
A/C recoil spectra. A response toy draws one multiplicity vector per production slice, in slice order. Rejected
toys consume their draws.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .corrections import correct_measured

REJECTIONS = ("exception", "abcd_branch_change")


def abcd_branches(correction: Mapping[str, Any]) -> tuple[str, ...]:
    """The ABCD solution method of every photon-pT bin of one correction."""

    return tuple(row["method"] for row in correction["solutions"])


@dataclass
class DataToys:
    """Accepted DATA toys: corrected photon yields (n, n_pt) and pair spectra (n, n_pt * n_xj), index-aligned."""

    photons: np.ndarray
    pairs: np.ndarray
    requested: int
    rejected: dict[str, int]


def poisson_data_toys(
    data: Mapping[str, np.ndarray],
    inputs,
    *,
    pt_edges: np.ndarray,
    rng: np.random.Generator,
    toys: int,
    subtract_combinatoric: bool = True,
    purity_strategy: str = "per_pt",
    purity_corrections=None,
) -> DataToys:
    """Independent-Poisson toys of ``abcd_counts`` and ``recoil_spectra``, corrected with the central ``inputs``.

    A toy's sumw2 equals its counts. ``data`` carries the keys of ``corrections.correct_measured``.
    """

    options = dict(pt_edges=pt_edges, purity_strategy=purity_strategy, purity_corrections=purity_corrections)
    central = correct_measured(data, inputs, subtract_combinatoric, **options)
    branches = abcd_branches(central)
    photons: list[np.ndarray] = []
    pairs: list[np.ndarray] = []
    rejected = dict.fromkeys(REJECTIONS, 0)
    for _ in range(toys):
        counts = rng.poisson(data["abcd_counts"]).astype(float)
        spectra = rng.poisson(data["recoil_spectra"]).astype(float)
        toy = {"abcd_counts": counts, "abcd_counts_sumw2": counts.copy(),
               "recoil_spectra": spectra, "recoil_spectra_sumw2": spectra.copy()}
        try:
            correction = correct_measured(toy, inputs, subtract_combinatoric, **options)
        except ValueError:
            rejected["exception"] += 1
            continue
        if abcd_branches(correction) != branches:
            rejected["abcd_branch_change"] += 1
            continue
        photons.append(correction["photons"].copy())
        pairs.append(correction["corrected"].ravel().copy())
    return DataToys(np.asarray(photons), np.asarray(pairs), toys, rejected)


# ----------------------------------------------------------------------------
# Source-occurrence bootstrap of the simulated response
# ----------------------------------------------------------------------------

def occurrence_counts(per_source: Mapping[str, Mapping[str, np.ndarray]]) -> dict[str, int]:
    """Number of source occurrences of every slice; all arrays of a slice share the leading axis."""

    counts = {}
    for name, arrays in per_source.items():
        lengths = {len(values) for values in arrays.values()}
        if len(lengths) != 1:
            raise ValueError(f"{name}: per-occurrence arrays disagree on the number of occurrences: {sorted(lengths)}")
        counts[name] = lengths.pop()
    return counts


def occurrence_weights(rng: np.random.Generator, counts: Mapping[str, int]) -> dict[str, np.ndarray]:
    """Bootstrap multiplicity of every occurrence: n draws with replacement within each slice, slices in order."""

    return {name: np.bincount(rng.integers(0, n, size=n), minlength=n).astype(float) for name, n in counts.items()}


def occurrence_sums(
    per_source: Mapping[str, Mapping[str, np.ndarray]],
    weights: Mapping[str, np.ndarray] | None = None,
    keys: Sequence[str] | None = None,
) -> dict[str, dict[str, np.ndarray]]:
    """Per-slice totals of per-occurrence accumulators; ``weights`` are bootstrap multiplicities (None: each once)."""

    totals = {}
    for name, arrays in per_source.items():
        selected = arrays if keys is None else {key: arrays[key] for key in keys}
        totals[name] = {key: (values.sum(axis=0) if weights is None else np.tensordot(weights[name], values, axes=(0, 0)))
                        for key, values in selected.items()}
    return totals


@dataclass
class ResponseToys:
    """One entry per requested toy: its window and corrected photon yields and pairs, or None when rejected."""

    responses: list[Any]
    photons: list[np.ndarray | None]
    pairs: list[np.ndarray | None]
    requested: int
    rejected: dict[str, int]

    def problem(self, dimension: str) -> list[tuple[Any, np.ndarray] | None]:
        """(window, measured) of every toy for the photon ("1D") or pair ("2D") problem; None where rejected."""

        measured = {"1D": self.photons, "2D": self.pairs}[dimension]
        return [None if window is None else (window, values) for window, values in zip(self.responses, measured)]

    def accepted(self, dimension: str) -> np.ndarray:
        """Corrected measured vectors of the accepted toys."""

        measured = {"1D": self.photons, "2D": self.pairs}[dimension]
        return np.asarray([values for values in measured if values is not None])


def source_occurrence_toys(
    per_source: Mapping[str, Mapping[str, np.ndarray]],
    build: Callable[[dict[str, dict[str, np.ndarray]]], tuple[Any, Any]],
    data: Mapping[str, np.ndarray],
    central_inputs,
    *,
    pt_edges: np.ndarray,
    rng: np.random.Generator,
    toys: int,
    keys: Sequence[str] | None = None,
    subtract_combinatoric: bool = True,
    purity_strategy: str = "per_pt",
    purity_corrections=None,
) -> ResponseToys:
    """Response toys from whole source occurrences resampled within each production slice.

    ``per_source[slice][name]`` holds one accumulator per occurrence along the leading axis, slices in their
    production order. ``build`` turns the resampled per-slice totals into (analysis window, correction inputs);
    a ``ValueError`` there rejects the toy. The central ``data`` are corrected with every toy's inputs, and the
    ABCD branches are compared with those of the central correction with ``central_inputs``.
    """

    options = dict(pt_edges=pt_edges, purity_strategy=purity_strategy, purity_corrections=purity_corrections)
    branches = abcd_branches(correct_measured(data, central_inputs, subtract_combinatoric, **options))
    counts = occurrence_counts(per_source)
    result = ResponseToys([], [], [], toys, dict.fromkeys(REJECTIONS, 0))
    for _ in range(toys):
        sums = occurrence_sums(per_source, occurrence_weights(rng, counts), keys)
        try:
            window, inputs = build(sums)
            correction = correct_measured(data, inputs, subtract_combinatoric, **options)
        except ValueError:
            reason = "exception"
        else:
            reason = None if abcd_branches(correction) == branches else "abcd_branch_change"
        if reason is not None:
            result.rejected[reason] += 1
            result.responses.append(None); result.photons.append(None); result.pairs.append(None)
            continue
        result.responses.append(window)
        result.photons.append(correction["photons"].copy())
        result.pairs.append(correction["corrected"].ravel().copy())
    return result


__all__ = ["DataToys", "ResponseToys", "abcd_branches", "occurrence_counts", "occurrence_sums", "occurrence_weights",
           "poisson_data_toys", "source_occurrence_toys"]
