"""Validation-only calibration and model-independent region classification.

Numerics are retained in numerics.py. BDT centrality fits use pooled quantiles;
isolation fits retain the historical constant-in-pT then linear-in-centrality
estimator. Neither estimator is rerun for test or full-population diagnostics.
"""

import numpy as np
from numerics import (
    weighted_threshold,
    threshold_uncertainty,
    constant_fit,
    linear_fit,
    HistogramData,
    quantile_from_histogram,
)

PT_EDGES = [15, 17, 19, 21, 23, 26, 35]
CENT_EDGES = list(range(0, 81, 5))
TARGETS = (0.7, 0.8, 0.9)
ID_REGIONS = {"INVALID": 0, "BELOW_FLOOR": 1, "NONTIGHT": 2, "EXCLUDED": 3, "TIGHT": 4}
ISO_REGIONS = {"INVALID": 0, "ISOLATED": 1, "GAP": 2, "NONISOLATED": 3}


def line(coefficients, x):
    return coefficients[0] + coefficients[1] * np.asarray(x, dtype=float)


def thresholds(package, et, centrality):
    x = np.asarray(centrality if package["axis"] == "centrality" else et)
    low, high = package["domain"]
    return {
        k: np.where((x >= low) & (x < high), line(v, x), np.nan)
        for k, v in package["curves"].items()
    }


def id_regions(score, tight, upper, lower=0.1):
    score, tight, upper, lower = np.broadcast_arrays(score, tight, upper, lower)
    good = (
        np.isfinite(score)
        & np.isfinite(tight)
        & np.isfinite(upper)
        & np.isfinite(lower)
        & (lower < upper)
        & (upper < tight)
    )
    r = np.zeros(score.shape, dtype=np.int32)
    r[good & (score <= lower)] = 1
    r[good & (score > lower) & (score < upper)] = 2
    r[good & (score >= upper) & (score <= tight)] = 3
    r[good & (score > tight)] = 4
    return r


def isolation_regions(value, isolated, nonisolated):
    value, isolated, nonisolated = np.broadcast_arrays(value, isolated, nonisolated)
    good = (
        np.isfinite(value)
        & np.isfinite(isolated)
        & np.isfinite(nonisolated)
        & (isolated <= nonisolated)
    )
    r = np.zeros(value.shape, dtype=np.int32)
    r[good] = 2  # strict boundaries put equality in the gap, even at zero gap
    r[good & (value < isolated)] = 1
    r[good & (value > nonisolated)] = 3
    return r


def point(values, weights, efficiency, isolation=False):
    if len(values) < 2:
        raise ValueError("calibration cell has fewer than two signal candidates")
    if isolation:
        edges = np.linspace(-20.0, 50.0, 701)
        h = np.histogram(values, bins=edges, weights=weights)[0]
        h2 = np.histogram(values, bins=edges, weights=weights**2)[0]
        under, over = values < edges[0], values >= edges[-1]
        # ROOT upper edge is overflow; numpy includes it in the final bin.
        at_end = values == edges[-1]
        h[-1] -= weights[at_end].sum()
        h2[-1] -= (weights[at_end] ** 2).sum()
        d = HistogramData(
            h,
            h2,
            edges,
            float(weights[under].sum()),
            float(weights[over].sum()),
            float((weights[under] ** 2).sum()),
            float((weights[over] ** 2).sum()),
            (),
        )
        p = quantile_from_histogram(d, efficiency)
    else:
        t = weighted_threshold(values, weights, efficiency)
        p = dict(
            threshold_gev=t,
            **threshold_uncertainty(
                values,
                weights,
                t,
                efficiency,
                density_window=0.01,
                score_histogram_bins=2000,
            )
        )
    return dict(
        threshold=float(p["threshold_gev"]),
        error=float(p["threshold_error_gev"]),
        n=len(values),
        sumw=float(weights.sum()),
        sumw2=float((weights**2).sum()),
    )


def fit_working_points(
    values, weights, et, centrality, signal, system, *, partition, isolation=False
):
    if partition != "validation":
        raise ValueError("working points may only be fitted on validation")
    arrays = [np.asarray(a) for a in (values, weights, et, centrality, signal)]
    values, weights, et, centrality, signal = arrays
    good = (
        signal.astype(bool)
        & np.isfinite(values)
        & np.isfinite(weights)
        & (weights > 0)
        & (et >= 15)
        & (et < 35)
    )
    if system == "auau":
        good &= np.isfinite(centrality) & (centrality >= 0) & (centrality < 80)
    cells, flat, pooled, curves, fit_records = [], [], [], {}, []
    cents = list(zip(CENT_EDGES[:-1], CENT_EDGES[1:])) if system == "auau" else [(0, 1)]
    for eff in TARGETS:
        for lo, hi in cents:
            region = good & (
                (centrality >= lo) & (centrality < hi) if system == "auau" else True
            )
            points = []
            for p0, p1 in zip(PT_EDGES[:-1], PT_EDGES[1:]):
                m = region & (et >= p0) & (et < p1)
                p = dict(
                    efficiency=eff,
                    cent_min=lo,
                    cent_max=hi,
                    pt_min=p0,
                    pt_max=p1,
                    **point(values[m], weights[m], eff, isolation)
                )
                cells.append(p)
                points.append(p)
            if system == "auau":
                f = constant_fit(
                    np.array([p["threshold"] for p in points]),
                    np.array([p["error"] for p in points]),
                )
                flat.append(
                    dict(
                        efficiency=eff,
                        axis=(lo + hi) / 2,
                        threshold=f["value_gev"],
                        error=f["error_gev"],
                        fit=f,
                    )
                )
                pooled.append(
                    dict(
                        efficiency=eff,
                        axis=(lo + hi) / 2,
                        **point(values[region], weights[region], eff, isolation)
                    )
                )
        if system == "auau":
            rows = [
                r for r in (flat if isolation else pooled) if r["efficiency"] == eff
            ]
        else:
            rows = [
                dict(axis=(r["pt_min"] + r["pt_max"]) / 2, **r)
                for r in cells
                if r["efficiency"] == eff
            ]
        f = linear_fit(
            np.array([r["axis"] for r in rows]),
            np.array([r["threshold"] for r in rows]),
            np.array([r["error"] for r in rows]),
        )
        name = ("ISO" if isolation else "WP") + str(round(100 * eff))
        curves[name] = [f["intercept_gev"], f["slope_gev_per_percent"]]
        fit_records.append(dict(name=name, inputs=rows, fit=f))
    p = dict(
        schema="PhotonIDWorkingPointsV1",
        system=system,
        axis="centrality" if system == "auau" else "et",
        domain=[0.0, 80.0] if system == "auau" else [15.0, 35.0],
        curves=curves,
        fit_form="linear",
        fit_population="validation",
        status="candidate",
        isolation=isolation,
        non_tight_lower=0.1,
        estimator=(
            "histogram_cdf_constant_then_linear"
            if isolation
            else "weighted_quantile_pooled_then_linear"
        ),
        tables=dict(
            cells=cells, constant_fits=flat, pooled=pooled, linear_fits=fit_records
        ),
    )
    validate_ordering(p)
    return p


def validate_ordering(p):
    x = np.asarray(p["domain"])  # linear curves: endpoints prove ordering throughout
    tag = "ISO" if p.get("isolation") else "WP"
    a, b, c = [line(p["curves"][tag + str(n)], x) for n in (70, 80, 90)]
    good = (
        (a < b) & (b < c)
        if p.get("isolation")
        else ((a > b) & (b > c) & (a <= 1) & (c >= 0) & (b > p["non_tight_lower"]))
    )
    if not np.all(np.isfinite([a, b, c])) or not np.all(good):
        raise ValueError(
            "fitted thresholds cross or leave their valid domain; review calibration, never clip silently"
        )


def efficiency_report(values, weights, et, centrality, signal, package):
    """Evaluate frozen curves only, including binomial effective-entry errors."""
    cuts = thresholds(package, et, centrality)
    rows = []
    for lo, hi in (
        zip(CENT_EDGES[:-1], CENT_EDGES[1:])
        if package["system"] == "auau"
        else zip(PT_EDGES[:-1], PT_EDGES[1:])
    ):
        axis = centrality if package["system"] == "auau" else et
        for name, t in cuts.items():
            m = (
                signal
                & np.isfinite(values)
                & np.isfinite(t)
                & np.isfinite(weights)
                & (weights > 0)
                & (axis >= lo)
                & (axis < hi)
            )
            if not np.any(m):
                rows.append(
                    dict(wp=name, low=lo, high=hi, efficiency=None, error=None, n=0)
                )
                continue
            passed = (values < t) if package.get("isolation") else (values > t)
            w = weights[m]
            e = float(weights[m & passed].sum() / w.sum())
            neff = float(w.sum() ** 2 / (w * w).sum())
            rows.append(
                dict(
                    wp=name,
                    low=lo,
                    high=hi,
                    efficiency=e,
                    error=float(np.sqrt(e * (1 - e) / neff)),
                    n=int(m.sum()),
                    effective_entries=neff,
                )
            )
    return rows
