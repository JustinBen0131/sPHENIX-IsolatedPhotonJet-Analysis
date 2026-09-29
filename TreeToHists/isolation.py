"""The Au+Au isolation package: reconstructed-isolation cutoffs of truth-signal photons in embedded simulation.

    python -m TreeToHists.isolation derive population_*.root --out isolation_package.json
    python -m TreeToHists.isolation freeze isolation_package.json --out CANONICAL_ISOLATION_FREEZE.json
    python -m TreeToHists.isolation check CANONICAL_ISOLATION_FREEZE.json population_*.root --out isolation_check.json

Input: FillIsolationPopulation outputs. For every photon-pT group and 5% centrality slice (native bin b, slice (b - 1) // 5)
the reconstructed isolation of the population is histogrammed (700 bins from -20 to 50 GeV, unweighted, finite values),
and the cutoff where the cumulative spectrum reaches each efficiency (70, 80, 90%) is read inside its bin by linear
interpolation. Its error is the binomial error of the efficiency divided by the local density, and at least half a bin.
Per slice a constant is fitted across the pT groups; per efficiency a straight line across the slice centres:

    threshold(c) = intercept + slope * c        (GeV; c the centrality in percent)

The package holds, per cone (R = 0.3 and 0.4), the cutoff points, the constant fits and the lines, with the sha256 of
every input. ``freeze`` takes the R = 0.3 line at the nominal efficiency (80%) with the sideband gap (0.8 GeV) into the
hash-bound package the eventTree conversion applies to every Au+Au photon (pho_iso_class, Format/isolation_class.py).
``check`` measures a frozen line on new simulation: its realized efficiency per slice against the target, and the line
refitted there against the frozen one.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import uproot

NBINS, XMIN, XMAX = 700, -20.0, 50.0
EDGES = np.linspace(XMIN, XMAX, NBINS + 1)
PT_GROUPS = (((15, 17), ((15, 17),)), ((17, 19), ((17, 19),)), ((19, 21), ((19, 21),)), ((21, 23), ((21, 23),)),
             ((23, 26), ((23, 26),)), ((26, 35), ((26, 30), (30, 35))))
CENTRALITY_BINS = tuple((lo, lo + 5) for lo in range(0, 80, 5))
EFFICIENCIES = (0.70, 0.80, 0.90)
CONES = {"R30": "iso_r03", "R40": "iso_r04"}


def histogram(values: np.ndarray) -> dict:
    """Counts on EDGES of the finite values, with under- and overflow; variances equal the counts (unweighted)."""
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    counts, _ = np.histogram(finite, bins=EDGES)
    under, over = float(np.sum(finite < XMIN)), float(np.sum(finite >= XMAX))
    return dict(values=counts.astype(float), variances=counts.astype(float), edges=EDGES.copy(), underflow=under, overflow=over,
                underflow_variance=under, overflow_variance=over)


def quantile_from_histogram(h: dict, probability: float) -> dict:
    """The cutoff at ``probability`` of the cumulative spectrum, with its error and the counts behind it."""
    values = np.asarray(h["values"], dtype=float)
    negative_tolerance = 1e-10 * max(1.0, float(np.max(np.abs(values))))
    if np.any(values < -negative_tolerance):
        raise ValueError("isolation histogram has negative weighted bins; the quantile is undefined")
    values = np.maximum(values, 0.0)
    in_range_total = float(np.sum(values))
    underflow, overflow = max(float(h["underflow"]), 0.0), max(float(h["overflow"]), 0.0)
    total = in_range_total + underflow + overflow
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("isolation histogram has no positive in-range weight")
    cumulative = np.cumsum(values)
    absolute_target = probability * total
    if absolute_target <= underflow:
        raise ValueError(f"quantile {probability:.3f} lies in the underflow below {h['edges'][0]:g} GeV")
    if absolute_target > underflow + in_range_total:
        raise ValueError(f"quantile {probability:.3f} lies in the overflow above {h['edges'][-1]:g} GeV")
    target = absolute_target - underflow
    index = min(max(int(np.searchsorted(cumulative, target, side="left")), 0), len(values) - 1)
    previous = float(cumulative[index - 1]) if index > 0 else 0.0
    bin_weight = float(values[index])
    fraction = min(max(0.5 if bin_weight <= 0.0 else (target - previous) / bin_weight, 0.0), 1.0)
    width = float(h["edges"][index + 1] - h["edges"][index])
    threshold = float(h["edges"][index] + fraction * width)
    sumw2 = float(np.sum(np.maximum(h["variances"], 0.0)) + h["underflow_variance"] + h["overflow_variance"])
    effective_entries = total * total / sumw2 if sumw2 > 0.0 else total
    density = bin_weight / (total * width) if width > 0.0 else 0.0
    if effective_entries > 0.0 and density > 0.0:
        statistical_error = math.sqrt(probability * (1.0 - probability) / effective_entries) / density
    else:
        statistical_error = width
    return {"threshold_gev": threshold, "threshold_error_gev": float(max(0.5 * width, statistical_error)), "sumw": total,
            "sumw2": sumw2, "effective_entries": float(effective_entries), "underflow": h["underflow"], "overflow": h["overflow"],
            "flow_fraction": float((underflow + overflow) / max(total, 1e-12))}


def constant_fit(y: np.ndarray, yerr: np.ndarray) -> dict:
    valid = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0.0)
    if int(np.count_nonzero(valid)) < 2:
        raise ValueError("a constant fit needs at least two valid points")
    y, yerr = y[valid], yerr[valid]
    weights = 1.0 / np.square(yerr)
    value = float(np.sum(weights * y) / np.sum(weights))
    return {"value_gev": value, "error_gev": float(math.sqrt(1.0 / np.sum(weights))),
            "chi2": float(np.sum(np.square((y - value) / yerr))), "ndf": int(y.size - 1), "n_points": int(y.size)}


def linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> dict:
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0.0)
    if int(np.count_nonzero(valid)) < 3:
        raise ValueError("a linear fit needs at least three valid points")
    x, y, yerr = x[valid], y[valid], yerr[valid]
    design = np.column_stack([np.ones_like(x), x])
    weights = 1.0 / np.square(yerr)
    covariance = np.linalg.inv(design.T @ (weights[:, None] * design))
    parameters = covariance @ (design.T @ (weights * y))
    chi2 = float(np.sum(np.square((y - design @ parameters) / yerr)))
    return {"intercept_gev": float(parameters[0]), "slope_gev_per_percent": float(parameters[1]),
            "intercept_error_gev": float(math.sqrt(covariance[0, 0])), "slope_error_gev_per_percent": float(math.sqrt(covariance[1, 1])),
            "covariance": covariance.tolist(), "chi2": chi2, "ndf": int(y.size - 2), "n_points": int(y.size)}


def derive(et: np.ndarray, native_bin: np.ndarray, iso: np.ndarray, groups=PT_GROUPS) -> tuple[list, list, list]:
    """(cutoff points, constant fits per slice, lines per efficiency) of one cone."""
    et, native, iso = np.asarray(et, dtype=float), np.asarray(native_bin, dtype=np.int64), np.asarray(iso, dtype=float)
    slice_index = (native - 1) // 5
    points = []
    for k, (cent_min, cent_max) in enumerate(CENTRALITY_BINS):
        in_cent = (slice_index == k) & (native >= 1)
        for target_pt, source_bins in groups:
            in_pt = np.zeros(len(et), dtype=bool)
            for lo, hi in source_bins:
                in_pt |= (et >= lo) & (et < hi)
            selected = in_cent & in_pt
            for efficiency in EFFICIENCIES:
                points.append({"efficiency": efficiency, "cent_min": cent_min, "cent_max": cent_max, "cent_center": 0.5 * (cent_min + cent_max),
                               "pt_min": target_pt[0], "pt_max": target_pt[1], "pt_center": 0.5 * (target_pt[0] + target_pt[1]),
                               "photons": int(selected.sum()), **quantile_from_histogram(histogram(iso[selected]), efficiency)})
    flat = []
    for cent_min, cent_max in CENTRALITY_BINS:
        for efficiency in EFFICIENCIES:
            chosen = [p for p in points if p["cent_min"] == cent_min and p["cent_max"] == cent_max and p["efficiency"] == efficiency]
            fit = constant_fit(np.asarray([p["threshold_gev"] for p in chosen]), np.asarray([p["threshold_error_gev"] for p in chosen]))
            flat.append({"efficiency": efficiency, "cent_min": cent_min, "cent_max": cent_max, "cent_center": 0.5 * (cent_min + cent_max), **fit})
    lines = []
    for efficiency in EFFICIENCIES:
        chosen = [f for f in flat if f["efficiency"] == efficiency]
        lines.append({"efficiency": efficiency, **linear_fit(np.asarray([f["cent_center"] for f in chosen]),
                                                             np.asarray([f["value_gev"] for f in chosen]),
                                                             np.asarray([f["error_gev"] for f in chosen]))})
    return points, flat, lines


def population(paths) -> dict:
    """The rows of every FillIsolationPopulation output, concatenated in order; one subtraction method."""
    columns = ("et", "eta", "centrality", "native_bin", "iso_r03", "iso_r04", "method", "weight", "file")
    parts = []
    for path in paths:
        with uproot.open(path) as f:
            parts.append(f["isolation_population"].arrays(columns, library="np"))
    rows = {c: np.concatenate([p[c] for p in parts]) if parts else np.zeros(0) for c in columns}
    methods = sorted(set(np.asarray(rows["method"]).tolist()))
    if len(methods) > 1:
        raise ValueError(f"more than one reconstructed isolation subtraction method: {methods}")
    return rows


def package(paths) -> dict:
    rows = population(paths)
    out = {"schema": "AuAuIsolationPackageV1", "method": "per pT group and 5% centrality slice: the 70/80/90% cutoff of the "
           "reconstructed-isolation spectrum of truth-signal photons; a constant across pT per slice; a line across centrality",
           "histogram": {"nbins": NBINS, "xmin": XMIN, "xmax": XMAX, "weighted": False},
           "pt_groups": [list(map(list, (g, b))) for g, b in PT_GROUPS], "efficiencies": list(EFFICIENCIES),
           "centrality_convention": "native bin b (1-based); slice k = bins 5k+1..5k+5",
           "photons": int(len(rows["et"])), "subtraction_method": int(rows["method"][0]) if len(rows["method"]) else None,
           "inputs": [{"path": str(p), "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in paths], "cones": {}}
    for cone, column in CONES.items():
        points, flat, lines = derive(rows["et"], rows["native_bin"], rows[column])
        out["cones"][cone] = {"lines": lines, "constant_fits": flat, "cutoff_points": points}
    return out


def freeze(derivation: Path, efficiency: float = 0.80, gap: float = 0.8) -> dict:
    """The hash-bound package of the nominal selection: the R = 0.3 line at ``efficiency`` and the sideband ``gap``."""
    doc = json.loads(Path(derivation).read_text())
    if doc.get("schema") != "AuAuIsolationPackageV1":
        raise ValueError("not an isolation derivation package")
    line = next((l for l in doc["cones"]["R30"]["lines"] if abs(l["efficiency"] - efficiency) < 1e-12), None)
    if line is None:
        raise ValueError(f"no R = 0.3 line at efficiency {efficiency}")
    return {"schema": "CanonicalIsolationFreezeV1", "system": "auau",
            "witness": {"radius": 0.3, "stride_index": 0, "subtraction_method": doc["subtraction_method"], "reconstructed": True},
            "threshold": {"intercept_gev": line["intercept_gev"], "slope_gev_per_percent": line["slope_gev_per_percent"],
                          "variable": "cent_value", "efficiency": efficiency},
            "sideband_gap_gev": gap, "domain": {"centrality_min": 0.0, "centrality_max_exclusive": 80.0},
            "derivation": {"file": Path(derivation).name, "package_sha256": hashlib.sha256(Path(derivation).read_bytes()).hexdigest()}}


def check(frozen: dict, et, centrality, native_bin, iso) -> dict:
    """A frozen line against a population of the derivation's kind (15 <= ET < 35 GeV): per 5% slice (native bin) and
    overall, the fraction of photons the line isolates, the realized efficiency against the freeze's target, with its
    binomial error and pull; and the line refitted on the population at the same efficiency against the frozen one,
    threshold by threshold at the slice centres, with the refit's error."""

    if frozen.get("schema") != "CanonicalIsolationFreezeV1":
        raise ValueError("not an isolation freeze")
    line = frozen["threshold"]
    target = float(line["efficiency"])
    et, c, iso = (np.asarray(a, dtype=float) for a in (et, centrality, iso))
    native = np.asarray(native_bin, dtype=np.int64)
    inside = (et >= 15.0) & (et < 35.0) & np.isfinite(iso) & np.isfinite(c) & (native >= 1)
    isolated = iso < line["intercept_gev"] + line["slope_gev_per_percent"] * c

    def fraction(selected):
        n = int(selected.sum())
        if n == 0:
            return {"photons": 0, "efficiency": None, "error": None, "pull": None}
        p = float(isolated[selected].mean())
        error = math.sqrt(max(p * (1.0 - p), 0.0) / n)
        return {"photons": n, "efficiency": p, "error": error, "pull": (p - target) / error if error > 0 else None}

    slices = [{"cent_min": lo, "cent_max": hi, **fraction(inside & ((native - 1) // 5 == k))}
              for k, (lo, hi) in enumerate(CENTRALITY_BINS)]
    out = {"schema": "AuAuIsolationCheckV1", "target_efficiency": target,
           "frozen_line": {"intercept_gev": line["intercept_gev"], "slope_gev_per_percent": line["slope_gev_per_percent"]},
           "slices": slices, "overall": fraction(inside)}
    try:
        _, _, lines = derive(et[inside], native[inside], iso[inside])
    except ValueError as error:  # too few photons in a slice and pT group to refit
        out["refit"] = {"error": str(error)}
        return out
    refit = next(l for l in lines if abs(l["efficiency"] - target) < 1e-12)
    covariance = np.asarray(refit["covariance"], dtype=float)
    centres = np.asarray([0.5 * (lo + hi) for lo, hi in CENTRALITY_BINS])
    design = np.column_stack([np.ones_like(centres), centres])
    difference = (refit["intercept_gev"] - line["intercept_gev"]) + (refit["slope_gev_per_percent"] - line["slope_gev_per_percent"]) * centres
    error = np.sqrt(np.einsum("ij,jk,ik->i", design, covariance, design))
    out["refit"] = {"line": refit, "threshold_minus_frozen_gev": difference.tolist(), "error_gev": error.tolist(),
                    "slice_centres": centres.tolist()}
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], fromfile_prefix_chars="@")
    sub = parser.add_subparsers(dest="command", required=True)
    k = sub.add_parser("check", help="a frozen line against FillIsolationPopulation outputs of new simulation")
    k.add_argument("freeze", type=Path)
    k.add_argument("populations", type=Path, nargs="+")
    k.add_argument("--out", type=Path, required=True)
    d = sub.add_parser("derive", help="the derivation package from FillIsolationPopulation outputs")
    d.add_argument("populations", type=Path, nargs="+")
    d.add_argument("--out", type=Path, required=True)
    f = sub.add_parser("freeze", help="the hash-bound package of the nominal selection")
    f.add_argument("derivation", type=Path)
    f.add_argument("--efficiency", type=float, default=0.80)
    f.add_argument("--gap", type=float, default=0.8)
    f.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.out.exists():
        raise FileExistsError(args.out)
    if args.command == "freeze":
        args.out.write_text(json.dumps(freeze(args.derivation, args.efficiency, args.gap), indent=1) + "\n")
        return 0
    if args.command == "check":
        frozen = json.loads(args.freeze.read_text())
        rows = population(args.populations)
        if len(rows["method"]) and int(rows["method"][0]) != frozen["witness"]["subtraction_method"]:
            raise ValueError("the population's isolation method differs from the freeze's witness")
        result = check(frozen, rows["et"], rows["centrality"], rows["native_bin"], rows["iso_r03"])
        result["freeze"] = {"path": str(args.freeze), "sha256": hashlib.sha256(args.freeze.read_bytes()).hexdigest()}
        result["inputs"] = [{"path": str(p), "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in args.populations]
        args.out.write_text(json.dumps(result, indent=1) + "\n")
        for s in result["slices"]:
            if s["photons"]:
                print(f"{s['cent_min']:>2}-{s['cent_max']:<2}% {s['photons']:>7} photons: {s['efficiency']:.3f} +- {s['error']:.3f}")
        return 0
    result = package(args.populations)
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for cone, block in result["cones"].items():
        for line in block["lines"]:
            print(f"{cone} {line['efficiency']:.2f}: {line['intercept_gev']:.4f} + {line['slope_gev_per_percent']:.5f} c "
                  f"(chi2/ndf {line['chi2']:.1f}/{line['ndf']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
