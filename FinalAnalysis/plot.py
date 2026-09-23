#!/usr/bin/env python3
"""Overlay final (1/N_gamma) dN/dxJgamma spectra, optionally with reference points.

Plotting is separate from the numbers: this reads the JSON written by
``run.py`` and draws it. Nothing here changes a value.

Examples
  python FinalAnalysis/plot.py --result "p+p=results/pp.json" \
      --result "Au+Au 0-20%=results/auau.json" --output results/xjgamma.png
  python FinalAnalysis/plot.py --result "p+p=results/pp.json" \
      --reference "ATLAS p+p=reference/atlas_xjgamma.csv:pp" --output results/xjgamma_vs_atlas.png

``--pt-bin all`` (default) draws the spectrum summed over the photon-pT
window; an index draws one photon-pT bin. Only the reported xJ bins of the
result are shown. Reference CSV columns: series, xj_low, xj_high, value,
stat, syst; stat and syst are added in quadrature. A reference at a
different photon pT is a shape comparison only; say so in its label.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def load_result(path: Path, pt_bin: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
    result = json.loads(Path(path).read_text(encoding="utf-8"))
    if result.get("stage") == "unfolded" and not result.get("analysis_gate", {}).get("passed"):
        raise ValueError("refusing to plot an unfolded result without a passed analysis gate")
    edges = np.asarray(result["xj_edges"], dtype=float)
    pt = result["pt_edges"]
    shown = np.asarray(result.get("reported_xj_mask", np.ones(len(edges) - 1, dtype=bool)), dtype=bool)
    if pt_bin == "all":
        if "density" not in result:
            raise SystemExit(f"{path}: no unfolded spectrum; run run.py with --response or pass --pt-bin")
        values = np.asarray(result["density"], dtype=float)
        errors = np.asarray(result["density_error"], dtype=float)
        label = f"{pt[0]:g} < E_T^gamma < {pt[-1]:g} GeV"
    else:
        index = int(pt_bin)
        values = np.asarray(result["density_per_ptgamma"], dtype=float)[index]
        errors = np.asarray(result["density_per_ptgamma_error"], dtype=float)[index]
        label = f"{pt[index]:g} < E_T^gamma < {pt[index + 1]:g} GeV"
    if np.any(shown & (~np.isfinite(values) | ~np.isfinite(errors) | (errors < 0))):
        raise ValueError("persisted reported values/errors are invalid")
    return edges, values, errors, shown, f"{label}, {result.get('stage', '')}"


def load_reference(spec: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path, _, series = spec.rpartition(":")
    rows = [r for r in csv.DictReader(Path(path).open(encoding="utf-8")) if r["series"] == series]
    if not rows:
        raise SystemExit(f"{path}: no rows for series {series!r}")
    low = np.asarray([float(r["xj_low"]) for r in rows]); high = np.asarray([float(r["xj_high"]) for r in rows])
    values = np.asarray([float(r["value"]) for r in rows])
    errors = np.hypot(np.asarray([float(r["stat"]) for r in rows]), np.asarray([float(r.get("syst", 0.0) or 0.0) for r in rows]))
    return np.append(low, high[-1]), values, errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--result", action="append", default=[], help="label=path/to/result.json; repeatable")
    parser.add_argument("--reference", action="append", default=[], help="label=path/to/points.csv:series; repeatable")
    parser.add_argument("--pt-bin", default="all")
    parser.add_argument("--title", default="sPHENIX work in progress")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    if not args.result and not args.reference:
        raise SystemExit("give at least one --result or --reference")

    figure, axis = plt.subplots(figsize=(6.4, 4.8))
    markers = ("o", "s", "^", "D", "v")
    subtitle = ""
    for index, entry in enumerate(args.result):
        label, _, path = entry.partition("=")
        edges, values, errors, shown, subtitle = load_result(Path(path), args.pt_bin)
        centres = 0.5 * (edges[:-1] + edges[1:]); half = 0.5 * np.diff(edges)
        axis.errorbar(centres[shown], values[shown], xerr=half[shown], yerr=errors[shown],
                      fmt=markers[index % len(markers)], capsize=2, label=label or Path(path).stem)
    for entry in args.reference:
        label, _, spec = entry.partition("=")
        edges, values, errors = load_reference(spec)
        centres = 0.5 * (edges[:-1] + edges[1:]); half = 0.5 * np.diff(edges); shown = np.isfinite(values)
        axis.errorbar(centres[shown], values[shown], xerr=half[shown], yerr=errors[shown], fmt="x", mfc="none",
                      capsize=2, alpha=0.8, label=label or spec)
    axis.set_xlabel(r"$x_{J\gamma} = p_T^{jet}/E_T^{\gamma}$")
    axis.set_ylabel(r"$(1/N_{\gamma})\,dN/dx_{J\gamma}$")
    axis.set_xlim(float(edges[0]), float(edges[-1]))
    axis.set_ylim(bottom=0.0)
    axis.set_title(f"{args.title}  ({subtitle})" if subtitle else args.title, fontsize=10)
    axis.legend(frameon=False, fontsize=8)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(args.output, dpi=150)
    print(f"[plot] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
