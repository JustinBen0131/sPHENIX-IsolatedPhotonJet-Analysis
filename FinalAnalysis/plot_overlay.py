#!/usr/bin/env python3
"""Overlay final (1/N_gamma) dN/dxJgamma spectra, optionally with published reference points.

Examples
  python FinalAnalysis/plot_overlay.py --result "p+p=results/pp.json" \\
      --result "Au+Au 0-20%=results/auau.json" --output results/xjgamma.png

  python FinalAnalysis/plot_overlay.py --result "p+p=results/pp.json" \\
      --reference "ATLAS p+p, 63-80 GeV=FinalAnalysis/reference/atlas_plb789_167_table1_xjgamma.csv:pp" \\
      --output results/xjgamma_vs_atlas.png

``--pt-bin all`` (default) draws the spectrum summed over the 15-35 GeV photon
window; an index draws one photon-pT bin.  Each series is drawn on its own xJ
bins with statistical error bars; only bins with low edge >= 0.3 and centre
<= 1.8 are shown.  Reference points are read from a CSV with the columns
series, xj_low, xj_high, value, stat, syst and drawn with stat and syst added
in quadrature.  A reference at a different photon pT is a shape comparison
only; the legend label should say so.
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

DISPLAY_XJ_LOW = 0.3
DISPLAY_XJ_HIGH = 1.8


def _shown(edges: np.ndarray, values: np.ndarray) -> np.ndarray:
    centres = 0.5 * (edges[:-1] + edges[1:])
    return (edges[:-1] >= DISPLAY_XJ_LOW) & (centres <= DISPLAY_XJ_HIGH) & np.isfinite(values)


def load_result(path: Path, pt_bin: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    result = json.loads(Path(path).read_text(encoding="utf-8"))
    edges = np.asarray(result["xj_edges"], dtype=float)
    pt = result["pt_edges"]
    if pt_bin == "all":
        if "density" not in result:
            raise SystemExit(f"{path}: no unfolded spectrum; run run_corrections.py with both responses or pass --pt-bin")
        values = np.asarray(result["density"], dtype=float)
        errors = np.asarray(result["density_error"], dtype=float)
        label = f"{pt[0]:g} < E_T^gamma < {pt[-1]:g} GeV"
    else:
        index = int(pt_bin)
        values = np.asarray(result["density_per_ptgamma"], dtype=float)[index]
        errors = np.asarray(result["density_per_ptgamma_error"], dtype=float)[index]
        label = f"{pt[index]:g} < E_T^gamma < {pt[index + 1]:g} GeV"
    stage = result.get("stage", "")
    return edges, values, errors, f"{label}, {stage}"


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
    parser.add_argument("--pt-bin", default="all", help="'all' or a photon-pT bin index on the 15/20/25/35 GeV grid")
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
        edges, values, errors, subtitle = load_result(Path(path), args.pt_bin)
        centres = 0.5 * (edges[:-1] + edges[1:]); half = 0.5 * np.diff(edges); shown = _shown(edges, values)
        axis.errorbar(centres[shown], values[shown], xerr=half[shown], yerr=errors[shown],
                      fmt=markers[index % len(markers)], capsize=2, label=label or Path(path).stem)
    for index, entry in enumerate(args.reference):
        label, _, spec = entry.partition("=")
        edges, values, errors = load_reference(spec)
        centres = 0.5 * (edges[:-1] + edges[1:]); half = 0.5 * np.diff(edges); shown = _shown(edges, values)
        axis.errorbar(centres[shown], values[shown], xerr=half[shown], yerr=errors[shown], fmt="x", mfc="none",
                      capsize=2, alpha=0.8, label=label or spec)
    axis.set_xlabel(r"$x_{J\gamma} = p_T^{jet}/E_T^{\gamma}$")
    axis.set_ylabel(r"$(1/N_{\gamma})\,dN/dx_{J\gamma}$")
    axis.set_xlim(0.0, 2.0)
    axis.set_ylim(bottom=0.0)
    axis.set_title(f"{args.title}  ({subtitle})" if subtitle else args.title, fontsize=10)
    axis.legend(frameon=False, fontsize=8)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(args.output, dpi=150)
    print(f"[plot_overlay] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
