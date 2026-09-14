#!/usr/bin/env python3
"""Overlay corrected or unfolded xJgamma spectra from run_corrections.py outputs.

Example
  python FinalAnalysis/plot_overlay.py --result "p+p=out/pp_result.json" \
      --result "Au+Au 0-20%=out/auau_result.json" --pt-bin 0 --stage unfolded \
      --output out/xjgamma_overlay_pt15_20.png

Each series is drawn on its own xJ bins with statistical error bars; the
display window is the reported 0.3 <= xJ <= 1.8 range.  The y axis is
(1/N_gamma) dN/dxJ as written by run_corrections.py.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DISPLAY_XJ_LOW = 0.3
DISPLAY_XJ_HIGH = 1.8
KEYS = {
    "corrected": ("corrected_per_photon_ptgamma_xj", "corrected_per_photon_stat_ptgamma_xj"),
    "unfolded": ("unfolded_per_photon_ptgamma_xj", "unfolded_per_photon_stat_ptgamma_xj"),
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--result", action="append", required=True, help="label=path/to/result.json; repeatable")
    parser.add_argument("--pt-bin", type=int, default=0, help="photon-pT bin index on the 15/20/25/35 GeV grid")
    parser.add_argument("--stage", choices=tuple(KEYS), default="unfolded")
    parser.add_argument("--title", default="sPHENIX work in progress")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    figure, axis = plt.subplots(figsize=(6.4, 4.8))
    markers = ("o", "s", "^", "D")
    pt_label = None
    for index, entry in enumerate(args.result):
        label, _, path = entry.partition("=")
        result = json.loads(Path(path).read_text(encoding="utf-8"))
        source = result["unfolding"] if args.stage == "unfolded" else result
        if args.stage == "unfolded" and "unfolding" not in result:
            raise SystemExit(f"{path}: no unfolding block; run run_corrections.py with --response")
        values = np.asarray(source[KEYS[args.stage][0]], dtype=float)[args.pt_bin]
        errors = np.asarray(source[KEYS[args.stage][1]], dtype=float)[args.pt_bin]
        edges = np.asarray(result["xj_edges"], dtype=float)
        centres = 0.5 * (edges[:-1] + edges[1:])
        half = 0.5 * np.diff(edges)
        shown = (edges[:-1] >= DISPLAY_XJ_LOW) & (centres <= DISPLAY_XJ_HIGH) & np.isfinite(values)
        axis.errorbar(centres[shown], values[shown], xerr=half[shown], yerr=errors[shown],
                      fmt=markers[index % len(markers)], capsize=2, label=label or Path(path).stem)
        pt = result["pt_edges"]
        pt_label = f"{pt[args.pt_bin]:g} < E_T^gamma < {pt[args.pt_bin + 1]:g} GeV"
    axis.set_xlabel(r"$x_{J\gamma} = p_T^{jet}/E_T^{\gamma}$")
    axis.set_ylabel(r"$(1/N_{\gamma})\,dN/dx_{J\gamma}$")
    axis.set_xlim(0.0, 2.0)
    axis.set_title(f"{args.title}  ({args.stage}, {pt_label})", fontsize=10)
    axis.legend(frameon=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(args.output, dpi=150)
    print(f"[plot_overlay] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
