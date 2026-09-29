"""The purity-nonclosure systematic's input: in simulation, the truth purity over the leakage-corrected ABCD purity,
per photon pT, from the QA fills (every simulation sample weighted as the measurement weights it; pseudo-data = all
candidates; truth purity = prompt-labelled Region-A leaders over all Region-A leaders; ABCD purity from FinalAnalysis
corrections.solve_leakage_abcd with the responses' signal leakage, the same as the closure figure). A straight line in
photon pT is fitted to the 2 GeV points (PPG12 Fig. 47 precedent) and averaged over each measurement bin with the
simulated Region-A spectrum as weight; the result is the per-bin `corrections.purity_scale` of the systematic variant.
A scale outside 0.5-1.5 is flagged implausible (plausible: false) and is not used without review.

    python -m FinalAnalysis.plots.purity_nonclosure --hists ROOT/PhotonJetHists_v1 --config FinalAnalysis/config/final_analysis.yaml --out N.json
        [--sample-weights W.json] [--unit-weights]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import yaml

from . import render_qa as Q


def nonclosure(store, system, c, measurement_edges):
    closure = Q.simulation_purity(store, system, c, Q.FINE_EDGES)
    if closure is None:
        return None
    fine, truth, corrected = Q.FINE_EDGES, closure['truth'], closure['corrected']
    ratio = truth / corrected
    x = (fine[:-1] + fine[1:]) / 2
    good = np.isfinite(ratio) & (closure['leaders'] > 0) & (corrected > 0.05) & (truth > 0)   # unstable ABCD points stay out
    if good.sum() < 2:
        return None
    slope, intercept = np.polyfit(x[good], ratio[good], 1, w=np.sqrt(closure['leaders'][good]))
    line = lambda z: intercept + slope * z  # noqa: E731
    # photon-weighted average of the line within each measurement bin (Region-A spectrum of the pseudo-data)
    native, edges = closure['native_leaders'], closure['native_edges']
    centers = (edges[:-1] + edges[1:]) / 2
    scale = []
    for lo, hi in zip(measurement_edges[:-1], measurement_edges[1:]):
        m = (centers >= lo) & (centers < hi)
        w = native[m]
        scale.append(float(np.sum(w * line(centers[m])) / w.sum()) if w.sum() > 0 else float(line(0.5 * (lo + hi))))
    return dict(points=dict(pt_edges=fine.tolist(), truth_purity=truth.tolist(), abcd_purity=corrected.tolist(), ratio=ratio.tolist()),
                fit=dict(form='a + b pT', intercept=float(intercept), slope=float(slope), points_used=int(good.sum())),
                purity_scale=scale, plausible=bool(all(0.5 < v < 1.5 for v in scale)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--hists', type=Path, required=True)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--sample-weights', type=Path)
    ap.add_argument('--unit-weights', action='store_true')
    args = ap.parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    store = Q.Store(args.hists, config, json.loads(args.sample_weights.read_text()) if args.sample_weights else None, args.unit_weights)
    edges = np.asarray(config['measurement']['photon_pt_edges_gev'], dtype=float)
    out = dict(schema='PurityNonclosureV1', hists=str(args.hists), unit_weights=args.unit_weights, measurement_pt_edges=edges.tolist(), entries={})
    for key, block in config['systems'].items():
        out['entries'][key] = nonclosure(store, block.get('system', key), block['class'], edges)
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: (v or {}).get('purity_scale') for k, v in out['entries'].items()}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
