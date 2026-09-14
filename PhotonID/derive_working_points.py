#!/usr/bin/env python3
"""Derive the tight score threshold from scored signal simulation.

Input:  scored photon+jet simulation trees (output of score_trees.py with a
        model) for one collision system.
Output: a YAML fragment with a binned ``tight_threshold`` giving the configured
        signal efficiency in every photon-pT bin (and centrality bin for Au+Au),
        to be pasted into ``config/nominal.yaml`` under
        ``systems.<system>.working_points``.

Signal candidates are reconstructed photons linked (``recoTruthLinks`` with
``link_class == 0``, photon to photon) to a truth photon that satisfies the
nominal truth-signal contract.  Efficiency is the weighted fraction of those
candidates with ``score > threshold``.  Only the tight threshold is derived
here; the non-tight band and the R=0.3 isolation thresholds are analysis
choices that this script does not invent.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

import numpy as np
import uproot
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from photon_selection import load_config, system_config, truth_signal_mask  # noqa: E402


LINK_MATCH = 0
TYPE_PHOTON = 1


def signal_candidates(path: Path, config: dict[str, Any]) -> dict[str, np.ndarray]:
    """Return score, photon_et, centrality, weight of truth-signal-linked candidates."""

    with uproot.open(path) as root:
        photons = root["photons"].arrays(
            ["source_file_index", "candidate_id_hi", "candidate_id_lo", "photon_et", "photon_eta",
             "centrality", "event_weight", "bdt_score"], library="np")
        truth = root["truthPhotons"].arrays(library="np")
        links = root["recoTruthLinks"].arrays(
            ["source_file_index", "reco_id_hi", "reco_id_lo", "truth_id_hi", "truth_id_lo",
             "reco_type", "truth_type", "link_class"], library="np")
    signal = truth_signal_mask(truth, config["truth_signal"])
    signal_truth = {
        (int(s), int(h), int(l))
        for s, h, l, ok in zip(truth["source_file_index"], truth["truth_photon_id_hi"], truth["truth_photon_id_lo"], signal)
        if ok
    }
    matched = set()
    for s, rh, rl, th, tl, rt, tt, lc in zip(
        links["source_file_index"], links["reco_id_hi"], links["reco_id_lo"], links["truth_id_hi"],
        links["truth_id_lo"], links["reco_type"], links["truth_type"], links["link_class"],
    ):
        if int(rt) == TYPE_PHOTON and int(tt) == TYPE_PHOTON and int(lc) == LINK_MATCH and (int(s), int(th), int(tl)) in signal_truth:
            matched.add((int(s), int(rh), int(rl)))
    keep = np.asarray([
        (int(s), int(h), int(l)) in matched
        for s, h, l in zip(photons["source_file_index"], photons["candidate_id_hi"], photons["candidate_id_lo"])
    ], dtype=bool)
    keep &= np.isfinite(photons["bdt_score"])
    keep &= np.abs(photons["photon_eta"]) < float(config["photon"]["abs_eta_max"])
    return {name: np.asarray(values)[keep] for name, values in photons.items()}


def weighted_threshold(score: np.ndarray, weight: np.ndarray, efficiency: float) -> float:
    """Largest threshold t with weighted fraction(score > t) >= efficiency."""

    order = np.argsort(score)[::-1]
    cumulative = np.cumsum(weight[order]) / weight.sum()
    index = int(np.searchsorted(cumulative, efficiency, side="left"))
    index = min(index, len(score) - 1)
    # Threshold just below the selected score so that ``score > t`` keeps it.
    return float(np.nextafter(score[order][index], -np.inf))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", action="append", required=True, type=Path, help="scored photon+jet simulation file; repeatable")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--target-efficiency", type=float, help="overrides tight_signal_efficiency_target from the configuration")
    parser.add_argument("--cent-edges", default="0,20,50,80", help="Au+Au centrality bin edges in percent")
    parser.add_argument("--output", required=True, type=Path, help="YAML fragment to write")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    block = system_config(config, args.system)
    target = args.target_efficiency
    if target is None:
        target = (block.get("working_points") or {}).get("tight_signal_efficiency_target")
    if target is None:
        raise SystemExit(
            f"no tight_signal_efficiency_target for {args.system}: set it in the configuration or pass --target-efficiency"
        )
    target = float(target)
    if not 0.0 < target < 1.0:
        raise SystemExit("target efficiency must be between 0 and 1")

    columns: dict[str, list[np.ndarray]] = {}
    for path in args.input:
        for name, values in signal_candidates(path, config).items():
            columns.setdefault(name, []).append(values)
    data = {name: np.concatenate(values) for name, values in columns.items()}
    if len(data["bdt_score"]) == 0:
        raise SystemExit("no truth-signal-linked scored candidates found; check the inputs were scored with a model")

    pt_edges = [float(x) for x in config["photon"]["pt_edges_gev"]]
    pt_bin = np.searchsorted(pt_edges, data["photon_et"], side="right") - 1
    inside = (pt_bin >= 0) & (pt_bin < len(pt_edges) - 1)
    result: dict[str, Any] = {"mode": "binned", "pt_edges": pt_edges}
    report: list[str] = []
    if args.system == "auau":
        cent_edges = [float(x) for x in args.cent_edges.split(",")]
        cent_bin = np.searchsorted(cent_edges, data["centrality"], side="right") - 1
        inside &= (cent_bin >= 0) & (cent_bin < len(cent_edges) - 1)
        values = []
        for i in range(len(pt_edges) - 1):
            row = []
            for j in range(len(cent_edges) - 1):
                mask = inside & (pt_bin == i) & (cent_bin == j)
                if mask.sum() < 20:
                    raise SystemExit(f"fewer than 20 signal candidates in pT bin {i} centrality bin {j}; cannot derive a threshold")
                row.append(weighted_threshold(data["bdt_score"][mask], data["event_weight"][mask], target))
                report.append(f"pt[{pt_edges[i]},{pt_edges[i+1]}) cent[{cent_edges[j]},{cent_edges[j+1]}): n={int(mask.sum())} t={row[-1]:.6f}")
            values.append(row)
        result["cent_edges"] = cent_edges
        result["values"] = values
    else:
        values = []
        for i in range(len(pt_edges) - 1):
            mask = inside & (pt_bin == i)
            if mask.sum() < 20:
                raise SystemExit(f"fewer than 20 signal candidates in pT bin {i}; cannot derive a threshold")
            values.append(weighted_threshold(data["bdt_score"][mask], data["event_weight"][mask], target))
            report.append(f"pt[{pt_edges[i]},{pt_edges[i+1]}): n={int(mask.sum())} t={values[-1]:.6f}")
        result["values"] = values

    fragment = {
        "working_points": {
            "tight_signal_efficiency_target": target,
            "tight_threshold": result,
            "derived_from": [str(p) for p in args.input],
            "signal_candidates": int(inside.sum()),
        }
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(fragment, sort_keys=False), encoding="utf-8")
    print("\n".join(report))
    print(f"[derive_working_points] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
