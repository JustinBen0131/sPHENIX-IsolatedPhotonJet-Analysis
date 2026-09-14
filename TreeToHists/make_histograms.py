#!/usr/bin/env python3
"""Nominal histograms from scored collaborator trees.

Input:  scored collaborator tree files for one collision system and
        ``config/nominal.yaml``.
Output: a JSON payload plus a ROOT file with

        * event-leading ABCD photon counts per photon-pT bin
          (A = tight & isolated, B = tight & non-isolated,
           C = non-tight & isolated, D = non-tight & non-isolated), and
        * the recoil xJgamma spectrum per photon-pT bin in regions A and C,

        each with Sumw2.  These are exactly the inputs of the purity
        correction in FinalAnalysis/run_corrections.py.

Selection (all from the configuration): event ``terminal_status``, vertex
window, Au+Au centrality window, photon ET in [15, 35) GeV and |eta| < 0.7,
isolation cone R = 0.3.  Working-point thresholds come from the configuration
(``--thresholds-from config``, default) or, for mechanism checks only, from the
per-candidate threshold branches stored in the trees (``--thresholds-from
tree``).  The event-leading photon is chosen per region among the photons
passing the kinematic selection: highest ET, ties to the lowest encounter
ordinal.  Weights are the stored ``event_weight`` applied exactly once.

With ``--truth-signal-only`` (simulation) only candidates linked to a truth
photon satisfying the nominal truth-signal contract are kept; the resulting
ABCD counts are the prompt-photon leakage inputs of the purity correction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "PhotonID"))
from photon_selection import (  # noqa: E402
    REGION_ORDER,
    classify,
    delta_phi_min,
    event_leading_index,
    load_config,
    missing_working_points,
    region_masks,
    system_config,
    truth_signal_mask,
    working_point_thresholds,
)

SCHEMA = "PhotonJetNominalHistogramsV1"
LINK_MATCH = 0
TYPE_PHOTON = 1
PHOTON_BRANCHES = [
    "source_file_index", "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
    "photon_encounter_ordinal", "photon_et", "photon_eta", "vertex_z", "centrality",
    "event_weight", "terminal_status", "bdt_score", "iso_r03",
    "bdt_tight_threshold", "bdt_nontight_low_threshold", "bdt_nontight_high_threshold",
    "iso_r03_threshold", "iso_r03_nonisolated_threshold",
]
PAIR_BRANCHES = [
    "source_file_index", "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
    "jet_pt", "jet_eta", "jet_radius", "delta_phi", "xjgamma", "event_weight",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_inputs(inputs: list[Path], lists: list[Path]) -> list[Path]:
    paths = [Path(p) for p in inputs]
    for listing in lists:
        for line in Path(listing).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(Path(line) if Path(line).is_absolute() else Path(listing).parent / line)
    if not paths:
        raise SystemExit("no input files given (use --input or --input-list)")
    return paths


def truth_signal_candidate_keys(root, config: dict[str, Any]) -> set[tuple[int, int, int]]:
    truth = root["truthPhotons"].arrays(library="np")
    signal = truth_signal_mask(truth, config["truth_signal"])
    signal_truth = {
        (int(s), int(h), int(l))
        for s, h, l, ok in zip(truth["source_file_index"], truth["truth_photon_id_hi"], truth["truth_photon_id_lo"], signal)
        if ok
    }
    links = root["recoTruthLinks"].arrays(library="np")
    keys = set()
    for s, rh, rl, th, tl, rt, tt, lc in zip(
        links["source_file_index"], links["reco_id_hi"], links["reco_id_lo"], links["truth_id_hi"], links["truth_id_lo"],
        links["reco_type"], links["truth_type"], links["link_class"],
    ):
        if int(rt) == TYPE_PHOTON and int(tt) == TYPE_PHOTON and int(lc) == LINK_MATCH and (int(s), int(th), int(tl)) in signal_truth:
            keys.add((int(s), int(rh), int(rl)))
    return keys


def build(
    paths: list[Path],
    *,
    system: str,
    config: dict[str, Any],
    thresholds_from: str,
    truth_signal_only: bool,
    non_tight_definition: str | None = None,
) -> dict[str, Any]:
    block = system_config(config, system)
    event_cfg = block["event"]
    non_tight = non_tight_definition or config["abcd"]["non_tight_definition"]
    pt_edges = np.asarray(config["photon"]["pt_edges_gev"], dtype=float)
    xj_edges = np.asarray(config["recoil"]["xj_edges"], dtype=float)
    n_pt, n_xj = len(pt_edges) - 1, len(xj_edges) - 1
    et_min, et_max = float(config["photon"]["et_min_gev"]), float(config["photon"]["et_max_gev"])
    abs_eta_max = float(config["photon"]["abs_eta_max"])
    vz_max = float(event_cfg["abs_vertex_z_max_cm"])
    terminal = int(event_cfg["terminal_status"])
    cent_window = event_cfg.get("centrality_percent")
    jet_radius = float(config["recoil"]["jet_radius"])
    jet_pt_min = float(config["recoil"]["jet_pt_min_gev"])
    jet_abs_eta_max = float(config["recoil"]["jet_abs_eta_max"])
    dphi_min = delta_phi_min(config)

    if thresholds_from == "config":
        missing = missing_working_points(block)
        if missing:
            raise SystemExit(
                f"working points for {system} are not derived in the configuration: {missing}. "
                "Run PhotonID/derive_working_points.py and complete config/nominal.yaml, or use "
                "--thresholds-from tree for a mechanism check with the stored per-candidate thresholds."
            )

    counts = np.zeros((4, n_pt))
    counts_w2 = np.zeros((4, n_pt))
    events = np.zeros((4, n_pt), dtype=np.int64)
    spectra = np.zeros((2, n_pt, n_xj))          # regions A and C
    spectra_w2 = np.zeros((2, n_pt, n_xj))
    pairs_used = 0
    selected_events: set[tuple[int, int, int]] = set()
    inputs: list[dict[str, Any]] = []

    for path in paths:
        with uproot.open(path) as root:
            photons = root["photons"].arrays(PHOTON_BRANCHES, library="np")
            pairs = root["photonJets"].arrays(PAIR_BRANCHES, library="np")
            signal_keys = truth_signal_candidate_keys(root, config) if truth_signal_only else None
        inputs.append({"path": str(path), "sha256": sha256(path), "photons": int(len(photons["photon_et"]))})

        if thresholds_from == "config":
            thresholds = working_point_thresholds(block, photons["photon_et"], photons["centrality"])
        else:
            thresholds = {
                "tight": photons["bdt_tight_threshold"],
                "nontight_low": photons["bdt_nontight_low_threshold"],
                "nontight_high": photons["bdt_nontight_high_threshold"],
                "isolated_max": photons["iso_r03_threshold"],
                "nonisolated_min": photons["iso_r03_nonisolated_threshold"],
            }
        flags = classify(photons["bdt_score"], photons["iso_r03"], thresholds)
        regions = region_masks(flags, non_tight)

        et, eta = photons["photon_et"], photons["photon_eta"]
        kinematic = (et >= et_min) & (et < et_max) & (np.abs(eta) < abs_eta_max)
        kinematic &= photons["terminal_status"] == terminal
        kinematic &= np.abs(photons["vertex_z"]) < vz_max
        if cent_window is not None:
            kinematic &= (photons["centrality"] >= float(cent_window[0])) & (photons["centrality"] < float(cent_window[1]))
        if signal_keys is not None:
            kinematic &= np.asarray([
                (int(s), int(h), int(l)) in signal_keys
                for s, h, l in zip(photons["source_file_index"], photons["candidate_id_hi"], photons["candidate_id_lo"])
            ], dtype=bool)

        # group candidates by event identity
        by_event: dict[tuple[int, int, int], list[int]] = {}
        for index, key in enumerate(zip(photons["source_file_index"], photons["event_id_hi"], photons["event_id_lo"])):
            by_event.setdefault((int(key[0]), int(key[1]), int(key[2])), []).append(index)

        # pairs indexed by candidate identity
        pairs_by_candidate: dict[tuple[int, int, int], list[int]] = {}
        for index, key in enumerate(zip(pairs["source_file_index"], pairs["candidate_id_hi"], pairs["candidate_id_lo"])):
            pairs_by_candidate.setdefault((int(key[0]), int(key[1]), int(key[2])), []).append(index)

        for event_key, rows in by_event.items():
            rows_array = np.asarray(rows)
            weights = photons["event_weight"][rows_array]
            if not np.all(weights == weights[0]):
                raise ValueError(f"{path}: event {event_key} carries more than one event_weight")
            weight = float(weights[0])
            for region_index, region in enumerate(REGION_ORDER):
                eligible = kinematic[rows_array] & regions[region][rows_array]
                leader_local = event_leading_index(et[rows_array], photons["photon_encounter_ordinal"][rows_array], eligible)
                if leader_local < 0:
                    continue
                leader = int(rows_array[leader_local])
                pt_bin = int(np.searchsorted(pt_edges, et[leader], side="right") - 1)
                if not 0 <= pt_bin < n_pt:
                    continue
                counts[region_index, pt_bin] += weight
                counts_w2[region_index, pt_bin] += weight * weight
                events[region_index, pt_bin] += 1
                selected_events.add(event_key)
                if region not in ("A", "C"):
                    continue
                spectrum_index = 0 if region == "A" else 1
                candidate_key = (int(photons["source_file_index"][leader]), int(photons["candidate_id_hi"][leader]), int(photons["candidate_id_lo"][leader]))
                for pair_index in pairs_by_candidate.get(candidate_key, []):
                    if not math.isclose(float(pairs["jet_radius"][pair_index]), jet_radius, abs_tol=1e-9):
                        continue
                    if not (pairs["jet_pt"][pair_index] > jet_pt_min and abs(pairs["jet_eta"][pair_index]) < jet_abs_eta_max):
                        continue
                    if not pairs["delta_phi"][pair_index] > dphi_min:
                        continue
                    xj = float(pairs["xjgamma"][pair_index])
                    xj_bin = int(np.searchsorted(xj_edges, xj, side="right") - 1)
                    if not 0 <= xj_bin < n_xj:
                        continue
                    pair_weight = float(pairs["event_weight"][pair_index])
                    spectra[spectrum_index, pt_bin, xj_bin] += pair_weight
                    spectra_w2[spectrum_index, pt_bin, xj_bin] += pair_weight * pair_weight
                    pairs_used += 1

    return {
        "schema": SCHEMA,
        "system": system,
        "selection_version": config.get("selection_version"),
        "thresholds_source": thresholds_from,
        "truth_signal_only": truth_signal_only,
        "selection": {
            "photon_et_gev": [et_min, et_max], "photon_abs_eta_max": abs_eta_max,
            "abs_vertex_z_max_cm": vz_max, "terminal_status": terminal,
            "centrality_percent": cent_window, "isolation_cone_radius": float(config["isolation"]["cone_radius"]),
            "non_tight_definition": non_tight, "jet_radius": jet_radius, "jet_pt_min_gev": jet_pt_min,
            "jet_abs_eta_max": jet_abs_eta_max, "delta_phi_min_rad": dphi_min,
            "working_points": block.get("working_points") if thresholds_from == "config" else "per_candidate_tree_branches",
        },
        "pt_edges": pt_edges.tolist(),
        "xj_edges": xj_edges.tolist(),
        "region_order": list(REGION_ORDER),
        "abcd_event_leading_counts_ptgamma": counts.tolist(),
        "abcd_event_leading_sumw2_ptgamma": counts_w2.tolist(),
        "abcd_event_leading_events_ptgamma": events.tolist(),
        "recoil_region_order": ["A", "C"],
        "inclusive_recoil_spectra_ptgamma_xj": spectra.tolist(),
        "inclusive_recoil_sumw2_ptgamma_xj": spectra_w2.tolist(),
        "selected_events": len(selected_events),
        "recoil_pairs": pairs_used,
        "inputs": inputs,
    }


def write_root(payload: dict[str, Any], output: Path) -> None:
    """Write the same histograms as ROOT objects (bin errors are sqrt(Sumw2))."""

    import ROOT

    pt_edges = np.asarray(payload["pt_edges"], dtype=float)
    xj_edges = np.asarray(payload["xj_edges"], dtype=float)
    output.parent.mkdir(parents=True, exist_ok=True)
    file = ROOT.TFile(str(output), "RECREATE")
    counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
    counts_w2 = np.asarray(payload["abcd_event_leading_sumw2_ptgamma"])
    for region_index, region in enumerate(payload["region_order"]):
        hist = ROOT.TH1D(f"h_photon_et_region_{region}", f"Event-leading photons, region {region};E_{{T}}^{{#gamma}} [GeV];Counts", len(pt_edges) - 1, pt_edges)
        hist.Sumw2()
        for pt in range(len(pt_edges) - 1):
            hist.SetBinContent(pt + 1, counts[region_index, pt])
            hist.SetBinError(pt + 1, math.sqrt(counts_w2[region_index, pt]))
        hist.Write()
    spectra = np.asarray(payload["inclusive_recoil_spectra_ptgamma_xj"])
    spectra_w2 = np.asarray(payload["inclusive_recoil_sumw2_ptgamma_xj"])
    for spectrum_index, region in enumerate(payload["recoil_region_order"]):
        hist2 = ROOT.TH2D(f"h_xj_vs_photon_et_region_{region}", f"Recoil x_{{J#gamma}}, region {region};E_{{T}}^{{#gamma}} [GeV];x_{{J#gamma}}", len(pt_edges) - 1, pt_edges, len(xj_edges) - 1, xj_edges)
        hist2.Sumw2()
        for pt in range(len(pt_edges) - 1):
            for xj in range(len(xj_edges) - 1):
                hist2.SetBinContent(pt + 1, xj + 1, spectra[spectrum_index, pt, xj])
                hist2.SetBinError(pt + 1, xj + 1, math.sqrt(spectra_w2[spectrum_index, pt, xj]))
        hist2.Write()
    meta = ROOT.TObjString(json.dumps({k: payload[k] for k in ("schema", "system", "selection_version", "selection", "inputs")}, sort_keys=True))
    meta.Write("histogram_provenance")
    file.Close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", action="append", default=[], type=Path, help="scored collaborator tree file; repeatable")
    parser.add_argument("--input-list", action="append", default=[], type=Path)
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path, help="config/nominal.yaml")
    parser.add_argument("--thresholds-from", choices=("config", "tree"), default="config")
    parser.add_argument("--non-tight-definition", choices=("bounded", "complement"))
    parser.add_argument("--truth-signal-only", action="store_true", help="simulation: keep only nominal truth-signal-linked candidates (leakage inputs)")
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args(argv)

    config = load_config(args.config)
    payload = build(
        read_inputs(args.input, args.input_list),
        system=args.system, config=config, thresholds_from=args.thresholds_from,
        truth_signal_only=args.truth_signal_only, non_tight_definition=args.non_tight_definition,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.output_root is not None:
        write_root(payload, args.output_root)
    counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
    print(f"[make_histograms] {args.output_json}: events={payload['selected_events']} pairs={payload['recoil_pairs']} "
          f"A/B/C/D={counts.sum(axis=1).tolist()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
