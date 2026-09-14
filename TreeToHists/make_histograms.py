#!/usr/bin/env python3
"""Nominal histograms from scored collaborator trees.

Input:  scored collaborator tree files of one sample, ``config/nominal.yaml``
        and ``config/samples.yaml``.
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
ordinal.

Weights: the complete analysis weight of every event from
``sample_weights.py`` (producer weight, source factor, centrality factor),
applied exactly once per photon count and per recoil pair; Sumw2 uses its
square.  Give the sample with ``--sample`` (all inputs) or as a second column
of ``--input-list`` lines.  Without a sample the stored event weight is used,
which is correct for data.

With ``--truth-signal-only`` (simulation) only reconstructed candidates that
are the best match (smallest match metric) of a truth photon satisfying the
nominal truth-signal contract are kept, one candidate per truth photon; the
resulting ABCD counts are the prompt-photon leakage inputs of the purity
correction and the Region-A row is the normalisation of the combinatoric
template.
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
sys.path.insert(0, str(HERE))
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
from sample_weights import complete_event_weights, leading_truth_jet_pt, load_manifest  # noqa: E402

SCHEMA = "PhotonJetNominalHistogramsV1"
LINK_MATCH = 0
TYPE_PHOTON = 1
EVENT_BRANCHES = ["source_file_index", "event_id_hi", "event_id_lo", "event_weight", "centrality"]
PHOTON_BRANCHES = [
    "source_file_index", "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
    "photon_encounter_ordinal", "photon_et", "photon_eta", "vertex_z", "centrality",
    "terminal_status", "bdt_score", "iso_r03",
    "bdt_tight_threshold", "bdt_nontight_low_threshold", "bdt_nontight_high_threshold",
    "iso_r03_threshold", "iso_r03_nonisolated_threshold",
]
PAIR_BRANCHES = [
    "source_file_index", "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
    "jet_pt", "jet_eta", "jet_radius", "delta_phi", "xjgamma",
]
TRUTH_JET_BRANCHES = ["source_file_index", "event_id_hi", "event_id_lo", "truth_jet_radius", "truth_jet_pt"]

EventKey = tuple[int, int, int]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_inputs(inputs: list[Path], lists: list[Path], sample: str | None) -> list[tuple[Path, str | None]]:
    """Return (path, sample) pairs; list lines are ``path`` or ``path sample``."""

    entries = [(Path(p), sample) for p in inputs]
    for listing in lists:
        for line in Path(listing).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            path = Path(parts[0]) if Path(parts[0]).is_absolute() else Path(listing).parent / parts[0]
            entries.append((path, parts[1] if len(parts) > 1 else sample))
    if not entries:
        raise SystemExit("no input files given (use --input or --input-list)")
    return entries


def truth_signal_best_candidates(root, config: dict[str, Any]) -> set[EventKey]:
    """Candidate identities that are the best match of a nominal truth-signal photon.

    One reconstructed candidate per truth photon: the photon-to-photon match
    link with the smallest match metric.  A candidate matched to two truth
    photons is kept once.
    """

    truth = root["truthPhotons"].arrays(list(root["truthPhotons"].keys()), library="np")
    signal = truth_signal_mask(truth, config["truth_signal"])
    signal_truth = {
        (int(s), int(h), int(l))
        for s, h, l, ok in zip(truth["source_file_index"], truth["truth_photon_id_hi"], truth["truth_photon_id_lo"], signal)
        if ok
    }
    links = root["recoTruthLinks"].arrays(list(root["recoTruthLinks"].keys()), library="np")
    best: dict[EventKey, tuple[float, EventKey]] = {}
    for s, rh, rl, th, tl, rt, tt, lc, metric in zip(
        links["source_file_index"], links["reco_id_hi"], links["reco_id_lo"], links["truth_id_hi"], links["truth_id_lo"],
        links["reco_type"], links["truth_type"], links["link_class"], links["match_metric"],
    ):
        truth_key = (int(s), int(th), int(tl))
        if int(rt) != TYPE_PHOTON or int(tt) != TYPE_PHOTON or int(lc) != LINK_MATCH or truth_key not in signal_truth:
            continue
        candidate = (int(s), int(rh), int(rl))
        score = float(metric) if math.isfinite(float(metric)) else math.inf
        if truth_key not in best or score < best[truth_key][0]:
            best[truth_key] = (score, candidate)
    return {candidate for _, candidate in best.values()}


def event_weights_for_file(root, *, sample: str | None, manifest, path: Path) -> tuple[dict[EventKey, float], dict[str, Any]]:
    """Complete weight per event key; events dropped by the contract are absent."""

    events = root["events"].arrays(EVENT_BRANCHES, library="np")
    keys = [(int(s), int(h), int(l)) for s, h, l in zip(events["source_file_index"], events["event_id_hi"], events["event_id_lo"])]
    if len(set(keys)) != len(keys):
        raise ValueError(f"{path}: event identities are not unique")
    if sample is None:
        weights = np.asarray(events["event_weight"], dtype=float)
        if np.any(~np.isfinite(weights)) or np.any(weights < 0):
            raise ValueError(f"{path}: stored event_weight must be finite and non-negative")
        return dict(zip(keys, weights.tolist())), {"sample": None, "components": ["producer_event_weight"], "dropped": {}}
    if manifest is None:
        raise ValueError("a sample manifest is required when --sample is given")
    spec = manifest.samples.get(sample)
    if spec is None:
        raise ValueError(f"sample {sample!r} is not in {manifest.path}")
    leading = None
    if spec.source_factor == "ownership_stitch":
        truth_jets = root["truthJets"].arrays(TRUTH_JET_BRANCHES, library="np")
        by_event = leading_truth_jet_pt(truth_jets)
        leading = np.asarray([by_event.get(key, math.nan) for key in keys], dtype=float)
    result = complete_event_weights(manifest, sample, event_weight=events["event_weight"], centrality=events["centrality"],
                                    leading_truth_jet_pt_gev=leading)
    weights = {key: float(w) for key, w, kept in zip(keys, result["weight"], result["kept"]) if kept}
    return weights, {"sample": sample, "family": result["family"], "components": list(result["components"]),
                     "dropped": result["dropped"], "campaign": result["campaign"], "manifest": result["manifest"]}


def build(
    entries: list[tuple[Path, str | None]],
    *,
    system: str,
    config: dict[str, Any],
    thresholds_from: str,
    truth_signal_only: bool,
    non_tight_definition: str | None = None,
    manifest=None,
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

    counts = np.zeros((4, n_pt)); counts_w2 = np.zeros((4, n_pt)); events_count = np.zeros((4, n_pt), dtype=np.int64)
    spectra = np.zeros((2, n_pt, n_xj)); spectra_w2 = np.zeros((2, n_pt, n_xj))
    pairs_used = 0
    selected_events: set[EventKey] = set()
    inputs: list[dict[str, Any]] = []
    weight_receipts: list[dict[str, Any]] = []
    weighted_events = 0
    sum_weights = 0.0

    for path, sample in entries:
        with uproot.open(path) as root:
            weights, receipt = event_weights_for_file(root, sample=sample, manifest=manifest, path=path)
            photons = root["photons"].arrays(PHOTON_BRANCHES, library="np")
            pairs = root["photonJets"].arrays(PAIR_BRANCHES, library="np")
            signal_keys = truth_signal_best_candidates(root, config) if truth_signal_only else None
        inputs.append({"path": str(path), "sha256": sha256(path), "sample": sample, "photons": int(len(photons["photon_et"])),
                       "events_weighted": len(weights)})
        weight_receipts.append({"path": str(path), **receipt})
        weighted_events += len(weights)
        sum_weights += float(sum(weights.values()))

        if thresholds_from == "config":
            thresholds = working_point_thresholds(block, photons["photon_et"], photons["centrality"])
        else:
            thresholds = {
                "tight": photons["bdt_tight_threshold"], "nontight_low": photons["bdt_nontight_low_threshold"],
                "nontight_high": photons["bdt_nontight_high_threshold"], "isolated_max": photons["iso_r03_threshold"],
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
            kinematic &= np.asarray([(int(s), int(h), int(l)) in signal_keys for s, h, l in
                                     zip(photons["source_file_index"], photons["candidate_id_hi"], photons["candidate_id_lo"])], dtype=bool)

        by_event: dict[EventKey, list[int]] = {}
        for index, key in enumerate(zip(photons["source_file_index"], photons["event_id_hi"], photons["event_id_lo"])):
            by_event.setdefault((int(key[0]), int(key[1]), int(key[2])), []).append(index)
        pairs_by_candidate: dict[EventKey, list[int]] = {}
        for index, key in enumerate(zip(pairs["source_file_index"], pairs["candidate_id_hi"], pairs["candidate_id_lo"])):
            pairs_by_candidate.setdefault((int(key[0]), int(key[1]), int(key[2])), []).append(index)

        for event_key, rows in by_event.items():
            if event_key not in weights:
                continue                                   # dropped by the weight contract
            weight = weights[event_key]
            rows_array = np.asarray(rows)
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
                events_count[region_index, pt_bin] += 1
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
                    xj_bin = int(np.searchsorted(xj_edges, float(pairs["xjgamma"][pair_index]), side="right") - 1)
                    if not 0 <= xj_bin < n_xj:
                        continue
                    spectra[spectrum_index, pt_bin, xj_bin] += weight
                    spectra_w2[spectrum_index, pt_bin, xj_bin] += weight * weight
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
        "abcd_event_leading_events_ptgamma": events_count.tolist(),
        "recoil_region_order": ["A", "C"],
        "inclusive_recoil_spectra_ptgamma_xj": spectra.tolist(),
        "inclusive_recoil_sumw2_ptgamma_xj": spectra_w2.tolist(),
        "selected_events": len(selected_events),
        "recoil_pairs": pairs_used,
        "weights": {"events_weighted": weighted_events, "sum_of_weights": sum_weights, "per_input": weight_receipts},
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
    meta = ROOT.TObjString(json.dumps({k: payload[k] for k in ("schema", "system", "selection_version", "selection", "weights", "inputs")}, sort_keys=True))
    meta.Write("histogram_provenance")
    file.Close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", action="append", default=[], type=Path, help="scored collaborator tree file; repeatable")
    parser.add_argument("--input-list", action="append", default=[], type=Path, help="lines: path [sample]")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path, help="config/nominal.yaml")
    parser.add_argument("--samples", type=Path, help="config/samples.yaml (default: next to --config)")
    parser.add_argument("--sample", help="sample name from the manifest applied to every --input")
    parser.add_argument("--thresholds-from", choices=("config", "tree"), default="config")
    parser.add_argument("--non-tight-definition", choices=("bounded", "complement"))
    parser.add_argument("--truth-signal-only", action="store_true", help="simulation: nominal truth-signal candidates only (leakage inputs)")
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args(argv)

    config = load_config(args.config)
    manifest_path = args.samples or (args.config.parent / "samples.yaml")
    manifest = load_manifest(manifest_path) if manifest_path.is_file() else None
    entries = read_inputs(args.input, args.input_list, args.sample)
    if any(sample is not None for _, sample in entries) and manifest is None:
        raise SystemExit(f"sample manifest not found: {manifest_path}")
    payload = build(entries, system=args.system, config=config, thresholds_from=args.thresholds_from,
                    truth_signal_only=args.truth_signal_only, non_tight_definition=args.non_tight_definition, manifest=manifest)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.output_root is not None:
        write_root(payload, args.output_root)
    counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
    print(f"[make_histograms] {args.output_json}: events={payload['selected_events']} pairs={payload['recoil_pairs']} "
          f"A/B/C/D={counts.sum(axis=1).tolist()} weighted_events={payload['weights']['events_weighted']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
