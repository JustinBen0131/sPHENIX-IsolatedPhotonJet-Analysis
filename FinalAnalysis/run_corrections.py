#!/usr/bin/env python3
"""Purity correction, unfolding and final points for one collision system.

Inputs
  --data             TreeToHists payload of the data sample
  --leakage          TreeToHists payload of the photon+jet simulation built
                     with --truth-signal-only (prompt-photon leakage)
  --pair-response    stem of the 2D response from
                     ``photonjet/cli.py response build --dimension 2D``
  --photon-response  stem of the per-event photon response from
                     ``photonjet/cli.py photon-response build``
  --iterations       fixed iteration count; when omitted, ``unfolding.iterations``
                     from the configuration is used, and when that is null the
                     maintained scan (2..12, gate and score) selects it

The chain itself lives in ``photonjet/analysis/chain.py``: leakage-aware ABCD
purity, optional combinatoric subtraction, joint (photon pT, xJ) unfolding,
photon-count unfolding, (1/N_gamma) dN/dxJ with toy statistical errors, and
the closure diagnostics.  Per-system defaults (combinatoric subtraction,
purity strategy) come from ``systems.<system>.corrections`` in the
configuration.

Without the two responses only the purity correction runs (useful for a
quick look at raw and corrected spectra).

Output: one JSON with every intermediate quantity and input hash, and a CSV of
the final points (xJ edges, value, statistical error) summed over the photon
pT window.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "PhotonID"))
from photonjet.analysis.background import PT_EDGES, CorrectionInputs, correct_measured  # noqa: E402
from photonjet.analysis.chain import ITERATION_CANDIDATES, load_bundle, restrict, run_chain, select_iterations  # noqa: E402
from photonjet.analysis.photon_response import PhotonResponse  # noqa: E402
from photon_selection import load_config, system_config  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_payload(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "PhotonJetNominalHistogramsV1":
        raise SystemExit(f"{path}: not a TreeToHists payload")
    if not np.allclose(payload["pt_edges"], PT_EDGES):
        raise SystemExit(f"{path}: photon-pT edges {payload['pt_edges']} differ from the correction grid {PT_EDGES.tolist()}")
    return payload


def purity_only(data: dict[str, Any], leakage: dict[str, Any], *, purity_strategy: str) -> dict[str, Any]:
    """The purity correction alone, without a response (no unfolding, no combinatoric subtraction)."""

    spectra = np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float)
    n_pt, n_xj = spectra.shape[1], spectra.shape[2]
    inputs = CorrectionInputs(
        leakage=np.asarray(leakage["abcd_event_leading_counts_ptgamma"], dtype=float),
        leakage_sumw2=np.asarray(leakage["abcd_event_leading_sumw2_ptgamma"], dtype=float),
        combinatoric_reco=np.zeros((n_pt, n_xj)), combinatoric_reco_sumw2=np.zeros((n_pt, n_xj)),
        combinatoric_normalization_reco=np.zeros(n_pt), combinatoric_normalization_reco_sumw2=np.zeros(n_pt),
    )
    measured = {key: np.asarray(data[key], dtype=float) for key in (
        "abcd_event_leading_counts_ptgamma", "abcd_event_leading_sumw2_ptgamma",
        "inclusive_recoil_spectra_ptgamma_xj", "inclusive_recoil_sumw2_ptgamma_xj")}
    result = correct_measured(measured, inputs, subtract_unmatched_recoil=False, purity_strategy=purity_strategy)
    xj_edges = np.asarray(data["xj_edges"], dtype=float)
    widths = np.diff(xj_edges)
    photons = np.asarray(result["photons"])
    inverse = np.divide(1.0, photons, out=np.full(photons.shape, np.nan), where=photons > 0)
    scale = inverse[:, None] / widths[None, :]
    return {
        "stage": "purity_correction_only",
        "purity": result["solutions"],
        "corrected_ptgamma_xj": np.asarray(result["corrected"]).tolist(),
        "corrected_variance_ptgamma_xj": np.asarray(result["variance"]).tolist(),
        "purity_corrected_photons_ptgamma": photons.tolist(),
        "density_per_ptgamma": (np.asarray(result["corrected"]) * scale).tolist(),
        "density_per_ptgamma_error": (np.sqrt(np.clip(result["variance"], 0, None)) * scale).tolist(),
        "xj_edges": xj_edges.tolist(), "pt_edges": data["pt_edges"],
    }


def write_points_csv(result: dict[str, Any], path: Path) -> None:
    edges = np.asarray(result["xj_edges"], dtype=float)
    values = np.asarray(result.get("density", []), dtype=float)
    errors = np.asarray(result.get("density_error", []), dtype=float)
    shown = np.asarray(result.get("reported_xj_mask", np.ones(len(edges) - 1, dtype=bool)), dtype=bool)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["xj_low", "xj_high", "one_over_ngamma_dn_dxj", "stat_error", "reported"])
        for index in range(len(values)):
            writer.writerow([edges[index], edges[index + 1], values[index], errors[index], int(shown[index])])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--leakage", required=True, type=Path)
    parser.add_argument("--pair-response", type=Path, help="stem of <stem>.npz/<stem>.json from response build")
    parser.add_argument("--photon-response", type=Path, help="stem of <stem>.npz/<stem>.json from photon-response build")
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--toys", type=int, help="final toy count (default: unfolding.final_toys)")
    parser.add_argument("--scan-toys", type=int, help="toys per scanned iteration count (default: unfolding.scan_toys)")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--subtract-combinatoric", choices=("auto", "on", "off"), default="auto")
    parser.add_argument("--purity-strategy", choices=("auto", "per_pt", "integrated_15_35_transfer"), default="auto")
    parser.add_argument("--label", default="")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    config = load_config(args.config)
    block = system_config(config, args.system)
    corrections = block.get("corrections", {}) or {}
    unfolding = config.get("unfolding", {}) or {}
    subtract = {"auto": bool(corrections.get("subtract_combinatoric", False)), "on": True, "off": False}[args.subtract_combinatoric]
    strategy = corrections.get("purity_strategy", "per_pt") if args.purity_strategy == "auto" else args.purity_strategy
    data = load_payload(args.data)
    leakage = load_payload(args.leakage)
    if data["system"] != args.system:
        raise SystemExit(f"{args.data} is a {data['system']} payload")
    if not leakage.get("truth_signal_only"):
        raise SystemExit(f"{args.leakage}: leakage payload must be built with --truth-signal-only")

    output: dict[str, Any] = {
        "schema": "PhotonJetFinalResultV1",
        "label": args.label or args.system,
        "system": args.system,
        "selection_version": data.get("selection_version"),
        "inputs": {"data": {"path": str(args.data), "sha256": sha256(args.data)},
                   "leakage": {"path": str(args.leakage), "sha256": sha256(args.leakage)}},
        "settings": {"subtract_combinatoric": subtract, "purity_strategy": strategy},
    }

    if args.pair_response is None or args.photon_response is None:
        output.update(purity_only(data, leakage, purity_strategy=strategy))
        output["note"] = "purity correction only: give --pair-response and --photon-response for unfolded results"
    else:
        bundle = load_bundle(args.pair_response)
        photon = PhotonResponse.load(args.photon_response)
        window = restrict(bundle, photon)
        output["inputs"]["pair_response"] = {"stem": str(args.pair_response), "sha256": sha256(args.pair_response.with_suffix(".npz"))}
        output["inputs"]["photon_response"] = {"stem": str(args.photon_response), "sha256": sha256(args.photon_response.with_suffix(".npz"))}
        final_toys = args.toys if args.toys is not None else int(unfolding.get("final_toys", 3000))
        iterations = args.iterations if args.iterations is not None else unfolding.get("iterations")
        if iterations is not None:
            result = run_chain(data, leakage, window, iterations=int(iterations), subtract_combinatoric=subtract,
                               purity_strategy=strategy, toys=final_toys, seed=args.seed)
        else:
            candidates = [int(k) for k in unfolding.get("candidates", ITERATION_CANDIDATES)]
            scan_toys = args.scan_toys if args.scan_toys is not None else int(unfolding.get("scan_toys", 300))
            result = select_iterations(data, leakage, window, candidates=candidates, scan_toys=scan_toys,
                                       final_toys=final_toys, seed=args.seed,
                                       subtract_combinatoric=subtract, purity_strategy=strategy)
        output.update(result)
        output["stage"] = "unfolded"
        write_points_csv(output, args.output.with_suffix(".csv"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, default=float) + "\n", encoding="utf-8")
    purities = [round(float(s["purity"]), 4) for s in output["purity"]]
    summary = f"[run_corrections] {args.output}: stage={output['stage']} purity={purities}"
    if "selection" in output:
        summary += f" iterations={output['selection']['selected_iterations']} gate={'pass' if output['selection']['selected_passed_gate'] else 'FAIL'}"
    elif "iterations" in output:
        summary += f" iterations={output['iterations']} refold_chi2_ndf={output['refold_chi2_ndf']:.3f}"
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
