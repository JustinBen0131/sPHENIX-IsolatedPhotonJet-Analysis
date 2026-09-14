#!/usr/bin/env python3
"""Purity correction and unfolding of the nominal recoil xJgamma spectra.

Input
  --data       TreeToHists payload for the data sample
  --leakage    TreeToHists payload for the photon+jet simulation built with
               ``--truth-signal-only`` (prompt-photon ABCD leakage)
  --response   stem of the response bundle written by
               ``python FinalAnalysis/photonjet/cli.py response build`` on the
               photon+jet simulation trees (``<stem>.npz`` and ``<stem>.json``)
  --iterations iterative-Bayes iteration count (required unless
               ``unfolding.iterations`` is set in the configuration)

Steps
  1. Leakage-aware ABCD: per photon-pT bin solve for the signal count in
     region A and the background transfer A/C from the data counts and the
     simulated prompt-photon leakage (photonjet.analysis.background).
  2. Purity-corrected recoil spectrum: region A minus transfer x region C,
     with the region-C prompt leakage restored.  Statistical variance from
     the Sumw2 of both regions.  No combinatoric (unmatched-recoil)
     subtraction is applied in this first pass; the kernel supports it once a
     matched/unmatched template is supplied.
  3. Unfolding: the corrected spectrum is placed on the response's reconstructed
     (pT, xJ) grid and unfolded with the D'Agostini iterative Bayes kernel
     including misses and fakes (photonjet.analysis.unfolding).  Statistical
     uncertainties are propagated by Gaussian toys of the corrected spectrum
     only; response statistics are not included.
  4. Normalisation: each pT bin is divided by its purity-corrected photon
     count from step 1 and by the xJ bin width, giving (1/N_gamma) dN/dxJ.
     The photon count is not corrected for reconstruction efficiency here.

Output: one JSON file with the corrected and unfolded spectra, uncertainties,
the ABCD solutions, the refolding chi2/ndf and the input hashes.
"""

from __future__ import annotations

import argparse
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
from photonjet.analysis.unfolding import chi2_ndf, iterative_bayes  # noqa: E402
from photon_selection import load_config  # noqa: E402


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


def correct(data: dict[str, Any], leakage: dict[str, Any]) -> dict[str, Any]:
    counts = np.asarray(data["abcd_event_leading_counts_ptgamma"], dtype=float)
    counts_w2 = np.asarray(data["abcd_event_leading_sumw2_ptgamma"], dtype=float)
    spectra = np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float)
    spectra_w2 = np.asarray(data["inclusive_recoil_sumw2_ptgamma_xj"], dtype=float)
    leak = np.asarray(leakage["abcd_event_leading_counts_ptgamma"], dtype=float)
    leak_w2 = np.asarray(leakage["abcd_event_leading_sumw2_ptgamma"], dtype=float)
    if not np.allclose(data["xj_edges"], leakage["xj_edges"]):
        raise SystemExit("data and leakage payloads use different xJ edges")
    n_pt, n_xj = spectra.shape[1], spectra.shape[2]
    inputs = CorrectionInputs(
        leakage=leak, leakage_sumw2=leak_w2,
        combinatoric_reco=np.zeros((n_pt, n_xj)), combinatoric_reco_sumw2=np.zeros((n_pt, n_xj)),
        combinatoric_normalization_reco=np.zeros(n_pt), combinatoric_normalization_reco_sumw2=np.zeros(n_pt),
    )
    measured = {
        "abcd_event_leading_counts_ptgamma": counts,
        "abcd_event_leading_sumw2_ptgamma": counts_w2,
        "inclusive_recoil_spectra_ptgamma_xj": spectra,
        "inclusive_recoil_sumw2_ptgamma_xj": spectra_w2,
    }
    return correct_measured(measured, inputs, subtract_unmatched_recoil=False)


def load_response(stem: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    npz = np.load(stem.with_suffix(".npz"))
    meta = json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
    if meta.get("schema") != "PhotonJetResponseBundleV1" or meta.get("dimension") != "2D":
        raise SystemExit(f"{stem}: a 2D PhotonJetResponseBundleV1 is required")
    return {name: npz[name] for name in npz.files}, meta


def unfold(
    corrected: np.ndarray,
    variance: np.ndarray,
    xj_edges: np.ndarray,
    response: dict[str, np.ndarray],
    *,
    iterations: int,
    toys: int,
    seed: int,
) -> dict[str, Any]:
    reco_edges = np.asarray(response["reco_ptgamma_edges"], dtype=float)
    truth_edges = np.asarray(response["truth_ptgamma_edges"], dtype=float)
    if not np.allclose(response["xj_edges"], xj_edges):
        raise SystemExit("response xJ edges differ from the histogram xJ edges")
    n_xj = len(xj_edges) - 1
    reco_bin = [int(np.flatnonzero(np.isclose(reco_edges, edge))[0]) for edge in PT_EDGES[:-1]]
    truth_bin = [int(np.flatnonzero(np.isclose(truth_edges, edge))[0]) for edge in PT_EDGES[:-1]]
    if any(np.isclose(reco_edges[r + 1], hi) is False for r, hi in zip(reco_bin, PT_EDGES[1:])):
        raise SystemExit("response reconstructed pT edges do not contain the measured 15/20/25/35 GeV bins")

    def place(values: np.ndarray) -> np.ndarray:
        vector = np.zeros(len(reco_edges) - 1) if False else np.zeros((len(reco_edges) - 1) * n_xj)
        for local, r in enumerate(reco_bin):
            vector[r * n_xj:(r + 1) * n_xj] = values[local]
        return vector

    def extract(vector: np.ndarray) -> np.ndarray:
        return np.stack([vector[t * n_xj:(t + 1) * n_xj] for t in truth_bin])

    measured = place(corrected)
    matrix = response["matrix"]
    misses = response["misses"]
    fakes = response["unfolding_fakes"]
    unfolded, refolded = iterative_bayes(measured, matrix, misses, iterations, fakes=fakes)
    rng = np.random.default_rng(seed)
    sigma = np.sqrt(np.clip(variance, 0.0, None))
    draws = []
    for _ in range(toys):
        toy = corrected + rng.normal(size=corrected.shape) * sigma
        draws.append(iterative_bayes(place(toy), matrix, misses, iterations, fakes=fakes)[0])
    toy_std = np.std(np.stack(draws), axis=0, ddof=1) if toys > 1 else np.full_like(unfolded, np.nan)
    measured_mask = measured != 0
    refold = chi2_ndf(refolded[measured_mask], measured[measured_mask], place(variance)[measured_mask]) if measured_mask.any() else float("nan")
    return {
        "iterations": iterations,
        "toys": toys,
        "seed": seed,
        "unfolded_ptgamma_xj": extract(unfolded).tolist(),
        "unfolded_stat_ptgamma_xj": extract(toy_std).tolist(),
        "refold_chi2_ndf": float(refold),
        "reco_grid_bins": reco_bin,
        "truth_grid_bins": truth_bin,
    }


def normalise(values: np.ndarray, errors: np.ndarray, photons: np.ndarray, xj_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    width = np.diff(xj_edges)
    inverse = np.divide(1.0, photons, out=np.full(photons.shape, np.nan), where=photons > 0)
    scale = inverse[:, None] / width[None, :]
    return values * scale, errors * scale


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--leakage", required=True, type=Path)
    parser.add_argument("--response", type=Path, help="response bundle stem; when omitted only the purity correction runs")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--toys", type=int, default=200)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--label", default="")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    config = load_config(args.config)
    data = load_payload(args.data)
    leakage = load_payload(args.leakage)
    if not leakage.get("truth_signal_only"):
        raise SystemExit(f"{args.leakage}: leakage payload must be built with --truth-signal-only")
    xj_edges = np.asarray(data["xj_edges"], dtype=float)

    result = correct(data, leakage)
    corrected = np.asarray(result["corrected"])
    variance = np.asarray(result["variance"])
    photons = np.asarray(result["photons"])
    corrected_norm, corrected_err = normalise(corrected, np.sqrt(np.clip(variance, 0, None)), photons, xj_edges)
    output: dict[str, Any] = {
        "schema": "PhotonJetCorrectedSpectraV1",
        "label": args.label or data["system"],
        "system": data["system"],
        "selection_version": data.get("selection_version"),
        "pt_edges": data["pt_edges"],
        "xj_edges": xj_edges.tolist(),
        "abcd_solutions": result["solutions"],
        "purity_corrected_photons_ptgamma": photons.tolist(),
        "corrected_ptgamma_xj": corrected.tolist(),
        "corrected_variance_ptgamma_xj": variance.tolist(),
        "corrected_per_photon_ptgamma_xj": corrected_norm.tolist(),
        "corrected_per_photon_stat_ptgamma_xj": corrected_err.tolist(),
        "combinatoric_subtraction": "not applied",
        "photon_efficiency_correction": "not applied",
        "inputs": {
            "data": {"path": str(args.data), "sha256": sha256(args.data)},
            "leakage": {"path": str(args.leakage), "sha256": sha256(args.leakage)},
        },
    }

    if args.response is not None:
        iterations = args.iterations if args.iterations is not None else config.get("unfolding", {}).get("iterations")
        if iterations is None:
            raise SystemExit("no iteration count: pass --iterations or set unfolding.iterations in the configuration "
                             "after running the maintained refold/toy scan")
        response, meta = load_response(args.response)
        unfolded = unfold(corrected, variance, xj_edges, response, iterations=int(iterations), toys=args.toys, seed=args.seed)
        values = np.asarray(unfolded["unfolded_ptgamma_xj"])
        errors = np.asarray(unfolded["unfolded_stat_ptgamma_xj"])
        norm, norm_err = normalise(values, errors, photons, xj_edges)
        unfolded["unfolded_per_photon_ptgamma_xj"] = norm.tolist()
        unfolded["unfolded_per_photon_stat_ptgamma_xj"] = norm_err.tolist()
        unfolded["response"] = {
            "stem": str(args.response), "sha256_npz": sha256(args.response.with_suffix(".npz")),
            "status": meta.get("status"), "system": meta.get("system"),
        }
        output["unfolding"] = unfolded

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    purities = [round(s["purity"], 4) for s in result["solutions"]]
    print(f"[run_corrections] {args.output}: purity per pT bin={purities} photons={photons.round(2).tolist()}"
          + (f" refold chi2/ndf={output['unfolding']['refold_chi2_ndf']:.3f}" if "unfolding" in output else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
