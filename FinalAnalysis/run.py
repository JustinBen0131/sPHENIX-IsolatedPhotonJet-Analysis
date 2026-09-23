#!/usr/bin/env python3
"""Purity correction, unfolding and the final xJgamma points for one system.

Inputs (all TreeToHists packages)
  --data       the data sample family
  --leakage    the photon+jet simulation built with --truth-signal-only
               (prompt-photon leakage into B, C, D; region-A row normalises
               the combinatoric template)
  --response   the photon+jet simulation built without --truth-signal-only
               (photon and pair responses); when omitted only the purity
               correction runs

The chain (ported from the validated collaboration implementation,
``photonjet/analysis/chain.py`` at main a7f15eb):

  1. leakage-aware ABCD purity per photon-pT bin and the purity-corrected
     recoil spectrum;
  2. optional combinatoric subtraction (Au+Au in the maintained analysis);
  3. joint iterative-Bayes unfolding of (photon pT, xJ) with the truth
     spectrum as prior and only detector and boundary fakes as fake cause;
  4. one-dimensional unfolding of the purity-corrected photon count;
  5. density = unfolded pairs summed over pT / (unfolded photons x bin width);
  6. toys: the event-bootstrap replicas stored in the data package
     (``--toys bootstrap``) or Gaussian fluctuations of the sufficient
     statistics (``--toys gaussian``), each rerunning 1-5; the error is the
     central-68 half width and the covariance the winsorised toy covariance;
  7. diagnostics: refolding, photon refolding, simulation closure, zero
     efficiency bins, negative-input fraction, toy success fraction.

With ``unfolding.iterations: null`` in the configuration the iteration count
is scanned over the candidates, gated and scored, and the minimum passing score
kept. An entirely failed scan aborts nominal analysis. Output: one JSON with every intermediate quantity and input hash, and
a CSV of the final points.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent


def _load(name: str, path: Path):
    """Import a sibling module by path so that ``statistics.py`` never shadows the standard library."""

    import importlib.util

    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_corrections = _load("photonjet_corrections", HERE / "corrections.py")
_unfolding = _load("photonjet_unfolding", HERE / "unfolding.py")
_statistics = _load("photonjet_statistics", HERE / "statistics.py")
_contract = _load("photonjet_histogram_contract", HERE.parent / "TreeToHists" / "histogram_contract.py")

CorrectionInputs = _corrections.CorrectionInputs
bootstrap_data_toy = _corrections.bootstrap_data_toy
correct_measured = _corrections.correct_measured
draw_combinatoric_template_toy = _corrections.draw_combinatoric_template_toy
draw_joint_data_toy = _corrections.draw_joint_data_toy
Package = _contract.Package
central_68 = _statistics.central_68
chi2_ndf = _statistics.chi2_ndf
gate_and_score = _statistics.gate_and_score
reported_mask = _statistics.reported_mask
toy_covariance = _statistics.toy_covariance
AnalysisResponse = _unfolding.AnalysisResponse
iterative_bayes = _unfolding.iterative_bayes
restrict = _unfolding.restrict
unfold_pairs = _unfolding.unfold_pairs
unfold_photons = _unfolding.unfold_photons


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def measured_arrays(package: Package) -> dict[str, np.ndarray]:
    a = package.arrays
    return {k: np.asarray(a[k], dtype=float) for k in ("abcd_counts", "abcd_counts_sumw2", "recoil_spectra", "recoil_spectra_sumw2",
                                                         "bootstrap_abcd_counts", "bootstrap_recoil_spectra")}


def correction_inputs(leakage: Package, window: AnalysisResponse | None, n_pt: int, n_xj: int) -> CorrectionInputs:
    counts = np.asarray(leakage.arrays["leakage_abcd_counts"], dtype=float)
    counts_w2 = np.asarray(leakage.arrays["leakage_abcd_counts_sumw2"], dtype=float)
    if counts.shape != (4, n_pt):
        raise ValueError(f"leakage package must carry 4 x {n_pt} event-leading counts")
    if window is None:
        template = np.zeros((n_pt, n_xj)); template_w2 = np.zeros((n_pt, n_xj))
    else:
        template, template_w2 = window.xj_combinatoric, window.xj_combinatoric_sumw2
    return CorrectionInputs(leakage=counts, leakage_sumw2=counts_w2, combinatoric_reco=template, combinatoric_reco_sumw2=template_w2,
                            combinatoric_normalization_reco=counts[0], combinatoric_normalization_reco_sumw2=counts_w2[0])


def densities(unfolded: np.ndarray, photons: np.ndarray, xj_edges: np.ndarray, n_pt: int) -> tuple[np.ndarray, np.ndarray]:
    nx = len(xj_edges) - 1
    widths = np.diff(xj_edges)
    blocks = np.asarray(unfolded).reshape(n_pt, nx)
    total = float(np.sum(photons))
    if not math.isfinite(total) or total <= 0 or np.any(~np.isfinite(photons)) or np.any(photons <= 0):
        raise ValueError("normalization requires positive finite unfolded photon denominators")
    summed = np.divide(blocks.sum(axis=0), total * widths, out=np.zeros(nx), where=(total * widths) > 0)
    per_pt = np.divide(blocks, photons[:, None] * widths[None, :], out=np.zeros((n_pt, nx)), where=(photons[:, None] * widths[None, :]) > 0)
    return summed, per_pt


def run_chain(data: Package, leakage: Package, window: AnalysisResponse, *, iterations: int, subtract_combinatoric: bool,
              purity_strategy: str, toys: int, toy_mode: str, seed: int, reported: tuple[float, float]) -> dict[str, Any]:
    g = data.grids
    if not np.allclose(g.xj_edges, window.xj_edges) or not np.allclose(g.pt_edges, window.pt_edges):
        raise ValueError("data grids differ from the response grids")
    n_pt, nx = g.n_pt, g.n_xj
    measured = measured_arrays(data)
    inputs = correction_inputs(leakage, window, n_pt, nx)

    correction = correct_measured(measured, inputs, subtract_unmatched_recoil=subtract_combinatoric, pt_edges=g.pt_edges, purity_strategy=purity_strategy)
    corrected = np.asarray(correction["corrected"]); variance = np.asarray(correction["variance"]); photons = np.asarray(correction["photons"])
    unfolded, refolded = unfold_pairs(corrected, window, iterations)
    unfolded_photons, refolded_photons = unfold_photons(photons, window, iterations)
    density, density_per_pt = densities(unfolded, unfolded_photons, g.xj_edges, n_pt)

    # --- toys ---------------------------------------------------------------
    rng = np.random.default_rng(seed)
    central_methods = tuple(row["method"] for row in correction["solutions"])
    replicas = int(measured["bootstrap_abcd_counts"].shape[0])
    if toy_mode == "bootstrap":
        if toys > replicas:
            raise ValueError("requested more independent bootstrap draws than the package retains")
        toys = replicas if toys <= 0 else toys
    if toys < 2:
        raise ValueError("at least two uncertainty toys are required")
    toy_density: list[np.ndarray] = []; toy_per_pt: list[np.ndarray] = []; toy_measured: list[np.ndarray] = []
    for k in range(toys):
        toy_data = bootstrap_data_toy(measured, k) if toy_mode == "bootstrap" else draw_joint_data_toy(measured, rng)
        if subtract_combinatoric:
            template, normalisation = draw_combinatoric_template_toy(inputs, rng)
        else:
            template, normalisation = inputs.combinatoric_reco, inputs.combinatoric_normalization_reco
        try:
            toy = correct_measured(toy_data, inputs, subtract_unmatched_recoil=subtract_combinatoric, pt_edges=g.pt_edges,
                                   combinatoric_reco=template, combinatoric_normalization_reco=normalisation, purity_strategy=purity_strategy)
            if tuple(row["method"] for row in toy["solutions"]) != central_methods:
                continue  # a different ABCD branch is another estimator, not a fluctuation
            toy_unfolded, _ = unfold_pairs(toy["corrected"], window, iterations)
            toy_photons, _ = unfold_photons(toy["photons"], window, iterations)
        except ValueError:
            continue
        if toy_photons.sum() <= 0:
            continue
        try:
            summed, per_pt = densities(toy_unfolded, toy_photons, g.xj_edges, n_pt)
        except ValueError:
            continue
        toy_density.append(summed); toy_per_pt.append(per_pt); toy_measured.append(np.asarray(toy["corrected"]).reshape(-1))
    successes = len(toy_density)

    if successes > 1:
        toy_array = np.asarray(toy_density)
        lower, upper, density_error = central_68(toy_array)
        density_covariance = toy_covariance(toy_array, density_error, lower, upper)
        _, _, per_pt_error = central_68(np.asarray(toy_per_pt))
        measured_variance = np.asarray(toy_measured).var(axis=0, ddof=1)
    else:
        density_error = np.full(nx, np.nan); density_covariance = np.full((nx, nx), np.nan)
        per_pt_error = np.full_like(density_per_pt, np.nan)
        measured_variance = np.clip(variance.reshape(-1), 1e-12, None)

    # --- diagnostics --------------------------------------------------------
    measured_vector = corrected.reshape(-1)
    live = measured_variance > 0
    refold = chi2_ndf(refolded[live], measured_vector[live], measured_variance[live]) if live.any() else math.nan
    photon_refold = chi2_ndf(refolded_photons, photons, np.maximum(photons, 1.0))
    closure_reco = window.xj_response.sum(axis=0) + window.unfolding_fakes
    mc_unfolded, mc_refolded = iterative_bayes(closure_reco, window.xj_response, window.xj_misses, iterations,
                                               window.xj_truth.reshape(-1), fakes=window.unfolding_fakes)
    closure_truth = window.xj_truth.reshape(-1)
    supported = closure_truth > 0
    mc_closure = chi2_ndf(mc_unfolded[supported], closure_truth[supported], np.clip(closure_truth + mc_unfolded, 1e-12, None)[supported]) if supported.any() else math.nan
    reco_live = closure_reco > 0
    mc_refold = chi2_ndf(mc_refolded[reco_live], closure_reco[reco_live], np.maximum(closure_reco, 1.0)[reco_live]) if reco_live.any() else math.nan
    xj_totals = window.xj_response.sum(axis=1) + window.xj_misses
    xj_efficiency = np.divide(window.xj_response.sum(axis=1), xj_totals, out=np.zeros_like(xj_totals), where=xj_totals > 0)
    shown = np.tile(reported_mask(g.xj_edges, *reported), n_pt)
    photon_totals = window.photon_response.sum(axis=1) + window.photon_misses
    photon_efficiency = np.divide(window.photon_response.sum(axis=1), photon_totals, out=np.zeros_like(photon_totals), where=photon_totals > 0)
    negative_fraction = float(np.mean(measured_vector[shown] < 0)) if shown.any() else 0.0

    return {
        "iterations": iterations, "subtract_combinatoric": subtract_combinatoric, "purity_strategy": purity_strategy,
        "xj_edges": g.xj_edges.tolist(), "pt_edges": g.pt_edges.tolist(),
        "purity": correction["solutions"],
        "corrected_ptgamma_xj": corrected.tolist(), "corrected_variance_ptgamma_xj": variance.tolist(),
        "purity_corrected_photons_ptgamma": photons.tolist(),
        "combinatoric_subtracted_ptgamma_xj": np.asarray(correction["combinatoric_subtracted"]).tolist(),
        "unfolded_ptgamma_xj": unfolded.reshape(n_pt, nx).tolist(), "refolded_ptgamma_xj": refolded.reshape(n_pt, nx).tolist(),
        "unfolded_photons_ptgamma": unfolded_photons.tolist(),
        "density": density.tolist(), "density_error": density_error.tolist(), "density_stat_covariance": density_covariance.tolist(),
        "density_per_ptgamma": density_per_pt.tolist(), "density_per_ptgamma_error": per_pt_error.tolist(),
        "reported_xj_mask": reported_mask(g.xj_edges, *reported).tolist(),
        "toy_mode": toy_mode, "toy_requested": toys, "toy_successes": successes,
        "refold_chi2_ndf": float(refold), "photon_refold_chi2_ndf": float(photon_refold),
        "mc_closure_chi2_ndf": float(mc_closure), "mc_refold_chi2_ndf": float(mc_refold),
        "negative_input_fraction": negative_fraction,
        "response_observability": {
            "xj_zero_efficiency_supported_bins": int(np.sum(supported & (xj_efficiency == 0))),
            "xj_zero_efficiency_reported_bins": int(np.sum(supported & shown & (xj_efficiency == 0))),
            "photon_zero_efficiency_supported_bins": int(np.sum((window.photon_truth > 0) & (photon_efficiency == 0))),
        },
        "uncertainty_note": "Data ABCD/recoil share stored event-bootstrap draws (Gaussian mode is an explicit approximation). "
                            "Response matrices and leakage are fixed. With combinatoric subtraction, template and denominator "
                            "fluctuate independently as Gaussians, following the predecessor; their mutual/response covariance is not retained.",
    }


def select_iterations(data: Package, leakage: Package, window: AnalysisResponse, *, candidates: Sequence[int], scan_toys: int,
                      final_toys: int, seed: int, reported: tuple[float, float], **chain_options: Any) -> dict[str, Any]:
    table = []
    for k in candidates:
        try:
            result = run_chain(data, leakage, window, iterations=k, toys=scan_toys, seed=seed, reported=reported, **chain_options)
            passed, score, metrics = gate_and_score(result, *reported)
        except ValueError as error:
            passed, score, metrics = False, math.inf, {"error": str(error)}
        table.append({"iterations": k, "passed": passed, "score": score, "metrics": metrics})
    passing = [row for row in table if row["passed"] and math.isfinite(row["score"])]
    if not passing:
        raise ValueError(f"FAILED_ANALYSIS: no iteration passed the gate: {table}")
    best = min(passing, key=lambda row: (row["score"], row["iterations"]))
    final = run_chain(data, leakage, window, iterations=best["iterations"], toys=final_toys, seed=seed + 1, reported=reported, **chain_options)
    passed, _, metrics = gate_and_score(final, *reported)
    if not passed:
        raise ValueError(f"FAILED_ANALYSIS: selected iteration failed the final-toy gate: {metrics}")
    final["selection"] = {"candidates": list(candidates), "scan_toys": scan_toys, "final_toys": final_toys, "table": table,
                          "selected_iterations": best["iterations"], "selected_passed_gate": best["passed"]}
    return final


def purity_only(data: Package, leakage: Package, *, purity_strategy: str) -> dict[str, Any]:
    g = data.grids
    measured = measured_arrays(data)
    inputs = correction_inputs(leakage, None, g.n_pt, g.n_xj)
    result = correct_measured(measured, inputs, subtract_unmatched_recoil=False, pt_edges=g.pt_edges, purity_strategy=purity_strategy)
    widths = np.diff(g.xj_edges)
    photons = np.asarray(result["photons"])
    inverse = np.divide(1.0, photons, out=np.full(photons.shape, np.nan), where=photons > 0)
    scale = inverse[:, None] / widths[None, :]
    return {
        "stage": "purity_correction_only", "purity": result["solutions"],
        "corrected_ptgamma_xj": np.asarray(result["corrected"]).tolist(), "corrected_variance_ptgamma_xj": np.asarray(result["variance"]).tolist(),
        "purity_corrected_photons_ptgamma": photons.tolist(),
        "density_per_ptgamma": (np.asarray(result["corrected"]) * scale).tolist(),
        "density_per_ptgamma_error": (np.sqrt(np.clip(result["variance"], 0, None)) * scale).tolist(),
        "xj_edges": g.xj_edges.tolist(), "pt_edges": g.pt_edges.tolist(),
    }


def write_points_csv(result: Mapping[str, Any], path: Path) -> None:
    edges = np.asarray(result["xj_edges"], dtype=float)
    values = np.asarray(result.get("density", []), dtype=float)
    errors = np.asarray(result.get("density_error", []), dtype=float)
    shown = np.asarray(result.get("reported_xj_mask", np.ones(len(edges) - 1, dtype=bool)), dtype=bool)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["xj_low", "xj_high", "one_over_ngamma_dn_dxj", "stat_error", "reported"])
        for index in range(len(values)):
            writer.writerow([edges[index], edges[index + 1], values[index], errors[index], int(shown[index])])


def require_compatible(package: Package, data: Package, config: Mapping[str, Any]) -> None:
    if package.metadata.get("status") != "analysis_ready":
        raise ValueError("FinalAnalysis refuses diagnostic/incomplete histogram packages")
    digest = hashlib.sha256(json.dumps(config, sort_keys=True, allow_nan=False).encode()).hexdigest()
    if package.metadata.get("measurement_sha256") != digest:
        raise ValueError("package was built with a different measurement configuration")
    for key in ("system", "model", "model_sha256", "feature_definition_sha256"):
        if not package.metadata.get(key) or package.metadata[key] != data.metadata.get(key):
            raise ValueError(f"incompatible package binding: {key}")
    if package.grids.to_dict() != data.grids.to_dict():
        raise ValueError("package grids differ")


def json_finite(value):
    """Diagnostic non-finites are JSON null; never silently serialize NaN tokens."""
    if isinstance(value, Mapping): return {key: json_finite(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)): return [json_finite(v) for v in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(value): return None
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path, help="config/measurement.yaml")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--data", required=True, type=Path, help="stem of the data package")
    parser.add_argument("--leakage", required=True, type=Path, help="stem of the truth-signal-only photon+jet package")
    parser.add_argument("--response", type=Path, help="stem of the photon+jet package with responses")
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--toys", choices=("bootstrap", "gaussian"), default="bootstrap")
    parser.add_argument("--final-toys", type=int, help="gaussian mode: final toy count (default unfolding.final_toys)")
    parser.add_argument("--scan-toys", type=int)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--subtract-combinatoric", choices=("auto", "on", "off"), default="auto")
    parser.add_argument("--purity-strategy", choices=("auto", "per_pt", "integrated_15_35_transfer"), default="auto")
    parser.add_argument("--label", default="")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    if args.output.exists() or args.output.with_suffix(".csv").exists():
        raise FileExistsError(args.output)
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if config.get("contract", {}).get("status") != "approved" or not config.get("contract", {}).get("evidence"):
        raise SystemExit("measurement/uncertainty contract is unresolved; no nominal analysis is permitted")
    block = config["systems"][args.system]
    corrections = block.get("corrections", {}) or {}
    unfolding = config.get("unfolding", {}) or {}
    reported = tuple(float(v) for v in config["recoil"]["reported_xj"])
    subtract = {"auto": bool(corrections.get("subtract_combinatoric", False)), "on": True, "off": False}[args.subtract_combinatoric]
    strategy = corrections.get("purity_strategy", "per_pt") if args.purity_strategy == "auto" else args.purity_strategy

    data = Package.load(args.data)
    leakage = Package.load(args.leakage)
    for name, package in (("data", data), ("leakage", leakage)):
        require_compatible(package, data, config)
        if package.metadata.get("system") != args.system:
            raise SystemExit(f"{name} package is a {package.metadata.get('system')} package")
    if data.metadata.get("family") != "data":
        raise ValueError("--data requires the data sample family")
    if not leakage.metadata.get("truth_signal_only"):
        raise SystemExit(f"{args.leakage}: leakage package must be built with --truth-signal-only")
    if leakage.metadata.get("family") != "photonjet":
        raise SystemExit(f"{args.leakage}: leakage package must come from the photonjet family")

    output: dict[str, Any] = {
        "schema": "PhotonJetFinalResultV2", "label": args.label or args.system, "system": args.system,
        "implementation_sha256": {name: sha256(HERE / name) for name in
                                  ("run.py", "corrections.py", "unfolding.py", "statistics.py")},
        "measurement_version": data.metadata.get("measurement_version"),
        "inputs": {"data": {"stem": str(args.data), "sha256": sha256(args.data.with_suffix(".npz"))},
                   "leakage": {"stem": str(args.leakage), "sha256": sha256(args.leakage.with_suffix(".npz"))}},
        "settings": {"subtract_combinatoric": subtract, "purity_strategy": strategy, "toys": args.toys},
    }

    if args.response is None:
        output.update(purity_only(data, leakage, purity_strategy=strategy))
        output["note"] = "purity correction only: give --response for unfolded results"
    else:
        response = Package.load(args.response)
        if response.metadata.get("truth_signal_only") or response.metadata.get("family") != "photonjet":
            raise SystemExit(f"{args.response}: response package must be the photonjet family without --truth-signal-only")
        require_compatible(response, data, config)
        if subtract:
            leak_inputs = {(r["sha256"], r["sample"]) for r in leakage.metadata["inputs"]}
            response_inputs = {(r["sha256"], r["sample"]) for r in response.metadata["inputs"]}
            if leak_inputs != response_inputs or leakage.metadata.get("sample_manifest_sha256") != response.metadata.get("sample_manifest_sha256"):
                raise ValueError("combinatoric template and denominator need the same normalized simulation population")
            if not np.allclose(response.arrays["combinatoric_photon_denominator"], leakage.arrays["leakage_abcd_counts"][0], rtol=1e-9, atol=1e-9):
                raise ValueError("combinatoric denominator differs from leakage A; resolve the population contract")
        window = restrict(response)
        output["inputs"]["response"] = {"stem": str(args.response), "sha256": sha256(args.response.with_suffix(".npz"))}
        final_toys = args.final_toys if args.final_toys is not None else int(unfolding.get("final_toys", 3000))
        if args.toys == "bootstrap":
            final_toys = 0  # every stored replica
        iterations = args.iterations if args.iterations is not None else unfolding.get("iterations")
        options = dict(subtract_combinatoric=subtract, purity_strategy=strategy, toy_mode=args.toys)
        if iterations is not None:
            result = run_chain(data, leakage, window, iterations=int(iterations), toys=final_toys, seed=args.seed, reported=reported, **options)
        else:
            candidates = [int(k) for k in unfolding.get("candidates", range(2, 13))]
            scan_toys = args.scan_toys if args.scan_toys is not None else int(unfolding.get("scan_toys", 300))
            if args.toys == "bootstrap":
                scan_toys = 0
            result = select_iterations(data, leakage, window, candidates=candidates, scan_toys=scan_toys, final_toys=final_toys,
                                       seed=args.seed, reported=reported, **options)
        passed, _, metrics = gate_and_score(result, *reported)
        if not passed:
            raise ValueError(f"FAILED_ANALYSIS: nominal result failed gate: {metrics}")
        output.update(result)
        output["analysis_gate"] = {"passed": passed, "metrics": metrics}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        output["stage"] = "unfolded"
        write_points_csv(output, args.output.with_suffix(".csv"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(json_finite(output), indent=2, sort_keys=True, default=float, allow_nan=False) + "\n", encoding="utf-8")
    purities = [round(float(s["purity"]), 4) for s in output["purity"]]
    summary = f"[run] {args.output}: stage={output['stage']} purity={purities}"
    if "selection" in output:
        summary += f" iterations={output['selection']['selected_iterations']} gate={'pass' if output['selection']['selected_passed_gate'] else 'FAIL'}"
    elif "iterations" in output:
        summary += f" iterations={output['iterations']} refold_chi2_ndf={output['refold_chi2_ndf']:.3f}"
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
