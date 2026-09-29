#!/usr/bin/env python3
"""The operating-region scan: every configuration explored, with the health of its whole chain.

    python -m FinalAnalysis.scan --config FinalAnalysis/config/final_analysis.yaml --system auau --variants variants.yaml \
        --data data/*.hists.root --simulation auau_photon12 sim12/*.root --simulation auau_photon20 sim20/*.root \
        --output landscape.json

A variant is a named set of Settings overrides (variants.yaml: {name: {field: value, ...}}); "nominal" is the
configuration itself. Every variant runs the complete measurement (measurement.measure) and records its health, never
only its final curve: prior-fit quality, data and response toy success, statistical precision, purity and its
factorization, the iteration gates, closure and half-closure pulls, the size of the combinatoric subtraction and its
simulation closure, the data/simulation gates, and the negative corrected bins. A variant that fails keeps its error in the record. Nothing is selected here: the choice is
made on these measurement-quality criteria and recorded with the whole landscape.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np
import yaml

from .. import analysis as measurement

HERE = Path(__file__).resolve().parent


def _median(values) -> float | None:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else None


def health(result: Mapping[str, Any]) -> dict[str, Any]:
    """The measurement-quality metrics of one measure() result."""

    entry = result["entry"]
    shown = np.asarray(result["reported_xj_mask"], dtype=bool)
    central, sigma = np.asarray(entry["central"], dtype=float), np.asarray(entry["sigma_total"], dtype=float)
    out: dict[str, Any] = {
        "iterations": result["iterations"],
        "prior_fit": {name: {"chi2_ndf": fit["chi2"] / fit["ndf"] if fit["ndf"] else None, "converged": bool(fit["converged"]),
                             "theta": np.asarray(fit["theta"], dtype=float).tolist()} for name, fit in entry["fits"].items()},
        "data_toys": {"requested": result["data_toys_requested"], "accepted": result["data_toys_accepted"],
                      "success": result["data_toys_accepted"] / max(result["data_toys_requested"], 1), "rejected": result["data_toys_rejected"],
                      "finite": int(entry["data_toys_finite"]), "used_after_trim": int(entry["data_toys_used"])},
        "response_toys": {"finite": int(entry["response_toys_finite"]), "used_after_trim": int(entry["response_toys_used"]),
                          "rejected": int(entry["response_toys_rejected"])},
        "relative_uncertainty_median_reported": _median(np.divide(sigma[shown], np.abs(central[shown]),
                                                                 out=np.full(int(shown.sum()), np.nan), where=central[shown] != 0)),
        "negative_reported_bins": int(np.sum(central[shown] < 0)),
        "purity": [{"method": row["method"], "purity": float(row["purity"]),
                    "factorization_residual": float(row["factorization_residual"])} for row in result["purity"]],
    }
    if result.get("gates"):
        g = result["gates"]
        out["gates"] = {"jets_per_photon_corrected_data_over_simulation": np.asarray(g["jets_per_photon"]["corrected_data_over_simulation"]).tolist(),
                        "mean_xj_data_over_simulation": np.asarray(g["mean_xj"]["data_over_simulation"]).tolist(),
                        "simulation_detector_fakes_over_matched": np.asarray(g["simulation"]["detector_fakes_over_matched"]).tolist()}
    if result.get("iteration_rule"):
        rule = result["iteration_rule"]
        out["iteration_rule"] = {"selected": rule["selected"], "photon_refold_eligible": rule["photon_refold"]["eligible_k"],
                                 "gate_failures": {d: {row["k"]: row["non_k_gate_failures"] for row in rows if row["k"] in (2, 3, 4)}
                                                   for d, rows in rule["rows"].items()}}
    if result.get("closure"):
        out["closure"] = {mode: {"max_abs_pull": r["max_abs_pull_reported"], "chi2_ndf_diagonal": r["chi2_ndf_diagonal_reported"],
                                 "toys": r["toys"], "failed": r["toys_failed"]} for mode, r in result["closure"]["modes"].items()}
    if result.get("unfolding_closure"):
        out["unfolding_closure"] = {key: {"chi2_ndf": r["chi2_ndf_bootstrap_covariance_reported_bins"], "max_abs_pull": r["max_abs_pull_reported"]}
                                    for key, r in result["unfolding_closure"]["results"].items()}
    if result.get("combinatoric"):
        c = result["combinatoric"]
        background = np.asarray(c["purity_corrected_background"], dtype=float)
        corrected = np.asarray(result["corrected_ptgamma_xj"], dtype=float)
        out["combinatoric"] = {"fraction_of_corrected_pairs": float(background.sum() / max(abs(corrected.sum() + background.sum()), 1e-12)),
                               "v2": [f.get("v2") for f in c["fits"]], "shape": c["shape"], "scope": c["scope"],
                               "window_occupancy": None if c.get("window_occupancy") is None else np.asarray(c["window_occupancy"]).tolist()}
        closure = c.get("simulation_closure")
        if closure:
            per_pt = closure.get("per_pt")
            out["combinatoric"]["simulation_closure"] = {"error": closure["error"]} if "error" in closure else {
                "estimate_over_truth": np.divide(per_pt["estimate"], per_pt["truth"], out=np.full(len(per_pt["truth"]), np.nan),
                                                 where=np.asarray(per_pt["truth"]) > 0).tolist(),
                "pull": np.asarray(per_pt["pull"]).tolist(),
                "sideband_share_with_truth_jet": np.asarray(closure["sideband_share_with_truth_jet"]).tolist()}
        halves = c.get("simulation_closure_occupancy_halves")
        if halves:
            out["combinatoric"]["occupancy_half_closure_pull"] = {
                mode: r.get("error") or np.asarray(r["per_pt"]["pull"]).tolist() for mode, r in halves.items()}
    return out


def scan(variants: Mapping[str, Mapping[str, Any]], base: measurement.Settings,
         run: Callable[[measurement.Settings], Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One record per variant: its overrides, settings, health, per-photon result, or the error that stopped it."""

    records = []
    for name, overrides in variants.items():
        record: dict[str, Any] = {"variant": name, "overrides": dict(overrides)}
        try:
            settings = dataclasses.replace(base, **overrides)
            record["settings"] = dataclasses.asdict(settings)
            result = run(settings)
        except ValueError as error:  # an invalid or failing configuration is part of the landscape
            record.update(status="FAILED", error=str(error))
        else:
            record.update(status="OK", health=health(result), per_photon=np.asarray(result["per_photon"]).tolist(),
                          sigma_total=np.asarray(result["sigma_total"]).tolist(), reported_xj_mask=np.asarray(result["reported_xj_mask"]).tolist())
        records.append(record)
        print(f"[scan] {name}: {record['status']}{'' if record['status'] == 'OK' else ' ' + record['error'][:160]}", flush=True)
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     fromfile_prefix_chars="@")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--system", required=True,
                        help="a key of the configuration's systems: pp, auau, or a class entry such as auau_cent20_50")
    parser.add_argument("--variants", required=True, type=Path, help="YAML mapping: variant name -> Settings overrides")
    parser.add_argument("--data", required=True, type=Path, nargs="+")
    parser.add_argument("--simulation", required=True, nargs="+", action="append", metavar=("SAMPLE", "FILE"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError(args.output)
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    measurement.physical_system(config, args.system)
    base = measurement.Settings.from_config(config, args.system)
    variants = {"nominal": {}, **(yaml.safe_load(args.variants.read_text(encoding="utf-8")) or {})}
    files = {group[0]: [Path(p) for p in group[1:]] for group in args.simulation}
    samples, identities = measurement.simulation_samples(files, config, args.system, base)
    labels = measurement.half_split(identities)

    def run(settings):
        dphi = measurement.simulation_dphi(files, config, args.system, settings)
        table = measurement.data_table(args.data, config, args.system, settings,
                                       occupancy=measurement.window_occupancy(dphi, settings))
        data = table.measured()
        toys = (measurement.event_bootstrap_data_toys(table, settings) if settings.data_toy_model == "event_bootstrap"
                else measurement.poisson_data_toys(data, settings))
        return measurement.measure(samples, data, toys, settings, half_labels=labels, table=table, simulation_dphi=dphi,
                                   simulation_dphi_halves=measurement.simulation_dphi_halves(files, config, args.system, settings, labels))

    records = scan(variants, base, run)
    output = {"schema": "PhotonJetOperatingRegionScanV1", "system": args.system, "config": {"path": str(args.config), "sha256": measurement.sha256(args.config)},
              "variants": {"path": str(args.variants), "sha256": measurement.sha256(args.variants)}, "records": records,
              "code_sha256": {name: measurement.sha256(HERE / name) for name in measurement.CODE + ("scan.py",)}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(measurement.jsonable(output), indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
