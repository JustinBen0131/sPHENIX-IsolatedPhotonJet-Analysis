"""Command-line entry points for tree checks and the two response builds.

    python FinalAnalysis/photonjet/cli.py trees validate --input trees/pp.root --model-input-count 11
    python FinalAnalysis/photonjet/cli.py purity --input trees/pp.root --non-tight-definition bounded --isolation-radius 0.3
    python FinalAnalysis/photonjet/cli.py response build --config config/nominal.yaml --system pp \\
        --sample pp_photonjet --input trees/pp_photonjet_sim.root --output-stem response/pp_pairs
    python FinalAnalysis/photonjet/cli.py photon-response build --config config/nominal.yaml --system pp \\
        --sample pp_photonjet --input trees/pp_photonjet_sim.root --output-stem response/pp_photons

``trees validate`` checks a file against the branch contract in ``contracts/``.
``purity`` is the independent event-leading ABCD counter used as a cross-check
of TreeToHists.  ``response build`` constructs the (photon pT, xJ) pair
response and ``photon-response build`` the per-event photon response, both
from photon+jet simulation trees with the nominal truth-signal contract and,
when ``--sample`` is given, the complete analysis weights of
``config/samples.yaml``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "PhotonID"))
sys.path.insert(0, str(REPO / "TreeToHists"))

import numpy as np
import uproot

from photonjet.analysis.photon_response import build_photon_response
from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.response_builder import ResponseBuildConfig, write_response_artifacts
from photonjet.io.tree_validation import validate as validate_trees
from photonjet.provenance import write_json


def _inputs(values: Iterable[Path], lists: Iterable[Path]) -> list[Path]:
    result = [Path(value) for value in values]
    for list_path in lists:
        result.extend(Path(line.strip().split()[0]) for line in Path(list_path).read_text(encoding="utf-8").splitlines()
                      if line.strip() and not line.lstrip().startswith("#"))
    if not result:
        raise ValueError("at least one --input or --input-list is required")
    return result


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", action="append", type=Path, default=[], help="collaborator tree file; repeatable")
    parser.add_argument("--input-list", action="append", type=Path, default=[], help="text file with one path per line")


def _add_nominal(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", required=True, type=Path, help="config/nominal.yaml")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--samples", type=Path, help="config/samples.yaml (default: next to --config)")
    parser.add_argument("--sample", help="sample name in the manifest; enables the complete analysis weights")
    parser.add_argument("--region", choices=("A", "B", "C", "D"), default="A")
    parser.add_argument("--truth-definition", choices=("nominal", "legacy_fixture"), default="nominal",
                        help="legacy_fixture only for trees that predate the R=0.3 truth-isolation branches")


def _selection_from_config(config: dict[str, Any], region: str) -> RecoilSelection:
    return RecoilSelection(
        photon_et_min=float(config["photon"]["et_min_gev"]), photon_et_max=float(config["photon"]["et_max_gev"]),
        photon_abs_eta_max=float(config["photon"]["abs_eta_max"]), jet_pt_min=float(config["recoil"]["jet_pt_min_gev"]),
        jet_abs_eta_max=float(config["recoil"]["jet_abs_eta_max"]),
        delta_phi_min=float(config["recoil"]["delta_phi_min_over_pi"]) * math.pi,
        jet_radius=float(config["recoil"]["jet_radius"]), region=region,
        non_tight_definition=str(config["abcd"]["non_tight_definition"]),
        isolation_radius=float(config["isolation"]["cone_radius"]),
    )


def _complete_weights(paths: list[Path], sample: str, manifest_path: Path) -> dict[tuple[int, int, int], float]:
    """Complete analysis weight per event key for the given sample (see TreeToHists/sample_weights.py)."""

    from sample_weights import complete_event_weights, leading_truth_jet_pt, load_manifest

    manifest = load_manifest(manifest_path)
    spec = manifest.samples.get(sample)
    if spec is None:
        raise SystemExit(f"sample {sample!r} is not in {manifest_path}")
    weights: dict[tuple[int, int, int], float] = {}
    for path in paths:
        with uproot.open(path) as root:
            events = root["events"].arrays(["source_file_index", "event_id_hi", "event_id_lo", "event_weight", "centrality"], library="np")
            leading = None
            if spec.source_factor == "ownership_stitch":
                truth_jets = root["truthJets"].arrays(["source_file_index", "event_id_hi", "event_id_lo", "truth_jet_radius", "truth_jet_pt"], library="np")
                by_event = leading_truth_jet_pt(truth_jets)
                leading = np.asarray([by_event.get((int(s), int(h), int(l)), math.nan) for s, h, l in
                                      zip(events["source_file_index"], events["event_id_hi"], events["event_id_lo"])], dtype=float)
        result = complete_event_weights(manifest, sample, event_weight=events["event_weight"], centrality=events["centrality"],
                                        leading_truth_jet_pt_gev=leading)
        for s, h, l, w, kept in zip(events["source_file_index"], events["event_id_hi"], events["event_id_lo"], result["weight"], result["kept"]):
            if kept:
                weights[(int(s), int(h), int(l))] = float(w)
    return weights


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="photonjet", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = root.add_subparsers(dest="command", required=True)

    trees = commands.add_parser("trees", help="collaborator tree checks")
    tree_commands = trees.add_subparsers(dest="tree_command", required=True)
    validate = tree_commands.add_parser("validate", help="validate a file against the tree contract")
    _add_inputs(validate)
    validate.add_argument("--model-input-count", type=int, choices=(11, 14), help="11 for p+p, 14 for Au+Au")
    validate.add_argument("--require-scaled-bit30", action="store_true", help="p+p data: require the scaled photon-trigger branch")
    validate.add_argument("--report", type=Path)

    purity = commands.add_parser("purity", help="event-leading ABCD counts (cross-check of TreeToHists)")
    _add_inputs(purity)
    purity.add_argument("--non-tight-definition", required=True, choices=("bounded", "complement"))
    purity.add_argument("--isolation-radius", required=True, type=float, choices=(0.3, 0.4))
    purity.add_argument("--output", type=Path)

    response = commands.add_parser("response", help="pair response from photon+jet simulation trees")
    response_commands = response.add_subparsers(dest="response_command", required=True)
    build = response_commands.add_parser("build", help="build the (photon pT x xJ) response bundle")
    _add_inputs(build)
    _add_nominal(build)
    build.add_argument("--dimension", choices=("1D", "2D"), default="2D")
    build.add_argument("--output-stem", required=True, type=Path, help="writes <stem>.npz, <stem>.json and <stem>.receipt.json")
    build.add_argument("--receipt", type=Path)

    photon = commands.add_parser("photon-response", help="per-event photon response for the N_gamma denominator")
    photon_commands = photon.add_subparsers(dest="photon_command", required=True)
    photon_build = photon_commands.add_parser("build")
    _add_inputs(photon_build)
    _add_nominal(photon_build)
    photon_build.add_argument("--output-stem", required=True, type=Path, help="writes <stem>.npz and <stem>.json")
    return root


def _nominal_context(args: argparse.Namespace, inputs: list[Path]):
    from photon_selection import load_config, truth_signal_mask

    config = load_config(args.config)
    contract = config["truth_signal"]
    truth_signal = (lambda columns: truth_signal_mask(columns, contract)) if args.truth_definition == "nominal" else None  # noqa: E731
    weights = None
    if args.sample:
        manifest_path = args.samples or (args.config.parent / "samples.yaml")
        weights = _complete_weights(inputs, args.sample, manifest_path)
    return config, truth_signal, weights


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    inputs = _inputs(args.input, args.input_list)
    if args.command == "trees":
        lines = validate_trees(inputs, model_input_count=args.model_input_count, require_scaled_bit30=args.require_scaled_bit30)
        text = "VALIDATION=PASS\n" + "\n".join(lines) + "\n"
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(text, encoding="utf-8")
        print(text, end="")
        return 0
    if args.command == "purity":
        payload = purity_counts(inputs, non_tight_definition=args.non_tight_definition, isolation_radius=args.isolation_radius)
        if args.output:
            write_json(args.output, payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command == "response":
        config, truth_signal, weights = _nominal_context(args, inputs)
        build_config = ResponseBuildConfig(system=args.system, dimension=args.dimension,
                                           selection=_selection_from_config(config, args.region),
                                           truth_definition=args.truth_definition, truth_signal=truth_signal, event_weights=weights)
        receipt = write_response_artifacts(inputs, args.output_stem, build_config, receipt_path=args.receipt)
        print(json.dumps({"status": receipt["state"], "outputs": [o["name"] for o in receipt["outputs"]],
                          "weights": "complete analysis weights" if weights is not None else "stored event_weight"}, indent=2))
        return 0
    if args.command == "photon-response":
        if args.truth_definition != "nominal":
            raise SystemExit("photon-response build needs the nominal truth definition")
        config, truth_signal, weights = _nominal_context(args, inputs)
        response = build_photon_response(inputs, system=args.system, selection=_selection_from_config(config, args.region),
                                         truth_signal=truth_signal, event_weights=weights)
        response.provenance["sample"] = args.sample
        npz, meta = response.save(args.output_stem)
        print(json.dumps({"status": "PASS", "outputs": [npz.name, meta.name], "counts": response.provenance["counts"]}, indent=2))
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
