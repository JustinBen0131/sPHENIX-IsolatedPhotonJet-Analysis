"""Transparent command-line dispatcher for directly callable stages."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection, recoil_histogram, write_recoil_skim
from photonjet.analysis.response_builder import (
    ResponseBuildConfig,
    write_response_artifacts,
)
from photonjet.io.tree_validation import validate as validate_trees
from photonjet.plotting.acceptance import accept_rendered_plot
from photonjet.plotting.contract import compile_plot_contract, write_histogram_receipt
from photonjet.plotting.render import render_histogram
from photonjet.plotting.root_bridge import (
    emit_root_annotation_header,
    verify_root_annotation_header,
)
from photonjet.provenance import canonical_json, write_json


def _inputs(values: Iterable[Path], lists: Iterable[Path]) -> list[Path]:
    result = [Path(value) for value in values]
    for list_path in lists:
        result.extend(
            Path(line.strip())
            for line in Path(list_path).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    if not result:
        raise ValueError("at least one --input or --input-list is required")
    return result


def _selection(args: argparse.Namespace) -> RecoilSelection:
    return RecoilSelection(
        photon_et_min=args.photon_et_min,
        photon_et_max=args.photon_et_max,
        photon_abs_eta_max=args.photon_abs_eta_max,
        jet_pt_min=args.jet_pt_min,
        jet_abs_eta_max=args.jet_abs_eta_max,
        delta_phi_min=args.delta_phi_min,
        jet_radius=args.jet_radius,
        region=args.region,
        non_tight_definition=args.non_tight_definition,
    )


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", action="append", type=Path, default=[])
    parser.add_argument("--input-list", action="append", type=Path, default=[])


def _add_recoil_selection(
    parser: argparse.ArgumentParser,
    *,
    event_leading_only: bool = False,
) -> None:
    parser.add_argument("--photon-et-min", type=float, default=15.0)
    parser.add_argument("--photon-et-max", type=float, default=35.0)
    parser.add_argument("--photon-abs-eta-max", type=float, default=0.7)
    parser.add_argument("--jet-pt-min", type=float, default=5.0)
    parser.add_argument("--jet-abs-eta-max", type=float, default=0.7)
    parser.add_argument("--delta-phi-min", type=float, default=7.0 * 3.141592653589793 / 8.0)
    parser.add_argument("--jet-radius", type=float, default=0.4)
    parser.add_argument(
        "--region",
        choices=(("A", "B", "C", "D") if event_leading_only else ("inclusive", "A", "B", "C", "D")),
        default="A",
    )
    parser.add_argument(
        "--non-tight-definition",
        choices=("bounded", "complement"),
        default="bounded",
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="photonjet", description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    trees = commands.add_parser("trees", help="validate PhotonJetTrees_v1")
    tree_commands = trees.add_subparsers(dest="tree_command", required=True)
    validate = tree_commands.add_parser("validate", help="validate the public tree contract")
    _add_inputs(validate)
    validate.add_argument("--model-input-count", type=int, choices=(11, 14))
    validate.add_argument("--require-scaled-bit30", action="store_true")
    validate.add_argument("--report", type=Path)

    histogram = commands.add_parser("histogram", help="make a canonical xJgamma payload")
    _add_inputs(histogram)
    _add_recoil_selection(histogram)
    histogram.add_argument("--output", required=True, type=Path)
    histogram.add_argument("--receipt", type=Path)
    histogram.add_argument(
        "--dataset",
        required=True,
        type=Path,
        help="typed dataset manifest bound into the histogram receipt",
    )
    histogram.add_argument("--edges", default=",".join(str(round(0.1 * i, 10)) for i in range(19)))

    skim = commands.add_parser("skim", help="write a documented event-leading recoil skim")
    _add_inputs(skim)
    _add_recoil_selection(skim)
    skim.add_argument("--output", required=True, type=Path)
    skim.add_argument("--receipt", type=Path)

    purity = commands.add_parser("purity", help="compute event-leading ABCD occupancy")
    _add_inputs(purity)
    purity.add_argument("--non-tight-definition", required=True, choices=("bounded", "complement"))
    purity.add_argument("--isolation-radius", required=True, type=float, choices=(0.3, 0.4))
    purity.add_argument("--output", type=Path)

    response = commands.add_parser(
        "response",
        help="construct a response bundle from validated PhotonJetTrees_v1",
    )
    response_commands = response.add_subparsers(dest="response_command", required=True)
    response_build = response_commands.add_parser(
        "build",
        help="build deterministic 1D or 2D response artifacts",
    )
    _add_inputs(response_build)
    _add_recoil_selection(response_build, event_leading_only=True)
    response_build.add_argument("--system", required=True, choices=("pp", "auau"))
    response_build.add_argument("--dimension", choices=("1D", "2D"), default="2D")
    response_build.add_argument("--output-stem", required=True, type=Path)
    response_build.add_argument("--receipt", type=Path)

    plot = commands.add_parser("plot", help="render provenance-compiled plots")
    plot_commands = plot.add_subparsers(dest="plot_command", required=True)
    render = plot_commands.add_parser("render", help="render a receipt-bound xJgamma histogram")
    render.add_argument("--histogram", required=True, type=Path)
    render.add_argument("--histogram-receipt", required=True, type=Path)
    render.add_argument("--dataset", required=True, type=Path)
    render.add_argument("--output", required=True, type=Path)
    render.add_argument("--receipt", type=Path)
    accept = plot_commands.add_parser(
        "accept",
        help="replay and issue the canonical promotion token for one exact render",
    )
    accept.add_argument("--histogram", required=True, type=Path)
    accept.add_argument("--histogram-receipt", required=True, type=Path)
    accept.add_argument("--dataset", required=True, type=Path)
    accept.add_argument("--image", required=True, type=Path)
    accept.add_argument("--render-receipt", required=True, type=Path)
    accept.add_argument("--output", required=True, type=Path)
    root_header = plot_commands.add_parser(
        "emit-root-header", help="emit generated TLatex annotations for a ROOT macro"
    )
    root_header.add_argument("--histogram", required=True, type=Path)
    root_header.add_argument("--histogram-receipt", required=True, type=Path)
    root_header.add_argument("--dataset", required=True, type=Path)
    root_header.add_argument("--output", required=True, type=Path)
    root_header.add_argument("--receipt", type=Path)
    verify_header = plot_commands.add_parser(
        "verify-root-header", help="fail if a generated ROOT annotation include is edited or stale"
    )
    verify_header.add_argument("--histogram", required=True, type=Path)
    verify_header.add_argument("--histogram-receipt", required=True, type=Path)
    verify_header.add_argument("--dataset", required=True, type=Path)
    verify_header.add_argument("--header", required=True, type=Path)
    verify_header.add_argument("--receipt", type=Path)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    inputs = (
        _inputs(args.input, args.input_list)
        if hasattr(args, "input") and hasattr(args, "input_list")
        else []
    )
    if args.command == "trees" and args.tree_command == "validate":
        lines = validate_trees(
            inputs,
            model_input_count=args.model_input_count,
            require_scaled_bit30=args.require_scaled_bit30,
        )
        text = "VALIDATION=PASS\n" + "\n".join(lines) + "\n"
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(text, encoding="utf-8")
        print(text, end="")
        return 0
    if args.command == "histogram":
        edges = [float(value) for value in args.edges.split(",") if value.strip()]
        selection = _selection(args)
        payload = recoil_histogram(inputs, selection, edges)
        write_json(args.output, payload)
        receipt_path = args.receipt or args.output.with_name(args.output.name + ".receipt.json")
        write_histogram_receipt(
            input_paths=inputs,
            histogram_path=args.output,
            selection=selection,
            dataset_manifest_path=args.dataset,
            receipt_path=receipt_path,
        )
        print(
            canonical_json(
                {
                    "status": "PASS",
                    "output": args.output.name,
                    "receipt": receipt_path.name,
                }
            ),
            end="",
        )
        return 0
    if args.command == "skim":
        receipt = write_recoil_skim(inputs, args.output, _selection(args))
        if args.receipt:
            write_json(args.receipt, receipt)
        print(canonical_json({"status": "PASS", **receipt}), end="")
        return 0
    if args.command == "purity":
        payload = purity_counts(
            inputs,
            non_tight_definition=args.non_tight_definition,
            isolation_radius=args.isolation_radius,
        )
        if args.output:
            write_json(args.output, payload)
        print(canonical_json(payload), end="")
        return 0
    if args.command == "response" and args.response_command == "build":
        config = ResponseBuildConfig(
            system=args.system,
            dimension=args.dimension,
            selection=_selection(args),
        )
        receipt = write_response_artifacts(
            inputs,
            args.output_stem,
            config,
            receipt_path=args.receipt,
        )
        print(
            canonical_json(
                {
                    "status": receipt["state"],
                    "outputs": [output["name"] for output in receipt["outputs"]],
                    "receipt": (
                        args.receipt.name
                        if args.receipt is not None
                        else args.output_stem.with_suffix(".receipt.json").name
                    ),
                }
            ),
            end="",
        )
        return 0
    if args.command == "plot":
        if args.plot_command == "accept":
            receipt = accept_rendered_plot(
                histogram_path=args.histogram,
                histogram_receipt_path=args.histogram_receipt,
                dataset_manifest_path=args.dataset,
                image_path=args.image,
                render_receipt_path=args.render_receipt,
                acceptance_receipt_path=args.output,
            )
            print(
                canonical_json(
                    {
                        "status": receipt["state"],
                        "output": args.output.name,
                        "image": args.image.name,
                        "acceptance_sha256": receipt["acceptance_sha256"],
                    }
                ),
                end="",
            )
            return 0
        contract = compile_plot_contract(
            histogram_path=args.histogram,
            histogram_receipt_path=args.histogram_receipt,
            dataset_manifest_path=args.dataset,
        )
        if args.plot_command == "render":
            receipt_path = args.receipt or args.output.with_name(args.output.name + ".receipt.json")
            render_histogram(
                contract=contract,
                histogram_path=args.histogram,
                output_path=args.output,
                receipt_path=receipt_path,
            )
            print(
                canonical_json(
                    {
                        "status": "PASS",
                        "artifact_state": "CANDIDATE_RAW",
                        "output": args.output.name,
                        "receipt": receipt_path.name,
                        "contract_sha256": contract["contract_sha256"],
                    }
                ),
                end="",
            )
            return 0
        if args.plot_command == "emit-root-header":
            output = emit_root_annotation_header(contract, args.output)
            receipt = verify_root_annotation_header(contract, output)
            receipt_path = args.receipt or args.output.with_name(args.output.name + ".receipt.json")
            write_json(receipt_path, receipt)
            print(
                canonical_json(
                    {
                        "status": "PASS",
                        "output": output.name,
                        "receipt": receipt_path.name,
                        "contract_sha256": contract["contract_sha256"],
                    }
                ),
                end="",
            )
            return 0
        if args.plot_command == "verify-root-header":
            receipt = verify_root_annotation_header(contract, args.header)
            receipt_path = args.receipt or args.header.with_name(args.header.name + ".receipt.json")
            write_json(receipt_path, receipt)
            print(canonical_json({"status": "PASS", **receipt}), end="")
            return 0
        raise AssertionError(args.plot_command)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
