#!/usr/bin/env python3
"""Score collaborator trees with the per-system model and regenerate selection flags.

Input:  collaborator tree files (eight trees, see contracts/) with or without a
        stored photon-identification score.
Output: one file per input with the same tree layout, where every photon row in
        ``photons``, ``photonJets`` and ``eventTree`` carries

        * ``bdt_score`` evaluated with the supplied TMVA model (or the stored
          score when no model is given),
        * the configured working-point thresholds copied into the
          ``bdt_*_threshold`` and ``iso_r03_*threshold`` branches, and
        * ``bdt_is_tight``, ``bdt_is_nontight``, ``bdt_is_not_tight`` and
          ``iso_r03_pass`` recomputed from those values.

        A JSON receipt next to each output records the model hash, the
        selection version, the thresholds used and the candidate counts.

Rows are joined by the stable candidate identity
(``source_file_index``, ``candidate_id_hi``, ``candidate_id_lo``), never by
entry order.  Every pair and every event-array entry must resolve to exactly
one photon row or the file is rejected.

Whole files are held in memory; production-size inputs should be split into
parts first (the tree builder already writes parts).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import awkward as ak
import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "TreeProduction"))
# The stored convenience leader indices are recomputed with the tree builder's
# own rule, so scoring can never leave them inconsistent with the new flags.
from build_photonjet_collaboration_tree import _leader_index  # noqa: E402
from photon_selection import (  # noqa: E402
    FEATURE_COUNT,
    classify,
    load_config,
    missing_working_points,
    system_config,
    working_point_thresholds,
)


THRESHOLD_BRANCHES = {
    "tight": "bdt_tight_threshold",
    "nontight_low": "bdt_nontight_low_threshold",
    "nontight_high": "bdt_nontight_high_threshold",
    "isolated_max": "iso_r03_threshold",
    "nonisolated_min": "iso_r03_nonisolated_threshold",
}
FLAG_BRANCHES = {
    "bdt_is_tight": "bdt_is_tight",
    "bdt_is_nontight": "bdt_is_nontight",
    "bdt_is_not_tight": "bdt_is_not_tight",
    "iso_pass": "iso_r03_pass",
}
UPDATED_PHOTON_BRANCHES = ("bdt_score", *THRESHOLD_BRANCHES.values(), *FLAG_BRANCHES.values())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def make_evaluator(model_path: Path, input_count: int):
    """TMVA RBDT evaluator over the ordered model inputs stored in the trees."""

    import ROOT  # PyROOT is required only when a model is supplied

    model = ROOT.TMVA.Experimental.RBDT("myBDT", str(model_path))

    def evaluate(rows: np.ndarray) -> np.ndarray:
        if rows.ndim != 2 or rows.shape[1] != input_count:
            raise ValueError(f"model expects {input_count} inputs, got shape {rows.shape}")
        scores = np.empty(len(rows), dtype=float)
        for index, row in enumerate(rows):
            values = ROOT.std.vector("float")()
            for value in row:
                values.push_back(float(value))
            result = model.Compute(values)
            if len(result) != 1:
                raise ValueError(f"model returned {len(result)} scores instead of one")
            scores[index] = float(result[0])
        return scores

    return evaluate


def candidate_keys(source_index: np.ndarray, hi: np.ndarray, lo: np.ndarray) -> list[tuple[int, int, int]]:
    return list(zip(
        np.asarray(source_index).astype(np.int64).tolist(),
        np.asarray(hi).astype(np.uint64).tolist(),
        np.asarray(lo).astype(np.uint64).tolist(),
    ))


def score_file(
    input_path: Path,
    output_path: Path,
    *,
    system: str,
    config: dict[str, Any],
    evaluator,
    model_path: Path | None,
) -> dict[str, Any]:
    block = system_config(config, system)
    n_inputs = FEATURE_COUNT[system]
    with uproot.open(input_path) as root:
        trees = {name: root[name].arrays(library="ak") for name in root.keys(cycle=False)}
    for required in ("photons", "photonJets", "eventTree"):
        if required not in trees:
            raise ValueError(f"{input_path}: missing tree {required}")

    photons = {name: ak.to_numpy(trees["photons"][name]) for name in ak.fields(trees["photons"])}
    keys = candidate_keys(photons["source_file_index"], photons["candidate_id_hi"], photons["candidate_id_lo"])
    if len(set(keys)) != len(keys):
        raise ValueError(f"{input_path}: photon candidate identities are not unique")
    row_of = {key: index for index, key in enumerate(keys)}

    # --- score -----------------------------------------------------------
    score = np.asarray(photons["bdt_score"], dtype=float).copy()
    state = photons["bdt_evaluation_state"].copy() if "bdt_evaluation_state" in photons else None
    if evaluator is not None:
        inputs = np.stack([np.asarray(photons[f"bdt_input_{j:02d}"], dtype=float) for j in range(n_inputs)], axis=1)
        complete = np.asarray(photons["bdt_input_count"]) == n_inputs
        complete &= np.all(np.isfinite(inputs), axis=1)
        score[:] = np.nan
        if complete.any():
            score[complete] = evaluator(inputs[complete])
        if state is not None:
            state = np.where(complete, 1, 0).astype(state.dtype)
    elif state is not None and np.any(state == 0):
        print(f"[score_trees] {input_path.name}: {int((state == 0).sum())} candidates have no stored score",
              file=sys.stderr)

    # --- thresholds and flags --------------------------------------------
    thresholds = working_point_thresholds(block, photons["photon_et"], photons["centrality"])
    flags = classify(score, photons["iso_r03"], thresholds)

    updated = dict(photons)
    updated["bdt_score"] = score
    if state is not None:
        updated["bdt_evaluation_state"] = state
    for key, branch in THRESHOLD_BRANCHES.items():
        updated[branch] = np.asarray(thresholds[key], dtype=float)
    for key, branch in FLAG_BRANCHES.items():
        updated[branch] = flags[key].astype(photons[branch].dtype)

    # --- propagate to photonJets by candidate identity --------------------
    pairs = {name: ak.to_numpy(trees["photonJets"][name]) for name in ak.fields(trees["photonJets"])}
    pair_rows = np.empty(len(pairs["candidate_id_hi"]), dtype=np.int64)
    for index, key in enumerate(candidate_keys(pairs["source_file_index"], pairs["candidate_id_hi"], pairs["candidate_id_lo"])):
        if key not in row_of:
            raise ValueError(f"{input_path}: photonJets row {index} references an unknown photon candidate")
        pair_rows[index] = row_of[key]
    for branch in UPDATED_PHOTON_BRANCHES:
        if branch in pairs:
            pairs[branch] = updated[branch][pair_rows].astype(pairs[branch].dtype)
    if state is not None and "bdt_evaluation_state" in pairs:
        pairs["bdt_evaluation_state"] = updated["bdt_evaluation_state"][pair_rows]

    # --- propagate to eventTree photon arrays -----------------------------
    event = trees["eventTree"]
    event_fields = {name: event[name] for name in ak.fields(event)}
    source = ak.to_numpy(event["source_file_index"])
    hi_arrays = event["photon_candidate_id_hi"]
    lo_arrays = event["photon_candidate_id_lo"]
    counts = np.asarray(ak.to_numpy(ak.num(hi_arrays)), dtype=np.int64)
    flat_rows = np.empty(int(counts.sum()), dtype=np.int64)
    cursor = 0
    for index in range(len(counts)):
        for hi, lo in zip(ak.to_numpy(hi_arrays[index]), ak.to_numpy(lo_arrays[index])):
            key = (int(source[index]), int(hi), int(lo))
            if key not in row_of:
                raise ValueError(f"{input_path}: eventTree row {index} references an unknown photon candidate")
            flat_rows[cursor] = row_of[key]
            cursor += 1
    for branch in UPDATED_PHOTON_BRANCHES:
        jagged = f"photon_{branch}"
        if jagged in event_fields:
            flat = ak.to_numpy(ak.flatten(event_fields[jagged]))
            values = updated[branch][flat_rows].astype(flat.dtype)
            event_fields[jagged] = ak.unflatten(values, counts)

    # Refresh the stored leader indices: they depend on the score, the
    # thresholds and the R=0.4 isolation witnesses, so new flags imply new
    # leaders.  The R=0.4 witnesses are untouched by this stage.
    leader_fields = [name for name in event_fields if name.startswith("leader_")]
    if leader_fields:
        leader_inputs = ("bdt_score", "bdt_tight_threshold", "bdt_nontight_low_threshold",
                         "bdt_nontight_high_threshold", "iso_r04", "iso_r04_threshold",
                         "iso_r04_nonisolated_threshold", "photon_et", "photon_encounter_ordinal")
        rows_by_event: list[list[dict[str, Any]]] = []
        cursor = 0
        for count in counts:
            rows = [
                {name: updated[name][flat_rows[cursor + offset]] for name in leader_inputs}
                for offset in range(int(count))
            ]
            rows_by_event.append(rows)
            cursor += int(count)
        for name in leader_fields:
            region = name.split("_")[1]
            complement = name.endswith("_complement_index")
            dtype = ak.to_numpy(event_fields[name]).dtype
            event_fields[name] = np.asarray(
                [_leader_index(rows, region, complement=complement) for rows in rows_by_event],
                dtype=dtype,
            )

    # --- write ------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as root:
        for name, arrays in trees.items():
            if name == "photons":
                write_tree(root, name, updated)
            elif name == "photonJets":
                write_tree(root, name, pairs)
            elif name == "eventTree":
                write_tree(root, name, event_fields)
            else:
                write_tree(root, name, {field: arrays[field] for field in ak.fields(arrays)})

    receipt = {
        "schema": "PhotonSelectionReceiptV1",
        "selection_version": config.get("selection_version"),
        "system": system,
        "input": {"path": str(input_path), "sha256": sha256(input_path)},
        "output": {"path": str(output_path), "sha256": sha256(output_path)},
        "model": None if model_path is None else {"path": str(model_path), "sha256": sha256(model_path), "input_count": n_inputs},
        "score_source": "model" if evaluator is not None else "stored",
        "working_points": block.get("working_points"),
        "working_points_missing": missing_working_points(block),
        "counts": {
            "photons": int(len(score)),
            "scored": int(np.isfinite(score).sum()),
            "id_valid": int(flags["id_valid"].sum()),
            "iso_valid": int(flags["iso_valid"].sum()),
            "tight": int(flags["bdt_is_tight"].sum()),
            "nontight": int(flags["bdt_is_nontight"].sum()),
            "isolated_r03": int(flags["iso_pass"].sum()),
            "nonisolated_r03": int(flags["iso_nonisolated"].sum()),
            "pairs": int(len(pair_rows)),
            "events": int(len(counts)),
        },
    }
    receipt_path = output_path.with_suffix(".selection.json")
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def write_tree(root, name: str, columns: dict[str, Any]) -> None:
    """Write one tree preserving the input branch order.

    ``root[name] = dict`` lets uproot choose the layout; the public tree
    contract fixes the branch order, so the branches are declared explicitly.
    """

    types = {}
    for field, values in columns.items():
        array = values if isinstance(values, ak.Array) else ak.Array(np.asarray(values))
        types[field] = str(ak.type(array)).split(" * ", 1)[1]
    tree = root.mktree(name, types)
    if len(next(iter(columns.values()))) > 0:
        tree.extend(columns)


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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", action="append", default=[], type=Path, help="collaborator tree file; repeatable")
    parser.add_argument("--input-list", action="append", default=[], type=Path, help="text file with one input path per line")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path, help="config/nominal.yaml")
    parser.add_argument("--model", type=Path, help="TMVA RBDT ROOT file; when omitted the stored score is kept")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    block = system_config(config, args.system)
    model_path = args.model
    if model_path is None and block.get("model_file"):
        model_path = (args.config.parent / block["model_file"]).resolve()
    evaluator = None
    if model_path is not None:
        if not model_path.is_file():
            raise SystemExit(f"model file not found: {model_path}")
        pinned = block.get("model_sha256")
        if pinned and sha256(model_path) != pinned:
            raise SystemExit(f"model {model_path} does not match model_sha256 in the configuration")
        evaluator = make_evaluator(model_path, FEATURE_COUNT[args.system])
    missing = missing_working_points(block)
    if missing:
        print(f"[score_trees] working points not derived for {args.system}: {missing}; flags are written as 0 "
              "and thresholds as NaN until config/nominal.yaml is completed", file=sys.stderr)

    for input_path in read_inputs(args.input, args.input_list):
        output_path = args.output_dir / input_path.name
        if output_path.resolve() == input_path.resolve():
            raise SystemExit("output directory must differ from the input location")
        receipt = score_file(input_path, output_path, system=args.system, config=config, evaluator=evaluator, model_path=model_path)
        counts = receipt["counts"]
        print(f"[score_trees] {output_path}: photons={counts['photons']} scored={counts['scored']} "
              f"tight={counts['tight']} nontight={counts['nontight']} isolated_r03={counts['isolated_r03']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
