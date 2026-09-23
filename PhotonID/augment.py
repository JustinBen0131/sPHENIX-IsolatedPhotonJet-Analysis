#!/usr/bin/env python3
"""Evaluate one registered model on canonical trees and write a score sidecar.

The base tree is never modified. For each input file this writes a sidecar
ROOT file with one table, ``PhotonScores``, keyed by the stable candidate
identity (``source``, ``event``, ``photon`` 128-bit words) so it joins to the
base tree by identity and never by row order.

Per candidate the sidecar records

    model            registry name
    score            the model output, NaN unless state is evaluated_finite
    state            one of
                       1 evaluated_finite     inputs complete, score finite
                       2 evaluated_nonfinite  inputs complete, score not finite
                       3 missing_inputs       at least one input is NaN
                       4 out_of_domain        outside the model's domain
                       5 model_unavailable    no bound model file
    complete, in_domain, shower_valid   the three facts behind the state

Working points, tight/non-tight classes and ABCD regions are not decided
here. They are measurement choices and live in ``config/measurement.yaml``,
applied by ``TreeToHists/build.py``.

Usage
    python PhotonID/augment.py --registry PhotonID/model_registry.yaml \
        --model canonical_pp_v1 --input trees/pp_data.root --output-dir scores/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from features import FeatureMatrix, ModelSpec, feature_identity, features_for_file, load_registry, sha256  # noqa: E402

STATE_EVALUATED_FINITE = 1
STATE_EVALUATED_NONFINITE = 2
STATE_MISSING_INPUTS = 3
STATE_OUT_OF_DOMAIN = 4
STATE_MODEL_UNAVAILABLE = 5
STATE_NAMES = {1: "evaluated_finite", 2: "evaluated_nonfinite", 3: "missing_inputs",
               4: "out_of_domain", 5: "model_unavailable"}


def make_evaluator(spec: ModelSpec) -> Callable[[np.ndarray], np.ndarray] | None:
    """Return a function rows -> scores, or None when the model is not bound."""

    if not spec.bound:
        return None
    if not spec.model_file.is_file():
        raise SystemExit(f"model file not found: {spec.model_file}")
    actual = sha256(spec.model_file)
    if actual != spec.model_sha256:
        raise SystemExit(f"{spec.model_file}: sha256 {actual} differs from the registry's {spec.model_sha256}")
    k = len(spec.features)
    if spec.format == "tmva_rbdt":
        import ROOT  # PyROOT is needed only when a TMVA model is bound

        model = ROOT.TMVA.Experimental.RBDT("myBDT", str(spec.model_file))

        def evaluate(rows: np.ndarray) -> np.ndarray:
            if rows.ndim != 2 or rows.shape[1] != k:
                raise ValueError(f"model expects {k} inputs, got {rows.shape}")
            out = np.empty(len(rows))
            for i, row in enumerate(rows):
                values = ROOT.std.vector("float")()
                for v in row:
                    values.push_back(float(v))
                result = model.Compute(values)
                if len(result) != 1:
                    raise ValueError("model returned more than one score")
                out[i] = float(result[0])
            return out

        return evaluate
    if spec.format == "xgboost_json":
        import xgboost

        booster = xgboost.Booster()
        booster.load_model(str(spec.model_file))

        def evaluate(rows: np.ndarray) -> np.ndarray:
            if rows.ndim != 2 or rows.shape[1] != k:
                raise ValueError(f"model expects {k} inputs, got {rows.shape}")
            return np.asarray(booster.predict(xgboost.DMatrix(rows.astype(np.float32))), dtype=float)

        return evaluate
    raise SystemExit(f"unknown model format {spec.format!r}")


def score_matrix(matrix: FeatureMatrix, evaluator) -> tuple[np.ndarray, np.ndarray]:
    """Return (score, state) for every candidate of the feature matrix."""

    n = len(matrix.keys)
    score = np.full(n, np.nan)
    state = np.full(n, STATE_MODEL_UNAVAILABLE, dtype=np.int32)
    if evaluator is None:
        return score, state
    state[:] = STATE_OUT_OF_DOMAIN
    candidates = matrix.in_domain
    state[candidates & ~matrix.complete] = STATE_MISSING_INPUTS
    ready = candidates & matrix.complete
    if ready.any():
        values = evaluator(matrix.values[ready])
        if np.asarray(values).shape != (int(ready.sum()),):
            raise ValueError("model did not return one score per candidate")
        score[ready] = values
        finite = np.isfinite(values)
        state[np.flatnonzero(ready)[finite]] = STATE_EVALUATED_FINITE
        state[np.flatnonzero(ready)[~finite]] = STATE_EVALUATED_NONFINITE
        score[np.flatnonzero(ready)[~finite]] = np.nan
    return score, state


def write_sidecar(output: Path, matrix: FeatureMatrix, spec: ModelSpec, score: np.ndarray, state: np.ndarray,
                  receipt: dict[str, Any]) -> None:
    keys = np.asarray(matrix.keys, dtype=np.uint64).reshape(-1, 4) if matrix.keys else np.zeros((0, 4), dtype=np.uint64)
    columns = {
        "source_hi": matrix.source[:, 0].astype(np.uint64) if len(matrix.keys) else np.zeros(0, np.uint64),
        "source_lo": matrix.source[:, 1].astype(np.uint64) if len(matrix.keys) else np.zeros(0, np.uint64),
        "event_hi": keys[:, 0], "event_lo": keys[:, 1], "photon_hi": keys[:, 2], "photon_lo": keys[:, 3],
        "model": np.asarray([spec.name] * len(matrix.keys)),
        "score": score.astype(np.float64),
        "state": state.astype(np.int32),
        "complete": matrix.complete.astype(np.int32),
        "in_domain": matrix.in_domain.astype(np.int32),
        "shower_valid": matrix.shower_valid.astype(np.int32),
    }
    for j, name in enumerate(matrix.names):
        columns[f"input_{j:02d}"] = matrix.values[:, j].astype(np.float64)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    with uproot.create(output) as root:
        types = {name: ("string" if name == "model" else value.dtype) for name, value in columns.items()}
        root.mktree("PhotonScores", types)
        if len(matrix.keys):
            root["PhotonScores"].extend(columns)
        root["metadata"] = json.dumps(receipt, indent=2, sort_keys=True)


def augment_file(input_path: Path, output_dir: Path, spec: ModelSpec, evaluator) -> dict[str, Any]:
    input_hash = sha256(input_path)
    matrix = features_for_file(input_path, spec)
    score, state = score_matrix(matrix, evaluator)
    counts = {STATE_NAMES[s]: int((state == s).sum()) for s in STATE_NAMES}
    receipt = {
        "schema": "PhotonScoresSidecarV1",
        "feature_definition_sha256": feature_identity(spec),
        "model": spec.name,
        "system": spec.system,
        "shower_definition": spec.shower_definition,
        "features": list(spec.features),
        "ratio_policy": spec.ratio_policy,
        "domain": dict(spec.domain),
        "model_file": None if spec.model_file is None else str(spec.model_file),
        "model_sha256": spec.model_sha256,
        "score_direction": spec.score_direction,
        "input": {"path": str(input_path), "sha256": input_hash},
        "candidates": len(matrix.keys),
        "states": counts,
    }
    output = output_dir / f"{input_path.stem}.{spec.name}.scores.root"
    if sha256(input_path) != input_hash:
        raise ValueError("base product changed during scoring")
    write_sidecar(output, matrix, spec, score, state, receipt)
    receipt["output"] = {"path": str(output), "sha256": sha256(output)}
    output.with_suffix(".json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def read_inputs(inputs: list[Path], lists: list[Path]) -> list[Path]:
    paths = list(inputs)
    for listing in lists:
        for line in Path(listing).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                p = Path(line.split()[0])
                paths.append(p if p.is_absolute() else Path(listing).parent / p)
    if not paths:
        raise SystemExit("no input files given (use --input or --input-list)")
    if len({p.resolve() for p in paths}) != len(paths):
        raise ValueError("duplicate base input")
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--registry", required=True, type=Path, help="PhotonID/model_registry.yaml")
    parser.add_argument("--model", required=True, help="registry entry to evaluate")
    parser.add_argument("--input", action="append", default=[], type=Path, help="canonical tree file; repeatable")
    parser.add_argument("--input-list", action="append", default=[], type=Path, help="text file with one path per line")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    registry = load_registry(args.registry)
    if args.model not in registry:
        raise SystemExit(f"model {args.model!r} is not in {args.registry}; known: {sorted(registry)}")
    spec = registry[args.model]
    evaluator = make_evaluator(spec)
    if evaluator is None:
        print(f"[augment] {spec.name}: no bound model file; every candidate is written as model_unavailable",
              file=sys.stderr)
    for path in read_inputs(args.input, args.input_list):
        receipt = augment_file(path, args.output_dir, spec, evaluator)
        print(f"[augment] {receipt['output']['path']}: {receipt['states']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
