#!/usr/bin/env python3
"""One immutable-base, multi-model ROOT sidecar writer.

Model inference remains shared with training. Model, WP and isolation identities
are separate: raw scores do not require a selection recipe. All relations carry
source/event/photon words, never input row numbers.
"""

from pathlib import Path
import hashlib
import json
import os
import sys
import tempfile
from typing import Callable
import numpy as np
import uproot
from features import (
    FeatureMatrix,
    ModelSpec,
    feature_identity,
    features_for_file,
    sha256,
    base_metadata,
)
from registry import digest, file_hash
from inputs import isolation_for_file
from calibration import thresholds, line, id_regions, isolation_regions

STATE_EVALUATED_FINITE = 1
STATE_EVALUATED_NONFINITE = 2
STATE_MISSING_INPUTS = 3
STATE_OUT_OF_DOMAIN = 4
STATE_MODEL_UNAVAILABLE = 5
STATE_NAMES = {
    1: "evaluated_finite",
    2: "evaluated_nonfinite",
    3: "missing_inputs",
    4: "out_of_domain",
    5: "model_unavailable",
}


def make_evaluator(spec: ModelSpec) -> Callable[[np.ndarray], np.ndarray] | None:
    """Return a function rows -> scores, or None when the model is not bound."""

    if not spec.bound:
        return None
    if not spec.model_file.is_file():
        raise ValueError(f"model file not found: {spec.model_file}")
    actual = sha256(spec.model_file)
    if actual != spec.model_sha256:
        raise ValueError(
            f"{spec.model_file}: sha256 {actual} differs from the registry's {spec.model_sha256}"
        )
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
            return np.asarray(
                booster.predict(
                    xgboost.DMatrix(
                        rows.astype(np.float32), feature_names=list(spec.features)
                    )
                ),
                dtype=float,
            )

        return evaluate
    raise ValueError(f"unknown model format {spec.format!r}")


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


KEY_TYPES = {
    k: "uint64"
    for k in (
        "source_hi",
        "source_lo",
        "event_hi",
        "event_lo",
        "photon_hi",
        "photon_lo",
    )
}
MODEL_TYPES = dict(
    model_id="uint64",
    model_name="string",
    model_version="string",
    collision_system="int32",
    model_sha256="string",
    feature_schema_sha256="string",
    training_receipt_sha256="string",
)
SCORE_TYPES = dict(
    KEY_TYPES,
    model_id="uint64",
    score="float64",
    score_valid="int32",
    state="int32",
    complete="int32",
    in_domain="int32",
    shower_valid="int32",
)
RECIPE_TYPES = dict(
    selection_id="uint64",
    selection_name="string",
    model_id="uint64",
    id_package_sha256="string",
    isolation_package_sha256="string",
)
SELECTION_TYPES = dict(
    KEY_TYPES,
    selection_id="uint64",
    model_id="uint64",
    bdt_score="float64",
    tight_threshold="float64",
    nontight_upper_threshold="float64",
    nontight_lower_threshold="float64",
    id_region="int32",
    isolation_value="float64",
    isolation_radius="float64",
    isolation_method="int32",
    isolated_threshold="float64",
    nonisolated_threshold="float64",
    isolation_region="int32",
    valid="int32",
)


def stable_id(kind, identity):
    return int.from_bytes(
        hashlib.sha256((kind + ":" + digest(identity)).encode()).digest()[:8], "big"
    )


def key_columns(matrix):
    a = np.column_stack(
        (matrix.source, np.asarray(matrix.keys, dtype=np.uint64).reshape(-1, 4))
    ).astype(np.uint64)
    if len(np.unique(a, axis=0)) != len(a):
        raise ValueError("duplicate sidecar candidate identity")
    return {key: a[:, i] for i, key in enumerate(KEY_TYPES)}


def selection_columns(matrix, scores, state, recipe, model_id):
    wp, iso = recipe["wp"], recipe["isolation"]
    cuts = thresholds(wp, matrix.et, matrix.centrality)
    identity = dict(
        name=recipe["name"],
        model_id=int(model_id),
        wp=recipe["wp_sha256"],
        isolation=recipe["isolation_sha256"],
    )
    sid = stable_id("selection", identity)
    n = len(scores)
    region = id_regions(scores, cuts["WP70"], cuts["WP80"], wp["non_tight_lower"])
    region[state != 1] = 0
    axis = matrix.centrality if iso["axis"] == "centrality" else matrix.et
    lo, hi = iso["domain"]
    valid_domain = np.isfinite(axis) & (axis >= lo) & (axis < hi)
    isolated = np.where(valid_domain, line(iso["isolated"], axis), np.nan)
    nonisolated = np.where(valid_domain, line(iso["nonisolated"], axis), np.nan)
    value = recipe["_isolation_values"]
    ir = isolation_regions(value, isolated, nonisolated)
    cols = dict(
        key_columns(matrix),
        selection_id=np.full(n, sid, dtype=np.uint64),
        model_id=np.full(n, model_id, dtype=np.uint64),
        bdt_score=scores,
        tight_threshold=cuts["WP70"],
        nontight_upper_threshold=cuts["WP80"],
        nontight_lower_threshold=np.full(n, wp["non_tight_lower"]),
        id_region=region,
        isolation_value=value,
        isolation_radius=np.full(n, iso["radius"]),
        isolation_method=np.full(n, iso["method"], dtype=np.int32),
        isolated_threshold=isolated,
        nonisolated_threshold=nonisolated,
        isolation_region=ir,
        valid=((region != 0) & (ir != 0)).astype(np.int32),
    )
    row = dict(
        selection_id=sid,
        selection_name=recipe["name"],
        model_id=model_id,
        id_package_sha256=recipe["wp_sha256"],
        isolation_package_sha256=recipe["isolation_sha256"],
    )
    return cols, row


def write_table(root, name, types, blocks):
    root.mktree(name, types)
    for block in blocks:
        if not block or not len(next(iter(block.values()))):
            continue
        root[name].extend(
            {
                k: np.asarray(v, dtype=str if types[k] == "string" else types[k])
                for k, v in block.items()
            }
        )


def write_sidecar(output, score_blocks, models, selection_blocks, recipes, receipt):
    """Explicit TTrees, exclusive atomic publication, count/schema readback."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    fd, name = tempfile.mkstemp(
        prefix="." + output.name + ".", suffix=".partial", dir=output.parent
    )
    os.close(fd)
    tmp = Path(name)
    try:
        with uproot.recreate(tmp) as f:
            write_table(
                f,
                "Models",
                MODEL_TYPES,
                [{k: [r[k] for r in models] for k in MODEL_TYPES}],
            )
            write_table(f, "PhotonScores", SCORE_TYPES, score_blocks)
            write_table(
                f,
                "SelectionRecipes",
                RECIPE_TYPES,
                [{k: [r[k] for r in recipes] for k in RECIPE_TYPES}],
            )
            write_table(f, "PhotonSelections", SELECTION_TYPES, selection_blocks)
            f["metadata"] = json.dumps(receipt, sort_keys=True, allow_nan=False)
        with uproot.open(tmp) as f:
            for table, types, n in [
                ("Models", MODEL_TYPES, len(models)),
                (
                    "PhotonScores",
                    SCORE_TYPES,
                    sum(len(b["score"]) for b in score_blocks),
                ),
                ("SelectionRecipes", RECIPE_TYPES, len(recipes)),
                (
                    "PhotonSelections",
                    SELECTION_TYPES,
                    sum(len(b["bdt_score"]) for b in selection_blocks),
                ),
            ]:
                if (
                    f[table].classname != "TTree"
                    or set(f[table].keys()) != set(types)
                    or f[table].num_entries != n
                ):
                    raise ValueError("sidecar schema/count readback failed")
            for table, blocks in [
                ("PhotonScores", score_blocks),
                ("PhotonSelections", selection_blocks),
            ]:
                for key in KEY_TYPES:
                    expected = (
                        np.concatenate([b[key] for b in blocks])
                        if blocks
                        else np.empty(0, dtype=np.uint64)
                    )
                    if not np.array_equal(f[table][key].array(library="np"), expected):
                        raise ValueError("identity changed during serialization")
        # Link creation is exclusive even if another process published meanwhile.
        os.link(tmp, output)
    finally:
        tmp.unlink(missing_ok=True)


def augment_file(
    source, output, specs, evaluators, recipes, registry_rows, campaign_identity
):
    path = source["_path"]
    before = file_hash(path)
    if before != source["sha256"]:
        raise ValueError("base changed since campaign preflight")
    with uproot.open(path) as f:
        meta = base_metadata(f)
        s = f["Sources"].arrays(["source_hi", "source_lo"], library="np")
        if [int(s["source_hi"][0]), int(s["source_lo"][0])] != source["source_id"]:
            raise ValueError("manifest/ROOT source mismatch")
        expected_kind = "1" if source["lane"].endswith("_data") else "2"
        if meta.get("data_kind") != expected_kind:
            raise ValueError("manifest/ROOT data kind mismatch")
        if expected_kind == "2" and meta.get("simulation_role") != (
            "1" if "photon" in source["lane"] else "2"
        ):
            raise ValueError("manifest/ROOT simulation role mismatch")
    scores_out, selections_out, model_rows, recipe_rows = [], [], [], []
    ids = set()
    selection_ids = set()
    common = None
    counts = {}
    for spec in specs:
        if not spec.bound or evaluators[spec.name] is None:
            raise ValueError(f"unbound model {spec.name}")
        m = features_for_file(path, spec)
        if not np.all(m.source == np.asarray(source["source_id"], dtype=np.uint64)):
            raise ValueError("event source identity differs from manifest/Source table")
        keys = key_columns(m)
        if common is not None and any(
            not np.array_equal(common[k], keys[k]) for k in keys
        ):
            raise ValueError("model feature adapters disagree on candidate identities")
        common = keys
        mid = stable_id(
            "model",
            dict(
                name=spec.name, sha=spec.model_sha256, features=feature_identity(spec)
            ),
        )
        if mid in ids:
            raise ValueError("duplicate/colliding model ID")
        ids.add(mid)
        score, state = score_matrix(m, evaluators[spec.name])
        scores_out.append(
            dict(
                keys,
                model_id=np.full(len(score), mid, dtype=np.uint64),
                score=score,
                score_valid=(state == 1).astype(np.int32),
                state=state,
                complete=m.complete.astype(np.int32),
                in_domain=m.in_domain.astype(np.int32),
                shower_valid=m.shower_valid.astype(np.int32),
            )
        )
        counts[spec.name] = {
            STATE_NAMES[s]: int((state == s).sum()) for s in STATE_NAMES
        }
        model_rows.append(
            dict(
                model_id=mid,
                model_name=spec.name,
                model_version=str(registry_rows[spec.name]["version"]),
                collision_system=1 if spec.system == "pp" else 2,
                model_sha256=spec.model_sha256,
                feature_schema_sha256=feature_identity(spec),
                training_receipt_sha256=file_hash(spec.training_manifest),
            )
        )
        for recipe in recipes:
            if recipe["model"] != spec.name:
                continue
            r = dict(
                recipe,
                _isolation_values=isolation_for_file(
                    path,
                    m,
                    recipe["isolation"]["method"],
                    recipe["isolation"]["radius"],
                ),
            )
            cols, row = selection_columns(m, score, state, r, mid)
            if row["selection_id"] in selection_ids:
                raise ValueError("duplicate/colliding selection ID")
            selection_ids.add(row["selection_id"])
            selections_out.append(cols)
            recipe_rows.append(row)
    receipt = dict(
        schema="PhotonIDSidecarV2",
        status="complete",
        scientific_acceptance="NOT_REVIEWED",
        campaign=campaign_identity,
        base_sha256=before,
        source_id=source["source_id"],
        lane=source["lane"],
        states=counts,
        models=model_rows,
        recipes=recipes,
        join_key=list(KEY_TYPES),
        id_regions={
            "INVALID": 0,
            "BELOW_FLOOR": 1,
            "NONTIGHT": 2,
            "EXCLUDED": 3,
            "TIGHT": 4,
        },
        isolation_regions={"INVALID": 0, "ISOLATED": 1, "GAP": 2, "NONISOLATED": 3},
    )
    if file_hash(path) != before:
        raise ValueError("base changed during scoring")
    write_sidecar(output, scores_out, model_rows, selections_out, recipe_rows, receipt)
    if file_hash(path) != before:
        Path(output).unlink()  # only our just-created, invalid sidecar; never the base
        raise ValueError("base changed during sidecar publication")
    return dict(output=str(output), sha256=file_hash(output), states=counts)


if __name__ == "__main__":
    from frontend import main

    raise SystemExit(main("augment", sys.argv[1:]))
