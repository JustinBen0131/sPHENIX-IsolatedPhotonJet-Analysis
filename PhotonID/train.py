#!/usr/bin/env python3
"""Train from canonical trees using a fully bound, explicitly approved recipe.

features.py owns both training and inference inputs. Source/event hash groups
keep related candidates in one partition. Alternate productions with different
source identities require an explicit common grouping map before combining.
Unknown truth relations are excluded. Weight fitting and thinning use training
rows only; held-out AUC uses unit weights. This is an explicit supported recipe,
not a claim that the final pp/AuAu training recipe has been approved.

Registry entries remain unresolved until labels, weighting, split, window and
hyperparameters are reconciled with the accepted training authority. No silent
fallback recipe exists. Outputs are new model files plus a manifest; this
script never publishes models or rewrites the registry/base trees.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from features import FeatureMatrix, ModelSpec, feature_identity, features_for_file, load_registry, sha256, base_metadata  # noqa: E402

ASSOCIATION_MATCHED = 3  # photonjet::AssociationState::Matched

# ----------------------------------------------------------------------------
# Weighting
# ----------------------------------------------------------------------------

def inverse_pdf_weights(values: np.ndarray, *, n_bins: int = 20, fixed_range: tuple[float, float] | None = None,
                        cap: float | None = None) -> np.ndarray:
    """1 / (spline through a normalised histogram), capped, mean one."""

    from scipy.interpolate import UnivariateSpline

    x = np.asarray(values, dtype=float)
    if len(x) < n_bins:
        raise ValueError(f"fewer rows ({len(x)}) than bins ({n_bins}) for inverse-pdf weights")
    low, high = fixed_range if fixed_range is not None else (float(x.min()), float(x.max()))
    if not high > low:
        raise ValueError("degenerate range for inverse-pdf weights")
    edges = np.linspace(low, high, n_bins + 1)
    density = np.histogram(x, bins=edges, density=True)[0] * n_bins
    centres = 0.5 * (edges[:-1] + edges[1:])
    spline = UnivariateSpline(centres, density, s=0.0)
    pdf = np.clip(spline(x), 1.0e-3, None)
    local = 1.0 / pdf
    if cap is not None:
        local = np.minimum(local, float(cap))
    return local / local.mean()


def training_weights(et: np.ndarray, eta: np.ndarray, labels: np.ndarray, weighting: dict[str, Any]) -> tuple[np.ndarray, list[str]]:
    labels = np.asarray(labels, dtype=int)
    weights = np.ones(len(labels))
    steps = []
    n_total = len(labels)
    for cls in (0, 1):
        mask = labels == cls
        factor = n_total / (2.0 * mask.sum())
        weights[mask] *= factor
        steps.append(f"class {cls}: factor n_total/(2 n_class) = {factor:.6g}")
    eta_range = tuple(float(v) for v in weighting["eta_range"])
    bins = int(weighting["bins"])
    cap = float(weighting["et_cap"]) if weighting.get("et_cap") is not None else None
    for cls in (0, 1):
        mask = labels == cls
        weights[mask] *= inverse_pdf_weights(eta[mask], n_bins=bins, fixed_range=eta_range)
        weights[mask] *= inverse_pdf_weights(et[mask], n_bins=bins, fixed_range=None, cap=cap)
        steps.append(f"class {cls}: eta inverse-pdf ({bins} bins on {list(eta_range)}, mean one) "
                     f"then ET inverse-pdf ({bins} bins on the data range, cap {cap}, mean one)")
    return weights, steps


def flatten_low_et_background(et: np.ndarray, labels: np.ndarray, spec: dict[str, Any]) -> np.ndarray:
    keep = np.ones(len(et), dtype=bool)
    if not spec or not spec.get("enabled", False):
        return keep
    et_max = float(spec["et_max_gev"])
    bins = int(spec["bins"])
    rng = np.random.RandomState(int(spec["seed"]))
    low_bkg = (labels == 0) & (et < et_max)
    if not low_bkg.any():
        return keep
    edges = np.linspace(float(et[low_bkg].min()), et_max, bins + 1)
    index = np.clip(np.searchsorted(edges, et, side="right") - 1, 0, bins - 1)
    counts = np.bincount(index[low_bkg], minlength=bins)
    target = int(counts[counts > 0].min())
    for b in range(bins):
        rows = np.flatnonzero(low_bkg & (index == b))
        if len(rows) > target:
            keep[rng.choice(rows, size=len(rows) - target, replace=False)] = False
    return keep


# ----------------------------------------------------------------------------
# Labels from the canonical relations
# ----------------------------------------------------------------------------

def signal_linked_candidates(path: Path) -> tuple[set, set]:
    """Candidate keys whose Matched photon-truth relation points at an analysis-signal truth photon."""

    with uproot.open(path) as root:
        complete = root["Events"]["truth_denominator_complete"].array(library="np")
        if np.any(complete != 1):
            raise ValueError("training source has incomplete truth census")
        truth = root["TruthPhotons"].arrays(["event_hi", "event_lo", "truth_photon_hi", "truth_photon_lo", "analysis_signal"], library="np")
        links = root["PhotonTruthLinks"].arrays(["event_hi", "event_lo", "photon_hi", "photon_lo",
                                                 "truth_photon_hi", "truth_photon_lo", "state"], library="np")
    signal_truth = {
        (int(eh), int(el), int(th), int(tl))
        for eh, el, th, tl, ok in zip(truth["event_hi"], truth["event_lo"], truth["truth_photon_hi"], truth["truth_photon_lo"], truth["analysis_signal"])
        if bool(ok)
    }
    linked, known = set(), set()
    for eh, el, ph, pl, th, tl, state in zip(links["event_hi"], links["event_lo"], links["photon_hi"], links["photon_lo"],
                                             links["truth_photon_hi"], links["truth_photon_lo"], links["state"]):
        if int(state) in (ASSOCIATION_MATCHED, 4):
            known.add((int(eh), int(el), int(ph), int(pl)))
        if int(state) == ASSOCIATION_MATCHED and (int(eh), int(el), int(th), int(tl)) in signal_truth:
            linked.add((int(eh), int(el), int(ph), int(pl)))
    return linked, known


def load_rows(path: Path, spec: ModelSpec, role: str) -> dict[str, Any]:
    matrix: FeatureMatrix = features_for_file(path, spec)
    with uproot.open(path) as root:
        metadata = base_metadata(root)
        if metadata.get("data_kind") != "2" or metadata.get("simulation_role") != ("1" if role == "signal" else "2"):
            raise ValueError(f"{path}: training role differs from source metadata")
    linked, known = signal_linked_candidates(path)
    is_linked = np.asarray([key in linked for key in matrix.keys], dtype=bool)
    complete = matrix.complete & matrix.in_domain & np.isfinite(matrix.et) & np.isfinite(matrix.eta)
    complete &= np.abs(matrix.eta) < float(spec.domain.get("max_abs_eta", 0.7))
    keep = complete & np.asarray([key in known for key in matrix.keys]) & (is_linked if role == "signal" else ~is_linked)
    return {
        "features": matrix.values[keep],
        "groups": np.asarray([f"{int(src[0])}:{int(src[1])}:{key[0]}:{key[1]}"
                               for src, key in zip(matrix.source, matrix.keys)])[keep],
        "keys": [key for key, take in zip(matrix.keys, keep) if take],
        "label": np.full(int(keep.sum()), 1 if role == "signal" else 0, dtype=np.int32),
        "et": matrix.et[keep], "eta": matrix.eta[keep],
        "n_total": len(matrix.keys), "n_complete": int(complete.sum()), "n_linked": int((complete & is_linked).sum()),
    }


def split_rows(groups: np.ndarray, spec: dict[str, Any]) -> np.ndarray:
    """Stable source/event groups; candidate order and file partition do not enter."""
    seed = int(spec["seed"])
    test_fraction = float(spec["test_fraction"])
    validation_fraction = float(spec["validation_fraction"])
    if not (0 < test_fraction < 1 and 0 <= validation_fraction < 1 - test_fraction):
        raise ValueError("invalid held-out fractions")
    def fraction(group):
        digest = hashlib.sha256(f"photonid-split-v1:{seed}:{group}".encode()).digest()
        return int.from_bytes(digest[:8], "big") / 2**64
    draws = np.asarray([fraction(group) for group in groups])
    return np.where(draws < test_fraction, "test",
                    np.where(draws < test_fraction + validation_fraction, "validation", "train"))


def training_contract(spec: ModelSpec) -> dict[str, Any]:
    training = dict(spec.training)
    required = {"status", "recipe_evidence", "label_rule", "et_window_gev", "weighting",
                "low_et_background_flattening", "split", "xgboost", "weight_fit_population", "heldout_weights"}
    if required - training.keys() or training.get("status") != "approved" or not training.get("recipe_evidence"):
        raise ValueError("training recipe is unresolved; bind the reviewed recipe before training")
    if training["label_rule"] != "source_role_known_truth" or training["weight_fit_population"] != "train_only" or training["heldout_weights"] != "unit":
        raise ValueError("requested label/weight recipe is not implemented; do not substitute another")
    if training["split"].get("rule") != "source_event_hash_v1":
        raise ValueError("training requires source_event_hash_v1 grouping")
    for key in ("eta_range", "bins", "et_cap"):
        if key not in training["weighting"]: raise ValueError(f"missing weighting.{key}")
    thinning = training["low_et_background_flattening"]
    if "enabled" not in thinning:
        raise ValueError("background thinning must be explicitly enabled or disabled")
    if thinning["enabled"] and any(k not in thinning for k in ("et_max_gev", "bins", "seed")):
        raise ValueError("thinning parameters must be explicit")
    if not training["xgboost"]: raise ValueError("explicit XGBoost hyperparameters required")
    return training


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--model", required=True, help="registry entry to train (its features and system)")
    parser.add_argument("--signal", action="append", required=True, type=Path, help="photon+jet simulation tree; repeatable")
    parser.add_argument("--background", action="append", required=True, type=Path, help="inclusive-jet simulation tree; repeatable")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args(argv)

    try:
        from sklearn.metrics import roc_auc_score
        from xgboost import XGBClassifier
    except ModuleNotFoundError as exc:
        raise SystemExit(f"training needs xgboost and scikit-learn ({exc.name} is not installed)") from exc
    import ROOT

    registry = load_registry(args.registry)
    if args.model not in registry:
        raise SystemExit(f"model {args.model!r} is not in {args.registry}")
    spec = registry[args.model]
    training = training_contract(spec)
    params = dict(training["xgboost"])

    parts = [load_rows(p, spec, "signal") for p in args.signal]
    parts += [load_rows(p, spec, "background") for p in args.background]
    X = np.concatenate([p["features"] for p in parts])
    y = np.concatenate([p["label"] for p in parts])
    et = np.concatenate([p["et"] for p in parts])
    eta = np.concatenate([p["eta"] for p in parts])

    groups = np.concatenate([p["groups"] for p in parts])
    keys = [key for part in parts for key in part["keys"]]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate candidate identity across training inputs")
    low, high = (float(v) for v in training["et_window_gev"])
    if not (np.isfinite(low) and np.isfinite(high) and low < high):
        raise ValueError("invalid training window")
    window = (et >= low) & (et < high)
    X, y, et, eta, groups = X[window], y[window], et[window], eta[window], groups[window]
    split = split_rows(groups, training["split"])
    train, valid, test = split == "train", split == "validation", split == "test"
    # Fit weighting and thinning on training observations only. Held-out AUC
    # uses unit weights; this explicit recipe needs separate science approval.
    keep = np.ones(len(y), dtype=bool)
    keep[train] = flatten_low_et_background(et[train], y[train], training["low_et_background_flattening"])
    X, y, et, eta, split = X[keep], y[keep], et[keep], eta[keep], split[keep]
    train, valid, test = split == "train", split == "validation", split == "test"
    for name, mask in (("train", train), ("test", test), ("validation", valid)):
        if (name != "validation" or valid.any()) and len(np.unique(y[mask])) != 2:
            raise ValueError(f"{name} partition lacks both classes; revise the declared split")
    weights = np.ones(len(y))
    weights[train], steps = training_weights(et[train], eta[train], y[train], training["weighting"])

    seed = int(params.pop("random_state", (training.get("split") or {}).get("seed", 42)))
    if args.output_dir.exists():
        raise FileExistsError("training output directory must be new")
    model = XGBClassifier(random_state=seed, n_jobs=args.n_jobs, objective="binary:logistic", **params)
    eval_set = [(X[valid], y[valid])] if valid.any() else None
    eval_weights = [weights[valid]] if valid.any() else None
    model.fit(X[train], y[train], sample_weight=weights[train], eval_set=eval_set,
              sample_weight_eval_set=eval_weights, verbose=False)
    scores = model.predict_proba(X)[:, 1]
    auc = {name: float(roc_auc_score(y[mask], scores[mask], sample_weight=weights[mask]))
           for name, mask in (("train", train), ("validation", valid), ("test", test)) if mask.any()}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / f"{spec.name}.json"
    tmva_path = args.output_dir / f"{spec.name}.root"
    model.save_model(str(json_path))
    booster = model.get_booster()
    booster.feature_names = [f"f{i}" for i in range(len(spec.features))]
    ROOT.TMVA.Experimental.SaveXGBoost(model, "myBDT", str(tmva_path), num_inputs=len(spec.features))

    manifest = {
        "schema": "PhotonIdModelManifestV2",
        "feature_definition_sha256": feature_identity(spec),
        "training_contract": training,
        "model": spec.name,
        "system": spec.system,
        "shower_definition": spec.shower_definition,
        "features": list(spec.features),
        "ratio_policy": spec.ratio_policy,
        "score_direction": "higher_is_signal",
        "tmva_model": {"path": str(tmva_path), "sha256": sha256(tmva_path)},
        "xgboost_model": {"path": str(json_path), "sha256": sha256(json_path)},
        "xgboost_params": {**params, "random_state": seed, "objective": "binary:logistic"},
        "label_rule": {
            "signal": "signal files: candidate with a Matched PhotonTruthLinks relation to a TruthPhotons row with analysis_signal",
            "background": "background files: known Matched non-signal or Fake relation; Unknown excluded",
            "dropped": "cross-role candidates are excluded, never relabelled; no isolation or working point enters",
        },
        "weighting_steps": steps,
        "producer_event_weight_used": False,
        "low_et_background_flattening": training.get("low_et_background_flattening"),
        "training_et_window_gev": [low, high],
        "split": {"rule": "source_event_hash_v1", **(training.get("split") or {}),
                  "train": int(train.sum()), "validation": int(valid.sum()), "test": int(test.sum())},
        "samples": [{"role": "signal" if i < len(args.signal) else "background", "path": str(p), "sha256": sha256(p),
                     "rows_used": int(len(part["label"])), "candidates_total": part["n_total"],
                     "candidates_complete": part["n_complete"], "candidates_linked_to_signal_truth": part["n_linked"]}
                    for i, (p, part) in enumerate(zip([*args.signal, *args.background], parts))],
        "class_counts": {"signal": int((y == 1).sum()), "background": int((y == 0).sum())},
        "auc": auc,
    }
    (args.output_dir / f"{spec.name}.manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[train] wrote {tmva_path} AUC={auc}")
    print(f"[train] bind model_file and model_sha256={manifest['tmva_model']['sha256']} for {spec.name} in the registry")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
