#!/usr/bin/env python3
"""Train the per-system photon-identification BDT from collaborator trees.

Input:  unscored or scored collaborator trees for one collision system:
        ``--signal`` photon+jet simulation files and ``--background``
        inclusive-jet simulation files.
Output: ``<output-dir>/<system>_photon_bdt.root`` (TMVA RBDT for score_trees.py),
        the XGBoost JSON model, and ``<output-dir>/<system>_model_manifest.json``
        with the ordered features, label mapping, sample lists, split counts,
        held-out AUC and file hashes.

Labels
  signal      reconstructed candidates in the signal files linked
              (recoTruthLinks, link_class 0, photon to photon) to a truth photon
              satisfying the nominal truth-signal contract in config/nominal.yaml
              (valid Geant photon, valid generator association, PPG12 class in
              {-1,0,1,2}, truth isolation R=0.3 below 4 GeV, |eta| < 0.7).
  background  every candidate with complete inputs in the background files that
              is not linked to such a truth photon.
  dropped     candidates in the signal files that are not linked to a signal
              truth photon.  They are not relabelled as background.

Features are the ordered model inputs stored per candidate
(``bdt_input_00`` ... ``bdt_input_NN``; 11 for p+p, 14 for Au+Au) and must be
complete (``bdt_input_count`` equal to the feature count, all finite).  Rows
are weighted by ``event_weight`` and the two classes are balanced to equal
total weight.  The split is by event identity (70/15/15), so candidates of
one event never straddle training and held-out partitions.

Requires ``xgboost`` and ``scikit-learn`` in addition to PyROOT.
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
sys.path.insert(0, str(HERE.parent / "FinalAnalysis"))
from photon_selection import FEATURE_COUNT, load_config, system_config, truth_signal_mask  # noqa: E402
from photonjet.training.split import assign_group_split  # noqa: E402

LINK_MATCH = 0
TYPE_PHOTON = 1
XGBOOST_DEFAULTS = {
    "n_estimators": 450, "max_depth": 4, "learning_rate": 0.035, "subsample": 0.85,
    "colsample_bytree": 0.85, "tree_method": "hist", "reg_alpha": 5.0, "reg_lambda": 0.3,
    "grow_policy": "lossguide", "max_bin": 256,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows(path: Path, *, system: str, config: dict[str, Any], role: str) -> dict[str, Any]:
    n_inputs = FEATURE_COUNT[system]
    branches = ["source_file_index", "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
                "photon_et", "photon_eta", "centrality", "event_weight", "bdt_input_count"]
    branches += [f"bdt_input_{j:02d}" for j in range(n_inputs)]
    with uproot.open(path) as root:
        photons = root["photons"].arrays(branches, library="np")
        truth = root["truthPhotons"].arrays(library="np")
        links = root["recoTruthLinks"].arrays(library="np")
    signal = truth_signal_mask(truth, config["truth_signal"])
    signal_truth = {
        (int(s), int(h), int(l))
        for s, h, l, ok in zip(truth["source_file_index"], truth["truth_photon_id_hi"], truth["truth_photon_id_lo"], signal)
        if ok
    }
    linked = set()
    for s, rh, rl, th, tl, rt, tt, lc in zip(
        links["source_file_index"], links["reco_id_hi"], links["reco_id_lo"], links["truth_id_hi"], links["truth_id_lo"],
        links["reco_type"], links["truth_type"], links["link_class"],
    ):
        if int(rt) == TYPE_PHOTON and int(tt) == TYPE_PHOTON and int(lc) == LINK_MATCH and (int(s), int(th), int(tl)) in signal_truth:
            linked.add((int(s), int(rh), int(rl)))
    is_linked = np.asarray([
        (int(s), int(h), int(l)) in linked
        for s, h, l in zip(photons["source_file_index"], photons["candidate_id_hi"], photons["candidate_id_lo"])
    ], dtype=bool)
    features = np.stack([np.asarray(photons[f"bdt_input_{j:02d}"], dtype=float) for j in range(n_inputs)], axis=1)
    complete = (np.asarray(photons["bdt_input_count"]) == n_inputs) & np.all(np.isfinite(features), axis=1)
    complete &= np.abs(np.asarray(photons["photon_eta"], dtype=float)) < float(config["photon"]["abs_eta_max"])
    keep = complete & (is_linked if role == "signal" else ~is_linked)
    label = np.full(int(keep.sum()), 1 if role == "signal" else 0, dtype=np.int32)
    group = np.asarray([f"{path.name}\x1f{s}\x1f{h}\x1f{l}" for s, h, l in zip(
        photons["source_file_index"][keep], photons["event_id_hi"][keep], photons["event_id_lo"][keep])], dtype=object)
    return {
        "features": features[keep], "label": label, "weight": np.asarray(photons["event_weight"], dtype=float)[keep],
        "group": group, "photon_et": np.asarray(photons["photon_et"], dtype=float)[keep],
        "n_complete": int(complete.sum()), "n_linked": int((complete & is_linked).sum()), "n_total": int(len(complete)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--signal", action="append", required=True, type=Path, help="photon+jet simulation tree file; repeatable")
    parser.add_argument("--background", action="append", required=True, type=Path, help="inclusive-jet simulation tree file; repeatable")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--n-jobs", type=int, default=4)
    for name, value in XGBOOST_DEFAULTS.items():
        parser.add_argument(f"--{name.replace('_', '-')}", type=type(value), default=value)
    args = parser.parse_args(argv)

    try:
        from sklearn.metrics import roc_auc_score
        from xgboost import XGBClassifier
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"training needs xgboost and scikit-learn ({exc.name} is not installed): "
            "pip install xgboost scikit-learn"
        ) from exc
    import ROOT

    config = load_config(args.config)
    block = system_config(config, args.system)
    features = list(block["features"])

    parts = [load_rows(p, system=args.system, config=config, role="signal") for p in args.signal]
    parts += [load_rows(p, system=args.system, config=config, role="background") for p in args.background]
    X = np.concatenate([p["features"] for p in parts])
    y = np.concatenate([p["label"] for p in parts])
    w = np.concatenate([p["weight"] for p in parts])
    groups = np.concatenate([p["group"] for p in parts])
    if (y == 1).sum() < 100 or (y == 0).sum() < 100:
        raise SystemExit(f"too few rows to train: signal={int((y == 1).sum())} background={int((y == 0).sum())}")
    if np.any(~np.isfinite(w)) or np.any(w <= 0):
        raise SystemExit("event_weight must be finite and positive for every training row")
    # class balance: equal total weight per class, mean weight one
    balanced = w.copy()
    for cls in (0, 1):
        balanced[y == cls] *= 0.5 * w.sum() / w[y == cls].sum()
    balanced /= balanced.mean()

    split = np.asarray(assign_group_split(groups.tolist(), seed=args.seed))
    train, valid, test = split == "train", split == "validation", split == "test"
    params = {name: getattr(args, name) for name in XGBOOST_DEFAULTS}
    model = XGBClassifier(random_state=args.seed, n_jobs=args.n_jobs, objective="binary:logistic", **params)
    model.fit(X[train], y[train], sample_weight=balanced[train],
              eval_set=[(X[valid], y[valid])], sample_weight_eval_set=[balanced[valid]], verbose=False)
    scores = model.predict_proba(X)[:, 1]
    auc = {name: float(roc_auc_score(y[mask], scores[mask], sample_weight=balanced[mask]))
           for name, mask in (("train", train), ("validation", valid), ("test", test))}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / f"{args.system}_photon_bdt.json"
    tmva_path = args.output_dir / f"{args.system}_photon_bdt.root"
    model.save_model(str(json_path))
    # TMVA export gives score_trees.py a PyROOT-only evaluator (same as the maintained analysis).
    ROOT.TMVA.Experimental.SaveXGBoost(model, "myBDT", str(tmva_path), num_inputs=len(features))

    manifest = {
        "schema": "PhotonIdModelManifestV1",
        "system": args.system,
        "selection_version": config.get("selection_version"),
        "features": features,
        "feature_count": len(features),
        "score_direction": "higher_is_signal",
        "tmva_model": {"path": str(tmva_path), "sha256": sha256(tmva_path)},
        "xgboost_model": {"path": str(json_path), "sha256": sha256(json_path)},
        "xgboost_params": params,
        "seed": args.seed,
        "label_mapping": {
            "signal": "signal files: candidate linked to a truth photon passing truth_signal in the configuration",
            "background": "background files: complete candidates not linked to such a truth photon",
            "dropped": "signal files: candidates not linked to a signal truth photon are excluded, not relabelled",
        },
        "truth_signal": config["truth_signal"],
        "weighting": "event_weight, classes balanced to equal total weight, mean one",
        "samples": [
            {"role": "signal" if i < len(args.signal) else "background",
             "path": str(p), "sha256": sha256(p), "rows_used": int(len(part["label"])),
             "candidates_total": part["n_total"], "candidates_complete": part["n_complete"],
             "candidates_linked_to_signal_truth": part["n_linked"]}
            for i, (p, part) in enumerate(zip([*args.signal, *args.background], parts))
        ],
        "split": {"train": int(train.sum()), "validation": int(valid.sum()), "test": int(test.sum())},
        "class_counts": {"signal": int((y == 1).sum()), "background": int((y == 0).sum())},
        "auc": auc,
    }
    (args.output_dir / f"{args.system}_model_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[train_photon_bdt] wrote {tmva_path} (AUC train={auc['train']:.4f} validation={auc['validation']:.4f} test={auc['test']:.4f})")
    print(f"[train_photon_bdt] set systems.{args.system}.model_file and model_sha256={manifest['tmva_model']['sha256']} in the configuration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
