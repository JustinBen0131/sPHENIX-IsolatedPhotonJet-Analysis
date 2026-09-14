#!/usr/bin/env python3
"""Train the per-system photon-identification BDT with the established recipe.

Input:  collaborator trees for one collision system: ``--signal`` photon+jet
        simulation files and ``--background`` inclusive-jet simulation files.
Output: ``<output-dir>/<system>_photon_bdt.root`` (TMVA RBDT read by
        score_trees.py), the XGBoost JSON model, and
        ``<output-dir>/<system>_model_manifest.json`` with the ordered
        features, the label mapping, the weighting steps, sample hashes, split
        counts and the held-out AUC.

Labels (config ``truth_signal``)
  signal      candidates in the signal files linked (photon-to-photon match)
              to a truth photon satisfying the nominal truth-signal contract;
  background  candidates in the background files not linked to such a photon;
  dropped     signal-file candidates without such a link.  They are never
              relabelled as background.

Weights (config ``systems.<system>.training.weighting``), computed on every
row with complete inputs before the training window is applied, as in the
maintained analysis:
  1. class factor n_total / (2 n_class), so both classes carry equal weight;
  2. per class, eta flattening: inverse of a spline through a 20-bin density
     on the fixed range [-0.7, 0.7], normalised to mean one;
  3. per class, ET flattening: the same with 20 bins over the data range and
     the weight capped at 800 before normalisation.
The producer event weight is not used in training.  p+p additionally thins
background rows below 15 GeV to equal counts per bin (the historical
low-ET background flattening); Au+Au does not.

Split: row-level, stratified on the label, with the per-system fractions and
seeds from the configuration (p+p 70/10/20 seeds 42/43; Au+Au 90/10 seed 13).

Hyperparameters: ``systems.<system>.training.xgboost``.  Requires ``xgboost``,
``scikit-learn`` and ``scipy`` in addition to PyROOT (for the TMVA export).
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
from photon_selection import FEATURE_COUNT, load_config, system_config, truth_signal_mask  # noqa: E402

LINK_MATCH = 0
TYPE_PHOTON = 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# Weighting (the maintained PPG12-style recipe)
# --------------------------------------------------------------------------
def inverse_pdf_weights(values: np.ndarray, *, n_bins: int = 20, fixed_range: tuple[float, float] | None = None,
                        cap: float | None = None) -> np.ndarray:
    """1 / (spline through a normalised histogram), capped, mean one.

    The histogram density is multiplied by the bin count so a flat sample has
    density one; the spline through the bin centres is evaluated at every
    value and clipped at 1e-3 before inversion.
    """

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
    """Class factor, then per-class eta and ET flattening; returns weights and a step log."""

    labels = np.asarray(labels, dtype=int)
    weights = np.ones(len(labels))
    steps = []
    n_total = len(labels)
    for cls in (0, 1):
        mask = labels == cls
        factor = n_total / (2.0 * mask.sum())
        weights[mask] *= factor
        steps.append(f"class {cls}: factor n_total/(2 n_class) = {factor:.6g}")
    eta_range = tuple(float(v) for v in weighting.get("eta_range", (-0.7, 0.7)))
    bins = int(weighting.get("bins", 20))
    cap = float(weighting["et_cap"]) if weighting.get("et_cap") is not None else None
    for cls in (0, 1):
        mask = labels == cls
        weights[mask] *= inverse_pdf_weights(eta[mask], n_bins=bins, fixed_range=eta_range)
        weights[mask] *= inverse_pdf_weights(et[mask], n_bins=bins, fixed_range=None, cap=cap)
        steps.append(f"class {cls}: eta inverse-pdf ({bins} bins on {list(eta_range)}, no cap, mean one) "
                     f"then ET inverse-pdf ({bins} bins on the data range, cap {cap}, mean one)")
    return weights, steps


def flatten_low_et_background(et: np.ndarray, labels: np.ndarray, spec: dict[str, Any]) -> np.ndarray:
    """Keep mask thinning background rows below ``et_max_gev`` to equal counts per bin."""

    keep = np.ones(len(et), dtype=bool)
    if not spec or not spec.get("enabled", False):
        return keep
    et_max = float(spec.get("et_max_gev", 15.0))
    bins = int(spec.get("bins", 20))
    rng = np.random.RandomState(int(spec.get("seed", 42)))
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


# --------------------------------------------------------------------------
# Rows
# --------------------------------------------------------------------------
def load_rows(path: Path, *, system: str, config: dict[str, Any], role: str) -> dict[str, Any]:
    n_inputs = FEATURE_COUNT[system]
    branches = ["source_file_index", "event_id_hi", "event_id_lo", "candidate_id_hi", "candidate_id_lo",
                "photon_et", "photon_eta", "centrality", "bdt_input_count"]
    branches += [f"bdt_input_{j:02d}" for j in range(n_inputs)]
    with uproot.open(path) as root:
        photons = root["photons"].arrays(branches, library="np")
        truth = root["truthPhotons"].arrays(list(root["truthPhotons"].keys()), library="np")
        links = root["recoTruthLinks"].arrays(list(root["recoTruthLinks"].keys()), library="np")
    signal = truth_signal_mask(truth, config["truth_signal"])
    signal_truth = {(int(s), int(h), int(l)) for s, h, l, ok in
                    zip(truth["source_file_index"], truth["truth_photon_id_hi"], truth["truth_photon_id_lo"], signal) if ok}
    linked = set()
    for s, rh, rl, th, tl, rt, tt, lc in zip(links["source_file_index"], links["reco_id_hi"], links["reco_id_lo"],
                                             links["truth_id_hi"], links["truth_id_lo"], links["reco_type"],
                                             links["truth_type"], links["link_class"]):
        if int(rt) == TYPE_PHOTON and int(tt) == TYPE_PHOTON and int(lc) == LINK_MATCH and (int(s), int(th), int(tl)) in signal_truth:
            linked.add((int(s), int(rh), int(rl)))
    is_linked = np.asarray([(int(s), int(h), int(l)) in linked for s, h, l in
                            zip(photons["source_file_index"], photons["candidate_id_hi"], photons["candidate_id_lo"])], dtype=bool)
    features = np.stack([np.asarray(photons[f"bdt_input_{j:02d}"], dtype=float) for j in range(n_inputs)], axis=1)
    complete = (np.asarray(photons["bdt_input_count"]) == n_inputs) & np.all(np.isfinite(features), axis=1)
    complete &= np.abs(np.asarray(photons["photon_eta"], dtype=float)) < float(config["photon"]["abs_eta_max"])
    keep = complete & (is_linked if role == "signal" else ~is_linked)
    return {
        "features": features[keep],
        "label": np.full(int(keep.sum()), 1 if role == "signal" else 0, dtype=np.int32),
        "photon_et": np.asarray(photons["photon_et"], dtype=float)[keep],
        "photon_eta": np.asarray(photons["photon_eta"], dtype=float)[keep],
        "n_total": int(len(complete)), "n_complete": int(complete.sum()), "n_linked": int((complete & is_linked).sum()),
    }


def split_rows(labels: np.ndarray, spec: dict[str, Any]) -> np.ndarray:
    """Row-level stratified split -> array of 'train' / 'validation' / 'test'."""

    from sklearn.model_selection import train_test_split

    seed = int(spec.get("seed", 42))
    test_fraction = float(spec.get("test_fraction", 0.2))
    validation_fraction = float(spec.get("validation_fraction", 0.0))
    index = np.arange(len(labels))
    parts = np.full(len(labels), "train", dtype=object)
    rest, test = train_test_split(index, test_size=test_fraction, random_state=seed, stratify=labels)
    parts[test] = "test"
    if validation_fraction > 0:
        adjusted = validation_fraction / (1.0 - test_fraction)
        _, validation = train_test_split(rest, test_size=adjusted, random_state=seed + 1, stratify=labels[rest])
        parts[validation] = "validation"
    return parts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--signal", action="append", required=True, type=Path, help="photon+jet simulation tree file; repeatable")
    parser.add_argument("--background", action="append", required=True, type=Path, help="inclusive-jet simulation tree file; repeatable")
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args(argv)

    try:
        from sklearn.metrics import roc_auc_score
        from xgboost import XGBClassifier
    except ModuleNotFoundError as exc:
        raise SystemExit(f"training needs xgboost and scikit-learn ({exc.name} is not installed): pip install xgboost scikit-learn scipy") from exc
    import ROOT

    config = load_config(args.config)
    block = system_config(config, args.system)
    training = block.get("training") or {}
    features = list(block["features"])
    params = dict(training.get("xgboost") or {})
    if not params:
        raise SystemExit(f"systems.{args.system}.training.xgboost is not set in the configuration")

    parts = [load_rows(p, system=args.system, config=config, role="signal") for p in args.signal]
    parts += [load_rows(p, system=args.system, config=config, role="background") for p in args.background]
    X = np.concatenate([p["features"] for p in parts])
    y = np.concatenate([p["label"] for p in parts])
    et = np.concatenate([p["photon_et"] for p in parts])
    eta = np.concatenate([p["photon_eta"] for p in parts])

    keep = flatten_low_et_background(et, y, training.get("low_et_background_flattening") or {})
    X, y, et, eta = X[keep], y[keep], et[keep], eta[keep]
    if (y == 1).sum() < 100 or (y == 0).sum() < 100:
        raise SystemExit(f"too few rows to train: signal={int((y == 1).sum())} background={int((y == 0).sum())}")
    weights, steps = training_weights(et, eta, y, training.get("weighting") or {})

    low, high = (float(v) for v in training.get("et_window_gev", (float(et.min()), float(et.max()) + 1.0)))
    window = (et >= low) & (et < high)
    X, y, weights = X[window], y[window], weights[window]
    split = split_rows(y, training.get("split") or {})
    train, valid, test = split == "train", split == "validation", split == "test"

    seed = int(params.pop("random_state", (training.get("split") or {}).get("seed", 42)))
    model = XGBClassifier(random_state=seed, n_jobs=args.n_jobs, objective="binary:logistic", **params)
    eval_set = [(X[valid], y[valid])] if valid.any() else [(X[test], y[test])]
    eval_weights = [weights[valid]] if valid.any() else [weights[test]]
    model.fit(X[train], y[train], sample_weight=weights[train], eval_set=eval_set, sample_weight_eval_set=eval_weights, verbose=False)
    scores = model.predict_proba(X)[:, 1]
    auc = {name: float(roc_auc_score(y[mask], scores[mask], sample_weight=weights[mask]))
           for name, mask in (("train", train), ("validation", valid), ("test", test)) if mask.any()}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / f"{args.system}_photon_bdt.json"
    tmva_path = args.output_dir / f"{args.system}_photon_bdt.root"
    model.save_model(str(json_path))
    booster = model.get_booster()
    booster.feature_names = [f"f{i}" for i in range(len(features))]
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
        "xgboost_params": {**params, "random_state": seed, "objective": "binary:logistic"},
        "label_mapping": {
            "signal": "signal files: candidate linked to a truth photon passing truth_signal in the configuration",
            "background": "background files: complete candidates not linked to such a truth photon",
            "dropped": "signal files: candidates not linked to a signal truth photon are excluded, not relabelled",
        },
        "truth_signal": config["truth_signal"],
        "weighting_steps": steps,
        "producer_event_weight_used": False,
        "low_et_background_flattening": training.get("low_et_background_flattening"),
        "training_et_window_gev": [low, high],
        "split": {"rule": "row-level stratified", **(training.get("split") or {}),
                  "train": int(train.sum()), "validation": int(valid.sum()), "test": int(test.sum())},
        "samples": [{"role": "signal" if i < len(args.signal) else "background", "path": str(p), "sha256": sha256(p),
                     "rows_used": int(len(part["label"])), "candidates_total": part["n_total"],
                     "candidates_complete": part["n_complete"], "candidates_linked_to_signal_truth": part["n_linked"]}
                    for i, (p, part) in enumerate(zip([*args.signal, *args.background], parts))],
        "class_counts": {"signal": int((y == 1).sum()), "background": int((y == 0).sum())},
        "auc": auc,
    }
    (args.output_dir / f"{args.system}_model_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[train_photon_bdt] wrote {tmva_path} AUC={auc}")
    print(f"[train_photon_bdt] set systems.{args.system}.model_file and model_sha256={manifest['tmva_model']['sha256']} in the configuration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
