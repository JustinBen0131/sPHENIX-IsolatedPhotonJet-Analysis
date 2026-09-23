#!/usr/bin/env python3
"""One candidate-package training path. No publication or scientific acceptance.

Canonical weights are fitted on TRAIN and then frozen. The parity profile
intentionally retains historical all-population weighting and row splitting.
Labels and population witnesses are adapted once in inputs.py. Calibration
never receives test rows; the on-disk freeze precedes their first prediction.
"""

from pathlib import Path
import hashlib
import json
import subprocess
import sys
import numpy as np
from features import feature_identity
from registry import digest, file_hash, write_json
from calibration import fit_working_points

PARTITIONS = ("train", "validation", "test")
PARITY_PARAMS = dict(
    n_estimators=750,
    max_depth=5,
    learning_rate=0.1,
    subsample=0.5,
    colsample_bytree=0.6,
    colsample_bylevel=1.0,
    reg_alpha=5.0,
    reg_lambda=0.3,
    tree_method="hist",
    grow_policy="lossguide",
    max_bin=256,
    random_state=42,
    n_jobs=4,
)


def check_profile(spec, profile):
    """Refuse silently different recipes under the same public profile identity."""
    if (
        list(spec.features) != profile["features"]
        or spec.system != profile["system"]
        or spec.shower_definition != profile["shower_definition"]
    ):
        raise ValueError(
            "registry and training profile disagree on ordered features/system/shower"
        )
    for field in (
        "feature_binding_evidence",
        "label_binding_evidence",
        "preselection_binding_evidence",
        "domain_binding_evidence",
    ):
        if not profile.get(field):
            raise ValueError(
                f"{spec.name}: unresolved {field} in training_profiles.yaml"
            )
    if profile.get("label_adapter") != "dominant_prompt_v1":
        raise ValueError("unknown label adapter")
    s = profile["split"]
    if spec.name == "ppg12_equivalent_v1":
        if (
            s
            != dict(
                rule="row_stratified_70_10_20",
                seed=42,
                train_fraction=0.7,
                validation_fraction=0.1,
                test_fraction=0.2,
            )
            or profile["xgboost"] != PARITY_PARAMS
            or profile["weight_fit_population"] != "all_before_split"
            or profile["low_et_background_flattening"] is not True
        ):
            raise ValueError(
                "PPG12-equivalent mechanics differ from recovered executable recipe"
            )
    else:
        if (
            s
            != dict(
                rule="physical_event_hash_v1",
                seed=13,
                train_fraction=0.8,
                validation_fraction=0.1,
                test_fraction=0.1,
            )
            or profile["weight_fit_population"] != "train_only"
            or profile["low_et_background_flattening"] is not False
            or profile["et_window_gev"] != [15.0, 35.0]
        ):
            raise ValueError(
                "canonical protocol must remain grouped 80/10/10, train-only weights, 15–35 GeV"
            )
    if profile["weighting"] != dict(
        bins=20, eta_range=[-0.7, 0.7], et_cap=800.0, class_balance=True
    ):
        raise ValueError("weighting differs from recovered ET/eta/class-balance recipe")
    if profile["working_points"]["fit_form"] != "linear":
        raise ValueError("unsupported WP fit form")
    window = profile["et_window_gev"]
    if not window or window != [spec.domain["min_et_gev"], spec.domain["max_et_gev"]]:
        raise ValueError(
            "bind the exact training domain consistently in registry and profile"
        )


def parity_split_indices(labels):
    from sklearn.model_selection import train_test_split

    rows = np.arange(len(labels))
    rest, test = train_test_split(rows, test_size=0.2, random_state=42, stratify=labels)
    train, val = train_test_split(
        rest, test_size=0.1 / (0.7 + 0.1), random_state=43, stratify=labels[rest]
    )
    return train, val, test


def split_rows(groups, spec, labels=None):
    if spec["rule"] == "physical_event_hash_v1":
        draws = np.array(
            [
                int.from_bytes(
                    hashlib.sha256(
                        f"photonid-physical-v1:{spec['seed']}:{g}".encode()
                    ).digest()[:8],
                    "big",
                )
                / 2**64
                for g in groups
            ]
        )
        return np.where(
            draws < 0.8, "train", np.where(draws < 0.9, "validation", "test")
        )
    if spec["rule"] != "row_stratified_70_10_20" or labels is None:
        raise ValueError("unsupported split")
    train, val, test = parity_split_indices(labels)
    result = np.full(len(labels), "train", dtype="<U10")
    result[test] = "test"
    result[val] = "validation"
    return result


def parity_row_order(et, labels, samples, ready, source_order):
    """Historical per-source low-ET flattening, including its sampled row order.

    A logical training source can span several base files. Its manifest order
    must reproduce the historical concatenation; no filesystem sorting occurs.
    """
    import pandas as pd

    rows = []
    if set(np.asarray(samples)[ready]) != set(source_order):
        raise ValueError("training source roster is incomplete or unexpected")
    for source in source_order:
        ix = np.flatnonzero(ready & (samples == source))
        low = ix[(labels[ix] == 0) & (et[ix] < 15)]
        if len(low):
            bins = pd.cut(
                et[low],
                np.linspace(et[low].min(), 15.0, 21),
                labels=False,
                include_lowest=True,
            )
            counts = np.bincount(bins.astype(int), minlength=20)
            target = int(counts[counts > 0].min())
            rng = np.random.RandomState(42)
            low = np.concatenate(
                [
                    rng.choice(low[bins == b], target, replace=False)
                    for b in range(20)
                    if counts[b]
                ]
            )
            ix = np.concatenate([low, ix[et[ix] >= 15]])
        rows.extend(ix.tolist())
    return np.array(rows, dtype=int)


def fit_weights(et, eta, labels, fit_mask, config):
    """Save spline knots/coefficients and TRAIN normalizations for later reuse.

    Formula retained from PPG12 KinematicReweighter. The fit population is the
    only deliberate methodological difference for the two canonical profiles.
    """
    from scipy.interpolate import UnivariateSpline

    fitted = {}
    total = int(fit_mask.sum())
    for cls in (0, 1):
        mask = fit_mask & (labels == cls)
        if mask.sum() < 20:
            raise ValueError("weight fit requires at least 20 candidates per class")
        record = {"class_balance": total / (2 * int(mask.sum()))}
        for name, values, limits, cap in [
            ("eta", eta, config["eta_range"], None),
            ("et", et, None, config["et_cap"]),
        ]:
            x = values[mask]
            lo, hi = limits if limits else (x.min(), x.max())
            if not hi > lo or not np.all(np.isfinite(x)):
                raise ValueError("degenerate/nonfinite weight-fit range")
            edges = np.linspace(lo, hi, 21)
            density = np.histogram(x, bins=edges, density=True)[0] * 20
            spline = UnivariateSpline((edges[:-1] + edges[1:]) / 2, density, s=0.0)
            t, c, k = spline._eval_args
            w = 1 / np.maximum(spline(x), 1e-3)
            if cap is not None:
                w = np.minimum(w, cap)
            record[name] = dict(
                knots=t.tolist(),
                coefficients=c.tolist(),
                degree=k,
                normalization=float(w.mean()),
                cap=cap,
            )
        fitted[str(cls)] = record
    return fitted


def apply_weights(et, eta, labels, fitted):
    from scipy.interpolate import BSpline

    result = np.full(len(labels), np.nan)
    for cls in (0, 1):
        m = labels == cls
        r = fitted[str(cls)]
        result[m] = r["class_balance"]
        for name, values in [("eta", eta), ("et", et)]:
            d = r[name]
            pdf = BSpline(d["knots"], d["coefficients"], d["degree"])(values[m])
            w = 1 / np.maximum(pdf, 1e-3)
            if d["cap"] is not None:
                w = np.minimum(w, d["cap"])
            result[m] *= w / d["normalization"]
    return result


def predict(model, X, features):
    import xgboost as xgb

    if not len(X):
        return np.empty(0)
    return model.get_booster().predict(xgb.DMatrix(X, feature_names=list(features)))


def train_package(spec, profile, data, output, provenance, reference=None):
    """Internal per-profile executor, also exercised by small synthetic fixtures."""
    from xgboost import XGBClassifier
    from diagnostics import save_diagnostics

    output = Path(output)
    if output.exists():
        raise FileExistsError(f"refusing existing package: {output}")
    eligible = data["ready"] & (data["y"] >= 0)
    canonical = profile["split"]["rule"] == "physical_event_hash_v1"
    if canonical:
        order = np.flatnonzero(eligible)
        partitions = split_rows(data["groups"], profile["split"])
    else:
        order = parity_row_order(
            data["et"], data["y"], data["sample"], eligible, profile["source_order"]
        )
        partitions = np.full(len(eligible), "excluded", dtype="<U10")
        partitions[order] = split_rows(
            data["groups"][order], profile["split"], data["y"][order]
        )
    if len(np.unique(data["keys"], axis=0)) != len(data["keys"]):
        raise ValueError("duplicate training candidate key")
    selected = np.zeros(len(eligible), bool)
    selected[order] = True
    masks = {p: selected & (partitions == p) for p in PARTITIONS}
    for p, m in masks.items():
        if set(data["y"][m]) != {0, 1}:
            raise ValueError(f"{p} lacks one class")
    if canonical:
        sets = {p: set(data["groups"][partitions == p]) for p in PARTITIONS}
        if any(
            sets[a] & sets[b]
            for a, b in [
                ("train", "validation"),
                ("train", "test"),
                ("validation", "test"),
            ]
        ):
            raise ValueError("physical-event partition overlap")
    fit_mask = masks["train"] if canonical else selected
    fitted = fit_weights(
        data["et"], data["eta"], data["y"], fit_mask, profile["weighting"]
    )
    weights = np.full(len(selected), np.nan)
    before_freeze = selected & (partitions != "test") if canonical else selected
    weights[before_freeze] = apply_weights(
        data["et"][before_freeze],
        data["eta"][before_freeze],
        data["y"][before_freeze],
        fitted,
    )
    if not np.all(np.isfinite(weights[before_freeze]) & (weights[before_freeze] > 0)):
        raise ValueError("invalid training/validation weights")
    output.mkdir(parents=True)
    for folder in ("model", "split", "working_points", "isolation", "diagnostics"):
        (output / folder).mkdir()
    write_json(
        output / "INCOMPLETE.json", {"status": "building candidate; not accepted"}
    )
    write_json(output / "training_profile.json", profile)
    write_json(
        output / "split/weighting.json",
        dict(fit_population=profile["weight_fit_population"], classes=fitted),
    )
    np.savez_compressed(
        output / "split/assignments.npz",
        keys=data["keys"],
        group=data["groups"],
        partition=partitions,
        training_eligible=selected,
        training_row_order=order,
    )
    # Preserve historical row order for parity; canonical has no row shuffling contract.
    tr = (
        order[partitions[order] == "train"]
        if canonical
        else order[parity_split_indices(data["y"][order])[0]]
    )
    va = masks["validation"]
    model = XGBClassifier(**profile["xgboost"])
    import pandas as pd

    model.fit(
        pd.DataFrame(data["X"][tr], columns=list(spec.features)),
        data["y"][tr],
        sample_weight=weights[tr],
    )
    model_path = output / "model/model.json"
    model.save_model(model_path)
    reloaded = XGBClassifier()
    reloaded.load_model(model_path)
    scores = np.full(len(selected), np.nan)
    scores[va] = predict(model, data["X"][va], spec.features)
    again = predict(reloaded, data["X"][va], spec.features)
    if not np.all(np.isfinite(scores[va])) or not np.array_equal(scores[va], again):
        raise ValueError("nonfinite score or serialization/reload inequivalence")
    binding = dict(
        model_name=spec.name,
        model_sha256=file_hash(model_path),
        feature_schema_sha256=feature_identity(spec),
    )
    wp = fit_working_points(
        scores[va],
        weights[va],
        data["et"][va],
        data["centrality"][va],
        data["y"][va] == 1,
        spec.system,
        partition="validation",
    )
    wp.update(binding)
    write_json(output / "working_points/id.json", wp)
    iso = None
    if spec.system == "auau" and profile.get("analysis_signal_binding_evidence"):
        iv = (partitions == "validation") & data["analysis"]
        iso = fit_working_points(
            data["isolation"][iv],
            np.ones(int(iv.sum())),
            data["et"][iv],
            data["centrality"][iv],
            np.ones(int(iv.sum()), bool),
            "auau",
            partition="validation",
            isolation=True,
        )
        iso.update(binding)
        write_json(output / "isolation/calibration.json", iso)
        # Calibration is useful before the nominal alias/sideband decision is bound.
        settings = profile["isolation"]
        if (
            settings.get("nominal_target")
            and settings.get("nonisolated_rule")
            and settings.get("binding_evidence")
        ):
            if settings["nonisolated_rule"] != "same_as_isolated":
                raise ValueError("unsupported bound sideband rule")
            curve = iso["curves"][settings["nominal_target"]]
            write_json(
                output / "isolation/selection.json",
                dict(
                    system="auau",
                    method=2,
                    radius=0.4,
                    axis="centrality",
                    domain=[0, 80],
                    isolated=curve,
                    nonisolated=curve,
                    calibration_sha256=file_hash(output / "isolation/calibration.json"),
                    status="candidate",
                    binding_evidence=settings["binding_evidence"],
                ),
            )
    elif spec.system == "auau":
        write_json(
            output / "isolation/UNBOUND.json",
            dict(
                status="unbound",
                reason="analysis_signal_binding_evidence is required for isolation calibration; model/ID remain independent",
            ),
        )
    freeze = {
        str(p.relative_to(output)): file_hash(p)
        for p in output.rglob("*")
        if p.is_file()
    }
    write_json(
        output / "FROZEN_CHOICES.json",
        dict(status="candidate choices frozen before test prediction", files=freeze),
    )
    # Only now may test rows be evaluated. Never call a fit function below here.
    if canonical:
        te = masks["test"]
        weights[te] = apply_weights(
            data["et"][te], data["eta"][te], data["y"][te], fitted
        )
    if not np.all(np.isfinite(weights[selected]) & (weights[selected] > 0)):
        raise ValueError("invalid frozen evaluation weights")
    rest = selected & ~va
    scores[rest] = predict(reloaded, data["X"][rest], spec.features)
    if not np.all(np.isfinite(scores[selected])):
        raise ValueError("nonfinite post-freeze scores")
    metrics = {}
    for population in ("validation", "test", "all"):
        mask = selected if population == "all" else masks[population]
        im = (
            data["analysis"]
            if population == "all"
            else data["analysis"] & (partitions == population)
        )
        payload = dict(
            score=scores[mask],
            y=data["y"][mask],
            weight=weights[mask],
            et=data["et"][mask],
            centrality=data["centrality"][mask],
            keys=data["keys"][mask],
            isolation=data["isolation"][im],
            isolation_et=data["et"][im],
            isolation_centrality=data["centrality"][im],
        )
        if reference is not None:
            payload["reference_score"] = reference(data["X"][mask])
        metrics[population] = save_diagnostics(
            output / "diagnostics" / population, payload, wp, iso, population, spec.name
        )
    for name, sha in freeze.items():
        if file_hash(output / name) != sha:
            raise ValueError("frozen choices changed during evaluation")
    for source in provenance.get("checked_inputs", []):
        if file_hash(source["path"]) != source["sha256"]:
            raise ValueError("training input changed during package construction")
    hashes = {
        str(p.relative_to(output)): file_hash(p)
        for p in output.rglob("*")
        if p.is_file() and p.name != "INCOMPLETE.json"
    }
    receipt = dict(
        schema="PhotonIDTrainingReceiptV1",
        status="candidate",
        mechanical_status="PASS",
        scientific_acceptance="NOT_REVIEWED",
        **binding,
        model_version=provenance.get("model_version", "1"),
        system=spec.system,
        features=list(spec.features),
        training_profile_sha256=digest(profile),
        provenance=provenance,
        split=profile["split"],
        split_counts={p: int(m.sum()) for p, m in masks.items()},
        split_group_counts={p: len(set(data["groups"][m])) for p, m in masks.items()},
        weights=profile["weight_fit_population"],
        hyperparameters=profile["xgboost"],
        metrics=metrics,
        files=hashes,
        reference_comparison="evaluated" if reference else "unbound; not evaluated",
        isolation_status=(
            "calibrated"
            if iso
            else (
                "unbound denominator; calibration/plots omitted"
                if spec.system == "auau"
                else "fixed pp package in registry"
            )
        ),
        id_package_sha256=hashes["working_points/id.json"],
        isolation_package_sha256=hashes.get("isolation/selection.json"),
        isolation_calibration_sha256=hashes.get("isolation/calibration.json"),
    )
    import importlib.metadata

    receipt["software_versions"] = {
        n: importlib.metadata.version(n)
        for n in (
            "numpy",
            "scipy",
            "pandas",
            "scikit-learn",
            "xgboost",
            "uproot",
            "matplotlib",
        )
    }
    write_json(output / "TRAINING_RECEIPT.json", receipt)
    (output / "INCOMPLETE.json").unlink()
    return receipt


if __name__ == "__main__":
    from frontend import main

    raise SystemExit(main("train", sys.argv[1:]))
