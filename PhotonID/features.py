"""The one feature-extraction implementation for photon identification.

Every model in ``model_registry.yaml`` takes its inputs from here, whether it
is being trained (``train.py``) or applied (``augment.py``). There is no
second copy of any feature definition, so a model can never be trained on one
definition and applied on another.

Inputs are the canonical trees written by ``TreeProduction`` (tables
``Events``, ``Photons``, ``PhotonShowerViews``); the base tree stores
primitives only. Each candidate is identified by ``(source, event, photon)``
128-bit identities, never by row order.

Feature vocabulary (names as the registries spell them)

    cluster_Et            Photons.et
    cluster_Eta           Photons.eta
    cluster_Phi           Photons.phi
    vertexz               Photons.producer_vertex_z
    centrality            Events.centrality_percent (NaN when not valid)
    cluster_<x>           PhotonShowerViews.<x> of the requested definition,
                          x in e11..e77, weta, wphi, weta_cog, wphi_cog,
                          weta_cogx, wphi_cogx, weta33_cogx, wphi33_cogx,
                          w32, w52, w72, et1..et4
    <a>_over_<b>          ratio of two rectangular sums of the same view;
                          the denominator rule is the registry's ratio_policy

A feature that cannot be formed is NaN. The caller decides what a NaN means
(``augment.py`` records it as missing_inputs); this module never substitutes
a plausible value.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import uproot
import yaml


EVENT_BRANCHES = ("source_hi", "source_lo", "event_hi", "event_lo",
                  "centrality_percent", "centrality_valid", "reco_vertex_z", "reco_vertex_valid")
PHOTON_BRANCHES = ("event_hi", "event_lo", "photon_hi", "photon_lo", "encounter_ordinal",
                   "energy", "et", "eta", "phi", "kinematics_finite", "producer_vertex_z")
SHOWER_SCALARS = ("e11", "e13", "e15", "e17", "e22", "e31", "e32", "e33", "e35", "e37",
                  "e51", "e52", "e53", "e55", "e57", "e71", "e72", "e73", "e75", "e77",
                  "weta", "wphi", "weta_cog", "wphi_cog", "weta_cogx", "wphi_cogx",
                  "weta33_cogx", "wphi33_cogx", "w32", "w52", "w72",
                  "et1", "et2", "et3", "et4", "valid")

RATIO_POLICIES = ("builder_zero", "nan")

CandidateKey = tuple[int, int, int, int]  # event_hi, event_lo, photon_hi, photon_lo


@dataclass(frozen=True)
class ModelSpec:
    """One registry entry, validated."""

    name: str
    system: str
    shower_definition: str
    features: tuple[str, ...]
    ratio_policy: str
    domain: Mapping[str, Any]
    format: str
    model_file: Path | None
    model_sha256: str | None
    cdb_key: str | None
    score_direction: str
    training_manifest: Path | None
    training: Mapping[str, Any] = field(default_factory=dict)

    @property
    def bound(self) -> bool:
        return self.model_file is not None and self.model_sha256 is not None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def feature_identity(spec: ModelSpec) -> str:
    """Bind feature order, arithmetic implementation, validity and domain."""
    definition = {"features": spec.features, "shower_definition": spec.shower_definition,
                  "ratio_policy": spec.ratio_policy, "domain": dict(spec.domain),
                  "numeric_type": "float32", "implementation_sha256": sha256(Path(__file__))}
    return hashlib.sha256(json.dumps(definition, sort_keys=True, allow_nan=False).encode()).hexdigest()


def base_metadata(root) -> dict[str, str]:
    """Read the producer's key=value metadata; incomplete sources are unusable."""
    if "metadata" not in root or "completion" not in root:
        raise ValueError("base product lacks metadata or completion")
    def parse(name):
        return dict(line.split("=", 1) for line in str(root[name]).splitlines() if "=" in line)
    metadata, completion = parse("metadata"), parse("completion")
    if metadata.get("contract") != "PhotonJetTrees" or completion.get("completion_status") != "complete":
        raise ValueError("base product is not a completed PhotonJetTrees file")
    if root["Events"].num_entries != int(completion["retained_events"]):
        raise ValueError("event rows disagree with completion accounting")
    if root["Sources"].num_entries != 1 or not bool(root["Sources"]["completed"].array(library="np")[0]):
        raise ValueError("base source is incomplete or crosses unsupported source boundaries")
    expected = int(root["Events"]["photon_count"].array(library="np").sum())
    if root["Photons"].num_entries != expected:
        raise ValueError("photon rows disagree with event accounting")
    return metadata


def load_registry(path: Path | str) -> dict[str, ModelSpec]:
    path = Path(path)
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or "models" not in raw:
        raise ValueError(f"{path}: not a model registry")
    models: dict[str, ModelSpec] = {}
    for name, spec in raw["models"].items():
        if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name):
            raise ValueError("unsafe model identifier")
        if spec.get("model_file") is not None and (spec.get("status") != "approved" or not spec.get("binding_evidence")):
            raise ValueError(f"{name}: model binding requires approval evidence")
        policy = spec.get("ratio_policy", "nan")
        if policy not in RATIO_POLICIES:
            raise ValueError(f"{name}: ratio_policy must be one of {RATIO_POLICIES}")
        features = tuple(str(f) for f in spec["features"])
        if not features or len(set(features)) != len(features):
            raise ValueError(f"{name}: duplicate feature names")
        for feature in features:
            parse_feature(feature)  # validates the vocabulary
        if spec.get("system") not in {"pp", "auau", "any"}:
            raise ValueError(f"{name}: invalid system")
        if spec.get("format") not in {"tmva_rbdt", "xgboost_json"}:
            raise ValueError(f"{name}: unsupported format")
        if spec.get("model_sha256") is not None and not re.fullmatch(r"[0-9a-f]{64}", spec["model_sha256"]):
            raise ValueError(f"{name}: malformed model hash")
        if (spec.get("model_file") is None) != (spec.get("model_sha256") is None):
            raise ValueError(f"{name}: model file and hash must be bound together")
        domain = spec.get("domain") or {}
        low, high, eta = (domain.get(k) for k in ("min_et_gev", "max_et_gev", "max_abs_eta"))
        if not all(isinstance(x, (int, float)) and np.isfinite(x) for x in (low, high, eta)) or not (low < high and eta > 0):
            raise ValueError(f"{name}: unresolved model domain")
        model_file = spec.get("model_file")
        manifest = spec.get("training_manifest")
        models[name] = ModelSpec(
            name=name,
            system=str(spec.get("system", "any")),
            shower_definition=str(spec.get("shower_definition", "H70")),
            features=features,
            ratio_policy=policy,
            domain=dict(spec.get("domain") or {}),
            format=str(spec.get("format", "tmva_rbdt")),
            model_file=None if model_file is None else (path.parent / model_file),
            model_sha256=spec.get("model_sha256"),
            cdb_key=spec.get("cdb_key"),
            score_direction=str(spec.get("score_direction", "higher_is_signal")),
            training_manifest=None if manifest is None else (path.parent / manifest),
            training=dict(spec.get("training") or {}),
        )
    return models


# ----------------------------------------------------------------------------
# Feature vocabulary
# ----------------------------------------------------------------------------

def parse_feature(name: str) -> tuple[str, tuple[str, ...]]:
    """Return (kind, arguments) for one feature name; raise for unknown names."""

    if name in ("cluster_Et", "cluster_Eta", "cluster_Phi", "vertexz", "vertex_z", "centrality"):
        return "scalar", (name,)
    if "_over_" in name:
        numerator, denominator = name.split("_over_", 1)
        for part in (numerator, denominator):
            if part not in SHOWER_SCALARS:
                raise ValueError(f"unknown ratio operand {part!r} in feature {name!r}")
        return "ratio", (numerator, denominator)
    if name.startswith("cluster_"):
        column = name[len("cluster_"):]
        if column in SHOWER_SCALARS:
            return "shower", (column,)
    if name in SHOWER_SCALARS:
        return "shower", (name,)
    raise ValueError(f"unknown feature name {name!r}")


# ----------------------------------------------------------------------------
# Reading the canonical tables
# ----------------------------------------------------------------------------

@dataclass
class CandidateTable:
    """Flat per-candidate columns joined from Events, Photons and one shower view."""

    keys: list[CandidateKey]
    source: np.ndarray            # (n, 2) source identity words
    columns: dict[str, np.ndarray]
    shower_definition: str

    def __len__(self) -> int:
        return len(self.keys)


def _identity_pairs(arrays: Mapping[str, np.ndarray], prefix: str) -> list[tuple[int, int]]:
    return list(zip(np.asarray(arrays[f"{prefix}_hi"]).astype(np.uint64).tolist(),
                    np.asarray(arrays[f"{prefix}_lo"]).astype(np.uint64).tolist()))


def read_candidates(path: Path | str, shower_definition: str = "H70") -> CandidateTable:
    """Join Photons with their event context and one shower-shape definition.

    Every photon row must have exactly one shower view of the requested
    definition and exactly one event row; anything else is an error, not a
    default.
    """

    path = Path(path)
    with uproot.open(path) as root:
        base_metadata(root)
        for required in ("Events", "Photons", "PhotonShowerViews"):
            if required not in root:
                raise ValueError(f"{path}: missing table {required}")
        events = root["Events"].arrays(list(EVENT_BRANCHES), library="np")
        photons = root["Photons"].arrays(list(PHOTON_BRANCHES), library="np")
        views = root["PhotonShowerViews"].arrays(
            ["event_hi", "event_lo", "photon_hi", "photon_lo", "definition", *SHOWER_SCALARS], library="np")

    event_index = {key: i for i, key in enumerate(_identity_pairs(events, "event"))}
    if len(event_index) != len(events["event_hi"]):
        raise ValueError(f"{path}: event identities are not unique")

    photon_keys: list[CandidateKey] = []
    for (eh, el), (ph, pl) in zip(_identity_pairs(photons, "event"), _identity_pairs(photons, "photon")):
        photon_keys.append((eh, el, ph, pl))
    if len(set(photon_keys)) != len(photon_keys):
        raise ValueError(f"{path}: photon identities are not unique")
    photon_index = {key: i for i, key in enumerate(photon_keys)}

    n = len(photon_keys)
    columns: dict[str, np.ndarray] = {}
    for name in ("energy", "et", "eta", "phi", "producer_vertex_z"):
        columns[name] = np.asarray(photons[name], dtype=float)
    columns["encounter_ordinal"] = np.asarray(photons["encounter_ordinal"], dtype=np.int64)
    columns["kinematics_finite"] = np.asarray(photons["kinematics_finite"], dtype=bool)

    # event context
    rows = np.empty(n, dtype=np.int64)
    for i, key in enumerate(photon_keys):
        try:
            rows[i] = event_index[(key[0], key[1])]
        except KeyError as error:
            raise ValueError(f"{path}: photon row {i} references an unknown event") from error
    centrality = np.asarray(events["centrality_percent"], dtype=float)[rows]
    centrality_valid = np.asarray(events["centrality_valid"], dtype=bool)[rows]
    columns["centrality"] = np.where(centrality_valid, centrality, np.nan)
    columns["reco_vertex_z"] = np.asarray(events["reco_vertex_z"], dtype=float)[rows]
    columns["reco_vertex_valid"] = np.asarray(events["reco_vertex_valid"], dtype=bool)[rows]
    source = np.stack([np.asarray(events["source_hi"]).astype(np.uint64)[rows],
                       np.asarray(events["source_lo"]).astype(np.uint64)[rows]], axis=1)

    # shower view of the requested definition
    definitions = np.asarray(views["definition"]).astype(str)
    selected = np.flatnonzero(definitions == shower_definition)
    if n and len(selected) == 0:
        raise ValueError(f"{path}: no shower view named {shower_definition!r}")
    view_rows = np.full(n, -1, dtype=np.int64)
    view_keys = list(zip(_identity_pairs(views, "event"), _identity_pairs(views, "photon")))
    for j in selected:
        (eh, el), (ph, pl) = view_keys[j]
        i = photon_index.get((eh, el, ph, pl))
        if i is None:
            raise ValueError(f"{path}: shower view row {j} references an unknown photon")
        if view_rows[i] >= 0:
            raise ValueError(f"{path}: photon row {i} has two {shower_definition!r} shower views")
        view_rows[i] = j
    if np.any(view_rows < 0):
        raise ValueError(f"{path}: {int((view_rows < 0).sum())} photons have no {shower_definition!r} shower view")
    for name in SHOWER_SCALARS:
        values = np.asarray(views[name])
        columns[f"shower_{name}"] = values[view_rows].astype(float if name != "valid" else bool)

    return CandidateTable(keys=photon_keys, source=source, columns=columns, shower_definition=shower_definition)


# ----------------------------------------------------------------------------
# Feature evaluation
# ----------------------------------------------------------------------------

def _ratio(numerator: np.ndarray, denominator: np.ndarray, policy: str) -> np.ndarray:
    out = np.full(numerator.shape, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        positive = np.isfinite(denominator) & (denominator > 0.0) & np.isfinite(numerator)
        out[positive] = numerator[positive] / denominator[positive]
    if policy == "builder_zero":
        # The reconstruction-resident evaluator returned 0 for a non-positive
        # denominator; historical models were trained and applied that way.
        nonpositive = np.isfinite(denominator) & ~(denominator > 0.0)
        out[nonpositive] = 0.0
    return out


def feature_column(table: CandidateTable, name: str, ratio_policy: str) -> np.ndarray:
    kind, args = parse_feature(name)
    if kind == "scalar":
        key = args[0]
        if key == "cluster_Et":
            return table.columns["et"]
        if key == "cluster_Eta":
            return table.columns["eta"]
        if key == "cluster_Phi":
            return table.columns["phi"]
        if key in ("vertexz", "vertex_z"):
            return table.columns["producer_vertex_z"]
        if key == "centrality":
            return table.columns["centrality"]
        raise AssertionError(key)
    if kind == "shower":
        return table.columns[f"shower_{args[0]}"]
    if kind == "ratio":
        return _ratio(table.columns[f"shower_{args[0]}"], table.columns[f"shower_{args[1]}"], ratio_policy)
    raise AssertionError(kind)


@dataclass
class FeatureMatrix:
    keys: list[CandidateKey]
    source: np.ndarray
    names: tuple[str, ...]
    values: np.ndarray          # (n, k) float32, shared training/inference precision
    complete: np.ndarray        # (n,) every input finite
    in_domain: np.ndarray       # (n,) inside the model's kinematic domain
    shower_valid: np.ndarray    # (n,) the shower view reported itself valid
    et: np.ndarray
    eta: np.ndarray
    centrality: np.ndarray
    encounter_ordinal: np.ndarray


def in_domain_mask(table: CandidateTable, domain: Mapping[str, Any]) -> np.ndarray:
    et = table.columns["et"]
    eta = table.columns["eta"]
    mask = np.isfinite(et) & np.isfinite(eta)
    if "min_et_gev" in domain:
        mask &= et >= float(domain["min_et_gev"])
    if "max_et_gev" in domain:
        mask &= et < float(domain["max_et_gev"])
    if "max_abs_eta" in domain:
        mask &= np.abs(eta) < float(domain["max_abs_eta"])
    if "centrality_percent" in domain and domain["centrality_percent"] is not None:
        low, high = (float(v) for v in domain["centrality_percent"])
        c = table.columns["centrality"]
        mask &= np.isfinite(c) & (c >= low) & (c < high)
    return mask


def build_features(table: CandidateTable, spec: ModelSpec) -> FeatureMatrix:
    """Evaluate the ordered feature vector of one model for every candidate."""

    if table.shower_definition != spec.shower_definition:
        raise ValueError(f"table carries {table.shower_definition!r} but model {spec.name} needs {spec.shower_definition!r}")
    columns = [feature_column(table, name, spec.ratio_policy) for name in spec.features]
    values = np.stack(columns, axis=1) if columns else np.zeros((len(table), 0))
    # TMVA and XGBoost consume float32. Check the actual evaluation precision,
    # including overflow introduced by narrowing, in both train and augment.
    with np.errstate(over="ignore", invalid="ignore"):
        values = values.astype(np.float32)
    return FeatureMatrix(
        keys=table.keys, source=table.source, names=spec.features, values=values,
        complete=(np.all(np.isfinite(values), axis=1) & table.columns["shower_valid"]
                  & table.columns["kinematics_finite"]),
        in_domain=in_domain_mask(table, spec.domain),
        shower_valid=table.columns["shower_valid"].astype(bool),
        et=table.columns["et"], eta=table.columns["eta"], centrality=table.columns["centrality"],
        encounter_ordinal=table.columns["encounter_ordinal"],
    )


def features_for_file(path: Path | str, spec: ModelSpec) -> FeatureMatrix:
    with uproot.open(path) as root:
        metadata = base_metadata(root)
    system = {"1": "pp", "2": "auau"}.get(metadata.get("collision_system"))
    if spec.system != "any" and system != spec.system:
        raise ValueError(f"{path}: collision system differs from model {spec.name}")
    return build_features(read_candidates(path, spec.shower_definition), spec)


__all__ = [
    "CandidateTable", "FeatureMatrix", "ModelSpec", "RATIO_POLICIES", "SHOWER_SCALARS",
    "build_features", "feature_column", "features_for_file", "in_domain_mask",
    "load_registry", "parse_feature", "read_candidates",
]
