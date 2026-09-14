"""Compile plot annotations from the exact executable selection lineage.

No public function accepts free-form physics labels.  Structured dataset
metadata and the receipt-bound selection program are the only semantic inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import uproot
import yaml

from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.selection import SelectionProgram, compile_recoil_selection
from photonjet.provenance import artifact, canonical_json, sha256_file, write_json


ROOT_EXPERIMENT_LABEL = "#it{#bf{sPHENIX}} Internal"
# Nested mathit/mathbf selects the inner font rather than combining styles.
# mathbfit requests the actual bold-oblique glyphs; Internal stays upright.
MATPLOTLIB_EXPERIMENT_LABEL = r"$\mathbfit{sPHENIX}$ Internal"
STYLE_ID = "sphenix_internal_provenance_compiled_v2"
SEMANTIC_COVERAGE_SCHEMA = "PhotonJetSemanticCoverageV1"

_PREDICATE_KEYS = {
    "event.accepted",
    "photon.et",
    "photon.abs_eta",
    "jet.pt",
    "jet.abs_eta",
    "jet.radius",
    "recoil.delta_phi",
}
_CLASSIFICATION_PREDICATE_KEYS = {
    "photon.bdt_class",
    "photon.isolation_class",
}
_FACT_KEYS = {
    "photon.candidate_scope",
    "photon.abcd_region",
    "photon.id_state",
    "photon.isolation_state",
    "photon.isolation_radius",
    "photon.non_tight_definition",
    "photon.score_source",
    "photon.bdt_threshold_source",
    "photon.isolation_threshold_source",
    "photon.leader_rule",
}
_DATASET_FIELDS = {
    "schema",
    "collision_system",
    "sample_kind",
    "sqrt_s_gev",
    "sqrt_s_nn_gev",
    "centrality_percent",
    "event_count",
    "experiment_status",
}


def _number(value: float) -> str:
    rounded = round(float(value))
    if math.isclose(float(value), rounded, rel_tol=0.0, abs_tol=1.0e-12):
        return str(rounded)
    return f"{float(value):.6g}"


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def _load_mapping(path: Path) -> dict[str, Any]:
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if source.suffix.lower() in {".yaml", ".yml"}:
        value = yaml.safe_load(source.read_text(encoding="utf-8"))
    else:
        value = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected an object in {source}")
    return value


def _selection_from_payload(payload: dict[str, Any]) -> RecoilSelection:
    if payload.get("schema") != "PhotonJetHistogramPayloadV1":
        raise ValueError("plot input must be PhotonJetHistogramPayloadV1")
    if payload.get("observable") != "xjgamma":
        raise ValueError("V1 renderer only accepts the xjgamma observable")
    selection = payload.get("selection")
    if not isinstance(selection, dict):
        raise ValueError("histogram payload lacks a structured selection")
    expected = set(RecoilSelection.__dataclass_fields__)
    if set(selection) != expected:
        missing = sorted(expected - set(selection))
        extra = sorted(set(selection) - expected)
        raise ValueError(f"selection fields differ from RecoilSelection; missing={missing}, extra={extra}")
    result = RecoilSelection(**selection)
    if payload.get("selection_program_sha256") != compile_recoil_selection(result).sha256:
        raise ValueError("histogram payload lacks its current executable selection identity")
    return result


def write_histogram_receipt(
    *,
    input_paths: Iterable[Path],
    histogram_path: Path,
    selection: RecoilSelection,
    dataset_manifest_path: Path,
    receipt_path: Path,
) -> dict[str, Any]:
    """Bind a histogram byte stream to the program that selected its rows."""

    histogram = Path(histogram_path).resolve()
    dataset_file = Path(dataset_manifest_path).resolve()
    dataset = DatasetDescriptor.from_mapping(_load_mapping(dataset_file))
    payload = _load_mapping(histogram)
    if _selection_from_payload(payload) != selection:
        raise ValueError("histogram selection differs from the requested receipt selection")
    program = compile_recoil_selection(selection)
    paths = [Path(path).resolve() for path in input_paths]
    total_events = 0
    accepted_events = 0
    observed_centrality: list[float] = []
    for path in paths:
        with uproot.open(path) as root:
            if "events" not in root:
                raise ValueError(f"dataset binding requires PhotonJetTrees_v1/events: {path}")
            total_events += int(root["events"].num_entries)
            accepted_events += sum(
                program.accepts_event(status)
                for status in root["events"]["terminal_status"].array(library="np")
            )
            observed_centrality.extend(
                np.asarray(root["events"]["centrality"].array(library="np"), dtype=float).tolist()
            )
    centrality_values = np.asarray(observed_centrality, dtype=float)
    pp_not_applicable = bool(np.all(centrality_values == -1.0))
    auau_percent = bool(
        np.all(np.isfinite(centrality_values))
        and np.all(centrality_values >= 0.0)
        and np.all(centrality_values < 100.0)
    )
    if dataset.collision_system == "pp" and not pp_not_applicable:
        raise ValueError("p+p dataset manifest requires the canonical centrality=-1 state")
    if dataset.collision_system == "auau" and not auau_percent:
        raise ValueError("Au+Au dataset manifest requires centrality within [0, 100)")
    if dataset.event_count is not None and dataset.event_count != total_events:
        raise ValueError(
            f"dataset event_count={dataset.event_count} differs from tree total {total_events}"
        )
    if dataset.centrality_percent is not None:
        low, high = dataset.centrality_percent
        if np.any(centrality_values < low) or np.any(centrality_values >= high):
            raise ValueError("tree events fall outside the declared centrality interval")
    receipt = {
        "schema": "PhotonJetHistogramReceiptV1",
        "histogram": artifact(histogram, "histogram_payload"),
        "inputs": [artifact(path, "PhotonJetTrees_v1") for path in paths],
        "selection_program": program.to_dict(),
        "selection_program_sha256": program.sha256,
        "dataset_manifest_sha256": sha256_file(dataset_file),
        "dataset": dataset.to_dict(),
        "dataset_observation": {
            "event_count": total_events,
            "accepted_events": accepted_events,
            "retained_terminal_events_excluded": total_events - accepted_events,
            "centrality_state": (
                "not_applicable_minus_one" if pp_not_applicable else "auau_percent"
            ),
        },
    }
    write_json(receipt_path, receipt)
    return receipt


@dataclass(frozen=True)
class DatasetDescriptor:
    collision_system: str
    sample_kind: str
    experiment_status: str
    sqrt_s_gev: float | None
    sqrt_s_nn_gev: float | None
    centrality_percent: tuple[float, float] | None
    event_count: int | None

    @classmethod
    def from_mapping(cls, value: dict[str, Any]) -> "DatasetDescriptor":
        if value.get("schema") != "PhotonJetDatasetManifestV1":
            raise ValueError("dataset schema must be PhotonJetDatasetManifestV1")
        extra = sorted(set(value) - _DATASET_FIELDS)
        if extra:
            raise ValueError(f"free-form or unknown dataset fields are forbidden: {extra}")
        system = str(value.get("collision_system"))
        if system not in {"pp", "auau"}:
            raise ValueError("collision_system must be pp or auau")
        kind = str(value.get("sample_kind"))
        if kind not in {"data", "simulation", "embedding", "engineering_fixture"}:
            raise ValueError("unsupported sample_kind")
        if "experiment_status" not in value:
            raise ValueError("experiment_status is required")
        status = str(value["experiment_status"])
        if status != "Internal":
            raise ValueError("V1 public renderer is fail-closed to the sPHENIX Internal status")
        sqrt_s = value.get("sqrt_s_gev")
        sqrt_s_nn = value.get("sqrt_s_nn_gev")
        if system == "pp":
            if sqrt_s is None or sqrt_s_nn is not None:
                raise ValueError("pp requires sqrt_s_gev and forbids sqrt_s_nn_gev")
        else:
            if sqrt_s_nn is None or sqrt_s is not None:
                raise ValueError("auau requires sqrt_s_nn_gev and forbids sqrt_s_gev")
        centrality = value.get("centrality_percent")
        parsed_centrality: tuple[float, float] | None = None
        if centrality is not None:
            if system != "auau" or not isinstance(centrality, list) or len(centrality) != 2:
                raise ValueError("centrality_percent is an Au+Au [low, high] interval")
            parsed_centrality = (float(centrality[0]), float(centrality[1]))
            if not 0 <= parsed_centrality[0] < parsed_centrality[1] <= 100:
                raise ValueError("centrality interval is outside [0, 100]")
        event_count = value.get("event_count")
        if event_count is not None and (not isinstance(event_count, int) or event_count <= 0):
            raise ValueError("event_count must be a positive integer")
        return cls(
            collision_system=system,
            sample_kind=kind,
            experiment_status=status,
            sqrt_s_gev=None if sqrt_s is None else float(sqrt_s),
            sqrt_s_nn_gev=None if sqrt_s_nn is None else float(sqrt_s_nn),
            centrality_percent=parsed_centrality,
            event_count=event_count,
        )

    def line(self, dialect: str) -> str:
        if dialect not in {"matplotlib", "root"}:
            raise ValueError(dialect)
        if self.collision_system == "pp":
            system = "p+p"
            energy = (
                rf"$\sqrt{{s}}={_number(float(self.sqrt_s_gev))}\ \mathrm{{GeV}}$"
                if dialect == "matplotlib"
                else f"#sqrt{{s}} = {_number(float(self.sqrt_s_gev))} GeV"
            )
        else:
            system = "Au+Au"
            energy = (
                rf"$\sqrt{{s_{{NN}}}}={_number(float(self.sqrt_s_nn_gev))}\ \mathrm{{GeV}}$"
                if dialect == "matplotlib"
                else f"#sqrt{{s_{{NN}}}} = {_number(float(self.sqrt_s_nn_gev))} GeV"
            )
        kind = {
            "data": "Data",
            "simulation": "Simulation",
            "embedding": "Embedded simulation",
            "engineering_fixture": "Engineering fixture",
        }[self.sample_kind]
        parts = [system, energy, kind]
        if self.centrality_percent is not None:
            low, high = self.centrality_percent
            parts.append(f"{_number(low)}-{_number(high)}% centrality")
        if self.event_count is not None and self.sample_kind == "engineering_fixture":
            parts.append(f"{self.event_count} input events")
        return ", ".join(parts)

    def to_dict(self) -> dict[str, Any]:
        return {
            "collision_system": self.collision_system,
            "sample_kind": self.sample_kind,
            "experiment_status": self.experiment_status,
            "sqrt_s_gev": self.sqrt_s_gev,
            "sqrt_s_nn_gev": self.sqrt_s_nn_gev,
            "centrality_percent": None if self.centrality_percent is None else list(self.centrality_percent),
            "event_count": self.event_count,
        }

    @classmethod
    def from_contract_mapping(cls, value: dict[str, Any]) -> "DatasetDescriptor":
        expected = _DATASET_FIELDS - {"schema"}
        if set(value) != expected:
            missing = sorted(expected - set(value))
            extra = sorted(set(value) - expected)
            raise ValueError(f"plot-contract dataset fields differ; missing={missing}, extra={extra}")
        return cls.from_mapping({"schema": "PhotonJetDatasetManifestV1", **value})


def _program_maps(program: SelectionProgram) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    predicates = {item["key"]: item for item in program.to_dict()["predicates"]}
    facts = {item["key"]: item["value"] for item in program.to_dict()["facts"]}
    scope = facts.get("photon.candidate_scope")
    expected_predicates = (
        _PREDICATE_KEYS
        if scope == "inclusive_pairs"
        else _PREDICATE_KEYS | _CLASSIFICATION_PREDICATE_KEYS
    )
    if set(predicates) != expected_predicates:
        raise ValueError(f"unsupported predicate inventory: {sorted(predicates)}")
    unknown_facts = sorted(set(facts) - _FACT_KEYS)
    if unknown_facts:
        raise ValueError(f"unsupported selection facts: {unknown_facts}")
    return predicates, facts


def _delta_phi(value: float, dialect: str) -> str:
    ratio = Fraction(value / math.pi).limit_denominator(16)
    if not math.isclose(float(ratio) * math.pi, value, rel_tol=0.0, abs_tol=1.0e-12):
        return _number(value)
    numerator = "" if ratio.numerator == 1 else str(ratio.numerator)
    fraction = f"{numerator}\\pi/{ratio.denominator}" if dialect == "matplotlib" else f"{numerator}#pi/{ratio.denominator}"
    return fraction


def _annotation_lines(program: SelectionProgram, dialect: str) -> list[str]:
    predicates, facts = _program_maps(program)
    photon_et = predicates["photon.et"]
    photon_eta = predicates["photon.abs_eta"]
    jet_pt = predicates["jet.pt"]
    jet_eta = predicates["jet.abs_eta"]
    radius = predicates["jet.radius"]
    delta_phi = predicates["recoil.delta_phi"]

    if dialect == "matplotlib":
        photon_line = (
            rf"${_number(photon_et['lower'])}\leq p_T^\gamma<{_number(photon_et['upper'])}\ "
            rf"\mathrm{{GeV}},\ |\eta^\gamma|<{_number(photon_eta['value'])}$"
        )
        jet_line = (
            rf"$R={_number(radius['value'])}\ \mathrm{{jets}},\ p_T^{{\mathrm{{jet}}}}>{_number(jet_pt['value'])}\ "
            rf"\mathrm{{GeV}},\ |\eta^{{\mathrm{{jet}}}}|<{_number(jet_eta['value'])},\ "
            rf"\Delta\phi>{_delta_phi(float(delta_phi['value']), dialect)}$"
        )
    elif dialect == "root":
        photon_line = (
            f"{_number(photon_et['lower'])} #leq p_{{T}}^{{#gamma}} < {_number(photon_et['upper'])} GeV, "
            f"|#eta^{{#gamma}}| < {_number(photon_eta['value'])}"
        )
        jet_line = (
            f"R = {_number(radius['value'])} jets, p_{{T}}^{{jet}} > {_number(jet_pt['value'])} GeV, "
            f"|#eta^{{jet}}| < {_number(jet_eta['value'])}, "
            f"#Delta#phi > {_delta_phi(float(delta_phi['value']), dialect)}"
        )
    else:
        raise ValueError(dialect)

    lines = [photon_line]
    scope = facts.get("photon.candidate_scope")
    if scope == "event_leading":
        region = str(facts["photon.abcd_region"])
        identity = "tight" if facts["photon.id_state"] == "tight" else "non-tight"
        iso_radius = _number(float(facts["photon.isolation_radius"]))
        score_source = str(facts["photon.score_source"])
        if score_source != "recorded_bdt_score":
            raise ValueError(f"unsupported photon score source: {score_source}")
        bdt_operator = predicates["photon.bdt_class"]["operator"]
        bdt_definition = {
            "gt_field": "score > per-candidate calibrated tight threshold",
            "open_between_fields": "per-candidate calibrated bounded score band",
            "not_gt_field": "per-candidate tight-threshold complement",
        }[bdt_operator]
        iso_operator = predicates["photon.isolation_class"]["operator"]
        iso_symbol = "<" if iso_operator == "lt_field" else ">"
        lines.append(
            f"Event-leading {identity} photon (ABCD {region}; recorded BDT score)"
        )
        lines.append(f"Photon ID: {bdt_definition}")
        lines.append(
            (
                f"E_{{T}}^{{iso}}(R = {iso_radius}) {iso_symbol} per-candidate isolation threshold"
                if dialect == "root"
                else (
                    rf"$E_T^{{\mathrm{{iso}}}}(R={iso_radius}) {iso_symbol} "
                    rf"T_{{\mathrm{{iso}}}}$ (per-candidate threshold)"
                )
            )
        )
    elif scope == "inclusive_pairs":
        lines.append("Inclusive photon-jet pairs")
    else:
        raise ValueError(f"unsupported candidate scope: {scope}")
    lines.append(jet_line)
    # The required event predicate is checked by SelectionProgram. This is
    # its audience-facing grammar, not a second configurable selection.
    lines.append("Accepted producer events")
    return lines


def _visible_coverage_entry(semantic_key: str, slot: str) -> dict[str, Any]:
    """Map one active semantic node to its exact renderer destinations."""

    return {
        "semantic_key": semantic_key,
        "disposition": "visible",
        "targets": [
            {"dialect": "matplotlib", "slot": slot},
            {"dialect": "root_tlatex", "slot": slot},
        ],
    }


def _audit_coverage_entry(semantic_key: str, reason: str) -> dict[str, Any]:
    """Record a bound non-display semantic instead of silently omitting it."""

    return {
        "semantic_key": semantic_key,
        "disposition": "audit_only",
        "reason": reason,
        "targets": [],
    }


def _semantic_coverage(
    program: SelectionProgram,
    dataset: DatasetDescriptor,
) -> dict[str, Any]:
    """Compile complete active semantic coverage for the V1 xJgamma plot.

    Each active selection or dataset node is represented exactly once.  A
    node may be visible in both render dialects or explicitly audit-only; an
    unclassified node is never accepted by :func:`verify_plot_contract`.
    """

    predicates, facts = _program_maps(program)
    scope = str(facts["photon.candidate_scope"])
    jet_line = 2 if scope == "inclusive_pairs" else 4
    predicate_slots = {
        "event.accepted": f"annotations.cuts[{jet_line + 1}]",
        "photon.et": "annotations.cuts[0]",
        "photon.abs_eta": "annotations.cuts[0]",
        "jet.pt": f"annotations.cuts[{jet_line}]",
        "jet.abs_eta": f"annotations.cuts[{jet_line}]",
        "jet.radius": f"annotations.cuts[{jet_line}]",
        "recoil.delta_phi": f"annotations.cuts[{jet_line}]",
        "photon.bdt_class": "annotations.cuts[2]",
        "photon.isolation_class": "annotations.cuts[3]",
    }
    fact_slots = {
        "photon.candidate_scope": "annotations.cuts[1]",
        "photon.abcd_region": "annotations.cuts[1]",
        "photon.id_state": "annotations.cuts[1]",
        "photon.non_tight_definition": "annotations.cuts[2]",
        "photon.score_source": "annotations.cuts[1]",
        "photon.bdt_threshold_source": "annotations.cuts[2]",
        "photon.isolation_state": "annotations.cuts[3]",
        "photon.isolation_radius": "annotations.cuts[3]",
        "photon.isolation_threshold_source": "annotations.cuts[3]",
    }

    entries: list[dict[str, Any]] = [
        _visible_coverage_entry("dataset.collision_system", "annotations.dataset"),
        _visible_coverage_entry("dataset.sample_kind", "annotations.dataset"),
        _visible_coverage_entry("dataset.experiment_status", "annotations.experiment"),
        _visible_coverage_entry("observable.key", "axes.x_label"),
        _visible_coverage_entry("observable.yield_accumulator", "axes.y_label"),
        _audit_coverage_entry("plot.style_id", "bound_renderer_style_not_physics_text"),
        _audit_coverage_entry(
            "selection.leader_branch",
            "recorded_tree_crosscheck_or_explicit_offline_recomputation",
        ),
    ]
    energy_key = (
        "dataset.sqrt_s_gev"
        if dataset.collision_system == "pp"
        else "dataset.sqrt_s_nn_gev"
    )
    entries.append(_visible_coverage_entry(energy_key, "annotations.dataset"))
    if dataset.centrality_percent is not None:
        entries.append(
            _visible_coverage_entry("dataset.centrality_percent", "annotations.dataset")
        )
    if dataset.event_count is not None:
        if dataset.sample_kind == "engineering_fixture":
            entries.append(_visible_coverage_entry("dataset.event_count", "annotations.dataset"))
        else:
            entries.append(
                _audit_coverage_entry(
                    "dataset.event_count",
                    "bound_sample_inventory_not_a_physics_normalization_label",
                )
            )
    for key in predicates:
        entries.append(_visible_coverage_entry(key, predicate_slots[key]))
    for key in facts:
        if key == "photon.leader_rule":
            entries.append(
                _audit_coverage_entry(
                    key,
                    "deterministic_tie_break_bound_in_selection_program",
                )
            )
        else:
            entries.append(_visible_coverage_entry(key, fact_slots[key]))

    entries.sort(key=lambda entry: str(entry["semantic_key"]))
    keys = [str(entry["semantic_key"]) for entry in entries]
    if len(keys) != len(set(keys)):
        raise ValueError("semantic coverage assigns an active key more than once")
    return {"schema": SEMANTIC_COVERAGE_SCHEMA, "entries": entries}


def _visible_semantic_keys(coverage: dict[str, Any]) -> list[str]:
    return sorted(
        str(entry["semantic_key"])
        for entry in coverage["entries"]
        if entry["disposition"] == "visible"
    )


def _validate_semantic_coverage_targets(
    coverage: dict[str, Any],
    annotations: dict[str, Any],
) -> None:
    """Fail if a semantic entry has no real, exact renderer destination."""

    if set(coverage) != {"schema", "entries"}:
        raise ValueError("semantic coverage fields differ")
    if coverage["schema"] != SEMANTIC_COVERAGE_SCHEMA:
        raise ValueError("semantic coverage schema differs")
    entries = coverage["entries"]
    if not isinstance(entries, list):
        raise ValueError("semantic coverage entries must be a list")
    keys: list[str] = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("semantic coverage entry must be an object")
        disposition = entry.get("disposition")
        expected_fields = (
            {"semantic_key", "disposition", "targets"}
            if disposition == "visible"
            else {"semantic_key", "disposition", "reason", "targets"}
        )
        if set(entry) != expected_fields:
            raise ValueError("semantic coverage entry fields differ")
        key = entry.get("semantic_key")
        if not isinstance(key, str) or not key:
            raise ValueError("semantic coverage key is invalid")
        keys.append(key)
        targets = entry.get("targets")
        if disposition == "audit_only":
            if targets != [] or not isinstance(entry.get("reason"), str) or not entry["reason"]:
                raise ValueError("audit-only semantic coverage is invalid")
            continue
        if disposition != "visible" or not isinstance(targets, list) or len(targets) != 2:
            raise ValueError("visible semantic coverage is invalid")
        if {target.get("dialect") for target in targets if isinstance(target, dict)} != {
            "matplotlib",
            "root_tlatex",
        }:
            raise ValueError("visible semantic coverage must target both render dialects")
        for target in targets:
            if not isinstance(target, dict) or set(target) != {"dialect", "slot"}:
                raise ValueError("semantic coverage target fields differ")
            dialect = str(target["dialect"])
            slot = str(target["slot"])
            if slot in {"annotations.experiment", "annotations.dataset"}:
                annotation_key = slot.split(".", 1)[1]
                if not isinstance(annotations[dialect].get(annotation_key), str):
                    raise ValueError("semantic coverage annotation target is absent")
            elif slot.startswith("annotations.cuts[") and slot.endswith("]"):
                index_text = slot[len("annotations.cuts["):-1]
                if not index_text.isdigit():
                    raise ValueError("semantic coverage cut target is invalid")
                index = int(index_text)
                if index >= len(annotations[dialect]["cuts"]):
                    raise ValueError("semantic coverage cut target is absent")
            elif slot not in {"axes.x_label", "axes.y_label"}:
                raise ValueError("semantic coverage target slot is unsupported")
    if keys != sorted(keys) or len(keys) != len(set(keys)):
        raise ValueError("semantic coverage keys are not unique and sorted")


def compile_plot_contract(
    *,
    histogram_path: Path,
    histogram_receipt_path: Path,
    dataset_manifest_path: Path,
) -> dict[str, Any]:
    """Return a hash-bound plot contract with no free-form physics prose."""

    histogram_file = Path(histogram_path).resolve()
    receipt_file = Path(histogram_receipt_path).resolve()
    dataset_file = Path(dataset_manifest_path).resolve()
    payload = _load_mapping(histogram_file)
    receipt = _load_mapping(receipt_file)
    dataset_mapping = _load_mapping(dataset_file)
    if receipt.get("schema") != "PhotonJetHistogramReceiptV1":
        raise ValueError("histogram receipt schema is not PhotonJetHistogramReceiptV1")
    recorded_histogram = receipt.get("histogram", {})
    if recorded_histogram.get("sha256") != sha256_file(histogram_file):
        raise ValueError("histogram bytes do not match their receipt")
    selection = _selection_from_payload(payload)
    program = compile_recoil_selection(selection)
    if receipt.get("selection_program_sha256") != program.sha256:
        raise ValueError("payload selection does not match the executable selection receipt")
    if receipt.get("selection_program") != program.to_dict():
        raise ValueError("selection program body differs from its canonical compilation")

    edges = np.asarray(payload.get("axis", {}).get("edges", []), dtype=float)
    sumw = np.asarray(payload.get("sumw", []), dtype=float)
    sumw2 = np.asarray(payload.get("sumw2", []), dtype=float)
    if len(edges) != len(sumw) + 1 or len(sumw) != len(sumw2):
        raise ValueError("histogram arrays and bin edges are inconsistent")
    if np.any(~np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError("histogram edges are not finite and increasing")
    if np.any(~np.isfinite(sumw)) or np.any(~np.isfinite(sumw2)) or np.any(sumw2 < 0):
        raise ValueError("histogram contents are not finite physical accumulators")

    dataset = DatasetDescriptor.from_mapping(dataset_mapping)
    if receipt.get("dataset_manifest_sha256") != sha256_file(dataset_file):
        raise ValueError("dataset manifest differs from the histogram producer receipt")
    if receipt.get("dataset") != dataset.to_dict():
        raise ValueError("dataset descriptor differs from the histogram producer receipt")
    semantic_coverage = _semantic_coverage(program, dataset)
    contract: dict[str, Any] = {
        "schema": "PhotonJetPlotContractV1",
        "style_id": STYLE_ID,
        "plot_kind": "xjgamma_histogram",
        "histogram_sha256": sha256_file(histogram_file),
        "histogram_receipt_sha256": sha256_file(receipt_file),
        "dataset_manifest_sha256": sha256_file(dataset_file),
        "dataset": dataset.to_dict(),
        "observable": {
            "key": "xjgamma",
            "x_label_matplotlib": r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$",
            "x_label_root": "x_{J#gamma} = p_{T}^{jet}/p_{T}^{#gamma}",
            "y_label": "Weighted photon-jet pairs",
        },
        "selection_program": program.to_dict(),
        "selection_program_sha256": program.sha256,
        "annotations": {
            "matplotlib": {
                "experiment": MATPLOTLIB_EXPERIMENT_LABEL,
                "dataset": dataset.line("matplotlib"),
                "cuts": _annotation_lines(program, "matplotlib"),
            },
            "root_tlatex": {
                "experiment": ROOT_EXPERIMENT_LABEL,
                "dataset": dataset.line("root"),
                "cuts": _annotation_lines(program, "root"),
            },
        },
        "semantic_coverage": semantic_coverage,
        "semantic_coverage_sha256": _sha256_json(semantic_coverage),
        "visible_semantic_keys": _visible_semantic_keys(semantic_coverage),
    }
    contract["contract_sha256"] = _sha256_json(contract)
    return verify_plot_contract(contract)


def verify_plot_contract(contract: dict[str, Any]) -> dict[str, Any]:
    """Recompile every semantic field and reject any post-compile mutation."""

    expected_keys = {
        "schema", "style_id", "plot_kind", "histogram_sha256",
        "histogram_receipt_sha256", "dataset_manifest_sha256", "dataset",
        "observable", "selection_program", "selection_program_sha256",
        "annotations", "semantic_coverage", "semantic_coverage_sha256",
        "visible_semantic_keys", "contract_sha256",
    }
    if set(contract) != expected_keys:
        missing = sorted(expected_keys - set(contract))
        extra = sorted(set(contract) - expected_keys)
        raise ValueError(f"plot contract fields differ; missing={missing}, extra={extra}")
    if contract["schema"] != "PhotonJetPlotContractV1":
        raise ValueError("renderer requires PhotonJetPlotContractV1")
    body = {key: value for key, value in contract.items() if key != "contract_sha256"}
    observed_hash = _sha256_json(body)
    if contract["contract_sha256"] != observed_hash:
        raise ValueError("plot contract changed after compilation")
    if contract["style_id"] != STYLE_ID or contract["plot_kind"] != "xjgamma_histogram":
        raise ValueError("unsupported plot style or plot kind")

    dataset = DatasetDescriptor.from_contract_mapping(dict(contract["dataset"]))
    program = SelectionProgram.from_dict(contract["selection_program"])
    if contract["selection_program_sha256"] != program.sha256:
        raise ValueError("plot-contract selection identity is invalid")
    observable = {
        "key": "xjgamma",
        "x_label_matplotlib": r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$",
        "x_label_root": "x_{J#gamma} = p_{T}^{jet}/p_{T}^{#gamma}",
        "y_label": "Weighted photon-jet pairs",
    }
    if contract["observable"] != observable:
        raise ValueError("observable labels or identity were mutated")
    annotations = {
        "matplotlib": {
            "experiment": MATPLOTLIB_EXPERIMENT_LABEL,
            "dataset": dataset.line("matplotlib"),
            "cuts": _annotation_lines(program, "matplotlib"),
        },
        "root_tlatex": {
            "experiment": ROOT_EXPERIMENT_LABEL,
            "dataset": dataset.line("root"),
            "cuts": _annotation_lines(program, "root"),
        },
    }
    if contract["annotations"] != annotations:
        raise ValueError("plot annotations differ from their compiled semantic sources")
    semantic_coverage = _semantic_coverage(program, dataset)
    _validate_semantic_coverage_targets(semantic_coverage, annotations)
    if contract["semantic_coverage"] != semantic_coverage:
        raise ValueError("semantic coverage differs from the active plot lineage")
    coverage_hash = _sha256_json(semantic_coverage)
    if contract["semantic_coverage_sha256"] != coverage_hash:
        raise ValueError("semantic coverage identity was mutated")
    visible = _visible_semantic_keys(semantic_coverage)
    if contract["visible_semantic_keys"] != visible:
        raise ValueError("visible semantic lineage was mutated")
    return contract
