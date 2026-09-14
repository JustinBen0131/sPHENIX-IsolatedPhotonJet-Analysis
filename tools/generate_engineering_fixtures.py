#!/usr/bin/env python3
"""Regenerate the public, engineering-only PhotonJetTrees fixtures.

The compact JSON seed is the human-auditable fixture authority.  The branch
contract is the schema authority.  This tool derives repeated event context,
object joins, eventTree arrays, leaders, and typed defaults instead of storing
thousands of expanded fields, then writes ROOT through a path-free handle.

ROOT container timestamps are intentionally not physics identity.  They and
the UUID are nevertheless fixed for these tiny fixtures so the pinned runtime
reproduces both exact bytes and exact ordered record semantics.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Iterator
import uuid
from unittest.mock import patch

import awkward as ak
import numpy as np
import uproot


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "contracts/photonjet_trees_v1_branches.json"
DEFAULT_SEED = ROOT / "tests/fixtures/engineering_fixture_seed_v1.json"
DEFAULT_FIXTURE_DIR = ROOT / "tests/fixtures"

SCALAR_DTYPES = {
    "int32_t": "int32",
    "int64_t": "int64",
    "uint64_t": "uint64",
    "double": "float64",
}
SYSTEM_ORDER = ("auau", "pp")
FIXED_CONTAINER_DATETIME = datetime.datetime(2024, 1, 1, 0, 0, 0)


class _FixedDateTime(datetime.datetime):
    """Give ROOT metadata a stable timestamp during fixture generation."""

    @classmethod
    def now(cls, timezone: datetime.tzinfo | None = None) -> "_FixedDateTime":
        value = cls(*FIXED_CONTAINER_DATETIME.timetuple()[:6])
        return value if timezone is None else value.replace(tzinfo=timezone)


@contextmanager
def _fixed_root_container_time() -> Iterator[None]:
    # Uproot writes TKey/TTree dates through ``datetime.datetime.now``.  The
    # pinned fixture runtime plus this narrow patch makes repeated generation
    # byte-exact without changing timestamps anywhere outside this context.
    with patch.object(datetime, "datetime", _FixedDateTime):
        yield


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact(path: Path, *, repo_root: Path = ROOT) -> dict[str, Any]:
    resolved = path.resolve()
    try:
        name = resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(f"artifact is outside the repository: {resolved}") from error
    return {
        "path": name,
        "sha256": _sha256(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def load_contract(path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    contract = _load_json(path)
    if contract.get("schema") != "PhotonJetTreeBranchContractV1":
        raise ValueError("unsupported branch-contract schema")
    tree_order = contract.get("tree_order")
    branches = contract.get("branches")
    if not isinstance(tree_order, list) or not isinstance(branches, dict):
        raise ValueError("branch contract lacks tree_order or branches")
    if list(branches) != tree_order:
        raise ValueError("branch map order differs from tree_order")
    for tree_name in tree_order:
        names: list[str] = []
        for branch in branches[tree_name]:
            name = branch.get("name")
            typename = branch.get("typename")
            base = typename[:-2] if isinstance(typename, str) and typename.endswith("[]") else typename
            if not isinstance(name, str) or base not in SCALAR_DTYPES:
                raise ValueError(f"unsupported branch declaration in {tree_name}: {branch}")
            names.append(name)
        if len(names) != len(set(names)):
            raise ValueError(f"duplicate branch name in {tree_name}")
    return contract


def _data_branch_names(branches: Iterable[dict[str, str]]) -> set[str]:
    """Return scalar counters that uproot generates for jagged branches."""

    return {
        f"n{branch['name']}"
        for branch in branches
        if branch["typename"].endswith("[]")
    }


def _ordered_record(
    contract: dict[str, Any],
    tree_name: str,
    values: dict[str, Any],
) -> dict[str, Any]:
    branches = contract["branches"][tree_name]
    generated_counters = _data_branch_names(branches)
    fields = [
        branch["name"]
        for branch in branches
        if branch["name"] not in generated_counters
    ]
    missing = [name for name in fields if name not in values]
    if missing:
        raise ValueError(f"{tree_name} materialization lacks fields: {', '.join(missing)}")
    return {name: values[name] for name in fields}


def _validate_seed(seed: dict[str, Any]) -> None:
    if seed.get("schema") != "PhotonJetEngineeringFixtureSeedV1":
        raise ValueError("unsupported engineering-fixture seed schema")
    if seed.get("physics_use") != "FORBIDDEN_ENGINEERING_ONLY":
        raise ValueError("engineering fixture seed must forbid physics use")
    systems = seed.get("systems")
    events = seed.get("events")
    if not isinstance(systems, dict) or set(systems) != {"pp", "auau"}:
        raise ValueError("fixture seed must define exactly pp and auau")
    if not isinstance(events, list) or len(events) != 2:
        raise ValueError("fixture seed must define exactly two events")
    photon_count = 0
    for event in events:
        if len(event.get("photons", [])) != 2 or len(event.get("jets", [])) != 2:
            raise ValueError("each fixture event must define two photons and two jets")
        for photon in event["photons"]:
            if len(photon.get("native", [])) != 7:
                raise ValueError("each fixture photon must define seven native observables")
            photon_count += 1
    bdt_fields = {"score", "tight", "nontight_low", "nontight_high"}
    for system, profile in systems.items():
        if len(profile.get("centrality", [])) != len(events):
            raise ValueError(f"{system} centrality count differs from events")
        if profile.get("model_input_count") not in (11, 14):
            raise ValueError(f"{system} model input count must be 11 or 14")
        bdt_rows = profile.get("bdt")
        if not isinstance(bdt_rows, list) or len(bdt_rows) != photon_count:
            raise ValueError(f"{system} BDT row count differs from photons")
        if any(not isinstance(row, dict) or set(row) != bdt_fields for row in bdt_rows):
            raise ValueError(f"{system} BDT row shape is invalid")


def _event_values(
    *,
    event_seed: dict[str, Any],
    event_index: int,
    system_profile: dict[str, Any],
) -> dict[str, Any]:
    sequence = event_index + 1
    return {
        "source_file_index": 0,
        "source_entry": event_index,
        "event_id_hi": 0,
        "event_id_lo": event_seed["event_id_lo"],
        "run": event_seed["run"],
        "event_sequence": sequence,
        "physical_event_sequence": sequence,
        "trigger_bits": 0,
        "live_trigger_bits": 0,
        "scaled_trigger_bits": system_profile["scaled_trigger_bits"],
        "scaled_bit30": system_profile["scaled_bit30"],
        "vertex_z": event_seed["vertex_z"],
        "centrality": system_profile["centrality"][event_index],
        "event_weight": 1.0,
        "total_calo_energy": event_seed["total_calo_energy"],
        "terminal_status": 0,
    }


def _bdt_inputs(global_photon_index: int, active_count: int) -> list[float | None]:
    step = 0.002 * global_photon_index
    raw = [
        0.02 + step,
        0.02 + step,
        0.01 + step,
        0.01 + step,
        0.03 + step,
        0.03 + step,
        0.02 + step,
        0.02 + step,
        0.90 + step,
        0.92 + step,
        0.25 + step,
        0.10 + step,
        0.15 + step,
        0.20 + step,
    ]
    return [
        float(np.float32(value)) if index < active_count else None
        for index, value in enumerate(raw)
    ]


def _leader_index(photons: list[dict[str, Any]], state: str) -> int:
    def accepts(row: dict[str, Any]) -> bool:
        isolated = bool(row["iso_r04_pass"])
        nonisolated = row["iso_r04"] > row["iso_r04_nonisolated_threshold"]
        if state == "A":
            return bool(row["bdt_is_tight"]) and isolated
        if state == "B":
            return bool(row["bdt_is_tight"]) and nonisolated
        if state == "C":
            return bool(row["bdt_is_nontight"]) and isolated
        if state == "D":
            return bool(row["bdt_is_nontight"]) and nonisolated
        if state == "C_complement":
            return bool(row["bdt_is_not_tight"]) and isolated
        if state == "D_complement":
            return bool(row["bdt_is_not_tight"]) and nonisolated
        raise ValueError(f"unsupported leader state: {state}")

    candidates = [index for index, row in enumerate(photons) if accepts(row)]
    return max(candidates, key=lambda index: photons[index]["photon_et"], default=-1)


def _event_tree_source(name: str) -> tuple[str, str]:
    prefix = name.split("_", 1)[0]
    if prefix == "photon":
        keep_prefix = {
            "photon_encounter_ordinal",
            "photon_et",
            "photon_eta",
            "photon_phi",
        }
        return prefix, name if name in keep_prefix else name.removeprefix("photon_")
    if prefix == "jet":
        return prefix, name
    if prefix == "pair":
        return prefix, name.removeprefix("pair_")
    raise ValueError(f"unsupported eventTree array branch: {name}")


def materialize_records(
    seed: dict[str, Any],
    contract: dict[str, Any],
) -> dict[str, Any]:
    """Expand the compact semantic seed into the exact eight-tree records."""

    _validate_seed(seed)
    bases = seed["id_bases"]
    fixed = seed["fixed_values"]
    event_seeds = seed["events"]
    materialized: dict[str, Any] = {
        "schema": "PhotonJetEngineeringFixtureRecordsV1",
        "branch_contract": "contracts/photonjet_trees_v1_branches.json",
        "fixture_class": seed["fixture_class"],
        "physics_use": "FORBIDDEN_ENGINEERING_ONLY",
        "null_float_semantics": "IEEE_NAN",
        "systems": {},
    }

    for system in ("pp", "auau"):
        profile = seed["systems"][system]
        offset = int(profile["id_offset"])
        trees: dict[str, list[dict[str, Any]]] = {
            name: [] for name in contract["tree_order"]
        }
        event_groups: list[
            tuple[
                dict[str, Any],
                list[dict[str, Any]],
                list[dict[str, Any]],
                list[dict[str, Any]],
            ]
        ] = []
        photon_ordinal = 0
        jet_ordinal = 0
        pair_ordinal = 0

        for event_index, event_seed in enumerate(event_seeds):
            event = _ordered_record(
                contract,
                "events",
                _event_values(
                    event_seed=event_seed,
                    event_index=event_index,
                    system_profile=profile,
                ),
            )
            trees["events"].append(event)
            event_photons: list[dict[str, Any]] = []
            event_jets: list[dict[str, Any]] = []
            event_pairs: list[dict[str, Any]] = []

            for local_index, photon_seed in enumerate(event_seed["photons"]):
                bdt = profile["bdt"][photon_ordinal]
                inputs = _bdt_inputs(photon_ordinal, profile["model_input_count"])
                isolated = photon_seed["isolation"] < fixed["isolation_threshold"]
                truth_matched = local_index == 0
                native = photon_seed["native"]
                values = {
                    **event,
                    "candidate_id_hi": 0,
                    "candidate_id_lo": bases["photon"] + offset + photon_ordinal + 1,
                    "photon_encounter_ordinal": local_index,
                    "photon_et": photon_seed["et"],
                    "photon_eta": photon_seed["eta"],
                    "photon_phi": photon_seed["phi"],
                    "bdt_score": bdt["score"],
                    "bdt_tight_threshold": bdt["tight"],
                    "bdt_nontight_low_threshold": bdt["nontight_low"],
                    "bdt_nontight_high_threshold": bdt["nontight_high"],
                    "bdt_is_tight": int(bdt["score"] > bdt["tight"]),
                    "bdt_is_nontight": int(
                        bdt["nontight_low"] < bdt["score"] < bdt["nontight_high"]
                    ),
                    "bdt_is_not_tight": int(bdt["score"] <= bdt["tight"]),
                    "bdt_input_count": profile["model_input_count"],
                    "iso_r03": photon_seed["isolation"],
                    "iso_r03_threshold": fixed["isolation_threshold"],
                    "iso_r03_nonisolated_threshold": fixed["nonisolated_threshold"],
                    "iso_r03_pass": int(isolated),
                    "iso_r04": photon_seed["isolation"],
                    "iso_r04_threshold": fixed["isolation_threshold"],
                    "iso_r04_nonisolated_threshold": fixed["nonisolated_threshold"],
                    "iso_r04_pass": int(isolated),
                    "truth_matched": int(truth_matched),
                    "truth_barcode": event_seed["truth_photon"]["barcode"] if truth_matched else -1,
                    "truth_generator_occurrence_embedding_id": 2 if truth_matched else -1,
                    "native_weta_cogx": native[0],
                    "native_wphi_cogx": native[1],
                    "native_weta33_cogx": native[2],
                    "native_wphi33_cogx": native[3],
                    "native_e11_over_e33": native[4],
                    "native_e32_over_e35": native[5],
                    "native_et1": native[6],
                    **{f"bdt_input_{index:02d}": value for index, value in enumerate(inputs)},
                }
                photon = _ordered_record(contract, "photons", values)
                trees["photons"].append(photon)
                event_photons.append(photon)
                photon_ordinal += 1

            for local_index, jet_seed in enumerate(event_seed["jets"]):
                values = {
                    **event,
                    "jet_id_hi": 0,
                    "jet_id_lo": bases["jet"] + offset + jet_ordinal + 1,
                    "jet_radius": fixed["jet_radius"],
                    "jet_raw_pt": jet_seed["raw_pt"],
                    "jet_pt": jet_seed["pt"],
                    "jet_eta": jet_seed["eta"],
                    "jet_phi": jet_seed["phi"],
                    "jet_mass": fixed["jet_mass"],
                    "jet_area": fixed["jet_area"],
                    "jet_quality_bitmask": 0,
                    "jet_order": local_index,
                }
                jet = _ordered_record(contract, "jets", values)
                trees["jets"].append(jet)
                event_jets.append(jet)
                jet_ordinal += 1

            for photon_index, photon in enumerate(event_photons):
                for jet_index, jet in enumerate(event_jets):
                    values = {
                        **photon,
                        **jet,
                        "pair_id_hi": 0,
                        "pair_id_lo": bases["pair"] + offset + pair_ordinal + 1,
                        "photon_index": photon_index,
                        "jet_index": jet_index,
                        "delta_phi": fixed["pair_delta_phi"],
                        "xjgamma": jet["jet_pt"] / photon["photon_et"],
                        "recoil_state": 1,
                        "photon_rank": photon_index,
                        "jet_rank": jet_index,
                        "wrong_photon_class": 0,
                        "wrong_recoil_class": 0,
                    }
                    pair = _ordered_record(contract, "photonJets", values)
                    trees["photonJets"].append(pair)
                    event_pairs.append(pair)
                    pair_ordinal += 1
            event_groups.append((event, event_photons, event_jets, event_pairs))

        for event_index, (event, photons, jets, pairs) in enumerate(event_groups):
            event_tree_values: dict[str, Any] = {
                **event,
                "nphotons": len(photons),
                "njets": len(jets),
                "npairs": len(pairs),
            }
            source_groups = {"photon": photons, "jet": jets, "pair": pairs}
            for branch in contract["branches"]["eventTree"]:
                name = branch["name"]
                if not branch["typename"].endswith("[]"):
                    continue
                prefix, source_field = _event_tree_source(name)
                event_tree_values[name] = [row[source_field] for row in source_groups[prefix]]
            event_tree_values.update(
                {
                    "leader_A_r04_index": _leader_index(photons, "A"),
                    "leader_B_r04_index": _leader_index(photons, "B"),
                    "leader_C_r04_index": _leader_index(photons, "C"),
                    "leader_D_r04_index": _leader_index(photons, "D"),
                    "leader_C_r04_complement_index": _leader_index(photons, "C_complement"),
                    "leader_D_r04_complement_index": _leader_index(photons, "D_complement"),
                }
            )
            trees["eventTree"].append(
                _ordered_record(contract, "eventTree", event_tree_values)
            )

            truth_photon_seed = event_seeds[event_index]["truth_photon"]
            trees["truthPhotons"].append(
                _ordered_record(
                    contract,
                    "truthPhotons",
                    {
                        "source_file_index": 0,
                        "event_id_hi": 0,
                        "event_id_lo": event["event_id_lo"],
                        "truth_photon_id_hi": 0,
                        "truth_photon_id_lo": bases["truth_photon"] + offset + event_index + 1,
                        "truth_photon_pt": truth_photon_seed["pt"],
                        "truth_photon_eta": truth_photon_seed["eta"],
                        "truth_photon_phi": truth_photon_seed["phi"],
                        "prompt_class": 1,
                        "source_role": 1,
                        "generator_barcode": truth_photon_seed["barcode"],
                        "generator_occurrence_embedding_id": 2,
                        "truth_isolation": truth_photon_seed["isolation"],
                    },
                )
            )
            leading_jet = jets[0]
            trees["truthJets"].append(
                _ordered_record(
                    contract,
                    "truthJets",
                    {
                        "source_file_index": 0,
                        "event_id_hi": 0,
                        "event_id_lo": event["event_id_lo"],
                        "truth_jet_id_hi": 0,
                        "truth_jet_id_lo": bases["truth_jet"] + offset + event_index + 1,
                        "truth_jet_radius": fixed["jet_radius"],
                        "truth_jet_pt": leading_jet["jet_raw_pt"],
                        "truth_jet_eta": leading_jet["jet_eta"],
                        "truth_jet_phi": leading_jet["jet_phi"],
                    },
                )
            )

        overlay = seed["overlay_truth_photon"]
        overlay_event = event_groups[overlay["event_index"]][0]
        overlay_barcode = event_seeds[overlay["event_index"]]["truth_photon"]["barcode"]
        trees["truthPhotons"].append(
            _ordered_record(
                contract,
                "truthPhotons",
                {
                    "source_file_index": 0,
                    "event_id_hi": 0,
                    "event_id_lo": overlay_event["event_id_lo"],
                    "truth_photon_id_hi": 0,
                    "truth_photon_id_lo": bases["truth_photon_overlay"] + offset,
                    "truth_photon_pt": overlay["pt"],
                    "truth_photon_eta": overlay["eta"],
                    "truth_photon_phi": overlay["phi"],
                    "prompt_class": 1,
                    "source_role": 0,
                    "generator_barcode": overlay_barcode,
                    "generator_occurrence_embedding_id": overlay["occurrence"],
                    "truth_isolation": overlay["isolation"],
                },
            )
        )

        link_ordinal = 0
        for reco_type, truth_type, reco_tree, truth_tree in (
            (1, 1, "photons", "truthPhotons"),
            (2, 2, "jets", "truthJets"),
        ):
            for event_index, (_, photons, jets, _) in enumerate(event_groups):
                reco = (photons if reco_tree == "photons" else jets)[0]
                truth = trees[truth_tree][event_index]
                reco_prefix = "candidate" if reco_tree == "photons" else "jet"
                truth_prefix = "truth_photon" if truth_tree == "truthPhotons" else "truth_jet"
                trees["recoTruthLinks"].append(
                    _ordered_record(
                        contract,
                        "recoTruthLinks",
                        {
                            "source_file_index": 0,
                            "event_id_hi": 0,
                            "event_id_lo": reco["event_id_lo"],
                            "link_id_hi": 0,
                            "link_id_lo": bases["link"] + offset + link_ordinal + 1,
                            "reco_type": reco_type,
                            "reco_id_hi": reco[f"{reco_prefix}_id_hi"],
                            "reco_id_lo": reco[f"{reco_prefix}_id_lo"],
                            "reco_index": 0,
                            "truth_type": truth_type,
                            "truth_id_hi": truth[f"{truth_prefix}_id_hi"],
                            "truth_id_lo": truth[f"{truth_prefix}_id_lo"],
                            "truth_index": 0,
                            "match_metric": 0.01 * (link_ordinal + 1),
                            "link_class": 0,
                        },
                    )
                )
                link_ordinal += 1

        materialized["systems"][system] = {
            "output_path": f"tests/fixtures/photonjet_trees_{system}.root",
            "collision_system": system,
            "trees": trees,
        }

    validate_records(materialized, contract)
    return materialized


def _writing_schema(branches: list[dict[str, str]]) -> dict[str, str]:
    generated_counters = _data_branch_names(branches)
    schema: dict[str, str] = {}
    for branch in branches:
        name = branch["name"]
        typename = branch["typename"]
        if name in generated_counters:
            continue
        if typename.endswith("[]"):
            schema[name] = f"var * {SCALAR_DTYPES[typename[:-2]]}"
        else:
            schema[name] = SCALAR_DTYPES[typename]
    return schema


def _validate_and_convert_column(
    *,
    tree_name: str,
    branch_name: str,
    dtype: str,
    rows: list[dict[str, Any]],
) -> np.ndarray | ak.Array:
    values = [row[branch_name] for row in rows]
    if dtype.startswith("var * "):
        if any(not isinstance(value, list) for value in values):
            raise ValueError(f"{tree_name}.{branch_name} must be a list in every record")
        if any(item is None for value in values for item in value):
            raise ValueError(f"{tree_name}.{branch_name} arrays may not contain null")
        return ak.values_astype(ak.Array(values), dtype.removeprefix("var * "))

    if any(value is None for value in values):
        if dtype != "float64":
            raise ValueError(f"null is only valid for floating fields: {tree_name}.{branch_name}")
        values = [math.nan if value is None else value for value in values]
    return np.asarray(values, dtype=dtype)


def validate_records(
    records: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    if records.get("schema") != "PhotonJetEngineeringFixtureRecordsV1":
        raise ValueError("unsupported engineering-fixture record schema")
    if records.get("physics_use") != "FORBIDDEN_ENGINEERING_ONLY":
        raise ValueError("engineering fixture must forbid physics use")
    if records.get("null_float_semantics") != "IEEE_NAN":
        raise ValueError("unsupported null-float semantics")
    systems = records.get("systems")
    if not isinstance(systems, dict) or set(systems) != {"pp", "auau"}:
        raise ValueError("fixture records must define exactly pp and auau")

    tree_order = contract["tree_order"]
    for system, specification in systems.items():
        if specification.get("collision_system") != system:
            raise ValueError(f"collision-system mismatch for {system}")
        expected_output = f"tests/fixtures/photonjet_trees_{system}.root"
        if specification.get("output_path") != expected_output:
            raise ValueError(f"unexpected fixture output path for {system}")
        trees = specification.get("trees")
        if not isinstance(trees, dict) or list(trees) != tree_order:
            raise ValueError(f"tree order differs in {system} record source")
        for tree_name in tree_order:
            rows = trees[tree_name]
            if not isinstance(rows, list):
                raise ValueError(f"{system}.{tree_name} must contain a record list")
            expected = [branch["name"] for branch in contract["branches"][tree_name]]
            generated = _data_branch_names(contract["branches"][tree_name])
            record_fields = [name for name in expected if name not in generated]
            for row_index, row in enumerate(rows):
                if not isinstance(row, dict) or list(row) != record_fields:
                    raise ValueError(
                        f"{system}.{tree_name} row {row_index} fields differ from the contract"
                    )


def write_fixture(
    *,
    output: Path,
    system: str,
    specification: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        temporary = tempfile.NamedTemporaryFile(
            mode="w+b",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".partial",
            delete=False,
        )
        temporary_name = temporary.name
        stable_uuid = uuid.uuid5(
            uuid.NAMESPACE_URL,
            "sphenix-photonjet-engineering-fixture:"
            f"{system}:{records_identity(specification)}",
        )
        with _fixed_root_container_time():
            with uproot.recreate(temporary, uuid_function=lambda: stable_uuid) as root:
                for tree_name in contract["tree_order"]:
                    branches = contract["branches"][tree_name]
                    rows = specification["trees"][tree_name]
                    schema = _writing_schema(branches)
                    tree = root.mktree(
                        tree_name,
                        schema,
                        initial_basket_capacity=max(2, min(10, len(rows))),
                    )
                    columns = {
                        name: _validate_and_convert_column(
                            tree_name=tree_name,
                            branch_name=name,
                            dtype=dtype,
                            rows=rows,
                        )
                        for name, dtype in schema.items()
                    }
                    tree.extend(columns)
        os.replace(temporary_name, output)
        output.chmod(0o644)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def records_identity(specification: dict[str, Any]) -> str:
    encoded = json.dumps(
        specification,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _normalize_root_value(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, list):
        return [_normalize_root_value(item) for item in value]
    if isinstance(value, dict):
        return {key: _normalize_root_value(item) for key, item in value.items()}
    return value


def check_fixture(
    *,
    path: Path,
    specification: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    with uproot.open(path) as root:
        observed_trees = root.keys(cycle=False)
        if observed_trees != contract["tree_order"]:
            raise ValueError(f"{path.name} tree order differs")
        for tree_name in contract["tree_order"]:
            branches = contract["branches"][tree_name]
            expected_schema = [
                (branch["name"], branch["typename"])
                for branch in branches
            ]
            if list(root[tree_name].typenames().items()) != expected_schema:
                raise ValueError(f"{path.name}:{tree_name} branch schema differs")
            generated_counters = _data_branch_names(branches)
            record_fields = [
                branch["name"]
                for branch in branches
                if branch["name"] not in generated_counters
            ]
            observed = _normalize_root_value(ak.to_list(
                root[tree_name].arrays(record_fields, library="ak")
            ))
            if observed != specification["trees"][tree_name]:
                raise ValueError(f"{path.name}:{tree_name} records differ")


def _fixture_rows(fixture_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for system in SYSTEM_ORDER:
        path = fixture_dir / f"photonjet_trees_{system}.root"
        row = _artifact(path)
        row.update(
            {
                "collision_system": system,
                "kind": "synthetic_public_tree_fixture",
            }
        )
        rows.append(row)
    return rows


def _write_receipts(
    *,
    contract_path: Path,
    seed_path: Path,
    fixture_dir: Path,
) -> None:
    fixture_rows = _fixture_rows(fixture_dir)
    source_seed = _artifact(seed_path)
    contract = _artifact(contract_path)
    generator = _artifact(Path(__file__))
    generation = {
        "schema": "PhotonJetEngineeringFixtureGenerationReceiptV2",
        "fixture_class": "sanitized_redistribution_safe_synthetic_records",
        "public_regeneration": "COMPACT_SEED_AND_GENERATOR_INCLUDED",
        "semantic_regeneration": "EXACT_MATERIALIZED_ORDERED_SCHEMA_AND_RECORD_VALUES",
        "root_container_byte_reproducibility": "BYTE_EXACT_WITH_PINNED_RUNTIME",
        "root_container_identity": {
            "fixed_datetime_utc": "2024-01-01T00:00:00Z",
            "uuid_derivation": "UUID5_NAMESPACE_URL_FROM_SYSTEM_RECORD_IDENTITY",
        },
        "physics_use": "FORBIDDEN_ENGINEERING_ONLY",
        "container_metadata_policy": "NO_MACHINE_PATHS",
        "negative_metadata_scan": "REQUIRED_BY_REPOSITORY_VALIDATOR",
        "branch_contract": contract,
        "source_seed": source_seed,
        "generator": generator,
        "runtime": {
            "uproot": uproot.__version__,
            "awkward": ak.__version__,
            "numpy": np.__version__,
        },
        "fixtures": fixture_rows,
    }
    manifest = {
        "schema": "PhotonJetEngineeringFixtureManifestV2",
        "generation_receipt": "tests/fixtures/GENERATION_RECEIPT.json",
        "source_seed": source_seed,
        "fixtures": fixture_rows,
    }
    (fixture_dir / "GENERATION_RECEIPT.json").write_text(
        json.dumps(generation, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    (fixture_dir / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def write_all(
    *,
    contract_path: Path = DEFAULT_CONTRACT,
    seed_path: Path = DEFAULT_SEED,
    fixture_dir: Path = DEFAULT_FIXTURE_DIR,
) -> None:
    contract = load_contract(contract_path)
    records = materialize_records(_load_json(seed_path), contract)
    for system in SYSTEM_ORDER:
        write_fixture(
            output=fixture_dir / f"photonjet_trees_{system}.root",
            system=system,
            specification=records["systems"][system],
            contract=contract,
        )
        check_fixture(
            path=fixture_dir / f"photonjet_trees_{system}.root",
            specification=records["systems"][system],
            contract=contract,
        )
    _write_receipts(
        contract_path=contract_path,
        seed_path=seed_path,
        fixture_dir=fixture_dir,
    )


def check_all(
    *,
    contract_path: Path = DEFAULT_CONTRACT,
    seed_path: Path = DEFAULT_SEED,
    fixture_dir: Path = DEFAULT_FIXTURE_DIR,
) -> None:
    contract = load_contract(contract_path)
    records = materialize_records(_load_json(seed_path), contract)
    for system in SYSTEM_ORDER:
        check_fixture(
            path=fixture_dir / f"photonjet_trees_{system}.root",
            specification=records["systems"][system],
            contract=contract,
        )
    with tempfile.TemporaryDirectory() as temporary:
        regeneration_root = Path(temporary)
        for system in SYSTEM_ORDER:
            regenerated = regeneration_root / f"photonjet_trees_{system}.root"
            write_fixture(
                output=regenerated,
                system=system,
                specification=records["systems"][system],
                contract=contract,
            )
            checked_in = fixture_dir / regenerated.name
            if _sha256(regenerated) != _sha256(checked_in):
                raise ValueError(
                    f"{regenerated.name} is not byte-exact under the pinned fixture runtime"
                )
    manifest = _load_json(fixture_dir / "MANIFEST.json")
    generation = _load_json(fixture_dir / "GENERATION_RECEIPT.json")
    if manifest.get("schema") != "PhotonJetEngineeringFixtureManifestV2":
        raise ValueError("unsupported fixture manifest")
    if generation.get("schema") != "PhotonJetEngineeringFixtureGenerationReceiptV2":
        raise ValueError("unsupported fixture generation receipt")
    if manifest.get("fixtures") != _fixture_rows(fixture_dir):
        raise ValueError("fixture manifest does not bind the checked-in ROOT bytes")
    if generation.get("fixtures") != manifest["fixtures"]:
        raise ValueError("generation receipt fixture rows differ from the manifest")
    if manifest.get("source_seed") != _artifact(seed_path):
        raise ValueError("fixture manifest does not bind the semantic seed")
    if generation.get("source_seed") != manifest["source_seed"]:
        raise ValueError("generation receipt semantic seed differs from the manifest")
    if generation.get("branch_contract") != _artifact(contract_path):
        raise ValueError("generation receipt does not bind the branch contract")
    if generation.get("generator") != _artifact(Path(__file__)):
        raise ValueError("generation receipt does not bind this generator")
    expected_runtime = {
        "uproot": uproot.__version__,
        "awkward": ak.__version__,
        "numpy": np.__version__,
    }
    if generation.get("runtime") != expected_runtime:
        raise ValueError("fixture runtime differs from the generation receipt")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true", help="regenerate fixtures and receipts")
    action.add_argument("--check", action="store_true", help="verify checked-in fixtures")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--seed", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--fixture-dir", type=Path, default=DEFAULT_FIXTURE_DIR)
    args = parser.parse_args(argv)
    if args.write:
        write_all(
            contract_path=args.contract,
            seed_path=args.seed,
            fixture_dir=args.fixture_dir,
        )
    else:
        check_all(
            contract_path=args.contract,
            seed_path=args.seed,
            fixture_dir=args.fixture_dir,
        )
    print("ENGINEERING_FIXTURES=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
