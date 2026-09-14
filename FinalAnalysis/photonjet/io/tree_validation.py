#!/usr/bin/env python3
"""Validate self-contained PhotonJetTrees_v1 files.

This validator starts at the public tree boundary. Migration adapters and
their private normalized inputs are deliberately not part of this package.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from importlib import resources
import json
from pathlib import Path
from typing import Any, Iterable

import awkward as ak
import numpy as np
import uproot

from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.selection import compile_recoil_selection


TREE_ORDER = (
    "events", "eventTree", "photons", "jets", "photonJets",
    "truthPhotons", "truthJets", "recoTruthLinks",
)
IDENTITIES = {
    "events": ("event_id_hi", "event_id_lo"),
    "photons": ("candidate_id_hi", "candidate_id_lo"),
    "jets": ("jet_id_hi", "jet_id_lo"),
    "photonJets": ("pair_id_hi", "pair_id_lo"),
    "truthPhotons": ("truth_photon_id_hi", "truth_photon_id_lo"),
    "truthJets": ("truth_jet_id_hi", "truth_jet_id_lo"),
    "recoTruthLinks": ("link_id_hi", "link_id_lo"),
}

_LINK_MATCH = 0
_LINK_RECO_FAKE = 1
_LINK_TRUTH_MISS = 2
_LINK_WRONG_PHOTON = 3
_LINK_WRONG_RECOIL = 4
_LINK_MATCH_CANDIDATE = 5
_LINK_CLASSES = {
    _LINK_MATCH,
    _LINK_RECO_FAKE,
    _LINK_TRUTH_MISS,
    _LINK_WRONG_PHOTON,
    _LINK_WRONG_RECOIL,
    _LINK_MATCH_CANDIDATE,
}


def _key(high: Any, low: Any) -> tuple[int, int]:
    return int(high), int(low)


def _contract() -> dict[str, Any]:
    try:
        contract_text = (
            resources.files("photonjet_contracts")
            .joinpath("photonjet_trees_v1_branches.json")
            .read_text(encoding="utf-8")
        )
    except ModuleNotFoundError:
        # Repository-local execution before installation. A built wheel takes
        # the importlib.resources path above and is tested separately.
        contract_text = (
            Path(__file__).resolve().parents[3]
            / "contracts"
            / "photonjet_trees_v1_branches.json"
        ).read_text(encoding="utf-8")
    payload = json.loads(contract_text)
    if payload.get("schema") != "PhotonJetTreeBranchContractV1":
        raise ValueError("unsupported branch-contract schema")
    if tuple(payload.get("tree_order", [])) != TREE_ORDER:
        raise ValueError("branch-contract tree order differs from PhotonJetTrees_v1")
    return payload


def _validate_schema(root: uproot.ReadOnlyDirectory) -> None:
    observed_trees = tuple(root.keys(cycle=False))
    if observed_trees != TREE_ORDER:
        raise ValueError(f"tree order differs: {observed_trees} != {TREE_ORDER}")
    contract = _contract()
    for tree_name in TREE_ORDER:
        expected = [
            (item["name"], item["typename"])
            for item in contract["branches"][tree_name]
        ]
        observed = list(root[tree_name].typenames().items())
        if observed != expected:
            first = next(
                (
                    index
                    for index in range(min(len(observed), len(expected)))
                    if observed[index] != expected[index]
                ),
                min(len(observed), len(expected)),
            )
            raise ValueError(
                f"{tree_name} branch contract differs at index {first}: "
                f"observed={observed[first:first + 1]} expected={expected[first:first + 1]}"
            )


def _identity_set(root: uproot.ReadOnlyDirectory, tree_name: str) -> set[tuple[int, int]]:
    high_name, low_name = IDENTITIES[tree_name]
    arrays = root[tree_name].arrays([high_name, low_name], library="np")
    identities = {
        _key(high, low)
        for high, low in zip(arrays[high_name], arrays[low_name])
    }
    if len(identities) != root[tree_name].num_entries:
        raise ValueError(f"duplicate identity in {tree_name}")
    return identities


def _validate_parent_events(root: uproot.ReadOnlyDirectory) -> None:
    events = root["events"].arrays(library="np")
    event_rows: dict[tuple[int, int], int] = {}
    for row, (source, high, low) in enumerate(
        zip(events["source_file_index"], events["event_id_hi"], events["event_id_lo"])
    ):
        identity = _key(high, low)
        if identity in event_rows:
            raise ValueError(f"duplicate event identity: {identity}")
        if int(source) < 0:
            raise ValueError("source_file_index must be nonnegative")
        event_rows[identity] = row

    for tree_name in TREE_ORDER[2:]:
        arrays = root[tree_name].arrays(library="np")
        parent_rows: list[int] = []
        for row, (high, low) in enumerate(
            zip(arrays["event_id_hi"], arrays["event_id_lo"])
        ):
            event = _key(high, low)
            if event not in event_rows:
                raise ValueError(f"{tree_name} row {row} has no parent event")
            parent_rows.append(event_rows[event])
        parent_indices = np.asarray(parent_rows, dtype="int64")
        for field in sorted(set(events).intersection(arrays)):
            if not _same(arrays[field], events[field][parent_indices]):
                raise ValueError(
                    f"{tree_name} event field disagrees with its parent: {field}"
                )


def _group_indices(arrays: Any) -> dict[tuple[int, int], list[int]]:
    result: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, (high, low) in enumerate(zip(arrays["event_id_hi"], arrays["event_id_lo"])):
        result[_key(high, low)].append(index)
    return result


def _same(left: Any, right: Any) -> bool:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    if left_array.dtype.kind == "f" or right_array.dtype.kind == "f":
        return bool(np.array_equal(left_array, right_array, equal_nan=True))
    return bool(np.array_equal(left_array, right_array))


def _validate_event_tree(root: uproot.ReadOnlyDirectory) -> None:
    event_tree = root["eventTree"].arrays(library="ak")
    events = root["events"].arrays(library="np")
    if len(event_tree) != len(events["event_id_hi"]):
        raise ValueError("eventTree count differs from events")
    for name in events:
        if not _same(event_tree[name], events[name]):
            raise ValueError(f"eventTree scalar branch differs from events: {name}")
    if int(ak.sum(event_tree["nphotons"])) != root["photons"].num_entries:
        raise ValueError("eventTree photon total differs from photons")
    if int(ak.sum(event_tree["njets"])) != root["jets"].num_entries:
        raise ValueError("eventTree jet total differs from jets")
    if int(ak.sum(event_tree["npairs"])) != root["photonJets"].num_entries:
        raise ValueError("eventTree pair total differs from photonJets")

    photons = root["photons"].arrays(library="np")
    jets = root["jets"].arrays(library="np")
    pairs = root["photonJets"].arrays(library="np")
    photon_groups = _group_indices(photons)
    jet_groups = _group_indices(jets)
    pair_groups = _group_indices(pairs)
    photon_fields = {
        "photon_candidate_id_hi": "candidate_id_hi",
        "photon_candidate_id_lo": "candidate_id_lo",
        "photon_encounter_ordinal": "photon_encounter_ordinal",
        "photon_et": "photon_et",
        "photon_eta": "photon_eta",
        "photon_phi": "photon_phi",
        "photon_bdt_score": "bdt_score",
        "photon_bdt_tight_threshold": "bdt_tight_threshold",
        "photon_bdt_nontight_low_threshold": "bdt_nontight_low_threshold",
        "photon_bdt_nontight_high_threshold": "bdt_nontight_high_threshold",
        "photon_iso_r04": "iso_r04",
        "photon_iso_r04_threshold": "iso_r04_threshold",
        "photon_iso_r04_nonisolated_threshold": "iso_r04_nonisolated_threshold",
    }
    jet_fields = {
        "jet_id_hi": "jet_id_hi", "jet_id_lo": "jet_id_lo", "jet_pt": "jet_pt",
        "jet_raw_pt": "jet_raw_pt", "jet_eta": "jet_eta", "jet_phi": "jet_phi",
        "jet_radius": "jet_radius", "jet_mass": "jet_mass", "jet_area": "jet_area",
    }
    pair_fields = {
        "pair_photon_index": "photon_index", "pair_jet_index": "jet_index",
        "pair_delta_phi": "delta_phi", "pair_xjgamma": "xjgamma",
        "pair_recoil_state": "recoil_state",
    }
    for event_index in range(len(event_tree)):
        event = _key(event_tree["event_id_hi"][event_index], event_tree["event_id_lo"][event_index])
        p_indices = photon_groups.get(event, [])
        j_indices = jet_groups.get(event, [])
        q_indices = pair_groups.get(event, [])
        if int(event_tree["nphotons"][event_index]) != len(p_indices):
            raise ValueError(f"eventTree photon count differs at row {event_index}")
        if int(event_tree["njets"][event_index]) != len(j_indices):
            raise ValueError(f"eventTree jet count differs at row {event_index}")
        if int(event_tree["npairs"][event_index]) != len(q_indices):
            raise ValueError(f"eventTree pair count differs at row {event_index}")
        for event_name, flat_name in photon_fields.items():
            if not _same(event_tree[event_name][event_index], photons[flat_name][p_indices]):
                raise ValueError(f"eventTree {event_name} differs at row {event_index}")
        for event_name, flat_name in jet_fields.items():
            if not _same(event_tree[event_name][event_index], jets[flat_name][j_indices]):
                raise ValueError(f"eventTree {event_name} differs at row {event_index}")
        for event_name, flat_name in pair_fields.items():
            if not _same(event_tree[event_name][event_index], pairs[flat_name][q_indices]):
                raise ValueError(f"eventTree {event_name} differs at row {event_index}")

        photon_indices = np.asarray(event_tree["pair_photon_index"][event_index], dtype="int32")
        jet_indices = np.asarray(event_tree["pair_jet_index"][event_index], dtype="int32")
        if np.any(photon_indices < 0) or np.any(photon_indices >= len(p_indices)):
            raise ValueError(f"photon local index is invalid at eventTree row {event_index}")
        if np.any(jet_indices < 0) or np.any(jet_indices >= len(j_indices)):
            raise ValueError(f"jet local index is invalid at eventTree row {event_index}")
        observed_pairs = list(zip(photon_indices.tolist(), jet_indices.tolist()))
        expected_pairs = {
            (photon_index, jet_index)
            for photon_index in range(len(p_indices))
            for jet_index in range(len(j_indices))
        }
        if len(observed_pairs) != len(set(observed_pairs)):
            raise ValueError(
                f"photonJets repeats a local photon/jet combination at row {event_index}"
            )
        if set(observed_pairs) != expected_pairs:
            raise ValueError(
                f"photonJets does not cover the local photon/jet Cartesian product "
                f"at row {event_index}"
            )
        if len(photon_indices):
            pair_rows = np.asarray(q_indices, dtype="int64")
            event_photon_rows = np.asarray(p_indices, dtype="int64")
            event_jet_rows = np.asarray(j_indices, dtype="int64")
            selected_photon_rows = event_photon_rows[photon_indices]
            selected_jet_rows = event_jet_rows[jet_indices]
            for side, target, selected_rows, fields in (
                (
                    "photon",
                    photons,
                    selected_photon_rows,
                    sorted((set(photons) - set(events)).intersection(pairs)),
                ),
                (
                    "jet",
                    jets,
                    selected_jet_rows,
                    sorted((set(jets) - set(events)).intersection(pairs)),
                ),
            ):
                for field in fields:
                    if not _same(pairs[field][pair_rows], target[field][selected_rows]):
                        raise ValueError(
                            f"photonJets {side} field disagrees with its local index "
                            f"at eventTree row {event_index}: {field}"
                        )
            expected_xj = (
                np.asarray(event_tree["jet_pt"][event_index], dtype="float64")[jet_indices]
                / np.asarray(event_tree["photon_et"][event_index], dtype="float64")[photon_indices]
            )
            if not np.array_equal(
                expected_xj,
                np.asarray(event_tree["pair_xjgamma"][event_index], dtype="float64"),
            ):
                raise ValueError(f"xJgamma local-index relation differs at row {event_index}")

        records = [
            {
                "photon_et": float(event_tree["photon_et"][event_index][local]),
                "photon_eta": float(event_tree["photon_eta"][event_index][local]),
                "photon_encounter_ordinal": int(
                    event_tree["photon_encounter_ordinal"][event_index][local]
                ),
                "photon_bdt_score": float(event_tree["photon_bdt_score"][event_index][local]),
                "photon_bdt_tight_threshold": float(
                    event_tree["photon_bdt_tight_threshold"][event_index][local]
                ),
                "photon_bdt_nontight_low_threshold": float(
                    event_tree["photon_bdt_nontight_low_threshold"][event_index][local]
                ),
                "photon_bdt_nontight_high_threshold": float(
                    event_tree["photon_bdt_nontight_high_threshold"][event_index][local]
                ),
                "photon_iso_r04": float(event_tree["photon_iso_r04"][event_index][local]),
                "photon_iso_r04_threshold": float(
                    event_tree["photon_iso_r04_threshold"][event_index][local]
                ),
                "photon_iso_r04_nonisolated_threshold": float(
                    event_tree["photon_iso_r04_nonisolated_threshold"][event_index][local]
                ),
            }
            for local in range(len(p_indices))
        ]
        for region in "ABCD":
            definitions = ("bounded", "complement") if region in "CD" else ("bounded",)
            for definition in definitions:
                program = compile_recoil_selection(
                    RecoilSelection(region=region, non_tight_definition=definition)
                )
                expected = program.choose_leader(records)
                branch = program.leader_branch
                if int(event_tree[branch][event_index]) != expected:
                    raise ValueError(f"{branch} differs at eventTree row {event_index}")


def _local_targets(
    root: uproot.ReadOnlyDirectory,
    tree_name: str,
    id_names: tuple[str, str],
) -> dict[tuple[int, int], tuple[int, tuple[int, int], int]]:
    arrays = root[tree_name].arrays(
        ["source_file_index", "event_id_hi", "event_id_lo", *id_names], library="np"
    )
    next_index: dict[tuple[int, int], int] = defaultdict(int)
    targets: dict[tuple[int, int], tuple[int, tuple[int, int], int]] = {}
    for source, event_hi, event_lo, object_hi, object_lo in zip(
        arrays["source_file_index"], arrays["event_id_hi"], arrays["event_id_lo"],
        arrays[id_names[0]], arrays[id_names[1]],
    ):
        event = _key(event_hi, event_lo)
        identity = _key(object_hi, object_lo)
        if identity in targets:
            raise ValueError(f"duplicate {tree_name} identity: {identity}")
        targets[identity] = (int(source), event, next_index[event])
        next_index[event] += 1
    return targets


def _validate_truth_links(root: uproot.ReadOnlyDirectory) -> None:
    target_maps = {
        "reco": {
            1: _local_targets(root, "photons", IDENTITIES["photons"]),
            2: _local_targets(root, "jets", IDENTITIES["jets"]),
        },
        "truth": {
            1: _local_targets(root, "truthPhotons", IDENTITIES["truthPhotons"]),
            2: _local_targets(root, "truthJets", IDENTITIES["truthJets"]),
        },
    }
    links = root["recoTruthLinks"].arrays(library="np")
    for row in range(len(links["link_id_hi"])):
        source = int(links["source_file_index"][row])
        event = _key(links["event_id_hi"][row], links["event_id_lo"][row])
        target_types: dict[str, int] = {}
        for side in ("reco", "truth"):
            target_type = int(links[f"{side}_type"][row])
            target_index = int(links[f"{side}_index"][row])
            target_identity = _key(
                links[f"{side}_id_hi"][row], links[f"{side}_id_lo"][row]
            )
            target_types[side] = target_type
            if target_type == 0:
                if target_index != -1 or target_identity != (0, 0):
                    raise ValueError(
                        f"absent {side} target is not the canonical null at link row {row}"
                    )
                continue
            if target_type not in target_maps[side]:
                raise ValueError(f"unsupported {side} target type at link row {row}")
            observed = target_maps[side][target_type].get(target_identity)
            if observed != (source, event, target_index):
                raise ValueError(f"wrong {side} target event/index at link row {row}")

        link_class = int(links["link_class"][row])
        if link_class not in _LINK_CLASSES:
            raise ValueError(f"unsupported link_class at link row {row}")
        reco_present = target_types["reco"] != 0
        truth_present = target_types["truth"] != 0
        if link_class == _LINK_RECO_FAKE:
            valid_topology = reco_present and not truth_present
        elif link_class == _LINK_TRUTH_MISS:
            valid_topology = not reco_present and truth_present
        else:
            valid_topology = reco_present and truth_present
        if not valid_topology:
            raise ValueError(f"link_class topology differs at link row {row}")


def _validate_truth_occurrence_identity(root: uproot.ReadOnlyDirectory) -> None:
    """Bind every barcode convenience field to an occurrence-safe truth row.

    HepMC barcodes can repeat across generator occurrences in the same event.
    A barcode is therefore never interpreted alone: the truth-row key is the
    pair ``(generator_occurrence_embedding_id, generator_barcode)`` and a
    reconstructed match must resolve through the typed reco-truth link to that
    exact truth row.
    """

    truth = root["truthPhotons"].arrays(
        [
            "event_id_hi",
            "event_id_lo",
            "truth_photon_id_hi",
            "truth_photon_id_lo",
            "generator_barcode",
            "generator_occurrence_embedding_id",
        ],
        library="np",
    )
    truth_by_event_and_generator: dict[
        tuple[tuple[int, int], tuple[int, int]], tuple[int, int]
    ] = {}
    for row in range(len(truth["generator_barcode"])):
        event = _key(truth["event_id_hi"][row], truth["event_id_lo"][row])
        barcode = int(truth["generator_barcode"][row])
        occurrence = int(truth["generator_occurrence_embedding_id"][row])
        if barcode < 0:
            raise ValueError(f"truthPhotons row {row} has a negative generator barcode")
        generator_key = (occurrence, barcode)
        composite = (event, generator_key)
        if composite in truth_by_event_and_generator:
            raise ValueError(
                "duplicate generator occurrence/barcode identity in truthPhotons "
                f"for event {event}: {generator_key}"
            )
        truth_by_event_and_generator[composite] = _key(
            truth["truth_photon_id_hi"][row], truth["truth_photon_id_lo"][row]
        )

    links = root["recoTruthLinks"].arrays(
        [
            "event_id_hi",
            "event_id_lo",
            "reco_type",
            "reco_id_hi",
            "reco_id_lo",
            "truth_type",
            "truth_id_hi",
            "truth_id_lo",
            "link_class",
        ],
        library="np",
    )
    linked_truth_by_photon: dict[
        tuple[tuple[int, int], tuple[int, int]], tuple[int, int]
    ] = {}
    for row in range(len(links["reco_type"])):
        if (
            int(links["reco_type"][row]) != 1
            or int(links["link_class"][row]) != _LINK_MATCH
        ):
            continue
        event = _key(links["event_id_hi"][row], links["event_id_lo"][row])
        reco_identity = _key(links["reco_id_hi"][row], links["reco_id_lo"][row])
        link_key = (event, reco_identity)
        truth_type = int(links["truth_type"][row])
        truth_identity = (
            _key(links["truth_id_hi"][row], links["truth_id_lo"][row])
            if truth_type == 1
            else (0, 0)
        )
        previous = linked_truth_by_photon.get(link_key)
        if previous is not None and previous != truth_identity:
            raise ValueError(f"photon has conflicting typed truth links: {reco_identity}")
        linked_truth_by_photon[link_key] = truth_identity

    photons = root["photons"].arrays(
        [
            "event_id_hi",
            "event_id_lo",
            "candidate_id_hi",
            "candidate_id_lo",
            "truth_matched",
            "truth_barcode",
            "truth_generator_occurrence_embedding_id",
        ],
        library="np",
    )
    photon_truth_state: dict[
        tuple[tuple[int, int], tuple[int, int]], tuple[int, int, int]
    ] = {}
    for row in range(len(photons["truth_matched"])):
        event = _key(photons["event_id_hi"][row], photons["event_id_lo"][row])
        photon_identity = _key(
            photons["candidate_id_hi"][row], photons["candidate_id_lo"][row]
        )
        matched = int(photons["truth_matched"][row])
        barcode = int(photons["truth_barcode"][row])
        occurrence = int(photons["truth_generator_occurrence_embedding_id"][row])
        if matched not in (0, 1):
            raise ValueError(f"photons row {row} truth_matched is not binary")
        if not matched:
            if (occurrence, barcode) != (-1, -1):
                raise ValueError(
                    f"unmatched photons row {row} must use the (-1, -1) truth sentinel"
                )
        else:
            expected_truth = truth_by_event_and_generator.get(
                (event, (occurrence, barcode))
            )
            if expected_truth is None:
                raise ValueError(
                    f"photons row {row} occurrence/barcode has no truthPhoton in its event"
                )
            linked_truth = linked_truth_by_photon.get((event, photon_identity))
            if linked_truth != expected_truth:
                raise ValueError(
                    f"photons row {row} occurrence/barcode disagrees with its typed truth link"
                )
        photon_truth_state[(event, photon_identity)] = (matched, occurrence, barcode)

    pairs = root["photonJets"].arrays(
        [
            "event_id_hi",
            "event_id_lo",
            "candidate_id_hi",
            "candidate_id_lo",
            "truth_matched",
            "truth_barcode",
            "truth_generator_occurrence_embedding_id",
        ],
        library="np",
    )
    for row in range(len(pairs["truth_matched"])):
        event = _key(pairs["event_id_hi"][row], pairs["event_id_lo"][row])
        photon_identity = _key(
            pairs["candidate_id_hi"][row], pairs["candidate_id_lo"][row]
        )
        observed = (
            int(pairs["truth_matched"][row]),
            int(pairs["truth_generator_occurrence_embedding_id"][row]),
            int(pairs["truth_barcode"][row]),
        )
        if photon_truth_state.get((event, photon_identity)) != observed:
            raise ValueError(
                f"photonJets row {row} truth identity differs from its photon row"
            )


def _validate_numerics(
    root: uproot.ReadOnlyDirectory,
    model_input_count: int | None,
    require_scaled_bit30: bool,
) -> None:
    event_arrays = root["events"].arrays(
        ["event_weight", "vertex_z", "centrality", "scaled_trigger_bits", "scaled_bit30"],
        library="np",
    )
    if not np.all(np.isfinite(event_arrays["event_weight"])):
        raise ValueError("event weights are not finite")
    if not np.all(np.isfinite(event_arrays["vertex_z"])):
        raise ValueError("event vertices are not finite")
    centrality = np.asarray(event_arrays["centrality"], dtype=float)
    pp_not_applicable = bool(np.all(centrality == -1.0))
    auau_percent = bool(
        np.all(np.isfinite(centrality))
        and np.all(centrality >= 0.0)
        and np.all(centrality < 100.0)
    )
    if not (pp_not_applicable or auau_percent):
        raise ValueError(
            "event centrality must be uniformly -1 for p+p or uniformly within [0, 100) for Au+Au"
        )
    if require_scaled_bit30:
        reconstructed = (
            (event_arrays["scaled_trigger_bits"] >> np.uint64(30)) & np.uint64(1)
        ).astype("int32")
        if not np.all(reconstructed == 1) or not np.all(event_arrays["scaled_bit30"] == 1):
            raise ValueError("scaled-bit-30 membership failed")

    fields = ["bdt_score", "iso_r03", "iso_r04"]
    if model_input_count is not None:
        if model_input_count not in (11, 14):
            raise ValueError("model_input_count must be 11 or 14")
        fields.extend(
            ["bdt_input_count", *[f"bdt_input_{index:02d}" for index in range(model_input_count)]]
        )
    arrays = root["photons"].arrays(fields, library="np")
    for field in fields:
        if field == "bdt_input_count":
            if not np.all(arrays[field] == model_input_count):
                raise ValueError("BDT input dimension differs")
        elif not np.all(np.isfinite(arrays[field])):
            raise ValueError(f"non-finite photon field: {field}")
    radii = root["jets"]["jet_radius"].array(library="np")
    if len(radii) and not np.allclose(radii, 0.4, rtol=0.0, atol=1.0e-12):
        raise ValueError("reconstructed jet radius differs from R=0.4")


def validate(
    input_paths: Iterable[Path],
    *,
    model_input_count: int | None = None,
    require_scaled_bit30: bool = False,
) -> list[str]:
    """Validate one or more public tree parts and return receipt-ready checks."""

    paths = [Path(path).resolve() for path in input_paths]
    if not paths:
        raise ValueError("at least one PhotonJetTrees_v1 input is required")
    reference_schema: dict[str, list[tuple[str, str]]] | None = None
    global_identities = {name: set() for name in IDENTITIES}
    totals = {name: 0 for name in TREE_ORDER}
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        with uproot.open(path) as root:
            _validate_schema(root)
            schema = {name: list(root[name].typenames().items()) for name in TREE_ORDER}
            if reference_schema is None:
                reference_schema = schema
            elif schema != reference_schema:
                raise ValueError(f"branch schema differs across parts: {path}")
            for name in IDENTITIES:
                current = _identity_set(root, name)
                overlap = global_identities[name].intersection(current)
                if overlap:
                    raise ValueError(f"identity repeats across parts in {name}: {next(iter(overlap))}")
                global_identities[name].update(current)
            _validate_parent_events(root)
            _validate_event_tree(root)
            _validate_truth_links(root)
            _validate_truth_occurrence_identity(root)
            _validate_numerics(root, model_input_count, require_scaled_bit30)
            for name in TREE_ORDER:
                totals[name] += int(root[name].num_entries)
    return [
        "tree_schema_exact=PASS",
        "identity_uniqueness=PASS",
        "event_parentage=PASS",
        "eventTree_flat_equivalence=PASS",
        "eventTree_executable_leaders=PASS",
        "truth_link_type_event_index=PASS",
        "truth_occurrence_identity=PASS",
        "numeric_contract=PASS",
        "counts=" + ",".join(f"{name}:{totals[name]}" for name in TREE_ORDER),
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--model-input-count", type=int, choices=(11, 14))
    parser.add_argument("--require-scaled-bit30", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    lines = validate(
        args.input,
        model_input_count=args.model_input_count,
        require_scaled_bit30=args.require_scaled_bit30,
    )
    text = "VALIDATION=PASS\n" + "\n".join(lines) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
