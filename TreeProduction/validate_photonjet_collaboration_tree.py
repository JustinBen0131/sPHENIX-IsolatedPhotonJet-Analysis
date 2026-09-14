#!/usr/bin/env python3
"""Validate a collaborator-facing photon+jet ROOT against its source files."""

from __future__ import annotations

import argparse
import gc
import json
import math
from pathlib import Path
import resource
import sys
import tempfile
from typing import Any

import awkward as ak
import numpy as np
import uproot
import trigger_scaler_interface
import gl1_trigger_interface


PREFIX = "ReplayFoundationV1/"
MODEL_ET_MIN_GEV = 15.0
MODEL_ET_MAX_GEV = 35.0
SOURCE_TREES = {
    "events": "RJEventV1",
    "photons": "RJPhotonCandidateV1",
    "jets": "RJJetV1",
    "photonJets": "RJPhotonJetPairV1",
    "truthPhotons": "RJTruthPhotonV1",
    "truthJets": "RJTruthJetV1",
    "recoTruthLinks": "RJRecoTruthLinkV1",
}
IDENTITY_BRANCHES = {
    "events": ("event_id_hi", "event_id_lo"),
    "photons": ("candidate_id_hi", "candidate_id_lo"),
    "jets": ("jet_id_hi", "jet_id_lo"),
    "photonJets": ("pair_id_hi", "pair_id_lo"),
    "truthPhotons": ("truth_photon_id_hi", "truth_photon_id_lo"),
    "truthJets": ("truth_jet_id_hi", "truth_jet_id_lo"),
    "recoTruthLinks": ("link_id_hi", "link_id_lo"),
}
VALIDATION_STEP_SIZE = 100_000
EVENT_VALIDATION_STEP_SIZE = 256
SCOPED_IDENTITY_DTYPE = np.dtype(
    [("source", ">u4"), ("hi", ">u8"), ("lo", ">u8")]
)
IDENTITY_PROOF_MAX_RECORDS = 250_000
IDENTITY_PROOF_FANOUT = 16


def _ids(source_index: np.ndarray, hi: np.ndarray, lo: np.ndarray) -> set[tuple[int, int, int]]:
    return {(int(source), int(high), int(low)) for source, high, low in zip(source_index, hi, lo)}


def _key(source: Any, hi: Any, lo: Any) -> tuple[int, int, int]:
    return int(source), int(hi), int(lo)


def _validate_truth_link_targets(root: Any, tree_name: str = "recoTruthLinks") -> None:
    """Require the exact typed target, event, and event-local encounter index."""
    if not root[tree_name].num_entries:
        return
    targets = {}
    event_fields = ("source_file_index", "event_id_hi", "event_id_lo")
    for side, names in (("reco", ("photons", "jets")),
                        ("truth", ("truthPhotons", "truthJets"))):
        for kind, name in enumerate(names, 1):
            per_event = {}
            lookup = targets[side, kind] = {}
            fields = (*event_fields, *IDENTITY_BRANCHES[name])
            for chunk in root[name].iterate(list(fields), step_size=VALIDATION_STEP_SIZE, library="np"):
                for source, event_hi, event_lo, high, low in zip(*(chunk[key] for key in fields)):
                    event = _key(source, event_hi, event_lo)
                    identity = _key(source, high, low)
                    if identity in lookup:
                        raise ValueError(f"recoTruthLinks duplicate typed target in {name}")
                    index = per_event.get(event, 0)
                    lookup[identity] = event, index
                    per_event[event] = index + 1
    fields = [*event_fields, *(f"{side}_{field}" for side in ("reco", "truth")
                              for field in ("type", "id_hi", "id_lo", "index"))]
    for links in root[tree_name].iterate(fields, step_size=VALIDATION_STEP_SIZE, library="np"):
        for values in zip(*(links[key] for key in fields)):
            source, event_hi, event_lo, *sides = map(int, values)
            event = _key(source, event_hi, event_lo)
            for side, values in (("reco", sides[:4]), ("truth", sides[4:])):
                kind, high, low, index = values
                if kind == 0:
                    if index != -1:
                        raise ValueError(f"recoTruthLinks {side} null local-index state is inconsistent")
                    continue
                if (side, kind) not in targets:
                    raise ValueError(f"recoTruthLinks unknown {side} target type")
                expected = targets[side, kind].get(_key(source, high, low))
                if expected is None:
                    raise ValueError(f"recoTruthLinks missing typed {side} target")
                if expected != (event, index):
                    raise ValueError(f"recoTruthLinks {side} target event or local index differs")


def _rss_mb() -> float:
    maximum = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return maximum / (1024.0 * 1024.0) if sys.platform == "darwin" else maximum / 1024.0


def _report_memory(stage: str) -> None:
    print(f"VALIDATOR_RSS_MB stage={stage} peak={_rss_mb():.1f}", flush=True)


def _scoped_identity_records(
    source: np.ndarray, high: np.ndarray, low: np.ndarray
) -> np.ndarray:
    records = np.empty(len(high), dtype=SCOPED_IDENTITY_DTYPE)
    records["source"] = np.asarray(source, dtype=">u4")
    records["hi"] = np.asarray(high, dtype=">u8")
    records["lo"] = np.asarray(low, dtype=">u8")
    return records


def _identity_nibbles(records: np.ndarray, depth: int) -> np.ndarray:
    raw = records.view(np.uint8).reshape(len(records), SCOPED_IDENTITY_DTYPE.itemsize)
    byte_values = raw[:, depth // 2]
    return (byte_values >> 4) if depth % 2 == 0 else (byte_values & 0x0F)


def _partition_identity_records(
    records: np.ndarray, depth: int, streams: dict[int, Any], paths: dict[int, Path]
) -> None:
    if not len(records):
        return
    buckets = _identity_nibbles(records, depth)
    order = np.argsort(buckets, kind="stable")
    ordered = buckets[order]
    boundaries = np.flatnonzero(np.diff(ordered)) + 1
    starts = np.concatenate((np.asarray([0]), boundaries))
    stops = np.concatenate((boundaries, np.asarray([len(order)])))
    for start, stop in zip(starts, stops):
        bucket = int(ordered[start])
        stream = streams.get(bucket)
        if stream is None:
            stream = paths[bucket].open("xb")
            streams[bucket] = stream
        stream.write(records[order[start:stop]].tobytes())


def _prove_unique_identity_file(
    path: Path,
    label: str,
    depth: int = 0,
    name_counter: list[int] | None = None,
) -> None:
    if name_counter is None:
        name_counter = [0]
    size = path.stat().st_size
    record_bytes = SCOPED_IDENTITY_DTYPE.itemsize
    if size % record_bytes:
        raise ValueError(f"identity-proof spool is truncated for {label}")
    count = size // record_bytes
    if count <= IDENTITY_PROOF_MAX_RECORDS:
        records = np.fromfile(path, dtype=SCOPED_IDENTITY_DTYPE)
        if len(records) > 1:
            records.sort(order=("source", "hi", "lo"))
            repeated = (
                (records["source"][1:] == records["source"][:-1])
                & (records["hi"][1:] == records["hi"][:-1])
                & (records["lo"][1:] == records["lo"][:-1])
            )
            indices = np.flatnonzero(repeated)
            if len(indices):
                index = int(indices[0])
                duplicate = (
                    int(records["source"][index]),
                    int(records["hi"][index]),
                    int(records["lo"][index]),
                )
                path.unlink()
                raise ValueError(f"duplicate source-scoped identity in {label}: {duplicate}")
        path.unlink()
        return
    if depth >= 2 * record_bytes:
        raise ValueError(f"duplicate full source-scoped identity in {label}")

    child_paths: dict[int, Path] = {}
    for bucket in range(IDENTITY_PROOF_FANOUT):
        name_counter[0] += 1
        child = path.parent / f".identity_{name_counter[0]:08x}.bin"
        if child.exists():
            raise FileExistsError(f"identity-proof child already exists: {child}")
        child_paths[bucket] = child
    streams: dict[int, Any] = {}
    try:
        with path.open("rb") as source:
            while True:
                records = np.fromfile(
                    source,
                    dtype=SCOPED_IDENTITY_DTYPE,
                    count=IDENTITY_PROOF_MAX_RECORDS,
                )
                if not len(records):
                    break
                _partition_identity_records(records, depth, streams, child_paths)
    finally:
        for stream in streams.values():
            stream.close()
    path.unlink()
    for child in child_paths.values():
        if child.exists():
            _prove_unique_identity_file(child, label, depth + 1, name_counter)


def _retention_mask(
    arrays: dict[str, np.ndarray],
    source_index: int,
    retained_photons: set[tuple[int, int, int]] | None,
    *,
    link_tree: bool = False,
) -> np.ndarray:
    if retained_photons is None:
        first = next(iter(arrays.values()))
        return np.ones(len(first), dtype=np.bool_)
    if link_tree:
        return np.fromiter(
            (
                int(reco_type) != 1
                or (source_index, int(high), int(low)) in retained_photons
                for reco_type, high, low in zip(
                    arrays["reco_type"], arrays["reco_id_hi"], arrays["reco_id_lo"]
                )
            ),
            dtype=np.bool_,
            count=len(arrays["reco_type"]),
        )
    return np.fromiter(
        (
            (source_index, int(high), int(low)) in retained_photons
            for high, low in zip(arrays["candidate_id_hi"], arrays["candidate_id_lo"])
        ),
        dtype=np.bool_,
        count=len(arrays["candidate_id_hi"]),
    )


def _iter_source_identities(
    inputs: list[Path],
    output_name: str,
    retained_photons: set[tuple[int, int, int]] | None,
):
    source_name = SOURCE_TREES[output_name]
    hi_name, lo_name = IDENTITY_BRANCHES[output_name]
    dependent = output_name in {"photons", "photonJets", "recoTruthLinks"}
    link_tree = output_name == "recoTruthLinks"
    for source_index, source in enumerate(inputs):
        branches = [hi_name, lo_name]
        if dependent and retained_photons is not None:
            branches += (
                ["reco_type", "reco_id_hi", "reco_id_lo"]
                if link_tree
                else ["candidate_id_hi", "candidate_id_lo"]
            )
        with uproot.open(source) as root:
            tree = root[PREFIX + source_name]
            for arrays in tree.iterate(
                branches, step_size=VALIDATION_STEP_SIZE, library="np"
            ):
                mask = _retention_mask(
                    arrays,
                    source_index,
                    retained_photons if dependent else None,
                    link_tree=link_tree,
                )
                count = int(np.count_nonzero(mask))
                if count:
                    yield (
                        np.full(count, source_index, dtype=np.int32),
                        np.asarray(arrays[hi_name][mask], dtype=np.uint64),
                        np.asarray(arrays[lo_name][mask], dtype=np.uint64),
                    )


def _iter_output_identities(root: uproot.ReadOnlyDirectory, output_name: str):
    hi_name, lo_name = IDENTITY_BRANCHES[output_name]
    for arrays in root[output_name].iterate(
        ["source_file_index", hi_name, lo_name],
        step_size=VALIDATION_STEP_SIZE,
        library="np",
    ):
        yield (
            np.asarray(arrays["source_file_index"], dtype=np.int32),
            np.asarray(arrays[hi_name], dtype=np.uint64),
            np.asarray(arrays[lo_name], dtype=np.uint64),
        )


def _compare_identity_streams(
    expected,
    observed,
    proof_path: Path,
    label: str,
) -> int:
    expected_iter = iter(expected)
    observed_iter = iter(observed)
    expected_chunk = observed_chunk = None
    expected_offset = observed_offset = 0
    total = 0
    with proof_path.open("xb") as proof:
        while True:
            if expected_chunk is None or expected_offset == len(expected_chunk[0]):
                expected_chunk = next(expected_iter, None)
                expected_offset = 0
            if observed_chunk is None or observed_offset == len(observed_chunk[0]):
                observed_chunk = next(observed_iter, None)
                observed_offset = 0
            if expected_chunk is None or observed_chunk is None:
                if expected_chunk is not None or observed_chunk is not None:
                    raise ValueError(f"{label} identity stream length differs")
                break
            count = min(
                len(expected_chunk[0]) - expected_offset,
                len(observed_chunk[0]) - observed_offset,
            )
            expected_slice = tuple(
                values[expected_offset : expected_offset + count]
                for values in expected_chunk
            )
            observed_slice = tuple(
                values[observed_offset : observed_offset + count]
                for values in observed_chunk
            )
            for field, expected_values, observed_values in zip(
                ("source_file_index", "identity_hi", "identity_lo"),
                expected_slice,
                observed_slice,
            ):
                mismatch = np.flatnonzero(expected_values != observed_values)
                if len(mismatch):
                    row = total + int(mismatch[0])
                    raise ValueError(
                        f"{label} ordered identity mismatch at row {row} field {field}: "
                        f"{int(observed_values[mismatch[0]])} != {int(expected_values[mismatch[0]])}"
                    )
            proof.write(_scoped_identity_records(*expected_slice).tobytes())
            total += count
            expected_offset += count
            observed_offset += count
    _prove_unique_identity_file(proof_path, label)
    return total


def _iter_source_pair_observable(
    inputs: list[Path],
    retained_photons: set[tuple[int, int, int]] | None,
    branch: str,
):
    if branch not in {"delta_phi", "xjgamma"}:
        raise ValueError(f"unsupported source pair observable: {branch}")
    for source_index, source in enumerate(inputs):
        with uproot.open(source) as root:
            tree = root[PREFIX + SOURCE_TREES["photonJets"]]
            for arrays in tree.iterate(
                ["candidate_id_hi", "candidate_id_lo", branch],
                step_size=VALIDATION_STEP_SIZE,
                library="np",
            ):
                mask = _retention_mask(arrays, source_index, retained_photons)
                if np.any(mask):
                    yield np.asarray(arrays[branch][mask], dtype=np.float64)


def _wrapped_abs_delta_phi(phi_a: np.ndarray, phi_b: np.ndarray) -> np.ndarray:
    delta = np.asarray(phi_a, dtype=np.float64) - np.asarray(phi_b, dtype=np.float64)
    return np.abs(np.arctan2(np.sin(delta), np.cos(delta)))


def _compare_float_streams(expected, observed, label: str) -> int:
    expected_iter = iter(expected)
    observed_iter = iter(observed)
    expected_chunk = observed_chunk = None
    expected_offset = observed_offset = 0
    total = 0
    while True:
        if expected_chunk is None or expected_offset == len(expected_chunk):
            expected_chunk = next(expected_iter, None)
            expected_offset = 0
        if observed_chunk is None or observed_offset == len(observed_chunk):
            observed_chunk = next(observed_iter, None)
            observed_offset = 0
        if expected_chunk is None or observed_chunk is None:
            if expected_chunk is not None or observed_chunk is not None:
                raise ValueError(f"{label} stream length differs")
            break
        count = min(
            len(expected_chunk) - expected_offset,
            len(observed_chunk) - observed_offset,
        )
        expected_values = expected_chunk[expected_offset : expected_offset + count]
        observed_values = observed_chunk[observed_offset : observed_offset + count]
        equal = (expected_values == observed_values) | (
            np.isnan(expected_values) & np.isnan(observed_values)
        )
        mismatch = np.flatnonzero(~equal)
        if len(mismatch):
            row = total + int(mismatch[0])
            raise ValueError(f"{label} mismatch at row {row}")
        total += count
        expected_offset += count
        observed_offset += count
    return total


def _source_identity_set(inputs: list[Path], tree_name: str, branches: tuple[str, str]) -> set[tuple[int, int, int]]:
    identities: set[tuple[int, int, int]] = set()
    for source_index, source in enumerate(inputs):
        with uproot.open(source) as root:
            arrays = root[PREFIX + tree_name].arrays(list(branches), library="np")
            current = _ids(
                np.full(len(arrays[branches[0]]), source_index, dtype="int32"),
                arrays[branches[0]],
                arrays[branches[1]],
            )
            if len(current) != len(arrays[branches[0]]):
                raise ValueError(f"duplicate source identity in {source}: {tree_name}")
            overlap = identities.intersection(current)
            if overlap:
                raise ValueError(f"duplicate source-scoped identity across inputs: {tree_name} {next(iter(overlap))}")
            identities.update(current)
    return identities


def _retained_photon_identities(
    inputs: list[Path], input_count: int | None
) -> set[tuple[int, int, int]] | None:
    if input_count is None:
        return None
    identities: set[tuple[int, int, int]] = set()
    for source_index, source in enumerate(inputs):
        with uproot.open(source) as root:
            candidates = root[PREFIX + "RJPhotonCandidateV1"].arrays(
                ["candidate_id_hi", "candidate_id_lo", "cluster_et"], library="np"
            )
            views = root[PREFIX + "RJShowerFeatureViewV1"].arrays(
                [
                    "candidate_id_hi", "candidate_id_lo", "definition_name",
                    "finite_feature_state", "ordered_features",
                ],
                library="ak",
                how=dict,
            )
        h70: dict[tuple[int, int], tuple[np.ndarray, int]] = {}
        for index in range(len(views["candidate_id_hi"])):
            if str(views["definition_name"][index]).upper() != "H70":
                continue
            candidate_key = (
                int(views["candidate_id_hi"][index]), int(views["candidate_id_lo"][index])
            )
            if candidate_key in h70:
                raise ValueError(f"duplicate H70 shower-feature view in {source}: {candidate_key}")
            h70[candidate_key] = (
                np.asarray(views["ordered_features"][index], dtype="float64"),
                int(views["finite_feature_state"][index]),
            )
        for high, low, photon_et in zip(
            candidates["candidate_id_hi"], candidates["candidate_id_lo"], candidates["cluster_et"]
        ):
            candidate_key = (int(high), int(low))
            view = h70.get(candidate_key)
            exact_inputs_valid = (
                view is not None
                and view[1] == 1
                and len(view[0]) == input_count
                and bool(np.all(np.isfinite(view[0])))
            )
            if not exact_inputs_valid:
                if MODEL_ET_MIN_GEV <= float(photon_et) < MODEL_ET_MAX_GEV:
                    raise ValueError(
                        f"candidate {candidate_key} in {source} lacks valid exact H70 inputs "
                        "inside the accepted 15--35 GeV model window"
                    )
                continue
            identity = (source_index, *candidate_key)
            if identity in identities:
                raise ValueError(f"duplicate source-scoped photon identity: {identity}")
            identities.add(identity)
    return identities


def _source_dependent_identity_set(
    inputs: list[Path],
    tree_name: str,
    identity_branches: tuple[str, str],
    retained_photons: set[tuple[int, int, int]],
    *,
    link_tree: bool = False,
) -> set[tuple[int, int, int]]:
    identities: set[tuple[int, int, int]] = set()
    for source_index, source in enumerate(inputs):
        dependent_branches = (
            ["reco_type", "reco_id_hi", "reco_id_lo"]
            if link_tree
            else ["candidate_id_hi", "candidate_id_lo"]
        )
        branches = list(identity_branches) + dependent_branches
        with uproot.open(source) as root:
            arrays = root[PREFIX + tree_name].arrays(branches, library="np")
        for index in range(len(arrays[identity_branches[0]])):
            if link_tree:
                retained = int(arrays["reco_type"][index]) != 1 or (
                    source_index,
                    int(arrays["reco_id_hi"][index]),
                    int(arrays["reco_id_lo"][index]),
                ) in retained_photons
            else:
                retained = (
                    source_index,
                    int(arrays["candidate_id_hi"][index]),
                    int(arrays["candidate_id_lo"][index]),
                ) in retained_photons
            if not retained:
                continue
            identity = (
                source_index,
                int(arrays[identity_branches[0]][index]),
                int(arrays[identity_branches[1]][index]),
            )
            if identity in identities:
                raise ValueError(f"duplicate retained source identity: {tree_name} {identity}")
            identities.add(identity)
    return identities


def _make_evaluator(model_path: Path | None, input_count: int | None):
    if model_path is None:
        return None
    if input_count not in (11, 14):
        raise ValueError("model validation requires --model-input-count 11 or 14")
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError("PyROOT is required for model validation") from exc
    model = ROOT.TMVA.Experimental.RBDT("myBDT", str(model_path.resolve()))

    def evaluate(values: np.ndarray) -> float:
        if len(values) != input_count:
            raise ValueError(f"BDT input dimension mismatch: {len(values)} != {input_count}")
        vector = ROOT.std.vector("float")()
        for value in values:
            vector.push_back(float(value))
        result = model.Compute(vector)
        if len(result) != 1:
            raise ValueError(f"BDT returned {len(result)} scores")
        return float(result[0])

    return evaluate


def _leader_index(
    scores: np.ndarray,
    tight_thresholds: np.ndarray,
    nontight_low: np.ndarray,
    nontight_high: np.ndarray,
    isolation: np.ndarray,
    isolated_thresholds: np.ndarray,
    nonisolated_thresholds: np.ndarray,
    photon_et: np.ndarray,
    encounter_ordinal: np.ndarray,
    region: str,
    *,
    complement: bool = False,
) -> int:
    tight = scores > tight_thresholds
    nontight = (~tight) if complement else (scores > nontight_low) & (scores < nontight_high)
    isolated = isolation < isolated_thresholds
    nonisolated = isolation > nonisolated_thresholds
    finite = (
        np.isfinite(scores) & np.isfinite(tight_thresholds)
        & np.isfinite(nontight_low) & np.isfinite(nontight_high) & np.isfinite(isolation)
        & np.isfinite(isolated_thresholds) & np.isfinite(nonisolated_thresholds)
    )
    mask = finite & {
        "A": tight & isolated,
        "B": tight & nonisolated,
        "C": nontight & isolated,
        "D": nontight & nonisolated,
    }[region]
    eligible = np.flatnonzero(mask)
    if len(eligible) == 0:
        return -1
    return int(min(eligible, key=lambda index: (-float(photon_et[index]), int(encounter_ordinal[index]), int(index))))


def _expected_thresholds(system: str, photon_et: np.ndarray, centrality: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if system == "pp":
        return (
            0.815625 - 0.0015625 * photon_et,
            0.7333333333333333 - 0.01333333333333333 * photon_et,
            0.684375 + 0.0015625 * photon_et,
        )
    if system == "auau":
        return (
            0.6927385742113304 + 0.0010476808935335573 * centrality,
            0.4634577166843198 + 0.0015466253667456042 * centrality,
            0.6071724560222583 + 0.0013019695025559648 * centrality,
        )
    raise ValueError(f"unsupported collision system: {system!r}")


def validate_interface(inputs, root, system, *, require_complete=False):
    """Compare every added field to its exact source and repeated event join."""
    import build_photonjet_collaboration_tree as builder

    if "collaborator_interface_version" not in root["events"].keys():
        if require_complete:
            raise ValueError("complete interface unavailable: legacy output")
        return "collaborator_interface=LEGACY_V1"
    keys = ["source_file_index", "event_id_hi", "event_id_lo"]
    branches = [*keys, *builder.INTERFACE_SCHEMA]
    actual = root["events"].arrays(branches, library="np")
    for source_index, source in enumerate(inputs):
        with uproot.open(source) as original:
            events = builder._arrays(original[PREFIX + "RJEventV1"], builder.EVENT_INPUT_BRANCHES)
            weights = (builder._arrays(original[PREFIX + "RJWeightComponentV1"], builder.WEIGHT_INPUT_BRANCHES)
                       if PREFIX + "RJWeightComponentV1" in original else None)
            expected = builder._interface_columns(events, weights, system,
                                                   require_complete=require_complete)
            selected = actual["source_file_index"] == source_index
            for name, values in expected.items():
                if not np.array_equal(actual[name][selected], values, equal_nan=True):
                    raise ValueError(f"interface source readback mismatch: {source_index}/{name}")
            expected_mbd = builder._mbd_columns(events)[1]
            expected_mbd.update(builder.centrality_replay.capture_columns(events)[1])
            expected_mbd.update(builder._truth_vertex_columns(events)[1])
            expected_mbd.update(builder._truth_jet_availability(events)[1])
            event_tree = root["eventTree" if "eventTree" in root else "events"].arrays(
                ["source_file_index", *builder.MBD_ARRAY_SCHEMA, *builder.centrality_replay.ARRAYS,
                 *builder.TRUTH_VERTEX_ARRAY_SCHEMA, *builder.TRUTH_JET_AVAILABILITY_SCHEMA], library="ak")
            mbd_mask = event_tree["source_file_index"] == source_index
            for name, values in expected_mbd.items():
                observed = event_tree[name][mbd_mask]
                if ak.to_list(ak.num(observed)) != ak.to_list(ak.num(values)) or not np.array_equal(
                    ak.to_numpy(ak.flatten(observed)), ak.to_numpy(ak.flatten(values)), equal_nan=True
                ):
                    raise ValueError(f"MBD per-PMT source readback mismatch: {source_index}/{name}")
    lookup = {_key(*identity): i for i, identity in enumerate(zip(*(actual[k] for k in keys)))}
    if len(lookup) != len(actual[keys[0]]):
        raise ValueError("duplicate interface event identity")
    for tree_name in ("eventTree", "photons", "jets", "photonJets"):
        if "storage_layout" in root and str(root["storage_layout"]) == "normalized_v1":
            # Normalized objects join event state by identity; it is not copied.
            continue
        for rows in root[tree_name].iterate(branches, library="np", step_size=VALIDATION_STEP_SIZE):
            indices = np.asarray([lookup[_key(*identity)] for identity in
                                  zip(*(rows[k] for k in keys))], dtype="int64")
            for name in builder.INTERFACE_SCHEMA:
                if not np.array_equal(rows[name], actual[name][indices], equal_nan=True):
                    raise ValueError(f"interface event join mismatch: {tree_name}/{name}")
    return "collaborator_interface_v2_source_and_event_joins=PASS"


def validate_dominant_truth_witnesses(inputs, root, *, require_complete=False,
                                     retained_photons=None):
    """Source-qualified transport proof, not acceptance of a training population."""
    defaults = {
        "dominant_truth_state": -2,
        "dominant_truth_evaluator_mode": 0,
        "dominant_truth_track_id": -1,
        "dominant_truth_pid": 0,
        "dominant_truth_barcode": -1,
        "dominant_truth_embedding_id": 0,
        "dominant_truth_energy_contribution": np.nan,
        "dominant_truth_vertex_id": -1,
    }
    fields = tuple(defaults)
    identity = ("candidate_id_hi", "candidate_id_lo", "event_id_hi", "event_id_lo")
    output = root["photons"]
    if not set(fields).issubset(output.keys()):
        raise ValueError("output omits dominant-primary witnesses")
    offset = 0
    for source_index, path in enumerate(inputs):
        with uproot.open(path) as source:
            tree = source[PREFIX + "RJPhotonCandidateV1"]
            missing = set(fields) - set(tree.keys())
            if require_complete and missing:
                raise ValueError("dominant-primary source witness unavailable")
            requested = [*identity, *[key for key in fields if key in tree]]
            for rows in tree.iterate(requested, step_size=VALIDATION_STEP_SIZE, library="np"):
                if retained_photons is not None:
                    mask = np.asarray([_key(source_index, hi, lo) in retained_photons
                        for hi, lo in zip(rows[identity[0]], rows[identity[1]])], dtype=bool)
                    rows = {key: values[mask] for key, values in rows.items()}
                size = len(rows[identity[0]])
                actual = output.arrays(["source_file_index", *identity, *fields],
                    entry_start=offset, entry_stop=offset + size, library="np")
                if len(actual["source_file_index"]) != size or np.any(actual["source_file_index"] != source_index):
                    raise ValueError("dominant-primary source/cardinality mismatch")
                for key in identity:
                    if not np.array_equal(actual[key], rows[key]):
                        raise ValueError("dominant-primary identity/order mismatch")
                for key, default in defaults.items():
                    expected = rows.get(key, np.full(size, default))
                    if not np.array_equal(actual[key], expected, equal_nan=True):
                        raise ValueError(f"dominant-primary source readback mismatch: {key}")
                state = actual["dominant_truth_state"]
                if not np.all(np.isin(state, [-2, -1, 0, 1, 2, 3])):
                    raise ValueError("unknown dominant-primary state")
                if require_complete and np.any(np.isin(state, [-2, 0, 3])):
                    raise ValueError("incomplete dominant-primary evidence; not a background label")
                valid = state == 2
                if np.any(valid & ((actual["dominant_truth_track_id"] <= 0) |
                                  (actual["dominant_truth_pid"] == 0) |
                                  ~np.isin(actual["dominant_truth_evaluator_mode"], [1, 2]) |
                                  ~np.isfinite(actual["dominant_truth_energy_contribution"]) |
                                  (actual["dominant_truth_energy_contribution"] < 0))):
                    raise ValueError("invalid declared dominant-primary association")
                offset += size
    if offset != output.num_entries:
        raise ValueError("unexpected dominant-primary rows")
    return "dominant_primary_source_readback=PASS"


def validate_jet_view_response(inputs, root, *, compact: bool) -> list[str]:
    """Prove exact view transport and the independent final jet partition."""
    import build_photonjet_collaboration_tree as builder

    output_jets = root["jets"].arrays(
        ["source_file_index", "jet_id_hi", "jet_id_lo", "jet_input_identity",
         "jet_subtraction_identity", "jet_view_provenance_state"], library="np"
    )
    output_response = root["jetResponseLinks"].arrays(library="np")
    output_availability = root["jetViewAvailability"].arrays(library="np")
    seen_response_ids: set[tuple[int, int, int]] = set()
    supported_sources = 0
    for source_index, path in enumerate(inputs):
        with uproot.open(path) as source:
            events = builder._arrays(source[PREFIX + "RJEventV1"], builder.EVENT_INPUT_BRANCHES)
            jets = builder._arrays(source[PREFIX + "RJJetV1"], builder.JET_INPUT_BRANCHES)
            truth_jets = builder._arrays(
                source[PREFIX + "RJTruthJetV1"], builder.TRUTH_JET_INPUT_BRANCHES
            )
            links = builder._arrays(
                source[PREFIX + "RJRecoTruthLinkV1"], builder.LINK_INPUT_BRANCHES
            )
        available_views, availability_rows = builder._reco_jet_view_availability(events)
        availability_mask = output_availability["source_file_index"] == source_index
        if int(np.count_nonzero(availability_mask)) != len(availability_rows):
            raise ValueError("jet view availability source/cardinality mismatch")
        for name, dtype in builder.JET_VIEW_AVAILABILITY_SCHEMA.items():
            if name == "source_file_index":
                expected = np.full(len(availability_rows), source_index, dtype="int32")
            elif name == "source_entry":
                expected = np.asarray([row["event_index"] for row in availability_rows], dtype="int64")
            elif name in ("event_id_hi", "event_id_lo"):
                expected = np.asarray([events[name][row["event_index"]] for row in availability_rows], dtype=dtype)
            else:
                expected = np.asarray([row[name] for row in availability_rows],
                                      dtype=object if dtype == "string" else dtype)
            if not np.array_equal(output_availability[name][availability_mask], expected):
                raise ValueError(f"jet view availability readback mismatch: {source_index}/{name}")
        expected_input, expected_subtraction, expected_state = builder._jet_view_columns(jets)
        jet_mask = output_jets["source_file_index"] == source_index
        if int(np.count_nonzero(jet_mask)) != len(jets["jet_id_hi"]):
            raise ValueError("jet view source/cardinality mismatch")
        for name, expected in (
            ("jet_id_hi", jets["jet_id_hi"]), ("jet_id_lo", jets["jet_id_lo"]),
            ("jet_input_identity", expected_input),
            ("jet_subtraction_identity", expected_subtraction),
            ("jet_view_provenance_state", expected_state),
        ):
            if not np.array_equal(output_jets[name][jet_mask], expected):
                raise ValueError(f"jet view source readback mismatch: {source_index}/{name}")
        reco_lookup = builder._local_identity_lookup(
            jets, "jet_id_hi", "jet_id_lo", "reco jet identity"
        )
        truth_lookup = builder._local_identity_lookup(
            truth_jets, "truth_jet_id_hi", "truth_jet_id_lo", "truth jet identity"
        )
        expected_rows = builder._jet_response_rows(
            source_index, jets, truth_jets, links, reco_lookup, truth_lookup, available_views
        )
        response_mask = output_response["source_file_index"] == source_index
        actual_count = int(np.count_nonzero(response_mask))
        if actual_count != len(expected_rows):
            raise ValueError("per-view jet response cardinality mismatch")
        for name, dtype in builder.JET_RESPONSE_SCHEMA.items():
            expected = np.asarray(
                [row[name] for row in expected_rows], dtype=object if dtype == "string" else dtype
            )
            actual = output_response[name][response_mask]
            equal = (np.array_equal(actual, expected) if dtype == "string" else
                     np.array_equal(actual, expected, equal_nan=True))
            if not equal:
                raise ValueError(f"per-view jet response mismatch: {source_index}/{name}")
        states = output_response["jet_response_state"][response_mask]
        if available_views is not None:
            if np.any(~np.isin(states, [0, 1, 2])):
                raise ValueError("jet view response contains an unknown authority state")
            if np.any(states == 0):
                unsupported = response_mask & (output_response["jet_response_state"] == 0)
                if (np.any(output_response["jet_input_identity"][unsupported] != "") or
                        np.any(output_response["jet_subtraction_identity"][unsupported] != "")):
                    raise ValueError("unsupported truth-only response inferred a reconstruction view")
            else:
                supported_sources += 1
        elif actual_count and np.any(states != 0):
            raise ValueError("legacy source jet response was silently promoted")
        for high, low in zip(output_response["link_id_hi"][response_mask],
                             output_response["link_id_lo"][response_mask]):
            identity = _key(source_index, high, low)
            if identity in seen_response_ids:
                raise ValueError(f"duplicate jetResponseLinks identity: {identity}")
            seen_response_ids.add(identity)

    if not compact:
        jets = root["jets"].arrays(
            ["source_file_index", "jet_id_hi", "jet_id_lo", "jet_input_identity",
             "jet_subtraction_identity", "jet_view_provenance_state"], library="np"
        )
        lookup = {_key(s, hi, lo): (str(inp), str(sub), int(state))
                  for s, hi, lo, inp, sub, state in zip(
                      jets["source_file_index"], jets["jet_id_hi"], jets["jet_id_lo"],
                      jets["jet_input_identity"], jets["jet_subtraction_identity"],
                      jets["jet_view_provenance_state"])}
        for pairs in root["photonJets"].iterate(
            ["source_file_index", "jet_id_hi", "jet_id_lo", "jet_input_identity",
             "jet_subtraction_identity", "jet_view_provenance_state"],
            step_size=VALIDATION_STEP_SIZE, library="np"
        ):
            for i, source_index in enumerate(pairs["source_file_index"]):
                expected = lookup[_key(source_index, pairs["jet_id_hi"][i], pairs["jet_id_lo"][i])]
                actual = (str(pairs["jet_input_identity"][i]),
                          str(pairs["jet_subtraction_identity"][i]),
                          int(pairs["jet_view_provenance_state"][i]))
                if actual != expected:
                    raise ValueError("expanded pair jet-view join mismatch")

    _validate_truth_link_targets(root, "jetResponseLinks")
    status = ("PASS" if supported_sources == len(inputs) else
              "UNSUPPORTED_MISSING_VIEW_AVAILABILITY_OR_LEGACY_PROVENANCE")
    lines = ["jet_view_source_and_pair_transport=PASS",
             f"per_view_jet_response={status}"]
    if status == "PASS":
        lines.append("jet_candidate_edges_preserved_and_final_partition_rebuilt=PASS")
    else:
        lines.extend([
            "jet_candidate_edges_preserved=PASS",
            "jet_final_partition=UNSUPPORTED_WITHOUT_EXACT_EVENT_VIEW_AVAILABILITY",
        ])
    return lines


def validate_truth_witnesses(inputs, root):
    """Prove new truth witnesses survive transport; absent legacy fields stay unknown."""
    integer_fields = ("truth_isolation_valid", "g4_photon_valid", "hepmc_association_valid",
                      "analysis_signal_r03", "native_track_id", "native_vertex_id", "embedding_id",
                      "sample_source_role", "generator_occurrence_embedding_id")
    fields = ("truth_isolation_r03", "truth_isolation_r04", *integer_fields)
    output = root["truthPhotons"]
    if not output.num_entries:
        return "truth_dual_cone_witnesses=EMPTY"
    if not set(fields).issubset(output.keys()):
        raise ValueError("output omits truth isolation/association witnesses")
    offset = 0
    for source_index, path in enumerate(inputs):
        with uproot.open(path) as source:
            tree = source[PREFIX + "RJTruthPhotonV1"]
            requested = ["truth_photon_id_hi", "truth_photon_id_lo", *[k for k in fields if k in tree]]
            for rows in tree.iterate(requested, step_size=VALIDATION_STEP_SIZE, library="np"):
                size = len(rows["truth_photon_id_hi"])
                actual = output.arrays(["source_file_index", *requested[:2], *fields],
                    entry_start=offset, entry_stop=offset + size, library="np")
                if len(actual["source_file_index"]) != size or np.any(actual["source_file_index"] != source_index):
                    raise ValueError("truth witness source/cardinality mismatch")
                for key in requested[:2]:
                    if not np.array_equal(actual[key], rows[key]):
                        raise ValueError("truth witness identity/order mismatch")
                for key in fields:
                    expected = rows.get(key, np.full(size, -1 if key in integer_fields else np.nan))
                    if not np.array_equal(actual[key], expected, equal_nan=True):
                        raise ValueError(f"truth witness source readback mismatch: {key}")
                offset += size
    if offset != output.num_entries:
        raise ValueError("unexpected truth witness rows")
    return "truth_dual_cone_source_readback=PASS"


def validate_normalized_layout(root, inputs, system, evaluator, input_count):
    """Validate compact transport without materializing an expanded ROOT copy.

    Tables are shard-sized. Only identity/kinematic columns are read for pair
    joins; no photon/event payload is replicated for each pair on disk.
    """
    import build_photonjet_collaboration_tree as builder
    schemas = builder.output_schemas("normalized_v1")
    if "eventTree" in root:
        raise ValueError("normalized layout contains a redundant expanded eventTree")
    for name, schema in schemas.items():
        if not set(schema).issubset(root[name].keys()):
            raise ValueError(f"normalized layout missing columns in {name}")
    for name in ("photons", "jets", "photonJets"):
        forbidden = set(builder.EVENT_SCHEMA) - set(builder.EVENT_JOIN_SCHEMA)
        if forbidden.intersection(root[name].keys()):
            raise ValueError(f"normalized layout repeats event payload in {name}")
    if (set(builder.PHOTON_FIELDS) | set(builder.JET_FIELDS)).difference(
            builder.PAIR_OBJECT_IDS).intersection(root["photonJets"].keys()):
        raise ValueError("normalized layout repeats object payload in photonJets")

    events = root["events"].arrays([*builder.EVENT_JOIN_SCHEMA, "centrality"], library="np")
    event_lookup = {_key(s, hi, lo): i for i, (s, hi, lo) in enumerate(zip(
        events["source_file_index"], events["event_id_hi"], events["event_id_lo"]))}
    photons = root["photons"].arrays(library="np")
    event_rows = np.asarray([event_lookup[_key(s, hi, lo)] for s, hi, lo in zip(
        photons["source_file_index"], photons["event_id_hi"], photons["event_id_lo"])], dtype="int64")
    c, et = events["centrality"][event_rows], photons["photon_et"]
    domain = (et >= MODEL_ET_MIN_GEV) & (et < MODEL_ET_MAX_GEV)
    if system == "auau":
        domain &= np.isfinite(c) & (c >= 0) & (c < 80)
    state, score = photons["bdt_evaluation_state"], photons["bdt_score"]
    scored = state != 0
    if (not np.all(np.isin(state, [0, 1, 2])) or
            np.any(scored != np.isfinite(score)) or np.any(scored & ~domain)):
        raise ValueError("normalized scored/unscored or model-domain state is inconsistent")
    expected_thresholds = _expected_thresholds(system, et, c)
    for name, expected in zip(("bdt_tight_threshold", "bdt_nontight_low_threshold",
                               "bdt_nontight_high_threshold"), expected_thresholds):
        expected[~domain] = np.nan
        if not np.allclose(photons[name], expected, rtol=0, atol=1e-12, equal_nan=True):
            raise ValueError(f"normalized threshold mismatch: {name}")
    tight, low, high = expected_thresholds
    for name, selected in (
        ("bdt_is_tight", score > tight),
        ("bdt_is_nontight", (score > low) & (score < high)),
        ("bdt_is_not_tight", ~(score > tight)),
    ):
        expected = np.where(scored, selected.astype("int32"), -1)
        if not np.array_equal(photons[name], expected):
            raise ValueError(f"unscored capture became a physics selection: {name}")
    if evaluator is not None:
        if np.any(scored & (state != 1)):
            raise ValueError("supplied model was not used for every scored row")
        for i in np.flatnonzero(scored):
            if photons["bdt_input_count"][i] != input_count:
                raise ValueError("scored input dimension mismatch")
            values = np.asarray([photons[f"bdt_input_{j:02d}"][i] for j in range(input_count)])
            if not np.all(np.isfinite(values)) or not math.isclose(
                    float(score[i]), evaluator(values), rel_tol=1e-6, abs_tol=1e-7):
                raise ValueError(f"normalized BDT score/input mismatch at row {i}")

    # Preserve the raw H70 features even for candidates that cannot be scored.
    for source_index, path in enumerate(inputs):
        with uproot.open(path) as source:
            views = source[PREFIX + "RJShowerFeatureViewV1"].arrays(
                ["candidate_id_hi", "candidate_id_lo", "definition_name", "ordered_features"],
                library="ak", how=dict)
            lookup = {(_key(source_index, hi, lo)): i for i, (hi, lo, s) in enumerate(zip(
                photons["candidate_id_hi"], photons["candidate_id_lo"], photons["source_file_index"]))
                if s == source_index}
            for i, definition in enumerate(views["definition_name"]):
                if str(definition).upper() != "H70":
                    continue
                key = _key(source_index, views["candidate_id_hi"][i], views["candidate_id_lo"][i])
                if key not in lookup:
                    raise ValueError("raw H70 candidate disappeared from normalized export")
                row = lookup[key]
                values = np.asarray(views["ordered_features"][i], dtype="float64")[:14]
                actual = np.asarray([photons[f"bdt_input_{j:02d}"][row] for j in range(len(values))])
                if photons["bdt_input_count"][row] != len(values) or not np.array_equal(
                        actual, values, equal_nan=True):
                    raise ValueError("normalized raw H70 input readback mismatch")

    parents = {}
    for name, prefix, fields in (("photons", "candidate", ("photon_et", "photon_phi")),
                                  ("jets", "jet", ("jet_pt", "jet_phi"))):
        rows = root[name].arrays([*builder.EVENT_JOIN_SCHEMA, f"{prefix}_id_hi",
                                  f"{prefix}_id_lo", *fields], library="np")
        per_event, lookup = {}, {}
        for i, (s, hi, lo) in enumerate(zip(rows["source_file_index"],
                rows[f"{prefix}_id_hi"], rows[f"{prefix}_id_lo"])):
            event = _key(s, rows["event_id_hi"][i], rows["event_id_lo"][i])
            key = _key(s, hi, lo)
            if key in lookup:
                raise ValueError("duplicate normalized parent identity")
            if rows["source_entry"][i] != events["source_entry"][event_lookup[event]]:
                raise ValueError("normalized source_entry/event join mismatch")
            lookup[key] = (event, per_event.get(event, 0), *(float(rows[k][i]) for k in fields))
            per_event[event] = per_event.get(event, 0) + 1
        parents[name] = lookup
    for pairs in root["photonJets"].iterate(library="np", step_size=VALIDATION_STEP_SIZE):
        for i, s in enumerate(pairs["source_file_index"]):
            event = _key(s, pairs["event_id_hi"][i], pairs["event_id_lo"][i])
            p = parents["photons"].get(_key(s, pairs["candidate_id_hi"][i], pairs["candidate_id_lo"][i]))
            j = parents["jets"].get(_key(s, pairs["jet_id_hi"][i], pairs["jet_id_lo"][i]))
            if (p is None or j is None or p[:2] != (event, pairs["photon_index"][i]) or
                    j[:2] != (event, pairs["jet_index"][i]) or
                    pairs["source_entry"][i] != events["source_entry"][event_lookup[event]]):
                raise ValueError("normalized pair parent identity/event/local index mismatch")
            if not all(math.isfinite(v) for v in (*p[2:], *j[2:])) or p[2] <= 0:
                raise ValueError("normalized pair has invalid kinematic inputs")
            dphi = abs(math.atan2(math.sin(p[3] - j[3]), math.cos(p[3] - j[3])))
            if (not math.isclose(pairs["xjgamma"][i], j[2] / p[2], rel_tol=1e-12, abs_tol=1e-12) or
                    not math.isclose(pairs["delta_phi"][i], dphi, rel_tol=1e-12, abs_tol=1e-12)):
                raise ValueError("normalized pair xjgamma/delta_phi mismatch")
    for observable in ("xjgamma", "delta_phi"):
        count = _compare_float_streams(_iter_source_pair_observable(inputs, None, observable),
            (np.asarray(chunk[observable], dtype=np.float64) for chunk in
             root["photonJets"].iterate([observable], library="np", step_size=VALIDATION_STEP_SIZE)),
            f"normalized source pair {observable}")
        if count != root["photonJets"].num_entries:
            raise ValueError("normalized source pair count mismatch")
    return ["normalized_single_payload_layout=PASS", "normalized_pair_joins_and_kinematics=PASS",
            "raw_H70_capture_readback=PASS",
            f"scored_candidates={int(np.sum(scored))}; unscored_candidates={int(np.sum(~scored))}",
            "raw_transport_only_not_model_or_physics_acceptance"]


def validate(
    inputs: list[Path],
    output: Path,
    system: str,
    model_path: Path | None = None,
    model_input_count: int | None = None,
    require_scaled_bit30: bool = False,
    require_complete_interface: bool = False,
    require_mbd_pmt: bool = False,
    require_centrality_replay: bool = False,
    require_trigger_scalers: bool = False,
) -> list[str]:
    import build_photonjet_collaboration_tree as builder
    inputs, output = builder.validate_io_paths(inputs, output)
    evaluator = _make_evaluator(model_path, model_input_count)
    if system not in {"pp", "auau"}:
        raise ValueError("system must be pp or auau")
    expected_model_input_count = 11 if system == "pp" else 14
    if (
        model_path is not None
        and model_input_count is not None
        and model_input_count != expected_model_input_count
    ):
        raise ValueError(
            f"BDT input dimension mismatch for {system}: "
            f"{model_input_count} != {expected_model_input_count}"
        )
    lines: list[str] = []
    with uproot.open(output) as root:
        layout = str(root["storage_layout"]) if "storage_layout" in root else "expanded"
        if layout not in builder.LAYOUTS:
            raise ValueError(f"unknown storage layout: {layout}")
        compact = layout == "normalized_v1"
        if compact and "source_files" not in root:
            raise ValueError("normalized package has no source-file locator/provenance map")
        if "source_files" in root:
            actual_sources = json.loads(str(root["source_files"]))
            expected_sources = []
            for source_index, path in enumerate(inputs):
                with uproot.open(path) as source:
                    expected_sources.append(builder.source_file_record(source, path, source_index))
            # The original path is a locator, not an identity. Exact copies may
            # move to an offline drive; UUID and stored upstream provenance stay bound.
            def identities(records):
                return [{key: value for key, value in record.items() if key != "input_path"}
                        for record in records]
            if identities(actual_sources) != identities(expected_sources):
                raise ValueError("source-file locator/provenance readback mismatch")
            lines.append("source_file_locator_and_upstream_provenance=PASS")
        retained_photons = (None if compact else _retained_photon_identities(
            inputs, model_input_count if model_path is not None else None))
        _report_memory("retained_photon_inventory")
        required = set(SOURCE_TREES) | {"jetResponseLinks", "jetViewAvailability"} | (set() if compact else {"eventTree"})
        missing = sorted(required.difference(root.keys(cycle=False)))
        if missing:
            raise ValueError("missing output trees: " + ", ".join(missing))
        trigger_scaler_interface.validate_copy(inputs, root, require=require_trigger_scalers)
        gl1_trigger_interface.validate_copy(inputs, root)
        lines.append("direct_gl1_decisions_and_run_configuration=PASS")
        lines.append("trigger_scaler_source_copy_and_event_joins=PASS")
        lines.append(validate_truth_witnesses(inputs, root))
        lines.append(validate_dominant_truth_witnesses(inputs, root,
                                                       retained_photons=retained_photons))
        lines.append(validate_interface(inputs, root, system,
                                        require_complete=require_complete_interface))
        if require_centrality_replay and ("centrality_replay_version" not in root["events"].keys() or
                                          np.any(root["events"]["centrality_replay_version"].array(library="np") != 1)):
            raise ValueError("complete centrality replay capture is required")
        if require_mbd_pmt and ("mbd_pmt_available" not in root["events"].keys() or
                               np.any(root["events"]["mbd_pmt_available"].array(library="np") != 1)):
            raise ValueError("per-PMT MBD unavailable")

        expected_counts: dict[str, int] = {}
        with tempfile.TemporaryDirectory(prefix=".photonjet_part_identity_proof_") as proof_dir:
            proof_root = Path(proof_dir)
            for name in SOURCE_TREES:
                expected_count = _compare_identity_streams(
                    _iter_source_identities(inputs, name, retained_photons),
                    _iter_output_identities(root, name),
                    proof_root / f"{name}.bin",
                    name,
                )
                observed = int(root[name].num_entries)
                if observed != expected_count:
                    raise ValueError(f"{name} count mismatch: {observed} != {expected_count}")
                expected_counts[name] = expected_count
                lines.append(f"{name}_identity_and_count=PASS ({observed})")
                gc.collect()
                _report_memory(f"identity_{name}")

        if not compact and root["eventTree"].num_entries != expected_counts["events"]:
            raise ValueError("eventTree count does not match events")

        events = root["events"].arrays(
            ["source_file_index", "event_id_hi", "event_id_lo", "scaled_trigger_bits", "scaled_bit30"],
            library="np",
        )
        event_ids = _ids(events["source_file_index"], events["event_id_hi"], events["event_id_lo"])
        if require_scaled_bit30:
            reconstructed = ((events["scaled_trigger_bits"] >> np.uint64(30)) & np.uint64(1)).astype("int32")
            if not np.all(reconstructed == 1) or not np.all(events["scaled_bit30"] == 1):
                raise ValueError("pp data scaled-bit-30 membership failed")
            lines.append("scaled_bit30_membership=PASS")

        for child_name in ("photons", "jets", "photonJets", "truthPhotons", "truthJets",
                           "recoTruthLinks", "jetResponseLinks", "jetViewAvailability"):
            for child in root[child_name].iterate(
                ["source_file_index", "event_id_hi", "event_id_lo"],
                step_size=VALIDATION_STEP_SIZE,
                library="np",
            ):
                for source, high, low in zip(
                    child["source_file_index"], child["event_id_hi"], child["event_id_lo"]
                ):
                    if _key(source, high, low) not in event_ids:
                        raise ValueError(f"{child_name} contains an event without a parent")
            gc.collect()
        lines.append("event_joins=PASS")
        _report_memory("event_joins")
        lines.extend(validate_jet_view_response(inputs, root, compact=compact))

        if compact:
            lines.extend(validate_normalized_layout(root, inputs, system, evaluator,
                                                     model_input_count))
            _validate_truth_link_targets(root)
            lines.append("truth_link_targets=PASS")
            return lines

        photons = root["photons"].arrays(library="np")
        if len(photons["bdt_score"]):
            if not np.all(np.isfinite(photons["bdt_score"])):
                raise ValueError("collaborator bdt_score contains a non-finite value")
            if model_input_count is not None and not np.all(photons["bdt_input_count"] == model_input_count):
                raise ValueError("stored BDT input dimension does not match the included model")
            if not np.all(np.isfinite(photons["iso_r03"])) or not np.all(np.isfinite(photons["iso_r04"])):
                raise ValueError("a photon lacks an R=0.3 or R=0.4 isolation witness")
            expected_tight, expected_low, expected_high = _expected_thresholds(
                system,
                np.asarray(photons["photon_et"], dtype="float64"),
                np.asarray(photons["centrality"], dtype="float64"),
            )
            for branch, expected in (
                ("bdt_tight_threshold", expected_tight),
                ("bdt_nontight_low_threshold", expected_low),
                ("bdt_nontight_high_threshold", expected_high),
            ):
                if not np.all(np.isfinite(photons[branch])) or not np.allclose(
                    photons[branch], expected, rtol=0.0, atol=1e-12
                ):
                    raise ValueError(f"{branch} does not match the established {system} selection surface")
            if not np.all(photons["bdt_nontight_low_threshold"] < photons["bdt_nontight_high_threshold"]):
                raise ValueError("bounded non-tight thresholds are not ordered")
            expected_tight_state = (photons["bdt_score"] > expected_tight).astype("int32")
            expected_nontight_state = (
                (photons["bdt_score"] > expected_low) & (photons["bdt_score"] < expected_high)
            ).astype("int32")
            if not np.array_equal(photons["bdt_is_tight"], expected_tight_state):
                raise ValueError("bdt_is_tight does not match score and threshold")
            if not np.array_equal(photons["bdt_is_nontight"], expected_nontight_state):
                raise ValueError("bdt_is_nontight does not match score and bounded thresholds")
            if not np.array_equal(photons["bdt_is_not_tight"], 1 - expected_tight_state):
                raise ValueError("bdt_is_not_tight does not match the tight-score complement")
            if evaluator is not None:
                for row_index in range(len(photons["bdt_score"])):
                    count = int(photons["bdt_input_count"][row_index])
                    values = np.asarray([photons[f"bdt_input_{index:02d}"][row_index] for index in range(count)], dtype="float64")
                    if not np.all(np.isfinite(values)):
                        raise ValueError(f"non-finite BDT input at photon row {row_index}")
                    expected_score = evaluator(values)
                    if not math.isclose(float(photons["bdt_score"][row_index]), expected_score, rel_tol=1e-6, abs_tol=1e-7):
                        raise ValueError(f"BDT score mismatch at photon row {row_index}")
            lines.append(f"uniform_bdt_score=PASS ({len(photons['bdt_score'])})")
            lines.append(f"{system}_selection_thresholds_and_states=PASS")
            lines.append("isolation_r03_r04=PASS")

        for pair_start, pairs in enumerate(
            root["photonJets"].iterate(
                [
                    "photon_et", "photon_phi", "jet_pt", "jet_phi",
                    "delta_phi", "xjgamma",
                ],
                step_size=VALIDATION_STEP_SIZE,
                library="np",
            )
        ):
            finite_xj_inputs = (
                np.isfinite(pairs["photon_et"])
                & np.isfinite(pairs["jet_pt"])
                & (pairs["photon_et"] != 0)
            )
            if not np.all(finite_xj_inputs) or not np.all(np.isfinite(pairs["xjgamma"])):
                raise ValueError(f"non-finite xjgamma input or value in chunk {pair_start}")
            expected_xj = pairs["jet_pt"] / pairs["photon_et"]
            if not np.allclose(
                pairs["xjgamma"], expected_xj, rtol=1e-12, atol=1e-12
            ):
                raise ValueError(f"xjgamma is not jet_pt/photon_et in chunk {pair_start}")
            if (
                not np.all(np.isfinite(pairs["photon_phi"]))
                or not np.all(np.isfinite(pairs["jet_phi"]))
                or not np.all(np.isfinite(pairs["delta_phi"]))
            ):
                raise ValueError(f"non-finite delta_phi input or value in chunk {pair_start}")
            expected_delta_phi = _wrapped_abs_delta_phi(
                pairs["photon_phi"], pairs["jet_phi"]
            )
            if not np.allclose(
                pairs["delta_phi"], expected_delta_phi, rtol=1e-12, atol=1e-12
            ):
                raise ValueError(
                    f"delta_phi does not match wrapped photon/jet phi in chunk {pair_start}"
                )
        lines.append("pair_xjgamma_and_delta_phi_numeric_parity=PASS")

        for observable in ("xjgamma", "delta_phi"):
            compared_pairs = _compare_float_streams(
                _iter_source_pair_observable(inputs, retained_photons, observable),
                (
                    np.asarray(chunk[observable], dtype=np.float64)
                    for chunk in root["photonJets"].iterate(
                        [observable], step_size=VALIDATION_STEP_SIZE, library="np"
                    )
                ),
                f"source pair {observable}",
            )
            if compared_pairs != expected_counts["photonJets"]:
                raise ValueError(
                    f"source pair {observable} count mismatch: {compared_pairs} != "
                    f"{expected_counts['photonJets']}"
                )
            lines.append(f"source_pair_{observable}=PASS")
        gc.collect()
        _report_memory("pair_numeric_parity")

        event_branches = [
            "nphotons", "njets", "npairs", "pair_photon_index", "pair_jet_index",
            "pair_delta_phi", "pair_xjgamma", "photon_et", "photon_phi",
            "jet_pt", "jet_phi", "photon_bdt_score",
            "photon_iso_r04", "photon_iso_r04_threshold",
            "photon_iso_r04_nonisolated_threshold", "photon_encounter_ordinal",
            "photon_bdt_tight_threshold", "photon_bdt_nontight_low_threshold",
            "photon_bdt_nontight_high_threshold",
        ] + [f"leader_{region}_r04_index" for region in "ABCD"] + [
            f"leader_{region}_r04_complement_index" for region in "CD"
        ]
        photon_total = jet_total = pair_total = 0
        event_offset = 0
        for event_tree in root["eventTree"].iterate(
            event_branches,
            step_size=EVENT_VALIDATION_STEP_SIZE,
            library="ak",
            how=dict,
        ):
            photon_total += int(ak.sum(event_tree["nphotons"]))
            jet_total += int(ak.sum(event_tree["njets"]))
            pair_total += int(ak.sum(event_tree["npairs"]))
            chunk_length = len(event_tree["nphotons"])
            for local_event_index in range(chunk_length):
                event_index = event_offset + local_event_index
                nphotons = int(event_tree["nphotons"][local_event_index])
                njets = int(event_tree["njets"][local_event_index])
                photon_indices = np.asarray(
                    event_tree["pair_photon_index"][local_event_index], dtype="int32"
                )
                jet_indices = np.asarray(
                    event_tree["pair_jet_index"][local_event_index], dtype="int32"
                )
                npairs = int(event_tree["npairs"][local_event_index])
                pair_delta_phi = np.asarray(
                    event_tree["pair_delta_phi"][local_event_index], dtype="float64"
                )
                pair_xjgamma = np.asarray(
                    event_tree["pair_xjgamma"][local_event_index], dtype="float64"
                )
                if not (
                    len(photon_indices)
                    == len(jet_indices)
                    == len(pair_delta_phi)
                    == len(pair_xjgamma)
                    == npairs
                ):
                    raise ValueError(
                        f"eventTree pair-array lengths differ at event {event_index}"
                    )
                if np.any(photon_indices < 0) or np.any(photon_indices >= nphotons):
                    raise ValueError(f"eventTree photon local index out of range at event {event_index}")
                if np.any(jet_indices < 0) or np.any(jet_indices >= njets):
                    raise ValueError(f"eventTree jet local index out of range at event {event_index}")
                if len(photon_indices):
                    photon_et = np.asarray(
                        event_tree["photon_et"][local_event_index], dtype="float64"
                    )[photon_indices]
                    jet_pt = np.asarray(
                        event_tree["jet_pt"][local_event_index], dtype="float64"
                    )[jet_indices]
                    photon_phi = np.asarray(
                        event_tree["photon_phi"][local_event_index], dtype="float64"
                    )[photon_indices]
                    jet_phi = np.asarray(
                        event_tree["jet_phi"][local_event_index], dtype="float64"
                    )[jet_indices]
                    if (
                        not np.all(np.isfinite(photon_et))
                        or not np.all(np.isfinite(jet_pt))
                        or np.any(photon_et == 0)
                        or not np.all(np.isfinite(pair_xjgamma))
                    ):
                        raise ValueError(
                            f"eventTree non-finite xjgamma input or value at event {event_index}"
                        )
                    if not np.allclose(
                        pair_xjgamma, jet_pt / photon_et, rtol=1e-12, atol=1e-12
                    ):
                        raise ValueError(
                            f"eventTree xjgamma local pair relation failed at event {event_index}"
                        )
                    if (
                        not np.all(np.isfinite(photon_phi))
                        or not np.all(np.isfinite(jet_phi))
                        or not np.all(np.isfinite(pair_delta_phi))
                    ):
                        raise ValueError(
                            f"eventTree non-finite delta_phi input or value at event {event_index}"
                        )
                    expected_delta_phi = _wrapped_abs_delta_phi(photon_phi, jet_phi)
                    if not np.allclose(
                        pair_delta_phi, expected_delta_phi, rtol=1e-12, atol=1e-12
                    ):
                        raise ValueError(
                            f"eventTree delta_phi local pair relation failed at event {event_index}"
                        )
                scores = np.asarray(
                    event_tree["photon_bdt_score"][local_event_index], dtype="float64"
                )
                isolation = np.asarray(
                    event_tree["photon_iso_r04"][local_event_index], dtype="float64"
                )
                isolated_thresholds = np.asarray(
                    event_tree["photon_iso_r04_threshold"][local_event_index], dtype="float64"
                )
                nonisolated_thresholds = np.asarray(
                    event_tree["photon_iso_r04_nonisolated_threshold"][local_event_index],
                    dtype="float64",
                )
                photon_et = np.asarray(
                    event_tree["photon_et"][local_event_index], dtype="float64"
                )
                encounter_ordinal = np.asarray(
                    event_tree["photon_encounter_ordinal"][local_event_index], dtype="int32"
                )
                tight_thresholds = np.asarray(
                    event_tree["photon_bdt_tight_threshold"][local_event_index],
                    dtype="float64",
                )
                nontight_low = np.asarray(
                    event_tree["photon_bdt_nontight_low_threshold"][local_event_index],
                    dtype="float64",
                )
                nontight_high = np.asarray(
                    event_tree["photon_bdt_nontight_high_threshold"][local_event_index],
                    dtype="float64",
                )
                for region in "ABCD":
                    branch = f"leader_{region}_r04_index"
                    expected_leader = _leader_index(
                        scores, tight_thresholds, nontight_low, nontight_high,
                        isolation, isolated_thresholds, nonisolated_thresholds,
                        photon_et, encounter_ordinal, region,
                    )
                    if int(event_tree[branch][local_event_index]) != expected_leader:
                        raise ValueError(f"leader mismatch: {branch} event {event_index}")
                for region in "CD":
                    branch = f"leader_{region}_r04_complement_index"
                    expected_leader = _leader_index(
                        scores, tight_thresholds, nontight_low, nontight_high,
                        isolation, isolated_thresholds, nonisolated_thresholds,
                        photon_et, encounter_ordinal, region, complement=True,
                    )
                    if int(event_tree[branch][local_event_index]) != expected_leader:
                        raise ValueError(f"leader mismatch: {branch} event {event_index}")
            event_offset += chunk_length
        if photon_total != expected_counts["photons"]:
            raise ValueError("eventTree photon-array total mismatch")
        if jet_total != expected_counts["jets"]:
            raise ValueError("eventTree jet-array total mismatch")
        if pair_total != expected_counts["photonJets"]:
            raise ValueError("eventTree pair-array total mismatch")
        lines.append("eventTree_arrays_local_indices_and_leaders=PASS")
        lines.append("eventTree_pair_xjgamma_and_delta_phi=PASS")
        gc.collect()
        _report_memory("event_tree")

        _validate_truth_link_targets(root)
        lines.append("truth_link_targets=PASS")

        radii = np.unique(root["jets"]["jet_radius"].array(library="np"))
        lines.append("jet_radius_inventory=" + (",".join(f"{value:.3g}" for value in radii) if len(radii) else "EMPTY"))
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--system", required=True, choices=("pp", "auau"))
    parser.add_argument("--model", type=Path)
    parser.add_argument("--model-input-count", type=int, choices=(11, 14))
    parser.add_argument("--require-scaled-bit30", action="store_true")
    parser.add_argument("--require-complete-interface", action="store_true")
    parser.add_argument("--require-mbd-pmt", action="store_true")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--require-centrality-replay", action="store_true")
    parser.add_argument("--require-trigger-scalers", action="store_true")
    args = parser.parse_args()
    lines = validate(
        args.input, args.output, args.system, model_path=args.model,
        model_input_count=args.model_input_count,
        require_scaled_bit30=args.require_scaled_bit30,
        require_complete_interface=args.require_complete_interface,
        require_mbd_pmt=args.require_mbd_pmt,
        require_centrality_replay=args.require_centrality_replay,
        require_trigger_scalers=args.require_trigger_scalers,
    )
    text = "PARITY=PASS\n" + "\n".join(lines) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
