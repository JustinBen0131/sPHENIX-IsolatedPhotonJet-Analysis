"""Lossless, source-scoped transport of cumulative GL1 counters; no luminosity inference."""
from __future__ import annotations

import re
import numpy as np

TREE = "RJTriggerScalersV1"
OUTPUT_TREE = "triggerScalers"
STEP = 2048
SCHEMA = {
    "source_file_index": "int32", "snapshot_id": "uint64",
    "run": "int32", "segment": "int32",
    "input_uri_hash": "string", "source_manifest_sha256": "string",
    "first_event_sequence": "int64", "last_event_sequence": "int64",
    "first_physical_event": "int64", "last_physical_event": "int64",
    "first_bco": "uint64", "last_bco": "uint64", "observed_events": "uint64",
    "state": "int32", "packet_version": "int32", "packet_status": "uint64",
    "discontinuity_flags": "uint32",
    "raw": "64 * uint64", "live": "64 * uint64", "scaled": "64 * uint64",
}
BRANCHES = tuple(key for key in SCHEMA if key != "source_file_index")
REFERENCE = "trigger_scaler_snapshot_id"


def _source_tree(root):
    if TREE not in root:
        return None
    tree = root[TREE]
    missing = set(BRANCHES) - set(tree.keys())
    if missing:
        raise ValueError(f"trigger scaler branches missing: {sorted(missing)}")
    return tree


def copy_scalers(root, output, events, source_index, *, require=False):
    """Copy fixed-size counter arrays in bounded batches and prove retained-event joins."""
    n = len(events["event_id_hi"])
    refs = np.asarray(events.get(REFERENCE, np.zeros(n, dtype="uint64")))
    if refs.dtype.kind != "u" or refs.dtype.itemsize != 8:
        raise ValueError("trigger scaler reference must be uint64")
    tree = _source_tree(root)
    if tree is None or tree.num_entries == 0:
        if require or np.any(refs != 0):
            raise ValueError("source has no trigger scaler snapshots for its required references")
        return 0
    if REFERENCE not in events or np.any(refs == 0):
        raise ValueError("source with trigger scalers has missing event snapshot references")

    metadata = {key: [] for key in ("snapshot_id", "run", "first_event_sequence", "last_event_sequence")}
    count = 0
    source_identity = None
    for batch in tree.iterate(list(BRANCHES), step_size=STEP, library="np"):
        size = len(batch["snapshot_id"])
        ids = np.asarray(batch["snapshot_id"])
        if ids.dtype.kind != "u" or ids.dtype.itemsize != 8:
            raise ValueError("snapshot ids must be uint64")
        if not np.array_equal(ids, np.arange(count + 1, count + size + 1, dtype="uint64")):
            raise ValueError("trigger scaler ids are not consecutive and unique within source")
        for field in ("raw", "live", "scaled"):
            values = np.asarray(batch[field])
            if values.shape != (size, 64) or values.dtype.kind != "u" or values.dtype.itemsize != 8:
                raise ValueError(f"{field} must contain 64 lossless uint64 counters per snapshot")
        if np.any(batch["observed_events"] == 0):
            raise ValueError("empty trigger scaler observation interval")
        if np.any(batch["first_event_sequence"] > batch["last_event_sequence"]):
            raise ValueError("reversed trigger scaler event interval")
        states = np.asarray(batch["state"])
        if np.any((states < 0) | (states > 4)):
            raise ValueError("unknown trigger scaler validity state")
        if require and np.any(states != 0):
            raise ValueError("required trigger scalers include missing, unsupported or bad packets")
        for run, segment, input_hash, manifest_hash in zip(
                batch["run"], batch["segment"], batch["input_uri_hash"], batch["source_manifest_sha256"]):
            identity = (int(run), int(segment), str(input_hash), str(manifest_hash))
            if identity[0] <= 0 or identity[1] < 0 or any(
                    not re.fullmatch(r"[0-9a-f]{64}", value) for value in identity[2:]):
                raise ValueError("invalid trigger scaler physical-source identity")
            if source_identity is None:
                source_identity = identity
            if source_identity != identity:
                raise ValueError("trigger scaler physical-source identity changes within input")
        for key in metadata:
            metadata[key].append(np.asarray(batch[key]))
        output.extend({"source_file_index": np.full(size, source_index, dtype="int32"), **batch})
        count += size

    meta = {key: np.concatenate(values) for key, values in metadata.items()}
    positions = np.searchsorted(meta["snapshot_id"], refs)
    if np.any(positions >= count) or not np.array_equal(meta["snapshot_id"][positions], refs):
        raise ValueError("event references a missing source-local trigger snapshot")
    event_sequence = np.asarray(events["event_sequence"])
    if np.any(meta["run"][positions] != np.asarray(events["run"])) or np.any(
            (event_sequence < meta["first_event_sequence"][positions]) |
            (event_sequence > meta["last_event_sequence"][positions])):
        raise ValueError("event run/sequence is outside its trigger snapshot interval")
    return count


def validate_copy(inputs, output_root, *, require=False):
    """Compare every saved counter/provenance field directly to its input, in batches."""
    import uproot
    if OUTPUT_TREE not in output_root:
        if require or any(_input_has_scalers(path) for path in inputs):
            raise ValueError("collaborator output omitted trigger scalers")
        return
    out = output_root[OUTPUT_TREE]
    offset = 0
    for source_index, path in enumerate(inputs):
        with uproot.open(path) as source:
            tree = _source_tree(source)
            if tree is None:
                if require:
                    raise ValueError("required original scaler tree is absent")
                continue
            if require and tree.num_entries == 0:
                raise ValueError("required original scaler tree is empty")
            for batch in tree.iterate(list(BRANCHES), step_size=STEP, library="np"):
                if require and np.any(batch["state"] != 0):
                    raise ValueError("required trigger scaler readback includes invalid packets")
                size = len(batch["snapshot_id"])
                actual = out.arrays(list(SCHEMA), entry_start=offset, entry_stop=offset + size, library="np")
                if len(actual["snapshot_id"]) != size or np.any(actual["source_file_index"] != source_index):
                    raise ValueError("trigger scaler source/cardinality mismatch")
                for field in BRANCHES:
                    if not np.array_equal(actual[field], batch[field]):
                        raise ValueError(f"trigger scaler copy differs in {field}")
                offset += size
    if offset != out.num_entries:
        raise ValueError("unexpected extra trigger scaler rows")
    if "events" in output_root:
        validate_event_references(inputs, output_root)


def validate_event_references(inputs, output_root):
    import uproot
    keys = ("source_file_index", "event_id_hi", "event_id_lo")
    if REFERENCE not in output_root["events"].keys():
        raise ValueError("collaborator events omitted scaler references")
    actual = output_root["events"].arrays([*keys, REFERENCE], library="np")
    for source_index, path in enumerate(inputs):
        with uproot.open(path) as original:
            event_tree = original["ReplayFoundationV1/RJEventV1"]
            expected = (event_tree[REFERENCE].array(library="np") if REFERENCE in event_tree
                        else np.zeros(event_tree.num_entries, dtype="uint64"))
            if not np.array_equal(actual[REFERENCE][actual["source_file_index"] == source_index], expected):
                raise ValueError("collaborator event scaler reference differs from source")
    lookup = {tuple(int(x) for x in identity): int(ref) for identity, ref in
              zip(zip(*(actual[key] for key in keys)), actual[REFERENCE])}
    normalized = ("storage_layout" in output_root and
                  str(output_root["storage_layout"]) == "normalized_v1")
    for name in ("eventTree", "photons", "jets", "photonJets"):
        if name not in output_root:
            continue
        branches = list(keys) if normalized else [*keys, REFERENCE]
        for rows in output_root[name].iterate(branches, step_size=STEP, library="np"):
            expected = np.asarray([lookup[tuple(int(x) for x in identity)] for identity in
                                   zip(*(rows[key] for key in keys))], dtype="uint64")
            if not normalized and not np.array_equal(rows[REFERENCE], expected):
                raise ValueError(f"derived event scaler join differs: {name}")


def _input_has_scalers(path):
    import uproot
    with uproot.open(path) as root:
        return TREE in root and root[TREE].num_entries > 0
