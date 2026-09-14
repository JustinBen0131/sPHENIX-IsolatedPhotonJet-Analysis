"""Direct GL1 contract: lossless event decisions and source-scoped RUN metadata."""
import numpy as np

EVENT_SCHEMA = {
    "trigger_capture_version": "int32", "trigger_decision_state": "int32",
    "trigger_decision_available": "uint32", "trigger_packet_version": "int32",
    "trigger_packet_status": "uint64",
}
WORDS = ("trigger_bits", "live_trigger_bits", "scaled_trigger_bits")
RUN_TREE = "RJTriggerRunInfoV1"
RUN_SCHEMA = {
    "source_file_index": "int32", "run": "int32", "segment": "int32",
    "input_uri_hash": "string", "source_manifest_sha256": "string",
    "bit": "int32", "name": "string", "state": "int32",
    "initial_prescale": "int32", "prescale": "float64",
    "run_raw": "uint64", "run_live": "uint64", "run_scaled": "uint64",
}


def event_columns(events):
    n = len(events["event_id_hi"])
    present = set(EVENT_SCHEMA) & set(events)
    if not present:
        # Legacy inputs remain readable, but absent evidence is not valid zero.
        return {key: np.full(n, -1 if key == "trigger_decision_state" else 0, dtype=dtype)
                for key, dtype in EVENT_SCHEMA.items()}
    if present != set(EVENT_SCHEMA) or not set(WORDS).issubset(events):
        raise ValueError("partial direct GL1 decision contract")
    out = {}
    for key, dtype in {**EVENT_SCHEMA, **{name: "uint64" for name in WORDS}}.items():
        value = np.asarray(events[key])
        if dtype == "uint64" and (value.dtype.kind != "u" or value.dtype.itemsize != 8):
            raise ValueError(f"{key} must be a lossless uint64 word")
        if value.shape != (n,):
            raise ValueError(f"invalid direct GL1 field shape: {key}")
        out[key] = value.astype(dtype, copy=False)
    version, state = out["trigger_capture_version"], out["trigger_decision_state"]
    mask, packet = out["trigger_decision_available"], out["trigger_packet_version"]
    if np.any(version != 1) or np.any(~np.isin(state, range(6))) or np.any(mask & ~np.uint32(7)):
        raise ValueError("unknown direct GL1 contract version/state/availability")
    valid = state == 0
    if np.any(valid & ((mask != 7) | ~np.isin(packet, [2, 3]) | (out["trigger_packet_status"] != 0))):
        raise ValueError("valid GL1 state contradicts packet evidence")
    modern = np.isin(packet, [2, 3]) & ((mask & 5) == 5)
    if np.any(modern & (out["trigger_bits"] != out["live_trigger_bits"])):
        raise ValueError("v2/v3 legacy trigger word is not its live alias")
    if np.any(((mask & 1) == 0) & (out["live_trigger_bits"] != 0)) or np.any(
            ((mask & 2) == 0) & (out["scaled_trigger_bits"] != 0)):
        raise ValueError("unavailable GL1 decision word contains a fabricated value")
    if np.any(((mask & 4) == 0) & (out["trigger_bits"] != 0)):
        raise ValueError("unavailable legacy GL1 word contains a fabricated value")
    empty = np.isin(state, [1, 3, 5])
    if np.any(empty & ((mask != 0) | (packet != 0) | (out["trigger_packet_status"] != 0))):
        raise ValueError("missing/error/not-DATA GL1 state contradicts availability")
    if np.any((state == 4) & (out["trigger_packet_status"] == 0)):
        raise ValueError("bad-status GL1 state lacks native bad status")
    if np.any((state == 5) & ((mask != 0) | (packet != 0))):
        raise ValueError("simulation GL1 state carries DATA trigger availability")
    return {key: out[key] for key in EVENT_SCHEMA}


def copy_run_info(root, output, events, source_index):
    contract = event_columns(events)
    new = np.any(contract["trigger_capture_version"] == 1)
    data = new and np.any(contract["trigger_decision_state"] != 5)
    if RUN_TREE not in root:
        if data:
            raise ValueError("new DATA trigger contract omitted native run configuration")
        return 0
    tree = root[RUN_TREE]
    keys = list(RUN_SCHEMA)[1:]
    if not set(keys).issubset(tree.keys()):
        raise ValueError("partial native trigger run configuration")
    rows = tree.arrays(keys, library="np")
    if tree.num_entries != 64 or not np.array_equal(rows["bit"], np.arange(64)):
        raise ValueError("trigger run configuration must enumerate all 64 bits exactly once")
    if np.any(~np.isin(rows["state"], [0, 1, 2, 3])):
        raise ValueError("unknown trigger run-info state")
    for key in ("run", "segment", "input_uri_hash", "source_manifest_sha256"):
        if len(set(rows[key])) != 1:
            raise ValueError("trigger run-info source identity changes within input")
    if len(events["run"]) and np.any(events["run"] != rows["run"][0]):
        raise ValueError("trigger run-info run differs from its events")
    output.extend({"source_file_index": np.full(64, source_index, dtype="int32"), **rows})
    return 64


def validate_copy(inputs, output):
    import uproot
    event_tree = output["events"]
    out_fields = [*EVENT_SCHEMA, *WORDS]
    if not set(out_fields).issubset(event_tree.keys()):
        # Old packages do not retroactively acquire a new contract.
        for path in inputs:
            with uproot.open(path) as source:
                if "trigger_capture_version" in source["ReplayFoundationV1/RJEventV1"]:
                    raise ValueError("collaborator output dropped direct GL1 fields")
        return
    actual = event_tree.arrays(["source_file_index", *out_fields], library="np")
    run_offset = 0
    for index, path in enumerate(inputs):
        with uproot.open(path) as source:
            tree = source["ReplayFoundationV1/RJEventV1"]
            fields = [key for key in ["event_id_hi", *EVENT_SCHEMA, *WORDS] if key in tree.keys()]
            original = tree.arrays(fields, library="np")
            expected = event_columns(original)
            use = actual["source_file_index"] == index
            if np.count_nonzero(use) != tree.num_entries:
                raise ValueError("direct GL1 event cardinality changed")
            for key in out_fields:
                value = (original.get(key, np.zeros(tree.num_entries, dtype="uint64"))
                         if key in WORDS else expected[key])
                if not np.array_equal(actual[key][use], value):
                    raise ValueError(f"direct GL1 event copy differs: {key}")
            if RUN_TREE in source:
                if "triggerRunInfo" not in output:
                    raise ValueError("collaborator output dropped trigger run configuration")
                keys = list(RUN_SCHEMA)[1:]
                expected_run = source[RUN_TREE].arrays(keys, library="np")
                size = source[RUN_TREE].num_entries
                actual_run = output["triggerRunInfo"].arrays(list(RUN_SCHEMA), entry_start=run_offset,
                                                          entry_stop=run_offset+size, library="np")
                if np.any(actual_run["source_file_index"] != index):
                    raise ValueError("trigger run configuration source join differs")
                for key in keys:
                    a, b = actual_run[key], expected_run[key]
                    same = (np.array_equal(a, b, equal_nan=True) if a.dtype.kind == "f"
                            else np.array_equal(a, b))
                    if not same:
                        raise ValueError(f"trigger run configuration differs: {key}")
                run_offset += size
    if "triggerRunInfo" in output and output["triggerRunInfo"].num_entries != run_offset:
        raise ValueError("extra trigger run configuration rows")
    # Expanded output repeats event fields: verify them in each consumer table.
    if str(output["storage_layout"]) != "normalized_v1":
        reference = {(int(s), int(eh), int(el)): i for i, (s, eh, el) in enumerate(zip(
            actual["source_file_index"], event_tree["event_id_hi"].array(library="np"),
            event_tree["event_id_lo"].array(library="np")))}
        for name in ("photons", "jets", "photonJets", "eventTree"):
            for rows in output[name].iterate(["source_file_index", "event_id_hi", "event_id_lo", *out_fields],
                                              step_size=2048, library="np"):
                pos = [reference[(int(s), int(h), int(l))] for s, h, l in zip(
                    rows["source_file_index"], rows["event_id_hi"], rows["event_id_lo"])]
                if any(not np.array_equal(rows[key], actual[key][pos]) for key in out_fields):
                    raise ValueError(f"expanded {name} GL1 event fields differ")
