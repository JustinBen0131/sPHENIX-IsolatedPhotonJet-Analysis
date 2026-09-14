#!/usr/bin/env python3
"""Local, occurrence-keyed AuAu centrality augmentation; never runs reconstruction.

The calibration manifest binds local CDB files to a run/tag and explicit node
names. Values are read from those files, not trusted from copied JSON numbers.
Only the pinned CentralityReco algorithm is supported. A friend is an event
calculation, not acceptance of reused objects, cuts, weights or physics results.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import struct
import tempfile

import awkward as ak
import numpy as np
import uproot

ALGORITHM = "CentralityReco_mbd_q0p5_time25_float_v1"
SUPPORTED_NATIVE_SHA256 = {
    "CentralityReco.cc": "1a05f62781d36ec2800d0140ea0dd91f2eeb2dbd8ff6246c26111d3e06a5cdc4",
    "CentralityReco.h": "74cbaf6a3f83043e9e81edf11809a1b841200871d1348a9a1a6b5a1f7aea68f4",
}
SCALARS = {
    **{name: "int32" for name in (
        "centrality_replay_version", "centrality_pmt_available", "centrality_inputs_valid",
        "centrality_mb_decision", "centrality_mbd_z_valid", "centrality_mbd_event",
        "centrality_mbd_clock", "centrality_mbd_femclock", "centrality_native_bin",
        "centrality_native_valid")},
    **{name: "float64" for name in (
        "centrality_mbd_z", "centrality_selected_charge", "centrality_native_centile")},
}
ARRAYS = {
    "centrality_pmt_id": "var * int32", "centrality_pmt_charge": "var * float64",
    "centrality_pmt_time": "var * float64", "centrality_pmt_valid": "var * int32",
    "centrality_pmt_selected": "var * int32",
}
KEYS = ("source_occurrence_id_hi", "source_occurrence_id_lo", "event_id_hi", "event_id_lo")

def event_layout(source, tree):
    """Compact IDs are source-local; collaborator files add their source index."""
    with uproot.open(source) as root:
        if tree is None:
            tree = "eventTree" if "eventTree" in root else "events"
        names = set(root[tree].keys())
        keys = KEYS + (("source_file_index",) if "source_file_index" in names else ())
        return tree, keys
NODES = ("pmt", "mbd_out", "minimum_bias", "centrality")
STATES = {0: "VALID", 1: "NON_MB", 2: "MISSING_OR_INVALID_INPUT",
          3: "VERTEX_OUTSIDE_PAYLOAD", 4: "NO_DIVISION_MATCH", 5: "NONFINITE_ESTIMATOR"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as src:
        for block in iter(lambda: src.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def capture_columns(events):
    """Preserve a complete capture group; legacy sources get explicit absence."""
    n = len(events["event_id_hi"])
    present = (set(SCALARS) | set(ARRAYS)) & set(events)
    if present and present != set(SCALARS) | set(ARRAYS):
        raise ValueError("partial centrality replay capture group")
    scalar = {k: np.asarray(events[k], dtype=v) if present else
              np.full(n, -1 if v == "int32" else math.nan, dtype=v) for k, v in SCALARS.items()}
    if not present:
        scalar["centrality_replay_version"][:] = 0
    arrays = {}
    for k, v in ARRAYS.items():
        dtype = v.split()[-1]
        rows = [np.asarray(row, dtype=dtype) for row in events[k]] if present else []
        if any(row.ndim != 1 for row in rows):
            raise ValueError("PMT arrays must have one channel dimension")
        counts = np.asarray([len(row) for row in rows], dtype="int64") if present else np.zeros(n, dtype="int64")
        flat = np.concatenate(rows) if rows else np.array([], dtype=dtype)
        arrays[k] = ak.unflatten(flat, counts)
    if present:
        if np.any(scalar["centrality_replay_version"] != 1):
            raise ValueError("unsupported centrality replay version")
        for i in range(n):
            lengths = {len(arrays[k][i]) for k in ARRAYS}
            available = scalar["centrality_pmt_available"][i]
            if available not in (-1, 0, 1) or lengths != ({128} if available == 1 else {0}):
                raise ValueError("malformed centrality PMT capture")
            if scalar["centrality_mb_decision"][i] not in (-1, 0, 1):
                raise ValueError("malformed MB decision")
    return scalar, arrays


def _cdb_columns(spec, names):
    path = Path(spec["path"])
    if digest(path) != spec["sha256"]:
        raise ValueError(f"CDB payload hash mismatch: {path}")
    with uproot.open(path) as src:
        rows = src["Multiple"].arrays(["IID", *names], library="np")
    ids = [int(x) for x in rows["IID"]]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate CDB IID")
    # CDB uses a type prefix in every field name: F, D, I, g.
    for name in names:
        kind = {"F": "f4", "D": "f8", "I": "i4"}[name[0]]
        if rows[name].dtype.kind != np.dtype(kind).kind or rows[name].dtype.itemsize != np.dtype(kind).itemsize:
            raise ValueError(f"CDB field type mismatch: {name}")
    return {i: {name: rows[name][j].item() for name in names} for j, i in enumerate(ids)}


def load_calibrations(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    if manifest.get("schema") != "AuAuCentralityCalibrationManifestV1" or manifest.get("algorithm") != ALGORITHM:
        raise ValueError("unsupported calibration/algorithm manifest")
    software = manifest["software"]
    if not software.get("release"):
        raise ValueError("software release is required")
    for name in ("CentralityReco.cc", "CentralityReco.h"):
        spec = software[name]
        if digest(spec["path"]) != spec["sha256"] or spec["sha256"] != SUPPORTED_NATIVE_SHA256[name]:
            raise ValueError("native source hash mismatch")
    result = {}
    for record in manifest["records"]:
        run = record["run"]
        if type(run) is not int or run <= 0 or run in result or not record.get("global_tag"):
            raise ValueError("invalid or repeated calibration run/tag")
        if not all(isinstance(record["nodes"].get(k), str) and record["nodes"][k] for k in NODES):
            raise ValueError("exact centrality node bindings required")
        domains = record["payloads"]
        divs = _cdb_columns(domains["Centrality"], ["Fcentralitydiv"])
        scale = _cdb_columns(domains["CentralityScale"], ["Dcentralityscale"])
        vtx = _cdb_columns(domains["CentralityVertexScale"],
                           ["Invertexbins", "Dscale", "Dlow_vertex", "Dhigh_vertex"])
        divisions = np.asarray([divs[i]["Fcentralitydiv"] for i in range(100)], dtype="float32")
        run_scale = float(scale[0]["Dcentralityscale"])
        count = vtx[0]["Invertexbins"]
        if type(count) is not int or not 0 < count <= 10000:
            raise ValueError("invalid vertex bin count")
        ranges = [[float(np.float32(vtx[i][k])) for k in
                   ("Dlow_vertex", "Dhigh_vertex", "Dscale")] for i in range(count)]
        if (not np.all(np.isfinite(divisions)) or np.any(np.diff(divisions) > 0)
                or not math.isfinite(run_scale) or run_scale <= 0):
            raise ValueError("invalid divisions or run scale")
        for i, (low, high, factor) in enumerate(ranges):
            if not all(map(math.isfinite, (low, high, factor))) or low >= high or factor <= 0:
                raise ValueError("invalid vertex scale interval")
            if i and low < ranges[i - 1][1]:
                raise ValueError("overlapping or unordered vertex intervals")
        result[run] = {**record, "divisions": divisions, "run_scale": run_scale, "vertex_ranges": ranges}
    if not result:
        raise ValueError("empty calibration manifest")
    return manifest, result


def calculate(charge, time, valid, z, mb_decision, calibration):
    """Mirror native float arithmetic for finite inputs; invalids stay explicit.

    Never replace the per-channel calculation by sum(q)*scale: rounding differs.
    NaN times pass the old native `fabs(t)>25` rejection accidentally. We refuse
    those events instead of claiming a supported, physically valid percentile.
    """
    result = dict(state=2, bin=0, centile=math.nan, percentile=math.nan,
                  selected_charge=math.nan, scaled_charge=math.nan,
                  vertex_scale=math.nan, run_scale=calibration["run_scale"])
    if mb_decision == 0:
        result["state"] = 1
        return result
    q, t = np.asarray(charge, dtype="float32"), np.asarray(time, dtype="float32")
    if (mb_decision != 1 or len(q) != 128 or len(t) != 128 or len(valid) != 128
            or not np.all(np.asarray(valid) == 1) or not np.all(np.isfinite(q))
            or not np.all(np.isfinite(t)) or not math.isfinite(z)):
        return result
    z = np.float32(z)
    factor = next((np.float32(v) for lo, hi, v in calibration["vertex_ranges"] if z > lo and z <= hi), None)
    if factor is None:
        result["state"], result["vertex_scale"] = 3, 0.0
        return result
    selected = np.float32(0)
    scaled = np.float32(0)
    with np.errstate(over="ignore", invalid="ignore"):
        for a, b in zip(q, t):
            if a < np.float32(0.5) or abs(float(b)) > 25:
                continue
            selected = np.float32(selected + a)
            # float * float first, then double scale and double addition.
            scaled = np.float32(float(scaled) + float(np.float32(a * factor)) * calibration["run_scale"])
    result.update(selected_charge=float(selected), scaled_charge=float(scaled), vertex_scale=float(factor))
    if not math.isfinite(float(scaled)):
        result["state"] = 5
        return result
    matches = np.flatnonzero(calibration["divisions"] < scaled)
    if not len(matches):
        result["state"] = 4
        return result
    index = int(matches[0]) + 1
    result.update(state=0, bin=index, centile=float(np.float32(index) / np.float32(100)), percentile=float(index))
    return result


def _parity(event, calibration):
    """Replay under the as-produced calibration and compare to the native witness.

    This is the only check that proves the offline arithmetic reproduces what
    CentralityReco actually wrote, so it needs the calibration the producer job
    used, not the calibration we are migrating to.
    """
    replayed = calculate(event["centrality_pmt_charge"], event["centrality_pmt_time"],
                         event["centrality_pmt_valid"], event["centrality_mbd_z"],
                         event["centrality_mb_decision"], calibration)
    native_bin = int(event["centrality_native_bin"])
    native_centile = float(event["centrality_native_centile"])
    agrees = (replayed["state"] == 0 and replayed["bin"] == native_bin
              and replayed["centile"] == native_centile)
    return agrees, replayed, native_bin, native_centile


def augment(source, manifest_path, binding_path, output, *, tree=None, step_size=10000,
            parity_manifest=None, allow_vertex_domain_divergence=False):
    """Write a new directory containing a ROOT friend plus source/hash receipt.

    Every event, including invalid and non-MB rows, gets exactly one keyed row.
    A disk-backed unique index bounds join memory. Inputs are rehashed at exit.
    No implicit friend-by-entry attachment and no overwriting existing output.

    `parity_manifest` is the calibration the producer job actually ran with. It
    is optional only because it does not always exist; without it the receipt
    records that native-versus-offline parity was NOT checked, and no consumer
    may read the friend as a verified reproduction of CentralityReco.
    """
    source, manifest_path, binding_path, output = map(Path, (source, manifest_path, binding_path, output))
    if output.exists():
        raise ValueError("output already exists; reuse or choose a new version")
    if step_size <= 0:
        raise ValueError("step size must be positive")
    _, calibrations = load_calibrations(manifest_path)
    parity_calibrations = None
    parity_hash = None
    if parity_manifest is not None:
        parity_manifest = Path(parity_manifest)
        parity_hash = digest(parity_manifest)
        _, parity_calibrations = load_calibrations(parity_manifest)
        if parity_hash == digest(manifest_path):
            raise ValueError("parity manifest must be the as-produced calibration, "
                             "not the calibration being migrated to")
    source_hash, manifest_hash, binding_hash = map(digest, (source, manifest_path, binding_path))
    binding = json.loads(binding_path.read_text())
    if binding.get("schema") != "AuAuCentralityInputBindingV1" or binding.get("source_sha256") != source_hash:
        raise ValueError("input binding/source hash mismatch")
    if binding.get("system") != "AuAu" or binding.get("sample_kind") not in ("data", "embedded"):
        raise ValueError("AuAu DATA or explicitly bound embedding required")
    if binding.get("input_measurement_definition") != ALGORITHM:
        raise ValueError("unproved input measurement/algorithm definition")
    for k in NODES:
        if not binding.get("nodes", {}).get(k):
            raise ValueError("source node provenance missing")
    if binding["sample_kind"] == "embedded" and not binding.get("background_identity_verified"):
        raise ValueError("embedding background-event identity is unverified")
    run_branch = binding.get("calibration_run_branch")
    if not isinstance(run_branch, str) or not run_branch:
        raise ValueError("explicit DATA/background calibration run branch required")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()  # Exclusive: a failed invocation leaves reviewable evidence.
    counts = Counter()
    rows_written = 0
    tree, key_fields = event_layout(source, tree)
    branches = [*key_fields, run_branch, "centrality", *SCALARS, *ARRAYS]
    schema = {**{k: "uint64" for k in key_fields}, "source_entry": "uint64", "calibration_run": "int32",
              "centrality_original": "float64", "native_centile_original": "float64",
              "native_bin_original": "int32", "native_valid_original": "int32",
              "centrality_state": "int32", "centrality_bin_new": "int32",
              **{k: "float64" for k in ("centrality_centile_new", "centrality_percentile_new",
                  "selected_charge", "scaled_charge", "vertex_scale", "run_scale")}}
    root_path = output / "centrality_friend.root"
    receipt = {"schema": "AuAuCentralityFriendReceiptV1", "status": "INCOMPLETE",
               "source": str(source.resolve()), "source_sha256": source_hash, "source_tree": tree,
               "manifest": str(manifest_path.resolve()), "manifest_sha256": manifest_hash,
               "binding": str(binding_path.resolve()), "binding_sha256": binding_hash,
               "key_fields": list(key_fields), "state_codes": STATES,
               "native_parity": "NOT_CHECKED" if parity_calibrations is None else "PENDING",
               "parity_manifest": None if parity_manifest is None else str(parity_manifest.resolve()),
               "parity_manifest_sha256": parity_hash,
               "downstream_analysis_requires_regeneration": True, "production_accepted": False}
    parity_counts = Counter()
    try:
        with tempfile.TemporaryDirectory(prefix="centrality-keys-") as tmp, \
                sqlite3.connect(str(Path(tmp) / "keys.sqlite")) as keys, \
                uproot.open(source) as src, uproot.create(root_path) as dest:
            keys.execute("CREATE TABLE keys (identity TEXT PRIMARY KEY)")
            event_tree = src[tree]
            missing = set(branches) - set(event_tree.keys())
            if missing:
                raise ValueError(f"source lacks replay branches: {sorted(missing)}")
            friend = dest.mktree("centralityFriend", schema)
            for batch in event_tree.iterate(branches, library="ak", step_size=step_size):
                columns = {k: [] for k in schema}
                for event in ak.to_list(batch):
                    identity = tuple(int(event[k]) for k in key_fields)
                    if any(v < 0 for v in identity) or identity[:2] == (0, 0) or identity[2:4] == (0, 0):
                        raise ValueError("missing occurrence/event identity")
                    keys.execute("INSERT INTO keys VALUES (?)", (":".join(map(str, identity)),))
                    run = int(event[run_branch])
                    if run != event[run_branch] or run not in calibrations:
                        raise ValueError(f"no exact calibration binding for run {run}")
                    calibration = calibrations[run]
                    if calibration["nodes"] != binding["nodes"]:
                        raise ValueError("source/calibration node identity mismatch")
                    if event["centrality_replay_version"] != 1:
                        raise ValueError("input was not captured for centrality replay")
                    pmt_present = event["centrality_pmt_available"] == 1
                    lengths = {len(event[k]) for k in ARRAYS}
                    if lengths != ({128} if pmt_present else {0}):
                        raise ValueError("malformed PMT arrays")
                    if pmt_present and event["centrality_pmt_id"] != list(range(128)):
                        raise ValueError("PMT channel identity/order mismatch")
                    if event["centrality_mb_decision"] not in (-1, 0, 1):
                        raise ValueError("invalid MB decision")
                    if event["centrality_mbd_z_valid"] not in (0, 1) or bool(event["centrality_mbd_z_valid"]) != math.isfinite(event["centrality_mbd_z"]):
                        raise ValueError("MBD vertex validity mismatch")
                    if pmt_present:
                        q = np.asarray(event["centrality_pmt_charge"], dtype="float32")
                        t = np.asarray(event["centrality_pmt_time"], dtype="float32")
                        valid = np.asarray(event["centrality_pmt_valid"])
                        if not np.all(np.isin(valid, (0, 1))) or np.any((valid == 1) & ~(np.isfinite(q) & np.isfinite(t))):
                            raise ValueError("invalid PMT validity witness")
                        selected = (valid == 1) & (q >= .5) & (np.abs(t) <= 25)
                        if not np.array_equal(selected.astype("int32"), event["centrality_pmt_selected"]):
                            raise ValueError("captured PMT selection disagrees with algorithm")
                    r = calculate(event["centrality_pmt_charge"], event["centrality_pmt_time"],
                                  event["centrality_pmt_valid"], event["centrality_mbd_z"],
                                  event["centrality_mb_decision"], calibration)
                    if r["state"] == 0 and event["centrality_inputs_valid"] != 1:
                        raise ValueError("valid replay contradicts captured input validity")
                    if math.isfinite(r["selected_charge"]) and r["selected_charge"] != event["centrality_selected_charge"]:
                        raise ValueError("captured selected charge disagrees with replay")
                    counts[STATES[r["state"]]] += 1
                    if parity_calibrations is not None:
                        if run not in parity_calibrations:
                            raise ValueError(f"no as-produced calibration binding for run {run}")
                        if int(event["centrality_native_valid"]) == 1:
                            agrees, replayed, native_bin, native_centile = _parity(
                                event, parity_calibrations[run])
                            parity_counts["checked"] += 1
                            detail = (f"{identity}: native bin={native_bin} "
                                      f"centile={native_centile!r}; replay state="
                                      f"{STATES[replayed['state']]} bin={replayed['bin']} "
                                      f"centile={replayed['centile']!r}")
                            if not agrees and replayed["state"] == 3:
                                # The as-produced payload did not cover this MBD
                                # vertex, so getVertexScale() returned 0 and the
                                # native pass still wrote a bin from a zero
                                # estimator. That is a real, documented native
                                # divergence, not an arithmetic disagreement.
                                if not allow_vertex_domain_divergence:
                                    raise ValueError(
                                        "native/offline centrality parity failed on an "
                                        "uncovered vertex domain at " + detail)
                                parity_counts["vertex_domain_divergence"] += 1
                            elif not agrees:
                                raise ValueError(
                                    "native/offline centrality parity failed at " + detail)
                            else:
                                parity_counts["agree"] += 1
                        else:
                            parity_counts["no_native_witness"] += 1
                    row = {**dict(zip(key_fields, identity)), "source_entry": rows_written, "calibration_run": run,
                           "centrality_original": event["centrality"],
                           "native_centile_original": event["centrality_native_centile"],
                           "native_bin_original": event["centrality_native_bin"],
                           "native_valid_original": event["centrality_native_valid"],
                           "centrality_state": r["state"], "centrality_bin_new": r["bin"],
                           "centrality_centile_new": r["centile"], "centrality_percentile_new": r["percentile"],
                           **{k: r[k] for k in ("selected_charge", "scaled_charge", "vertex_scale", "run_scale")}}
                    for k in columns:
                        columns[k].append(row[k])
                    rows_written += 1
                friend.extend({k: np.asarray(v, dtype=schema[k]) for k, v in columns.items()})
            if rows_written != event_tree.num_entries:
                raise ValueError("friend/source row count mismatch")
            if rows_written == 0:
                raise ValueError("empty source cannot certify replay")
        if [digest(p) for p in (source, manifest_path, binding_path)] != [source_hash, manifest_hash, binding_hash]:
            raise ValueError("inputs changed during augmentation")
        load_calibrations(manifest_path)  # Recheck payload and software hashes.
        with uproot.open(root_path) as readback:
            if readback["centralityFriend"].num_entries != rows_written:
                raise ValueError("ROOT friend readback row count mismatch")
        if parity_calibrations is not None:
            load_calibrations(parity_manifest)  # Recheck the as-produced payloads too.
            if not parity_counts["checked"]:
                raise ValueError("parity was requested but no event carried a native witness")
            receipt["native_parity"] = ("PASS_WITH_VERTEX_DOMAIN_DIVERGENCE"
                                        if parity_counts["vertex_domain_divergence"] else "PASS")
            receipt["parity_counts"] = dict(parity_counts)
        receipt.update(status="pass", rows=rows_written, states=dict(counts), friend_sha256=digest(root_path))
    except Exception as exc:
        receipt.update(status="FAIL", error=f"{type(exc).__name__}: {exc}", rows_written=rows_written)
        if parity_calibrations is not None:
            receipt["native_parity"] = "FAIL"
            receipt["parity_counts"] = dict(parity_counts)
        (output / "RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
        raise
    (output / "RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def verify(source, output, *, tree=None, step_size=10000):
    """Prove the friend joins its source by identity, not by entry alignment.

    Equal row counts are not a join.  Every source event is looked up in a
    disk-backed index of the friend's own key columns, the claimed source entry
    must equal the real source entry, and every preserved original column is
    compared against an independent read of the source.  The report also states
    how many events a consumer would move across the 0-80% analysis window, so
    a centrality migration is a measured change and not an assumption.
    """
    source, output = Path(source), Path(output)
    receipt = json.loads((output / "RECEIPT.json").read_text())
    if receipt.get("schema") != "AuAuCentralityFriendReceiptV1" or receipt.get("status") != "pass":
        raise ValueError("friend receipt is absent, unsupported or not a passing receipt")
    root_path = output / "centrality_friend.root"
    if digest(root_path) != receipt["friend_sha256"]:
        raise ValueError("friend ROOT file changed since its receipt")
    if digest(source) != receipt["source_sha256"]:
        raise ValueError("source changed since the friend was written")
    if str(source.resolve()) != receipt["source"]:
        raise ValueError("this friend was not written for this source")
    tree, key_fields = event_layout(source, tree or receipt["source_tree"])
    if tree != receipt["source_tree"] or list(key_fields) != receipt["key_fields"]:
        raise ValueError("friend/source event tree or identity key contract differs")
    key_width = len(key_fields)

    original = {"centrality": "centrality_original",
                "centrality_native_centile": "native_centile_original",
                "centrality_native_bin": "native_bin_original",
                "centrality_native_valid": "native_valid_original"}
    report = {"schema": "AuAuCentralityFriendVerificationV1", "rows": 0,
              "states": Counter(), "boundary": Counter(), "bin_changed": 0}
    with tempfile.TemporaryDirectory(prefix="centrality-verify-") as tmp, \
            sqlite3.connect(str(Path(tmp) / "friend.sqlite")) as index, \
            uproot.open(source) as src, uproot.open(root_path) as friend_file:
        friend = friend_file["centralityFriend"]
        if friend.num_entries != src[tree].num_entries or friend.num_entries != receipt["rows"]:
            raise ValueError("friend/source/receipt row counts differ")
        # Floating point columns travel as raw IEEE bytes: sqlite silently maps
        # a NaN REAL to NULL, and the invalid rows are exactly the ones whose
        # preserved originals must stay NaN.
        index.execute("CREATE TABLE rows (identity TEXT PRIMARY KEY, entry INTEGER, "
                      "state INTEGER, bin INTEGER, centile BLOB, "
                      "centrality BLOB, centile0 BLOB, bin0 INTEGER, valid0 INTEGER)")
        columns = [*key_fields, "source_entry", "centrality_state", "centrality_bin_new",
                   "centrality_centile_new", *original.values()]
        for batch in friend.iterate(columns, library="np", step_size=step_size):
            index.executemany(
                "INSERT INTO rows VALUES (?,?,?,?,?,?,?,?,?)",
                [(":".join(str(int(row[i])) for i in range(key_width)),
                  int(row[key_width]), int(row[key_width+1]), int(row[key_width+2]),
                  struct.pack("<d", float(row[key_width+3])),
                  struct.pack("<d", float(row[key_width+4])), struct.pack("<d", float(row[key_width+5])),
                  int(row[key_width+6]), int(row[key_width+7]))
                 for row in zip(*(batch[name] for name in columns))])
        entry = 0
        source_columns = [*key_fields, *original]
        for batch in src[tree].iterate(source_columns, library="np", step_size=step_size):
            for values in zip(*(batch[name] for name in source_columns)):
                identity = ":".join(str(int(v)) for v in values[:key_width])
                found = index.execute(
                    "SELECT entry, state, bin, centile, centrality, centile0, bin0, valid0 "
                    "FROM rows WHERE identity = ?", (identity,)).fetchone()
                if found is None:
                    raise ValueError(f"source event {identity} has no friend row")
                if found[0] != entry:
                    raise ValueError(
                        f"friend row for {identity} claims source entry {found[0]}, not {entry}")
                observed = dict(zip(original, values[key_width:]))
                for offset, name in enumerate(original):
                    stored = found[4 + offset]
                    claimed = struct.unpack("<d", stored)[0] if isinstance(stored, bytes) \
                        else float(stored)
                    actual = float(observed[name])
                    if not (claimed == actual or (math.isnan(claimed) and math.isnan(actual))):
                        raise ValueError(f"preserved original differs for {identity}: {name}")
                state, new_bin = found[1], found[2]
                new_centile = struct.unpack("<d", found[3])[0]
                if state == 0:
                    if not 1 <= new_bin <= 100 or new_centile != float(np.float32(
                            np.float32(new_bin) / np.float32(100))):
                        raise ValueError(f"valid friend row is internally inconsistent: {identity}")
                elif new_bin != 0 or not math.isnan(new_centile):
                    raise ValueError(f"invalid friend row carries a value: {identity}")
                report["states"][STATES[state]] += 1
                old_bin = int(observed["centrality_native_bin"])
                if int(observed["centrality_native_valid"]) == 1:
                    if new_bin != old_bin:
                        report["bin_changed"] += 1
                    was_in, now_in = old_bin <= 80, state == 0 and new_bin <= 80
                    if was_in and not now_in:
                        report["boundary"]["left_analysis_window"] += 1
                    elif now_in and not was_in:
                        report["boundary"]["entered_analysis_window"] += 1
                entry += 1
        if entry != friend.num_entries:
            raise ValueError("source iteration did not cover every friend row")
    report["rows"] = entry
    report["states"] = dict(report["states"])
    report["boundary"] = dict(report["boundary"])
    report["friend_sha256"] = receipt["friend_sha256"]
    report["native_parity"] = receipt.get("native_parity", "NOT_CHECKED")
    report["safe_to_consume_as_native_reproduction"] = report["native_parity"] == "PASS"
    if report["native_parity"] == "PASS_WITH_VERTEX_DOMAIN_DIVERGENCE":
        report["divergence"] = (
            "events whose MBD vertex the as-produced payload did not cover still "
            "carry a native bin computed from a zero estimator; the replay refuses "
            "to reproduce that value and marks them VERTEX_OUTSIDE_PAYLOAD")
    (output / "VERIFICATION.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--calibration-manifest", type=Path)
    parser.add_argument("--input-binding", type=Path)
    parser.add_argument("--as-produced-calibration-manifest", type=Path,
                        help="calibration the producer job ran with; enables native parity")
    parser.add_argument("--verify-only", action="store_true",
                        help="re-verify an existing friend directory against its source")
    parser.add_argument("--allow-vertex-domain-divergence", action="store_true",
                        help="accept native bins written from a zero estimator when the "
                             "as-produced vertex payload did not cover the MBD vertex")
    parser.add_argument("--tree", help="event tree; auto-select eventTree or normalized events")
    parser.add_argument("--step-size", type=int, default=10000)
    args = parser.parse_args()
    if args.verify_only:
        result = verify(args.source, args.output, tree=args.tree, step_size=args.step_size)
    else:
        if not args.calibration_manifest or not args.input_binding:
            parser.error("--calibration-manifest and --input-binding are required")
        augment(args.source, args.calibration_manifest, args.input_binding, args.output,
                tree=args.tree, step_size=args.step_size,
                parity_manifest=args.as_produced_calibration_manifest,
                allow_vertex_domain_divergence=args.allow_vertex_domain_divergence)
        result = verify(args.source, args.output, tree=args.tree, step_size=args.step_size)
    print(json.dumps(result, indent=2))
