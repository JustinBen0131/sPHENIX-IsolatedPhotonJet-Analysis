#!/usr/bin/env python3
"""Compare normalized canonical content by exact row identity, not ROOT bytes.

All tables, fields, ROOT types, integer/enum states and jagged shapes must
agree. Floating values use explicit tolerances, preserving NaN versus signed
infinity. Exit 0 means the implemented content comparison passed; it is not
scientific acceptance or a replacement for the reviewed migration ledger.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import uproot

IDENTITY_COLUMNS = {
    "Sources": ("source_hi", "source_lo"),
    "UpstreamRejectedEvents": ("source_hi", "source_lo", "physical_event_sequence", "source_entry"),
    "Events": ("event_hi", "event_lo"),
    "Photons": ("event_hi", "event_lo", "photon_hi", "photon_lo"),
    "PhotonShowerViews": ("event_hi", "event_lo", "photon_hi", "photon_lo", "definition"),
    "PhotonCells": ("event_hi", "event_lo", "photon_hi", "photon_lo", "subsystem", "tower_key"),
    "Isolation": ("event_hi", "event_lo", "photon_hi", "photon_lo", "radius", "method"),
    "IsolationConstituents": ("event_hi", "event_lo", "photon_hi", "photon_lo", "radius", "source", "subsystem", "native_key"),
    "Jets": ("event_hi", "event_lo", "jet_hi", "jet_lo"),
    "PhotonJetPairs": ("event_hi", "event_lo", "pair_hi", "pair_lo"),
    "TruthVertices": ("event_hi", "event_lo", "vertex_id"),
    "TruthPhotons": ("event_hi", "event_lo", "truth_photon_hi", "truth_photon_lo"),
    "TruthJets": ("event_hi", "event_lo", "truth_jet_hi", "truth_jet_lo"),
    "PhotonTruthLinks": ("link_hi", "link_lo"),
    "JetTruthLinks": ("link_hi", "link_lo"),
    "WeightComponents": ("event_hi", "event_lo", "component_type"),
    "TriggerScalers": ("source_hi", "source_lo", "snapshot_id"),
    "TriggerRunInfo": ("source_hi", "source_lo", "run", "bit"),
}

DEFAULT_RTOL = 1e-6
DEFAULT_ATOL = 1e-9

# Producer/build identities legitimately differ between converted and direct
# files. They are reported verbatim, not used as scientific equality tests.
PROVENANCE_ONLY = {"producer_git_commit", "producer_source_sha256", "producer_library_sha256",
                   "macro_sha256", "production_tag", "input_uri_sha256", "input_file_sha256"}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""): h.update(block)
    return h.hexdigest()


def read(path):
    tables, types = {}, {}
    with uproot.open(path) as root:
        for name in root.keys(cycle=False):
            obj = root[name]
            if hasattr(obj, "arrays"):
                tables[name] = obj.arrays(library="np")
                types[name] = obj.typenames()
        metadata = dict(line.split("=", 1) for line in str(root["metadata"]).splitlines() if "=" in line)
        completion = dict(line.split("=", 1) for line in str(root["completion"]).splitlines() if "=" in line)
        if metadata.get("contract") != "PhotonJetTrees" or completion.get("completion_status") != "complete":
            raise ValueError(f"{path}: incomplete or wrong-contract product")
    return tables, types, metadata, completion


def keys_of(table, columns):
    names = IDENTITY_COLUMNS[table]
    if not all(name in columns for name in names):
        raise ValueError(f"{table}: missing required identity columns")
    n = len(columns[names[0]])
    if any(len(v) != n for v in columns.values()):
        raise ValueError(f"{table}: inconsistent column lengths")
    keys = [tuple(columns[c][i].item() if hasattr(columns[c][i], "item") else columns[c][i] for c in names)
            for i in range(n)]
    if any(isinstance(v, float) and not math.isfinite(v) for key in keys for v in key):
        raise ValueError(f"{table}: non-finite identity")
    if len(set(keys)) != n:
        raise ValueError(f"{table}: duplicate identity")
    return keys


def equal_value(left, right, rtol, atol):
    """Recursive arrays; never cast uint64 identities/counters through float."""
    a, b = np.asarray(left), np.asarray(right)
    if a.shape != b.shape: return False
    if a.dtype.kind == "O" or b.dtype.kind == "O":
        if a.ndim == 0 and b.ndim == 0:
            return type(left) is type(right) and left == right
        return all(equal_value(x, y, rtol, atol) for x, y in zip(a, b))
    if a.dtype.kind in "biu" or b.dtype.kind in "biu":
        return a.dtype.kind == b.dtype.kind and bool(np.array_equal(a, b))
    if a.dtype.kind in "SU" or b.dtype.kind in "SU":
        return bool(np.array_equal(a, b))
    # NaN equals NaN; +inf only equals +inf and -inf only equals -inf.
    same_special = (np.isnan(a) & np.isnan(b)) | (np.isposinf(a) & np.isposinf(b)) | (np.isneginf(a) & np.isneginf(b))
    finite = np.isfinite(a) & np.isfinite(b)
    return bool(np.all(same_special | (finite & np.isclose(a, b, rtol=rtol, atol=atol))))


def compare_table(name, left, right, left_types, right_types, rtol, atol):
    lk, rk = keys_of(name, left), keys_of(name, right)
    li, ri = {k: i for i, k in enumerate(lk)}, {k: i for i, k in enumerate(rk)}
    common = [k for k in lk if k in ri]
    result = {"left_rows": len(lk), "right_rows": len(rk), "left_only": len(set(lk)-set(rk)),
              "right_only": len(set(rk)-set(lk)), "left_only_fields": sorted(set(left)-set(right)),
              "right_only_fields": sorted(set(right)-set(left)), "fields": {}}
    for field in sorted(set(left) & set(right)):
        differences = [key for key in common if not equal_value(left[field][li[key]], right[field][ri[key]], rtol, atol)]
        result["fields"][field] = {"different_rows": len(differences), "examples": differences[:5],
                                  "left_type": left_types[field], "right_type": right_types[field],
                                  "type_equal": left_types[field] == right_types[field]}
    result["passed"] = not any(result[k] for k in ("left_only", "right_only", "left_only_fields", "right_only_fields")) and all(
        row["type_equal"] and row["different_rows"] == 0 for row in result["fields"].values())
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", required=True, type=Path)
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    parser.add_argument("--atol", type=float, default=DEFAULT_ATOL)
    args = parser.parse_args(argv)
    if not all(math.isfinite(v) and v >= 0 for v in (args.rtol, args.atol)):
        parser.error("tolerances must be finite and nonnegative")
    hashes = {"legacy": sha256(args.legacy), "canonical": sha256(args.canonical)}
    left, lt, lm, lc = read(args.legacy)
    right, rt, rm, rc = read(args.canonical)
    report = {"schema": "CanonicalContentComparisonV2", "file_sha256": hashes,
              "file_bytes_equal": hashes["legacy"] == hashes["canonical"],
              "tolerances": {"rtol": args.rtol, "atol": args.atol}, "results": {},
              "left_only_tables": sorted(set(left)-set(right)), "right_only_tables": sorted(set(right)-set(left)),
              "unsupported_tables": sorted((set(left)|set(right))-set(IDENTITY_COLUMNS)),
              "metadata_left": lm, "metadata_right": rm,
              "metadata_differences": sorted(k for k in (set(lm)|set(rm))-PROVENANCE_ONLY if lm.get(k) != rm.get(k)),
              "completion_equal": lc == rc, "scientific_acceptance": False}
    for name in sorted(set(left) & set(right) & set(IDENTITY_COLUMNS)):
        try:
            report["results"][name] = compare_table(name, left[name], right[name], lt[name], rt[name], args.rtol, args.atol)
        except (ValueError, TypeError, KeyError) as error:
            report["results"][name] = {"passed": False, "error": str(error)}
    report["content_passed"] = (all(name in left and name in right for name in ("Sources", "Events", "Photons"))
        and not any(report[k] for k in ("left_only_tables", "right_only_tables", "unsupported_tables", "metadata_differences"))
        and report["completion_equal"] and all(r["passed"] for r in report["results"].values()))
    if sha256(args.legacy) != hashes["legacy"] or sha256(args.canonical) != hashes["canonical"]:
        raise ValueError("input changed during comparison")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(f"content_passed={report['content_passed']}; scientific acceptance remains a separate closure gate")
    return 0 if report["content_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
