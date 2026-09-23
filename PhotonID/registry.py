"""Declarative package binding and completed TreeProduction campaign discovery.

This module owns operational validation, never reconstruction or class labels.
Paths in a manifest/registry are relative to that document, not the caller's cwd.
"""

from pathlib import Path
import hashlib
import json
import re
import yaml

LANES = {
    "runPP": "pp_data",
    "runAuAu": "auau_data",
    "runPhotonJetSim": "pp_photon_sim",
    "runInclusiveJetSim": "pp_inclusive_sim",
    "runEmbeddedPhotonJetSim": "auau_photon_embedded",
    "runEmbeddedInclusiveJetSim": "auau_inclusive_embedded",
}
AGGREGATES = {
    "runAllData": ["runPP", "runAuAu"],
    "runAllSim": list(LANES)[2:],
    "runAll": ["runAllData", "runAllSim"],
}
TRAIN = {
    "trainAuAu": "canonical_auau_v1",
    "trainPP": "canonical_pp_v1",
    "trainPPG12Equivalent": "ppg12_equivalent_v1",
}


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for b in iter(lambda: stream.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read(path):
    return yaml.safe_load(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def safe_name(name):
    if not isinstance(name, str) or not re.fullmatch(
        r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name
    ):
        raise ValueError(f"unsafe identity: {name!r}")
    return name


def expand(command):
    if command in LANES:
        return [command]
    if command not in AGGREGATES:
        raise ValueError(f"unknown lane command: {command}")
    return [item for c in AGGREGATES[command] for item in expand(c)]


def bound_file(entry, base, field="path", hash_field="sha256"):
    if not entry.get(field) or not entry.get(hash_field):
        raise ValueError(f"unresolved {field}/{hash_field} binding")
    p = (Path(base) / entry[field]).resolve()
    if not p.is_file() or file_hash(p) != entry[hash_field]:
        raise ValueError(f"missing file or SHA256 mismatch: {p}")
    return p


def campaign(path, expected_hash):
    """Read an explicit completed-output manifest, not a population/DST list.

    TreeProduction does not emit this contract yet. No directory globbing,
    nearby-file fallback, inferred sample roster, or row-number identity.
    """
    if not path or not expected_hash:
        raise ValueError(
            "TreeProduction completed campaign manifest is unbound; configure input_campaign and input_campaign_sha256"
        )
    path = Path(path).resolve()
    if file_hash(path) != expected_hash:
        raise ValueError("campaign manifest SHA256 mismatch")
    d = read(path)
    if d.get("schema") != "PhotonJetTreeCampaignV1" or d.get("status") != "complete":
        raise ValueError("requires completed PhotonJetTreeCampaignV1 output manifest")
    safe_name(d.get("campaign"))
    if not d.get("sources"):
        raise ValueError("campaign has no sources")
    seen, files = set(), set()
    for s in d["sources"]:
        if s.get("lane") not in LANES.values():
            raise ValueError("unknown campaign lane")
        p = bound_file(s, path.parent)
        if not isinstance(s.get("source_id"), list) or len(s["source_id"]) != 2:
            raise ValueError("source_id must contain exact hi/lo words")
        key = tuple(s["source_id"])
        if any(
            not isinstance(v, int) or isinstance(v, bool) or not 0 <= v < 2**64
            for v in key
        ):
            raise ValueError("source identity is not uint64")
        if key in seen or p in files:
            raise ValueError("duplicate source identity or base file")
        seen.add(key)
        files.add(p)
        s["_path"] = p
    d["_manifest"] = path
    d["_sha256"] = expected_hash
    return d


def load_package(entry, base):
    return read(bound_file(entry, base)) if "path" in entry else dict(entry)


def resolve_selection(raw, name, specs, base):
    from features import feature_identity

    if name not in raw["selection_recipes"]:
        raise ValueError(f"unknown selection: {name}")
    r = dict(raw["selection_recipes"][name])
    spec = specs[r["model"]]
    if not spec.bound:
        raise ValueError(f"{name}: model {spec.name} is unbound")
    wp = load_package(raw["working_points"][r["id_package"]], base)
    if wp.get("model_sha256") != spec.model_sha256 or wp.get(
        "feature_schema_sha256"
    ) != feature_identity(spec):
        raise ValueError(
            f"{name}: ID working points belong to a different model/feature schema"
        )
    iso = load_package(raw["isolation_packages"][r["isolation_package"]], base)
    if wp.get("system") != spec.system or iso.get("system") != spec.system:
        raise ValueError(f"{name}: incompatible collision system")
    if not iso.get("isolated") or not iso.get("nonisolated"):
        raise ValueError(f"{name}: isolated/non-isolated boundary is unbound")
    from calibration import validate_ordering, line
    import numpy as np

    validate_ordering(wp)
    if wp.get("fit_population") != "validation":
        raise ValueError("this recipe requires validation-derived model-specific WPs")
    if (
        iso.get("axis") not in ("et", "centrality")
        or iso.get("method") not in (2, 3)
        or iso.get("radius") != 0.4
    ):
        raise ValueError("unsupported isolation observable/package")
    x = np.array(iso["domain"])
    a = line(iso["isolated"], x)
    b = line(iso["nonisolated"], x)
    if not np.all(np.isfinite([a, b])) or np.any(a > b):
        raise ValueError("isolation boundaries cross or are nonfinite")
    wp_entry = raw["working_points"][r["id_package"]]
    iso_entry = raw["isolation_packages"][r["isolation_package"]]
    return dict(
        name=name,
        model=spec.name,
        wp=wp,
        isolation=iso,
        wp_sha256=wp_entry["sha256"] if "path" in wp_entry else digest(wp),
        isolation_sha256=iso_entry["sha256"] if "path" in iso_entry else digest(iso),
    )
