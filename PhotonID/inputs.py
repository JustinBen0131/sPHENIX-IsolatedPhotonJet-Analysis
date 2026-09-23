"""One explicit adapter from canonical facts to training rows and isolation.

Prompt labels deliberately do not use TruthPhotons.analysis_signal, which
includes truth isolation and has different class semantics. Production bindings
must certify dominant-primary ancestry completeness before using this adapter.
"""

from pathlib import Path
import numpy as np
import uproot
from features import features_for_file, base_metadata
from registry import bound_file, read

KEY = ["event_hi", "event_lo", "photon_hi", "photon_lo"]


def keys(a, names=KEY):
    return list(zip(*(np.asarray(a[n]).tolist() for n in names)))


def isolation_for_file(path, matrix, method, radius=0.4):
    with uproot.open(path) as f:
        a = f["Isolation"].arrays(
            KEY + ["method", "radius", "cone_sum", "valid", "candidate_removed"],
            library="np",
        )
    values = {}
    for i, key in enumerate(keys(a)):
        if int(a["method"][i]) != method or not np.isclose(
            a["radius"][i], radius, atol=1e-8, rtol=0
        ):
            continue
        if key in values:
            raise ValueError("duplicate isolation identity/method/radius")
        values[key] = (
            float(a["cone_sum"][i])
            if a["valid"][i] and a["candidate_removed"][i]
            else np.nan
        )
    return np.array([values.get(key, np.nan) for key in matrix.keys])


def dominant_prompt_labels(path, matrix):
    """Return prompt {-1 unknown,0 known nonprompt,1 direct/fragmentation}.

    Track lookup includes event and embedding identity. A valid non-photon
    dominant primary or the measured NoPrimary state is nonprompt. Unknown
    ancestry/evaluator remains unknown; no Fake-link shortcut invents a label.
    This adapter is implemented but its THE358 equivalence remains unbound.
    """
    with uproot.open(path) as f:
        if np.any(f["Events"]["truth_denominator_complete"].array(library="np") != 1):
            raise ValueError("training requires a complete truth census")
        p = f["Photons"].arrays(
            KEY
            + [
                "dominant_truth_state",
                "dominant_truth_pid",
                "dominant_truth_track_id",
                "dominant_truth_embedding_id",
            ],
            library="np",
        )
        t = f["TruthPhotons"].arrays(
            [
                "event_hi",
                "event_lo",
                "track_id",
                "embedding_id",
                "prompt_class",
                "generator_association_valid",
                "analysis_signal",
            ],
            library="np",
        )
    truth = {}
    for i, key in enumerate(
        keys(t, ["event_hi", "event_lo", "track_id", "embedding_id"])
    ):
        if key in truth:
            raise ValueError("duplicate native truth identity")
        truth[key] = (
            int(t["prompt_class"][i]),
            bool(t["generator_association_valid"][i]),
            bool(t["analysis_signal"][i]),
        )
    rows = {}
    for i, key in enumerate(keys(p)):
        label, analysis = -1, False
        state = int(p["dominant_truth_state"][i])
        if state == 1:
            label = 0
        elif state == 2:
            if int(p["dominant_truth_pid"][i]) != 22:
                label = 0
            else:
                q = truth.get(
                    (
                        *key[:2],
                        int(p["dominant_truth_track_id"][i]),
                        int(p["dominant_truth_embedding_id"][i]),
                    )
                )
                if q is not None:
                    cls, gen_valid, analysis = q
                    if gen_valid and cls in (1, 2, 3):
                        label = int(cls in (1, 2))
        if key in rows:
            raise ValueError("duplicate photon identity")
        rows[key] = (label, analysis)
    return np.array([rows[k][0] for k in matrix.keys]), np.array(
        [rows[k][1] for k in matrix.keys]
    )


def grouped_events(records):
    """Connected components prevent shared signal OR embedded background reuse.

    Manifest witnesses supply physical-event component identities. Event/source
    IDs remain output join keys; they are never mistaken for physical ancestry.
    """
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    for components in records:
        if not components or any(not isinstance(v, str) or not v for v in components):
            raise ValueError("missing physical-event component identity")
        roots = [find(x) for x in components]
        low = min(roots)
        for r in roots:
            parent[r] = low
    return np.asarray([find(c[0]) for c in records])


def load_training(campaign, spec, profile):
    if (
        not profile.get("label_binding_evidence")
        or profile.get("label_adapter") != "dominant_prompt_v1"
    ):
        raise ValueError(
            "dominant_prompt_v1 truth-label mapping is unbound; do not substitute analysis_signal"
        )
    if not profile.get("preselection_binding_evidence"):
        raise ValueError(
            "training population/event/preselection equivalence is unbound"
        )
    required = (
        ["auau_photon_embedded", "auau_inclusive_embedded"]
        if spec.system == "auau"
        else ["pp_photon_sim", "pp_inclusive_sim"]
    )
    rows = [s for s in campaign["sources"] if s["lane"] in required]
    if set(s["lane"] for s in rows) != set(required):
        raise ValueError("campaign lacks a required simulation lane")
    if set(s.get("sample") for s in rows) != set(profile["source_order"]):
        raise ValueError(
            "canonical simulation source/slice roster is incomplete or unexpected"
        )
    rows = sorted(rows, key=lambda s: profile["source_order"].index(s["sample"]))
    pieces, component_rows, all_keys = [], [], []
    for s in rows:
        path = s["_path"]
        m = features_for_file(path, spec)
        with uproot.open(path) as f:
            meta = base_metadata(f)
            source = f["Sources"].arrays(["source_hi", "source_lo"], library="np")
            if [int(source["source_hi"][0]), int(source["source_lo"][0])] != s[
                "source_id"
            ]:
                raise ValueError("campaign source identity differs from ROOT")
            expected_role = "1" if "photon" in s["lane"] else "2"
            if (
                meta.get("data_kind") != "2"
                or meta.get("simulation_role") != expected_role
            ):
                raise ValueError("source role differs from ROOT metadata")
        prompt, analysis = dominant_prompt_labels(path, m)
        role = "signal" if "photon" in s["lane"] else "background"
        y = np.where(
            prompt == (1 if role == "signal" else 0), 1 if role == "signal" else 0, -1
        )
        # Population witnesses explicitly bind reviewed event/preselection policy.
        # They are temporary adapter inputs until those predicates are frozen.
        witness = read(bound_file(s["training_witness"], campaign["_manifest"].parent))
        if (
            witness.get("base_sha256") != s["sha256"]
            or witness.get("label_adapter") != profile["label_adapter"]
        ):
            raise ValueError("training witness belongs to another base/adapter")
        events = {tuple(v["event_id"]): v for v in witness.get("events", [])}
        if len(events) != len(witness.get("events", [])):
            raise ValueError("duplicate event witness")
        policy = {tuple(v["photon_id"]): v for v in witness["photons"]}
        if len(policy) != len(witness["photons"]):
            raise ValueError("duplicate candidate witness")
        eligible, iso_eligible, comp = [], [], []
        for key in m.keys:
            p = policy[key]
            if profile["split"]["rule"] == "physical_event_hash_v1":
                e = events[key[:2]]
                if (
                    spec.system == "auau"
                    and e.get("embedding_components_complete") is not True
                ):
                    raise ValueError("embedded physical-event grouping is incomplete")
                comp.append(e["physical_components"])
            else:
                comp.append(
                    [
                        "row-parity-event:"
                        + ":".join(str(v) for v in (*s["source_id"], *key[:2]))
                    ]
                )
            if any(
                type(p.get(field)) is not bool
                for field in ("training_eligible", "isolation_eligible")
            ):
                raise ValueError(
                    "training/isolation eligibility witnesses must be explicit booleans"
                )
            eligible.append(p["training_eligible"])
            iso_eligible.append(p["isolation_eligible"])
        if not np.all(m.source == np.asarray(s["source_id"], dtype=np.uint64)):
            raise ValueError("candidate source differs from manifest")
        component_rows.extend(comp)
        fullkeys = [(*src, *key) for src, key in zip(m.source.tolist(), m.keys)]
        all_keys.extend(fullkeys)
        pieces.append(
            dict(
                X=m.values,
                y=y,
                et=m.et,
                eta=m.eta,
                centrality=m.centrality,
                ready=m.complete & m.in_domain & np.asarray(eligible),
                analysis=analysis
                & bool(profile.get("analysis_signal_binding_evidence"))
                & np.asarray(iso_eligible)
                & m.in_domain
                & (role == "signal"),
                isolation=isolation_for_file(
                    path, m, 2 if spec.system == "auau" else 3
                ),
                sample=np.full(len(y), s["sample"]),
                keys=np.asarray(fullkeys, dtype=np.uint64),
            )
        )
    if len(set(all_keys)) != len(all_keys):
        raise ValueError("duplicate source/event/photon training identity")
    out = {k: np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    out["groups"] = grouped_events(component_rows)
    return out
