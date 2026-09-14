"""Strict ``PhotonJetTrees_v1`` to response-bundle construction.

This module is deliberately an engineering-fixture implementation boundary.
It consumes only the public normalized trees, validates the complete tree
contract before reading physics rows, and builds the same ``ResponseBundle``
shape for p+p and Au+Au.  It does not claim a mini-DST golden differential,
full-statistics response validation, or a production-ready model lineage.

The construction unit is one event-leading photon plus every accepted recoil
jet.  One authoritative outgoing typed match is allowed per reconstructed
object; deterministic jet matches are also one-to-one, while source-preserved
photon matches may be many-reco-to-one-truth before event-leading selection.
Link absence is interpreted as a fake or miss; conflicting links are errors.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
import fcntl
import hashlib
from importlib import resources
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np
import uproot

from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.response import (
    CLASSIFICATION_RECO_PTGAMMA_EDGES,
    CLASSIFICATION_TRUTH_PTGAMMA_EDGES,
    COMMON_XJ_EDGES,
    Category,
    GlobalBinMap,
    ResponseBundle,
    classify_state,
)
from photonjet.analysis.selection import SelectionProgram, compile_recoil_selection
from photonjet.io.tree_validation import validate as validate_trees
from photonjet.provenance import artifact, write_json


Identity = tuple[int, int]
EventKey = tuple[int, int, int]

_PHOTON = 1
_JET = 2
_MATCH_LINK_CLASS = 0
_RECO_FAKE_LINK_CLASS = 1
_TRUTH_MISS_LINK_CLASS = 2
_WRONG_PHOTON_LINK_CLASS = 3
_WRONG_RECOIL_LINK_CLASS = 4
_MATCH_CANDIDATE_LINK_CLASS = 5
_STATUS = "TREE_CONSTRUCTION_COMPLETE__GOLDEN_EQUIVALENCE_NOT_RUN"


@dataclass(frozen=True)
class ResponseBuildConfig:
    """Closed configuration for the first public tree-to-response stage.

    Truth signal ownership is intentionally explicit.  The defaults match the
    redistribution-safe fixture contract; changing them changes the receipted
    configuration but does not create a new hidden mode.
    """

    system: str
    dimension: str = "2D"
    selection: RecoilSelection = field(default_factory=RecoilSelection)
    truth_prompt_class: int = 1
    truth_source_role: int = 1
    truth_isolation_max: float = 4.0

    def __post_init__(self) -> None:
        if self.system not in {"pp", "auau"}:
            raise ValueError("response system must be pp or auau")
        if self.dimension not in {"1D", "2D"}:
            raise ValueError("response dimension must be 1D or 2D")
        if self.selection.region == "inclusive":
            raise ValueError("response construction requires an event-leading ABCD region")
        if not isinstance(self.truth_prompt_class, int):
            raise ValueError("truth_prompt_class must be an integer")
        if not isinstance(self.truth_source_role, int):
            raise ValueError("truth_source_role must be an integer")
        if not math.isfinite(self.truth_isolation_max) or self.truth_isolation_max <= 0:
            raise ValueError("truth_isolation_max must be finite and positive")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "PhotonJetResponseBuildConfigV1",
            "system": self.system,
            "dimension": self.dimension,
            "selection": asdict(self.selection),
            "truth_selection": {
                "prompt_class": self.truth_prompt_class,
                "source_role": self.truth_source_role,
                "isolation": {
                    "operator": "lt",
                    "value": self.truth_isolation_max,
                    "unit": "GeV",
                },
                "photon_abs_eta_max": self.selection.photon_abs_eta_max,
                "jet_pt_min": self.selection.jet_pt_min,
                "jet_abs_eta_max": self.selection.jet_abs_eta_max,
                "jet_radius": self.selection.jet_radius,
                "delta_phi_min": self.selection.delta_phi_min,
            },
        }


@dataclass(frozen=True)
class _Pair:
    photon_id: Identity
    jet_id: Identity
    ptgamma: float
    xjgamma: float


@dataclass
class _Accumulators:
    matrix: np.ndarray
    matrix_sumw2: np.ndarray
    truth: np.ndarray
    truth_sumw2: np.ndarray
    reco: np.ndarray
    reco_sumw2: np.ndarray
    misses: np.ndarray
    misses_sumw2: np.ndarray
    fakes: np.ndarray
    fakes_sumw2: np.ndarray
    fake_causes: dict[str, np.ndarray]
    fake_causes_sumw2: dict[str, np.ndarray]
    boundary_reco: dict[str, np.ndarray]
    boundary_reco_sumw2: dict[str, np.ndarray]
    boundary_truth: dict[str, np.ndarray]
    boundary_truth_sumw2: dict[str, np.ndarray]


def _identity(row: Mapping[str, Any], prefix: str) -> Identity:
    return int(row[f"{prefix}_hi"]), int(row[f"{prefix}_lo"])


def _event_key(row: Mapping[str, Any]) -> EventKey:
    return (
        int(row["source_file_index"]),
        int(row["event_id_hi"]),
        int(row["event_id_lo"]),
    )


def _python_scalar(value: Any) -> Any:
    return value.item() if isinstance(value, np.generic) else value


def _rows(tree: uproot.behaviors.TTree.TTree, branches: Sequence[str]) -> list[dict[str, Any]]:
    missing = sorted(set(branches) - set(tree.keys()))
    if missing:
        raise ValueError(f"{tree.name} lacks response branches: {', '.join(missing)}")
    arrays = tree.arrays(list(branches), library="np")
    return [
        {name: _python_scalar(arrays[name][index]) for name in branches}
        for index in range(tree.num_entries)
    ]


def _group(rows: Sequence[dict[str, Any]]) -> dict[EventKey, list[dict[str, Any]]]:
    grouped: dict[EventKey, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[_event_key(row)].append(row)
    return dict(grouped)


def _finite(value: Any, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def _delta_phi(left: float, right: float) -> float:
    return abs(math.atan2(math.sin(left - right), math.cos(left - right)))


def _axis_bin(value: float, edges: np.ndarray) -> int | None:
    if not math.isfinite(value):
        raise ValueError("response coordinate is not finite")
    index = int(np.searchsorted(edges, value, side="right") - 1)
    return index if 0 <= index < len(edges) - 1 else None


def _pair_bin(pair: _Pair, config: ResponseBuildConfig, *, truth: bool) -> int | None:
    pt_edges = (
        CLASSIFICATION_TRUTH_PTGAMMA_EDGES
        if truth
        else CLASSIFICATION_RECO_PTGAMMA_EDGES
    )
    pt_bin = _axis_bin(pair.ptgamma, pt_edges)
    if pt_bin is None:
        return None
    if config.dimension == "1D":
        return pt_bin
    xj_bin = _axis_bin(pair.xjgamma, COMMON_XJ_EDGES)
    if xj_bin is None:
        return None
    return GlobalBinMap(pt_edges, COMMON_XJ_EDGES).flatten(pt_bin, xj_bin)


def _empty_accumulators(config: ResponseBuildConfig) -> _Accumulators:
    truth_size = len(CLASSIFICATION_TRUTH_PTGAMMA_EDGES) - 1
    reco_size = len(CLASSIFICATION_RECO_PTGAMMA_EDGES) - 1
    if config.dimension == "2D":
        xj_size = len(COMMON_XJ_EDGES) - 1
        truth_size *= xj_size
        reco_size *= xj_size

    def zeros(size: int) -> np.ndarray:
        return np.zeros(size, dtype=np.float64)

    return _Accumulators(
        matrix=np.zeros((truth_size, reco_size), dtype=np.float64),
        matrix_sumw2=np.zeros((truth_size, reco_size), dtype=np.float64),
        truth=zeros(truth_size),
        truth_sumw2=zeros(truth_size),
        reco=zeros(reco_size),
        reco_sumw2=zeros(reco_size),
        misses=zeros(truth_size),
        misses_sumw2=zeros(truth_size),
        fakes=zeros(reco_size),
        fakes_sumw2=zeros(reco_size),
        fake_causes={},
        fake_causes_sumw2={},
        boundary_reco={},
        boundary_reco_sumw2={},
        boundary_truth={},
        boundary_truth_sumw2={},
    )


def _component(
    values: dict[str, np.ndarray],
    key: str,
    shape_like: np.ndarray,
) -> np.ndarray:
    if key not in values:
        values[key] = np.zeros_like(shape_like)
    return values[key]


def _fill(values: np.ndarray, sumw2: np.ndarray, index: int, weight: float) -> None:
    values[index] += weight
    sumw2[index] += weight * weight


def _truth_photon(
    rows: Sequence[dict[str, Any]],
    config: ResponseBuildConfig,
    event: EventKey,
) -> dict[str, Any] | None:
    eligible = []
    for row in rows:
        pt = _finite(row["truth_photon_pt"], "truth photon pT")
        eta = _finite(row["truth_photon_eta"], "truth photon eta")
        phi = _finite(row["truth_photon_phi"], "truth photon phi")
        isolation = _finite(row["truth_isolation"], "truth photon isolation")
        if pt <= 0:
            raise ValueError(f"truth photon pT is not positive in event {event}")
        if (
            int(row["prompt_class"]) == config.truth_prompt_class
            and int(row["source_role"]) == config.truth_source_role
            and isolation < config.truth_isolation_max
            and abs(eta) < config.selection.photon_abs_eta_max
        ):
            selected = dict(row)
            selected.update(
                truth_photon_pt=pt,
                truth_photon_eta=eta,
                truth_photon_phi=phi,
                truth_isolation=isolation,
            )
            eligible.append(selected)
    if len(eligible) > 1:
        identities = sorted(_identity(row, "truth_photon_id") for row in eligible)
        raise ValueError(f"event {event} has ambiguous truth-signal photons: {identities}")
    return eligible[0] if eligible else None


def _truth_pairs(
    photon: dict[str, Any] | None,
    jets: Sequence[dict[str, Any]],
    config: ResponseBuildConfig,
    event: EventKey,
) -> list[_Pair]:
    if photon is None:
        return []
    photon_id = _identity(photon, "truth_photon_id")
    photon_pt = float(photon["truth_photon_pt"])
    photon_phi = float(photon["truth_photon_phi"])
    result: list[_Pair] = []
    seen: set[Identity] = set()
    for row in jets:
        jet_id = _identity(row, "truth_jet_id")
        if jet_id in seen:
            raise ValueError(f"duplicate truth-jet identity in event {event}: {jet_id}")
        seen.add(jet_id)
        pt = _finite(row["truth_jet_pt"], "truth jet pT")
        eta = _finite(row["truth_jet_eta"], "truth jet eta")
        phi = _finite(row["truth_jet_phi"], "truth jet phi")
        radius = _finite(row["truth_jet_radius"], "truth jet radius")
        if pt <= 0:
            raise ValueError(f"truth jet pT is not positive in event {event}")
        if not (
            pt > config.selection.jet_pt_min
            and abs(eta) < config.selection.jet_abs_eta_max
            and abs(radius - config.selection.jet_radius) < 1.0e-9
            and _delta_phi(phi, photon_phi) > config.selection.delta_phi_min
        ):
            continue
        result.append(
            _Pair(
                photon_id=photon_id,
                jet_id=jet_id,
                ptgamma=photon_pt,
                xjgamma=pt / photon_pt,
            )
        )
    return sorted(result, key=lambda pair: pair.jet_id)


def _candidate_record(row: Mapping[str, Any]) -> dict[str, float]:
    return {
        "photon_et": float(row["photon_et"]),
        "photon_eta": float(row["photon_eta"]),
        "photon_encounter_ordinal": float(row["photon_encounter_ordinal"]),
        "photon_bdt_score": float(row["bdt_score"]),
        "photon_bdt_tight_threshold": float(row["bdt_tight_threshold"]),
        "photon_bdt_nontight_low_threshold": float(row["bdt_nontight_low_threshold"]),
        "photon_bdt_nontight_high_threshold": float(row["bdt_nontight_high_threshold"]),
        "photon_iso_r04": float(row["iso_r04"]),
        "photon_iso_r04_threshold": float(row["iso_r04_threshold"]),
        "photon_iso_r04_nonisolated_threshold": float(
            row["iso_r04_nonisolated_threshold"]
        ),
    }


def _reco_pairs(
    photons: Sequence[dict[str, Any]],
    jets: Sequence[dict[str, Any]],
    pairs: Sequence[dict[str, Any]],
    event_tree: Mapping[str, Any],
    program: SelectionProgram,
    event: EventKey,
) -> list[_Pair]:
    candidates = [_candidate_record(row) for row in photons]
    leader = program.choose_leader(candidates)
    if program.leader_branch is not None:
        recorded = int(event_tree[program.leader_branch])
        if recorded != leader:
            raise ValueError(
                f"event {event} recorded {program.leader_branch}={recorded} "
                f"but response selection computed {leader}"
            )
    if leader >= len(photons):
        raise ValueError(f"event {event} response leader is outside photons")
    leader_id = (
        _identity(photons[leader], "candidate_id")
        if leader >= 0
        else None
    )
    observed_pair_indices = [
        (int(row["photon_index"]), int(row["jet_index"]))
        for row in pairs
    ]
    expected_pair_indices = {
        (photon_index, jet_index)
        for photon_index in range(len(photons))
        for jet_index in range(len(jets))
    }
    if (
        len(observed_pair_indices) != len(set(observed_pair_indices))
        or set(observed_pair_indices) != expected_pair_indices
    ):
        raise ValueError(
            f"event {event} photonJets does not exactly cover the photon/jet "
            "Cartesian product"
        )
    result: list[_Pair] = []
    seen_jets: set[Identity] = set()
    for row in pairs:
        photon_index = int(row["photon_index"])
        jet_index = int(row["jet_index"])
        if not 0 <= photon_index < len(photons):
            raise ValueError(f"event {event} pair photon index is outside photons")
        if not 0 <= jet_index < len(jets):
            raise ValueError(f"event {event} pair jet index is outside jets")
        if _identity(row, "candidate_id") != _identity(
            photons[photon_index], "candidate_id"
        ):
            raise ValueError(f"event {event} pair photon index and identity disagree")
        if _identity(row, "jet_id") != _identity(jets[jet_index], "jet_id"):
            raise ValueError(f"event {event} pair jet index and identity disagree")

        photon = photons[photon_index]
        jet = jets[jet_index]
        photon_et = _finite(photon["photon_et"], "reconstructed photon ET")
        photon_eta = _finite(photon["photon_eta"], "reconstructed photon eta")
        photon_phi = _finite(photon["photon_phi"], "reconstructed photon phi")
        jet_pt = _finite(jet["jet_pt"], "reconstructed jet pT")
        jet_eta = _finite(jet["jet_eta"], "reconstructed jet eta")
        jet_phi = _finite(jet["jet_phi"], "reconstructed jet phi")
        jet_radius = _finite(jet["jet_radius"], "reconstructed jet radius")
        pair_delta_phi = _finite(row["delta_phi"], "photonJets delta_phi")
        if photon_et <= 0:
            raise ValueError(f"event {event} has a nonpositive reconstructed photon ET")
        expected_xjgamma = jet_pt / photon_et
        copied_witnesses = {
            "photon_et": photon_et,
            "photon_eta": photon_eta,
            "photon_phi": photon_phi,
            "jet_pt": jet_pt,
            "jet_eta": jet_eta,
            "jet_phi": jet_phi,
            "jet_radius": jet_radius,
            "xjgamma": expected_xjgamma,
        }
        for name, expected in copied_witnesses.items():
            if _finite(row[name], f"photonJets {name}") != expected:
                raise ValueError(
                    f"event {event} photonJets {name} disagrees with normalized objects"
                )
        if int(row["wrong_photon_class"]) != 0 or int(row["wrong_recoil_class"]) != 0:
            raise ValueError(
                "nonzero producer wrong-object witnesses require the golden response "
                "semantics and are not accepted by the engineering builder"
            )
        if leader < 0:
            continue
        if photon_index != leader:
            continue
        if _identity(row, "candidate_id") != leader_id:
            raise ValueError(f"event {event} pair photon index and identity disagree")
        if not program.accepts(
            "recoil",
            {
                "jet_pt": jet_pt,
                "jet_eta": jet_eta,
                "jet_radius": jet_radius,
                "delta_phi": pair_delta_phi,
            },
        ):
            continue
        jet_id = _identity(row, "jet_id")
        if jet_id in seen_jets:
            raise ValueError(f"event {event} has duplicate selected photon-jet rows")
        seen_jets.add(jet_id)
        if expected_xjgamma < 0:
            raise ValueError(f"event {event} has a nonphysical selected recoil pair")
        result.append(
            _Pair(
                photon_id=leader_id,
                jet_id=jet_id,
                ptgamma=photon_et,
                xjgamma=expected_xjgamma,
            )
        )
    return sorted(result, key=lambda pair: pair.jet_id)


def _typed_links(
    rows: Sequence[dict[str, Any]],
    event: EventKey,
) -> tuple[dict[Identity, Identity], dict[Identity, Identity]]:
    forward: dict[int, dict[Identity, Identity]] = {_PHOTON: {}, _JET: {}}
    matched_truth: dict[int, set[Identity]] = {_PHOTON: set(), _JET: set()}
    explicit_fakes: dict[int, set[Identity]] = {_PHOTON: set(), _JET: set()}
    explicit_misses: dict[int, set[Identity]] = {_PHOTON: set(), _JET: set()}
    for row in rows:
        reco_type = int(row["reco_type"])
        truth_type = int(row["truth_type"])
        link_class = int(row["link_class"])
        if link_class == _RECO_FAKE_LINK_CLASS:
            if reco_type not in {_PHOTON, _JET} or truth_type != 0:
                raise ValueError(f"event {event} has invalid explicit reco-fake topology")
            explicit_fakes[reco_type].add(_identity(row, "reco_id"))
            continue
        if link_class == _TRUTH_MISS_LINK_CLASS:
            if reco_type != 0 or truth_type not in {_PHOTON, _JET}:
                raise ValueError(f"event {event} has invalid explicit truth-miss topology")
            explicit_misses[truth_type].add(_identity(row, "truth_id"))
            continue
        if link_class == _MATCH_CANDIDATE_LINK_CLASS:
            if reco_type not in {_PHOTON, _JET} or truth_type not in {_PHOTON, _JET}:
                raise ValueError(f"event {event} has invalid match-candidate topology")
            continue
        if link_class in {_WRONG_PHOTON_LINK_CLASS, _WRONG_RECOIL_LINK_CLASS}:
            raise ValueError(
                "nonzero producer wrong-object link witnesses require the golden "
                "response semantics and are not accepted by the engineering builder"
            )
        if link_class != _MATCH_LINK_CLASS:
            raise ValueError(f"event {event} has unsupported response link_class")
        if reco_type not in {_PHOTON, _JET} or truth_type != reco_type:
            raise ValueError(f"event {event} has unsupported response match topology")
        reco_id = _identity(row, "reco_id")
        truth_id = _identity(row, "truth_id")
        if reco_id in forward[reco_type]:
            raise ValueError(f"event {event} has an ambiguous typed reco link")
        if reco_type == _JET and truth_id in matched_truth[_JET]:
            raise ValueError(f"event {event} has a many-to-one typed truth link")
        forward[reco_type][reco_id] = truth_id
        matched_truth[reco_type].add(truth_id)
    for target_type in (_PHOTON, _JET):
        if set(forward[target_type]) & explicit_fakes[target_type]:
            raise ValueError(f"event {event} marks a matched reco object as an explicit fake")
        if matched_truth[target_type] & explicit_misses[target_type]:
            raise ValueError(f"event {event} marks a matched truth object as an explicit miss")
    return forward[_PHOTON], forward[_JET]


def _boundary_key(
    truth_pair: _Pair,
    reco_pair: _Pair,
    config: ResponseBuildConfig,
) -> str:
    categories = classify_state(
        truth_ptgamma=truth_pair.ptgamma,
        reco_ptgamma=reco_pair.ptgamma,
        truth_xj=(truth_pair.xjgamma if config.dimension == "2D" else None),
        reco_xj=(reco_pair.xjgamma if config.dimension == "2D" else None),
    )
    return "+".join(category.value for category in categories)


def _accumulate_event(
    accumulators: _Accumulators,
    *,
    truth_pairs: Sequence[_Pair],
    reco_pairs: Sequence[_Pair],
    photon_links: Mapping[Identity, Identity],
    jet_links: Mapping[Identity, Identity],
    config: ResponseBuildConfig,
    weight: float,
    category_counts: Counter[str],
    category_sumw: Counter[str],
) -> tuple[int, int, int]:
    truth_by_id = {
        (pair.photon_id, pair.jet_id): pair
        for pair in truth_pairs
    }
    if len(truth_by_id) != len(truth_pairs):
        raise ValueError("selected truth response pairs are not unique")

    matched: dict[tuple[Identity, Identity], _Pair] = {}
    matched_reco: set[_Pair] = set()
    boundary_reco_count = 0
    fake_count = 0
    for reco_pair in reco_pairs:
        reco_bin = _pair_bin(reco_pair, config, truth=False)
        if reco_bin is not None:
            _fill(accumulators.reco, accumulators.reco_sumw2, reco_bin, weight)

        linked_photon = photon_links.get(reco_pair.photon_id)
        linked_jet = jet_links.get(reco_pair.jet_id)
        target = (
            (linked_photon, linked_jet)
            if linked_photon is not None and linked_jet is not None
            else None
        )
        truth_pair = truth_by_id.get(target) if target is not None else None
        if truth_pair is not None:
            if target in matched:
                raise ValueError("multiple selected reco pairs resolve to one truth pair")
            matched[target] = reco_pair
            matched_reco.add(reco_pair)
            truth_bin = _pair_bin(truth_pair, config, truth=True)
            categories = classify_state(
                truth_ptgamma=truth_pair.ptgamma,
                reco_ptgamma=reco_pair.ptgamma,
                truth_xj=(truth_pair.xjgamma if config.dimension == "2D" else None),
                reco_xj=(reco_pair.xjgamma if config.dimension == "2D" else None),
            )
            for category in categories:
                category_counts[category.value] += 1
                category_sumw[category.value] += weight
            if truth_bin is not None and reco_bin is not None:
                accumulators.matrix[truth_bin, reco_bin] += weight
                accumulators.matrix_sumw2[truth_bin, reco_bin] += weight * weight
            elif reco_bin is not None:
                key = _boundary_key(truth_pair, reco_pair, config)
                _fill(
                    _component(accumulators.boundary_reco, key, accumulators.reco),
                    _component(
                        accumulators.boundary_reco_sumw2,
                        key,
                        accumulators.reco_sumw2,
                    ),
                    reco_bin,
                    weight,
                )
                boundary_reco_count += 1
            elif truth_bin is not None:
                key = _boundary_key(truth_pair, reco_pair, config)
                _fill(
                    _component(accumulators.boundary_truth, key, accumulators.truth),
                    _component(
                        accumulators.boundary_truth_sumw2,
                        key,
                        accumulators.truth_sumw2,
                    ),
                    truth_bin,
                    weight,
                )
            continue

        if reco_bin is None:
            continue
        if linked_photon is not None and linked_jet is not None:
            boundary_key = Category.ASSIGNED_HARD_NONFIDUCIAL.value
            _fill(
                _component(accumulators.boundary_reco, boundary_key, accumulators.reco),
                _component(
                    accumulators.boundary_reco_sumw2,
                    boundary_key,
                    accumulators.reco_sumw2,
                ),
                reco_bin,
                weight,
            )
            category_counts[boundary_key] += 1
            category_sumw[boundary_key] += weight
            boundary_reco_count += 1
            continue
        cause = (
            Category.UNMATCHED_RECO.value
            if linked_photon is None
            else Category.COMBINATORIC.value
        )
        _fill(accumulators.fakes, accumulators.fakes_sumw2, reco_bin, weight)
        _fill(
            _component(accumulators.fake_causes, cause, accumulators.fakes),
            _component(
                accumulators.fake_causes_sumw2,
                cause,
                accumulators.fakes_sumw2,
            ),
            reco_bin,
            weight,
        )
        category_counts[cause] += 1
        category_sumw[cause] += weight
        fake_count += 1

    miss_count = 0
    for truth_pair in truth_pairs:
        truth_bin = _pair_bin(truth_pair, config, truth=True)
        if truth_bin is not None:
            _fill(accumulators.truth, accumulators.truth_sumw2, truth_bin, weight)
        target = (truth_pair.photon_id, truth_pair.jet_id)
        if target in matched:
            continue
        if truth_bin is not None:
            _fill(accumulators.misses, accumulators.misses_sumw2, truth_bin, weight)
            category_counts[Category.PHOTON_RECO_MISS.value] += 1
            category_sumw[Category.PHOTON_RECO_MISS.value] += weight
            miss_count += 1
    return len(matched_reco), fake_count + boundary_reco_count, miss_count


_EVENT_BRANCHES = (
    "source_file_index",
    "event_id_hi",
    "event_id_lo",
    "event_weight",
    "centrality",
    "terminal_status",
)
_PHOTON_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "candidate_id_hi", "candidate_id_lo", "photon_encounter_ordinal",
    "photon_et", "photon_eta", "photon_phi", "bdt_score", "bdt_tight_threshold",
    "bdt_nontight_low_threshold", "bdt_nontight_high_threshold",
    "iso_r04", "iso_r04_threshold", "iso_r04_nonisolated_threshold",
)
_JET_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "jet_id_hi", "jet_id_lo", "jet_pt", "jet_eta", "jet_phi", "jet_radius",
)
_PAIR_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "candidate_id_hi", "candidate_id_lo", "jet_id_hi", "jet_id_lo",
    "photon_index", "jet_index", "jet_pt", "jet_eta", "jet_phi", "jet_radius",
    "delta_phi", "photon_et", "photon_eta", "photon_phi", "xjgamma",
    "wrong_photon_class", "wrong_recoil_class",
)
_TRUTH_PHOTON_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "truth_photon_id_hi", "truth_photon_id_lo", "truth_photon_pt",
    "truth_photon_eta", "truth_photon_phi", "prompt_class", "source_role",
    "truth_isolation",
)
_TRUTH_JET_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "truth_jet_id_hi", "truth_jet_id_lo", "truth_jet_radius",
    "truth_jet_pt", "truth_jet_eta", "truth_jet_phi",
)
_LINK_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "reco_type", "reco_id_hi", "reco_id_lo", "truth_type",
    "truth_id_hi", "truth_id_lo", "link_class",
)


def _canonical_inputs(input_paths: Iterable[Path]) -> tuple[list[Path], list[dict[str, Any]]]:
    paths = [Path(path).resolve() for path in input_paths]
    if not paths:
        raise ValueError("at least one PhotonJetTrees_v1 input is required")
    if len(paths) != len(set(paths)):
        raise ValueError("duplicate response input path")
    for path in paths:
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"response input is not a regular file: {path}")
    entries = [(path, artifact(path, "PhotonJetTrees_v1")) for path in paths]
    names = [entry[1]["name"] for entry in entries]
    hashes = [entry[1]["sha256"] for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError("response inputs have ambiguous duplicate basenames")
    if len(hashes) != len(set(hashes)):
        raise ValueError("response inputs repeat identical content")
    entries.sort(key=lambda entry: (entry[1]["sha256"], entry[1]["name"]))
    return [entry[0] for entry in entries], [entry[1] for entry in entries]


def _tree_contract_sha256() -> str:
    try:
        payload = (
            resources.files("photonjet_contracts")
            .joinpath("photonjet_trees_v1_branches.json")
            .read_bytes()
        )
    except ModuleNotFoundError:
        payload = (
            Path(__file__).resolve().parents[3]
            / "contracts"
            / "photonjet_trees_v1_branches.json"
        ).read_bytes()
    return hashlib.sha256(payload).hexdigest()


def _precheck_collision_system(path: Path, system: str) -> None:
    """Reject an obvious system mismatch before selecting a model dimension."""

    with uproot.open(path) as root:
        if "events" not in root or "centrality" not in root["events"]:
            return
        centrality = np.asarray(
            root["events"]["centrality"].array(library="np"),
            dtype=float,
        )
    if system == "pp" and not np.all(centrality == -1.0):
        raise ValueError("p+p response input does not use centrality=-1")
    if system == "auau" and not (
        np.all(np.isfinite(centrality))
        and np.all(centrality >= 0.0)
        and np.all(centrality < 100.0)
    ):
        raise ValueError("Au+Au response input centrality is outside [0, 100)")


def build_response(
    input_paths: Iterable[Path],
    config: ResponseBuildConfig,
) -> ResponseBundle:
    """Build one deterministic response bundle from validated public trees."""

    paths, input_artifacts = _canonical_inputs(input_paths)
    model_input_count = 11 if config.system == "pp" else 14
    program = compile_recoil_selection(config.selection)
    accumulators = _empty_accumulators(config)
    category_counts: Counter[str] = Counter()
    category_sumw: Counter[str] = Counter()
    seen_events: set[EventKey] = set()
    event_payloads: list[
        tuple[
            EventKey,
            float,
            list[_Pair],
            list[_Pair],
            dict[Identity, Identity],
            dict[Identity, Identity],
        ]
    ] = []

    for path in paths:
        _precheck_collision_system(path, config.system)
    # Validate the collection as one logical input. In addition to validating
    # every file, this rejects event and object identities repeated across
    # parts, which per-file validation alone cannot see.
    validate_trees(paths, model_input_count=model_input_count)

    for path in paths:
        with uproot.open(path) as root:
            events = _rows(root["events"], _EVENT_BRANCHES)
            event_tree_branches = ["source_file_index", "event_id_hi", "event_id_lo"]
            if program.leader_branch is not None:
                event_tree_branches.append(program.leader_branch)
            event_trees = {
                _event_key(row): row
                for row in _rows(root["eventTree"], event_tree_branches)
            }
            photon_groups = _group(_rows(root["photons"], _PHOTON_BRANCHES))
            jet_groups = _group(_rows(root["jets"], _JET_BRANCHES))
            pair_groups = _group(_rows(root["photonJets"], _PAIR_BRANCHES))
            truth_photon_groups = _group(
                _rows(root["truthPhotons"], _TRUTH_PHOTON_BRANCHES)
            )
            truth_jet_groups = _group(_rows(root["truthJets"], _TRUTH_JET_BRANCHES))
            link_groups = _group(_rows(root["recoTruthLinks"], _LINK_BRANCHES))

        for event_row in events:
            event = _event_key(event_row)
            if event in seen_events:
                raise ValueError(f"event/source identity is repeated across inputs: {event}")
            seen_events.add(event)
            if event not in event_trees:
                raise ValueError(f"event {event} has no eventTree row")
            if not program.accepts_event(event_row["terminal_status"]):
                continue
            centrality = _finite(event_row["centrality"], "event centrality")
            if config.system == "pp" and centrality != -1.0:
                raise ValueError("p+p response input does not use centrality=-1")
            if config.system == "auau" and not 0.0 <= centrality < 100.0:
                raise ValueError("Au+Au response input centrality is outside [0, 100)")
            weight = _finite(event_row["event_weight"], "event weight")
            if weight < 0:
                raise ValueError("negative response weights are not supported by ResponseBundleV1")
            truth_photon = _truth_photon(
                truth_photon_groups.get(event, []), config, event
            )
            truth_pairs = _truth_pairs(
                truth_photon,
                truth_jet_groups.get(event, []),
                config,
                event,
            )
            reco_pairs = _reco_pairs(
                photon_groups.get(event, []),
                jet_groups.get(event, []),
                pair_groups.get(event, []),
                event_trees[event],
                program,
                event,
            )
            photon_links, jet_links = _typed_links(link_groups.get(event, []), event)
            event_payloads.append(
                (event, weight, truth_pairs, reco_pairs, photon_links, jet_links)
            )

        known = {_event_key(row) for row in events}
        for label, grouped in (
            ("eventTree", event_trees),
            ("photons", photon_groups),
            ("jets", jet_groups),
            ("photonJets", pair_groups),
            ("truthPhotons", truth_photon_groups),
            ("truthJets", truth_jet_groups),
            ("recoTruthLinks", link_groups),
        ):
            foreign = sorted(set(grouped) - known)
            if foreign:
                raise ValueError(f"{label} has rows outside the event/source set: {foreign[0]}")

    consumed_input_artifacts = [artifact(path, "PhotonJetTrees_v1") for path in paths]
    if consumed_input_artifacts != input_artifacts:
        raise ValueError("response input changed while it was being validated or read")

    matched_total = 0
    fake_total = 0
    miss_total = 0
    truth_pair_total = 0
    reco_pair_total = 0
    for _, weight, truth_pairs, reco_pairs, photon_links, jet_links in sorted(
        event_payloads,
        key=lambda payload: payload[0],
    ):
        truth_pair_total += len(truth_pairs)
        reco_pair_total += len(reco_pairs)
        matched, fakes, misses = _accumulate_event(
            accumulators,
            truth_pairs=truth_pairs,
            reco_pairs=reco_pairs,
            photon_links=photon_links,
            jet_links=jet_links,
            config=config,
            weight=weight,
            category_counts=category_counts,
            category_sumw=category_sumw,
        )
        matched_total += matched
        fake_total += fakes
        miss_total += misses

    provenance = {
        "schema": "PhotonJetResponseBuildProvenanceV1",
        "state": _STATUS,
        "input_ordering": "sha256_then_name",
        "inputs": input_artifacts,
        "tree_contract": "PhotonJetTrees_v1",
        "tree_contract_sha256": _tree_contract_sha256(),
        "configuration": config.to_dict(),
        "selection_program": program.to_dict(),
        "selection_program_sha256": program.sha256,
        "observed": {
            "input_events": len(seen_events),
            "accepted_events": len(event_payloads),
            "retained_terminal_events_excluded": len(seen_events) - len(event_payloads),
            "truth_pairs": truth_pair_total,
            "reco_pairs": reco_pair_total,
            "matched_pairs": matched_total,
            "fake_or_boundary_reco_pairs": fake_total,
            "missed_truth_pairs": miss_total,
            "category_counts": dict(sorted(category_counts.items())),
            "category_sumw": dict(sorted(category_sumw.items())),
        },
    }
    bundle = ResponseBundle(
        system=config.system,
        dimension=config.dimension,
        truth_ptgamma_edges=CLASSIFICATION_TRUTH_PTGAMMA_EDGES.copy(),
        reco_ptgamma_edges=CLASSIFICATION_RECO_PTGAMMA_EDGES.copy(),
        xj_edges=(COMMON_XJ_EDGES.copy() if config.dimension == "2D" else None),
        matrix=accumulators.matrix,
        matrix_sumw2=accumulators.matrix_sumw2,
        truth=accumulators.truth,
        truth_sumw2=accumulators.truth_sumw2,
        reco=accumulators.reco,
        reco_sumw2=accumulators.reco_sumw2,
        misses=accumulators.misses,
        misses_sumw2=accumulators.misses_sumw2,
        fakes=accumulators.fakes,
        fakes_sumw2=accumulators.fakes_sumw2,
        fake_causes=accumulators.fake_causes,
        fake_causes_sumw2=accumulators.fake_causes_sumw2,
        boundary_reco=accumulators.boundary_reco,
        boundary_reco_sumw2=accumulators.boundary_reco_sumw2,
        boundary_truth=accumulators.boundary_truth,
        boundary_truth_sumw2=accumulators.boundary_truth_sumw2,
        unresolved={
            "golden_mini_dst_differential": "NOT_RUN",
            "full_statistics_response_equivalence": "NOT_RUN",
            "producer_wrong_object_witness_semantics": "NOT_RUN",
        },
        provenance=provenance,
        status=_STATUS,
    )
    conservation = bundle.conservation()
    if not conservation["closes_rtol_1e-10_atol_1e-10"]:
        raise ValueError(f"response partition does not conserve exactly: {conservation}")
    return bundle


@contextmanager
def _exclusive_output_locks(outputs: Iterable[Path]) -> Iterator[None]:
    """Fail closed when another cooperating writer owns any output path."""

    handles = []
    try:
        for output in sorted({Path(path).resolve() for path in outputs}, key=str):
            lock_path = output.with_name(f".{output.name}.lock")
            lock_path.parent.mkdir(parents=True, exist_ok=True)
            handle = lock_path.open("a+b")
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                handle.close()
                raise RuntimeError(
                    f"another response writer owns output path: {output}"
                ) from error
            handles.append(handle)
        yield
    finally:
        for handle in reversed(handles):
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()


def write_response_artifacts(
    input_paths: Iterable[Path],
    output_stem: Path,
    config: ResponseBuildConfig,
    *,
    receipt_path: Path | None = None,
) -> dict[str, Any]:
    """Build deterministic artifacts and publish the receipt as commit marker.

    The two payloads are replaced before the receipt. Any prior receipt is
    invalidated before that replacement starts, so an interrupted publication
    can leave payload files but can never leave a stale receipt certifying
    them. A custom receipt directory is supported without relying on a
    cross-filesystem rename.
    """

    paths = [Path(path).resolve() for path in input_paths]
    stem = Path(output_stem).resolve()
    npz_path = stem.with_suffix(".npz")
    metadata_path = stem.with_suffix(".json")
    receipt = (
        Path(receipt_path).resolve()
        if receipt_path is not None
        else stem.with_suffix(".receipt.json")
    )
    outputs = {npz_path, metadata_path, receipt}
    if len(outputs) != 3:
        raise ValueError("response output paths collide")
    if outputs & set(paths):
        raise ValueError("response output cannot overwrite an input")
    stem.parent.mkdir(parents=True, exist_ok=True)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with _exclusive_output_locks(outputs):
        bundle = build_response(paths, config)
        with tempfile.TemporaryDirectory(dir=stem.parent, prefix=".response-build-") as temporary:
            temporary_stem = Path(temporary) / "response"
            temporary_npz, temporary_metadata = bundle.save(temporary_stem)

            def destination_artifact(
                source: Path,
                destination: Path,
                role: str,
            ) -> dict[str, Any]:
                record = artifact(source, role)
                record["name"] = destination.name
                return record

            payload = {
                "schema": "PhotonJetResponseBuildReceiptV1",
                "state": _STATUS,
                "configuration": config.to_dict(),
                "tree_contract_sha256": bundle.provenance["tree_contract_sha256"],
                "selection_program_sha256": bundle.provenance[
                    "selection_program_sha256"
                ],
                "inputs": bundle.provenance["inputs"],
                "outputs": [
                    destination_artifact(
                        temporary_npz,
                        npz_path,
                        "response_arrays",
                    ),
                    destination_artifact(
                        temporary_metadata,
                        metadata_path,
                        "response_metadata",
                    ),
                ],
                "diagnostics": bundle.diagnostics(),
                "unresolved": bundle.unresolved,
            }
            with tempfile.TemporaryDirectory(
                dir=receipt.parent,
                prefix=".response-receipt-",
            ) as receipt_temporary:
                temporary_receipt = Path(receipt_temporary) / receipt.name
                write_json(temporary_receipt, payload)

                # The receipt is the publication commit marker. Remove any old
                # marker before replacing either payload so a partial failure is
                # visibly uncertified rather than falsely certified.
                receipt.unlink(missing_ok=True)
                os.replace(temporary_npz, npz_path)
                os.replace(temporary_metadata, metadata_path)
                os.replace(temporary_receipt, receipt)
    return payload


__all__ = [
    "ResponseBuildConfig",
    "build_response",
    "write_response_artifacts",
]
