"""Tree-only recoil reduction and an intentionally small offline skim.

The event-oriented tree stores event-local photon, jet, and pair indices.
Using those indices avoids accidental cross-event joins and makes the
event-leading photon definition explicit.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Iterator

import awkward as ak
import numpy as np
import uproot

from photonjet.analysis.selection import SelectionProgram, compile_recoil_selection
from photonjet.provenance import artifact


DEFAULT_XJ_EDGES = tuple(round(0.1 * index, 10) for index in range(19))


@dataclass(frozen=True)
class RecoilSelection:
    photon_et_min: float = 15.0
    photon_et_max: float = 35.0
    photon_abs_eta_max: float = 0.7
    jet_pt_min: float = 5.0
    jet_abs_eta_max: float = 0.7
    delta_phi_min: float = 7.0 * math.pi / 8.0
    jet_radius: float = 0.4
    region: str = "A"
    non_tight_definition: str = "bounded"
    isolation_radius: float = 0.4

    def __post_init__(self) -> None:
        if self.isolation_radius not in {0.3, 0.4}:
            raise ValueError("isolation_radius must be 0.3 or 0.4")
        if self.region not in {"inclusive", "A", "B", "C", "D"}:
            raise ValueError("region must be inclusive, A, B, C, or D")
        if self.non_tight_definition not in {"bounded", "complement"}:
            raise ValueError("non_tight_definition must be bounded or complement")
        if not self.photon_et_min < self.photon_et_max:
            raise ValueError("photon pT interval must be increasing")
        if min(self.photon_abs_eta_max, self.jet_pt_min, self.jet_abs_eta_max, self.jet_radius) <= 0:
            raise ValueError("selection scales must be positive")

    @property
    def recorded_leader_branch(self) -> str | None:
        """Optional V1 witness for the exact canonical photon eligibility."""

        return compile_recoil_selection(self).leader_branch

    @property
    def leader_branch(self) -> str:
        """Return the applicable stored witness or reject an inapplicable use."""

        branch = self.recorded_leader_branch
        if branch is None:
            raise ValueError(
                "selection has no recorded leader branch; use offline recomputation"
            )
        return branch


EVENT_BRANCHES = (
    "source_file_index",
    "event_id_hi",
    "event_id_lo",
    "event_weight",
    "terminal_status",
    "photon_candidate_id_hi",
    "photon_candidate_id_lo",
    "photon_encounter_ordinal",
    "photon_et",
    "photon_eta",
    "photon_bdt_score",
    "photon_bdt_tight_threshold",
    "photon_bdt_nontight_low_threshold",
    "photon_bdt_nontight_high_threshold",
    "photon_iso_r03",
    "photon_iso_r03_threshold",
    "photon_iso_r03_nonisolated_threshold",
    "photon_iso_r04",
    "photon_iso_r04_threshold",
    "photon_iso_r04_nonisolated_threshold",
    "jet_id_hi",
    "jet_id_lo",
    "jet_pt",
    "jet_eta",
    "jet_radius",
    "pair_photon_index",
    "pair_jet_index",
    "pair_delta_phi",
    "pair_xjgamma",
)


def _event_arrays(path: Path, leader_branch: str | None) -> ak.Array:
    with uproot.open(path) as root:
        if "eventTree" not in root:
            raise ValueError(f"{path} does not contain PhotonJetTrees_v1/eventTree")
        tree = root["eventTree"]
        required = [*EVENT_BRANCHES]
        if leader_branch is not None:
            required.append(leader_branch)
        missing = sorted(set(required) - set(tree.keys()))
        if missing:
            raise ValueError(f"{path} eventTree lacks branches: {', '.join(missing)}")
        return tree.arrays(required, library="ak")


def _selected_rows(
    input_paths: Iterable[Path],
    program: SelectionProgram,
) -> Iterator[dict[str, Any]]:
    scope = next(fact.value for fact in program.facts if fact.key == "photon.candidate_scope")
    if scope == "inclusive_pairs":
        yield from _selected_inclusive_rows(input_paths, program)
        return
    for raw_path in input_paths:
        path = Path(raw_path).resolve()
        arrays = _event_arrays(path, program.leader_branch)
        for index in range(len(arrays)):
            if not program.accepts_event(arrays["terminal_status"][index]):
                continue
            candidate_records = [
                {
                    "photon_et": float(arrays["photon_et"][index][local]),
                    "photon_eta": float(arrays["photon_eta"][index][local]),
                    "photon_encounter_ordinal": int(
                        arrays["photon_encounter_ordinal"][index][local]
                    ),
                    "photon_bdt_score": float(arrays["photon_bdt_score"][index][local]),
                    "photon_bdt_tight_threshold": float(
                        arrays["photon_bdt_tight_threshold"][index][local]
                    ),
                    "photon_bdt_nontight_low_threshold": float(
                        arrays["photon_bdt_nontight_low_threshold"][index][local]
                    ),
                    "photon_bdt_nontight_high_threshold": float(
                        arrays["photon_bdt_nontight_high_threshold"][index][local]
                    ),
                    "photon_iso_r03": float(arrays["photon_iso_r03"][index][local]),
                    "photon_iso_r03_threshold": float(
                        arrays["photon_iso_r03_threshold"][index][local]
                    ),
                    "photon_iso_r03_nonisolated_threshold": float(
                        arrays["photon_iso_r03_nonisolated_threshold"][index][local]
                    ),
                    "photon_iso_r04": float(arrays["photon_iso_r04"][index][local]),
                    "photon_iso_r04_threshold": float(
                        arrays["photon_iso_r04_threshold"][index][local]
                    ),
                    "photon_iso_r04_nonisolated_threshold": float(
                        arrays["photon_iso_r04_nonisolated_threshold"][index][local]
                    ),
                }
                for local in range(len(arrays["photon_et"][index]))
            ]
            leader = program.choose_leader(candidate_records)
            if program.leader_branch is not None:
                recorded_leader = int(arrays[program.leader_branch][index])
                if leader != recorded_leader:
                    raise ValueError(
                        f"eventTree {program.leader_branch}={recorded_leader} differs from "
                        f"the executable selection program result {leader} in {path} row {index}"
                    )
            if leader < 0:
                continue
            photon_et = float(arrays["photon_et"][index][leader])
            photon_eta = float(arrays["photon_eta"][index][leader])
            if not program.accepts(
                "photon", {"photon_et": photon_et, "photon_eta": photon_eta}
            ):
                continue
            photon_hi = int(arrays["photon_candidate_id_hi"][index][leader])
            photon_lo = int(arrays["photon_candidate_id_lo"][index][leader])
            for pair_index, (photon_index, jet_index) in enumerate(
                zip(
                    arrays["pair_photon_index"][index],
                    arrays["pair_jet_index"][index],
                )
            ):
                if int(photon_index) != leader:
                    continue
                jet_local = int(jet_index)
                jet_pt = float(arrays["jet_pt"][index][jet_local])
                jet_eta = float(arrays["jet_eta"][index][jet_local])
                jet_radius = float(arrays["jet_radius"][index][jet_local])
                delta_phi = float(arrays["pair_delta_phi"][index][pair_index])
                if not program.accepts(
                    "recoil",
                    {
                        "jet_pt": jet_pt,
                        "jet_eta": jet_eta,
                        "jet_radius": jet_radius,
                        "delta_phi": delta_phi,
                    },
                ):
                    continue
                yield {
                    "source_file_index": int(arrays["source_file_index"][index]),
                    "event_id_hi": int(arrays["event_id_hi"][index]),
                    "event_id_lo": int(arrays["event_id_lo"][index]),
                    "photon_id_hi": photon_hi,
                    "photon_id_lo": photon_lo,
                    "jet_id_hi": int(arrays["jet_id_hi"][index][jet_local]),
                    "jet_id_lo": int(arrays["jet_id_lo"][index][jet_local]),
                    "photon_index": leader,
                    "jet_index": jet_local,
                    "photon_et": photon_et,
                    "photon_eta": photon_eta,
                    "photon_bdt_score": candidate_records[leader]["photon_bdt_score"],
                    "photon_bdt_tight_threshold": candidate_records[leader][
                        "photon_bdt_tight_threshold"
                    ],
                    "photon_bdt_nontight_low_threshold": candidate_records[leader][
                        "photon_bdt_nontight_low_threshold"
                    ],
                    "photon_bdt_nontight_high_threshold": candidate_records[leader][
                        "photon_bdt_nontight_high_threshold"
                    ],
                    "photon_iso_r04": candidate_records[leader]["photon_iso_r04"],
                    "photon_iso_r04_threshold": candidate_records[leader][
                        "photon_iso_r04_threshold"
                    ],
                    "photon_iso_r04_nonisolated_threshold": candidate_records[leader][
                        "photon_iso_r04_nonisolated_threshold"
                    ],
                    "jet_pt": jet_pt,
                    "jet_eta": jet_eta,
                    "jet_radius": jet_radius,
                    "delta_phi": delta_phi,
                    "xjgamma": float(arrays["pair_xjgamma"][index][pair_index]),
                    "event_weight": float(arrays["event_weight"][index]),
                }


def _selected_inclusive_rows(
    input_paths: Iterable[Path],
    program: SelectionProgram,
) -> Iterator[dict[str, Any]]:
    """Match the direct flat photon-jet loop without choosing a photon leader."""

    branches = [
        "source_file_index", "event_id_hi", "event_id_lo",
        "candidate_id_hi", "candidate_id_lo", "jet_id_hi", "jet_id_lo",
        "photon_index", "jet_index", "photon_et", "photon_eta",
        "jet_pt", "jet_eta", "jet_radius", "delta_phi", "xjgamma",
        "event_weight", "terminal_status",
    ]
    for raw_path in input_paths:
        path = Path(raw_path).resolve()
        with uproot.open(path) as root:
            if "photonJets" not in root:
                raise ValueError(f"{path} does not contain PhotonJetTrees_v1/photonJets")
            tree = root["photonJets"]
            missing = sorted(set(branches) - set(tree.keys()))
            if missing:
                raise ValueError(f"{path} photonJets lacks branches: {', '.join(missing)}")
            arrays = tree.arrays(branches, library="np")
        for index in range(len(arrays["xjgamma"])):
            if not program.accepts_event(arrays["terminal_status"][index]):
                continue
            if not program.accepts(
                "photon",
                {
                    "photon_et": float(arrays["photon_et"][index]),
                    "photon_eta": float(arrays["photon_eta"][index]),
                },
            ):
                continue
            if not program.accepts(
                "recoil",
                {
                    "jet_pt": float(arrays["jet_pt"][index]),
                    "jet_eta": float(arrays["jet_eta"][index]),
                    "jet_radius": float(arrays["jet_radius"][index]),
                    "delta_phi": float(arrays["delta_phi"][index]),
                },
            ):
                continue
            yield {
                "source_file_index": int(arrays["source_file_index"][index]),
                "event_id_hi": int(arrays["event_id_hi"][index]),
                "event_id_lo": int(arrays["event_id_lo"][index]),
                "photon_id_hi": int(arrays["candidate_id_hi"][index]),
                "photon_id_lo": int(arrays["candidate_id_lo"][index]),
                "jet_id_hi": int(arrays["jet_id_hi"][index]),
                "jet_id_lo": int(arrays["jet_id_lo"][index]),
                "photon_index": int(arrays["photon_index"][index]),
                "jet_index": int(arrays["jet_index"][index]),
                "photon_et": float(arrays["photon_et"][index]),
                "photon_eta": float(arrays["photon_eta"][index]),
                "jet_pt": float(arrays["jet_pt"][index]),
                "jet_eta": float(arrays["jet_eta"][index]),
                "jet_radius": float(arrays["jet_radius"][index]),
                "delta_phi": float(arrays["delta_phi"][index]),
                "xjgamma": float(arrays["xjgamma"][index]),
                "event_weight": float(arrays["event_weight"][index]),
            }


def recoil_histogram(
    input_paths: Iterable[Path],
    selection: RecoilSelection = RecoilSelection(),
    edges: Iterable[float] = DEFAULT_XJ_EDGES,
) -> dict[str, Any]:
    paths = [Path(path).resolve() for path in input_paths]
    if not paths:
        raise ValueError("at least one PhotonJetTrees_v1 input is required")
    axis = np.asarray(tuple(edges), dtype=float)
    if axis.ndim != 1 or len(axis) < 2 or np.any(~np.isfinite(axis)) or np.any(np.diff(axis) <= 0):
        raise ValueError("histogram edges must be finite and strictly increasing")
    program = compile_recoil_selection(selection)
    rows = list(_selected_rows(paths, program))
    values = np.asarray([row["xjgamma"] for row in rows], dtype=float)
    weights = np.asarray([row["event_weight"] for row in rows], dtype=float)
    counts, _ = np.histogram(values, bins=axis, weights=weights)
    sumw2, _ = np.histogram(values, bins=axis, weights=weights * weights)
    events = {
        (row["source_file_index"], row["event_id_hi"], row["event_id_lo"])
        for row in rows
    }
    underflow = float(weights[values < axis[0]].sum()) if len(values) else 0.0
    overflow = float(weights[values >= axis[-1]].sum()) if len(values) else 0.0
    return {
        "schema": "PhotonJetHistogramPayloadV1",
        "observable": "xjgamma",
        "selection": asdict(selection),
        "selection_program_sha256": program.sha256,
        "axis": {"edges": axis.tolist(), "underflow": underflow, "overflow": overflow},
        "sumw": counts.tolist(),
        "sumw2": sumw2.tolist(),
        "selected_events": len(events),
        "selected_pairs": len(rows),
    }


SKIM_SCHEMA = {
    "source_file_index": "int32",
    "event_id_hi": "uint64",
    "event_id_lo": "uint64",
    "photon_id_hi": "uint64",
    "photon_id_lo": "uint64",
    "jet_id_hi": "uint64",
    "jet_id_lo": "uint64",
    "photon_index": "int32",
    "jet_index": "int32",
    "photon_et": "float64",
    "photon_eta": "float64",
    "jet_pt": "float64",
    "jet_eta": "float64",
    "jet_radius": "float64",
    "delta_phi": "float64",
    "xjgamma": "float64",
    "event_weight": "float64",
}


def write_recoil_skim(
    input_paths: Iterable[Path],
    output_path: Path,
    selection: RecoilSelection = RecoilSelection(),
) -> dict[str, Any]:
    paths = [Path(path).resolve() for path in input_paths]
    program = compile_recoil_selection(selection)
    rows = list(_selected_rows(paths, program))
    columns = {
        name: np.asarray([row[name] for row in rows], dtype=dtype)
        for name, dtype in SKIM_SCHEMA.items()
    }
    output = Path(output_path).resolve()
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
        # Passing a file-like sink makes uproot serialize no filesystem path in
        # the ROOT header while streaming to disk instead of buffering the skim.
        with uproot.recreate(temporary) as root:
            tree = root.mktree(
                "recoilPairs",
                SKIM_SCHEMA,
                title="Documented event-leading photon and recoil-jet offline skim",
            )
            if rows:
                tree.extend(columns)
        sync_fd = os.open(temporary_name, os.O_RDONLY)
        try:
            os.fsync(sync_fd)
        finally:
            os.close(sync_fd)
        os.replace(temporary_name, output)
        temporary_name = None
    finally:
        if temporary_name is not None and Path(temporary_name).exists():
            Path(temporary_name).unlink()
    return {
        "schema": "PhotonJetRecoilSkimReceiptV1",
        "output": artifact(output, "recoil_skim"),
        "selection": asdict(selection),
        "selection_program_sha256": program.sha256,
        "selected_pairs": len(rows),
        "selected_events": len(
            {
                (row["source_file_index"], row["event_id_hi"], row["event_id_lo"])
                for row in rows
            }
        ),
    }
