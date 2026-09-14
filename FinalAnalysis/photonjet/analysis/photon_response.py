"""Per-event photon response for the photon denominator.

The recoil spectrum is normalised by the number of prompt photons, and that
number is unfolded with its own one-dimensional response: truth photon pT
against the reconstructed ET of the event-leading Region-A photon.  This
module builds that response from photon+jet simulation trees, one entry per
event, on the classification grids used by the pair response
(truth 5-40 GeV, reconstructed 10-40 GeV); the analysis window 15-35 GeV is
cut out later by ``chain.restrict``.

Per event
  truth photon   the unique truth photon passing the nominal truth-signal
                 contract (an event with two is an error);
  reco photon    the event-leading Region-A candidate chosen by the same
                 selection program as the pair response;
  matched        the reconstruction link of that candidate points to that
                 truth photon.

Fills (weight w = complete event weight)
  matched, both on grid          matrix[truth, reco]
  matched, reco off grid         misses[truth]
  matched, truth off grid        boundary_fakes[reco]
  truth without a matched reco   misses[truth]
  reco without a match           fakes[reco]      (removed from data by the ABCD purity)

so that truth = matrix.sum(1) + misses and reco = matrix.sum(0) + fakes + boundary_fakes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import numpy as np
import uproot

from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.response import (
    CLASSIFICATION_RECO_PTGAMMA_EDGES,
    CLASSIFICATION_TRUTH_PTGAMMA_EDGES,
)
from photonjet.analysis.response_builder import (
    _EVENT_BRANCHES,
    _PHOTON_BRANCHES,
    _candidate_record,
    _event_key,
    _finite,
    _group,
    _identity,
    _rows,
    _typed_links,
)
from photonjet.analysis.selection import compile_recoil_selection
from photonjet.provenance import artifact

EventKey = tuple[int, int, int]
LINK_BRANCHES = (
    "source_file_index", "event_id_hi", "event_id_lo",
    "reco_id_hi", "reco_id_lo", "truth_id_hi", "truth_id_lo",
    "reco_type", "truth_type", "link_class",
)


@dataclass
class PhotonResponse:
    system: str
    truth_edges: np.ndarray
    reco_edges: np.ndarray
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
    boundary_fakes: np.ndarray
    boundary_fakes_sumw2: np.ndarray
    provenance: dict[str, Any] = field(default_factory=dict)

    ARRAYS = ("truth_edges", "reco_edges", "matrix", "matrix_sumw2", "truth", "truth_sumw2",
              "reco", "reco_sumw2", "misses", "misses_sumw2", "fakes", "fakes_sumw2",
              "boundary_fakes", "boundary_fakes_sumw2")

    def check_conservation(self) -> None:
        if not np.allclose(self.truth, self.matrix.sum(axis=1) + self.misses, rtol=1e-9, atol=1e-9):
            raise ValueError("photon response: truth is not matrix + misses")
        if not np.allclose(self.reco, self.matrix.sum(axis=0) + self.fakes + self.boundary_fakes, rtol=1e-9, atol=1e-9):
            raise ValueError("photon response: reco is not matrix + fakes + boundary fakes")

    def save(self, stem: Path) -> tuple[Path, Path]:
        stem = Path(stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        npz, meta = stem.with_suffix(".npz"), stem.with_suffix(".json")
        np.savez_compressed(npz, **{name: getattr(self, name) for name in self.ARRAYS})
        meta.write_text(json.dumps({"schema": "PhotonResponseV1", "system": self.system,
                                    "provenance": self.provenance}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz, meta

    @classmethod
    def load(cls, stem: Path) -> "PhotonResponse":
        stem = Path(stem)
        meta = json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
        if meta.get("schema") != "PhotonResponseV1":
            raise ValueError(f"{stem}: not a PhotonResponseV1 file")
        with np.load(stem.with_suffix(".npz")) as data:
            arrays = {name: np.asarray(data[name]) for name in cls.ARRAYS}
        response = cls(system=meta["system"], provenance=meta.get("provenance", {}), **arrays)
        response.check_conservation()
        return response


def _axis_bin(value: float, edges: np.ndarray) -> int | None:
    index = int(np.searchsorted(edges, value, side="right") - 1)
    return index if 0 <= index < len(edges) - 1 else None


def build_photon_response(
    input_paths: Iterable[Path],
    *,
    system: str,
    selection: RecoilSelection,
    truth_signal: Callable[[Mapping[str, np.ndarray]], np.ndarray],
    event_weights: Mapping[EventKey, float] | None = None,
    truth_edges: np.ndarray = CLASSIFICATION_TRUTH_PTGAMMA_EDGES,
    reco_edges: np.ndarray = CLASSIFICATION_RECO_PTGAMMA_EDGES,
) -> PhotonResponse:
    """Build the per-event photon response from photon+jet simulation trees.

    ``truth_signal`` receives the truthPhotons columns of one event (a mapping
    of branch name to array) and returns a boolean mask; the nominal contract
    is ``PhotonID.photon_selection.truth_signal_mask``.  ``event_weights``
    optionally replaces the stored event weight by the complete analysis
    weight; events absent from it are dropped.
    """

    if system not in {"pp", "auau"}:
        raise ValueError("system must be pp or auau")
    truth_edges = np.asarray(truth_edges, dtype=float)
    reco_edges = np.asarray(reco_edges, dtype=float)
    # Leader over the whole reconstructed grid, so migration across the
    # analysis window edges stays visible to the restriction step.
    program = compile_recoil_selection(RecoilSelection(
        photon_et_min=float(reco_edges[0]), photon_et_max=float(reco_edges[-1]),
        photon_abs_eta_max=selection.photon_abs_eta_max, jet_pt_min=selection.jet_pt_min,
        jet_abs_eta_max=selection.jet_abs_eta_max, delta_phi_min=selection.delta_phi_min,
        jet_radius=selection.jet_radius, region="A", non_tight_definition=selection.non_tight_definition,
        isolation_radius=selection.isolation_radius,
    ))
    nt, nr = len(truth_edges) - 1, len(reco_edges) - 1
    matrix = np.zeros((nt, nr)); matrix_w2 = np.zeros((nt, nr))
    truth = np.zeros(nt); truth_w2 = np.zeros(nt)
    reco = np.zeros(nr); reco_w2 = np.zeros(nr)
    misses = np.zeros(nt); misses_w2 = np.zeros(nt)
    fakes = np.zeros(nr); fakes_w2 = np.zeros(nr)
    boundary = np.zeros(nr); boundary_w2 = np.zeros(nr)
    counts = {"events": 0, "accepted": 0, "with_truth_photon": 0, "with_reco_leader": 0, "matched": 0, "dropped_no_weight": 0}
    inputs = []

    for raw_path in input_paths:
        path = Path(raw_path).resolve()
        inputs.append(artifact(path, "photon_response_input"))
        with uproot.open(path) as root:
            events = _rows(root["events"], _EVENT_BRANCHES)
            photon_groups = _group(_rows(root["photons"], _PHOTON_BRANCHES))
            truth_tree = root["truthPhotons"]
            truth_rows = _rows(truth_tree, tuple(truth_tree.keys()))
            link_groups = _group(_rows(root["recoTruthLinks"], LINK_BRANCHES))
        truth_groups = _group(truth_rows)

        for event_row in events:
            key = _event_key(event_row)
            counts["events"] += 1
            if not program.accepts_event(event_row["terminal_status"]):
                continue
            if event_weights is not None:
                if key not in event_weights:
                    counts["dropped_no_weight"] += 1
                    continue
                weight = float(event_weights[key])
            else:
                weight = _finite(event_row["event_weight"], "event weight")
            if weight < 0:
                raise ValueError(f"event {key} has a negative weight")
            counts["accepted"] += 1

            rows = truth_groups.get(key, [])
            signal_rows: list[dict[str, Any]] = []
            if rows:
                columns = {name: np.asarray([row[name] for row in rows]) for name in rows[0]}
                mask = np.asarray(truth_signal(columns), dtype=bool)
                for row, ok in zip(rows, mask):
                    if ok and abs(_finite(row["truth_photon_eta"], "truth eta")) < selection.photon_abs_eta_max:
                        signal_rows.append(row)
            if len(signal_rows) > 1:
                raise ValueError(f"event {key} has more than one truth-signal photon")
            truth_row = signal_rows[0] if signal_rows else None
            truth_bin = _axis_bin(_finite(truth_row["truth_photon_pt"], "truth pT"), truth_edges) if truth_row else None
            if truth_row is not None:
                counts["with_truth_photon"] += 1

            candidates = photon_groups.get(key, [])
            leader = program.choose_leader([_candidate_record(row) for row in candidates]) if candidates else -1
            reco_bin = _axis_bin(_finite(candidates[leader]["photon_et"], "photon ET"), reco_edges) if leader >= 0 else None
            if leader >= 0:
                counts["with_reco_leader"] += 1
            photon_links, _ = _typed_links(link_groups.get(key, []), key)
            matched = (
                leader >= 0 and truth_row is not None
                and photon_links.get(_identity(candidates[leader], "candidate_id")) == _identity(truth_row, "truth_photon_id")
            )
            w2 = weight * weight
            if truth_bin is not None:
                truth[truth_bin] += weight; truth_w2[truth_bin] += w2
            if reco_bin is not None:
                reco[reco_bin] += weight; reco_w2[reco_bin] += w2
            if matched:
                counts["matched"] += 1
                if truth_bin is not None and reco_bin is not None:
                    matrix[truth_bin, reco_bin] += weight; matrix_w2[truth_bin, reco_bin] += w2
                elif truth_bin is not None:
                    misses[truth_bin] += weight; misses_w2[truth_bin] += w2
                elif reco_bin is not None:
                    boundary[reco_bin] += weight; boundary_w2[reco_bin] += w2
            else:
                if truth_bin is not None:
                    misses[truth_bin] += weight; misses_w2[truth_bin] += w2
                if reco_bin is not None:
                    fakes[reco_bin] += weight; fakes_w2[reco_bin] += w2

    response = PhotonResponse(
        system=system, truth_edges=truth_edges, reco_edges=reco_edges,
        matrix=matrix, matrix_sumw2=matrix_w2, truth=truth, truth_sumw2=truth_w2,
        reco=reco, reco_sumw2=reco_w2, misses=misses, misses_sumw2=misses_w2,
        fakes=fakes, fakes_sumw2=fakes_w2, boundary_fakes=boundary, boundary_fakes_sumw2=boundary_w2,
        provenance={"inputs": inputs, "counts": counts, "selection_program_sha256": program.sha256,
                    "weights": "complete analysis weights" if event_weights is not None else "stored event_weight"},
    )
    response.check_conservation()
    return response


__all__ = ["PhotonResponse", "build_photon_response", "LINK_BRANCHES"]
