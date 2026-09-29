"""The histogram package: the named arrays of the simulated response that FinalAnalysis builds from the
FillPhotonJetResponse outputs (closure.package_window, weighted_package) and unfolds (unfolding.restrict). The data side
is the event tables (package.event_table_arrays). Nothing here computes physics; this module only names the arrays and
fixes their shapes.

Grids (the measurement's, FinalAnalysis Settings.grids())
    pt_edges          measurement photon-pT bins, n_pt
    xj_edges          xJ bins, n_xj
    reco_pt_edges     response classification grid, reconstructed, n_reco
    truth_pt_edges    response classification grid, truth, n_truth

Arrays (truth-conditioned simulation; pairs flattened g = i_pt * n_xj + i_xj)
    truth_photons              (n_truth,)        truth-signal photons on the truth grid
    photon_reco                (n_reco,)         region-A leaders
    photon_response            (n_truth, n_reco) truth pT x reco ET of the region-A leader, matched
    photon_misses              (n_truth,)
    photon_fakes               (n_reco,)         region-A leader without a truth match
    photon_boundary_fakes      (n_reco,)         matched, truth off the truth grid
    pair_truth                 (n_truth * n_xj,) truth photon-jet pairs
    pair_reco                  (n_reco * n_xj,)
    pair_response              (n_truth * n_xj, n_reco * n_xj)
    pair_misses                (n_truth * n_xj,)
    pair_fakes_<cause>         (n_reco * n_xj,)  UNMATCHED_RECO, COMBINATORIC (with pair_fakes_COMBINATORIC_sumw2)
    pair_boundary_reco_<key>   (n_reco * n_xj,)  ASSIGNED_HARD_NONFIDUCIAL: detector fakes
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

# Response categories (names shared with FinalAnalysis).
UNMATCHED_RECO = "UNMATCHED_RECO"
COMBINATORIC = "COMBINATORIC"
ASSIGNED_HARD_NONFIDUCIAL = "ASSIGNED_HARD_NONFIDUCIAL"


@dataclass
class Grids:
    pt_edges: np.ndarray
    xj_edges: np.ndarray
    reco_pt_edges: np.ndarray
    truth_pt_edges: np.ndarray

    def __post_init__(self) -> None:
        for name in ("pt_edges", "xj_edges", "reco_pt_edges", "truth_pt_edges"):
            values = np.asarray(getattr(self, name), dtype=float)
            if values.ndim != 1 or len(values) < 2 or np.any(~np.isfinite(values)) or np.any(np.diff(values) <= 0):
                raise ValueError(f"{name} must be strictly increasing with at least two edges")
            setattr(self, name, values)

    @property
    def n_pt(self) -> int:
        return len(self.pt_edges) - 1

    @property
    def n_xj(self) -> int:
        return len(self.xj_edges) - 1

    @property
    def n_reco(self) -> int:
        return len(self.reco_pt_edges) - 1

    @property
    def n_truth(self) -> int:
        return len(self.truth_pt_edges) - 1


@dataclass
class Package:
    grids: Grids
    arrays: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def empty(cls, grids: Grids) -> "Package":
        """Zero arrays of every response quantity on the grids."""

        g = grids
        nt, nr, nx = g.n_truth, g.n_reco, g.n_xj
        package = cls(grids=g)
        for name, shape in (
            ("truth_photons", (nt,)), ("photon_response", (nt, nr)), ("photon_misses", (nt,)),
            ("photon_reco", (nr,)), ("photon_fakes", (nr,)), ("photon_boundary_fakes", (nr,)),
            ("pair_truth", (nt * nx,)), ("pair_reco", (nr * nx,)), ("pair_response", (nt * nx, nr * nx)),
            ("pair_misses", (nt * nx,)),
            (f"pair_fakes_{UNMATCHED_RECO}", (nr * nx,)), (f"pair_fakes_{COMBINATORIC}", (nr * nx,)),
            (f"pair_fakes_{COMBINATORIC}_sumw2", (nr * nx,)),
            (f"pair_boundary_reco_{ASSIGNED_HARD_NONFIDUCIAL}", (nr * nx,)),
        ):
            package.arrays[name] = np.zeros(shape)
        return package

    def all_pair_misses(self) -> np.ndarray:
        total = self.arrays["pair_misses"].copy()
        for name, values in self.arrays.items():
            if name.startswith("pair_boundary_truth_") and not name.endswith("_sumw2"):
                total = total + values
        return total

    def fake_causes(self) -> dict[str, np.ndarray]:
        return {name[len("pair_fakes_"):]: values for name, values in self.arrays.items()
                if name.startswith("pair_fakes_") and not name.endswith("_sumw2")}

    def boundary_reco(self) -> dict[str, np.ndarray]:
        return {name[len("pair_boundary_reco_"):]: values for name, values in self.arrays.items()
                if name.startswith("pair_boundary_reco_") and not name.endswith("_sumw2")}


__all__ = ["ASSIGNED_HARD_NONFIDUCIAL", "COMBINATORIC", "Grids", "Package", "UNMATCHED_RECO"]


