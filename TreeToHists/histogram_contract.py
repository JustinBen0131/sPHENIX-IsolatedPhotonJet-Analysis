"""The histogram package: what TreeToHists writes and FinalAnalysis reads.

One package holds the sufficient statistics of the implemented analysis
contract. Unsupported covariance/measurement extensions must fail explicitly;
FinalAnalysis never returns to a tree. Arrays are numpy, saved as one ``.npz``
beside a ``.json`` of metadata. Nothing here computes physics; this module
only names the arrays, fixes their shapes and checks them.

Grids (from config/measurement.yaml)
    pt_edges          measurement photon-pT bins, n_pt
    xj_edges          xJ bins, n_xj
    reco_pt_edges     response classification grid, reconstructed, n_reco
    truth_pt_edges    response classification grid, truth, n_truth

Arrays (every count array has a matching ``<name>_sumw2``)

  Data and simulation
    abcd_counts                (4, n_pt)         event-leading photons per region A,B,C,D
    abcd_events                (4, n_pt) int     unweighted event counts
    recoil_spectra             (2, n_pt, n_xj)   recoil xJ of the region-A and region-C leader
    photon_et_denominator      (n_pt,)           every kinematically selected photon (no ID), for efficiencies
    bootstrap_abcd_counts      (R, 4, n_pt)      event-bootstrap replicas
    bootstrap_recoil_spectra   (R, 2, n_pt, n_xj)

  Simulation only (truth-conditioned, same selection)
    leakage_abcd_counts        (4, n_pt)         event-leading truth-signal photons per region
    truth_photons              (n_truth,)        truth-signal photons on the truth grid
    photon_response            (n_truth, n_reco) truth pT x reco ET of the region-A leader, matched
    photon_misses              (n_truth,)
    photon_fakes               (n_reco,)         region-A leader without a truth match
    photon_boundary_fakes      (n_reco,)         matched, truth off the truth grid
    pair_truth                 (n_truth * n_xj,) truth photon-jet pairs, flattened g = i_pt * n_xj + i_xj
    pair_reco                  (n_reco * n_xj,)
    pair_response              (n_truth * n_xj, n_reco * n_xj)
    pair_misses                (n_truth * n_xj,)
    pair_fakes_<cause>         (n_reco * n_xj,)  UNMATCHED_RECO, COMBINATORIC
    pair_boundary_reco_<key>   (n_reco * n_xj,)  matched with truth off grid, and ASSIGNED_HARD_NONFIDUCIAL
    pair_boundary_truth_<key>  (n_truth * n_xj,) matched with reco off grid

Metadata (json)
    schema, family, system, measurement_version, selection, grids, inputs
    (paths, hashes, sample, weight receipts), counts, model, working points,
    truth_signal_only flag, bootstrap seed and replica count.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA = "PhotonJetHistogramPackageV1"
REGION_ORDER = ("A", "B", "C", "D")
RECOIL_REGIONS = ("A", "C")

# Response categories (names shared with FinalAnalysis).
UNMATCHED_RECO = "UNMATCHED_RECO"
COMBINATORIC = "COMBINATORIC"
ASSIGNED_HARD_NONFIDUCIAL = "ASSIGNED_HARD_NONFIDUCIAL"
PHOTON_RECO_MISS = "PHOTON_RECO_MISS"


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

    def flatten(self, i_pt: int, i_xj: int) -> int:
        return i_pt * self.n_xj + i_xj

    def to_dict(self) -> dict[str, list[float]]:
        return {name: getattr(self, name).tolist() for name in ("pt_edges", "xj_edges", "reco_pt_edges", "truth_pt_edges")}

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> "Grids":
        return cls(
            pt_edges=config["photon"]["pt_edges_gev"], xj_edges=config["recoil"]["xj_edges"],
            reco_pt_edges=config["response"]["reco_ptgamma_edges"], truth_pt_edges=config["response"]["truth_ptgamma_edges"],
        )


def axis_bin(value: float, edges: np.ndarray) -> int | None:
    """Bin index of value on edges, None when outside or non-finite."""

    if not np.isfinite(value):
        return None
    index = int(np.searchsorted(edges, value, side="right") - 1)
    return index if 0 <= index < len(edges) - 1 else None


@dataclass
class Package:
    grids: Grids
    arrays: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    sample_arrays: dict[str, dict[str, np.ndarray]] = field(default_factory=dict)

    # ---- construction --------------------------------------------------------------

    @classmethod
    def empty(cls, grids: Grids, *, simulation: bool, replicas: int) -> "Package":
        if replicas < 2 or replicas > 10000:
            raise ValueError("bootstrap replicas must be in [2, 10000]")
        g = grids
        cells = replicas * (4*g.n_pt + 2*g.n_pt*g.n_xj) + 2*g.n_truth*g.n_reco*g.n_xj**2
        if cells * 8 > 256 * 1024**2:
            raise ValueError("single package would exceed 256 MiB; reduce explicit fanout")
        package = cls(grids=g)
        a = package.arrays
        a["abcd_counts"] = np.zeros((4, g.n_pt)); a["abcd_counts_sumw2"] = np.zeros((4, g.n_pt))
        a["abcd_events"] = np.zeros((4, g.n_pt), dtype=np.int64)
        a["recoil_spectra"] = np.zeros((2, g.n_pt, g.n_xj)); a["recoil_spectra_sumw2"] = np.zeros((2, g.n_pt, g.n_xj))
        a["photon_et_denominator"] = np.zeros(g.n_pt); a["photon_et_denominator_sumw2"] = np.zeros(g.n_pt)
        a["bootstrap_abcd_counts"] = np.zeros((replicas, 4, g.n_pt))
        a["bootstrap_recoil_spectra"] = np.zeros((replicas, 2, g.n_pt, g.n_xj))
        for name, size in (("accepted_event_weight", 1), ("combinatoric_photon_denominator", g.n_pt)):
            a[name] = np.zeros(size); a[f"{name}_sumw2"] = np.zeros(size)
        if simulation:
            nt, nr, nx = g.n_truth, g.n_reco, g.n_xj
            for name, shape in (
                ("leakage_abcd_counts", (4, g.n_pt)),
                ("truth_photons", (nt,)), ("photon_response", (nt, nr)), ("photon_misses", (nt,)),
                ("photon_reco", (nr,)), ("photon_fakes", (nr,)), ("photon_boundary_fakes", (nr,)),
                ("pair_truth", (nt * nx,)), ("pair_reco", (nr * nx,)), ("pair_response", (nt * nx, nr * nx)),
                ("pair_misses", (nt * nx,)),
                (f"pair_fakes_{UNMATCHED_RECO}", (nr * nx,)), (f"pair_fakes_{COMBINATORIC}", (nr * nx,)),
                (f"pair_boundary_reco_{ASSIGNED_HARD_NONFIDUCIAL}", (nr * nx,)),
            ):
                a[name] = np.zeros(shape); a[f"{name}_sumw2"] = np.zeros(shape)
        return package

    def component(self, name: str, shape: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
        """Return (values, sumw2) of a named array, creating it when absent."""

        if name not in self.arrays:
            self.arrays[name] = np.zeros(shape); self.arrays[f"{name}_sumw2"] = np.zeros(shape)
        return self.arrays[name], self.arrays[f"{name}_sumw2"]

    def fill(self, name: str, index, weight: float) -> None:
        if not np.isfinite(weight) or weight < 0:
            raise ValueError("this count contract requires finite nonnegative weights")
        self.arrays[name][index] += weight
        self.arrays[f"{name}_sumw2"][index] += weight * weight

    # ---- checks --------------------------------------------------------------------

    def check(self) -> dict[str, Any]:
        a = self.arrays
        report: dict[str, Any] = {}
        if not a or "abcd_counts" not in a or "bootstrap_abcd_counts" not in a:
            raise ValueError("incomplete histogram package")
        expected = Package.empty(self.grids, simulation="pair_response" in a,
                                 replicas=a["bootstrap_abcd_counts"].shape[0]).arrays
        for name, values in expected.items():
            if name not in a or a[name].shape != values.shape:
                raise ValueError(f"missing or malformed required array {name}")
        for name, values in a.items():
            if np.any(~np.isfinite(values)) or np.any(values < -1.0e-12):
                raise ValueError(f"{name} has a non-finite or negative entry")
            if name.endswith("_sumw2") or name.startswith("bootstrap_") or name == "abcd_events":
                continue
            if f"{name}_sumw2" not in a:
                raise ValueError(f"{name} has no sumw2 companion")
            if a[f"{name}_sumw2"].shape != values.shape:
                raise ValueError(f"{name}_sumw2 shape differs from {name}")
        if "pair_response" in a:
            g = self.grids
            reco_partition = a["pair_response"].sum(axis=0) + self.all_pair_fakes()
            truth_partition = a["pair_response"].sum(axis=1) + self.all_pair_misses()
            report["pair_reco_closes"] = bool(np.allclose(a["pair_reco"], reco_partition, rtol=1e-9, atol=1e-9))
            report["pair_truth_closes"] = bool(np.allclose(a["pair_truth"], truth_partition, rtol=1e-9, atol=1e-9))
            report["photon_truth_closes"] = bool(np.allclose(a["truth_photons"], a["photon_response"].sum(axis=1) + a["photon_misses"], rtol=1e-9, atol=1e-9))
            report["photon_reco_closes"] = bool(np.allclose(a["photon_reco"], a["photon_response"].sum(axis=0) + a["photon_fakes"] + a["photon_boundary_fakes"], rtol=1e-9, atol=1e-9))
            if not all(report.values()):
                raise ValueError(f"response partition does not close: {report}")
            report["pair_response_shape"] = [g.n_truth * g.n_xj, g.n_reco * g.n_xj]
        return report

    def all_pair_fakes(self) -> np.ndarray:
        total = np.zeros_like(self.arrays["pair_reco"])
        for name, values in self.arrays.items():
            if (name.startswith("pair_fakes_") or name.startswith("pair_boundary_reco_")) and not name.endswith("_sumw2"):
                total = total + values
        return total

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

    # ---- persistence ---------------------------------------------------------------

    def save(self, stem: Path) -> tuple[Path, Path]:
        """Write additive arrays first; a hash-bound JSON receipt completes the package."""
        stem = Path(stem)
        npz, meta = stem.with_suffix(".npz"), stem.with_suffix(".json")
        if npz.exists() or meta.exists():
            raise FileExistsError(f"refusing to overwrite {stem}")
        check = self.check()
        stem.parent.mkdir(parents=True, exist_ok=True)
        arrays = {**self.arrays, **{f"grid_{k}": v for k, v in self.grids.to_dict().items()}}
        for sample, subset in self.sample_arrays.items():
            for name, values in subset.items():
                arrays[f"sample::{sample}::{name}"] = values
        np.savez_compressed(npz, **{k: np.asarray(v) for k, v in arrays.items()})
        payload = {**self.metadata, "schema": SCHEMA, "grids": self.grids.to_dict(),
                   "arrays": sorted(self.arrays), "samples": sorted(self.sample_arrays),
                   "check": check, "npz_sha256": file_sha256(npz)}
        meta.write_text(json.dumps(payload, indent=2, sort_keys=True, default=float, allow_nan=False) + "\n", encoding="utf-8")
        return npz, meta

    @classmethod
    def load(cls, stem: Path) -> "Package":
        stem = Path(stem)
        meta = json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
        if meta.get("schema") != SCHEMA or meta.get("npz_sha256") != file_sha256(stem.with_suffix(".npz")):
            raise ValueError(f"{stem}: package schema or array hash differs")
        grids = Grids(**meta["grids"])
        samples = {}
        with np.load(stem.with_suffix(".npz"), allow_pickle=False) as data:
            arrays = {name: np.asarray(data[name]) for name in data.files
                      if not name.startswith(("grid_", "sample::"))}
            for name, edges in grids.to_dict().items():
                if not np.array_equal(data[f"grid_{name}"], edges):
                    raise ValueError(f"{stem}: JSON and array grids disagree")
            for name in data.files:
                if name.startswith("sample::"):
                    _, sample, field_name = name.split("::", 2)
                    samples.setdefault(sample, {})[field_name] = np.asarray(data[name])
        if sorted(arrays) != meta["arrays"] or sorted(samples) != meta["samples"]:
            raise ValueError("package inventory differs from its receipt")
        package = cls(grids=grids, arrays=arrays, sample_arrays=samples,
                      metadata={k: v for k, v in meta.items() if k not in
                                ("schema", "grids", "arrays", "samples", "check", "npz_sha256")})
        package.check()
        if samples:
            combined = {}
            for sample, subset in samples.items():
                cls(grids, subset).check()
                for name, values in subset.items():
                    combined[name] = combined.get(name, np.zeros_like(values)) + values
            if set(combined) != set(arrays) or any(not np.allclose(v, arrays[k], rtol=1e-12, atol=1e-12)
                                                  for k, v in combined.items()):
                raise ValueError("sample partitions do not reproduce family totals")
        return package


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


__all__ = ["ASSIGNED_HARD_NONFIDUCIAL", "COMBINATORIC", "Grids", "PHOTON_RECO_MISS", "Package",
           "RECOIL_REGIONS", "REGION_ORDER", "SCHEMA", "UNMATCHED_RECO", "axis_bin"]
