"""The correction and unfolding chain on synthetic, self-consistent responses."""

from __future__ import annotations

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from photonjet.analysis.background import PT_EDGES
from photonjet.analysis.chain import gate_and_score, load_bundle, restrict, run_chain, select_iterations
from photonjet.analysis.photon_response import PhotonResponse
from photonjet.analysis.response import (
    CLASSIFICATION_RECO_PTGAMMA_EDGES,
    CLASSIFICATION_TRUTH_PTGAMMA_EDGES,
    COMMON_XJ_EDGES,
    Category,
    ResponseBundle,
)

TRUTH_EDGES = np.asarray(CLASSIFICATION_TRUTH_PTGAMMA_EDGES, dtype=float)
RECO_EDGES = np.asarray(CLASSIFICATION_RECO_PTGAMMA_EDGES, dtype=float)
XJ_EDGES = np.asarray(COMMON_XJ_EDGES, dtype=float)
NX = len(XJ_EDGES) - 1


def truth_shape() -> np.ndarray:
    """A smooth xJ spectrum peaked at 0.85, zero above 2."""

    centres = 0.5 * (XJ_EDGES[:-1] + XJ_EDGES[1:])
    shape = np.exp(-0.5 * ((centres - 0.85) / 0.25) ** 2)
    shape[centres > 2.0] = 0.0
    return shape


def synthetic_bundle(scale: float = 2000.0) -> ResponseBundle:
    nt, nr = len(TRUTH_EDGES) - 1, len(RECO_EDGES) - 1
    matrix = np.zeros((nt, NX, nr, NX))
    shape = truth_shape()
    for t in range(nt):
        same = [r for r in range(nr) if math.isclose(RECO_EDGES[r], TRUTH_EDGES[t])]
        if not same:
            continue
        r = same[0]
        for x in range(NX):
            total = scale * shape[x] / (1.0 + 0.3 * t)
            matrix[t, x, r, x] += 0.7 * total
            if x > 0:
                matrix[t, x, r, x - 1] += 0.08 * total
            if x < NX - 1:
                matrix[t, x, r, x + 1] += 0.08 * total
            # photon-pT migration into the neighbouring reconstructed bins
            if r > 0:
                matrix[t, x, r - 1, x] += 0.07 * total
            if r < nr - 1:
                matrix[t, x, r + 1, x] += 0.07 * total
    flat = matrix.reshape(nt * NX, nr * NX)
    misses = 0.1 * flat.sum(axis=1)
    truth = flat.sum(axis=1) + misses
    fake_causes = {
        Category.UNMATCHED_RECO.value: 0.05 * flat.sum(axis=0),
        Category.COMBINATORIC.value: 0.03 * flat.sum(axis=0),
    }
    fakes = sum(fake_causes.values())
    boundary_reco = {Category.ASSIGNED_HARD_NONFIDUCIAL.value: 0.02 * flat.sum(axis=0)}
    reco = flat.sum(axis=0) + fakes + sum(boundary_reco.values())
    return ResponseBundle(
        system="pp", dimension="2D", truth_ptgamma_edges=TRUTH_EDGES, reco_ptgamma_edges=RECO_EDGES, xj_edges=XJ_EDGES,
        matrix=flat, matrix_sumw2=flat, truth=truth, truth_sumw2=truth, reco=reco, reco_sumw2=reco,
        misses=misses, misses_sumw2=misses, fakes=fakes, fakes_sumw2=fakes,
        fake_causes=fake_causes, fake_causes_sumw2=fake_causes,
        boundary_reco=boundary_reco, boundary_reco_sumw2=boundary_reco,
        provenance={"synthetic": True}, status="SYNTHETIC",
    )


def synthetic_photon_response(scale: float = 5000.0) -> PhotonResponse:
    nt, nr = len(TRUTH_EDGES) - 1, len(RECO_EDGES) - 1
    matrix = np.zeros((nt, nr))
    for t in range(nt):
        same = [r for r in range(nr) if math.isclose(RECO_EDGES[r], TRUTH_EDGES[t])]
        if same:
            matrix[t, same[0]] = scale / (1.0 + 0.3 * t)
            if same[0] + 1 < nr:
                matrix[t, same[0] + 1] = 0.05 * scale / (1.0 + 0.3 * t)
    misses = 0.15 * matrix.sum(axis=1)
    fakes = 0.05 * matrix.sum(axis=0)
    boundary = np.zeros(nr)
    return PhotonResponse(
        system="pp", truth_edges=TRUTH_EDGES, reco_edges=RECO_EDGES, matrix=matrix, matrix_sumw2=matrix,
        truth=matrix.sum(axis=1) + misses, truth_sumw2=matrix.sum(axis=1) + misses,
        reco=matrix.sum(axis=0) + fakes + boundary, reco_sumw2=matrix.sum(axis=0) + fakes + boundary,
        misses=misses, misses_sumw2=misses, fakes=fakes, fakes_sumw2=fakes,
        boundary_fakes=boundary, boundary_fakes_sumw2=boundary,
    )


def synthetic_payloads(window) -> tuple[dict, dict]:
    """Data as a reco-like spectrum with background; leakage from a pure-signal simulation."""

    reco = window.xj_reco
    counts_a = window.photon_reco.copy()
    data = {
        "schema": "PhotonJetNominalHistogramsV1", "system": "pp", "pt_edges": PT_EDGES.tolist(), "xj_edges": XJ_EDGES.tolist(),
        "abcd_event_leading_counts_ptgamma": np.vstack([counts_a, 0.25 * counts_a, 0.6 * counts_a, 0.9 * counts_a]).tolist(),
        "abcd_event_leading_sumw2_ptgamma": np.vstack([counts_a, 0.25 * counts_a, 0.6 * counts_a, 0.9 * counts_a]).tolist(),
        "inclusive_recoil_spectra_ptgamma_xj": np.stack([1.1 * reco, 0.5 * reco]).tolist(),
        "inclusive_recoil_sumw2_ptgamma_xj": np.stack([1.1 * reco, 0.5 * reco]).tolist(),
    }
    leak_a = 10.0 * counts_a
    leakage = {
        **data, "truth_signal_only": True,
        "abcd_event_leading_counts_ptgamma": np.vstack([leak_a, 0.05 * leak_a, 0.02 * leak_a, 0.001 * leak_a]).tolist(),
        "abcd_event_leading_sumw2_ptgamma": np.vstack([leak_a, 0.05 * leak_a, 0.02 * leak_a, 0.001 * leak_a]).tolist(),
    }
    return data, leakage


class RestrictTest(unittest.TestCase):
    def test_window_is_conserved(self) -> None:
        window = restrict(synthetic_bundle(), synthetic_photon_response())
        self.assertEqual(window.xj_response.shape, (3 * NX, 3 * NX))
        np.testing.assert_allclose(window.xj_truth.reshape(-1), window.xj_response.sum(axis=1) + window.xj_misses, rtol=1e-9)
        np.testing.assert_allclose(window.photon_truth, window.photon_response.sum(axis=1) + window.photon_misses, rtol=1e-9)
        self.assertTrue(np.all(window.xj_combinatoric >= 0) and np.all(window.xj_detector_fakes >= 0))
        # feed-in from the 10-15 and 35-40 GeV truth rows becomes boundary fakes
        self.assertGreater(window.xj_boundary_fakes.sum(), 0.0)

    def test_saved_bundle_round_trips(self) -> None:
        bundle = synthetic_bundle()
        with tempfile.TemporaryDirectory() as temporary:
            stem = Path(temporary) / "pairs"
            bundle.save(stem)
            loaded = load_bundle(stem)
            np.testing.assert_allclose(loaded.matrix, bundle.matrix)
            self.assertEqual(set(loaded.fake_causes), set(bundle.fake_causes))
            photon = synthetic_photon_response()
            photon.save(Path(temporary) / "photons")
            reloaded = PhotonResponse.load(Path(temporary) / "photons")
            np.testing.assert_allclose(reloaded.matrix, photon.matrix)


class ChainTest(unittest.TestCase):
    def setUp(self) -> None:
        self.window = restrict(synthetic_bundle(), synthetic_photon_response())
        self.data, self.leakage = synthetic_payloads(self.window)

    def test_run_chain_is_finite_and_closes(self) -> None:
        result = run_chain(self.data, self.leakage, self.window, iterations=3, subtract_combinatoric=True,
                           purity_strategy="per_pt", toys=25, seed=3)
        density = np.asarray(result["density"])
        self.assertEqual(density.shape, (NX,))
        self.assertTrue(np.all(np.isfinite(density)))
        self.assertGreater(result["toy_successes"], 0)
        self.assertTrue(math.isfinite(result["refold_chi2_ndf"]))
        self.assertTrue(math.isfinite(result["mc_closure_chi2_ndf"]))
        self.assertLess(result["mc_closure_chi2_ndf"], 1.0)          # closure of a self-consistent response
        self.assertEqual(len(result["purity"]), 3)
        self.assertEqual(np.asarray(result["density_stat_covariance"]).shape, (NX, NX))
        self.assertEqual(np.asarray(result["density_per_ptgamma"]).shape, (3, NX))

    def test_gate_and_score(self) -> None:
        result = run_chain(self.data, self.leakage, self.window, iterations=3, subtract_combinatoric=False,
                           purity_strategy="per_pt", toys=25, seed=5)
        passed, score, metrics = gate_and_score(result)
        self.assertTrue(math.isfinite(score))
        self.assertIn("toy_success_fraction", metrics)
        self.assertIsInstance(passed, bool)

    def test_select_iterations_returns_a_table(self) -> None:
        result = select_iterations(self.data, self.leakage, self.window, candidates=(2, 3), scan_toys=8, final_toys=12,
                                   seed=1, subtract_combinatoric=False, purity_strategy="per_pt")
        self.assertEqual([row["iterations"] for row in result["selection"]["table"]], [2, 3])
        self.assertIn(result["selection"]["selected_iterations"], (2, 3))
        self.assertEqual(result["toy_requested"], 12)


if __name__ == "__main__":
    unittest.main()
