"""Training recipe: weighting, low-ET flattening, split, and label mapping."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

import train_photon_bdt as trainer
from photon_selection import load_config
from synthetic_trees import match_link, signal_truth_photon, write_tree_file

ROOT = Path(__file__).resolve().parents[1]
CONFIG = load_config(ROOT / "config" / "nominal.yaml")
HAVE_SKLEARN = importlib.util.find_spec("sklearn") is not None


class WeightingTest(unittest.TestCase):
    def test_inverse_pdf_weights_are_mean_one_and_flatten(self) -> None:
        rng = np.random.default_rng(1)
        values = rng.normal(0.0, 0.3, size=5000)
        weights = trainer.inverse_pdf_weights(values, n_bins=20, fixed_range=(-0.7, 0.7))
        self.assertAlmostEqual(weights.mean(), 1.0, places=12)
        self.assertGreater(weights[np.abs(values) > 0.5].mean(), weights[np.abs(values) < 0.1].mean())
        capped = trainer.inverse_pdf_weights(rng.exponential(5.0, size=5000) + 6.0, n_bins=20, cap=800.0)
        self.assertAlmostEqual(capped.mean(), 1.0, places=12)
        with self.assertRaises(ValueError):
            trainer.inverse_pdf_weights(np.arange(5.0), n_bins=20)

    def test_class_balance_survives_flattening(self) -> None:
        rng = np.random.default_rng(2)
        labels = np.concatenate([np.ones(3000, dtype=int), np.zeros(1000, dtype=int)])
        et = rng.uniform(6.0, 40.0, size=4000)
        eta = rng.uniform(-0.7, 0.7, size=4000)
        weights, steps = trainer.training_weights(et, eta, labels, CONFIG["systems"]["pp"]["training"]["weighting"])
        # Each flattening factor has mean one within its class, so the class
        # balance of the first step survives up to the correlation of the two
        # factors (as in the maintained recipe).
        ratio = weights[labels == 1].sum() / weights[labels == 0].sum()
        self.assertLess(abs(ratio - 1.0), 0.05)
        self.assertLess(abs(weights.sum() / 4000.0 - 1.0), 0.05)
        self.assertEqual(len(steps), 4)

    def test_low_et_background_flattening(self) -> None:
        rng = np.random.default_rng(3)
        et = np.concatenate([rng.uniform(6.0, 14.9, size=2000), rng.uniform(15.0, 40.0, size=500)])
        labels = np.concatenate([np.zeros(2000, dtype=int), np.ones(500, dtype=int)])
        keep = trainer.flatten_low_et_background(et, labels, CONFIG["systems"]["pp"]["training"]["low_et_background_flattening"])
        self.assertTrue(keep[labels == 1].all())
        edges = np.linspace(et[labels == 0].min(), 15.0, 21)
        counts = np.bincount(np.clip(np.searchsorted(edges, et[(labels == 0) & keep], side="right") - 1, 0, 19), minlength=20)
        self.assertEqual(counts.min(), counts.max())
        self.assertTrue(trainer.flatten_low_et_background(et, labels, {"enabled": False}).all())

    @unittest.skipUnless(HAVE_SKLEARN, "scikit-learn not installed")
    def test_split_fractions(self) -> None:
        labels = np.repeat([0, 1], 500)
        parts = trainer.split_rows(labels, {"test_fraction": 0.2, "validation_fraction": 0.1, "seed": 42})
        self.assertEqual(int((parts == "test").sum()), 200)
        self.assertEqual(int((parts == "validation").sum()), 100)
        self.assertEqual(int((parts == "train").sum()), 700)


class LabelMappingTest(unittest.TestCase):
    def test_signal_file_candidates_without_link_are_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = {"source_file_index": 0, "event_id_hi": 1, "event_id_lo": 1, "run": 1, "event_weight": 1.0,
                    "centrality": -1.0, "vertex_z": 0.0, "terminal_status": 0}
            inputs = {f"bdt_input_{i:02d}": 0.5 for i in range(11)}
            photons = [
                {**base, "candidate_id_hi": 5, "candidate_id_lo": 11, "photon_et": 21.0, "photon_eta": 0.1, "bdt_input_count": 11, **inputs},
                {**base, "candidate_id_hi": 5, "candidate_id_lo": 12, "photon_et": 18.0, "photon_eta": 0.2, "bdt_input_count": 11, **inputs},
                {**base, "candidate_id_hi": 5, "candidate_id_lo": 13, "photon_et": 18.0, "photon_eta": 0.2, "bdt_input_count": 10, **inputs},
            ]
            path = write_tree_file(Path(temporary) / "signal.root", events=[base], photons=photons,
                                   truth_photons=[signal_truth_photon(0, 1, 1, 9, 1, 22.0)],
                                   links=[match_link(0, 1, 1, 5, 11, 9, 1)])
            signal = trainer.load_rows(path, system="pp", config=CONFIG, role="signal")
            background = trainer.load_rows(path, system="pp", config=CONFIG, role="background")
            self.assertEqual(signal["label"].tolist(), [1])
            self.assertEqual(background["label"].tolist(), [0])          # the unlinked complete candidate only
            self.assertEqual(signal["n_complete"], 2)                    # the 10-input candidate is incomplete
            self.assertEqual(signal["n_linked"], 1)


if __name__ == "__main__":
    unittest.main()
