"""The complete per-event weight: order, ownership windows, centrality support."""

from __future__ import annotations

import math
from pathlib import Path
import unittest

import numpy as np

from sample_weights import complete_event_weights, leading_truth_jet_pt, load_manifest

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = load_manifest(ROOT / "config" / "samples.yaml")


class ManifestTest(unittest.TestCase):
    def test_source_factors(self) -> None:
        self.assertAlmostEqual(MANIFEST.samples["auau_jet20"].factor_value, 38811.785 / 6564888, places=18)
        self.assertAlmostEqual(MANIFEST.samples["auau_photon12"].factor_value, 2598.12425 / 9991946, places=18)
        self.assertEqual(MANIFEST.samples["auau_jet40"].ownership_window, (41.0, None))
        self.assertIsNone(MANIFEST.samples["pp_photonjet"].factor_value)      # generated_events not recorded
        self.assertEqual(MANIFEST.samples["pp_data"].kind, "data")

    def test_centrality_maps(self) -> None:
        self.assertEqual(len(MANIFEST.centrality.edges), 17)
        for family in ("photonjet", "inclusivejet"):
            self.assertEqual(MANIFEST.centrality.values[family].shape, (16,))
            self.assertTrue(np.all(MANIFEST.centrality.values[family] > 0))


class CompleteWeightTest(unittest.TestCase):
    def test_data_keeps_the_producer_weight(self) -> None:
        result = complete_event_weights(MANIFEST, "auau_data", event_weight=np.array([1.0, 2.5]), centrality=np.array([5.0, 95.0]))
        np.testing.assert_allclose(result["weight"], [1.0, 2.5])
        self.assertEqual(result["components"], ("producer_event_weight",))

    def test_ownership_window_and_centrality_in_order(self) -> None:
        sample = MANIFEST.samples["auau_jet20"]
        result = complete_event_weights(
            MANIFEST, "auau_jet20",
            event_weight=np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
            centrality=np.array([7.0, 7.0, 7.0, 85.0, 7.0]),
            leading_truth_jet_pt_gev=np.array([25.0, 31.0, 20.999, 25.0, math.nan]),
        )
        factor = MANIFEST.centrality.values["inclusivejet"][1]     # 5-10 percent
        expected = [sample.factor_value * factor, math.nan, math.nan, math.nan, math.nan]
        np.testing.assert_allclose(result["weight"], expected, rtol=1e-12)
        self.assertEqual(result["dropped"], {"outside_ownership_window": 3, "outside_centrality_support": 1})
        self.assertEqual(result["components"][0], "producer_event_weight")
        self.assertTrue(result["components"][1].startswith("source_stitch:auau_jet20"))
        self.assertTrue(result["components"][2].startswith("centrality:inclusivejet"))

    def test_photon_sample_uses_the_photonjet_map(self) -> None:
        result = complete_event_weights(MANIFEST, "auau_photon12", event_weight=np.array([2.0]), centrality=np.array([79.9]))
        expected = 2.0 * MANIFEST.samples["auau_photon12"].factor_value * MANIFEST.centrality.values["photonjet"][15]
        np.testing.assert_allclose(result["weight"], [expected], rtol=1e-12)

    def test_bad_inputs_are_errors_not_clamps(self) -> None:
        with self.assertRaises(ValueError):
            complete_event_weights(MANIFEST, "auau_photon12", event_weight=np.array([1.0]), centrality=np.array([math.nan]))
        with self.assertRaises(ValueError):
            complete_event_weights(MANIFEST, "auau_photon12", event_weight=np.array([-1.0]), centrality=np.array([5.0]))
        with self.assertRaises(ValueError):
            complete_event_weights(MANIFEST, "auau_jet12", event_weight=np.array([1.0]), centrality=np.array([5.0]))
        with self.assertRaises(KeyError):
            complete_event_weights(MANIFEST, "no_such_sample", event_weight=np.array([1.0]), centrality=np.array([5.0]))

    def test_pp_without_generated_events_applies_event_weight_only(self) -> None:
        result = complete_event_weights(MANIFEST, "pp_photonjet", event_weight=np.array([0.3]), centrality=np.array([-1.0]))
        np.testing.assert_allclose(result["weight"], [0.3])
        self.assertIn("not_applied", result["components"][1])


class LeadingTruthJetTest(unittest.TestCase):
    def test_radius_filter_and_maximum(self) -> None:
        truth_jets = {
            "source_file_index": np.array([0, 0, 0, 1]), "event_id_hi": np.array([1, 1, 1, 2]), "event_id_lo": np.array([1, 1, 1, 2]),
            "truth_jet_radius": np.array([0.4, 0.4, 0.2, 0.4]), "truth_jet_pt": np.array([12.0, 30.0, 50.0, math.nan]),
        }
        self.assertEqual(leading_truth_jet_pt(truth_jets), {(0, 1, 1): 30.0})


if __name__ == "__main__":
    unittest.main()
