"""Boundary and contract tests for PhotonID/photon_selection.py."""

from __future__ import annotations

import unittest

import numpy as np

from photon_selection import (
    classify,
    event_leading_index,
    region_masks,
    threshold_values,
    truth_signal_mask,
)


TRUTH_CONTRACT = {
    "require_geant_photon": True,
    "require_generator_association": True,
    "valid_prompt_classes": [-1, 0, 1, 2],
    "truth_isolation_max_gev": 4.0,
    "abs_eta_max": 0.7,
}


def thresholds(n: int, tight=0.8, low=0.4, high=0.7, iso_max=4.0, noniso_min=7.0) -> dict[str, np.ndarray]:
    return {
        "tight": np.full(n, tight), "nontight_low": np.full(n, low), "nontight_high": np.full(n, high),
        "isolated_max": np.full(n, iso_max), "nonisolated_min": np.full(n, noniso_min),
    }


class FlagBoundaryTest(unittest.TestCase):
    def test_tight_and_bounded_non_tight_boundaries(self) -> None:
        score = np.array([0.8, 0.81, 0.4, 0.41, 0.7, 0.75, np.nan])
        flags = classify(score, np.full(7, 1.0), thresholds(7))
        self.assertEqual(flags["bdt_is_tight"].tolist(), [0, 1, 0, 0, 0, 0, 0])      # score == tight is not tight
        self.assertEqual(flags["bdt_is_nontight"].tolist(), [0, 0, 0, 1, 1, 0, 0])   # (low, high]; 0.75 is in the gap
        self.assertEqual(flags["bdt_is_not_tight"].tolist(), [1, 0, 1, 1, 1, 1, 0])
        self.assertEqual(flags["id_valid"].tolist(), [True] * 6 + [False])

    def test_isolation_boundaries_leave_a_gap(self) -> None:
        iso = np.array([3.9, 4.0, 5.0, 7.0, 7.1])
        flags = classify(np.full(5, 0.9), iso, thresholds(5))
        self.assertEqual(flags["iso_pass"].tolist(), [1, 0, 0, 0, 0])
        self.assertEqual(flags["iso_nonisolated"].tolist(), [0, 0, 0, 0, 1])

    def test_missing_thresholds_make_candidates_invalid_not_background(self) -> None:
        t = thresholds(2)
        t["tight"] = np.full(2, np.nan)
        flags = classify(np.array([0.9, 0.1]), np.array([1.0, 1.0]), t)
        self.assertEqual(flags["id_valid"].tolist(), [False, False])
        self.assertEqual(flags["bdt_is_not_tight"].tolist(), [0, 0])
        regions = region_masks(flags, "complement")
        self.assertFalse(regions["C"].any())

    def test_regions_partition_by_flags(self) -> None:
        score = np.array([0.9, 0.9, 0.5, 0.5])
        iso = np.array([1.0, 8.0, 1.0, 8.0])
        regions = region_masks(classify(score, iso, thresholds(4)), "bounded")
        self.assertEqual([regions[r].tolist() for r in "ABCD"],
                         [[True, False, False, False], [False, True, False, False],
                          [False, False, True, False], [False, False, False, True]])


class ThresholdSpecificationTest(unittest.TestCase):
    def test_constant_linear_and_binned(self) -> None:
        et = np.array([16.0, 22.0, 30.0, 40.0])
        cent = np.array([5.0, 15.0, 30.0, 60.0])
        np.testing.assert_allclose(threshold_values(0.5, et, cent), [0.5] * 4)
        np.testing.assert_allclose(threshold_values({"intercept": 1.0, "slope_photon_et": 0.1}, et, cent), 1.0 + 0.1 * et)
        binned = {"mode": "binned", "pt_edges": [15, 20, 25, 35], "values": [0.1, 0.2, 0.3]}
        out = threshold_values(binned, et, cent)
        np.testing.assert_allclose(out[:3], [0.1, 0.2, 0.3])
        self.assertTrue(np.isnan(out[3]))                       # outside the pT grid: not derived
        grid = {"mode": "binned", "pt_edges": [15, 20, 25, 35], "cent_edges": [0, 20, 50],
                "values": [[1, 2], [3, 4], [5, 6]]}
        np.testing.assert_allclose(threshold_values(grid, et, cent)[:3], [1, 3, 6])

    def test_null_threshold_is_nan(self) -> None:
        self.assertTrue(np.all(np.isnan(threshold_values(None, np.array([20.0]), np.array([0.0])))))


class TruthSignalContractTest(unittest.TestCase):
    def truth(self, **overrides):
        base = {
            "prompt_class": np.array([1, 2, -1, -1, 3, 0]),
            "g4_photon_valid": np.array([1, 1, 1, 1, 1, 1]),
            "hepmc_association_valid": np.array([1, 1, 1, 0, 1, 1]),
            "truth_isolation_valid": np.array([1, 1, 1, 1, 1, 1]),
            "truth_isolation_r03": np.array([1.0, 3.9, 2.0, 2.0, 1.0, 4.0]),
            "truth_photon_eta": np.array([0.1, -0.6, 0.0, 0.0, 0.2, 0.3]),
        }
        base.update(overrides)
        return base

    def test_contract(self) -> None:
        mask = truth_signal_mask(self.truth(), TRUTH_CONTRACT)
        # class 3 rejected; class -1 accepted only with a valid association; iso == 4 GeV rejected
        self.assertEqual(mask.tolist(), [True, True, True, False, False, False])

    def test_missing_fields_are_an_error(self) -> None:
        truth = self.truth()
        del truth["truth_isolation_r03"]
        with self.assertRaises(KeyError):
            truth_signal_mask(truth, TRUTH_CONTRACT)


class LeaderTest(unittest.TestCase):
    def test_highest_et_then_lowest_ordinal(self) -> None:
        et = np.array([20.0, 25.0, 25.0, 30.0])
        ordinal = np.array([3, 2, 1, 0])
        eligible = np.array([True, True, True, False])
        self.assertEqual(event_leading_index(et, ordinal, eligible), 2)
        self.assertEqual(event_leading_index(et, ordinal, np.zeros(4, dtype=bool)), -1)


if __name__ == "__main__":
    unittest.main()
