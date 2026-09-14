from __future__ import annotations

import unittest

import numpy as np

from photonjet.analysis.response import (
    CLASSIFICATION_RECO_PTGAMMA_EDGES,
    CLASSIFICATION_TRUTH_PTGAMMA_EDGES,
    COMMON_XJ_EDGES,
    MEASURED_PTGAMMA_EDGES,
    Category,
    GlobalBinMap,
    classify_state,
)
from photonjet.analysis.unfolding import iterative_bayes, scan_problem


class ResponseUnfoldingTest(unittest.TestCase):
    def test_axes_and_global_index_are_exact(self) -> None:
        self.assertEqual(MEASURED_PTGAMMA_EDGES.tolist(), [15, 20, 25, 35])
        self.assertEqual(
            CLASSIFICATION_RECO_PTGAMMA_EDGES.tolist(),
            [10, 15, 20, 25, 35, 40],
        )
        self.assertEqual(
            CLASSIFICATION_TRUTH_PTGAMMA_EDGES.tolist(),
            [5, 10, 15, 20, 25, 35, 40],
        )
        mapping = GlobalBinMap(MEASURED_PTGAMMA_EDGES, COMMON_XJ_EDGES)
        self.assertEqual(mapping.n_global, 63)
        for index in range(mapping.n_global):
            self.assertEqual(mapping.flatten(*mapping.invert(index)), index)

    def test_response_boundary_categories(self) -> None:
        cases = (
            (14, 20, Category.PHOTON_LOW_FEED_IN),
            (35, 20, Category.PHOTON_HIGH_FEED_IN),
            (20, 14, Category.PHOTON_LOW_FEED_OUT),
            (20, 35, Category.PHOTON_HIGH_FEED_OUT),
            (20, None, Category.PHOTON_RECO_MISS),
            (None, 20, Category.UNMATCHED_RECO),
        )
        for truth, reco, expected in cases:
            self.assertEqual(
                classify_state(truth_ptgamma=truth, reco_ptgamma=reco)[0],
                expected,
            )

    def test_iterative_bayes_is_reproducible(self) -> None:
        matrix = np.asarray([[80.0, 10.0], [5.0, 70.0]])
        misses = np.asarray([10.0, 25.0])
        measured = np.asarray([90.0, 85.0])
        first = iterative_bayes(measured, matrix, misses, 4)
        second = iterative_bayes(measured, matrix, misses, 4)
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[1], second[1])

    def test_iteration_selection_requires_response_covariance_toys(self) -> None:
        with self.assertRaisesRegex(ValueError, "response/leakage/Tcomb toys"):
            scan_problem(
                None,
                np.asarray([1.0]),
                np.asarray([1.0]),
                np.asarray([[1.0], [1.0]]),
                response_toys=None,
            )


if __name__ == "__main__":
    unittest.main()
