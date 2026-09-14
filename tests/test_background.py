from __future__ import annotations

from dataclasses import replace
import unittest

import numpy as np

from photonjet.analysis import background


def correction_fixture(n_xj: int = 3):
    """Synthetic factorized background with known 80-photon signal per bin."""
    leakage = np.tile(np.asarray([1.0, 0.10, 0.05, 0.02])[:, None], (1, 3))
    counts = np.tile(np.asarray([20.0, 40.0, 10.0, 20.0])[:, None], (1, 3)) + 80 * leakage
    spectra = np.zeros((2, 3, n_xj))
    spectra[0] = np.linspace(10.0, 25.0, n_xj)
    spectra[1] = np.linspace(1.0, 8.0, n_xj)
    data = {
        "abcd_event_leading_counts_ptgamma": counts,
        "abcd_event_leading_sumw2_ptgamma": counts.copy(),
        "inclusive_recoil_spectra_ptgamma_xj": spectra,
        "inclusive_recoil_sumw2_ptgamma_xj": spectra.copy(),
    }
    response = background.CorrectionInputs(
        leakage=leakage,
        leakage_sumw2=leakage.copy(),
        combinatoric_reco=np.full((3, n_xj), 5.0),
        combinatoric_reco_sumw2=np.full((3, n_xj), 5.0),
        combinatoric_normalization_reco=np.full(3, 50.0),
        combinatoric_normalization_reco_sumw2=np.full(3, 50.0),
    )
    return data, response


class BackgroundTest(unittest.TestCase):
    def test_analytic_solution_and_region_c_restoration(self):
        data, response = correction_fixture()
        out = background.correct_measured(data, response, False)
        np.testing.assert_allclose(out["photons"], [80.0] * 3, rtol=0, atol=1e-12)
        a, c = data["inclusive_recoil_spectra_ptgamma_xj"]
        np.testing.assert_allclose(out["corrected"], (a - 2 * c) / 0.9)
        self.assertTrue(all(row["method"] == "analytic_physical_root" for row in out["solutions"]))

    def test_raw_abcd_fallback_remains_visible(self):
        n = np.asarray([1930., 320., 330., 89.])
        result = background.solve_leakage_abcd(n, [1., .24823406, .10408785, .02072144])
        self.assertEqual(result.method, "raw_abcd_fallback_no_physical_leakage_root")
        self.assertEqual(result.signal_a, n[0] - n[1] * n[2] / n[3])

    def test_truth_conditioned_normalization_and_signed_recoil(self):
        data, response = correction_fixture()
        uncorrected = background.correct_measured(data, response, False)
        corrected = background.correct_measured(data, response, True)
        np.testing.assert_allclose(corrected["combinatoric_subtracted"], np.full((3, 3), 8.))
        np.testing.assert_allclose(corrected["corrected"], uncorrected["corrected"] - 8.)
        response = replace(response, combinatoric_reco=np.full((3, 3), 50.))
        signed = background.correct_measured(data, response, True)
        self.assertTrue(np.any(signed["corrected"] < 0))
        self.assertGreater(signed["negative_input_sum"], 0.)

    def test_template_and_normalization_variance_are_preserved(self):
        data, response = correction_fixture()
        out = background.correct_measured(data, response, True)
        expected = (80 / 50) ** 2 * 5 + (-80 * 5 / 50 ** 2) ** 2 * 50
        np.testing.assert_allclose(out["combinatoric_variance"], np.full((3, 3), expected))

    def test_pooled_transfer_and_integrated_view_are_distinct(self):
        data, response = correction_fixture()
        pooled = background.correct_measured(data, response, True, purity_strategy="integrated_15_35_transfer")
        integrated = background.correct_integrated_measured(data, response)
        self.assertEqual(pooled["corrected"].shape, (3, 3))
        self.assertEqual(integrated["corrected"].shape, (3,))
        self.assertTrue(all(x["method"] == "integrated_15_35_leakage_transfer" for x in pooled["solutions"]))
        self.assertAlmostEqual(integrated["photons"], 240.)

    def test_explicit_replica_inputs_are_used_together(self):
        data, response = correction_fixture()
        first = background.correct_measured(data, response, True)
        paired = background.correct_measured(
            data, response, True, combinatoric_reco=2 * response.combinatoric_reco,
            combinatoric_normalization_reco=2 * response.combinatoric_normalization_reco,
        )
        np.testing.assert_array_equal(first["corrected"], paired["corrected"])

    def test_invalid_strategy_shapes_and_leakage_fail(self):
        data, response = correction_fixture()
        cases = [
            {"purity_strategy": "from_plot_label"},
            {"combinatoric_reco": np.zeros((2, 3))},
            {"combinatoric_normalization_reco": np.zeros(2)},
        ]
        for kwargs in cases:
            with self.subTest(kwargs=list(kwargs)), self.assertRaises(ValueError):
                background.correct_measured(data, response, True, **kwargs)
        with self.assertRaises(ValueError):
            background.correct_measured(data, replace(response, leakage=np.zeros((4, 3))), True)

    def test_toy_helpers_are_seeded_and_preserve_sumw2(self):
        data, response = correction_fixture()
        first = background.draw_joint_data_toy(data, np.random.default_rng(817))
        second = background.draw_joint_data_toy(data, np.random.default_rng(817))
        for key in first:
            np.testing.assert_array_equal(first[key], second[key])
        np.testing.assert_array_equal(first["inclusive_recoil_sumw2_ptgamma_xj"], data["inclusive_recoil_sumw2_ptgamma_xj"])
        for a, b in zip(
            background.draw_combinatoric_template_toy(response, np.random.default_rng(911)),
            background.draw_combinatoric_template_toy(response, np.random.default_rng(911)),
        ):
            np.testing.assert_array_equal(a, b)
            self.assertTrue(np.all(a >= 0))


if __name__ == "__main__":
    unittest.main()
