"""A borderline line-fit error must preserve cuts and produce review evidence."""
import importlib.util
import pathlib
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from src import Models, WorkingPoints
from src import Training, Validation

spec = importlib.util.spec_from_file_location('photonid_quality_cli', HERE / 'scripts/photon_id.py')
CLI = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CLI)


class WorkingPointQualityTests(unittest.TestCase):
    def population(self, curvature=0., low=.3):
        recipe = Models.recipe(Models.load_config(HERE / 'config/training.yaml'), 'PP', 'nominal_retrained_v1')
        edges = recipe['working_points']['et_edges']
        n = 1000
        et = np.concatenate([np.linspace(a + 1e-4, b - 1e-4, n) for a, b in zip(edges[:-1], edges[1:])])
        score = np.concatenate([low + curvature * (i - 3.5) ** 2 + np.linspace(0., .25, n)
                                for i in range(len(edges) - 1)])
        population = types.SimpleNamespace(y=np.ones(len(et), np.uint8), et=et, weights=np.ones(len(et)))
        return population, np.ones(len(et), bool), score, dict(bounds={}), recipe

    def test_passing_calibration_is_identical_in_both_modes(self):
        args = self.population()
        self.assertEqual(WorkingPoints.derive(*args), WorkingPoints.derive(*args, strict_fit_checks=False))

    def test_borderline_error_records_all_working_points_without_changing_line(self):
        args = self.population(curvature=.001)
        with self.assertRaises(WorkingPoints.CalibrationError) as caught:
            WorkingPoints.derive(*args)
        original = caught.exception.report
        self.assertGreater(original['max_efficiency_error'], .02)
        self.assertLess(original['max_efficiency_error'], .04)
        result = WorkingPoints.derive(*args, strict_fit_checks=False)
        self.assertEqual(result['fit_quality_status'], 'REVIEW_REQUIRED')
        self.assertEqual(set(result['cuts']), {'WP70', 'WP80', 'WP90'})
        self.assertEqual(result['cuts']['WP70']['intercept'], original['intercept'])
        self.assertEqual(result['cuts']['WP70']['slope'], original['slope'])
        self.assertEqual(result['cuts']['WP70']['bins'], original['bins'])
        self.assertTrue(result['fit_warnings'])

    def test_advisory_precision_does_not_allow_empty_non_tight_region(self):
        args = self.population(low=0.)
        with self.assertRaisesRegex(ValueError, 'overlap or are empty'):
            WorkingPoints.derive(*args, strict_fit_checks=False)

    def test_training_records_warning_and_reuses_cell_without_another_fit(self):
        population, _, score, cell, recipe = self.population(curvature=.001)
        population.y = np.arange(len(score)) % 2
        population.train = np.arange(len(score)) % 4 < 2
        cell['id'] = 'cell_000'
        artifact = dict(native_export_validation=dict(status='PASS'))
        def fit(pop, mask, r, path, cell_id, **kwargs):
            Models.save_json(path, artifact)
            for suffix in ('.xgboost.json', '.tmva.root', '.history.json'):
                path.with_suffix(suffix).write_bytes(b'unchanged fitted artifact')
            return artifact, lambda _: score
        with tempfile.TemporaryDirectory() as tmp, \
             mock.patch.object(Training, 'fit', side_effect=fit) as fitter, \
             mock.patch.object(Training, 'evaluator', return_value=lambda _: score), \
             mock.patch.object(Validation, 'check', return_value=(score, dict(status='PASS'))):
            folder = pathlib.Path(tmp)
            first = CLI.validated_cell(population, recipe, cell, folder, 'request', None)
            self.assertEqual(first[2]['fit_quality_status'], 'REVIEW_REQUIRED')
            receipt = folder / 'artifacts/cell_000.validation.json'
            before = receipt.read_bytes()
            second = CLI.validated_cell(population, recipe, cell, folder, 'request', None, allow_fit=False)
            self.assertEqual(first[2], second[2])
            self.assertEqual(receipt.read_bytes(), before)
            self.assertEqual(fitter.call_count, 1)
            self.assertEqual(fitter.call_args.args[2], recipe)


if __name__ == '__main__':
    unittest.main()
