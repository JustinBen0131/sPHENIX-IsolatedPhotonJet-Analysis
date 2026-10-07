"""A borderline line-fit error must preserve cuts and produce review evidence."""
import importlib.util
import json
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
from src import Diagnostics, Training, Validation

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


class RetainedCompletionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.folder = pathlib.Path(self.tmp.name); (self.folder / 'artifacts').mkdir()
        self.population, _, self.score, _, self.recipe = WorkingPointQualityTests().population(curvature=.001)
        self.population.y = np.arange(len(self.score)) % 2
        self.population.train = np.arange(len(self.score)) % 4 < 2
        self.population.receipt = dict(input_recipe_sha256='prepared', split_sha256='split', split_reference=None)
        self.cell = list(Training.cells(self.recipe))[0]
        self.model = self.folder / 'artifacts' / (self.cell['id'] + '.json')
        Models.save_json(self.model, dict(native_export_validation=dict(status='PASS')))
        for suffix in ('.xgboost.json', '.tmva.root'):
            self.model.with_suffix(suffix).write_bytes(b'retained original model')
        Models.save_json(self.model.with_suffix('.history.json'), dict(rounds=450))
        self.fit_code = dict(Models.code_identity(), **{'WorkingPoints.py': 'a' * 64, 'photon_id.py': 'b' * 64})
        self.original = dict(recipe=self.recipe, inputs=[], environment=dict(runtime='known'),
                             code_sha256=self.fit_code, split_reference=None)
        self.save_original()
        self.diagnosis = dict(schema='PhotonIDRetainedFitDiagnosisV1', status='REPORT_WRITTEN', model_fits=0,
            scientific_acceptance=False, request_sha256=self.original['request_sha256'], fit_code_sha256=self.fit_code,
            population_sha256='prepared', cells=[dict(cell=self.cell, validation=dict(status='PASS'),
            artifacts={str(p.relative_to(self.folder)): Models.file_hash(p) for p in (self.folder / 'artifacts').iterdir()})])
        Models.save_json(self.folder / 'validation_diagnosis.json', self.diagnosis)
        self.source_manifest = self.folder / 'fit_CODE_MANIFEST.json'
        self.save_fit_manifest()
        self.before = {p.name: p.read_bytes() for p in (self.folder / 'artifacts').iterdir()}
        self.request_before = (self.folder / 'request.json').read_bytes()
        self.args = types.SimpleNamespace(command='finalizeCampaign', check=False, progress_identity=None,
                                        fit_code_manifest=self.source_manifest, assemble_only=False)
        self.plan = dict(model=self.recipe['name'], split_reference=None, release=str(self.folder))

    def save_original(self):
        self.original['request_sha256'] = Models.digest(dict(recipe=self.recipe, inputs=[], code=self.fit_code,
                                                             environment=self.original['environment']))
        (self.folder / 'request.json').write_text(json.dumps(self.original, indent=2) + '\n')

    def save_fit_manifest(self):
        self.source_manifest.write_text(json.dumps(dict(schema='PhotonIDPublicCodeSnapshotV1', members=[
            dict(path='PhotonID/scripts/photon_id.py' if name == 'photon_id.py' else 'PhotonID/src/' + name, sha256=h)
            for name, h in self.fit_code.items()])) + '\n')

    def complete(self):
        def save_qa(folder, *args, **kwargs):
            Models.save_json(folder / 'qa.json', dict(status='numerical test fixture'))
        with mock.patch.object(Models, 'environment', return_value=dict(runtime='known')), \
             mock.patch.object(CLI, 'release_request', return_value=(self.recipe, self.folder, 'current-request',
                               dict(status='PASS', missing=[]), self.plan, [])), \
             mock.patch.object(CLI, 'prepared_population', return_value=(self.population, None)) as prepared, \
             mock.patch.object(Training, 'fit', side_effect=AssertionError('completion must never fit')) as fit, \
             mock.patch.object(Training, 'evaluator', return_value=lambda _: self.score), \
             mock.patch.object(Validation, 'check', return_value=(self.score, dict(status='PASS'))), \
             mock.patch.object(Diagnostics, 'save', side_effect=save_qa) as qa:
            result = CLI.train_release([], {}, 'PP', self.recipe['name'], self.args, simulation=[])
            self.assertEqual(fit.call_count, 0)
            return result, prepared.call_count, qa.call_count

    def test_completion_preserves_original_fit_and_request_and_seals_review_evidence(self):
        result, prepared, qa = self.complete()
        self.assertEqual((prepared, qa), (1, 1))
        self.assertEqual(result['calibration_quality']['status'], 'REVIEW_REQUIRED')
        manifest = Models.load_release(self.folder, request='current-request')
        self.assertEqual(manifest['retained_fit']['fit_request_sha256'], self.original['request_sha256'])
        self.assertFalse(manifest['calibration_quality']['scientific_acceptance'])
        self.assertEqual((self.folder / 'request.json').read_bytes(), self.request_before)
        self.assertEqual(self.before, {p.name: p.read_bytes() for p in (self.folder / 'artifacts').iterdir()
                                      if p.name in self.before})
        self.assertTrue((self.folder / 'finalization_request.json').is_file())
        self.assertEqual(self.complete()[0]['status'], 'REUSED')

    def test_missing_diagnosis_never_replaces_original_request_or_starts_training(self):
        (self.folder / 'validation_diagnosis.json').unlink()
        with self.assertRaises(FileNotFoundError): self.complete()
        self.assertEqual((self.folder / 'request.json').read_bytes(), self.request_before)
        self.assertFalse((self.folder / 'manifest.json').exists())

    def test_changed_fitter_is_refused_before_completion(self):
        self.fit_code['Training.py'] = 'c' * 64
        self.save_original(); self.save_fit_manifest()
        with self.assertRaisesRegex(ValueError, 'unchanged fitter'): self.complete()
        self.assertFalse((self.folder / 'finalization_request.json').exists())

    def test_changed_fitted_artifact_is_refused(self):
        self.model.with_suffix('.xgboost.json').write_bytes(b'changed model')
        with self.assertRaisesRegex(ValueError, 'changed after diagnosis'): self.complete()
        self.assertFalse((self.folder / 'manifest.json').exists())

    def test_different_population_cannot_be_finalized(self):
        self.population.receipt['input_recipe_sha256'] = 'different'
        with self.assertRaisesRegex(ValueError, 'another prepared population'): self.complete()
        self.assertFalse((self.folder / 'manifest.json').exists())

    def test_original_runtime_must_match(self):
        self.original['environment'] = dict(runtime='different'); self.save_original()
        with self.assertRaisesRegex(ValueError, 'another recipe, source population or runtime'): self.complete()

    def test_incomplete_fit_cannot_be_completed(self):
        path = self.model.with_suffix('.history.json'); path.write_text(json.dumps(dict(rounds=449)))
        self.diagnosis['cells'][0]['artifacts'][str(path.relative_to(self.folder))] = Models.file_hash(path)
        (self.folder / 'validation_diagnosis.json').write_text(json.dumps(self.diagnosis))
        with self.assertRaisesRegex(ValueError, 'complete training'): self.complete()

    def test_frontend_requires_one_original_model_and_catalogue(self):
        args = CLI.arguments(['finalizeCampaign', 'PP', '--model', 'nominal_retrained_v1', '--inputs', 'inputs.yaml'])
        self.assertTrue(args.no_score)
        with self.assertRaises(SystemExit): CLI.arguments(['finalizeCampaign', 'PP'])


if __name__ == '__main__':
    unittest.main()
