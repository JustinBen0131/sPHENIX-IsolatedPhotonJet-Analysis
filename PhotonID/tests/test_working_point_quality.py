"""Quality warnings preserve fitted models; numerical failures still stop release."""
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
from src import CellInputs, Diagnostics, Training, Validation

spec = importlib.util.spec_from_file_location('photonid_quality_cli', HERE / 'scripts/photon_id.py')
CLI = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CLI)


class PerformanceQualityTests(unittest.TestCase):
    def fixture(self):
        rows = 80
        y = np.arange(rows) % 2
        score = y.astype(float)
        score[40:50] = 1. - score[40:50]
        population = types.SimpleNamespace(x=np.arange(rows)[:, None], y=y, weights=np.ones(rows),
            train=np.arange(rows) < 40, holdout=np.arange(rows) >= 40, partition=None,
            receipt=dict(split_sha256='fixed-split'))
        recipe = dict(validation=dict(saved_score_atol=1e-8, min_rows_per_class_per_split=10,
                      min_holdout_auc=.55, max_train_holdout_auc_gap=.05))
        evaluator = lambda x: score[x[:, 0]]
        return population, np.ones(rows, bool), evaluator, evaluator, recipe

    def test_gap_warning_keeps_scores_metrics_and_recipe_limits(self):
        args = self.fixture()
        with self.assertRaisesRegex(ValueError, 'discrimination gap failed'):
            Validation.check(*args)
        score, report = Validation.check(*args, strict_performance=False)
        self.assertEqual(report['status'], 'PASS')
        self.assertEqual(report['train']['auc'], 1.)
        self.assertEqual(report['holdout']['auc'], .75)
        quality = report['performance_quality']
        self.assertEqual(quality['status'], 'REVIEW_REQUIRED')
        self.assertFalse(quality['scientific_acceptance'])
        self.assertEqual(quality['warnings'], [dict(criterion='max_train_holdout_auc_gap', measured=.25, limit=.05)])
        np.testing.assert_array_equal(score, args[2](args[0].x))

    def test_advisory_mode_still_rejects_changed_saved_scores(self):
        population, mask, evaluator, _, recipe = self.fixture()
        with self.assertRaisesRegex(ValueError, 'differs from trainer'):
            Validation.check(population, mask, evaluator, lambda x: evaluator(x) + .01,
                             recipe, strict_performance=False)

    def test_advisory_mode_still_rejects_nonfinite_saved_scores(self):
        population, mask, evaluator, _, recipe = self.fixture()
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            Validation.check(population, mask, evaluator, lambda x: np.full(len(x), np.nan),
                             recipe, strict_performance=False)


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

    def complete(self, validation_report=None):
        def save_qa(folder, *args, **kwargs):
            Models.save_json(folder / 'qa.json', dict(status='numerical test fixture'))
        with mock.patch.object(Models, 'environment', return_value=dict(runtime='known')), \
             mock.patch.object(CLI, 'release_request', return_value=(self.recipe, self.folder, 'current-request',
                               dict(status='PASS', missing=[]), self.plan, [])), \
             mock.patch.object(CLI, 'prepared_population', return_value=(self.population, None)) as prepared, \
             mock.patch.object(Training, 'fit', side_effect=AssertionError('completion must never fit')) as fit, \
             mock.patch.object(Training, 'evaluator', return_value=lambda _: self.score), \
             mock.patch.object(Validation, 'check', return_value=(self.score, validation_report or dict(status='PASS'))), \
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

    def test_compact_completion_uses_original_fit_request(self):
        (self.folder / '.work').mkdir()
        (self.folder / '.work/cell_inputs.json').write_text('{}')
        self.args.models = self.folder.parent
        with mock.patch.object(CellInputs, 'load', return_value=(self.population,
                np.arange(len(self.score)), None, None)) as loader, \
             mock.patch.object(CellInputs, 'verify_origin'), \
             mock.patch.object(CellInputs, 'population_witness', return_value={}):
            self.assertEqual(self.complete()[0]['status'], 'PASS')
            self.assertEqual(loader.call_args.args[1], self.original['request_sha256'])
            self.assertNotEqual(loader.call_args.args[1], 'current-request')

    def test_completion_preserves_performance_warning_in_manifest(self):
        warning = dict(criterion='max_train_holdout_auc_gap', measured=.06, limit=.05)
        report = dict(status='PASS', performance_quality=dict(status='REVIEW_REQUIRED',
                      warnings=[warning], scientific_acceptance=False))
        self.diagnosis['cells'][0]['validation'] = report
        (self.folder / 'validation_diagnosis.json').write_text(json.dumps(self.diagnosis) + '\n')
        self.assertEqual(self.complete(validation_report=report)[0]['performance_quality']['status'], 'REVIEW_REQUIRED')
        manifest = Models.read_json(self.folder / 'manifest.json')
        self.assertEqual(manifest['performance_quality']['warnings'], [dict(cell=self.cell['id'], **warning)])
        self.assertFalse(manifest['performance_quality']['scientific_acceptance'])
        self.assertEqual(self.before, {p.name: p.read_bytes() for p in (self.folder / 'artifacts').iterdir()
                                      if p.name in self.before})

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

    def test_preparation_update_requires_explicit_original_population(self):
        self.fit_code.update({'Inputs.py': 'd' * 64, 'Models.py': 'e' * 64})
        self.save_original(); self.save_fit_manifest()
        self.diagnosis['request_sha256'] = self.original['request_sha256']
        (self.folder / 'validation_diagnosis.json').write_text(json.dumps(self.diagnosis) + '\n')
        with self.assertRaisesRegex(ValueError, 'unchanged fitter'): self.complete()
        self.args.prepared_population = self.folder / 'original-population'
        self.assertEqual(self.complete()[0]['status'], 'PASS')

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


class RetainedPopulationTests(unittest.TestCase):
    def setUp(self):
        from src import Inputs
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.work = pathlib.Path(self.tmp.name).resolve()
        self.recipe = WorkingPointQualityTests().population(curvature=.001)[4]
        self.sources = [dict(logical_id='sealed-source', entries=20)]
        n = 20
        arrays = dict(x=np.zeros((n, len(self.recipe['features'])), np.float32),
            y=np.arange(n, dtype=np.int8) % 2, et=np.full(n, 20.), centrality=np.zeros(n),
            weights=np.ones(n), train=np.arange(n) % 10 != 0, partition=(np.arange(n) % 10 == 0).astype(np.uint8) * 2,
            split_rank=np.arange(n, dtype=np.uint64), source_index=np.zeros(n, np.uint32),
            entry=np.arange(n), photon=np.zeros(n, np.int32), event_key=np.arange(n).astype('S32'))
        for name, values in arrays.items(): np.save(self.work / (name + '.npy'), values)
        definition = Inputs.population_definition(self.recipe)
        original_code = dict(Models.code_identity(), **{'Inputs.py': 'a' * 64, 'Models.py': 'b' * 64})
        self.manifest = self.work / 'fit_CODE_MANIFEST.json'
        Models.save_json(self.manifest, dict(schema='PhotonIDPublicCodeSnapshotV1', members=[
            dict(path='PhotonID/scripts/photon_id.py' if key == 'photon_id.py' else 'PhotonID/src/' + key, sha256=value)
            for key, value in original_code.items()]))
        self.receipt = dict(schema='PhotonIDPopulationV1', rows=n, sources=self.sources, definition=definition,
            split_reference=None, array_files={key:key + '.npy' for key in arrays},
            arrays={key + '.npy':Models.file_hash(self.work / (key + '.npy')) for key in arrays},
            input_recipe_sha256=Models.digest(dict(inputs=self.sources, **definition, preparation_code=original_code['Inputs.py'])),
            event_split_sha256='legacy retired event-store witness')
        self.save_receipt()
        self.args = types.SimpleNamespace(command='finalizeCampaign', fit_code_manifest=self.manifest,
                                         prepared_population=self.work)

    def save_receipt(self):
        (self.work / 'population.json').write_text(json.dumps(self.receipt, indent=2) + '\n')

    def test_reuse_opens_original_arrays_without_database_or_new_files(self):
        before = {p.name:p.read_bytes() for p in self.work.iterdir()}
        population, work = CLI.retained_population(self.sources, self.recipe, self.args)
        self.assertEqual(work, self.work)
        self.assertEqual(population.receipt, self.receipt)
        self.assertEqual(population.x.shape, (20, len(self.recipe['features'])))
        self.assertEqual(before, {p.name:p.read_bytes() for p in self.work.iterdir()})
        self.assertEqual(self.args.population_reuse['arrays'], self.receipt['arrays'])

    def test_changed_array_is_refused(self):
        np.save(self.work / 'weights.npy', np.full(20, 2.))
        with self.assertRaisesRegex(ValueError, 'array changed'): CLI.retained_population(self.sources, self.recipe, self.args)

    def test_wrong_sources_recipe_or_original_code_are_refused(self):
        with self.assertRaisesRegex(ValueError, 'original sources'): CLI.retained_population([], self.recipe, self.args)
        self.receipt['input_recipe_sha256'] = 'wrong'; self.save_receipt()
        with self.assertRaisesRegex(ValueError, 'original preparation code'): CLI.retained_population(self.sources, self.recipe, self.args)

    def test_missing_array_and_borrowed_split_are_refused(self):
        original = dict(self.receipt['array_files'])
        self.receipt['array_files'].pop('event_key'); self.save_receipt()
        with self.assertRaisesRegex(ValueError, 'array roster'): CLI.retained_population(self.sources, self.recipe, self.args)
        self.receipt['array_files'] = original; self.receipt['split_reference'] = dict(input_recipe_sha256='other'); self.save_receipt()
        with self.assertRaisesRegex(ValueError, 'split reference'): CLI.retained_population(self.sources, self.recipe, self.args)

    def test_recovery_option_is_completion_only(self):
        with self.assertRaises(SystemExit): CLI.arguments(['trainCampaign', 'PP', '--prepared-population', str(self.work)])
        with self.assertRaises(SystemExit): CLI.arguments(['finalizeCampaign', 'PP', '--model', 'nominal_retrained_v1',
                                                        '--inputs', 'inputs.yaml', '--prepared-population', str(self.work)])


if __name__ == '__main__':
    unittest.main()
