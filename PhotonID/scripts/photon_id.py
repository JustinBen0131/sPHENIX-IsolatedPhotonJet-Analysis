#!/usr/bin/env python3
"""Public PhotonID driver: project, verify caches, prepare, fit and score.

All algorithms share one eventTree input contract and saved split. Existing
validated releases are reused; an interrupted candidate is never promoted by
filename alone. Completed groups are published by TreeProduction's installed
eventTree library (MYINSTALL). The shell entry point selects an already
configured Python environment.
"""
from __future__ import annotations

import argparse
import datetime
import fcntl
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))


def arguments(argv=None):
    parser = argparse.ArgumentParser(description='Train and score named photon models from completed eventTrees.')
    parser.add_argument('command', choices=['reTrainNominal', 'trainNew', 'inspectCampaign', 'projectCampaign', 'verifyProjections',
                                         'prepareCampaign', 'trainCampaign', 'collectQA', 'score', 'listModels', 'planCampaign'])
    parser.add_argument('target', nargs='?', help='AuAu/PP/AuAuAndPP, model recipe or saved model name')
    parser.add_argument('system', nargs='?', choices=['AuAu', 'PP'])
    parser.add_argument('--config', type=Path, default=HERE / 'config/training.yaml')
    parser.add_argument('--models', type=Path, default=HERE / 'models', help='model release root')
    parser.add_argument('--group', type=Path, action='append', default=[], help='complete published group; repeat as needed')
    parser.add_argument('--code-root', type=Path, default=HERE.parent, help='repository containing TreeProduction and its eventTree schema')
    parser.add_argument('--release-version', default='v1', help='explicit version for a new training recipe/input population')
    parser.add_argument('--check', action='store_true', help='read-only input, recipe and dependency preview')
    parser.add_argument('--no-score', action='store_true', help='finish model release without publishing dataset augmentation')
    parser.add_argument('--progress-context', type=Path, help='YAML containing task_id, workstream_id, campaign_tag and executor')
    parser.add_argument('--inputs', type=Path, help='sealed training-input YAML from inspectCampaign')
    parser.add_argument('--require-roster', action='store_true', help='inspectCampaign must use a saved producer roster; never fall back to tree headers')
    parser.add_argument('--projection-root', type=Path, help='shared source-row cache; defaults to <models>/.projections')
    parser.add_argument('--shard-index', type=int, help='one projectCampaign shard from its --check plan')
    parser.add_argument('--max-events', type=int, default=100000, help='maximum events per projection shard (at most 100000)')
    parser.add_argument('--cache-layout', choices=('packed', 'source'), default='packed', help='projection storage; packed is one row bundle per shard')
    parser.add_argument('--scratch-root', type=Path, help='existing physical data-staging directory; otherwise use _CONDOR_SCRATCH_DIR')
    parser.add_argument('--auau-reference', type=Path, help='sealed canonical AuAu release for the pp/AuAu heldout comparison')
    parser.add_argument('--model', action='append', default=[], help='one trainCampaign recipe; repeat to select several, then run collectQA')
    args = parser.parse_args(argv)
    if args.command == 'reTrainNominal' and args.target not in ('AuAu', 'PP', 'AuAuAndPP'):
        parser.error('reTrainNominal requires AuAu, PP or AuAuAndPP')
    if args.command == 'trainNew' and not ((args.target in ('centAndETdepModels', 'routingStudy', 'rangeStudy', 'depthStudy', 'MLP', 'logReg') and args.system == 'AuAu') or (args.target == 'PPG12' and args.system == 'PP')):
        parser.error('trainNew requires an AuAu study or PPG12 PP')
    if args.command in ('planCampaign', 'inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign', 'trainCampaign', 'collectQA') and (args.target not in ('AuAu', 'PP') or args.system):
        parser.error(args.command + ' requires AuAu or PP')
    if args.command in ('inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign', 'trainCampaign') or args.inputs:
        args.no_score = True
    # Training produces portable releases. Publishing scored tree copies is a
    # separate explicit command; reduction normally evaluates models in place.
    if args.command in ('reTrainNominal', 'trainNew'): args.no_score = True
    if args.command == 'score' and (not args.target or not args.system):
        parser.error('score requires a saved model name and AuAu/PP')
    if args.command == 'listModels' and args.target not in (None, 'AuAu', 'PP'):
        parser.error('listModels accepts AuAu or PP')
    if args.inputs and (args.group or args.target == 'AuAuAndPP' or args.command in ('score', 'listModels', 'planCampaign')):
        parser.error('--inputs is a training catalogue; do not combine it with --group, paired systems or scoring')
    if args.command == 'verifyProjections' and not args.inputs:
        parser.error('verifyProjections requires --inputs from inspectCampaign')
    if args.command != 'projectCampaign' and (args.shard_index is not None or args.max_events != 100000):
        parser.error('--shard-index and --max-events apply only to projectCampaign')
    if args.command != 'projectCampaign' and (args.scratch_root is not None or args.cache_layout != 'packed'):
        parser.error('--scratch-root and --cache-layout apply only to projectCampaign')
    if args.scratch_root is not None and args.cache_layout == 'source':
        parser.error('--scratch-root requires packed projection storage')
    if not 0 < args.max_events <= 100000 or (args.shard_index is not None and args.shard_index < 0):
        parser.error('projection event bound must be 1..100000 and shard index must be nonnegative')
    if args.auau_reference and (args.command not in ('trainCampaign', 'collectQA') or args.target != 'PP'):
        parser.error('--auau-reference applies to trainCampaign PP or collectQA PP')
    if args.model and args.command != 'trainCampaign':
        parser.error('--model applies only to trainCampaign')
    if args.require_roster and args.command != 'inspectCampaign':
        parser.error('--require-roster applies only to inspectCampaign')
    return args


def training_baseline(inputs, *, verify_payload=True):
    """Pin the stage-2 ancestor so later score/fraction additions do not retrain."""
    from src.Models import digest, file_hash, read_json, require
    result = []
    cache = {}
    for item in inputs:
        group = Path(item['group']); version = item['version']
        if (str(group), version) not in cache:
            visited = set()
            while True:
                require(version not in visited and len(visited) < 256, 'invalid dataset version lineage')
                visited.add(version)
                directory = group / 'versions' / version
                m = read_json(directory / 'MANIFEST.json')
                require(digest({k: v for k, v in m.items() if k != 'version'}) == m['version'] == version,
                        'dataset ancestor manifest changed')
                if m['stage'] == 2:
                    break
                require(m['stage'] == 3 and m.get('base_version'), 'stage-2 training ancestor unavailable')
                version = m['base_version']
            cache[(str(group), item['version'])] = directory, m, {r['logical_id']: r for r in m['files']}
        directory, m, roster = cache[(str(group), item['version'])]
        require(item['logical_id'] in roster, 'augmentation lineage lost a logical source')
        baseline = roster[item['logical_id']]
        path = directory / baseline['relative_path']
        if verify_payload: require(file_hash(path) == baseline['sha256'], 'stage-2 training ancestor changed')
        result.append(dict(item, path=str(path), sha256=baseline['sha256'], version=m['version']))
    return result


def file_snapshot(paths):
    """Cheap in-process witnesses after full hashing, including replacement detection."""
    result = []
    for path in paths:
        path = Path(path).resolve(); stat = path.stat()
        result.append((str(path), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))
    return result


def prepared_population(simulation, recipe, config, args, cache):
    """Reuse verified matrices within one invocation; widen a frozen nominal split.

    Cached reuse checks source and matrix file metadata without rereading all
    ROOT bytes for every learner. A later invocation performs full hash checks.
    """
    from src import Inputs, Models
    projections = getattr(args, 'projections', None)
    if projections is None and getattr(args, 'use_projection_cache', False):
        projections = args.projections = campaign_projections(simulation, config, args, system=recipe['system'])
    reference = None; reference_pin = None
    if recipe['name'].startswith('wideET_'):
        nominal = Models.recipe(config, recipe['system'], 'nominal_retrained_v1')
        _, reference = prepared_population(simulation, nominal, config, args, cache)
        reference_pin = Inputs.population_identity(simulation, nominal)
    pin = Inputs.population_identity(simulation, recipe, reference_pin)
    if pin in cache:
        population, work, paths, snapshot = cache[pin]
        Models.require(file_snapshot(paths) == snapshot, 'cached source or prepared population changed')
        return population, work
    work = args.models.resolve() / '.populations' / pin
    work.mkdir(parents=True, exist_ok=True)
    with (work / '.preparation.lock').open('a') as preparation_lock:
        fcntl.flock(preparation_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        population = Inputs.prepare(simulation, recipe, work, split_reference=reference,
                                    progress_identity=getattr(args, 'progress_identity', None), projections=projections)
    paths = [work / name for name in population.receipt['arrays']]
    if projections is None: paths += [i['path'] for i in simulation]
    paths += [work / 'events.sqlite', work / 'population.json']
    cache[pin] = population, work, paths, file_snapshot(paths)
    return population, work


def prepare_campaign(simulation, config, args, populations):
    """Prepare the two system-specific populations without fitting or scoring."""
    from src import Inputs, Models
    nominal = Models.recipe(config, args.target, 'nominal_retrained_v1')
    required = {sample for samples in nominal['sources'].values() for sample in samples}
    observed = {str(item['source']['sample']).lower() for item in simulation}
    Models.require(observed == required, 'preparation requires the complete configured training sample population')
    result = []; prepared = {}
    reference_pin = Inputs.population_identity(simulation, nominal)
    for name in ('nominal_retrained_v1', 'wideET_v1' if args.target == 'AuAu' else 'ppg12_retrained_v1'):
        recipe = Models.recipe(config, args.target, name)
        pin = Inputs.population_identity(simulation, recipe, reference_pin if name == 'wideET_v1' else None)
        record = dict(status='CHECKED', model=name, input_recipe_sha256=pin,
                      input_events=sum(item['entries'] for item in simulation),
                      event_visits=0, cached_photon_visits=None,
                      et_ranges_gev=recipe['et_ranges_gev'])
        if not args.check:
            population, work = prepared_population(simulation, recipe, config, args, populations)
            record.update(status='PREPARED', population=str(work), rows=population.receipt['rows'],
                          event_split_sha256=population.receipt['event_split_sha256'],
                          cached_photon_visits=population.receipt.get('cached_photon_visits'))
            prepared[name] = population, work
        result.append(record)
    if not args.check:
        preparation_census(simulation, config, args, prepared)
    return result


def preparation_census(simulation, config, args, prepared):
    """Save actual fit-cell sizes from prepared arrays, without fitting a model.

    Count in bounded chunks and reuse identical routing cells across depth
    variants. The report supports batch sizing and exposes sparse cells before
    dispatch; the learner still enforces its ordinary class/split requirements.
    """
    from types import SimpleNamespace
    import numpy as np
    from src import Inputs, Models, Training
    recipes = [Models.recipe(config, args.target, name) for name in Models.campaign_models(args.target)]
    total = sum(len(list(Training.cells(recipe))) for recipe in recipes)
    tracker = Models.WorkProgress(args.models / '.progress' / ('preparation-census-' + args.target.lower() + '.json'),
        total, 'routing cells', 'preparation-census-' + args.target.lower(), getattr(args, 'progress_identity', None))
    counts = {}; fit_models = []; populations = {}
    tracker.write('count-prepared-cells', 0)
    try:
        for name, (population, work) in prepared.items():
            populations[name] = dict(input_recipe_sha256=population.receipt['input_recipe_sha256'],
                population=str(work), population_receipt_sha256=Models.file_hash(work / 'population.json'),
                rows=len(population.y), features=population.x.shape[1],
                array_bytes=sum((work / file).stat().st_size for file in population.receipt['arrays']),
                partition_counts=population.receipt['partition_counts'])
        for recipe in recipes:
            base = recipe['name'] if recipe['name'] in prepared else 'nominal_retrained_v1'
            population, _ = prepared[base]; cells = []
            for cell in Training.cells(recipe):
                key = (base, Models.digest(cell['bounds']))
                if key not in counts:
                    tally = np.zeros(6, np.int64)
                    for begin in range(0, len(population.y), 65536):
                        end = min(begin + 65536, len(population.y)); section = slice(begin, end)
                        view = SimpleNamespace(y=population.y[section], et=population.et[section],
                                               centrality=population.centrality[section])
                        mask = Training.cell_mask(view, cell)
                        part = population.partition[section] if population.partition is not None else np.where(
                            population.train[section], 0, 2)
                        tally += np.bincount((2 * part + view.y)[mask], minlength=6)
                        tracker.stage('count-' + base + '-' + cell['id'], end, len(population.y), 'prepared photons')
                    counts[key] = dict(rows=int(tally.sum()), partition_class_counts=tally.reshape(3, 2).tolist())
                cells.append(dict(cell, **counts[key]))
                tracker.write('count-prepared-cells', tracker.processed + 1)
            fit_models.append(dict(model=recipe['name'], recipe_sha256=Models.digest(recipe),
                                   population=base, cells=cells))
        body = dict(schema='PhotonIDPreparationCensusV1', status='PASS_PREPARED_CELL_CENSUS', system=args.target,
            inputs_identity=Models.digest(Inputs.catalogue_body(simulation, args.target)),
            source_files=len(simulation), source_events=sum(item['entries'] for item in simulation),
            populations=populations, models=fit_models, model_fits=0,
            partition_order=['train', 'validation', 'holdout'], class_order=['background', 'signal'])
        Models.save_json(args.models / args.target / 'preparation_census.json', dict(body, census_sha256=Models.digest(body)))
        tracker.write('complete', total, status='PASS')
    except BaseException as error:
        tracker.fail(error)
        raise


def campaign_projections(simulation, config, args, *, system):
    """Verify complete wide-support cache coverage once for this invocation."""
    from src import Inputs, Models, ProjectionCache
    recipe = Models.projection_recipe(config, system)
    root = args.projection_root or args.models / '.projections'
    return ProjectionCache.load_campaign(simulation, recipe, root,
        progress_path=args.models / '.progress' / ('verify-projections-' + system.lower() + '.json'),
        progress_identity=getattr(args, 'progress_identity', None))


def verify_projections(simulation, config, args, *, system):
    """Verify complete cached rows and save their census without opening trees.

    This uses the preparation verifier, including member hashes and exact source
    coverage. The resulting photon count describes the raw campaign support;
    individual model populations still apply their own domain and weights.
    """
    from src import Inputs, Models
    recipe = Models.projection_recipe(config, system)
    body = dict(schema='PhotonIDProjectionCensusV1', status='CHECKED', system=system,
                inputs_identity=Models.digest(Inputs.catalogue_body(simulation, system)),
                source_files=len(simulation), source_events=sum(i['entries'] for i in simulation),
                projection_definition_sha256=Models.digest(Inputs.population_definition(recipe)),
                et_ranges_gev=recipe['et_ranges_gev'])
    if args.check:
        return body
    projections = campaign_projections(simulation, config, args, system=system)
    receipts = [(key, part['receipt']) for key, part in sorted(projections.items())]
    rows = sum(receipt['rows'] for _, receipt in receipts)
    folders = {(str(part['folder']), part.get('index_file', 'projection.json'))
               for part in projections.values()}
    body.update(status='PASS_VERIFIED_PROJECTION_CENSUS', photon_rows=rows,
                row_bytes=rows * Inputs.projection_dtype(recipe).itemsize,
                packed_bundles=sum(index == 'bundle.json' for _, index in folders),
                source_caches=sum(index == 'projection.json' for _, index in folders),
                projection_receipts_sha256=Models.digest(
                    [(key, receipt['projection_sha256']) for key, receipt in receipts]),
                verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
    report = args.models / system / 'projection_census.json'
    Models.save_json(report, dict(body, census_sha256=Models.digest(body)))
    return dict(body, census=str(report), census_file_sha256=Models.file_hash(report))


def calibration_population(population, recipe):
    """Calibrate and validate over this model's full training/application domain."""
    from src import Models
    Models.require(recipe.get('training_domain', recipe['domain']) == recipe['domain'],
                   'training, working-point calibration and application domains must match')
    return population


def train_release(inputs, config, system, name, args, *, simulation=None, populations=None):
    """Prepare once, fit/check every cell, then seal the complete model release.

    --check returns the resolved request and dependency status before fitting.
    Reuse requires the same request; scoring is a separate phase in main.
    """
    from src import Diagnostics, Inputs, Models, Training, Validation, WorkingPoints
    import numpy as np
    r = Models.recipe(config, system, name)
    if args.release_version != 'v1':
        import re
        Models.require(re.fullmatch(r'v[1-9][0-9]*', args.release_version) is not None, 'invalid release version')
        r['name'] = r['name'].removesuffix('v1') + args.release_version
        r['score_branch'] = r['score_branch'].removesuffix('v1') + args.release_version
    folder = Models.release_dir(args.models, system, r['name'])
    required = {sample for samples in r['sources'].values() for sample in samples}
    Models.require(all(r['sources'][role] for role in ('signal', 'background')),
                   system + ' training source identities are not yet bound in config')
    selected = [i for i in inputs if i['source']['si_di_role'] in ('EMBEDDED', 'SI') and
                str(i['source']['sample']).lower() in required]
    simulation = training_baseline(selected) if simulation is None else simulation
    observed = {str(i['source']['sample']).lower() for i in simulation}
    Models.require(observed == required, 'training sample population differs; missing=' + str(sorted(required - observed)) +
                   ', unexpected=' + str(sorted(observed - required)))
    reference_pin = Inputs.population_identity(simulation, Models.recipe(config, system, 'nominal_retrained_v1')) \
        if name == 'wideET_v1' else None
    request = Models.request_identity(r, simulation, split_reference=reference_pin)
    runtime = Training.runtime_status(r['algorithm'])
    dependency = runtime['status'] == 'PASS'
    plan = dict(model=r['name'], system=system, score_branch=r['score_branch'],
                request_sha256=request, sources=len(simulation), release=str(folder),
                split_reference=reference_pin)
    if args.check:
        return dict(status='CHECKED' if dependency else 'BLOCKED_DEPENDENCY',
                    runtime=runtime, **plan)
    Models.require(dependency, 'training runtime is incomplete: ' + '; '.join(runtime['missing']))
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.training.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (folder / 'manifest.json').exists():
            Models.load_release(folder, request=request)
            return dict(status='REUSED', **plan)
        Models.save_json(folder / 'request.json', dict(request_sha256=request, recipe=r, inputs=simulation,
                                                      split_reference=reference_pin, code_sha256=Models.code_identity(),
                                                      environment=Models.environment()))
        population_cache = {} if populations is None else populations
        population, _ = prepared_population(simulation, r, config, args, population_cache)
        calibration = calibration_population(population, r)
        qa_scores = np.full(len(calibration.y), np.nan, np.float64)
        qa_coverage = np.zeros(len(calibration.y), np.uint8)
        cells = []; all_validation = {}; all_wp = {}
        identity = getattr(args, 'progress_identity', None)
        for cell in Training.cells(r):
            model_path = folder / 'artifacts' / (cell['id'] + '.json')
            receipt_path = folder / 'artifacts' / (cell['id'] + '.validation.json')
            model_path.parent.mkdir(parents=True, exist_ok=True)
            mask = Training.cell_mask(population, cell)
            tracker = Models.WorkProgress(folder / 'artifacts' / 'progress.json', 1,
                                          'validated cells', r['name'] + '/' + cell['id'] + '/validation', identity)
            tracker.write('validate-cell', 0)
            minimum = r['validation']['min_rows_per_class_per_split']
            Models.require(all(int((mask & (population.train == split) & (population.y == k)).sum()) >= minimum
                               for split in (True, False) for k in (0, 1)),
                           'insufficient class/split coverage: ' + cell['id'])
            if receipt_path.exists():
                receipt = Models.read_json(receipt_path)
                Models.require(receipt['request_sha256'] == request, 'cell receipt binds another request')
                for relative, h in receipt['artifacts'].items():
                    Models.require(Models.file_hash(folder / relative) == h, 'completed cell artifact changed')
                calibration_mask = Training.cell_mask(calibration, cell)
                score, validation = Validation.check(calibration, calibration_mask, Training.evaluator(model_path),
                                                     Training.evaluator(model_path), r, progress=tracker.stage)
                wp = WorkingPoints.derive(calibration, calibration_mask, score, cell, r, progress=tracker.stage)
                artifact = Models.read_json(model_path)
                if 'native_export_validation' in artifact:
                    validation['native_export'] = artifact['native_export_validation']
                Models.require(validation == receipt['validation'] and wp == receipt['working_points'], 'cell revalidation differs')
            else:
                Models.require(not model_path.exists(), 'interrupted fit lacks a receipt; retain it and use a new release version')
                artifact, trainer = Training.fit(population, mask, r, model_path, cell['id'], progress_identity=identity)
                calibration_mask = Training.cell_mask(calibration, cell)
                score, validation = Validation.check(calibration, calibration_mask, trainer, Training.evaluator(model_path), r, progress=tracker.stage)
                wp = WorkingPoints.derive(calibration, calibration_mask, score, cell, r, progress=tracker.stage)
                if 'native_export_validation' in artifact:
                    validation['native_export'] = artifact['native_export_validation']
                files = [model_path]
                for suffix in ('.xgboost.json', '.tmva.root', '.history.json'):
                    product = model_path.with_suffix(suffix)
                    if product.exists(): files.append(product)
                receipt = dict(status='PASS', request_sha256=request, validation=validation, working_points=wp,
                               artifacts={str(p.relative_to(folder)): Models.file_hash(p) for p in files})
                Models.save_json(receipt_path, receipt)
            cells.append(dict(cell, artifact=str(model_path.relative_to(folder))))
            all_validation[cell['id']] = validation; all_wp[cell['id']] = wp
            Models.require(not np.any(qa_coverage[calibration_mask]), 'QA model cells overlap')
            qa_scores[calibration_mask] = score
            qa_coverage[calibration_mask] += 1
            tracker.write('cell-validated', 1)
        Models.require(np.all(qa_coverage == 1), 'QA model cells leave uncovered population rows')
        Diagnostics.save(folder, calibration, r, qa_scores, cells, config, request, progress_identity=identity)
        Models.save_json(folder / 'working_points.json', dict(schema='PhotonIDWorkingPointsV1',
                         calibration_domain=r['domain'], cells=all_wp))
        Models.save_json(folder / 'validation.json', dict(schema='PhotonIDValidationV1', status='PASS', cells=all_validation,
                         population=population.receipt, calibration_population=calibration.receipt,
                         evaluation_domain=r['domain'], qualification='local numerical; farm and physics acceptance separate'))
        Models.save_json(folder / 'selection.json', dict(schema='PhotonIDSelectionV1', system=system,
                         name=r['name'], score_branch=r['score_branch'], state_branch=None,
                          working_points_file='working_points.json', domain=r['domain'],
                          et_ranges_gev=r['et_ranges_gev'], cells=cells,
                         non_tight_lower=0.1, tight='score > WP70', non_tight='0.1 < score < WP80',
                         downstream_adapter='PhotonJetSelectionV1: photon_id_selection'))
        artifacts = {str(p.relative_to(folder)): Models.file_hash(p) for p in folder.rglob('*')
                     if p.is_file() and '.work' not in p.parts and p.name not in
                     ('manifest.json', '.training.lock', 'progress.json', 'score_progress.json')}
        manifest = dict(schema='PhotonIDReleaseV1', status='PASS', name=r['name'], system=system,
                        score_branch=r['score_branch'], request_sha256=request, recipe=r, cells=cells,
                        artifacts=artifacts, split_sha256=population.receipt['split_sha256'],
                        calibration_population_sha256=calibration.receipt['input_recipe_sha256'],
                        calibration_split_sha256=calibration.receipt['split_sha256'],
                        code_sha256=Models.code_identity(), environment=Models.environment(),
                        feature_missing_value_rule='reject invalid/non-finite H70 inputs in model domain; unscored outside domain')
        manifest['split_reference'] = population.receipt['split_reference']
        manifest['release_sha256'] = Models.digest(manifest)
        Models.save_json(folder / 'manifest.json', manifest)
        Models.load_release(folder, request=request)
    return dict(status='PASS', **plan)


def main(argv=None):
    """Dispatch one model request and score complete groups after release validation."""
    args = arguments(argv)
    from src import Models
    campaign_names = Models.campaign_models(args.target) if args.command == 'trainCampaign' else []
    if args.model:
        Models.require(len(args.model) == len(set(args.model)), 'duplicate --model recipe')
        Models.require(set(args.model) <= set(campaign_names),
                       'unknown campaign recipe: ' + ', '.join(sorted(set(args.model) - set(campaign_names))))
        campaign_names = [name for name in campaign_names if name in args.model]
    config = Models.load_config(args.config)
    args.progress_identity = Models.progress_identity(Models.yaml.safe_load(args.progress_context.read_text())) \
        if args.progress_context else None
    if args.command == 'planCampaign':
        result = Models.campaign_plan(config, args.group, args.target)
    elif args.command == 'listModels':
        result = Models.list_models(args.models, args.target)
    elif args.command == 'collectQA':
        from src import Diagnostics
        version = args.release_version
        releases = [Models.release_dir(args.models, args.target, name.removesuffix('v1') + version)
                    for name in Models.campaign_models(args.target)]
        result = dict(status='CHECKED' if args.check else 'INDEXED',
                      model_qa_collection=str(Diagnostics.collection(args.models, releases, system=args.target,
                                                                   auau_reference=args.auau_reference, check_only=args.check,
                                                                   progress_identity=args.progress_identity)))
    else:
        from src import Diagnostics, Inputs, Scoring, Training
        systems = ['AuAu', 'PP'] if args.target == 'AuAuAndPP' else [
            args.target if args.command in ('reTrainNominal', 'inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign', 'trainCampaign') else args.system]
        names = {'reTrainNominal': ['nominal_retrained_v1'], 'prepareCampaign': ['nominal_retrained_v1', 'wideET_v1'],
                 'inspectCampaign': [], 'projectCampaign': [], 'verifyProjections': [], 'trainCampaign': campaign_names, 'trainNew': {
            'centAndETdepModels': ['centCoarse_v1', 'centFine_v1', 'ET_v1', 'ETcentCoarse_v1', 'ETcentFine_v1'],
            'routingStudy': Models.ROUTING_MODELS,
            'rangeStudy': ['wideET_v1'],
            'depthStudy': Models.DEPTH_STUDY,
            'MLP': ['MLP_v1'], 'logReg': ['logReg_v1'], 'PPG12': ['ppg12_retrained_v1']}.get(args.target, []), 'score': [args.target]}[args.command]
        result = []
        for system in systems:
            if args.command != 'score':
                algorithms = {'prepare'} if args.command in ('inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign') else {
                    Models.recipe(config, system, name)['algorithm'] for name in names}
                runtime = [Training.runtime_status(algorithm) for algorithm in sorted(algorithms)]
                unavailable = [report for report in runtime if report['status'] != 'PASS']
                if unavailable:
                    Models.require(args.check, 'training runtime is incomplete: ' +
                                   '; '.join(message for report in unavailable for message in report['missing']))
                    result.append(dict(status='BLOCKED_DEPENDENCY', system=system, runtime=unavailable,
                                       input_status='NOT_READ'))
                    continue
            groups = args.group or [Path(p) for p in config['datasets'][system]]
            # For a paired command, configuration owns the separate system collections.
            Models.require(not (len(systems) > 1 and args.group), 'paired systems use separate config dataset groups')
            inspection = {} if args.check or args.command == 'score' else dict(
                progress_path=args.models / '.progress' / 'inputs.json', progress_identity=args.progress_identity)
            # Source hashes are checked by the bounded projection invocation.
            # Subsequent training verifies the GPFS cache, not the TG payload.
            if args.command != 'score': inspection['verify_payload'] = False
            if args.require_roster: inspection['require_population'] = True
            inputs = Inputs.load_catalogue(args.inputs, system) if args.inputs else Inputs.published_inputs(groups, system, **inspection)
            simulation = None; populations = {}
            if args.command != 'score':
                required = {s for samples in config['sources'][system].values() for s in samples}
                selected = [i for i in inputs if
                    i['source']['si_di_role'] in ('EMBEDDED', 'SI') and
                     str(i['source']['sample']).lower() in required]
                simulation = selected if args.inputs else training_baseline(selected, verify_payload=False)
                Models.require({str(i['source']['sample']).lower() for i in simulation} == required,
                               'campaign requires the complete configured training sample population')
            if args.command == 'inspectCampaign':
                body = Inputs.catalogue_body(simulation, system)
                record = dict(status='CHECKED' if args.check else 'INSPECTED',
                              input_events=body['input_events'], files=body['files'], inputs_sha256=Models.digest(body))
                if not args.check: record['inputs'] = str(Inputs.save_catalogue(simulation, system, args.models))
                result.append(record)
                continue
            if args.command == 'projectCampaign':
                recipe = Models.projection_recipe(config, system)
                shards = Inputs.projection_shards(simulation, args.max_events)
                if args.shard_index is not None:
                    Models.require(args.shard_index < len(shards), 'projection shard index is outside the plan')
                    selected_shards = [shards[args.shard_index]]
                else: selected_shards = shards
                if args.check:
                    result.append(dict(status='CHECKED', input_events=sum(i['entries'] for i in simulation),
                                       shards=selected_shards, shard_count=len(shards), max_events=args.max_events,
                                       cache_layout=args.cache_layout, persistent_files_per_shard=4 if args.cache_layout == 'packed' else None,
                                       et_ranges_gev=recipe['et_ranges_gev']))
                else:
                    Models.require(args.shard_index is not None or len(shards) == 1,
                                   'choose one --shard-index from the projection plan; full-source fanout is separate')
                    root = args.projection_root or args.models / '.projections'
                    if args.cache_layout == 'packed':
                        from src import ProjectionCache
                        result.append(ProjectionCache.project_shard(simulation, recipe, root, selected_shards[0],
                            scratch_root=args.scratch_root, progress_identity=args.progress_identity))
                    else:
                        result.append(Inputs.project_shard(simulation, recipe, root, selected_shards[0],
                                                          progress_identity=args.progress_identity))
                continue
            if args.command == 'verifyProjections':
                result.append(verify_projections(simulation, config, args, system=system))
                continue
            args.projections = None
            args.use_projection_cache = not args.check and args.command != 'score'
            if args.command == 'prepareCampaign':
                result.extend(prepare_campaign(simulation, config, args, populations))
                continue
            for name in names:
                if args.command == 'score':
                    folder = Models.release_dir(args.models, system, name)
                    manifest = Models.load_release(folder)
                    record = dict(status='CHECKED', model=name, groups=len(groups)) if args.check else dict(
                        status='PASS', model=name, groups=Scoring.score_groups(inputs, folder, code_root=args.code_root))
                else:
                    record = train_release(inputs, config, system, name, args,
                                           simulation=simulation, populations=populations)
                    if not args.check and not args.no_score:
                        # Resolve again after each family: a previous family may have advanced CURRENT.
                        current = Inputs.published_inputs(groups, system)
                        record['augmentation'] = Scoring.score_groups(current, record['release'], code_root=args.code_root)
                result.append(record)
            # Individually scheduled recipes leave collection to collectQA,
            # after every campaign release has been validated and sealed.
            if args.command == 'trainCampaign' and not args.check and not args.model:
                collection = Diagnostics.collection(args.models, [record['release'] for record in result], system=system,
                                                    auau_reference=args.auau_reference, progress_identity=args.progress_identity)
                for record in result: record['model_qa_collection'] = str(collection)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, ModuleNotFoundError, BlockingIOError) as exc:
        print('PhotonID blocked: ' + str(exc), file=sys.stderr)
        raise SystemExit(2)
