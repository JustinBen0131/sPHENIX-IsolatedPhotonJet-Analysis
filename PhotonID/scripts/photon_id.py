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
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))


def arguments(argv=None):
    parser = argparse.ArgumentParser(description='Train and score named photon models from completed eventTrees.')
    parser.add_argument('command', choices=['reTrainNominal', 'trainNew', 'inspectCampaign', 'projectCampaign', 'verifyProjections',
                                         'prepareCampaign', 'prepareCells', 'trainCampaign', 'diagnoseCampaign', 'finalizeCampaign', 'collectQA', 'score', 'listModels', 'planCampaign'])
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
    parser.add_argument('--model', action='append', default=[], help='trainCampaign/prepareCells recipe; repeat to select several')
    parser.add_argument('--cell', help='trainCampaign: one exact cell from prepared compact inputs')
    parser.add_argument('--assemble-only', action='store_true', help='trainCampaign: seal completed cells; refuse any missing fit')
    parser.add_argument('--summary', action='store_true', help='planCampaign: compact human-readable preview')
    parser.add_argument('--execution-plan', action='store_true', help='planCampaign: preview independent cells and release barriers; no dispatch')
    parser.add_argument('--census', type=Path, help='planCampaign --execution-plan: saved preparation_census.json for exact cell counts')
    parser.add_argument('--fit-code-manifest', type=Path, help='diagnoseCampaign: original code manifest when the diagnostic source differs from the fit')
    parser.add_argument('--prepared-population', type=Path, help='finalizeCampaign: reuse the original completed numeric population after a preparation-code update')
    parser.add_argument('--projection-code-manifest', type=Path, help='original extraction code manifest for packed caches predating recorded code provenance')
    args = parser.parse_args(argv)
    if args.command == 'reTrainNominal' and args.target not in ('AuAu', 'PP', 'AuAuAndPP'):
        parser.error('reTrainNominal requires AuAu, PP or AuAuAndPP')
    if args.command == 'trainNew' and not ((args.target in ('centAndETdepModels', 'routingStudy', 'rangeStudy', 'depthStudy', 'MLP', 'logReg') and args.system == 'AuAu') or (args.target == 'PPG12' and args.system == 'PP')):
        parser.error('trainNew requires an AuAu study or PPG12 PP')
    if args.command in ('planCampaign', 'inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign', 'prepareCells', 'trainCampaign', 'diagnoseCampaign', 'finalizeCampaign', 'collectQA') and (args.target not in ('AuAu', 'PP') or args.system):
        parser.error(args.command + ' requires AuAu or PP')
    if args.command in ('inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign', 'prepareCells', 'trainCampaign') or args.inputs:
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
    if args.command != 'projectCampaign' and args.cache_layout != 'packed':
        parser.error('--cache-layout applies only to projectCampaign')
    if args.scratch_root is not None and args.command not in (
            'projectCampaign', 'prepareCampaign', 'prepareCells', 'trainCampaign', 'reTrainNominal',
            'trainNew', 'diagnoseCampaign', 'finalizeCampaign', 'collectQA'):
        parser.error('--scratch-root applies only to extraction or population-consuming commands')
    if args.scratch_root is not None and args.cache_layout == 'source':
        parser.error('--scratch-root requires packed projection storage')
    if not 0 < args.max_events <= 100000 or (args.shard_index is not None and args.shard_index < 0):
        parser.error('projection event bound must be 1..100000 and shard index must be nonnegative')
    if args.auau_reference and (args.command not in ('trainCampaign', 'collectQA') or args.target != 'PP'):
        parser.error('--auau-reference applies to trainCampaign PP or collectQA PP')
    if args.model and args.command not in ('trainCampaign', 'prepareCells', 'diagnoseCampaign', 'finalizeCampaign'):
        parser.error('--model applies only to trainCampaign, prepareCells, diagnoseCampaign or finalizeCampaign')
    if args.command in ('diagnoseCampaign', 'finalizeCampaign') and (len(args.model) != 1 or not args.inputs):
        parser.error(args.command + ' requires one --model and the original --inputs catalogue')
    if args.fit_code_manifest and args.command not in ('diagnoseCampaign', 'finalizeCampaign'):
        parser.error('--fit-code-manifest applies only to diagnoseCampaign or finalizeCampaign')
    if args.prepared_population and (args.command != 'finalizeCampaign' or not args.fit_code_manifest):
        parser.error('--prepared-population requires finalizeCampaign and the original --fit-code-manifest')
    if args.projection_code_manifest and args.command not in (
            'verifyProjections', 'prepareCampaign', 'prepareCells', 'trainCampaign', 'reTrainNominal', 'trainNew', 'diagnoseCampaign'):
        parser.error('--projection-code-manifest applies only to packed-cache consumers')
    if (args.cell or args.assemble_only) and args.command != 'trainCampaign':
        parser.error('--cell and --assemble-only apply only to trainCampaign')
    if args.cell and (len(args.model) != 1 or not args.inputs or args.assemble_only):
        parser.error('--cell requires one --model and --inputs; do not combine with --assemble-only')
    if args.require_roster and args.command != 'inspectCampaign':
        parser.error('--require-roster applies only to inspectCampaign')
    if args.command != 'planCampaign' and (args.summary or args.execution_plan or args.census):
        parser.error('--summary, --execution-plan and --census apply only to planCampaign')
    if args.census and not args.execution_plan:
        parser.error('--census requires --execution-plan')
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


def retained_population(simulation, recipe, args):
    """Verify the original numeric arrays once; require no retired event store.

    This recovery path preserves the original preparation identity. It opens
    only its receipt and twelve arrays, and is available only for completion
    of an already fitted model whose original source manifest is supplied.
    """
    from src import Inputs, Models
    Models.require(args.command == 'finalizeCampaign' and args.fit_code_manifest,
                   'retained population requires calibration-only completion')
    work = Path(args.prepared_population).absolute()
    Models.require(work.resolve() == work and work.is_dir(), 'retained population directory is aliased or missing')
    path = work / 'population.json'
    Models.require(path.is_file() and not path.is_symlink(), 'original population receipt is missing or aliased')
    receipt_hash = Models.file_hash(path)
    receipt = Models.read_json(path)
    definition = Inputs.population_definition(recipe)
    Models.require(receipt.get('schema') == 'PhotonIDPopulationV1' and receipt.get('sources') == simulation and
                   receipt.get('definition') == definition and receipt.get('split_reference') is None,
                   'retained population differs from the original sources, recipe or split reference')
    manifest = Models.read_json(args.fit_code_manifest)
    Models.require(manifest.get('schema') == 'PhotonIDPublicCodeSnapshotV1', 'original fit code manifest is invalid')
    original = {Path(row['path']).name: row['sha256'] for row in manifest['members']
                if row['path'].startswith('PhotonID/src/') or row['path'] == 'PhotonID/scripts/photon_id.py'}
    expected = Models.digest(dict(inputs=simulation, **definition, preparation_code=original['Inputs.py']))
    Models.require(receipt['input_recipe_sha256'] == expected,
                   'retained population does not bind the original preparation code')
    names = {'x', 'y', 'et', 'centrality', 'weights', 'train', 'partition',
             'split_rank', 'source_index', 'entry', 'photon', 'event_key'}
    files = receipt['array_files']
    Models.require(set(files) == names and set(files.values()) == set(receipt['arrays']) and
                   len(set(files.values())) == len(names), 'retained population array roster is incomplete')
    for filename, checksum in receipt['arrays'].items():
        array = work / filename
        Models.require(Path(filename).name == filename and array.is_file() and not array.is_symlink() and
                       Models.file_hash(array) == checksum, 'retained population array changed: ' + filename)
    population = Inputs.load_population(work, receipt)
    Models.require(population.x.shape == (receipt['rows'], len(recipe['features'])) and
                   all(len(getattr(population, key)) == receipt['rows'] for key in
                       ('y', 'et', 'centrality', 'weights', 'train', 'partition', 'split_rank')),
                   'retained population array dimensions differ')
    Models.require(Models.file_hash(path) == receipt_hash, 'original population receipt changed during reuse')
    args.population_reuse = dict(population_receipt_sha256=receipt_hash,
                                original_population_sha256=expected,
                                numeric_split_sha256=Inputs.split_witness(receipt),
                                arrays=receipt['arrays'])
    return population, work


def prepared_population(simulation, recipe, config, args, cache):
    """Reuse verified matrices within one invocation; widen a frozen nominal split.

    Cached reuse checks source and matrix file metadata without rereading all
    ROOT bytes for every learner. A later invocation performs full hash checks.
    """
    from src import Inputs, Models
    if getattr(args, 'prepared_population', None) is not None:
        return retained_population(simulation, recipe, args)
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
                                    progress_identity=getattr(args, 'progress_identity', None), projections=projections,
                                    scratch_root=getattr(args, 'scratch_root', None) or os.environ.get('_CONDOR_SCRATCH_DIR'))
    paths = [work / name for name in population.receipt['arrays']]
    if projections is None: paths += [i['path'] for i in simulation]
    paths += [work / 'population.json']
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
        progress_identity=getattr(args, 'progress_identity', None),
        code_manifest=getattr(args, 'projection_code_manifest', None))


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


def release_recipe(config, system, name, args):
    """Resolve the named scientific recipe and explicit output version."""
    from src import Models
    r = Models.recipe(config, system, name)
    if args.release_version != 'v1':
        import re
        Models.require(re.fullmatch(r'v[1-9][0-9]*', args.release_version) is not None, 'invalid release version')
        r['name'] = r['name'].removesuffix('v1') + args.release_version
        r['score_branch'] = r['score_branch'].removesuffix('v1') + args.release_version
    return r


def release_request(inputs, config, system, name, args, *, simulation=None):
    """Resolve the one authoritative recipe, input identity and runtime request."""
    from src import Inputs, Models, Training
    r = release_recipe(config, system, name, args)
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
    runtime = Training.runtime_status(r['algorithm'])
    request = Models.request_identity(r, simulation, split_reference=reference_pin)
    plan = dict(model=r['name'], system=system, score_branch=r['score_branch'],
                request_sha256=request, sources=len(simulation), release=str(folder),
                split_reference=reference_pin)
    return r, folder, request, runtime, plan, simulation


def validated_cell(population, r, cell, folder, request, identity, *, cell_input_sha=None, allow_fit=True,
                   retained_artifacts=None):
    """Fit or revalidate one complete cell using the established scientific code.

    Existing receipts must bind this exact request and input. Assembly sets
    allow_fit=False: a missing cell is an error, never an implicit training job.
    """
    from src import Models, Training, Validation, WorkingPoints
    model_path = folder / 'artifacts' / (cell['id'] + '.json')
    receipt_path = folder / 'artifacts' / (cell['id'] + '.validation.json')
    model_path.parent.mkdir(parents=True, exist_ok=True)
    mask = Training.cell_mask(population, cell)
    tracker = Models.WorkProgress(model_path.with_suffix('.validation.progress.json'), 1,
                                  'validated cells', r['name'] + '/' + cell['id'] + '/validation', identity)
    tracker.write('validate-cell', 0)
    minimum = r['validation']['min_rows_per_class_per_split']
    Models.require(all(int((mask & (population.train == split) & (population.y == k)).sum()) >= minimum
                       for split in (True, False) for k in (0, 1)),
                   'insufficient class/split coverage: ' + cell['id'])
    if receipt_path.exists():
        receipt = Models.read_json(receipt_path)
        Models.require(receipt.get('status') == 'PASS' and receipt['request_sha256'] == request and
                       receipt.get('cell') == cell and receipt.get('cell_input_sha256') == cell_input_sha,
                       'cell receipt binds another request, routing cell or compact population')
        products = [model_path] + ([model_path.with_suffix(s) for s in
                    ('.xgboost.json', '.tmva.root', '.history.json')] if r['algorithm'] == 'bdt' else [])
        Models.require(set(receipt['artifacts']) == {str(p.relative_to(folder)) for p in products},
                       'completed cell has missing or unexpected artifacts')
        for relative, h in receipt['artifacts'].items():
            Models.require(Models.file_hash(folder / relative) == h, 'completed cell artifact changed')
        calibration_mask = mask
        score, validation = Validation.check(population, calibration_mask, Training.evaluator(model_path),
                                             Training.evaluator(model_path), r, progress=tracker.stage,
                                             strict_performance=False)
        wp = WorkingPoints.derive(population, calibration_mask, score, cell, r, progress=tracker.stage,
                                  strict_fit_checks=False)
        artifact = Models.read_json(model_path)
        if 'native_export_validation' in artifact:
            validation['native_export'] = artifact['native_export_validation']
        # Older receipts already passed the same performance limits. Preserve
        # their bytes while adding the separately reported quality to this run.
        compared = validation if 'performance_quality' in receipt['validation'] else {
            k: v for k, v in validation.items() if k != 'performance_quality'}
        Models.require(compared == receipt['validation'] and wp == receipt['working_points'], 'cell revalidation differs')
    else:
        if retained_artifacts is not None:
            Models.require(not allow_fit, 'retained calibration must disable fitting')
            Models.require(all(Models.file_hash(folder / p) == h for p, h in retained_artifacts.items()),
                           'retained fit differs from its diagnostic hashes')
            artifact = Models.read_json(model_path)
            trainer = Training.evaluator(model_path)
        else:
            Models.require(allow_fit, 'release assembly requires every validated cell: ' + cell['id'])
            Models.require(not model_path.exists(), 'interrupted fit lacks a receipt; retain it and use a new release version')
            artifact, trainer = Training.fit(population, mask, r, model_path, cell['id'], progress_identity=identity,
                                            progress_path=model_path.with_suffix('.progress.json'))
        calibration_mask = mask
        score, validation = Validation.check(population, calibration_mask, trainer, Training.evaluator(model_path), r,
                                             progress=tracker.stage, strict_performance=False)
        wp = WorkingPoints.derive(population, calibration_mask, score, cell, r, progress=tracker.stage,
                                  strict_fit_checks=False)
        if 'native_export_validation' in artifact:
            validation['native_export'] = artifact['native_export_validation']
        files = [model_path]
        for suffix in ('.xgboost.json', '.tmva.root', '.history.json'):
            product = model_path.with_suffix(suffix)
            if product.exists(): files.append(product)
        receipt = dict(schema='PhotonIDCellValidationV1', status='PASS', request_sha256=request,
                       cell=cell, cell_input_sha256=cell_input_sha, validation=validation, working_points=wp,
                       artifacts={str(p.relative_to(folder)): Models.file_hash(p) for p in files})
        Models.save_json(receipt_path, receipt)
    tracker.write('complete', 1, status='PASS')
    return score, validation, wp


def verify_retained_request(folder, r, simulation, plan, args, *, finalizing=False):
    """Bind an original fit independently of later calibration/driver changes."""
    from src import Models
    original = Models.read_json(folder / 'request.json')
    Models.require(original.get('recipe') == r and original.get('inputs') == simulation and
                   original.get('environment') == Models.environment(),
                   'retained fit binds another recipe, source population or runtime')
    fit_code = original.get('code_sha256')
    Models.require(isinstance(fit_code, dict) and fit_code, 'retained fit lacks its source hashes')
    current_code = Models.code_identity()
    if fit_code != current_code:
        source_manifest = getattr(args, 'fit_code_manifest', None)
        Models.require(source_manifest is not None, 'changed diagnostic source requires the original --fit-code-manifest')
        manifest = Models.read_json(source_manifest)
        Models.require(manifest.get('schema') == 'PhotonIDPublicCodeSnapshotV1', 'original fit code manifest is invalid')
        source_hashes = {Path(row['path']).name: row['sha256'] for row in manifest['members']
                         if row['path'].startswith('PhotonID/src/') or row['path'] == 'PhotonID/scripts/photon_id.py'}
        Models.require(source_hashes == fit_code, 'original fit code manifest differs from the retained request')
    if finalizing:
        changed_coordination = {'WorkingPoints.py', 'Validation.py', 'photon_id.py'}
        if getattr(args, 'prepared_population', None) is not None:
            # The original arrays bypass preparation. Recipe/runtime, original
            # request, diagnosed models and every array are still checked.
            changed_coordination.update(('Inputs.py', 'Models.py', 'ProjectionCache.py'))
        Models.require(set(fit_code) == set(current_code) and all(
            fit_code[key] == current_code[key] for key in fit_code if key not in changed_coordination),
            'retained completion requires unchanged fitter, evaluator and inputs')
    fit_request = dict(recipe=r, inputs=simulation, code=fit_code, environment=original['environment'])
    if plan['split_reference'] is not None:
        Models.require(original.get('split_reference') == plan['split_reference'], 'retained fit split reference differs')
        fit_request['split_reference'] = plan['split_reference']
    Models.require(Models.digest(fit_request) == original.get('request_sha256'), 'retained fit request hash is inconsistent')
    return original


def retained_completion(folder, r, simulation, plan, args):
    """Require a complete diagnosed fit before calibration-only completion."""
    from src import Models, Training
    original = verify_retained_request(folder, r, simulation, plan, args, finalizing=True)
    diagnosis_path = folder / 'validation_diagnosis.json'
    diagnosis = Models.read_json(diagnosis_path)
    Models.require(diagnosis.get('schema') == 'PhotonIDRetainedFitDiagnosisV1' and
                   diagnosis.get('status') == 'REPORT_WRITTEN' and diagnosis.get('model_fits') == 0 and
                   diagnosis.get('scientific_acceptance') is False and
                   diagnosis.get('request_sha256') == original['request_sha256'] and
                   diagnosis.get('fit_code_sha256') == original['code_sha256'],
                   'retained completion requires its original completed diagnosis')
    expected = list(Training.cells(r))
    cells = diagnosis['cells']
    Models.require([c['cell'] for c in cells] == expected and r['algorithm'] == 'bdt',
                   'retained completion requires every expected BDT cell')
    for row in cells:
        path = folder / 'artifacts' / (row['cell']['id'] + '.json')
        products = [path] + [path.with_suffix(s) for s in ('.xgboost.json', '.tmva.root', '.history.json')]
        Models.require(set(row['artifacts']) == {str(p.relative_to(folder)) for p in products},
                       'retained diagnosis has missing or unexpected fitted artifacts')
        Models.require(all(Models.file_hash(folder / p) == h for p, h in row['artifacts'].items()),
                       'retained artifacts changed after diagnosis')
        Models.require(Models.read_json(path.with_suffix('.history.json')).get('rounds') == r['learner']['n_estimators'] and
                       Models.read_json(path).get('native_export_validation', {}).get('status') == 'PASS' and
                       row.get('validation', {}).get('status') == 'PASS',
                       'retained fit lacks complete training and numerical-validation evidence')
    return dict(fit_request_sha256=original['request_sha256'], fit_code_sha256=original['code_sha256'],
                diagnosis_sha256=Models.file_hash(diagnosis_path), population_sha256=diagnosis['population_sha256'],
                cells={c['cell']['id']: c for c in cells})


def diagnose_fit(inputs, config, system, name, args, *, simulation=None, populations=None):
    """Inspect retained fitted cells without fitting, changing cuts or sealing a release.

    The original request and prepared population must match. A failed working
    point records its measured bins; the diagnostic itself grants no acceptance.
    """
    from src import Models, Training, Validation, WorkingPoints
    r, folder, request, runtime, plan, simulation = release_request(inputs, config, system, name, args, simulation=simulation)
    if args.check:
        return dict(status='CHECKED', **plan)
    Models.require(runtime['status'] == 'PASS', 'diagnostic runtime is incomplete')
    Models.require(folder.is_dir() and not (folder / 'manifest.json').exists(),
                   'diagnoseCampaign requires an existing unsealed fitted release')
    original = verify_retained_request(folder, r, simulation, plan, args)
    fit_code = original['code_sha256']
    plan.update(diagnostic_request_sha256=request, request_sha256=original['request_sha256'])
    output = folder / 'validation_diagnosis.json'
    with (folder / '.training.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        Models.require(not output.exists(), 'retained-fit diagnosis already exists; inspect its saved report')
        population, _ = prepared_population(simulation, r, config, args, {} if populations is None else populations)
        tracker = Models.WorkProgress(folder / 'diagnosis.progress.json', 1, 'diagnostic reports',
                                      'diagnose-' + r['name'], getattr(args, 'progress_identity', None))
        tracker.write('inspect-retained-fit', 0)
        reports = []
        for cell in Training.cells(r):
            path = folder / 'artifacts' / (cell['id'] + '.json')
            artifact = Models.read_json(path)
            products = [path]
            if r['algorithm'] == 'bdt':
                products += [path.with_suffix(s) for s in ('.xgboost.json', '.tmva.root', '.history.json')]
                history = Models.read_json(path.with_suffix('.history.json'))
                Models.require(history.get('rounds') == r['learner']['n_estimators'] and
                               artifact.get('native_export_validation', {}).get('status') == 'PASS',
                               'retained BDT has no complete fit/native-export witness')
            hashes = {str(p.relative_to(folder)): Models.file_hash(p) for p in products}
            mask = Training.cell_mask(population, cell)
            evaluator = Training.evaluator(path)
            score, metrics = Validation.check(population, mask, evaluator, evaluator, r,
                                              progress=tracker.stage, strict_performance=False)
            try:
                wp = WorkingPoints.derive(population, mask, score, cell, r, progress=tracker.stage)
                result = dict(status='WORKING_POINTS_PASS', working_points=wp)
            except WorkingPoints.CalibrationError as error:
                result = dict(status='WORKING_POINTS_FAILED', error=str(error), measured=error.report)
            Models.require(all(Models.file_hash(folder / relative) == h for relative, h in hashes.items()),
                           'retained fit changed during diagnosis')
            reports.append(dict(cell=cell, artifacts=hashes, validation=metrics, **result))
        report = dict(schema='PhotonIDRetainedFitDiagnosisV1', status='REPORT_WRITTEN',
                      scientific_acceptance=False, model_fits=0, request_sha256=original['request_sha256'],
                      diagnostic_request_sha256=request, fit_code_sha256=fit_code,
                      diagnostic_code_sha256=Models.code_identity(),
                      population_sha256=population.receipt['input_recipe_sha256'], cells=reports)
        Models.save_json(output, report)
        tracker.write('complete', 1, status='REPORT_WRITTEN')
    return dict(status='REPORT_WRITTEN', scientific_acceptance=False, model_fits=0,
                report=str(output), report_sha256=Models.file_hash(output), **plan)


def prepare_cells(inputs, config, system, name, args, *, simulation=None, populations=None):
    """Copy a routed recipe once; give cell workers a small request pointer."""
    from src import CellInputs, Inputs, Models, Training
    r, folder, request, runtime, plan, simulation = release_request(inputs, config, system, name, args, simulation=simulation)
    Models.require(len(list(Training.cells(r))) > 1, 'prepareCells requires a routed campaign model')
    if args.check: return dict(status='CHECKED', **plan)
    Models.require(runtime['status'] == 'PASS', 'training runtime is incomplete')
    Models.require(args.inputs is not None, 'prepareCells requires the saved --inputs catalogue')
    census = Models.read_json(args.models / system / 'preparation_census.json')
    Models.execution_plan(config, system=system, census=census)
    Models.require(census['inputs_identity'] == Models.digest(Inputs.catalogue_body(simulation, system)),
                   'prepared cell census binds another input catalogue')
    counts = next(c['cells'] for c in census['models'] if c['model'] == name)
    population, work = prepared_population(simulation, r, config, args, {} if populations is None else populations)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.training.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (folder / 'manifest.json').exists():
            Models.load_release(folder, request=request)
            return dict(status='REUSED', **plan)
        Models.save_json(folder / 'request.json', dict(request_sha256=request, recipe=r, inputs=simulation,
            split_reference=plan['split_reference'], code_sha256=Models.code_identity(), environment=Models.environment()))
        result = CellInputs.prepare(args.models, request, r, population, work, counts,
                                   progress_identity=getattr(args, 'progress_identity', None))
        Models.save_json(folder / '.work' / 'cell_inputs.json', dict(schema='PhotonIDCellRequestV1',
            request_sha256=request, recipe=r, code_sha256=Models.code_identity(), environment=Models.environment(),
            catalogue_sha256=Models.file_hash(args.inputs), population=CellInputs.population_witness(population.receipt),
            sources=len(simulation), split_reference=plan['split_reference']))
    return dict(**result, model=r['name'], request_sha256=request)


def train_cell(config, system, name, args):
    """Fit one prepared cell without loading the full source/population arrays."""
    from src import CellInputs, Models, Training
    r = release_recipe(config, system, name, args)
    cells = list(Training.cells(r)); selected = [c for c in cells if c['id'] == args.cell]
    Models.require(len(cells) > 1 and len(selected) == 1, 'choose an exact routed cell from planCampaign')
    runtime = Training.runtime_status(r['algorithm'])
    Models.require(runtime['status'] == 'PASS', 'training runtime is incomplete')
    folder = Models.release_dir(args.models, system, r['name'])
    pin = Models.read_json(folder / '.work' / 'cell_inputs.json')
    Models.require(pin.get('schema') == 'PhotonIDCellRequestV1' and pin['recipe'] == r and
                   pin['code_sha256'] == Models.code_identity() and pin['environment'] == Models.environment(),
                   'prepared cell request has different recipe, source code or runtime')
    Models.require(args.inputs and Models.file_hash(args.inputs) == pin['catalogue_sha256'],
                   'cell input catalogue changed')
    request = pin['request_sha256']; cell = selected[0]
    plan = dict(model=r['name'], system=system, cell=cell['id'], request_sha256=request, release=str(folder))
    if args.check: return dict(status='CHECKED', **plan)
    with (folder / '.training.lock').open('a') as release_lock:
        # Shared for independent cells, exclusive for sequential training/assembly.
        fcntl.flock(release_lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
        if (folder / 'manifest.json').exists():
            Models.load_release(folder, request=request)
            return dict(status='REUSED', **plan)
        with (folder / '.work' / (cell['id'] + '.lock')).open('a') as cell_lock:
            fcntl.flock(cell_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            population, _, input_sha, _ = CellInputs.load(args.models, request, r, cell,
                                                        expected_population=pin['population'])
            validated_cell(population, r, cell, folder, request, getattr(args, 'progress_identity', None),
                           cell_input_sha=input_sha)
    return dict(status='CELL_VALIDATED', release_sealed=False, **plan)


def train_release(inputs, config, system, name, args, *, simulation=None, populations=None):
    """Fit/check a whole recipe, or assemble its exact set of completed cells.

    Sequential and independent-cell execution share the fitter, validation and
    working-point code. A manifest is sealed only after complete QA coverage.
    """
    from src import CellInputs, Diagnostics, Models, Training
    import numpy as np
    r, folder, request, runtime, plan, simulation = release_request(inputs, config, system, name, args, simulation=simulation)
    if args.check:
        return dict(status='CHECKED' if runtime['status'] == 'PASS' else 'BLOCKED_DEPENDENCY', runtime=runtime, **plan)
    Models.require(runtime['status'] == 'PASS', 'training runtime is incomplete: ' + '; '.join(runtime['missing']))
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.training.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (folder / 'manifest.json').exists():
            Models.load_release(folder, request=request)
            return dict(status='REUSED', **plan)
        completion = retained_completion(folder, r, simulation, plan, args) \
            if getattr(args, 'command', None) == 'finalizeCampaign' else None
        current_request = dict(request_sha256=request, recipe=r, inputs=simulation,
            split_reference=plan['split_reference'], code_sha256=Models.code_identity(), environment=Models.environment())
        # Preserve the original fitted request when only calibration/QA changes.
        Models.save_json(folder / ('finalization_request.json' if completion else 'request.json'), current_request)
        expected_cells = list(Training.cells(r))
        observed_receipts = set((folder / 'artifacts').glob('*.validation.json'))
        expected_receipts = {folder / 'artifacts' / (c['id'] + '.validation.json') for c in expected_cells}
        Models.require(observed_receipts <= expected_receipts, 'release has unexpected cell receipts')
        if getattr(args, 'assemble_only', False):
            Models.require(observed_receipts == expected_receipts, 'release assembly requires every validated cell')
        population, work = prepared_population(simulation, r, config, args, {} if populations is None else populations)
        if completion:
            Models.require(population.receipt['input_recipe_sha256'] == completion['population_sha256'],
                           'retained diagnosis used another prepared population')
            if getattr(args, 'population_reuse', None) is not None:
                completion['population_reuse'] = args.population_reuse
        calibration = calibration_population(population, r)
        qa_scores = np.full(len(calibration.y), np.nan, np.float64)
        qa_coverage = np.zeros(len(calibration.y), np.uint8)
        cells = []; all_validation = {}; all_wp = {}
        identity = getattr(args, 'progress_identity', None)
        compact = (folder / '.work' / 'cell_inputs.json').exists()
        artifact_request = completion['fit_request_sha256'] if completion else request
        for cell in expected_cells:
            if compact:
                cell_population, original, input_sha, data = CellInputs.load(args.models, artifact_request, r, cell,
                    expected_population=CellInputs.population_witness(population.receipt))
                CellInputs.verify_origin(data, work, population.receipt)
            else:
                cell_population, original, input_sha = population, np.flatnonzero(Training.cell_mask(population, cell)), None
            score, validation, wp = validated_cell(cell_population, r, cell, folder, artifact_request, identity,
                cell_input_sha=input_sha, allow_fit=not completion and not getattr(args, 'assemble_only', False),
                retained_artifacts=completion['cells'][cell['id']]['artifacts'] if completion else None)
            if completion:
                expected_validation = completion['cells'][cell['id']]['validation']
                compared = {k: v for k, v in validation.items() if k != 'native_export' and
                            (k != 'performance_quality' or k in expected_validation)}
                Models.require(compared == expected_validation,
                               'retained model numerical validation changed since diagnosis')
            cells.append(dict(cell, artifact='artifacts/' + cell['id'] + '.json'))
            all_validation[cell['id']] = validation; all_wp[cell['id']] = wp
            Models.require(not np.any(qa_coverage[original]), 'QA model cells overlap')
            qa_scores[original] = score; qa_coverage[original] += 1
        Models.require(np.all(qa_coverage == 1), 'QA model cells leave uncovered population rows')
        calibration_warnings = [dict(cell=cell_id, **warning) for cell_id, wp in all_wp.items()
                                for warning in wp.get('fit_warnings', [])]
        calibration_quality = dict(status='REVIEW_REQUIRED' if calibration_warnings else 'PASS',
                                   warnings=calibration_warnings, scientific_acceptance=False)
        performance_warnings = [dict(cell=cell_id, **warning) for cell_id, report in all_validation.items()
                                for warning in report.get('performance_quality', {}).get('warnings', [])]
        performance_quality = dict(status='REVIEW_REQUIRED' if performance_warnings else 'PASS',
                                   warnings=performance_warnings, scientific_acceptance=False)
        Diagnostics.save(folder, calibration, r, qa_scores, cells, config, request, progress_identity=identity)
        Models.save_json(folder / 'working_points.json', dict(schema='PhotonIDWorkingPointsV1',
                         calibration_domain=r['domain'], cells=all_wp))
        Models.save_json(folder / 'validation.json', dict(schema='PhotonIDValidationV1', status='PASS', cells=all_validation,
                         population=population.receipt, calibration_population=calibration.receipt,
                         calibration_quality=calibration_quality,
                         performance_quality=performance_quality,
                         evaluation_domain=r['domain'], qualification='local numerical; farm and physics acceptance separate'))
        Models.save_json(folder / 'selection.json', dict(schema='PhotonIDSelectionV1', system=system,
                         name=r['name'], score_branch=r['score_branch'], state_branch=None,
                          working_points_file='working_points.json', domain=r['domain'],
                          et_ranges_gev=r['et_ranges_gev'], cells=cells,
                         non_tight_lower=0.1, tight='score > WP70', non_tight='0.1 < score < WP80',
                         downstream_adapter='PhotonJetSelectionV1: photon_id_selection'))
        artifacts = {str(p.relative_to(folder)): Models.file_hash(p) for p in folder.rglob('*')
                     if p.is_file() and '.work' not in p.parts and p.name not in
                     ('manifest.json', '.training.lock', 'progress.json', 'score_progress.json') and
                     not p.name.endswith('.progress.json')}
        manifest = dict(schema='PhotonIDReleaseV1', status='PASS', name=r['name'], system=system,
                        score_branch=r['score_branch'], request_sha256=request, recipe=r, cells=cells,
                        artifacts=artifacts, split_sha256=population.receipt['split_sha256'],
                        calibration_population_sha256=calibration.receipt['input_recipe_sha256'],
                        calibration_split_sha256=calibration.receipt['split_sha256'],
                        calibration_quality=calibration_quality,
                        performance_quality=performance_quality,
                        code_sha256=Models.code_identity(), environment=Models.environment(),
                        feature_missing_value_rule='reject invalid/non-finite H70 inputs in model domain; unscored outside domain')
        manifest['split_reference'] = population.receipt['split_reference']
        if completion:
            Models.require(all(Models.file_hash(folder / p) == h for row in completion['cells'].values()
                               for p, h in row['artifacts'].items()), 'fitted artifacts changed during completion')
            manifest['retained_fit'] = {k: v for k, v in completion.items() if k != 'cells'}
        manifest['release_sha256'] = Models.digest(manifest)
        Models.save_json(folder / 'manifest.json', manifest)
        Models.load_release(folder, request=request)
    return dict(status='PASS', calibration_quality=calibration_quality,
                performance_quality=performance_quality, **plan)


def main(argv=None):
    """Dispatch one model request and score complete groups after release validation."""
    args = arguments(argv)
    from src import Models
    campaign_names = Models.campaign_models(args.target) if args.command in ('trainCampaign', 'prepareCells', 'diagnoseCampaign', 'finalizeCampaign') else []
    if args.model:
        Models.require(len(args.model) == len(set(args.model)), 'duplicate --model recipe')
        Models.require(set(args.model) <= set(campaign_names),
                       'unknown campaign recipe: ' + ', '.join(sorted(set(args.model) - set(campaign_names))))
        campaign_names = [name for name in campaign_names if name in args.model]
    config = Models.load_config(args.config)
    if args.command == 'prepareCells':
        from src import Training
        routed = [name for name in campaign_names if len(list(Training.cells(Models.recipe(config, args.target, name)))) > 1]
        Models.require(not args.model or routed == campaign_names, 'prepareCells accepts routed models only')
        Models.require(routed, 'no routed campaign models for this system')
        campaign_names = routed
    args.progress_identity = Models.progress_identity(Models.yaml.safe_load(args.progress_context.read_text())) \
        if args.progress_context else None
    if args.command == 'planCampaign':
        result = Models.execution_plan(config, args.group, args.target,
            census=Models.read_json(args.census) if args.census else None) if args.execution_plan else \
            Models.campaign_plan(config, args.group, args.target)
    elif args.command == 'listModels':
        result = Models.list_models(args.models, args.target)
    elif args.command == 'trainCampaign' and args.cell:
        # Compact workers stream-check the catalogue hash without parsing its
        # full roster or opening global population/projection arrays.
        result = [train_cell(config, args.target, campaign_names[0], args)]
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
            args.target if args.command in ('reTrainNominal', 'inspectCampaign', 'projectCampaign', 'verifyProjections', 'prepareCampaign', 'prepareCells', 'trainCampaign', 'diagnoseCampaign', 'finalizeCampaign') else args.system]
        names = {'reTrainNominal': ['nominal_retrained_v1'], 'prepareCampaign': ['nominal_retrained_v1', 'wideET_v1'],
                 'inspectCampaign': [], 'projectCampaign': [], 'verifyProjections': [], 'prepareCells': campaign_names, 'trainCampaign': campaign_names, 'diagnoseCampaign': args.model, 'finalizeCampaign': args.model, 'trainNew': {
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
            if args.command == 'prepareCells':
                result.extend(prepare_cells(inputs, config, system, name, args,
                                           simulation=simulation, populations=populations) for name in names)
                continue
            if args.command == 'diagnoseCampaign':
                result.extend(diagnose_fit(inputs, config, system, name, args,
                                          simulation=simulation, populations=populations) for name in names)
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
    print(Models.campaign_summary(result) if args.summary else json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, ModuleNotFoundError, BlockingIOError) as exc:
        print('PhotonID blocked: ' + str(exc), file=sys.stderr)
        raise SystemExit(2)
