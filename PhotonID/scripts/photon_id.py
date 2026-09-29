#!/usr/bin/env python3
"""Public PhotonID command driver: prepare, fit, validate, calibrate, score.

All algorithms share one eventTree input contract and saved split. Existing
validated releases are reused; an interrupted candidate is never promoted by
filename alone. Completed groups are published by the maintained eventTree API.
The shell entry point selects an already configured Python environment.
"""
from __future__ import annotations

import argparse
import fcntl
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))


def arguments(argv=None):
    parser = argparse.ArgumentParser(description='Train and score named photon models from completed eventTrees.')
    parser.add_argument('command', choices=['reTrainNominal', 'trainNew', 'score', 'listModels'])
    parser.add_argument('target', nargs='?', help='AuAu/PP/AuAuAndPP, model recipe or saved model name')
    parser.add_argument('system', nargs='?', choices=['AuAu', 'PP'])
    parser.add_argument('--config', type=Path, default=HERE / 'config/training.yaml')
    parser.add_argument('--models', type=Path, default=HERE / 'models', help='model release root')
    parser.add_argument('--group', type=Path, action='append', default=[], help='complete published group; repeat as needed')
    parser.add_argument('--code-root', type=Path, default=HERE.parent, help='repository containing shared TreeProduction eventTree module')
    parser.add_argument('--release-version', default='v1', help='explicit version for a new training recipe/input population')
    parser.add_argument('--check', action='store_true', help='read-only input, recipe and dependency preview')
    parser.add_argument('--no-score', action='store_true', help='finish model release without publishing dataset augmentation')
    args = parser.parse_args(argv)
    if args.command == 'reTrainNominal' and args.target not in ('AuAu', 'PP', 'AuAuAndPP'):
        parser.error('reTrainNominal requires AuAu, PP or AuAuAndPP')
    if args.command == 'trainNew' and (args.target not in ('centAndETdepModels', 'MLP', 'logReg') or args.system != 'AuAu'):
        parser.error('trainNew requires centAndETdepModels, MLP or logReg followed by AuAu')
    if args.command == 'score' and (not args.target or not args.system):
        parser.error('score requires a saved model name and AuAu/PP')
    if args.command == 'listModels' and args.target not in (None, 'AuAu', 'PP'):
        parser.error('listModels accepts AuAu or PP')
    return args


def training_baseline(inputs):
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
            cache[(str(group), item['version'])] = directory, m
        directory, m = cache[(str(group), item['version'])]
        roster = {r['logical_id']: r for r in m['files']}
        require(item['logical_id'] in roster, 'augmentation lineage lost a logical source')
        baseline = roster[item['logical_id']]
        path = directory / baseline['relative_path']
        require(file_hash(path) == baseline['sha256'], 'stage-2 training ancestor changed')
        result.append(dict(item, path=str(path), sha256=baseline['sha256'], version=m['version']))
    return result


def train_release(inputs, config, system, name, args):
    from src import Inputs, Models, Training, Validation, WorkingPoints
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
    simulation = training_baseline(selected)
    observed = {str(i['source']['sample']).lower() for i in simulation}
    Models.require(observed == required, 'training sample population differs; missing=' + str(sorted(required - observed)) +
                   ', unexpected=' + str(sorted(observed - required)))
    request = Models.request_identity(r, simulation)
    if r['algorithm'] == 'bdt':
        import importlib.util
        dependency = importlib.util.find_spec('xgboost') is not None
    else:
        dependency = True
    plan = dict(model=r['name'], system=system, score_branch=r['score_branch'],
                request_sha256=request, sources=len(simulation), release=str(folder))
    if args.check:
        return dict(status='CHECKED' if dependency else 'BLOCKED_DEPENDENCY',
                    dependency='available' if dependency else 'XGBoost missing', **plan)
    Models.require(dependency, 'BDT training requires XGBoost; no environment or packages were changed')
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.training.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (folder / 'manifest.json').exists():
            Models.load_release(folder, request=request)
            return dict(status='REUSED', **plan)
        Models.save_json(folder / 'request.json', dict(request_sha256=request, recipe=r, inputs=simulation))
        work = args.models.resolve() / '.populations' / Inputs.population_identity(simulation, r)
        work.mkdir(parents=True, exist_ok=True)
        with (work / '.preparation.lock').open('a') as preparation_lock:
            fcntl.flock(preparation_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            population = Inputs.prepare(simulation, r, work)
        cells = []; all_validation = {}; all_wp = {}
        for cell in Training.cells(r):
            model_path = folder / 'artifacts' / (cell['id'] + '.json')
            receipt_path = folder / 'artifacts' / (cell['id'] + '.validation.json')
            model_path.parent.mkdir(parents=True, exist_ok=True)
            mask = Training.cell_mask(population, cell)
            minimum = r['validation']['min_rows_per_class_per_split']
            Models.require(all(int((mask & (population.train == split) & (population.y == k)).sum()) >= minimum
                               for split in (True, False) for k in (0, 1)),
                           'insufficient class/split coverage: ' + cell['id'])
            if receipt_path.exists():
                receipt = Models.read_json(receipt_path)
                Models.require(receipt['request_sha256'] == request, 'cell receipt binds another request')
                for relative, h in receipt['artifacts'].items():
                    Models.require(Models.file_hash(folder / relative) == h, 'completed cell artifact changed')
                score, validation = Validation.check(population, mask, Training.evaluator(model_path),
                                                     Training.evaluator(model_path), r)
                wp = WorkingPoints.derive(population, mask, score, cell, r)
                Models.require(validation == receipt['validation'] and wp == receipt['working_points'], 'cell revalidation differs')
            else:
                Models.require(not model_path.exists(), 'interrupted fit lacks a receipt; retain it and use a new release version')
                _, trainer = Training.fit(population, mask, r, model_path, cell['id'])
                score, validation = Validation.check(population, mask, trainer, Training.evaluator(model_path), r)
                wp = WorkingPoints.derive(population, mask, score, cell, r)
                files = [model_path]
                booster = model_path.with_suffix('.xgboost.json')
                if booster.exists(): files.append(booster)
                receipt = dict(status='PASS', request_sha256=request, validation=validation, working_points=wp,
                               artifacts={str(p.relative_to(folder)): Models.file_hash(p) for p in files})
                Models.save_json(receipt_path, receipt)
            cells.append(dict(cell, artifact=str(model_path.relative_to(folder))))
            all_validation[cell['id']] = validation; all_wp[cell['id']] = wp
        Models.save_json(folder / 'working_points.json', dict(schema='PhotonIDWorkingPointsV1', cells=all_wp))
        Models.save_json(folder / 'validation.json', dict(schema='PhotonIDValidationV1', status='PASS', cells=all_validation,
                         population=population.receipt, qualification='local numerical; farm and physics acceptance separate'))
        Models.save_json(folder / 'selection.json', dict(schema='PhotonIDSelectionV1', system=system,
                         name=r['name'], score_branch=r['score_branch'], state_branch=None,
                         working_points_file='working_points.json', domain=r['domain'], cells=cells,
                         non_tight_lower=0.1, tight='score > WP70', non_tight='0.1 < score < WP80',
                         downstream_adapter='PhotonJetSelectionV1: photon_id_selection'))
        artifacts = {str(p.relative_to(folder)): Models.file_hash(p) for p in folder.rglob('*')
                     if p.is_file() and '.work' not in p.parts and p.name not in
                     ('manifest.json', '.training.lock', 'progress.json', 'score_progress.json')}
        manifest = dict(schema='PhotonIDReleaseV1', status='PASS', name=r['name'], system=system,
                        score_branch=r['score_branch'], request_sha256=request, recipe=r, cells=cells,
                        artifacts=artifacts, split_sha256=population.receipt['split_sha256'])
        manifest['release_sha256'] = Models.digest(manifest)
        Models.save_json(folder / 'manifest.json', manifest)
        Models.load_release(folder, request=request)
    return dict(status='PASS', **plan)


def main(argv=None):
    args = arguments(argv)
    from src import Inputs, Models, Scoring
    config = Models.load_config(args.config)
    if args.command == 'listModels':
        result = Models.list_models(args.models, args.target)
    else:
        systems = ['AuAu', 'PP'] if args.target == 'AuAuAndPP' else [args.target if args.command == 'reTrainNominal' else args.system]
        names = {'reTrainNominal': ['nominal_retrained_v1'], 'trainNew': {
            'centAndETdepModels': ['centCoarse_v1', 'centFine_v1', 'ET_v1', 'ETcentCoarse_v1', 'ETcentFine_v1'],
            'MLP': ['MLP_v1'], 'logReg': ['logReg_v1']}.get(args.target, []), 'score': [args.target]}[args.command]
        result = []
        for system in systems:
            groups = args.group or [Path(p) for p in config['datasets'][system]]
            # For a paired command, configuration owns the separate system collections.
            Models.require(not (len(systems) > 1 and args.group), 'paired systems use separate config dataset groups')
            inputs = Inputs.published_inputs(groups, system, code_root=args.code_root)
            for name in names:
                if args.command == 'score':
                    folder = Models.release_dir(args.models, system, name)
                    manifest = Models.load_release(folder)
                    record = dict(status='CHECKED', model=name, groups=len(groups)) if args.check else dict(
                        status='PASS', model=name, groups=Scoring.score_groups(inputs, folder, code_root=args.code_root))
                else:
                    record = train_release(inputs, config, system, name, args)
                    if not args.check and not args.no_score:
                        # Resolve again after each family: a previous family may have advanced CURRENT.
                        current = Inputs.published_inputs(groups, system, code_root=args.code_root)
                        record['augmentation'] = Scoring.score_groups(current, record['release'], code_root=args.code_root)
                result.append(record)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, ModuleNotFoundError, BlockingIOError) as exc:
        print('PhotonID blocked: ' + str(exc), file=sys.stderr)
        raise SystemExit(2)
