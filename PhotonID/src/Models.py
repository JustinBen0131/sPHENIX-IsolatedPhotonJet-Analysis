"""Models: immutable photon-ID releases and their scientific definitions.

The command driver uses this module to resolve recipes, verify model bytes and
resume only work bound to the same inputs, code and environment. Generated
artifacts never replace the bundled nominal models. A release becomes usable
only after its numerical and working-point validation has passed.
"""
from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys

import yaml

HERE = Path(__file__).resolve().parents[1]
SYSTEMS = {'AuAu': 'auau', 'PP': 'pp'}
NAMES = {
    'nominal_retrained_v1': ('bdt', 'pho_bdt_score_nominal_retrained_v1', ()),
    'centCoarse_v1': ('bdt', 'pho_bdt_score_cent_coarse_v1', ('centrality_coarse',)),
    'centFine_v1': ('bdt', 'pho_bdt_score_cent_fine_v1', ('centrality_fine',)),
    'ET_v1': ('bdt', 'pho_bdt_score_et_v1', ('et',)),
    'ETcentCoarse_v1': ('bdt', 'pho_bdt_score_et_cent_coarse_v1', ('et', 'centrality_coarse')),
    'ETcentFine_v1': ('bdt', 'pho_bdt_score_et_cent_fine_v1', ('et', 'centrality_fine')),
    'MLP_v1': ('mlp', 'pho_mlp_score_v1', ()),
    'logReg_v1': ('logreg', 'pho_logreg_score_v1', ()),
    'wideET_v1': ('bdt', 'pho_bdt_score_et12to40_v1', ()),
    'ppg12_retrained_v1': ('bdt', 'pho_bdt_score_ppg12_retrained_v1', ()),
}
DEPTH_MODELS = tuple(f'depth{depth}_v1' for depth in range(2, 11) if depth != 4)
NAMES.update({name: ('bdt', f'pho_bdt_score_depth_{name[5:-3]}_v1', ())
              for name in DEPTH_MODELS})
ROUTING_MODELS = ('centCoarse_v1', 'centFine_v1', 'ETcentCoarse_v1', 'ETcentFine_v1')
# The depth-4 point is the nominal release, so the scan never fits it twice.
DEPTH_STUDY = ('nominal_retrained_v1',) + DEPTH_MODELS
CAMPAIGN_MODELS = ('nominal_retrained_v1', 'wideET_v1') + ROUTING_MODELS + DEPTH_MODELS
PP_CAMPAIGN_MODELS = ('nominal_retrained_v1', 'ppg12_retrained_v1')


def campaign_models(system):
    require(system in SYSTEMS, 'unknown collision system')
    return CAMPAIGN_MODELS if system == 'AuAu' else PP_CAMPAIGN_MODELS


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_hash(path, *, progress=None):
    """Stream SHA-256; an optional callback reports bytes read of a fixed size."""
    h = hashlib.sha256()
    total = Path(path).stat().st_size if progress else None
    processed = 0
    if progress: progress(0, total)
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            h.update(block)
            if progress:
                processed += len(block)
                require(processed <= total, 'file grew during hashing')
                progress(processed, total)
    if progress: require(processed == total, 'file size changed during hashing')
    return h.hexdigest()


def save_json(path, value):
    """Create once: an existing different receipt is an error, not an overwrite."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = canonical(value) + '\n'
    try:
        with path.open('x') as out:
            out.write(raw)
            out.flush()
            os.fsync(out.fileno())
    except FileExistsError:
        require(path.read_text() == raw, 'immutable artifact differs: ' + str(path))


def read_json(path):
    return json.loads(Path(path).read_text())


def environment():
    versions = {}
    for name in ('numpy', 'scipy', 'uproot', 'awkward', 'xgboost', 'scikit-learn', 'PyYAML', 'ROOT'):
        module = sys.modules.get({'scikit-learn': 'sklearn', 'PyYAML': 'yaml'}.get(name, name))
        if module is not None and hasattr(module, '__version__'):
            versions[name] = module.__version__
            continue
        if name == 'ROOT' and module is not None and hasattr(module, 'gROOT'):
            versions[name] = module.gROOT.GetVersion()
            continue
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            try:
                versions[name] = importlib.metadata.version('xgboost-cpu') if name == 'xgboost' else None
            except importlib.metadata.PackageNotFoundError:
                versions[name] = None
    return dict(python=sys.version, executable=sys.executable, packages=versions)


def load_config(path=None):
    """Read the maintained recipe/deployment configuration and reject unknown structure."""
    path = Path(path or HERE / 'config/training.yaml').resolve()
    config = yaml.safe_load(path.read_text())
    require(config.get('version') == 1, 'unsupported training configuration')
    require(config['training']['train_fraction'] == .9, 'nominal training fraction must be 0.9')
    require(config['training']['seed'] == 13, 'nominal split seed must be 13')
    require(config['training']['minimum_bias_branch'] == 'cent_mb_decision' and
            config['training']['minimum_bias_pass'] == 1, 'corrected minimum-bias witness changed')
    require(len(config['features']['AuAu']) == 14 and len(config['features']['PP']) == 11,
            'nominal feature counts differ')
    expected = [n for n in config['features']['AuAu'] if n not in
                ('cluster_weta33_cogx', 'cluster_wphi33_cogx', 'centrality')]
    require(config['features']['PP'] == expected, 'pp feature order differs from nominal projection')
    require(config['weights'] == dict(bins=20, eta_range=[-.7, .7], et_cap=800),
            'canonical weighting settings changed')
    for edges in config['bins'].values():
        require(len(edges) >= 2 and all(a < b for a, b in zip(edges, edges[1:])), 'invalid bin edges')
    return config


def recipe(config, system, name):
    """Resolve one named model's features, domain, routing and learner settings."""
    require(system in SYSTEMS and name in NAMES, 'unknown system/model')
    require((system == 'AuAu' and name != 'ppg12_retrained_v1') or
            (system == 'PP' and name in PP_CAMPAIGN_MODELS), 'model is not defined for this system')
    algorithm, branch, axes = NAMES[name]
    domain = dict(et_min=15., et_max_exclusive=35., eta_abs_max=.7, vertex_abs_max_cm=10.)
    if system == 'AuAu':
        domain.update(centrality_min=0., centrality_max_exclusive=80.)
    result = copy.deepcopy(dict(schema='PhotonIDRecipeV1', name=name, system=system,
                algorithm=algorithm, score_branch=branch, features=config['features'][system],
                domain=domain, training_domain=dict(domain), axes={a: config['bins'][a] for a in axes},
                training=config['training'], weights=config['weights'],
                learner=config['learners'][algorithm], validation=config['validation'],
                working_points=config['working_points'],
                sources=config['sources'][system],
                feature_view='H70', isolation_in_training=False))
    if name == 'ppg12_retrained_v1':
        # The historical global fit used all accepted 6--40 GeV rows, including
        # the upper/acceptance boundaries. nextafter represents those closed
        # limits through the common half-open reader and selection interface.
        from math import nextafter, inf
        result['domain'].update(et_min=6., et_max_exclusive=nextafter(40., inf),
                                eta_abs_max=nextafter(.7, inf), vertex_abs_max_cm=nextafter(10., inf))
        result['training_domain'] = dict(result['domain'])
        result['training'] = dict(train_fraction=.7, validation_fraction=.1, test_fraction=.2,
                                  seed=42, split_identity='ppg12_row_stratified_70_10_20',
                                  split_method='row_stratified', weight_dtype='float32',
                                  background_flattening=dict(et_threshold=15., bins=20, seed=42))
        result['learner'] = copy.deepcopy(config['learners']['ppg12'])
        result['working_points']['et_edges'] = [6., 10., 15., 20., 25., 35.]
        reference = HERE / 'models/PP/nominal'
        frozen = read_json(reference / 'manifest.json')
        require(frozen['definition']['features']['order'] == result['features'], 'frozen PPG12 feature order differs')
        result['reference_model_sha256'] = frozen['artifacts']['model.tmva.root']
    elif name == 'wideET_v1':
        result['domain'].update(et_min=12., et_max_exclusive=40.)
        result['training_domain'].update(et_min=12., et_max_exclusive=40.)
    elif name in DEPTH_MODELS:
        result['learner']['max_depth'] = int(name[5:-3])
    low, high = result['domain']['et_min'], result['domain']['et_max_exclusive']
    result['et_ranges_gev'] = dict(
        training=[result['training_domain']['et_min'], result['training_domain']['et_max_exclusive']],
        working_point_calibration=[low, high], application=[low, high],
        analysis=[low + 2., high - 2.], underflow=[low, low + 2.], overflow=[high - 2., high])
    if name == 'ppg12_retrained_v1':
        # This is the historical-model comparison recipe. It does not declare
        # a new measurement interval or invent guard bins for the pp analysis.
        result['et_ranges_gev'].update(analysis=None, underflow=None, overflow=None)
    return result


def projection_recipe(config, system):
    """One raw extraction covers every model; no split or flattening is applied."""
    return recipe(config, system, 'wideET_v1' if system == 'AuAu' else 'ppg12_retrained_v1')


def campaign_plan(config, groups=(), system='AuAu'):
    """Describe the finite BDT matrix without opening trees or fitting models.

    Group paths are declarations, not input qualification. Routing and depth
    recipes share a preparation identity; the wider range has its own population.
    """
    from math import prod
    models = []
    for name in campaign_models(system):
        r = recipe(config, system, name)
        models.append(dict(name=name, score_branch=r['score_branch'],
                           training_et_gev=[r['training_domain']['et_min'], r['training_domain']['et_max_exclusive']],
                           evaluation_et_gev=[r['domain']['et_min'], r['domain']['et_max_exclusive']],
                           et_ranges_gev=r['et_ranges_gev'],
                           max_depth=r['learner']['max_depth'], routing=r['axes'],
                           training=r['training'], learner=r['learner'],
                           component_models=prod(len(edges) - 1 for edges in r['axes'].values()),
                           recipe_sha256=digest(r)))
    selected = list(groups) or config['datasets'][system]
    return dict(schema='PhotonIDCampaignPlanV1', status='PLAN_ONLY', system=system,
                groups=[str(path) for path in selected],
                input_status='DECLARED_UNVERIFIED' if selected else 'UNBOUND',
                sources=config['sources'][system], models=models,
                configurations=len(models), component_models=sum(m['component_models'] for m in models),
                depth4_release='nominal_retrained_v1',
                common_event_split='wideET_v1 inherits nominal assignments; new events use its hash cutoff' if system == 'AuAu' else
                    'canonical pp: event-grouped 90/10; PPG12: separate row-stratified 70/10/20',
                regions=dict(tight='score > WP70(C)', non_tight='0.1 < score < WP80(C)',
                             gap='WP80(C) <= score <= WP70(C), excluded'))


def execution_plan(config, groups=(), system='AuAu', census=None):
    """Preview exact independent fits and release barriers; never dispatch.

    Cell parallelism uses the existing Cartesian cells and frozen population
    weights/splits. Compact cell inputs and concurrent release assembly still
    need qualification. Memory targets are measurement goals, never requests.
    A supplied census is hash/recipe/count checked, not an array revalidation.
    """
    from .Training import cells
    campaign = campaign_plan(config, groups, system)
    names = campaign_models(system)
    recipes = {name: recipe(config, system, name) for name in names}
    populations = ('nominal_retrained_v1', 'wideET_v1' if system == 'AuAu' else 'ppg12_retrained_v1')
    measured_cells = {}
    if census is not None:
        require(isinstance(census, dict), 'execution census must be a JSON object')
        body = {k: v for k, v in census.items() if k != 'census_sha256'}
        require(census.get('schema') == 'PhotonIDPreparationCensusV1' and
                census.get('status') == 'PASS_PREPARED_CELL_CENSUS' and
                census.get('system') == system and census.get('census_sha256') == digest(body),
                'execution census hash, system or status differs')
        require(isinstance(census.get('populations'), dict) and isinstance(census.get('models'), list) and
                all(isinstance(m, dict) and {'model', 'recipe_sha256', 'population', 'cells'} <= set(m)
                    for m in census['models']), 'malformed execution census populations or models')
        require(set(census['populations']) == set(populations) and
                [m['model'] for m in census['models']] == list(names), 'execution census population/model set differs')
        require(all(isinstance(p, dict) and type(p.get('rows')) is int and p['rows'] > 0 and
                isinstance(p.get('partition_counts'), list) and len(p['partition_counts']) == 3 and
                all(type(n) is int and n >= 0 for n in p['partition_counts'])
                for p in census['populations'].values()), 'malformed execution census population totals')
        for model in census['models']:
            r = recipes[model['model']]
            base = model['model'] if model['model'] in populations else populations[0]
            expected = list(cells(r)); observed = model['cells']; total = [0, 0, 0]
            require(isinstance(observed, list) and model['recipe_sha256'] == digest(r) and model['population'] == base and
                    len(observed) == len(expected), 'execution census recipe or cell count differs')
            for saved, cell in zip(observed, expected):
                require(isinstance(saved, dict) and {'id', 'bounds', 'rows', 'partition_class_counts'} <= set(saved),
                        'malformed execution census cell')
                counts = saved['partition_class_counts']
                require(saved['id'] == cell['id'] and saved['bounds'] == cell['bounds'] and
                        isinstance(counts, list) and len(counts) == 3 and all(isinstance(pair, list) and len(pair) == 2 and
                        all(type(n) is int and n >= 0 for n in pair) for pair in counts),
                        'execution census cell identity or class counts differ')
                require(type(saved['rows']) is int and saved['rows'] == sum(map(sum, counts)), 'execution census cell total differs')
                minimum = r['validation']['min_rows_per_class_per_split']
                require(all(counts[0][k] >= minimum and counts[1][k] + counts[2][k] >= minimum and
                            counts[2][k] > 0 for k in (0, 1)), 'execution census cell has insufficient split coverage')
                measured_cells[(model['model'], cell['id'])] = saved
                total = [a + sum(b) for a, b in zip(total, counts)]
            require(total == census['populations'][base]['partition_counts'] and
                    sum(total) == census['populations'][base]['rows'], 'execution census routing coverage differs')
    prepare_id = 'prepare:' + system
    tasks = [dict(id=prepare_id, phase='prepare', depends_on=[], populations=list(populations),
                  reuse_existing=True, status='CENSUS_CHECKED' if census else 'WAITING_FOR_PREPARATION')]
    for name, r in recipes.items():
        routed = bool(r['axes']); base = name if name in populations else populations[0]
        dependency = prepare_id
        if routed:
            dependency = 'cell-inputs:' + name
            tasks.append(dict(id=dependency, phase='cell_inputs', depends_on=[prepare_id], model=name,
                              population=base, preserve=['row_order', 'global_weights', 'event_split', 'original_row_index'],
                              scans_per_population=1, qualification='IMPLEMENTATION_AND_EQUIVALENCE_PENDING'))
        fit_ids = []
        for cell in cells(r):
            fit_id = 'fit:' + name + ':' + cell['id']; fit_ids.append(fit_id)
            saved = measured_cells.get((name, cell['id']))
            tasks.append(dict(id=fit_id, phase='fit', depends_on=[dependency], model=name, cell=cell,
                recipe_sha256=digest(r), population=base, cpus=r['learner']['n_jobs'],
                priority=0 if name == 'nominal_retrained_v1' else 1 if name in populations else 2 if routed else 3,
                rows=saved['rows'] if saved else None,
                partition_class_counts=saved['partition_class_counts'] if saved else None,
                memory=dict(request_mb=None, target_mb=4096 if routed else 8192, ceiling_mb=12288,
                            status='MEASUREMENT_REQUIRED', basis='Exact cell input plus measured fit/export/validation peak'),
                qualification='CELL_EXECUTOR_EQUIVALENCE_REQUIRED' if routed else 'EXISTING_WHOLE_MODEL_PATH',
                task_sha256=digest(dict(system=system, recipe=digest(r), cell=cell))))
        tasks.append(dict(id='release:' + name, phase='assemble', depends_on=fit_ids, model=name,
            expected_cells=list(cells(r)), qualification='ASSEMBLY_EQUIVALENCE_REQUIRED',
            checks=['exact_recipe_and_population', 'one_receipt_per_cell', 'no_overlap_or_missing_rows',
                    'saved_artifact_hashes', 'canonical_row_order_for_diagnostics', 'working_point_validation']))
    tasks.append(dict(id='collect:' + system, phase='collect', depends_on=['release:' + n for n in names],
                      external_dependencies=['accepted_AuAu_nominal_reference'] if system == 'PP' else []))
    fits = [t for t in tasks if t['phase'] == 'fit']
    ids = {t['id'] for t in tasks}
    require(len(ids) == len(tasks) and all(d in ids for t in tasks for d in t['depends_on']),
            'execution plan has duplicate tasks or a missing dependency')
    return dict(schema='PhotonIDExecutionPlanV1', status='DESIGN_ONLY', submit_ready=False, system=system,
        campaign=campaign, tasks=tasks, counts=dict(fits=len(fits), routed_fits=sum(bool(t['cell']['bounds']) for t in fits),
            global_fits=sum(not t['cell']['bounds'] for t in fits), releases=len(names),
            cell_input_passes=sum(t['phase'] == 'cell_inputs' for t in tasks)),
        census_sha256=census['census_sha256'] if census else None,
        input_status='CENSUS_HASH_CHECKED_ARRAYS_NOT_REVERIFIED' if census else 'UNMEASURED',
        supported_execution='One trainCampaign --model invocation per recipe',
        reuse_existing_preparation=True, independent_release_consumption=True,
        invariants=['sources', 'features', 'row_order', 'weights', 'split', 'seed', 'learner', 'threads', 'working_points'],
        qualification_required=['compact_cell_input_equivalence', 'isolated_cell_fit_equivalence',
                                'release_assembly_equivalence', 'full_scale_memory_and_walltime'])


def campaign_summary(plan):
    """Small collaborator preview; planned tasks never imply execution readiness."""
    execution = plan.get('schema') == 'PhotonIDExecutionPlanV1'
    campaign = plan['campaign'] if execution else plan
    labels = dict(nominal_retrained_v1='Canonical', wideET_v1='Wide range',
        centCoarse_v1='Coarse centrality', centFine_v1='Fine centrality',
        ETcentCoarse_v1='ET x coarse centrality', ETcentFine_v1='ET x fine centrality',
        ppg12_retrained_v1='Historical PPG12')
    lines = [f"PhotonID {campaign['system']}: {campaign['configurations']} model recipes, "
             f"{campaign['component_models']} BDT fits", "", 'Model                         Range (GeV)    Fits    CPUs']
    for model in campaign['models']:
        lo, hi = model['training_et_gev']
        label = labels.get(model['name'], 'Depth ' + str(model['max_depth']))
        lines.append(f"{label:<29} {lo:g}-{hi:g}".ljust(46) +
                     f"{model['component_models']:>4}    {model['learner']['n_jobs']:>4}")
    lines += ['', 'Extract once; prepare two populations; fit, validate and calibrate; collect QA.',
              'Depth and routing variants reuse canonical weights and split. Depth 4 is the nominal model.']
    if execution:
        c = plan['counts']
        lines += [f"Candidate scheduling: {c['routed_fits']} independent routing-cell fits and "
                  f"{c['global_fits']} full-population fits; {c['releases']} release barriers.",
                  'Small-cell memory target: 4 GiB. Actual requests remain unset until measured.',
                  'Cell dispatch and assembly await equivalence qualification. This plan cannot submit jobs.']
    else:
        lines += [f"Supported batch scheduling: {campaign['configurations']} independent recipe invocations."]
    return '\n'.join(lines)


def code_identity():
    """Source hashes for the fitting, evaluator and command implementations."""
    code = {p.name: file_hash(p) for p in sorted((HERE / 'src').glob('*.py'))}
    code['photon_id.py'] = file_hash(HERE / 'scripts/photon_id.py')
    return code


def request_identity(recipe, inputs, *, split_reference=None):
    """Bind semantics and backend source bytes; code changes require a new release."""
    request = dict(recipe=recipe, inputs=inputs, code=code_identity(), environment=environment())
    if split_reference is not None:
        request['split_reference'] = split_reference
    return digest(request)


def release_dir(root, system, name):
    require(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', name) is not None, 'invalid release name')
    require(system in SYSTEMS, 'invalid collision system')
    return Path(root).resolve() / system / name


def load_release(folder, *, request=None):
    """Verify a sealed release and every artifact hash before reuse or scoring.

    When supplied, request must also match the input/recipe/code identity.
    A directory name or an isolated model file is not a completed release.
    """
    folder = Path(folder).resolve()
    manifest = read_json(folder / 'manifest.json')
    require(manifest.get('schema') == 'PhotonIDReleaseV1' and manifest.get('status') == 'PASS',
            'model release is not validated')
    body = {k: v for k, v in manifest.items() if k != 'release_sha256'}
    require(digest(body) == manifest['release_sha256'], 'model manifest changed')
    if request is not None:
        require(manifest['request_sha256'] == request,
                'existing model has different inputs/recipe/environment; choose a new release version')
    for relative, expected in manifest['artifacts'].items():
        p = folder / relative
        require(p.resolve().is_relative_to(folder), 'artifact path escapes model release')
        require(file_hash(p) == expected, 'model artifact changed: ' + relative)
    return manifest


def list_models(root, system=None):
    """List verified generated releases and the separately pinned nominal snapshots."""
    result = []
    for collision in ([system] if system else SYSTEMS):
        for p in sorted((Path(root) / collision).glob('*/manifest.json')):
            m = read_json(p)
            if m.get('schema') == 'PhotonIDReleaseV1':
                load_release(p.parent)
            elif m.get('schema') == 'PhotonIDFrozenNominalV1':
                for relative, h in m['artifacts'].items():
                    artifact = (p.parent / relative).resolve()
                    require(artifact.is_relative_to(p.parent.resolve()), 'nominal artifact path escapes package')
                    require(file_hash(artifact) == h, 'bundled nominal payload changed')
            else:
                raise ValueError('unknown model manifest: ' + str(p))
            result.append(dict(system=collision, name=p.parent.name,
                               score_branch=m['score_branch'], status=m['status'], path=str(p.parent)))
    return result


def progress(path, phase, processed, total, **extra):
    """Atomic in-process progress; terminal receipts remain the validation authority."""
    import datetime
    value = dict(phase=phase, processed=int(processed), total=int(total),
                 remaining=int(total - processed), updated_at=datetime.datetime.now(
                     datetime.timezone.utc).isoformat(), **extra)
    import resource
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    value['peak_memory_bytes'] = int(maximum if sys.platform == 'darwin' else maximum * 1024)
    require(0 <= processed <= total, 'invalid progress count')
    if os.environ.get('PHOTONID_PROGRESS_MODE') == 'stdout':
        # Batch controllers mirror these records locally. Shared storage gets
        # one terminal record; an earlier record is retained on verified reuse.
        if phase in ('complete', 'error'):
            path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
            try:
                with path.open('x') as stream: stream.write(canonical(value) + '\n')
            except FileExistsError:
                pass
        print('PHOTONID_PROGRESS\t' + canonical(value), flush=True)
        return
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.pending')
    temp.write_text(canonical(value) + '\n')
    os.replace(temp, path)


def progress_identity(value):
    """Optional job identity from YAML, without changing a scientific recipe."""
    if value is None: return None
    fields = ('task_id', 'workstream_id', 'campaign_tag', 'executor')
    require(isinstance(value, dict) and set(value) == set(fields) and
            all(isinstance(value[k], str) and value[k].strip() for k in fields),
            'progress identity requires task_id, workstream_id, campaign_tag and executor')
    return dict(value)


class WorkProgress:
    """One in-process work counter, with a fixed unit and total across phases.

    Input inspection counts files; population preparation counts event visits
    over both passes. Selected photons are separate witnesses, never the event
    denominator. Identity can be supplied by an analyst's controlling job.
    """
    def __init__(self, path, total, unit, row_id, identity=None, *, skip_total=0):
        import socket
        import time
        self.path = Path(path)
        self.total = int(total)
        require(self.total >= 0, 'negative work total')
        self.skip_total = int(skip_total); self.skipped = 0
        require(0 <= self.skip_total <= self.total, 'invalid deterministic skip total')
        defaults = dict(task_id='PhotonID', workstream_id='preparation',
                        campaign_tag=row_id, executor='photon_id.py')
        defaults.update(progress_identity(identity) or {})
        self.identity = dict(defaults, host=socket.gethostname(), row_id=row_id)
        self.started = time.monotonic()
        self.processed = 0
        self.extra = dict(unit=unit, latest_error=None)
        self.last_write = self.started
        self.phase_started = self.started
        self.phase_name = None

    def write(self, phase, processed=None, *, skip_processed=None, **extra):
        import time
        processed = self.processed if processed is None else int(processed)
        require(self.processed <= processed <= self.total, 'work progress regressed or exceeded its total')
        skipped = self.skipped if skip_processed is None else int(skip_processed)
        require(self.skipped <= skipped <= self.skip_total and
                self.processed - self.skipped <= processed - skipped <= self.total - self.skip_total,
                'processed/skip progress regressed or exceeded its total')
        self.skipped = skipped
        self.processed = processed
        self.extra.update(extra)
        elapsed = time.monotonic() - self.started
        self.last_write = time.monotonic()
        self.phase_name = phase
        rate = processed / elapsed if elapsed else 0.
        row_rate = (processed - skipped) / elapsed if elapsed else 0.
        progress(self.path, phase, processed, self.total,
                 schema='LongRunningCanaryProgressV1', producer='PhotonID', **self.identity,
                 row_processed=processed - skipped, row_total=self.total - self.skip_total,
                 row_remaining=self.total - self.skip_total - processed + skipped,
                 deterministic_skip_processed=skipped, deterministic_skip_total=self.skip_total,
                 deterministic_skip_remaining=self.skip_total - skipped,
                 skip_processed=skipped, skip_total=self.skip_total, skip_remaining=self.skip_total - skipped,
                 rate_per_second=rate, row_rate_per_second=row_rate,
                 eta_seconds=(self.total - self.skip_total - processed + skipped) / row_rate
                     if row_rate else (0. if processed == self.total else None),
                 elapsed_seconds=elapsed, **self.extra)

    def stage(self, phase, processed, total, unit):
        """Report auxiliary work without moving the event/photon denominator."""
        import time
        require(type(processed) is int and type(total) is int and 0 <= processed <= total,
                'invalid preparation stage count')
        if phase != self.phase_name or unit != self.extra.get('phase_unit') or 'phase_total' not in self.extra:
            self.phase_started = time.monotonic()
        else:
            require(total == self.extra['phase_total'] and processed >= self.extra['phase_processed'],
                    'preparation stage count regressed or total changed')
        elapsed = time.monotonic() - self.phase_started
        rate = processed / elapsed if elapsed else 0.
        self.write(phase, phase_processed=processed, phase_total=total, phase_remaining=total-processed,
                   phase_unit=unit, phase_rate_per_second=rate,
                   phase_eta_seconds=(total-processed)/rate if rate else (0. if processed == total else None))

    def beat(self):
        """Refresh from the processing invocation during a long SQLite call."""
        import time
        if time.monotonic() - self.last_write >= 30:
            self.write(self.phase_name)

    def fail(self, error):
        self.write('error', latest_error=type(error).__name__ + ': ' + str(error), status='ERROR')
