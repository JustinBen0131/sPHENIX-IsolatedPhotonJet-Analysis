"""Diagnostics: saved model-training QA for offline FinalAnalysis rendering.

Use the prepared population, its training weights and already validated cell
scores. ROOT owns histogram accumulation and exact feature moments; Validation
owns the weighted performance definitions. No eventTree reading, fitting,
score publication or plotting occurs here. Each model's full support supplies
its distributions and global metrics; a narrower QA slice is labelled explicitly.
"""
from __future__ import annotations

from array import array
import hashlib
from pathlib import Path

import numpy as np

from . import Models, Validation


FEATURE_FIELDS = {
    'index': 'I', 'signal_entries': 'L', 'background_entries': 'L',
    'signal_mean': 'D', 'background_mean': 'D', 'signal_std': 'D', 'background_std': 'D',
    'standardized_mean_difference': 'D', 'separation_valid': 'I',
    'weighted_signal_sumw': 'D', 'weighted_background_sumw': 'D',
    'weighted_signal_sumw2': 'D', 'weighted_background_sumw2': 'D',
    'weighted_signal_mean': 'D', 'weighted_background_mean': 'D',
    'weighted_signal_std': 'D', 'weighted_background_std': 'D',
    'weighted_standardized_mean_difference': 'D', 'weighted_separation_valid': 'I',
    'total_gain': 'D', 'total_gain_fraction': 'D', 'gain_valid': 'I',
}
METRIC_FIELDS = {'rows': 'L', 'background_rows': 'L', 'signal_rows': 'L',
                 'sumw': 'D', 'auc': 'D', 'logloss': 'D', 'brier': 'D', 'valid': 'I'}


def histogram(name, values, weights, bins, low, high):
    """Accumulate finite values/positive weights in bounded ROOT FillN batches.

    No observation may leave the supplied half-open axis. ROOT stores the
    actual-value moments and Sumw2, rather than approximating moments with
    bin centers. The detached histogram remains owned by the caller.
    """
    import ROOT
    Models.require(len(values) == len(weights) and high > low, 'invalid QA histogram inputs')
    hist = ROOT.TH1D(name, '', bins, low, high)
    hist.SetDirectory(0); hist.Sumw2()
    for start in range(0, len(values), 65536):
        x = np.ascontiguousarray(values[start:start + 65536], dtype=np.float64)
        w = np.ascontiguousarray(weights[start:start + 65536], dtype=np.float64)
        Models.require(np.isfinite(x).all() and np.isfinite(w).all() and (w > 0.).all() and
                       ((x >= low) & (x < high)).all(), 'QA histogram has invalid or out-of-axis rows')
        hist.FillN(len(x), x, w)
    return hist


def moments(values, weights, name):
    """Population mean/standard deviation and weight sums from ROOT statistics."""
    Models.require(len(values) > 0 and np.isfinite(values).all(), 'feature QA needs finite class rows')
    low, high = float(np.min(values)), float(np.max(values))
    if high == low:
        low -= max(1., abs(low) * .01); high += max(1., abs(high) * .01)
    else:
        high = float(np.nextafter(high, np.inf))
    hist = histogram(name, values, weights, 64, low, high)
    stats = array('d', [0.] * 4); hist.GetStats(stats)
    return dict(entries=int(hist.GetEntries()), sumw=stats[0], sumw2=stats[1],
                mean=float(hist.GetMean()), std=float(hist.GetStdDev()))


def separation(signal, background):
    """Signed standardized mean difference; zero dispersion is unavailable."""
    scale = np.sqrt((signal['std'] ** 2 + background['std'] ** 2) / 2.)
    return float((signal['mean'] - background['mean']) / scale) if scale > 0. else None


def gain(folder, recipe, cells):
    """Sum saved XGBoost total_gain in recipe order across all routing cells.

    The fraction normalizes that sum over features. Unused BDT features have
    zero gain; non-BDT algorithms have no tree gain. Unknown feature names,
    changed feature order and escaped booster paths refuse the product.
    """
    folder = Path(folder).resolve(); features = recipe['features']
    totals = np.zeros(len(features), np.float64); by_cell = {}
    if recipe['algorithm'] != 'bdt':
        return dict(status='NOT_APPLICABLE', total=totals.tolist(), fraction=totals.tolist(), cells={})
    from xgboost import Booster
    for cell in cells:
        path = (folder / cell['artifact']).resolve()
        Models.require(path.is_relative_to(folder), 'QA model path escapes release')
        artifact = Models.read_json(path)
        Models.require(artifact['features'] == features, 'gain feature order differs')
        booster_path = (path.parent / artifact['booster_file']).resolve()
        Models.require(booster_path.is_relative_to(folder), 'QA booster path escapes release')
        booster = Booster(); booster.load_model(booster_path)
        values = np.zeros(len(features), np.float64)
        for name, value in booster.get_score(importance_type='total_gain').items():
            if name in features:
                index = features.index(name)
            else:
                Models.require(name.startswith('f') and name[1:].isdigit(), 'unknown gain feature: ' + name)
                index = int(name[1:])
            Models.require(0 <= index < len(features) and np.isfinite(value) and value >= 0., 'invalid feature gain')
            values[index] += value
        by_cell[cell['id']] = values.tolist(); totals += values
    fraction = totals / totals.sum() if totals.sum() > 0. else totals.copy()
    return dict(status='PASS', total=totals.tolist(), fraction=fraction.tolist(), cells=by_cell)


def feature_table(population, recipe, gains, *, tracker=None):
    """Save raw and weighted class moments; raw separation is Figure 41's default."""
    indices = {k: np.flatnonzero(population.y == k) for k in (0, 1)}
    rows = []
    for index, name in enumerate(recipe['features']):
        raw, weighted = {}, {}
        for klass, label in ((1, 'signal'), (0, 'background')):
            selected = indices[klass]
            values = population.x[selected, index]
            raw[label] = moments(values, np.ones(len(selected)), f'feature_{index}_{label}_raw')
            weighted[label] = moments(values, population.weights[selected], f'feature_{index}_{label}_weighted')
        difference = separation(raw['signal'], raw['background'])
        weighted_difference = separation(weighted['signal'], weighted['background'])
        row = dict(index=index, feature=name, standardized_mean_difference=difference,
                   separation_valid=int(difference is not None),
                   weighted_standardized_mean_difference=weighted_difference,
                   weighted_separation_valid=int(weighted_difference is not None),
                   total_gain=gains['total'][index], total_gain_fraction=gains['fraction'][index],
                   gain_valid=int(gains['status'] == 'PASS'))
        for label in ('signal', 'background'):
            row.update({label + '_' + key: raw[label][key] for key in ('entries', 'mean', 'std')})
            row.update({'weighted_' + label + '_' + key: weighted[label][key]
                        for key in ('sumw', 'sumw2', 'mean', 'std')})
        rows.append(row)
        if tracker: tracker.stage('qa-feature-moments', index + 1, len(recipe['features']), 'features')
    return rows


def performance(population, scores, selected):
    """Canonical weighted metrics; empty/single-class bins remain explicitly invalid."""
    indices = np.flatnonzero(selected); y = population.y[indices]; w = population.weights[indices]
    counts = [int((y == k).sum()) for k in (0, 1)]
    row = dict(rows=len(indices), background_rows=counts[0], signal_rows=counts[1],
               sumw=float(w.sum()), valid=int(all(counts)), auc=None, logloss=None, brier=None)
    if row['valid']:
        metrics = Validation.metrics(y, scores[indices], w)
        row.update({key: metrics[key] for key in ('auc', 'logloss', 'brier')})
    return row


def grid_edges(recipe, config):
    """Extend the canonical QA grid to model support and wide guard boundaries.

    These are diagnostic bins, not a declaration of analysis/unfolding bins.
    Nominal routing comparisons retain their existing eight ET bins.
    """
    low, high = recipe['domain']['et_min'], recipe['domain']['et_max_exclusive']
    edges = [low, high] + [e for e in config['bins']['et'] if low < e < high]
    if low != config['bins']['et'][0] or high != config['bins']['et'][-1]:
        edges += recipe['et_ranges_gev']['analysis'] or []
    return sorted(set(edges))


def holdout_grid(population, scores, recipe, config):
    """Full-support ET × C3 holdout QA, with explicit unavailable-class bins."""
    rows = []
    edges = grid_edges(recipe, config)
    for et_low, et_high in zip(edges[:-1], edges[1:]):
        centralities = list(zip(config['bins']['centrality_coarse'][:-1], config['bins']['centrality_coarse'][1:])) \
            if recipe['system'] == 'AuAu' else [(0., 80.)]
        for cent_low, cent_high in centralities:
            mask = population.holdout & (population.et >= et_low) & (population.et < et_high)
            if recipe['system'] == 'AuAu':
                mask &= (population.centrality >= cent_low) & (population.centrality < cent_high)
            rows.append(dict(et_low=et_low, et_high=et_high, centrality_low=cent_low,
                             centrality_high=cent_high, **performance(population, scores, mask)))
    return rows


def write_table(root, name, fields, rows, *, tracker=None, total=None):
    """Write scalar ROOT columns; unavailable float metrics use NaN plus a validity flag."""
    import ROOT
    root.cd(); tree = ROOT.TTree(name, name)
    buffers = {key: array({'D': 'd', 'I': 'i', 'L': 'q'}[kind], [0]) for key, kind in fields.items()}
    for key, kind in fields.items(): tree.Branch(key, buffers[key], key + '/' + kind)
    if tracker: tracker.stage('qa-table-' + name, 0, total, 'rows')
    count = 0
    for row in rows:
        for key, kind in fields.items():
            value = row[key]
            buffers[key][0] = np.nan if value is None else float(value) if kind == 'D' else int(value)
        tree.Fill()
        count += 1
        if tracker and count % 65536 == 0: tracker.stage('qa-table-' + name, count, total, 'rows')
    if tracker:
        Models.require(count == total, 'QA table row count changed')
        tracker.stage('qa-table-' + name, count, total, 'rows')
    tree.Write(); tree.SetDirectory(0)


def score_diagnostics(root, population, scores, selected, directory):
    """Class score densities and weighted ROC from the same labelled holdout."""
    from sklearn.metrics import roc_curve
    folder = root.mkdir(directory); rows = np.flatnonzero(selected)
    report = performance(population, scores, selected)
    if not report['valid']: return report
    for klass, label in ((1, 'signal'), (0, 'background')):
        indices = rows[population.y[rows] == klass]
        for weighting, weights in [('raw', np.ones(len(indices))), ('weighted', population.weights[indices])]:
            h = histogram(label + '_' + weighting, scores[indices], weights, 50, 0., float(np.nextafter(1., np.inf)))
            folder.cd(); h.Write()
    fpr, tpr, thresholds = roc_curve(population.y[rows], scores[rows], sample_weight=population.weights[rows], drop_intermediate=False)
    write_table(folder, 'roc', dict(false_positive_rate='D', true_positive_rate='D', threshold='D'),
                (dict(false_positive_rate=f, true_positive_rate=t, threshold=s) for f, t, s in zip(fpr, tpr, thresholds)))
    return report


def original_comparison(root, folder, population, recipe, scores, tracker):
    """Pair frozen PPG12 and retrained predictions on exactly the new test rows.

    Both scores use the same H70 float32 vector. No stored producer sentinel or
    inference-domain gate is substituted for direct evaluation of the frozen
    reference model. Rows retain their original source-entry/photon identity.
    """
    from . import Training
    reference = Models.HERE / 'models/PP/nominal/model.tmva.root'
    Models.require(Models.file_hash(reference) == recipe['reference_model_sha256'], 'frozen PPG12 model bytes changed')
    native = Training.native_evaluator(reference, dict(type='ROOT_RBDT', file=reference.name, key='myBDT', inputs=11))
    rows = np.flatnonzero(population.holdout); originals = np.empty(len(rows), np.float64)
    tracker.stage('qa-original-predictions', 0, len(rows), 'photons')
    for begin in range(0, len(rows), 65536):
        part = rows[begin:begin + 65536]; originals[begin:begin + len(part)] = native(population.x[part])
        tracker.stage('qa-original-predictions', begin + len(part), len(rows), 'photons')
    Models.require(np.isfinite(originals).all() and ((originals >= 0.) & (originals <= 1.)).all(), 'invalid frozen PPG12 scores')
    files = population.receipt['array_files']; work = Path(population.x.filename).parent
    identities = {key: np.load(work / files[key], mmap_mode='r') for key in ('source_index', 'entry', 'photon')}
    comparison = root.mkdir('original_comparison')
    fields = dict(source_index='I', entry='L', photon='I', label='I', et='D', weight='D', original_score='D', retrained_score='D', difference='D')
    def paired_rows():
        for index, row in enumerate(rows):
            yield dict(source_index=identities['source_index'][row], entry=identities['entry'][row], photon=identities['photon'][row],
                       label=population.y[row], et=population.et[row], weight=population.weights[row],
                       original_score=originals[index], retrained_score=scores[row], difference=scores[row]-originals[index])
    write_table(comparison, 'paired_scores', fields, paired_rows(), tracker=tracker, total=len(rows))
    for klass, label in ((1, 'signal'), (0, 'background')):
        selected = population.y[rows] == klass; weights = population.weights[rows[selected]]
        delta = scores[rows[selected]] - originals[selected]
        h = histogram(label + '_difference', delta, weights, 100, -1., float(np.nextafter(1., np.inf)))
        comparison.cd(); h.Write()
    return dict(reference_model_sha256=recipe['reference_model_sha256'], reference='bundled frozen PPG12 RBDT',
                rows=len(rows), split='test only; validation excluded', score_difference='retrained minus original',
                source_roster=population.receipt['sources'], paired_scores_fields=fields,
                original_metrics=Validation.metrics(population.y[rows], originals, population.weights[rows]),
                retrained_metrics=Validation.metrics(population.y[rows], scores[rows], population.weights[rows]))


def save(folder, population, recipe, scores, cells, config, request_sha256):
    """Seal ROOT/JSON QA from the full population and already validated scores.

    Existing matching products are hash-verified and reused. An interrupted
    ROOT product is retained and reported, never silently overwritten.
    The release manifest seals both output files after this function succeeds.
    """
    import ROOT
    folder = Path(folder).resolve(); path = folder / 'diagnostics.root'; json_path = folder / 'diagnostics.json'
    scores = np.asarray(scores, np.float64)
    n = len(population.y)
    Models.require(population.x.shape == (n, len(recipe['features'])) and
                   all(np.shape(getattr(population, key)) == (n,) for key in ('et', 'centrality', 'weights', 'train')) and
                   np.isin(population.y, [0, 1]).all() and np.isfinite(population.weights).all() and
                   (population.weights > 0.).all(), 'QA population arrays or training weights are invalid')
    Models.require(scores.shape == population.y.shape and np.isfinite(scores).all() and
                   ((scores >= 0.) & (scores <= 1.)).all(), 'QA scores must cover every population row exactly once')
    Models.require(all(key in population.receipt for key in
                       ('input_recipe_sha256', 'split_sha256', 'event_split_sha256', 'sources')),
                   'QA population provenance is incomplete')
    identity = dict(request_sha256=request_sha256, recipe_sha256=Models.digest(recipe),
                    population_sha256=population.receipt['input_recipe_sha256'],
                    split_sha256=population.receipt['split_sha256'],
                    event_split_sha256=population.receipt['event_split_sha256'],
                    population_receipt_sha256=Models.digest(population.receipt),
                    training_weights_sha256=hashlib.sha256(np.ascontiguousarray(population.weights, dtype=np.float64).view(np.uint8)).hexdigest(),
                    score_sha256=hashlib.sha256(np.ascontiguousarray(scores).view(np.uint8)).hexdigest())
    model_artifacts = {}
    for cell in cells:
        artifact_path = (folder / cell['artifact']).resolve()
        Models.require(artifact_path.is_relative_to(folder.resolve()), 'QA artifact path escapes release')
        artifact = Models.read_json(artifact_path)
        paths = [artifact_path]
        if recipe['algorithm'] == 'bdt':
            paths += [artifact_path.parent / artifact['booster_file'],
                      artifact_path.parent / artifact['native_evaluator']['file']]
        for product in paths:
            Models.require(product.resolve().is_relative_to(folder.resolve()), 'QA model product escapes release')
            model_artifacts[str(product.relative_to(folder.resolve()))] = Models.file_hash(product)
    identity['model_artifacts_sha256'] = Models.digest(model_artifacts)
    if json_path.exists():
        saved = Models.read_json(json_path)
        Models.require(saved['identity'] == identity, 'existing model QA binds different scores or population')
        Models.require(Models.file_hash(path) == saved['root_sha256'], 'saved model QA ROOT artifact changed')
        return saved
    Models.require(not path.exists(), 'interrupted QA ROOT product retained: ' + str(path))
    tracker = Models.WorkProgress(folder / 'progress.json', 1, 'QA products', request_sha256)
    tracker.write('qa-feature-moments', input_photons=n)
    gains = gain(folder, recipe, cells)
    features = feature_table(population, recipe, gains, tracker=tracker)
    metrics = [dict(split=index, **performance(population, scores, mask))
               for index, mask in enumerate((population.train, population.holdout))]
    grid = holdout_grid(population, scores, recipe, config)
    Models.require(sum(row['rows'] for row in grid) == int(population.holdout.sum()),
                   'primary QA grid does not cover the full held-out support')
    metadata = dict(schema='PhotonIDDiagnosticsV1', status='PASS', name=recipe['name'], system=recipe['system'],
                    identity=identity, sources=population.receipt['sources'], code_sha256=Models.code_identity(),
                    model_artifacts=model_artifacts, environment=Models.environment(), ROOT=ROOT.gROOT.GetVersion(),
                    et_ranges_gev=recipe['et_ranges_gev'], features=recipe['features'],
                    feature_order=recipe['features'], max_depth=recipe['learner'].get('max_depth'),
                    weight_definition='class balance times per-class eta/ET flattening; no analysis event weights',
                    score_origin='validated saved ROOT_RBDT cell scores' if recipe['algorithm'] == 'bdt' else 'validated saved model scores',
                    feature_separation='(signal_mean-background_mean)/sqrt((signal_variance+background_variance)/2), population variances',
                    gain_definition='sum of saved XGBoost total_gain across cells; fraction normalized over all features',
                    global_metric_scope='full model application support',
                    grid_scope='full model support; model-specific training weights; diagnostic bins, not analysis bins',
                    grid_et_edges=grid_edges(recipe, config), grid_centrality_edges=config['bins']['centrality_coarse'] if recipe['system'] == 'AuAu' else [],
                    split_method=population.receipt.get('split_method', 'event_grouped_90_10'),
                    qualification='saved training QA; downstream physics acceptance is separate')
    root = ROOT.TFile.Open(str(path), 'CREATE')
    Models.require(root and not root.IsZombie(), 'cannot create model QA ROOT output')
    try:
        training = root.mkdir('training')
        for klass, label in ((1, 'signal'), (0, 'background')):
            selected = np.flatnonzero(population.y == klass)
            variables = [('et', population.et[selected], recipe['domain']['et_min'], recipe['domain']['et_max_exclusive']),
                         ('eta', population.x[selected, recipe['features'].index('cluster_Eta')], -.7, .7)]
            for variable, values, low, high in variables:
                for weighting, weights in [('raw', np.ones(len(selected))), ('weighted', population.weights[selected])]:
                    name = f'{variable}_{label}_{weighting}'
                    hist = histogram(name, values, weights, recipe['weights']['bins'], low, high)
                    training.cd(); hist.Write()
        write_table(root, 'features', FEATURE_FIELDS, features)
        write_table(root, 'performance', dict(split='I', **METRIC_FIELDS), metrics)
        fields = dict(et_low='D', et_high='D', centrality_low='D', centrality_high='D', **METRIC_FIELDS)
        write_table(root, 'holdout_et_c3', fields, grid)
        metadata['heldout_scores'] = score_diagnostics(root, population, scores, population.holdout, 'holdout_scores')
        if recipe['system'] == 'AuAu':
            peripheral = population.holdout & (population.centrality >= 50.) & (population.centrality < 80.) & (population.et >= 15.) & (population.et < 35.)
            if all(np.any(peripheral & (population.y == k)) for k in (0, 1)):
                metadata['peripheral_15to35_scores'] = score_diagnostics(root, population, scores, peripheral, 'peripheral_15to35_scores')
        if recipe['name'].startswith('ppg12_retrained_'):
            metadata['original_comparison'] = original_comparison(root, folder, population, recipe, scores, tracker)
        root.mkdir('meta').cd(); ROOT.TObjString(Models.canonical(metadata)).Write('diagnostics')
        root.Write()
    finally:
        root.Close()
    report = dict(metadata, root_file=path.name, root_sha256=Models.file_hash(path),
                  feature_rows=features, gain=gains, performance=metrics, holdout_et_c3=grid,
                  root_tables=dict(features=FEATURE_FIELDS, performance=dict(split='I', **METRIC_FIELDS),
                                   holdout_et_c3=fields))
    Models.save_json(json_path, report)
    tracker.write('complete', 1, status='PASS')
    return report


def collection(root, releases, *, system='AuAu', auau_reference=None, check_only=False):
    """Create an immutable, portable YAML index of sealed model releases.

    FinalAnalysis reads this index for its selector and saved routing/depth QA.
    Comparison groups require identical source/population/event-split/weight
    pins. The wider model is indexed with its own support and calibration; it
    is never admitted to a same-population routing or depth comparison.
    """
    import os
    import yaml
    root = Path(root).resolve(); entries = []; recipes = {}
    for release in releases:
        folder = Path(release).resolve()
        Models.require(folder.is_relative_to(root), 'QA collection release escapes models root')
        manifest = Models.load_release(folder)
        Models.require(manifest['system'] == system and manifest['recipe']['algorithm'] == 'bdt',
                       'QA collection collision system or learner differs')
        diagnostics = Models.read_json(folder / 'diagnostics.json')
        Models.require(Models.file_hash(folder / 'diagnostics.root') == diagnostics['root_sha256'] and
                       diagnostics['identity']['request_sha256'] == manifest['request_sha256'] and
                       diagnostics['identity']['recipe_sha256'] == Models.digest(manifest['recipe']),
                       'model QA is not bound to its sealed release')
        pins = {key: diagnostics['identity'][key] for key in
                ('population_sha256', 'split_sha256', 'event_split_sha256', 'training_weights_sha256')}
        pins['sources_sha256'] = Models.digest(diagnostics['sources'])
        pins['implementation_sha256'] = Models.digest(dict(code=manifest['code_sha256'],
                                                         environment=manifest['environment'], ROOT=diagnostics['ROOT']))
        ranges = manifest['recipe']['et_ranges_gev']
        Models.require(diagnostics['et_ranges_gev'] == ranges, 'QA collection support differs from model')
        selection = Models.read_json(folder / 'selection.json')
        Models.require(selection['et_ranges_gev'] == ranges and selection['domain'] == manifest['recipe']['domain'],
                       'selection range differs from model QA')
        entries.append(dict(name=manifest['name'], system=system, path=str(folder.relative_to(root)),
                            release_sha256=manifest['release_sha256'], et_ranges_gev=ranges, pins=pins,
                            diagnostics_json_sha256=Models.file_hash(folder / 'diagnostics.json'),
                            diagnostics_root_sha256=diagnostics['root_sha256'],
                            axes=manifest['recipe']['axes'], max_depth=manifest['recipe']['learner']['max_depth']))
        recipes[manifest['name']] = manifest['recipe']
    Models.require(entries and len({entry['name'] for entry in entries}) == len(entries), 'empty or duplicate model collection')
    nominal = [entry for entry in entries if entry['name'].startswith('nominal_retrained_')]
    Models.require(len(nominal) == 1, 'QA collection needs exactly one nominal reference')
    baseline = nominal[0]; routed = []; depth = []
    for entry in entries:
        Models.require(all(entry['pins'][key] == baseline['pins'][key] for key in
                           ('sources_sha256', 'implementation_sha256')),
                       'model QA collection has different input sources or implementations')
        if entry['name'].startswith('wideET_'):
            ignored = ('name', 'score_branch', 'domain', 'training_domain', 'et_ranges_gev')
            reference = {key: value for key, value in recipes[baseline['name']].items() if key not in ignored}
            variant = {key: value for key, value in recipes[entry['name']].items() if key not in ignored}
            shared_domain = {key: value for key, value in recipes[baseline['name']]['domain'].items()
                             if key not in ('et_min', 'et_max_exclusive')}
            wide_domain = {key: value for key, value in recipes[entry['name']]['domain'].items()
                           if key not in ('et_min', 'et_max_exclusive')}
            Models.require(variant == reference and wide_domain == shared_domain,
                           'wider-range model changes additional training settings')
        if entry['axes'] or entry['name'].startswith('depth') or entry is baseline:
            Models.require(entry['pins'] == baseline['pins'] and entry['et_ranges_gev'] == baseline['et_ranges_gev'],
                           'routing/depth QA comparison has different sources, population, split or weights')
            reference = {key: value for key, value in recipes[baseline['name']].items()
                         if key not in ('name', 'score_branch', 'axes', 'learner')}
            variant = {key: value for key, value in recipes[entry['name']].items()
                       if key not in ('name', 'score_branch', 'axes', 'learner')}
            learner, expected = dict(recipes[entry['name']]['learner']), dict(recipes[baseline['name']]['learner'])
            if not entry['axes']:
                learner.pop('max_depth'); expected.pop('max_depth')
            Models.require(variant == reference and learner == expected,
                           'routing/depth comparison changes additional training settings')
            (routed if entry['axes'] else depth).append(entry['name'])
    body = dict(schema='PhotonIDModelQACollectionV1', system=system, default_model=baseline['name'], models=entries,
                comparisons=dict(routing=dict(reference=baseline['name'], models=routed),
                                 depth=dict(reference=baseline['name'], models=depth)),
                qualification='sealed numerical model releases and saved QA; physics acceptance separate')
    if system == 'PP':
        ppg12 = [entry for entry in entries if entry['name'].startswith('ppg12_retrained_')]
        Models.require(len(ppg12) == 1, 'pp campaign needs exactly one PPG12-equivalent model')
        body['comparisons'] = dict(ppg12=dict(model=ppg12[0]['name'], reference='bundled frozen PPG12',
                                            table='original_comparison/paired_scores'),
                                   auau_settings=dict(pp_model=baseline['name'], et_range_gev=[15., 35.],
                                                      auau_centrality_range=[50., 80.], status='UNBOUND'))
        if auau_reference is not None:
            reference_path = Path(auau_reference).resolve()
            Models.require(reference_path.is_relative_to(root), 'AuAu reference escapes models root')
            reference = Models.load_release(reference_path)
            Models.require(reference['system'] == 'AuAu' and reference['name'].startswith('nominal_retrained_'),
                           'pp comparison requires the canonical AuAu release')
            aa_recipe = reference['recipe']; pp_recipe = recipes[baseline['name']]
            Models.require(aa_recipe['learner'] == pp_recipe['learner'] and
                           aa_recipe['training'] == pp_recipe['training'] and aa_recipe['weights'] == pp_recipe['weights'] and
                           [f for f in aa_recipe['features'] if f not in ('cluster_weta33_cogx', 'cluster_wphi33_cogx', 'centrality')] == pp_recipe['features'] and
                           [aa_recipe['domain'][k] for k in ('et_min', 'et_max_exclusive', 'eta_abs_max', 'vertex_abs_max_cm')] ==
                           [pp_recipe['domain'][k] for k in ('et_min', 'et_max_exclusive', 'eta_abs_max', 'vertex_abs_max_cm')],
                           'AuAu/pp comparison changes the canonical training settings')
            qa = Models.read_json(reference_path / 'diagnostics.json')
            Models.require('peripheral_15to35_scores' in qa, 'AuAu reference has no saved peripheral heldout comparison')
            body['comparisons']['auau_settings'].update(status='BOUND', auau_reference=dict(
                path=str(reference_path.relative_to(root)), name=reference['name'], release_sha256=reference['release_sha256'],
                diagnostics_json_sha256=Models.file_hash(reference_path / 'diagnostics.json'), diagnostics_root_sha256=qa['root_sha256']))
    pin = Models.digest(body); folder = root / '.campaigns' / pin
    path = folder / 'models.yaml'
    # Paths resolve against the YAML file, so copying the models root preserves the collection.
    for entry in body['models']:
        entry['path'] = os.path.relpath(root / entry['path'], folder)
    if system == 'PP' and auau_reference is not None:
        body['comparisons']['auau_settings']['auau_reference']['path'] = os.path.relpath(reference_path, folder)
    body['collection_sha256'] = Models.digest(body)
    raw = yaml.safe_dump(body, sort_keys=True)
    if check_only: return path
    folder.mkdir(parents=True, exist_ok=True)
    try:
        with path.open('x') as stream: stream.write(raw)
    except FileExistsError:
        Models.require(path.read_text() == raw, 'immutable model QA collection differs')
    return path
