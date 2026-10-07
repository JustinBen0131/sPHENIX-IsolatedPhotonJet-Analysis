"""Training: one population and split, with three small fitting backends.

BDT uses the established XGBoost settings and exports ROOT RBDT for scoring
and working-point calibration. The saved XGBoost model remains its checked
training reference. NumPy MLP and linear sigmoid backends retain the existing
tanh/Adam/weighted-loss numerical conventions.
Fixed epochs avoid selecting a model against the outer holdout. Training
returns portable artifacts; saved-model and working-point checks live elsewhere.
"""
from __future__ import annotations

import itertools
from pathlib import Path
import time

import numpy as np

from .Models import progress, require, save_json


def runtime_status(algorithm):
    """Read-only import/API preflight; never install or switch environments.

    ROOT is needed for saved QA and native BDT export. API availability does
    not prove that another ROOT version can read this runtime's artifacts.
    """
    import importlib
    missing = []; versions = {}
    modules = ('ROOT', 'xgboost', 'sklearn') if algorithm == 'bdt' else ('ROOT',)
    for name in modules:
        try:
            module = importlib.import_module(name)
            versions[name] = module.gROOT.GetVersion() if name == 'ROOT' else module.__version__
            if name == 'xgboost' and module.__version__ != '3.2.0':
                missing.append('reference BDT fitting requires xgboost==3.2.0; observed ' + module.__version__)
            if name == 'ROOT':
                # ROOT can import despite a failed Cling/bootstrap environment.
                # Prove basic C++ bindings before opening large input groups.
                vector = module.std.vector('float')()
                if int(vector.size()) != 0: missing.append('ROOT C++ binding probe failed')
            if name == 'ROOT' and algorithm == 'bdt':
                namespace = getattr(getattr(module, 'TMVA', None), 'Experimental', None)
                for api in ('SaveXGBoost', 'RBDT'):
                    if not hasattr(namespace, api): missing.append('ROOT.TMVA.Experimental.' + api)
        except (ImportError, OSError, AttributeError, RuntimeError) as exc:
            missing.append(name + ': ' + str(exc))
    return dict(schema='PhotonIDRuntimeCheckV1', status='PASS' if not missing else 'BLOCKED_DEPENDENCY',
                algorithm=algorithm, versions=versions, missing=missing,
                qualification='imports/export APIs only; cross-runtime artifact compatibility requires validation')


def cells(recipe):
    """Yield the Cartesian routing cells; an empty axis map gives one global model."""
    axes = list(recipe['axes'])
    bounds = [list(zip(recipe['axes'][a][:-1], recipe['axes'][a][1:])) for a in axes]
    for index, values in enumerate(itertools.product(*bounds) if axes else [()]):
        yield dict(id=f'cell_{index:03d}', bounds={a: list(v) for a, v in zip(axes, values)})


def cell_mask(population, cell):
    """Select a cell using half-open ET (GeV) and centrality (percent) intervals."""
    mask = np.ones(len(population.y), dtype=bool)
    for axis, (lo, hi) in cell['bounds'].items():
        values = population.et if axis == 'et' else population.centrality
        mask &= (values >= lo) & (values < hi)
    return mask


def sigmoid(values):
    values = np.asarray(values, np.float64)
    out = np.empty_like(values)
    positive = values >= 0
    out[positive] = 1. / (1. + np.exp(-values[positive]))
    exponential = np.exp(values[~positive])
    out[~positive] = exponential / (1. + exponential)
    return out


def numpy_predict(artifact, x):
    """Evaluate saved preprocessing and layers in recipe feature order.

    Artifact weights are stored output-by-input; the evaluation transposes
    them for row-wise matrix multiplication and returns one sigmoid per row.
    """
    preprocessing = artifact['preprocessing']
    values = (np.asarray(x, np.float64) - preprocessing['mean']) / preprocessing['scale']
    if preprocessing['clip'] > 0:
        values = np.clip(values, -preprocessing['clip'], preprocessing['clip'])
    for index, layer in enumerate(artifact['layers']):
        values = values @ np.asarray(layer['weights'], np.float64).T + layer['bias']
        if index < len(artifact['layers']) - 1:
            values = np.tanh(values)
    return sigmoid(values[:, 0])


def preprocessing(population, indices, clip):
    """Fit scaling on the selected training rows, in bounded batches."""
    sums = np.zeros(population.x.shape[1]); squares = sums.copy()
    for start in range(0, len(indices), 65536):
        x = np.asarray(population.x[indices[start:start + 65536]], np.float64)
        sums += x.sum(axis=0); squares += (x * x).sum(axis=0)
    mean = sums / len(indices)
    variance = np.maximum(squares / len(indices) - mean * mean, 0.)
    scale = np.sqrt(variance)
    scale[scale <= 1e-9] = 1.
    return dict(mean=mean.tolist(), scale=scale.tolist(), clip=float(clip))


def fit_numpy(population, mask, recipe, output, cell_id):
    """Fit one MLP/logistic cell for fixed epochs and save its portable artifact.

    Only mask AND population.train enters fitting or scaling. The outer
    holdout is checked later by Validation, without choosing an epoch here.
    """
    settings = recipe['learner']; seed = recipe['training']['seed']
    indices = np.flatnonzero(mask & population.train)
    transform = preprocessing(population, indices, settings['input_clip'])
    hidden = settings.get('hidden_layers', [])
    sizes = [population.x.shape[1]] + hidden + [1]
    rng = np.random.default_rng(seed)
    layers = []
    for n_in, n_out in zip(sizes[:-1], sizes[1:]):
        if recipe['algorithm'] == 'logreg':
            weights = rng.normal(0., .01, (n_in, n_out))
        else:
            limit = np.sqrt(6. / (n_in + n_out))
            weights = rng.uniform(-limit, limit, (n_in, n_out))
        layers.append(dict(w=weights, b=np.zeros(n_out)))
    moments = [{k: np.zeros_like(layer[k]) for k in ('w', 'b')} for layer in layers]
    squares = [{k: np.zeros_like(layer[k]) for k in ('w', 'b')} for layer in layers]
    mean_weight = float(np.mean(population.weights[indices])); step = 0
    epochs = int(settings['epochs'])
    batches = (len(indices) + int(settings['batch_size']) - 1) // int(settings['batch_size'])
    started = time.monotonic(); last_update = started
    progress(Path(output).parent / 'progress.json', 'fit-' + cell_id, 0, epochs * batches,
             learner=recipe['algorithm'], training_photons=len(indices), latest_error=None)
    for epoch in range(epochs):
        rng.shuffle(indices)
        for start in range(0, len(indices), int(settings['batch_size'])):
            idx = indices[start:start + int(settings['batch_size'])]
            x = (population.x[idx].astype(np.float64) - transform['mean']) / transform['scale']
            if transform['clip'] > 0:
                x = np.clip(x, -transform['clip'], transform['clip'])
            activations = [x]
            for number, layer in enumerate(layers):
                z = activations[-1] @ layer['w'] + layer['b']
                activations.append(z if number == len(layers) - 1 else np.tanh(z))
            score = sigmoid(activations[-1][:, 0])
            # Global mean-weight normalization gives an unbiased stochastic
            # gradient of the full weighted objective, including short batches.
            dz = ((score - population.y[idx]) * population.weights[idx] /
                  (len(idx) * mean_weight))[:, None]
            gradients = [None] * len(layers)
            for number in reversed(range(len(layers))):
                gradients[number] = dict(w=activations[number].T @ dz + settings['l2'] * layers[number]['w'],
                                         b=dz.sum(axis=0))
                if number:
                    dz = (dz @ layers[number]['w'].T) * (1. - activations[number] ** 2)
            step += 1
            for layer, gradient, moment, square in zip(layers, gradients, moments, squares):
                for name in ('w', 'b'):
                    moment[name] *= .9; moment[name] += .1 * gradient[name]
                    square[name] *= .999; square[name] += .001 * gradient[name] ** 2
                    layer[name] -= settings['learning_rate'] * (moment[name] / (1. - .9 ** step)) / (
                        np.sqrt(square[name] / (1. - .999 ** step)) + 1e-8)
            require(all(np.isfinite(layer['w']).all() and np.isfinite(layer['b']).all() for layer in layers),
                    'nonfinite learner parameters')
            now = time.monotonic()
            if now - last_update >= 30. or start + len(idx) == len(indices):
                elapsed = now - started; rate = step / elapsed if elapsed else 0.
                progress(Path(output).parent / 'progress.json', 'fit-' + cell_id, step, epochs * batches,
                         learner=recipe['algorithm'], training_photons=len(indices), epoch=epoch + 1,
                         unit='optimizer batches', elapsed_seconds=elapsed, rate=rate,
                         eta_seconds=(epochs * batches - step) / rate if rate else None, latest_error=None)
                last_update = now
    artifact = dict(schema='PhotonIDNumpyModelV1', algorithm=recipe['algorithm'],
                    features=recipe['features'], preprocessing=transform,
                    hidden_activation='tanh', output='sigmoid',
                    layers=[dict(weights=p['w'].T.tolist(), bias=p['b'].tolist()) for p in layers])
    save_json(output, artifact)
    return artifact, lambda x: numpy_predict(artifact, x)


def fit_bdt(population, mask, recipe, output, cell_id):
    """Fit on training rows, export ROOT RBDT and check all applicable predictions.

    Only the exported native evaluator supplies downstream scores and working
    points. The XGBoost artifact is kept for training/export provenance.
    """
    try:
        from xgboost import XGBClassifier
        from xgboost.callback import TrainingCallback
    except ImportError as exc:
        raise ValueError('XGBoost is unavailable in the selected environment; no packages were changed') from exc
    indices = np.flatnonzero(mask & population.train)
    holdout = np.flatnonzero(mask & population.holdout)
    if recipe['training'].get('split_method') == 'row_stratified':
        require(population.split_rank is not None, 'historical training row order is missing')
        indices = indices[np.argsort(population.split_rank[indices], kind='stable')]
    started = time.monotonic(); total = recipe['learner']['n_estimators']
    class ProgressCallback(TrainingCallback):
        def after_iteration(self, model, epoch, evals_log):
            elapsed = time.monotonic() - started; rate = (epoch + 1) / elapsed if elapsed else 0.
            progress(Path(output).parent / 'progress.json', 'fit-' + cell_id, epoch + 1, total,
                     unit='boosting rounds', elapsed_seconds=elapsed, rate=rate,
                     eta_seconds=(total - epoch - 1) / rate if rate else None, latest_error=None)
            return False
    import inspect
    fit_callbacks = 'callbacks' in inspect.signature(XGBClassifier.fit).parameters
    model = XGBClassifier(**recipe['learner'], random_state=recipe['training']['seed'],
                          objective='binary:logistic', eval_metric=['auc', 'logloss'],
                          **({} if fit_callbacks else {'callbacks':[ProgressCallback()]}))
    progress(Path(output).parent / 'progress.json', 'fit-' + cell_id, 0, recipe['learner']['n_estimators'])
    model.fit(population.x[indices], population.y[indices], sample_weight=population.weights[indices],
              eval_set=[(population.x[indices], population.y[indices]),
                        (population.x[holdout], population.y[holdout])],
              sample_weight_eval_set=[population.weights[indices], population.weights[holdout]], verbose=False,
              **({'callbacks':[ProgressCallback()]} if fit_callbacks else {}))
    booster_path = Path(output).with_suffix('.xgboost.json')
    model.save_model(booster_path)
    history_path = Path(output).with_suffix('.history.json')
    save_json(history_path, dict(schema='PhotonIDTrainingHistoryV1',
              populations=['train', 'holdout'], metrics=model.evals_result(),
              rounds=total, early_stopping=False))
    import ROOT
    ROOT.gROOT.SetBatch(True)
    require(hasattr(ROOT.TMVA.Experimental, 'SaveXGBoost'), 'ROOT does not provide the XGBoost-to-RBDT exporter')
    native_path = Path(output).with_suffix('.tmva.root')
    require(not native_path.exists(), 'native model artifact already exists')
    export_native_bdt(model, native_path, len(recipe['features']))
    specification = dict(type='ROOT_RBDT', file=native_path.name, key='myBDT', inputs=len(recipe['features']))
    native = native_evaluator(native_path, specification)
    from .Validation import native_export
    reference = lambda x: model.predict_proba(x)[:, 1].astype(np.float64)
    agreement = native_export(population, mask, reference, native, recipe)
    agreement.update(ROOT=ROOT.gROOT.GetVersion())
    artifact = dict(schema='PhotonIDXGBoostModelV1', algorithm='bdt', features=recipe['features'],
                    booster_file=booster_path.name, native_evaluator=specification,
                    training_history=history_path.name,
                    native_export_validation=agreement)
    save_json(output, artifact)
    progress(Path(output).parent / 'progress.json', 'fit-' + cell_id,
             recipe['learner']['n_estimators'], recipe['learner']['n_estimators'])
    return artifact, native


def export_native_bdt(model, path, inputs):
    """Use the established binary XGBoost3-to-ROOT6.32 export adapter.

    ROOT expects scalar base_score syntax and numbered feature names. Only
    its Python configuration view changes during export; saved booster bytes,
    fitted parameters and subsequent XGBoost predictions remain unchanged.
    """
    import re
    import ROOT
    booster = model.get_booster()
    original_config, original_names = booster.save_config, booster.feature_names
    try:
        booster.feature_names = [f'f{i}' for i in range(inputs)]
        booster.save_config = lambda: re.sub(r'"base_score":"\[([0-9eE+\-.]+)\]"',
                                             r'"base_score":"\1"', original_config())
        ROOT.TMVA.Experimental.SaveXGBoost(model, 'myBDT', str(path), num_inputs=inputs)
    finally:
        booster.save_config = original_config
        booster.feature_names = original_names


def fit(population, mask, recipe, output, cell_id):
    """Create one cell artifact; retain failure progress instead of overwriting a fit."""
    require(not Path(output).exists(), 'model artifact already exists; validate/reuse its saved release')
    try:
        return (fit_bdt if recipe['algorithm'] == 'bdt' else fit_numpy)(population, mask, recipe, output, cell_id)
    except Exception as exc:
        from .Models import read_json
        path = Path(output).parent / 'progress.json'
        previous = read_json(path) if path.exists() else {'processed':0,'total':0}
        progress(path, 'blocked-' + cell_id, previous['processed'], previous['total'], latest_error=str(exc))
        raise


def native_evaluator(path, specification):
    """Load one ROOT RBDT cell; accept finite float32 rows in released feature order."""
    require(specification['type'] == 'ROOT_RBDT', 'unsupported native evaluator')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    model = ROOT.TMVA.Experimental.RBDT(specification['key'], str(path))
    count = specification['inputs']
    def predict(x):
        values = np.ascontiguousarray(x, dtype=np.float32)
        require(values.ndim == 2 and values.shape[1] == count and np.isfinite(values).all(),
                'invalid native feature matrix')
        if not len(values):
            return np.empty(0, np.float64)
        result = np.asarray(model.Compute(values), np.float64)
        require(result.shape == (len(values), 1), 'native score shape differs')
        return result[:, 0]
    return predict


def evaluator(path):
    """Reload a saved artifact: ROOT RBDT for new BDTs, NumPy for other learners.

    Older BDT artifacts without a native export retain their XGBoost reader;
    they cannot supply the native on-read interface until separately exported.
    """
    from .Models import read_json
    path = Path(path); artifact = read_json(path)
    if artifact['schema'] == 'PhotonIDNumpyModelV1':
        return lambda x: numpy_predict(artifact, x)
    require(artifact['schema'] == 'PhotonIDXGBoostModelV1', 'unsupported model artifact')
    if 'native_evaluator' in artifact:
        specification = artifact['native_evaluator']
        require(specification['inputs'] == len(artifact['features']), 'native feature count differs')
        native_path = path.parent / specification['file']
        require(native_path.resolve().is_relative_to(path.parent.resolve()), 'native model path escapes its artifact folder')
        require(artifact['native_export_validation']['status'] == 'PASS', 'native export is not validated')
        return native_evaluator(native_path, specification)
    try:
        from xgboost import XGBClassifier
    except ImportError as exc:
        raise ValueError('this BDT release requires XGBoost in the selected environment') from exc
    model = XGBClassifier(n_jobs=1)
    model.load_model(path.parent / artifact['booster_file'])
    return lambda x: model.predict_proba(x)[:, 1].astype(np.float64)
