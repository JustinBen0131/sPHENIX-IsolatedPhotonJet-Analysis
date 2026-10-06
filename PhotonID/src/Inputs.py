"""Inputs: source-bound H70 training rows from completed eventTree groups.

Reuse TreeProduction's publication manifests and photon feature map through its
installed ROOT-only eventTree library. Two bounded ROOT passes write a
disk-backed matrix; SQLite checks identities and assigns an event-grouped
split once for every learner. Isolation never enters the labels. Weights
reproduce the existing class balance and spline ET/eta flattening.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3

import awkward as ak
import numpy as np

from .Models import SYSTEMS, digest, file_hash, progress, read_json, require, save_json

_EVENTTREE = None


def eventtree():
    """TreeProduction's installed eventTree library, bound through PyROOT.

    MYINSTALL names the installation prefix (lib/ and include/TreeProduction/).
    Only the ROOT product library is loaded, never a detector library.
    """
    global _EVENTTREE
    if _EVENTTREE is None:
        prefix = Path(os.environ.get('MYINSTALL', ''))
        require(os.environ.get('MYINSTALL') and (prefix / 'include/TreeProduction/EventTree.h').is_file(),
                'install TreeProduction (libEventTree) and set MYINSTALL to its installation prefix')
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.gInterpreter.AddIncludePath(str(prefix / 'include'))
        require(ROOT.gSystem.Load(str(prefix / 'lib' / 'libEventTree')) >= 0, 'cannot load libEventTree from ' + str(prefix))
        require(ROOT.gInterpreter.Declare('#include "TreeProduction/EventTree.h"'), 'cannot declare the eventTree interface')
        _EVENTTREE = ROOT.eventtree
    return _EVENTTREE


def run(function, *arguments):
    """One library call; its refusal becomes the ValueError this package reports."""
    try:
        return function(*arguments)
    except ValueError:
        raise
    except Exception as error:
        raise ValueError(str(error).strip().splitlines()[-1].split('runtime_error: ', 1)[-1]) from error


def metadata(path):
    """The saved file's direct meta/ records."""
    return json.loads(str(run(eventtree().metadataJson, str(path))))


def feature_branches():
    """Model input name -> eventTree branch, as TreeProduction stores it."""
    return {str(name): str(branch) for name, branch in eventtree().featureBranches()}


def published_inputs(groups, system):
    """Resolve each CURRENT once and verify its exact logical-source roster."""
    require(groups, 'no completed dataset groups configured; use --group or config datasets')
    inputs, seen = [], set()
    for group in groups:
        resolved = json.loads(str(run(eventtree().resolveGroup, str(group))))
        version, manifest = Path(resolved['directory']), resolved['manifest']
        require(manifest['stage'] in (2, 3), 'training/scoring requires a completed group')
        for item in manifest['files']:
            path = (version / item['relative_path']).resolve()
            require(path.is_relative_to(version) and path.is_file(), 'invalid group file path')
            require(file_hash(path) == item['sha256'], 'published input hash differs: ' + str(path))
            require(item['logical_id'] not in seen, 'duplicate logical source across selected groups')
            seen.add(item['logical_id'])
            meta = metadata(path)
            require(meta['schema']['system'] == SYSTEMS[system], 'mixed collision systems')
            require(meta['source']['logical_source_id'] == item['logical_id'], 'logical identity mismatch')
            require('consolidation' in meta, 'input is not a completed logical-source file')
            inputs.append(dict(path=str(path), sha256=item['sha256'], logical_id=item['logical_id'],
                               relative_path=item['relative_path'], group=str(Path(group).resolve()),
                               group_id=manifest['group_id'], version=manifest['version'],
                               source=meta['source']))
    return sorted(inputs, key=lambda item: item['logical_id'])


def event_chunks(path, *, all_branches=False, branches=None, step_size='16 MB'):
    """Chunk-local metadata keeps view/code dictionaries valid after consolidation."""
    import uproot
    meta = metadata(path)
    with uproot.open(path) as f:
        tree = f['eventTree']
        chunks = meta.get('consolidation', {}).get('chunks')
        require(chunks is not None, 'consolidated provenance missing')
        previous = 0
        for chunk in chunks:
            begin, end = chunk['entry_begin'], chunk['entry_end']
            require(begin == previous and begin <= end, 'chunk coverage overlaps or has gaps')
            previous = end
            for arrays in tree.iterate(None if all_branches else branches, entry_start=begin,
                                       entry_stop=end, step_size=step_size, library='ak'):
                yield arrays, chunk['original_meta']
        require(previous == tree.num_entries, 'chunk coverage is incomplete')


def feature_rows(arrays, meta, recipe):
    """Return one float32 feature vector per photon; prove named H70 parity."""
    BRANCH = feature_branches()
    counts = ak.to_numpy(ak.num(arrays['pho_et'])).astype(np.int64)
    n = int(counts.sum())
    if not n:
        # Photon-free chunks have no feature stride or H70 payload to verify.
        # Retain their events and emit aligned empty score vectors.
        empty = np.empty(0, np.float64)
        return np.empty((0, len(recipe['features'])), np.float32), counts, empty, empty, np.zeros(0, bool), np.zeros(0, bool)
    names = meta['constants']['showerViews.definition_name']
    require(names.count('H70') == 1 and len(names) == 7, 'H70 shower-view dictionary missing or ambiguous')
    view = names.index('H70')
    expected = 'fe6633650503f23cf3cfa662a2884d29efbcd7ca465b5b02cfe31ec220db3987'
    require(meta['constants']['showerViews.semantic_sha256'][view] == expected,
            'H70 feature semantics changed')
    floor = meta['constants']['showerViews.floor_gev'][view]['f64']
    require(np.uint64(floor).view(np.float64) == .07, 'H70 tower floor changed')
    columns = []
    for name in recipe['features']:
        branch = BRANCH[name]
        if branch in ('vtx_z', 'cent_value'):
            values = np.repeat(ak.to_numpy(arrays[branch]), counts)
        else:
            values = ak.to_numpy(ak.flatten(arrays[branch]))
            if branch.startswith('shv_'):
                require(len(values) == 7 * n, 'shower-view photon alignment differs')
                values = values.reshape(n, 7)[:, view]
        require(len(values) == n, 'feature photon alignment differs: ' + name)
        columns.append(values)
    x = np.column_stack(columns).astype(np.float32)
    # Stored H70 contains fourteen values in AuAu, eleven in pp.
    stride = meta.get('strides', meta.get('schema', {}).get('strides', {})).get('shv_features')
    require(stride == len(recipe['features']), 'H70 feature order/count differs')
    stored = ak.to_numpy(ak.flatten(arrays['shv_features'])).reshape(n, 7, stride)[:, view]
    require(np.array_equal(x, stored.astype(np.float32), equal_nan=True), 'named features differ from H70 vector')
    state = ak.to_numpy(ak.flatten(arrays['shv_finite_feature_state'])).reshape(n, 7)[:, view]
    valid = (state == 1) & np.isfinite(x).all(axis=1)
    et = ak.to_numpy(ak.flatten(arrays['pho_et']))
    eta = ak.to_numpy(ak.flatten(arrays['pho_eta']))
    cent = np.repeat(ak.to_numpy(arrays['cent_value']), counts)
    if recipe['system'] == 'AuAu':
        finite_cent = cent[np.isfinite(cent)]
        require(np.all(finite_cent == np.round(finite_cent)), 'AuAu centrality is not an integer bin')
    vertex = np.repeat(ak.to_numpy(arrays['vtx_z']), counts)
    domain = recipe['domain']
    inside = ((et >= domain['et_min']) & (et < domain['et_max_exclusive']) &
              (np.abs(eta) < domain['eta_abs_max']) & (np.abs(vertex) < domain['vertex_abs_max_cm']))
    if recipe['system'] == 'AuAu':
        inside &= (cent >= domain['centrality_min']) & (cent < domain['centrality_max_exclusive'])
    event_ok = ((ak.to_numpy(arrays['evt_record']) == 0) &
                (ak.to_numpy(arrays['evt_nodes_ready']) == 1) &
                (ak.to_numpy(arrays['vtx_valid']) == 1) &
                (ak.to_numpy(arrays['evt_terminal_status']) == 0))
    inside &= np.repeat(event_ok, counts)
    return x, counts, et, cent, inside, valid


def labels(arrays, role):
    """Dominant-primary barcode join, independent of isolated-analysis matches."""
    result = []
    for dominant, state, pid, barcodes, classes in zip(
            ak.to_list(arrays['pho_dom_barcode']), ak.to_list(arrays['pho_dom_state']),
            ak.to_list(arrays['pho_dom_pid']),
            ak.to_list(arrays['truthpho_generator_barcode']), ak.to_list(arrays['truthpho_prompt_class'])):
        require(len(dominant) == len(state) == len(pid) and len(barcodes) == len(classes), 'truth witness alignment differs')
        inventory = {}
        for barcode, klass in zip(barcodes, classes):
            inventory.setdefault(barcode, set()).add(klass)
        for barcode, validity, particle in zip(dominant, state, pid):
            matches = inventory.get(barcode, set())
            if validity != 2 or len(matches) > 1 or (particle == 22 and not matches):
                result.append(-1); continue
            prompt = bool(matches and next(iter(matches)) in (1, 2))
            if matches and next(iter(matches)) < 1:
                result.append(-1); continue
            result.append(1 if role == 'signal' and prompt else
                          0 if role == 'background' and not prompt else -1)
    return np.asarray(result, np.int8)


def source_role(source, recipe):
    sample = str(source.get('sample', '')).lower()
    matches = [role for role, samples in recipe['sources'].items() if sample in samples]
    require(len(matches) == 1, 'SIM sample not bound to the training recipe: ' + sample)
    require(source.get('si_di_role') in ('EMBEDDED', 'SI'), 'nominal training excludes DATA/DI')
    if recipe['system'] == 'AuAu':
        require(source['si_di_role'] == 'EMBEDDED', 'AuAu training requires embedded SIM')
    return matches[0]


def training_rows(item, recipe):
    """Read only training witnesses; exclude classes or events outside the recipe."""
    BRANCH = feature_branches()
    role = source_role(item['source'], recipe)
    branches = sorted(set([BRANCH[n] for n in recipe['features']] + [
        'pho_et','pho_eta','cent_value','vtx_z','vtx_valid','evt_record','evt_nodes_ready',
        'evt_terminal_status','evt_src_entry','evt_embedded_mb_decision',
        'shv_features','shv_finite_feature_state','pho_dom_barcode','pho_dom_state','pho_dom_pid',
        'truthpho_generator_barcode','truthpho_prompt_class']))
    if recipe['system'] == 'AuAu':
        branches = sorted(set(branches + [recipe['training']['minimum_bias_branch']]))
    for arrays, meta in event_chunks(item['path'], branches=branches):
        x, counts, et, cent, inside, valid = feature_rows(arrays, meta, recipe)
        y = labels(arrays, role)
        if recipe['system'] == 'AuAu':
            witness = recipe['training']['minimum_bias_branch']
            inside &= np.repeat(ak.to_numpy(arrays[witness]) == recipe['training']['minimum_bias_pass'], counts)
        keep = inside & valid & (y >= 0)
        entries = np.repeat(ak.to_numpy(arrays['evt_src_entry']), counts)
        photon_indices = np.concatenate([np.arange(c) for c in counts]) if counts.sum() else np.array([], np.int64)
        yield x[keep], y[keep], et[keep], cent[keep], entries[keep], photon_indices[keep]


@dataclass
class Population:
    x: np.ndarray
    y: np.ndarray
    et: np.ndarray
    centrality: np.ndarray
    weights: np.ndarray
    train: np.ndarray
    receipt: dict


def inverse_pdf(values, *, fixed_range=None, cap=None):
    """Port of the existing PPG12 exact twenty-bin spline weighting helper."""
    from scipy.interpolate import UnivariateSpline
    values = np.asarray(values, np.float64)
    require(np.isfinite(values).all() and len(values) >= 20, 'insufficient finite rows for exact weights')
    lo, hi = fixed_range or (float(values.min()), float(values.max()))
    require(hi > lo, 'zero-width weighting population')
    histogram, edges = np.histogram(values, bins=np.linspace(lo, hi, 21), density=True)
    spline = UnivariateSpline((edges[:-1] + edges[1:]) / 2, histogram * 20., s=0.)
    factors = 1. / np.clip(spline(values), 1e-3, None)
    if cap is not None:
        factors = np.minimum(factors, cap)
    return factors / factors.mean()


def canonical_weights(x, y, features):
    weights = np.ones(len(y), np.float64)
    for klass in (0, 1):
        mask = y == klass; n = int(mask.sum())
        require(n >= 20, 'both training classes need at least twenty eligible rows')
        weights[mask] = len(y) / (2. * n)
        weights[mask] *= inverse_pdf(x[mask, features.index('cluster_Eta')], fixed_range=(-.7, .7))
        weights[mask] *= inverse_pdf(x[mask, features.index('cluster_Et')], cap=800.)
    require(np.isfinite(weights).all() and (weights > 0).all(), 'invalid canonical training weights')
    return weights


def population_identity(inputs, recipe):
    """Learner and routing variations share preparation, weights and event split."""
    return digest(dict(inputs=inputs, features=recipe['features'], training=recipe['training'],
                       sources=recipe['sources'], domain=recipe['domain'], weights=recipe['weights'],
                       preparation_code=file_hash(__file__)))


def prepare(inputs, recipe, work):
    """Create/reuse one disk-backed population, with duplicate-free source keys."""
    work = Path(work); work.mkdir(parents=True, exist_ok=True)
    pin = population_identity(inputs, recipe)
    done = work / 'population.json'
    if done.exists():
        receipt = read_json(done)
        require(receipt['input_recipe_sha256'] == pin, 'prepared population belongs to different inputs')
        for item in inputs:
            require(file_hash(item['path']) == item['sha256'], 'training input changed')
        for filename, h in receipt['arrays'].items():
            require(file_hash(work / filename) == h, 'prepared matrix changed')
        return load_population(work, receipt)
    require(not (work / 'events.sqlite').exists(), 'incomplete population: retain evidence and use a fresh work directory')
    db = sqlite3.connect(work / 'events.sqlite')
    db.execute('CREATE TABLE events(source TEXT,entry INTEGER,hash TEXT,holdout INTEGER DEFAULT 0,PRIMARY KEY(source,entry))')
    db.execute('CREATE TABLE photons(source TEXT,entry INTEGER,photon INTEGER,PRIMARY KEY(source,entry,photon))')
    total = 0; counts = [0, 0]
    try:
        for index, item in enumerate(inputs):
            require(file_hash(item['path']) == item['sha256'], 'training input changed')
            for x, y, et, cent, entries, photons in training_rows(item, recipe):
                for entry, photon in zip(entries, photons):
                    key = f"{item['logical_id']}:{int(entry)}"
                    hash_text = '|'.join(('event50', str(recipe['training']['seed']),
                                          recipe['training']['split_identity'], key))
                    h = hashlib.sha256(hash_text.encode()).hexdigest()[:8]
                    db.execute('INSERT OR IGNORE INTO events(source,entry,hash) VALUES(?,?,?)',
                               (item['logical_id'], int(entry), h))
                    db.execute('INSERT INTO photons VALUES(?,?,?)', (item['logical_id'], int(entry), int(photon)))
                total += len(y)
                for klass in (0, 1): counts[klass] += int((y == klass).sum())
            db.commit()
            progress(work / 'progress.json', 'prepare-count', index + 1, len(inputs), eligible_photons=total)
        require(total > 0 and min(counts) >= 20, 'nominal training population is empty or missing a class')
        nevents = db.execute('SELECT COUNT(*) FROM events').fetchone()[0]
        nholdout = max(1, min(nevents - 1, round(nevents * .1)))
        db.execute('UPDATE events SET holdout=1 WHERE rowid IN (SELECT rowid FROM events ORDER BY hash,source,entry LIMIT ?)', (nholdout,))
        db.commit()
        shapes = {'x': (total, len(recipe['features'])), 'y': (total,), 'et': (total,),
                  'centrality': (total,), 'weights': (total,), 'train': (total,)}
        types = {'x': np.float32, 'y': np.int8, 'et': np.float64, 'centrality': np.float64,
                 'weights': np.float64, 'train': np.bool_}
        data = {name: np.lib.format.open_memmap(work / (name + '.npy'), mode='w+',
                                              dtype=types[name], shape=shape) for name, shape in shapes.items()}
        offset = 0
        for item in inputs:
            for x, y, et, cent, entries, photons in training_rows(item, recipe):
                n = len(y); sl = slice(offset, offset + n)
                data['x'][sl] = x; data['y'][sl] = y; data['et'][sl] = et; data['centrality'][sl] = cent
                data['train'][sl] = [not db.execute('SELECT holdout FROM events WHERE source=? AND entry=?',
                                                  (item['logical_id'], int(e))).fetchone()[0] for e in entries]
                offset += n
            require(file_hash(item['path']) == item['sha256'], 'input changed during population preparation')
            progress(work / 'progress.json', 'prepare-fill', offset, total)
        require(offset == total, 'training population changed between passes')
        data['weights'][:] = canonical_weights(data['x'], data['y'], recipe['features'])
        for array in data.values(): array.flush()
        receipt = dict(schema='PhotonIDPopulationV1', input_recipe_sha256=pin, rows=total,
                       class_counts=counts, events=nevents, holdout_events=nholdout,
                       arrays={name + '.npy': file_hash(work / (name + '.npy')) for name in data},
                       split_sha256=file_hash(work / 'train.npy'), isolation_used=False,
                       event_weights_used=False, sources=inputs)
        save_json(done, receipt)
    finally:
        db.close()
    return load_population(work, receipt)


def load_population(work, receipt):
    return Population(**{n: np.load(Path(work) / (n + '.npy'), mmap_mode='r') for n in
                         ('x', 'y', 'et', 'centrality', 'weights', 'train')}, receipt=receipt)
