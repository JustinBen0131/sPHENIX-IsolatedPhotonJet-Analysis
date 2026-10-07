"""Compact, source-bound inputs for independently scheduled routing cells.

Preparation copies complete cells from verified population arrays in their
original row order. Global weights and frozen partitions are copied verbatim;
this module computes no features, labels, splits or weights. Each fit verifies
only its compact input. Release assembly uses original zero-based row indices
to restore the full population order before diagnostics.
"""
from __future__ import annotations

import fcntl
import os
from pathlib import Path
import re
from types import SimpleNamespace

import numpy as np

from . import Inputs, Models, Training


def directory(root, request):
    """One immutable input family per resolved recipe/input/runtime request."""
    Models.require(isinstance(request, str) and re.fullmatch('[0-9a-f]{64}', request),
                   'invalid compact cell request identity')
    root = Path(root).resolve(); folder = root / '.cells' / request
    Models.require(folder.resolve().is_relative_to(root), 'compact cell directory escapes model root')
    return folder


def record_dtype(columns):
    """Retain each source array's dtype and trailing dimensions exactly."""
    return np.dtype([(c['name'], c['dtype'], tuple(c['shape'])) if c['shape'] else
                     (c['name'], c['dtype']) for c in columns])


def population_witness(receipt):
    """Small provenance record; the full source roster stays in preparation."""
    return {key: receipt[key] for key in ('input_recipe_sha256', 'rows',
            'split_sha256', 'event_split_sha256', 'split_method')}


def prepare(root, request, recipe, population, work, counts, *, progress_identity=None):
    """Copy one routing family in bounded chunks, checking actual census counts.

    counts comes from the verified preparation census. An interrupted family
    is retained, never overwritten; use a new release version for a new attempt.
    """
    cells = list(Training.cells(recipe)); folder = directory(root, request)
    expected = [(c['id'], c['bounds']) for c in cells]
    Models.require([(c['id'], c['bounds']) for c in counts] == expected,
                   'compact input census has different routing cells')
    Models.require(all(type(c['rows']) is int and c['rows'] > 0 for c in counts) and
                   sum(c['rows'] for c in counts) == len(population.y),
                   'compact input census does not cover the population')
    folder.mkdir(parents=True, exist_ok=True)
    tracker = Models.WorkProgress(folder / 'progress.json', len(population.y),
        'prepared photons', 'prepare-cells-' + recipe['name'], progress_identity)
    with (folder / '.preparation.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (folder / 'family.json').exists():
            for cell in cells:
                load(root, request, recipe, cell, expected_population=population_witness(population.receipt))
            return dict(status='REUSED', family=str(folder / 'family.json'), cells=len(cells))
        Models.require(not any(folder.glob('cell_*')),
                       'interrupted compact preparation retained; choose a new release version')
        arrays = {name: np.load(Path(work) / filename, mmap_mode='r', allow_pickle=False)
                  for name, filename in population.receipt['array_files'].items()}
        Models.require(set(('x', 'y', 'et', 'centrality', 'weights', 'train', 'partition',
                            'split_rank', 'source_index', 'entry', 'photon', 'event_key')) <= set(arrays),
                       'prepared population lacks cell identity witnesses')
        Models.require(all(len(a) == len(population.y) for a in arrays.values()),
                       'prepared cell arrays are not aligned')
        columns = [dict(name='original_index', dtype=np.dtype(np.int64).str, shape=[])] + [
            dict(name=name, dtype=a.dtype.str, shape=list(a.shape[1:])) for name, a in arrays.items()]
        dtype = record_dtype(columns)
        # Stream output chunks instead of keeping every destination mmap resident.
        # The ordinary NumPy header retains the same dtype/shape contract.
        targets = {}
        for c in counts:
            stream = (folder / (c['id'] + '.npy')).open('xb')
            targets[c['id']] = stream
            np.lib.format.write_array_header_1_0(stream, dict(
                descr=np.lib.format.dtype_to_descr(dtype), fortran_order=False, shape=(c['rows'],)))
        positions = {c['id']: 0 for c in cells}; tallies = {c['id']: np.zeros(6, np.int64) for c in cells}
        source_paths = [Path(work) / name for name in population.receipt['arrays']]
        before = [(p.stat().st_dev, p.stat().st_ino, p.stat().st_size,
                   p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in source_paths]
        tracker.write('copy-cells', 0)
        population_receipt_sha = Models.file_hash(Path(work) / 'population.json')
        try:
            for begin in range(0, len(population.y), 65536):
                end = min(begin + 65536, len(population.y)); section = slice(begin, end)
                view = SimpleNamespace(y=population.y[section], et=population.et[section],
                                       centrality=population.centrality[section])
                coverage = np.zeros(end - begin, np.uint8)
                for cell, count in zip(cells, counts):
                    selected = np.flatnonzero(Training.cell_mask(view, cell))
                    coverage[selected] += 1
                    offset = positions[cell['id']]; stop = offset + len(selected)
                    Models.require(stop <= count['rows'], 'compact cell exceeds its census')
                    target = np.empty(len(selected), dtype=dtype)
                    target['original_index'] = begin + selected
                    for name, array in arrays.items():
                        target[name] = array[section][selected]
                    tallies[cell['id']] += np.bincount(
                        2 * target['partition'].astype(np.int64) + target['y'], minlength=6)
                    targets[cell['id']].write(target.tobytes())
                    positions[cell['id']] = stop
                Models.require(np.all(coverage == 1), 'compact routing overlaps or misses population rows')
                tracker.write('copy-cells', end)
            Models.require(before == [(p.stat().st_dev, p.stat().st_ino, p.stat().st_size,
                           p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in source_paths],
                           'prepared population changed during compact copy')
            entries = []
            for cell, count in zip(cells, counts):
                Models.require(positions[cell['id']] == count['rows'] and
                               tallies[cell['id']].reshape(3, 2).tolist() == count['partition_class_counts'],
                               'compact cell count or partition differs from census')
                targets[cell['id']].flush(); os.fsync(targets[cell['id']].fileno())
                targets[cell['id']].close()
                path = folder / (cell['id'] + '.npy')
                body = dict(schema='PhotonIDCellInputV1', status='PASS', request_sha256=request,
                    recipe_sha256=Models.digest(recipe), cell=cell, rows=count['rows'],
                    partition_class_counts=count['partition_class_counts'], columns=columns,
                    population=population_witness(population.receipt),
                    population_receipt_sha256=population_receipt_sha,
                    data_file=path.name, data_sha256=Models.file_hash(path))
                receipt = folder / (cell['id'] + '.json'); Models.save_json(receipt, body)
                entries.append(dict(cell, rows=count['rows'], receipt_sha256=Models.file_hash(receipt)))
            body = dict(schema='PhotonIDCellFamilyV1', status='PASS', request_sha256=request,
                recipe_sha256=Models.digest(recipe), population=population_witness(population.receipt),
                population_receipt_sha256=population_receipt_sha, cells=entries)
            Models.save_json(folder / 'family.json', dict(body, family_sha256=Models.digest(body)))
            tracker.write('complete', len(population.y), status='PASS')
        except BaseException as error:
            tracker.fail(error)
            raise
        finally:
            for stream in targets.values():
                stream.close()
    return dict(status='PREPARED', family=str(folder / 'family.json'), cells=len(cells))


def load(root, request, recipe, cell, *, expected_population=None):
    """Verify one finite cell and return its population, original indices and pin.

    Original indices are strictly increasing and refer to the globally prepared
    row array. They are neither a new split nor a reordering for fitting.
    """
    folder = directory(root, request); family = Models.read_json(folder / 'family.json')
    Models.require(family.get('schema') == 'PhotonIDCellFamilyV1' and family.get('status') == 'PASS' and
                   Models.digest({k: v for k, v in family.items() if k != 'family_sha256'}) == family['family_sha256'],
                   'compact cell family is incomplete or changed')
    Models.require(family['request_sha256'] == request and family['recipe_sha256'] == Models.digest(recipe),
                   'compact cell family binds another recipe/input/runtime request')
    expected = list(Training.cells(recipe))
    Models.require([(c['id'], c['bounds']) for c in family['cells']] ==
                   [(c['id'], c['bounds']) for c in expected] and
                   sum(c['rows'] for c in family['cells']) == family['population']['rows'],
                   'compact cell family has missing, duplicate or inconsistent cells')
    Models.require(cell in expected, 'unknown compact routing cell')
    if expected_population is not None:
        Models.require(family['population'] == expected_population,
                       'compact cells bind another prepared population or split')
    entry = next(c for c in family['cells'] if c['id'] == cell['id'])
    receipt = folder / (cell['id'] + '.json')
    Models.require(Models.file_hash(receipt) == entry['receipt_sha256'], 'compact cell receipt changed')
    body = Models.read_json(receipt)
    Models.require(body.get('schema') == 'PhotonIDCellInputV1' and body.get('status') == 'PASS' and
                   body['request_sha256'] == request and body['recipe_sha256'] == Models.digest(recipe) and
                   body['cell'] == cell and body['rows'] == entry['rows'] and
                   body['population'] == family['population'] and
                   body['population_receipt_sha256'] == family['population_receipt_sha256'],
                   'compact cell receipt identity differs')
    path = folder / (cell['id'] + '.npy')
    Models.require(body['data_file'] == path.name and Models.file_hash(path) == body['data_sha256'],
                   'compact cell data changed')
    data = np.load(path, mmap_mode='r', allow_pickle=False)
    Models.require(data.dtype == record_dtype(body['columns']) and data.shape == (body['rows'],),
                   'compact cell data layout differs')
    Models.require(set(('original_index', 'x', 'y', 'et', 'centrality', 'weights', 'train',
                        'partition', 'split_rank', 'source_index', 'entry', 'photon', 'event_key')) <= set(data.dtype.names),
                   'compact cell data lacks row identity witnesses')
    Models.require(data['x'].shape == (body['rows'], len(recipe['features'])) and data['x'].dtype == np.float32,
                   'compact cell feature shape or precision differs')
    indices = data['original_index']; previous = -1; tally = np.zeros(6, np.int64)
    for begin in range(0, len(data), 65536):
        part = data[begin:begin + 65536]; index = part['original_index']
        Models.require(len(index) and previous < int(index[0]) and np.all(index[1:] > index[:-1]) and
                       int(index[-1]) < family['population']['rows'], 'compact cell original row order differs')
        previous = int(index[-1])
        Models.require(np.all((part['y'] == 0) | (part['y'] == 1)) and
                       np.all((part['partition'] >= 0) & (part['partition'] <= 2)) and
                       np.array_equal(part['train'], part['partition'] == 0) and
                       np.all(Training.cell_mask(SimpleNamespace(y=part['y'], et=part['et'],
                                      centrality=part['centrality']), cell)), 'compact cell routing or partition differs')
        tally += np.bincount(2 * part['partition'].astype(np.int64) + part['y'], minlength=6)
    Models.require(len(data) > 0 and tally.reshape(3, 2).tolist() == body['partition_class_counts'],
                   'compact cell class counts differ')
    population = Inputs.Population(**{name: data[name] for name in
        ('x', 'y', 'et', 'centrality', 'weights', 'train', 'partition', 'split_rank')},
        receipt=body['population'])
    return population, indices, Models.digest(body), data


def verify_origin(data, work, receipt):
    """Check compact bytes against their original prepared rows before sealing.

    This bounded assembly check includes the source/event/photon witnesses as
    well as learner inputs, global weights and partitions. No recomputation or
    approximate numeric comparison is used.
    """
    arrays = {name: np.load(Path(work) / filename, mmap_mode='r', allow_pickle=False)
              for name, filename in receipt['array_files'].items()}
    for begin in range(0, len(data), 65536):
        part = data[begin:begin + 65536]; indices = part['original_index']
        for name, original in arrays.items():
            Models.require(np.ascontiguousarray(part[name]).tobytes() ==
                           np.ascontiguousarray(original[indices]).tobytes(),
                           'compact bytes differ from original population: ' + name)
