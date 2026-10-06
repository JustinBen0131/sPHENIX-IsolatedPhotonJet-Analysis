"""Scoring: photon-aligned optional scores and conserved dataset publication.

TreeProduction's installed eventTree library writes each candidate and
publishes complete groups. Each model adds exactly one vector<double> branch;
outside-domain values are NaN. Saved outputs are reopened, all previous
payloads and metadata are compared, and scores are recomputed before a group
pointer may advance. Fraction and other optional scores remain intact when a
further model is added.
"""
from __future__ import annotations

import copy
import fcntl
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from .Inputs import event_chunks, eventtree, feature_rows, metadata, run
from .Models import SYSTEMS, canonical, file_hash, load_release, progress, read_json, require, save_json
from .Training import evaluator


def predictors(folder, manifest):
    return [(cell, evaluator(Path(folder) / cell['artifact'])) for cell in manifest['cells']]


def score_arrays(arrays, meta, recipe, models):
    x, counts, et, cent, inside, valid = feature_rows(arrays, meta, recipe)
    require(not np.any(inside & ~valid), 'invalid required features for an in-domain photon')
    scores = np.full(len(x), np.nan, np.float64)
    assigned = np.zeros(len(x), np.int8)
    for cell, predict in models:
        selected = inside.copy()
        for axis, (low, high) in cell['bounds'].items():
            values = et if axis == 'et' else cent
            selected &= (values >= low) & (values < high)
        indices = np.flatnonzero(selected)
        for begin in range(0, len(indices), 65536):
            part = indices[begin:begin + 65536]
            values = np.asarray(predict(x[part]), np.float64)
            require(values.shape == (len(part),) and np.isfinite(values).all() and
                    ((values >= 0) & (values <= 1)).all(), 'model scoring failed inside its domain')
            scores[part] = values
        assigned[selected] += 1
    require(np.all(assigned[inside] == 1) and np.all(assigned[~inside] == 0), 'model cells overlap or leave gaps')
    return ak.unflatten(scores, counts), int(inside.sum())


def same_scores(expected, saved):
    """Photon alignment and bit-identical values (NaN outside the domain included)."""
    if not np.array_equal(ak.to_numpy(ak.num(expected)), ak.to_numpy(ak.num(saved))):
        return False
    a = np.ascontiguousarray(ak.to_numpy(ak.flatten(expected)), dtype=np.float64)
    b = np.ascontiguousarray(ak.to_numpy(ak.flatten(saved)), dtype=np.float64)
    return a.view(np.uint64).tolist() == b.view(np.uint64).tolist()


def schema_path(code_root):
    path = Path(code_root).resolve() / 'TreeProduction/config/eventTree_v1.yaml'
    require(path.is_file(), 'the TreeProduction eventTree schema is required: ' + str(path))
    return path


def score_metadata(old, item, manifest):
    """Declare the one permitted metadata addition before writing."""
    meta = copy.deepcopy(old)
    meta.setdefault('photon_models', {})[manifest['name']] = dict(
        schema='PhotonIDScoreReceiptV1', status='PASS', release_sha256=manifest['release_sha256'],
        request_sha256=manifest['request_sha256'], score_branch=manifest['score_branch'],
        model_name=manifest['name'], domain=manifest['recipe']['domain'], outside_domain='NaN',
        input_sha256=item['sha256'], input_version=item['version'], logical_source_id=item['logical_id'])
    return meta


def prove(original, output, folder, manifest):
    """The library's preservation proof, then an independent rescore of every photon."""
    branch = manifest['score_branch']; models = predictors(folder, manifest)
    offset = 0; photons = 0
    old_meta, new_meta = metadata(original), metadata(output)
    receipt = new_meta.get('photon_models', {}).get(manifest['name'], {})
    require(receipt.get('release_sha256') == manifest['release_sha256'], 'saved score metadata binds another release')
    expected = copy.deepcopy(old_meta)
    expected.setdefault('photon_models', {})[manifest['name']] = receipt
    expected['schema']['counters'] = new_meta['schema']['counters']
    preservation = json.loads(str(run(eventtree().provePreserved, json.dumps(dict(
        original=str(original), candidate=str(output), added=[branch], metadata=expected)))))
    with uproot.open(output) as after:
        tree = after['eventTree']
        for arrays, chunk_meta in event_chunks(original, all_branches=True):
            saved = tree.arrays([branch], entry_start=offset, entry_stop=offset + len(arrays), library='ak')
            expected_score, n = score_arrays(arrays, chunk_meta, manifest['recipe'], models)
            require(same_scores(expected_score, saved[branch]), 'saved score differs from model recompute')
            photons += n; offset += len(arrays)
        require(offset == tree.num_entries, 'score proof coverage incomplete')
    return dict(ok=True, **preservation, scored_photons=photons, score_branch=branch)


def augment_file(item, folder, manifest, candidate, *, code_root):
    """Write one candidate with the added score branch through the library.

    An existing candidate is checked, never overwritten or regenerated."""
    candidate = Path(candidate); receipt_path = candidate.with_suffix('.receipt.json')
    meta = metadata(item['path'])
    require(meta['schema']['system'] == SYSTEMS[manifest['system']], 'wrong collision system for model')
    if not candidate.exists():
        models = predictors(folder, manifest)
        addition = eventtree().BranchAddition()
        run(addition.open, json.dumps(dict(original=item['path'], candidate=str(candidate), input_sha256=item['sha256'],
                                           schema_path=str(schema_path(code_root)),
                                           added=[dict(name=manifest['score_branch'], type='double',
                                                       title='validated optional photon score')],
                                           metadata=score_metadata(meta, item, manifest))))
        total = json.loads(str(run(addition.status)))['entries']; done = 0
        for arrays, chunk_meta in event_chunks(item['path'], all_branches=True):
            if len(arrays):
                scores, _ = score_arrays(arrays, chunk_meta, manifest['recipe'], models)
            else:
                scores = ak.unflatten(np.empty(0, np.float64), np.empty(0, np.int64))
            counts = np.ascontiguousarray(ak.to_numpy(ak.num(scores)), dtype=np.longlong)
            values = np.ascontiguousarray(ak.to_numpy(ak.flatten(scores)), dtype=np.float64)
            run(addition.append, manifest['score_branch'], len(counts), values, len(values), counts)
            run(addition.commit, len(counts))
            done += len(counts)
            progress(candidate.parent / 'progress.json', 'score', done, total,
                     logical_source_id=item['logical_id'], model=manifest['name'])
        run(addition.finish)
    proof = prove(item['path'], candidate, folder, manifest)
    require(file_hash(item['path']) == item['sha256'], 'input changed during scoring')
    receipt = dict(schema='PhotonIDAugmentationReceiptV1', status='PASS', logical_source_id=item['logical_id'],
                   input_sha256=item['sha256'], input_version=item['version'], output=str(candidate),
                   output_sha256=file_hash(candidate), release_sha256=manifest['release_sha256'], proof=proof)
    save_json(receipt_path, receipt)
    return dict(logical_id=item['logical_id'], relative_path=item['relative_path'], path=str(candidate),
                sha256=receipt['output_sha256'], receipt_path=str(receipt_path), receipt_sha256=file_hash(receipt_path))


def score_groups(inputs, folder, *, code_root):
    manifest = load_release(folder)
    groups = {}
    for item in inputs:
        groups.setdefault(item['group'], []).append(item)
    results = []
    for group, items in groups.items():
        present = []
        for item in items:
            old = metadata(item['path']).get('photon_models', {}).get(manifest['name'])
            if old:
                require(old['release_sha256'] == manifest['release_sha256'], 'existing score binds a different release')
            present.append(old is not None)
        require(all(present) or not any(present), 'published group has partial optional-model coverage')
        if all(present):
            # Re-evaluate saved scores before declaring idempotent reuse.
            for item in items:
                models = predictors(folder, manifest)
                for arrays, meta in event_chunks(item['path'], all_branches=True):
                    scores, _ = score_arrays(arrays, meta, manifest['recipe'], models)
                    require(same_scores(scores, arrays[manifest['score_branch']]), 'existing saved score differs')
                require(file_hash(item['path']) == item['sha256'], 'reused group input changed')
            results.append(dict(group=group, status='REUSED', version=items[0]['version'])); continue
        base_version = items[0]['version']
        require(all(i['version'] == base_version for i in items), 'group input version differs')
        work = Path(group) / '.photonid_candidates' / manifest['release_sha256'] / base_version
        work.mkdir(parents=True, exist_ok=True)
        with (work / '.scoring.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                # Every candidate is proven when written and bound to its receipt by
                # digest; publication rehashes each file before the pointer moves.
                records = [augment_file(i, folder, manifest, work / i['relative_path'], code_root=code_root) for i in items]
                version = json.loads(str(run(eventtree().publishGroup, json.dumps(dict(
                    group=str(group), files=records, group_id=items[0]['group_id'],
                    expected_sources=[i['logical_id'] for i in items],
                    recipe=dict(schema='PhotonIDPublicationRecipeV1', model=manifest['name'],
                                release_sha256=manifest['release_sha256']),
                    base_version=base_version, stage=3)))))
            except Exception as exc:
                progress(work / 'progress.json', 'blocked', 0, len(items), model=manifest['name'],
                         group=str(group), unit='logical-source files', latest_error=str(exc))
                raise
        results.append(dict(group=group, status='PASS', version=version['version'], files=len(items)))
    return results
