"""Scoring: photon-aligned optional scores and conserved dataset publication.

Use the existing eventTree streaming writer and complete-group publisher. Each
model adds exactly one vector<double> branch; outside-domain values are NaN.
Saved outputs are reopened, all previous payloads and metadata are compared,
and scores are recomputed before a group pointer may advance. Fraction and
other optional scores remain intact when a further model is added.
"""
from __future__ import annotations

import copy
import fcntl
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from .Inputs import event_chunks, feature_rows, use_format
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


def extended_schema(meta, new_branch):
    from eventtree.schema import Branch, Schema
    schema = Schema()
    for spec in meta.get('photon_models', {}).values():
        name = spec['score_branch']
        require(name.startswith('pho_') and name not in schema.branches, 'optional branch collides with base schema')
        schema.branches[name] = Branch(name, 'pho', 'vector', np.float64, 0, 'optional:model',
                                       'validated optional photon score', '')
    require(new_branch not in schema.branches, 'score branch already exists')
    schema.branches[new_branch] = Branch(new_branch, 'pho', 'vector', np.float64, 0,
                                         'optional:model', 'validated optional photon score', '')
    return schema


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
    """Use the tree updater's preservation proof, then independently rescore."""
    from eventtree.integrity import metadata, equal_arrays
    from eventtree.update import prove_preserved
    branch = manifest['score_branch']; models = predictors(folder, manifest)
    offset = 0; photons = 0
    with uproot.open(original) as before, uproot.open(output) as after:
        old_meta, new_meta = metadata(before), metadata(after)
        receipt = new_meta.get('photon_models', {}).get(manifest['name'], {})
        require(receipt.get('release_sha256') == manifest['release_sha256'],
                'saved score metadata binds another release')
        expected = copy.deepcopy(old_meta)
        expected.setdefault('photon_models', {})[manifest['name']] = receipt
        expected['schema']['counters'] = new_meta['schema']['counters']
        preservation = prove_preserved(original, output, added={branch}, expected_metadata=expected)
        tree = after['eventTree']
        for arrays, chunk_meta in event_chunks(original, all_branches=True):
            saved = tree.arrays([branch], entry_start=offset, entry_stop=offset + len(arrays), library='ak')
            expected_score, n = score_arrays(arrays, chunk_meta, manifest['recipe'], models)
            require(equal_arrays(expected_score, saved[branch]), 'saved score differs from model recompute')
            require(bool(ak.all(ak.num(saved[branch]) == ak.num(arrays['pho_et']))), 'score/photon alignment differs')
            photons += n; offset += len(arrays)
        require(offset == tree.num_entries, 'score proof coverage incomplete')
    return dict(ok=True, **preservation, scored_photons=photons, score_branch=branch)


def augment_file(item, folder, manifest, candidate):
    from eventtree.integrity import metadata
    from eventtree.update import add_branches
    candidate = Path(candidate); receipt_path = candidate.with_suffix('.receipt.json')
    models = predictors(folder, manifest)
    with uproot.open(item['path']) as source:
        meta = metadata(source)
        require(meta['schema']['system'] == SYSTEMS[manifest['system']], 'wrong collision system for model')
        schema = extended_schema(meta, manifest['score_branch'])
    def transform(arrays, chunk_meta):
        if len(arrays):
            scores, _ = score_arrays(arrays, chunk_meta, manifest['recipe'], models)
        else:
            scores = ak.unflatten(np.empty(0, np.float64), np.empty(0, np.int64))
        return ak.with_field(arrays, scores, manifest['score_branch'])
    add_branches(item['path'], candidate, input_sha256=item['sha256'], schema=schema,
                 added={manifest['score_branch']},
                 chunks=lambda: event_chunks(item['path'], all_branches=True), transform=transform,
                 new_metadata=score_metadata(meta, item, manifest),
                 progress=lambda done, total: progress(candidate.parent / 'progress.json', 'score', done, total,
                     logical_source_id=item['logical_id'], model=manifest['name']))
    proof = prove(item['path'], candidate, folder, manifest)
    require(file_hash(item['path']) == item['sha256'], 'input changed during scoring')
    receipt = dict(schema='PhotonIDAugmentationReceiptV1', status='PASS', logical_source_id=item['logical_id'],
                   input_sha256=item['sha256'], input_version=item['version'], output=str(candidate),
                   output_sha256=file_hash(candidate), release_sha256=manifest['release_sha256'], proof=proof)
    save_json(receipt_path, receipt)
    return dict(logical_id=item['logical_id'], relative_path=item['relative_path'], path=str(candidate),
                sha256=receipt['output_sha256'], receipt_path=str(receipt_path), receipt_sha256=file_hash(receipt_path))


def score_groups(inputs, folder, *, code_root):
    use_format(code_root)
    from eventtree.publication import publish
    from eventtree.integrity import metadata
    manifest = load_release(folder)
    groups = {}
    for item in inputs:
        groups.setdefault(item['group'], []).append(item)
    results = []
    for group, items in groups.items():
        present = []
        for item in items:
            with uproot.open(item['path']) as f:
                meta = metadata(f)
                old = meta.get('photon_models', {}).get(manifest['name'])
                if old:
                    require(old['release_sha256'] == manifest['release_sha256'], 'existing score binds a different release')
                present.append(old is not None)
        require(all(present) or not any(present), 'published group has partial optional-model coverage')
        if all(present):
            # Re-evaluate saved scores before declaring idempotent reuse.
            for item in items:
                models = predictors(folder, manifest)
                with uproot.open(item['path']) as f:
                    for arrays, meta in event_chunks(item['path'], all_branches=True):
                        scores, _ = score_arrays(arrays, meta, manifest['recipe'], models)
                        from eventtree.integrity import equal_arrays
                        require(equal_arrays(scores, arrays[manifest['score_branch']]), 'existing saved score differs')
                require(file_hash(item['path']) == item['sha256'], 'reused group input changed')
            results.append(dict(group=group, status='REUSED', version=items[0]['version'])); continue
        base_version = items[0]['version']
        require(all(i['version'] == base_version for i in items), 'group input version differs')
        work = Path(group) / '.photonid_candidates' / manifest['release_sha256'] / base_version
        work.mkdir(parents=True, exist_ok=True)
        with (work / '.scoring.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                records = [augment_file(i, folder, manifest, work / i['relative_path']) for i in items]
                originals = {str(Path(r['path']).resolve()): i for r, i in zip(records, items)}
                version = publish(group, records, group_id=items[0]['group_id'],
                                  expected_sources=[i['logical_id'] for i in items],
                                  recipe=dict(schema='PhotonIDPublicationRecipeV1', model=manifest['name'],
                                              release_sha256=manifest['release_sha256']),
                                  base_version=base_version, stage=3,
                                  validator=lambda path: prove(originals[str(Path(path).resolve())]['path'], path, folder, manifest))
            except Exception as exc:
                progress(work / 'progress.json', 'blocked', 0, len(items), model=manifest['name'],
                         group=str(group), unit='logical-source files', latest_error=str(exc))
                raise
        results.append(dict(group=group, status='PASS', version=version['version'], files=len(items)))
    return results
