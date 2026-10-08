"""Scoring: photon-aligned optional scores and conserved dataset publication.

TreeProduction's installed eventTree library writes each candidate and
publishes complete groups. All selected models share one source augmentation;
each model adds one photon-aligned double array. Outside-domain values are
NaN. Saved outputs are reopened, all previous
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

from .Inputs import event_chunks, eventtree, feature_branches, feature_rows, metadata, run
from .Models import SYSTEMS, WorkProgress, canonical, digest, file_hash, load_release, progress, require, save_json
from .Training import evaluator


def predictors(folder, manifest):
    return [(cell, evaluator(Path(folder) / cell['artifact'])) for cell in manifest['cells']]


def score_arrays(arrays, meta, recipe, models, *, rows=None):
    """Return photon-aligned jagged scores and the number inside the model domain.

    Outside-domain photons keep NaN. Every inside-domain photon needs valid
    inputs, exactly one routing cell and a finite prediction, or the file fails.
    """
    x, counts, et, cent, inside, valid = rows if rows is not None else feature_rows(arrays, meta, recipe)
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
    return np.array_equal(a.view(np.uint64), b.view(np.uint64))


def selected_releases(folder, manifest=None):
    """Normalize a single release or a collection without changing its recipe."""
    folders = [folder] if isinstance(folder, (str, Path)) else list(folder)
    manifests = ([manifest] if isinstance(manifest, dict) else list(manifest)) if manifest is not None else [
        load_release(path) for path in folders]
    require(folders and len(folders) == len(manifests), 'select at least one model release')
    require(len({m['system'] for m in manifests}) == 1, 'score one collision system per augmentation')
    require(len({m['name'] for m in manifests}) == len(manifests), 'duplicate selected model name')
    require(len({m['score_branch'] for m in manifests}) == len(manifests), 'selected model score branches collide')
    return list(zip(map(Path, folders), manifests))


def score_chunk(arrays, meta, releases, models):
    """Read H70 rows once per distinct feature/domain contract, then route every model."""
    require(len(releases) == len(models), 'selected model evaluator coverage differs')
    rows, scores, counts = {}, {}, {}
    for (_, manifest), predictors_ in zip(releases, models):
        recipe = manifest['recipe']
        key = canonical({k: recipe[k] for k in ('system', 'features', 'domain')})
        if key not in rows:
            rows[key] = feature_rows(arrays, meta, recipe)
        values, count = score_arrays(arrays, meta, recipe, predictors_, rows=rows[key])
        scores[manifest['score_branch']] = values
        counts[manifest['name']] = count
    return scores, counts


def input_branches(releases):
    """Only the H70 witnesses and application domain enter model evaluation.

    The native writer/proof preserves the other branches without loading jets,
    truth payloads or earlier ID decisions into Python scoring batches.
    """
    names = feature_branches()
    branches = {'pho_et', 'pho_eta', 'cent_value', 'vtx_z', 'vtx_valid',
                'evt_record', 'evt_nodes_ready', 'evt_terminal_status',
                'shv_features', 'shv_finite_feature_state'}
    branches.update(names[name] for _, model in releases for name in model['recipe']['features'])
    return sorted(branches)


def schema_path(code_root):
    path = Path(code_root).resolve() / 'TreeProduction/config/eventTree_v1.yaml'
    require(path.is_file(), 'the TreeProduction eventTree schema is required: ' + str(path))
    return path


def score_metadata(old, item, manifest):
    """Declare only the selected model receipts, preserving every earlier addition."""
    manifests = [manifest] if isinstance(manifest, dict) else manifest
    meta = copy.deepcopy(old)
    for model in manifests:
        require(model['name'] not in meta.get('photon_models', {}), 'model score already exists in augmentation base')
        meta.setdefault('photon_models', {})[model['name']] = dict(
            schema='PhotonIDScoreReceiptV1', status='PASS', release_sha256=model['release_sha256'],
            request_sha256=model['request_sha256'], score_branch=model['score_branch'],
            model_name=model['name'], domain=model['recipe']['domain'], outside_domain='NaN',
            input_sha256=item['sha256'], input_version=item['version'], logical_source_id=item['logical_id'])
    return meta


def prove(original, output, folder, manifest, *, model_predictors=None, input_item=None, tracker=None):
    """One native preservation proof and one saved-score read for all selected models."""
    releases = selected_releases(folder, manifest)
    branches = [m['score_branch'] for _, m in releases]
    models = model_predictors if model_predictors is not None else [predictors(path, m) for path, m in releases]
    offset = 0; photons = {m['name']: 0 for _, m in releases}
    old_meta, new_meta = metadata(original), metadata(output)
    expected = score_metadata(old_meta, input_item, [m for _, m in releases]) if input_item else copy.deepcopy(old_meta)
    for _, model in releases:
        receipt = new_meta.get('photon_models', {}).get(model['name'], {})
        require(receipt.get('release_sha256') == model['release_sha256'] and
                receipt.get('score_branch') == model['score_branch'] and
                receipt.get('domain') == model['recipe']['domain'], 'saved score metadata binds another release/domain')
        if input_item:
            require(receipt == expected['photon_models'][model['name']], 'saved score source/version receipt differs')
        expected.setdefault('photon_models', {})[model['name']] = receipt
    expected['schema']['counters'] = new_meta['schema']['counters']
    preservation = json.loads(str(run(eventtree().provePreserved, json.dumps(dict(
        original=str(original), candidate=str(output), added=branches, metadata=expected)))))
    with uproot.open(output) as after:
        tree = after['eventTree']
        for arrays, chunk_meta in event_chunks(original, branches=input_branches(releases)):
            saved = tree.arrays(branches, entry_start=offset, entry_stop=offset + len(arrays), library='ak')
            expected_scores, n = score_chunk(arrays, chunk_meta, releases, models)
            for branch in branches:
                require(same_scores(expected_scores[branch], saved[branch]), 'saved score differs from model recompute: ' + branch)
            for name in photons: photons[name] += n[name]
            offset += len(arrays)
            if tracker: tracker.write('verify-scores', tracker.total // 2 + offset)
        require(offset == tree.num_entries, 'score proof coverage incomplete')
    result = dict(ok=True, **preservation, scored_photons=photons, score_branches=branches)
    if len(releases) == 1:
        result.update(scored_photons=next(iter(photons.values())), score_branch=branches[0])
    return result


def augment_file(item, folder, manifest, candidate, *, code_root, model_predictors=None, progress_identity=None):
    """Write one candidate with all selected score branches through the library.

    An existing candidate is checked, never overwritten or regenerated."""
    candidate = Path(candidate); receipt_path = candidate.with_suffix('.receipt.json')
    releases = selected_releases(folder, manifest)
    manifests = [m for _, m in releases]
    meta = metadata(item['path'])
    require(meta['schema']['system'] == SYSTEMS[manifests[0]['system']], 'wrong collision system for model')
    require(meta['source']['logical_source_id'] == item['logical_id'], 'augmentation source identity differs')
    total = meta['schema']['entries']; reused = candidate.exists()
    tracker = WorkProgress(candidate.parent / 'progress.json', 2 * total, 'scoring event visits',
                           item['logical_id'], progress_identity, skip_total=total if reused else 0)
    models = model_predictors if model_predictors is not None else [predictors(path, m) for path, m in releases]
    try:
        tracker.write('reuse-candidate' if reused else 'copy-source', total if reused else 0,
                      skip_processed=total if reused else 0, models=[m['name'] for m in manifests])
        if not reused:
            addition = eventtree().BranchAddition()
            run(addition.open, json.dumps(dict(original=item['path'], candidate=str(candidate), input_sha256=item['sha256'],
                                               schema_path=str(schema_path(code_root)),
                                               added=[dict(name=m['score_branch'], type='double',
                                                           title='validated optional photon score') for m in manifests],
                                               metadata=score_metadata(meta, item, manifests))))
            require(json.loads(str(run(addition.status)))['entries'] == total, 'source event total changed')
            done = 0
            for arrays, chunk_meta in event_chunks(item['path'], branches=input_branches(releases)):
                scores, _ = score_chunk(arrays, chunk_meta, releases, models)
                for branch, values in scores.items():
                    counts = np.ascontiguousarray(ak.to_numpy(ak.num(values)), dtype=np.longlong)
                    flat = np.ascontiguousarray(ak.to_numpy(ak.flatten(values)), dtype=np.float64)
                    run(addition.append, branch, len(counts), flat, len(flat), counts)
                run(addition.commit, len(arrays))
                done += len(arrays); tracker.write('score', done)
            run(addition.finish)
        tracker.write('preserve-fields', total)
        proof = prove(item['path'], candidate, folder, manifest, model_predictors=models,
                      input_item=item, tracker=tracker)
        tracker.write('seal-artifact')
        require(file_hash(item['path']) == item['sha256'], 'input changed during scoring')
        receipt = dict(schema='PhotonIDAugmentationReceiptV1', status='PASS', logical_source_id=item['logical_id'],
                       input_sha256=item['sha256'], input_version=item['version'], output=str(candidate),
                       output_sha256=file_hash(candidate),
                       models={m['name']: m['release_sha256'] for m in manifests}, proof=proof)
        if len(manifests) == 1: receipt['release_sha256'] = manifests[0]['release_sha256']
        if receipt_path.exists():
            require(json.loads(receipt_path.read_text()) == receipt, 'existing augmentation receipt differs')
        else:
            save_json(receipt_path, receipt)
        tracker.write('complete', tracker.total, status='PASS')
    except BaseException as error:
        tracker.fail(error)
        raise
    return dict(logical_id=item['logical_id'], relative_path=item['relative_path'], path=str(candidate),
                sha256=receipt['output_sha256'], receipt_path=str(receipt_path), receipt_sha256=file_hash(receipt_path))


def score_groups(inputs, folder, *, code_root, candidate_root=None, progress_identity=None):
    """Augment every file of a resolved group before advancing its CURRENT pointer.

    Candidates preserve the logical-source roster and original payloads.
    The native publication compare-and-swap catches a changed base version;
    failed or incomplete candidates remain available for inspection. Supplying
    candidate_root stages proofs only; the shared publisher owns any later
    whole-group commit and never uses stale candidates as a new base.
    """
    releases = selected_releases(folder)
    manifests = [m for _, m in releases]
    groups = {}
    for item in inputs:
        groups.setdefault(item['group'], []).append(item)
    results = []
    for group, items in groups.items():
        saved_metadata = [metadata(item['path']).get('photon_models', {}) for item in items]
        missing = []
        for path, manifest in releases:
            present = []
            for old_models in saved_metadata:
                old = old_models.get(manifest['name'])
                if old:
                    require(old['release_sha256'] == manifest['release_sha256'] and
                            old['score_branch'] == manifest['score_branch'] and
                            old['domain'] == manifest['recipe']['domain'], 'existing score binds a different release/domain')
                present.append(old is not None)
            require(all(present) or not any(present), 'published group has partial optional-model coverage')
            if not all(present): missing.append((path, manifest))
        if not missing:
            # Re-evaluate saved scores before declaring idempotent reuse.
            models = [predictors(path, m) for path, m in releases]
            for item in items:
                for arrays, meta in event_chunks(item['path'], branches=input_branches(releases) +
                                                [m['score_branch'] for _, m in releases]):
                    scores, _ = score_chunk(arrays, meta, releases, models)
                    for branch, values in scores.items():
                        require(same_scores(values, arrays[branch]), 'existing saved score differs: ' + branch)
                require(file_hash(item['path']) == item['sha256'], 'reused group input changed')
            results.append(dict(group=group, status='REUSED', version=items[0]['version'])); continue
        base_version = items[0]['version']
        require(all(i['version'] == base_version for i in items), 'group input version differs')
        recipe = dict(schema='PhotonIDPublicationRecipeV1',
                      models={m['name']: m['release_sha256'] for m in manifests})
        identity = manifests[0]['release_sha256'] if len(manifests) == 1 else digest(recipe)
        root = Path(candidate_root) / digest(group) if candidate_root else Path(group) / '.photonid_candidates'
        work = root / identity / base_version
        work.mkdir(parents=True, exist_ok=True)
        with (work / '.scoring.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                # Every candidate is proven when written and bound to its receipt by
                # digest; publication rehashes each file before the pointer moves.
                folders, selected = zip(*missing)
                models = [predictors(path, m) for path, m in missing]
                records = [augment_file(i, folders, selected, work / i['relative_path'], code_root=code_root,
                                        model_predictors=models, progress_identity=progress_identity) for i in items]
                if candidate_root:
                    results.append(dict(group=group, status='CANDIDATES_READY', base_version=base_version,
                                        group_id=items[0]['group_id'], recipe=recipe, files=records))
                    continue
                version = json.loads(str(run(eventtree().publishGroup, json.dumps(dict(
                    group=str(group), files=records, group_id=items[0]['group_id'],
                    expected_sources=[i['logical_id'] for i in items],
                    recipe=recipe,
                    base_version=base_version, stage=3)))))
            except Exception as exc:
                progress(work / 'progress.json', 'blocked', 0, len(items), models=[m['name'] for m in manifests],
                         group=str(group), unit='logical-source files', latest_error=str(exc))
                raise
        results.append(dict(group=group, status='PASS', version=version['version'], files=len(items)))
    return results
