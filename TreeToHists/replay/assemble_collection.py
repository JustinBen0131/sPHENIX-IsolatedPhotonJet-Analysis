"""Assemble already validated C++ replay products; never calculate physics.

The complete nominal manifest fixes the sample roster and parent hashes.
No ROOT bytes, source identities, selection or normalization are rewritten.
"""
from pathlib import Path
import argparse
import copy
import datetime as dt
import hashlib
import json
import os
import uproot
import yaml


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def publish(path, text):
    with path.open('x') as stream:
        stream.write(text)


def assemble(parent_path, plan_path, output):
    parent_hash = sha(parent_path)
    parent = yaml.safe_load(parent_path.read_text())
    plan = json.loads(plan_path.read_text())
    assert parent['schema'] == 'PhotonJetHistogramCollectionV1' and parent['status'] == 'VALIDATED'
    assert 'measurement_replay' not in parent and 'compact_of' not in parent
    recipe = plan['histograms']['measurement_replay']
    if 'selection_change' in recipe:
        assert recipe['selection_change']['parent_selection_sha256'] == parent['selection_sha256']
        assert recipe['selection_change']['effective_selection_sha256'] == plan['selection']['selection_sha256']
    else:
        assert parent['selection_sha256'] == plan['selection']['selection_sha256']
    assert recipe['admission'] == 'MECHANICAL_FIXTURE_ONLY'
    originals = parent['samples']
    if 'centrality_weight' in recipe:
        centrality = recipe['centrality_weight']
        assert centrality['collection_subset'] == 'all_auau_parent_samples'
        assert centrality['parent_collection']['sha256'] == parent_hash
        originals = [row for row in originals if row['system'] == 'AuAu']
        assert centrality['sample_families'] == {row['key']: row['family'] for row in originals}
        assert centrality['parent_packages'] == {row['key']: row['sha256'] for row in originals}
    parents = {row['key']: row['sha256'] for row in originals}
    assert len(parents) == len(originals)
    assert set(parents) == {row['key'] for row in plan['samples']} == set(recipe['sample_effects'])
    pending, rows, files = {}, [], []
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    for original in originals:
        key = original['key']
        assert key == original['system'] + '/' + original['name'] and '..' not in Path(key).parts
        root, side = output / (key + '.root'), output / (key + '.json')
        receipt = json.loads(side.read_text())
        with uproot.open(root) as source:
            meta = json.loads(str(source['meta/package']))
        assert receipt['sha256'] == sha(root)
        assert receipt['plan_sha256'] == meta['plan_sha256'] == plan['plan_sha256']
        assert meta['histograms'] == plan['histograms'] and meta['selection'] == plan['selection']
        assert meta['source_count'] == original['source_count'] and meta['coverage'] == original['coverage']
        assert meta['system'] == original['system'] and meta['name'] == original['name']
        assert meta['family'] == original['family'] and meta['package_kind'] == 'FINAL'
        assert meta.get('run') == original.get('run') and meta.get('run_packages') == original.get('run_packages')
        provenance = meta['replay_provenance']
        assert provenance['recipe'] == recipe and provenance['parent_package_sha256'] == parents[key]
        assert provenance['normalization_applied'] is False and provenance['nominal_diagnostics'] == 'parent_only'
        validation = receipt['validation']
        assert validation['status'] == 'PASS' and validation['scope'] == 'measurement_response_purity_inputs'
        row = copy.deepcopy(original)
        row.update(sha256=receipt['sha256'], event_records=meta['event_records'], histograms=meta['histograms_written'],
                   package_plan_sha256=meta['plan_sha256'], package_binding_sha256=meta['package_binding_sha256'], validated_at=now)
        rows.append(row)
        files.extend([dict(path=key + '.root', sha256=receipt['sha256']), dict(path=key + '.json', sha256=sha(side))])
    for name, value in [('selection.yaml', plan['selection']), ('histograms.yaml', plan['histograms'])]:
        pending[name] = yaml.safe_dump(value, sort_keys=True)
        files.append(dict(path=name, sha256=hashlib.sha256(pending[name].encode()).hexdigest()))
    recipe_hash = hashlib.sha256(json.dumps(recipe, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    manifest = dict(schema='PhotonJetHistogramCollectionV1', status='VALIDATED', created_at=now,
                    input_snapshot_sha256=plan['plan_sha256'], selection_sha256=plan['selection']['selection_sha256'], samples=rows, files=files,
                    measurement_replay=dict(schema='PhotonJetMeasurementReplayCollectionV1', recipe_sha256=recipe_hash,
                        admission=recipe['admission'], parent_collection_sha256=parent_hash, parent_packages=parents,
                        nominal_diagnostics='parent_collection'))
    if 'centrality_weight' in recipe:
        manifest['measurement_replay']['system_subset'] = 'AuAu'
    pending['collection.yaml'] = yaml.safe_dump(manifest, sort_keys=True)
    for name in pending:
        assert not (output / name).exists(), 'Preserve existing declarations: ' + name
    assert sha(parent_path) == parent_hash
    for name, text in pending.items():
        publish(output / name, text)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = assemble(args.parent, args.plan, args.out)
    print('ASSEMBLED_VALIDATED_REPLAY_PACKAGES', len(result['samples']))
