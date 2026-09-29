"""Consistency check of one eventTree v1 file on its own (no source part needed).

Checks the structure a reader relies on: equal-length groups, strides, counters, every event-local
index inside its target collection, empty upstream-rejected entries, absent branches empty, and integer
codes inside their dictionaries. Returns a report; `ok` is False on any violation.
"""
import json

import awkward as ak
import numpy as np
import uproot

from .schema import Schema

# (index branch, target group count branch or list of them, allows -1)
INDICES = (('cell_pho', ['pho_et'], False), ('isoc_pho', ['pho_et'], False), ('pair_pho', ['pho_et'], False),
           ('pho_assoc_truthpho', ['truthpho_pt'], True), ('plink_reco', ['pho_et'], True),
           ('plink_truth', ['truthpho_pt'], True))
FIXED_PER_PHOTON = {'shv': 7, 'isow': 2}


def load(path):
    with uproot.open(path) as f:
        meta = {key: json.loads(str(f['meta/' + key])) for key in ('schema', 'codes', 'collections', 'constants') if 'meta/' + key in f}
        return f['eventTree'].arrays(library='ak'), meta


# Coded fields that index a dictionary entry (reasons, origins, evidence and
# scoring bindings) use -1 for "no entry", the augmentation convention.
INDEX_SUFFIXES = ('_reason', '_reason_index', '_origin_index', '_evidence_index', '_scoring_node', '_scoring_evidence')


def tree_keys(path):
    """Every TTree key of the file with its cycle: a valid eventTree file holds exactly one, eventTree;1 or any single
    cycle (an autosave cycle left beside the final write would be a second tree for every reader that lists keys)."""
    with uproot.open(path) as f:
        return sorted(k for k, c in f.classnames(recursive=True).items() if c in ('TTree', 'TNtuple'))


def check(path, schema=None):
    with uproot.open(path) as f:
        consolidated = 'meta/consolidation' in f
    if consolidated:
        from .consolidate import check_consolidated
        return dict(check_consolidated(path),path=str(path))
    arrays, meta = load(path)
    report = check_arrays(arrays, meta, schema)
    trees = tree_keys(path)
    if len(trees) != 1 or not trees[0].startswith('eventTree;'):
        report['problems'] = list(report.get('problems', [])) + ['the file holds %d TTree keys, not one eventTree: %s' % (len(trees), trees)]
        report['ok'] = False
    return dict(report, path=str(path))


def check_arrays(arrays, meta, schema=None):
    schema = schema or Schema()
    problems, notes = [], []
    info = meta['schema']
    codes = meta.get('codes', {})
    collections = meta.get('collections', {'reco': []})
    n = len(arrays)
    if info.get('schema') != 'RJEventTreeV1':
        problems.append('meta/schema is not RJEventTreeV1')
    if info.get('schema_sha256') != schema.sha256:
        notes.append('written with schema %s; this checker has %s' % (info.get('schema_sha256', '?')[:12], schema.sha256[:12]))
    if n != info.get('entries'):
        problems.append('entry count %d differs from meta %s' % (n, info.get('entries')))
    record = ak.to_numpy(arrays['evt_record'])
    if not np.isin(record, (0, 1)).all():
        problems.append('evt_record outside {0, 1}')
    if int((record == 0).sum()) != info.get('retained_events'):
        problems.append('retained entries differ from meta')
    missing = sorted(set(schema.branches) - set(arrays.fields))
    if missing:
        problems.append('branches missing: %s' % missing[:10])
    counts = {}
    for name in arrays.fields:
        if name.startswith('n_'):
            continue
        if 'var' in str(ak.type(arrays[name])):
            counts[name] = ak.to_numpy(ak.num(arrays[name], axis=1))
    strides = info.get('strides', {})
    absent = set(info.get('absent_branches', []))
    # Equal-length groups (strided branches carry stride x the group length).
    by_group = {}
    for name, c in counts.items():
        branch = schema.branches.get(name)
        if branch is None or name in absent:
            continue
        by_group.setdefault(branch.group, []).append(name)
    group_len = {}
    for group, names in by_group.items():
        plain = [x for x in names if x not in strides]
        if not plain:
            continue
        reference = counts[plain[0]]
        group_len[group] = reference
        for name in plain[1:]:
            if not np.array_equal(counts[name], reference):
                problems.append('%s length differs from %s' % (name, plain[0]))
        for name in names:
            if name in strides and not np.array_equal(counts[name], reference * strides[name]):
                problems.append('%s is not %d x the %s length' % (name, strides[name], group))
    for group, per in FIXED_PER_PHOTON.items():
        if group in group_len and 'pho' in group_len and not np.array_equal(group_len[group], per * group_len['pho']):
            problems.append('%s is not %d per photon' % (group, per))
    # Counters equal the lengths they count.
    for name, counter in (info.get('counters') or {}).items():
        if counter in arrays.fields and name in counts and not np.array_equal(ak.to_numpy(arrays[counter]), counts[name]):
            problems.append('%s differs from the length of %s' % (counter, name))
    # Event-local indices stay inside their target collections.
    zero = np.zeros(n, np.int64)
    jets = sum((counts.get(e['name'] + '_pt', zero) for e in collections['reco']), zero)
    truthjets = sum((counts.get('truthjet%s_pt' % r, zero) for r in ('02', '03', '04')), zero)
    checks = list(INDICES) + [('pair_jet', None, False), ('jlink_reco', None, True), ('jlink_truth', None, True)]
    for name, target, negative in checks:
        if name not in counts or name in absent:
            continue
        size = jets if name in ('pair_jet', 'jlink_reco') else truthjets if name == 'jlink_truth' else counts[target[0]]
        values = arrays[name]
        limit = ak.broadcast_arrays(ak.from_numpy(np.asarray(size, dtype=np.int64)), values)[0]
        bad = ak.any((values < (-1 if negative else 0)) | (values >= limit), axis=1)
        if ak.any(bad):
            problems.append('%s out of range in %d entries' % (name, int(ak.sum(bad))))
    # Upstream-rejected entries carry no collection values.
    rejected = record == 1
    if rejected.any():
        for name, c in counts.items():
            if c[rejected].any():
                problems.append('%s has values on upstream-rejected entries' % name)
                break
    for name in absent:
        if name in counts and counts[name].any():
            problems.append('absent branch %s holds values' % name)
    for name, dictionary in codes.items():
        if name in arrays.fields:
            values = arrays[name]
            flat = ak.to_numpy(ak.flatten(values, axis=1) if values.ndim > 1 else values)
            # Index fields hold -1 for "no entry"; category codes never do.
            lowest = -1 if name.endswith(INDEX_SUFFIXES) else 0
            if len(flat) and (flat.min() < lowest or flat.max() >= len(dictionary)):
                problems.append('%s code outside its dictionary' % name)
    if problems:
        notes.append('photon-ID inputs not checked: the file has structure problems')
    else:
        problems += photon_id_inputs(arrays, meta, notes)
    return dict(ok=not problems, entries=n, problems=problems, notes=notes,
                encoding=info.get('encoding'), groups_checked=len(group_len))


def photon_id_inputs(arrays, meta, notes):
    """The photon-ID model inputs: one view per shower-view slot, the model's view with the frozen semantics and
    floor, its stored features equal to the named branches, and (Au+Au) an integer centrality."""
    from .photon_id import feature_parity, load
    constants = meta.get('constants') or {}
    system = meta['schema'].get('system')
    per_photon = ak.num(arrays['pho_et'], axis=1)
    n = int(ak.sum(per_photon))
    if not n or 'showerViews.definition_name' not in constants:
        return []
    spec = load()[system]['features']
    names = constants['showerViews.definition_name']
    if names.count(spec['view']) != 1:
        return ['shower view %s appears %d times in the view order' % (spec['view'], names.count(spec['view']))]
    position = names.index(spec['view'])
    problems = []
    if constants['showerViews.semantic_sha256'][position] != spec['view_semantic_sha256']:
        problems.append('%s view semantic hash differs from the frozen model contract' % spec['view'])
    if np.uint64(constants['showerViews.floor_gev'][position]['f64']).view(np.float64) != spec['floor_gev']:
        problems.append('%s view floor differs from the frozen model contract' % spec['view'])
    stride = meta['schema']['strides']['shv_features']
    features = ak.to_numpy(ak.flatten(arrays['shv_features'])).reshape(n, 7, stride)[:, position, :]
    if stride != spec['count']:
        return problems + ['%s view carries %d features, the model %d' % (spec['view'], stride, spec['count'])]

    def column(branch):
        if branch.startswith('shv_'):
            return ak.to_numpy(ak.flatten(arrays[branch])).reshape(n, 7)[:, position]
        if branch in ('vtx_z', 'cent_value'):
            return np.repeat(ak.to_numpy(arrays[branch]), ak.to_numpy(per_photon))
        return ak.to_numpy(ak.flatten(arrays[branch]))
    mismatches = feature_parity(spec['order'], column, features)
    if mismatches:
        problems.append('%d photons: %s features differ from the named shower branches' % (mismatches, spec['view']))
    if system == 'auau':
        centrality = column('cent_value')
        finite = centrality[np.isfinite(centrality)]
        if np.any(finite != np.round(finite)):
            problems.append('Au+Au centrality is not an integer bin for every photon')
    notes.append('photon-ID inputs checked on %d photons (%s view)' % (n, spec['view']))
    return problems
