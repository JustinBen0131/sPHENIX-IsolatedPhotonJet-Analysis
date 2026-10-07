"""Inputs: source-bound H70 training rows from completed eventTree groups.

Reuse TreeProduction's publication manifests and photon feature map through its
installed ROOT-only eventTree library. A bounded source projection supplies
both training ranges; cached rows build disk-backed matrices without another
tree pass. SQLite checks identities and assigns an event-grouped
split once for every learner. Isolation never enters the labels. Weights
reproduce the existing class balance and spline ET/eta flattening.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sqlite3

import awkward as ak
import numpy as np
import yaml

from .Models import SYSTEMS, WorkProgress, digest, file_hash, read_json, require, save_json

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


def published_population(group, manifest, system):
    """Read the producer's pinned SIM roster without opening any eventTree.

    Older publications without a saved roster retain metadata inspection.
    Actual payload, metadata and event counts are checked during projection.
    """
    pin = manifest.get('recipe', {}).get('population_manifest_sha256')
    if manifest['stage'] != 2 or not pin:
        return None
    group = Path(group)
    candidates = [group.parent / 'population.json',
                  group.parent / 'populations' / (group.name + '.sources.tsv')]
    available = [path for path in candidates if path.is_file()]
    if not available:
        return None
    matches = [path for path in available if file_hash(path) == pin]
    require(len(matches) == 1, 'saved producer population hash is missing or ambiguous')
    path = matches[0]
    identity_keys = ('logical_source_id', 'sample', 'period', 'si_di_role',
                     'run', 'segment', 'source_file_ordinal')
    if path.suffix == '.tsv':
        with path.open(newline='') as stream:
            table = csv.DictReader(stream, delimiter='\t')
            require(table.fieldnames and len(table.fieldnames) == len(set(table.fieldnames)) and
                    set(identity_keys + ('event_count',)) <= set(table.fieldnames),
                    'producer source table identity columns are incomplete')
            sources = list(table)
        for source in sources:
            for key in ('run', 'segment', 'source_file_ordinal', 'event_count'):
                source[key] = int(source[key])
    else:
        population = read_json(path)
        require(population.get('schema_version') == 2 and isinstance(population.get('sources'), list),
                'unsupported saved producer population')
        sources = population['sources']
    sources = [source for source in sources if source.get('sample') == group.name]
    if not sources:
        return None
    expected_role = 'EMBEDDED' if system == 'AuAu' else 'SI'
    if any(source.get('si_di_role') != expected_role for source in sources):
        return None
    roster = {}
    for source in sources:
        require(set(identity_keys) <= set(source) and
                all(type(source[key]) is int and source[key] >= 0
                    for key in ('run', 'segment', 'source_file_ordinal', 'event_count')),
                'producer population identity or event count is invalid')
        identity = source['logical_source_id']
        require(identity not in roster, 'duplicate producer population source')
        roster[identity] = dict(source={key: source[key] for key in identity_keys},
                                entries=source['event_count'])
    require(set(roster) == {item['logical_id'] for item in manifest['files']} and
            len(roster) == len(manifest['files']) and
            manifest['population']['status'] == 'PASS' and
            sum(item['entries'] for item in roster.values()) == manifest['population']['entries'],
            'producer population differs from published sources or event total')
    require(file_hash(path) == pin, 'producer population changed during inspection')
    return roster, dict(path=str(path.resolve()), sha256=pin)


def published_inputs(groups, system, *, progress_path=None, progress_identity=None, verify_payload=True):
    """Resolve each CURRENT once and verify its exact logical-source roster."""
    require(groups, 'no completed dataset groups configured; use --group or config datasets')
    resolved_groups = []
    for group in groups:
        resolved = json.loads(str(run(eventtree().resolveGroup, str(group))))
        version, manifest = Path(resolved['directory']), resolved['manifest']
        require(manifest['stage'] in (2, 3), 'training/scoring requires a completed group')
        resolved_groups.append((group, version, manifest))
    total = sum(len(m['files']) for _, _, m in resolved_groups)
    tracker = WorkProgress(progress_path, total, 'input files', 'input-inspection', progress_identity) if progress_path else None
    inputs, seen = [], set()
    if tracker: tracker.write('inspect-inputs', input_events=0)
    try:
        for group, version, manifest in resolved_groups:
            population = None if verify_payload else published_population(group, manifest, system)
            group_path = str(Path(group).resolve())
            for item in manifest['files']:
                relative = Path(item['relative_path'])
                require(not relative.is_absolute() and '..' not in relative.parts, 'invalid group file path')
                path = version / relative
                if population is None:
                    path = path.resolve()
                    require(path.is_relative_to(version) and path.is_file(), 'invalid group file path')
                if verify_payload:
                    require(file_hash(path) == item['sha256'], 'published input hash differs: ' + str(path))
                require(item['logical_id'] not in seen, 'duplicate logical source across selected groups')
                seen.add(item['logical_id'])
                if population is None:
                    meta = metadata(path)
                    require(meta['schema']['system'] == SYSTEMS[system], 'mixed collision systems')
                    require(meta['source']['logical_source_id'] == item['logical_id'], 'logical identity mismatch')
                    record = dict(source=meta['source'], entries=source_entries(meta))
                else:
                    record = dict(population[0][item['logical_id']], producer_population=population[1])
                inputs.append(dict(path=str(path), sha256=item['sha256'], logical_id=item['logical_id'],
                                   relative_path=item['relative_path'], group=group_path,
                                   group_id=manifest['group_id'], version=manifest['version'],
                                   **record))
                if tracker:
                    tracker.write('inspect-inputs', len(inputs), source_logical_id=item['logical_id'],
                                  input_events=tracker.extra['input_events'] + inputs[-1]['entries'])
        if tracker: tracker.write('complete', total, status='PASS')
    except BaseException as error:
        if tracker: tracker.fail(error)
        raise
    return sorted(inputs, key=lambda item: item['logical_id'])


def source_entries(meta):
    """Count saved eventTree rows from exact contiguous consolidation metadata."""
    require('consolidation' in meta, 'input is not a completed logical-source file')
    total = meta['schema']['entries']
    require(type(total) is int and total >= 0, 'invalid source event count')
    previous = 0
    for chunk in meta['consolidation']['chunks']:
        begin, end = chunk['entry_begin'], chunk['entry_end']
        require(type(begin) is int and type(end) is int and begin == previous and end >= begin,
                'source event metadata overlaps or has gaps')
        previous = end
    require(previous == total == meta['consolidation']['expected_events'], 'source event metadata is incomplete')
    return total


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
    """Resolve a SIM sample to its recipe role; refuse DATA, DI and ambiguous roles."""
    sample = str(source.get('sample', '')).lower()
    matches = [role for role, samples in recipe['sources'].items() if sample in samples]
    require(len(matches) == 1, 'SIM sample not bound to the training recipe: ' + sample)
    require(source.get('si_di_role') in ('EMBEDDED', 'SI'), 'nominal training excludes DATA/DI')
    if recipe['system'] == 'AuAu':
        require(source['si_di_role'] == 'EMBEDDED', 'AuAu training requires embedded SIM')
    return matches[0]


def training_rows(item, recipe, *, event_progress=None):
    """Read only training witnesses; exclude classes or events outside the recipe."""
    BRANCH = feature_branches()
    role = source_role(item['source'], recipe)
    branches = sorted(set([BRANCH[n] for n in recipe['features']] + [
        'pho_et','pho_eta','cent_value','vtx_z','vtx_valid','evt_record','evt_nodes_ready',
        'evt_terminal_status','evt_src_entry','evt_run','evt_seq','evt_embedded_mb_decision',
        'shv_features','shv_finite_feature_state','pho_dom_barcode','pho_dom_state','pho_dom_pid',
        'truthpho_generator_barcode','truthpho_prompt_class']))
    if recipe['system'] == 'AuAu':
        branches = sorted(set(branches + [recipe['training']['minimum_bias_branch']]))
    training_recipe = dict(recipe, domain=recipe.get('training_domain', recipe['domain']))
    for arrays, meta in event_chunks(item['path'], branches=branches,
                                    step_size=100 if event_progress else '16 MB'):
        x, counts, et, cent, inside, valid = feature_rows(arrays, meta, training_recipe)
        y = labels(arrays, role)
        if recipe['system'] == 'AuAu':
            witness = recipe['training']['minimum_bias_branch']
            inside &= np.repeat(ak.to_numpy(arrays[witness]) == recipe['training']['minimum_bias_pass'], counts)
        keep = inside & valid & (y >= 0)
        entries = np.repeat(ak.to_numpy(arrays['evt_src_entry']), counts)
        photon_indices = np.concatenate([np.arange(c) for c in counts]) if counts.sum() else np.array([], np.int64)
        keys = np.repeat(event_keys(meta['source'], ak.to_numpy(arrays['evt_run']),
                                    ak.to_numpy(arrays['evt_seq'])), counts)
        yield x[keep], y[keep], et[keep], cent[keep], entries[keep], photon_indices[keep], keys[keep]
        if event_progress: event_progress(len(arrays), int(keep.sum()))


def event_keys(source, runs, sequences):
    """Historical 128-bit event identity from lane/sample/run/segment/sequence.

    These original provenance fields survive publication and file reordering.
    The producer sequence is one-based; source-entry indices remain zero-based.
    Never substitute the publication hash or the current file-list position.
    """
    require(all(k in source for k in ('lane', 'sample', 'segment')),
            'training source lacks historical event identity fields')
    require(len(runs) == len(sequences), 'event identity columns are not aligned')
    keys = []
    for run, sequence in zip(runs, sequences):
        require(int(sequence) >= 1, 'invalid one-based producer event sequence')
        text = '|'.join((str(source['lane']), str(source['sample']), str(int(run)),
                         str(int(source['segment'])), str(int(sequence))))
        key = hashlib.sha256(text.encode()).hexdigest()[:32]
        keys.append(key if int(key, 16) else '0' * 31 + '1')
    return np.asarray(keys, dtype='S32')


# Source inspection and bounded projection. Raw rows carry no training weights;
# those are derived from each complete range, never independently per shard.
def catalogue_body(inputs, system):
    require(system in SYSTEMS and inputs, 'empty or unsupported training catalogue')
    sources = sorted(inputs, key=lambda item: item['logical_id'])
    require(len({item['logical_id'] for item in sources}) == len(sources), 'duplicate catalogue source')
    for item in sources:
        require(Path(item['path']).is_absolute() and isinstance(item['sha256'], str) and
                len(item['sha256']) == 64 and all(c in '0123456789abcdef' for c in item['sha256']),
                'invalid catalogue source path/hash')
        require(type(item['entries']) is int and item['entries'] >= 0, 'invalid catalogue event count')
        require(item['source']['si_di_role'] in ('EMBEDDED', 'SI'), 'training catalogue excludes DATA/DI')
    return dict(schema='PhotonIDTrainingInputsV1', system=system, stage=2, sources=sources,
                files=len(sources), input_events=sum(item['entries'] for item in sources))


def save_catalogue(inputs, system, root):
    """Save inspected stage-2 sources as one immutable analyst-facing YAML file."""
    body = catalogue_body(inputs, system); pin = digest(body)
    path = Path(root) / '.inputs' / pin / 'inputs.yaml'
    path.parent.mkdir(parents=True, exist_ok=True)
    records = path.with_name('inputs.json')
    save_json(records, body)
    value = dict(schema=body['schema'], system=system, stage=2, files=body['files'], input_events=body['input_events'],
                 inputs_sha256=pin, sources_file=records.name, sources_sha256=file_hash(records))
    if path.exists():
        require(yaml.safe_load(path.read_text()) == value, 'training catalogue changed')
    else:
        with path.open('x') as stream: stream.write(yaml.safe_dump(value, sort_keys=False))
    return path


def load_catalogue(path, system):
    """Verify the roster itself without reopening its TG source files."""
    path = Path(path)
    value = yaml.safe_load(path.read_text())
    require(isinstance(value, dict) and value.get('schema') == 'PhotonIDTrainingInputsV1' and
            value.get('system') == system and value.get('stage') == 2, 'unsupported training catalogue')
    require(value.get('sources_file') == 'inputs.json', 'catalogue source record path differs')
    records = path.with_name('inputs.json')
    require(file_hash(records) == value['sources_sha256'], 'catalogue source records changed')
    record_body = read_json(records)
    body = catalogue_body(record_body['sources'], system)
    require(record_body == body and value == dict(schema=body['schema'], system=system, stage=2,
            files=body['files'], input_events=body['input_events'], inputs_sha256=digest(body),
            sources_file=records.name, sources_sha256=file_hash(records)), 'training catalogue hash/counts differ')
    return body['sources']


def projection_dtype(recipe):
    """Little-endian fixed records; entry/photon are original source identities."""
    return np.dtype([('x', '<f4', (len(recipe['features']),)), ('y', 'i1'), ('et', '<f8'),
                     ('centrality', '<f8'), ('entry', '<i8'), ('photon', '<i8'), ('event_key', 'S32')])


def projection_key(item, recipe):
    return digest(dict(source=item, definition=population_definition(recipe), code_sha256=file_hash(__file__)))


def projection_directory(item, recipe, root):
    pin = projection_key(item, recipe)
    return Path(root) / pin[:2] / pin


def projection_witness(folder, *, index_file='projection.json'):
    result = []
    for name in ('rows.bin', index_file):
        stat = (Path(folder) / name).stat()
        result.append((stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))
    return result


def projection_shards(inputs, max_events=100000):
    """Contiguous, duplicate-free whole-source shards with exact event totals."""
    require(type(max_events) is int and 0 < max_events <= 100000, 'projection event bound must be 1..100000')
    require(inputs == sorted(inputs, key=lambda item: item['logical_id']) and
            len({item['logical_id'] for item in inputs}) == len(inputs), 'projection roster must be sorted and unique')
    shards = []; begin = 0; events = 0
    for index, item in enumerate(inputs):
        n = item['entries']
        require(type(n) is int and 0 <= n <= max_events,
                'one source exceeds the projection event bound; use qualified source-range partitioning')
        if index > begin and events + n > max_events:
            shards.append(dict(index=len(shards), begin=begin, end=index, events=events))
            begin = index; events = 0
        events += n
    if inputs: shards.append(dict(index=len(shards), begin=begin, end=len(inputs), events=events))
    return shards


def load_projection(item, recipe, root):
    """Admit only a sealed, byte-verified raw-row cache for this exact source."""
    folder = projection_directory(item, recipe, root)
    witness = projection_witness(folder)
    receipt = read_json(folder / 'projection.json')
    body = {k: v for k, v in receipt.items() if k != 'projection_sha256'}
    require(receipt['projection_sha256'] == digest(body) and receipt['status'] == 'PASS' and
            receipt['schema'] == 'PhotonIDSourceProjectionV1', 'projection receipt changed or incomplete')
    require(receipt['source'] == item and receipt['definition'] == population_definition(recipe) and
            receipt['cache_identity'] == projection_key(item, recipe), 'projection source/definition differs')
    dtype = projection_dtype(recipe)
    require(digest(receipt['record_dtype']) == digest(dtype.descr) and
            type(receipt['rows']) is int and receipt['rows'] >= 0 and
            receipt['source_events'] == item['entries'], 'projection counts/layout differ')
    path = folder / 'rows.bin'
    require(path.stat().st_size == receipt['rows'] * dtype.itemsize and
            file_hash(path) == receipt['rows_sha256'], 'projection rows changed')
    require(projection_witness(folder) == witness, 'projection changed during verification')
    return dict(folder=folder, recipe=recipe, receipt=receipt, witness=witness)


def project_source(item, recipe, root, *, progress_identity=None, event_progress=None):
    """Read one complete TG source once; append bounded chunks to an immutable cache."""
    require(type(item['entries']) is int and 0 <= item['entries'] <= 100000, 'source projection exceeds its event bound')
    folder = projection_directory(item, recipe, root); folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.projection.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (folder / 'projection.json').exists(): return load_projection(item, recipe, root)
        require(not (folder / 'rows.bin').exists(), 'incomplete projection: retain evidence and use a fresh cache root')
        tracker = WorkProgress(folder / 'progress.json', item['entries'], 'tree events',
                               item['logical_id'], progress_identity)
        tracker.write('project-verify', eligible_photons=0)
        def witness():
            stat = Path(item['path']).stat()
            return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns
        try:
            original = witness()
            require(file_hash(item['path']) == item['sha256'], 'projection input changed')
            meta = metadata(item['path'])
            source_matches = (all(meta['source'].get(key) == value for key, value in item['source'].items())
                              if 'producer_population' in item else meta['source'] == item['source'])
            require(source_matches and meta['schema']['system'] == SYSTEMS[recipe['system']] and
                    meta['source']['logical_source_id'] == item['logical_id'] and source_entries(meta) == item['entries'],
                    'projection source metadata differs')
            dtype = projection_dtype(recipe); rows = 0; classes = [0, 0]
            def observed(events, selected):
                tracker.write('project-source', tracker.processed + events, eligible_photons=rows)
                if event_progress: event_progress(events, selected)
            with (folder / 'rows.bin').open('xb') as stream:
                for x, y, et, cent, entries, photons, keys in training_rows(item, recipe, event_progress=observed):
                    records = np.empty(len(y), dtype=dtype)
                    for key, values in zip(dtype.names, (x, y, et, cent, entries, photons, keys)): records[key] = values
                    records.tofile(stream); rows += len(y)
                    for klass in (0, 1): classes[klass] += int((y == klass).sum())
                stream.flush(); os.fsync(stream.fileno())
            require(tracker.processed == item['entries'] and witness() == original,
                    'projection event count or input changed during reading')
            tracker.write('project-seal')
            body = dict(schema='PhotonIDSourceProjectionV1', status='PASS', cache_identity=projection_key(item, recipe),
                        source=item, source_events=item['entries'], definition=population_definition(recipe),
                        rows=rows, class_counts=classes, record_dtype=dtype.descr,
                        rows_sha256=file_hash(folder / 'rows.bin'), isolation_used=False, weights_used=False)
            save_json(folder / 'projection.json', dict(body, projection_sha256=digest(body)))
            tracker.write('complete', tracker.total, status='PASS')
        except BaseException as error:
            tracker.fail(error)
            raise
    return load_projection(item, recipe, root)


def project_shard(inputs, recipe, root, shard, *, progress_identity=None):
    """Execute one admitted whole-source shard; verified parts are deterministic skips."""
    require(all(type(shard[k]) is int for k in ('index', 'begin', 'end', 'events')) and
            shard['index'] >= 0 and 0 <= shard['begin'] < shard['end'] <= len(inputs), 'invalid projection shard interval')
    selected = inputs[shard['begin']:shard['end']]
    require(selected and sum(item['entries'] for item in selected) == shard['events'] <= 100000,
            'projection shard event bounds differ')
    pin = digest(dict(inputs=inputs, definition=population_definition(recipe), shard=shard, code=file_hash(__file__)))
    folder = Path(root) / '.shards' / pin
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.shard.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        cached = {item['logical_id'] for item in selected
                  if (projection_directory(item, recipe, root) / 'projection.json').exists()}
        skips = sum(item['entries'] for item in selected if item['logical_id'] in cached)
        tracker = WorkProgress(folder / 'progress.json', shard['events'], 'source events',
                               'projection-shard-' + pin, progress_identity, skip_total=skips)
        tracker.write('project-start', sources=len(selected), eligible_photons=0)
        parts = {}
        try:
            for item in selected:
                tracker.write('project-verify', source_logical_id=item['logical_id'])
                if item['logical_id'] in cached:
                    part = load_projection(item, recipe, root)
                    tracker.write('project-reuse', tracker.processed + item['entries'],
                                  skip_processed=tracker.skipped + item['entries'],
                                  eligible_photons=tracker.extra['eligible_photons'] + part['receipt']['rows'])
                else:
                    def observed(events, selected_rows):
                        tracker.write('project-source', tracker.processed + events,
                                      eligible_photons=tracker.extra['eligible_photons'] + selected_rows)
                    part = project_source(item, recipe, root, progress_identity=progress_identity, event_progress=observed)
                parts[item['logical_id']] = part['receipt']['projection_sha256']
            require(tracker.processed == tracker.total and tracker.skipped == skips, 'projection shard coverage differs')
            save_json(folder / 'shard.json', dict(schema='PhotonIDProjectionShardV1', status='PASS',
                      shard=shard, shard_sha256=pin, parts=parts, input_events=shard['events']))
            tracker.write('complete', tracker.total, status='PASS')
        except BaseException as error:
            tracker.fail(error)
            raise
    return dict(status='REUSED' if len(cached) == len(selected) else 'PROJECTED', shard=shard,
                receipt=str(folder / 'shard.json'), progress=str(folder / 'progress.json'))


def require_projection_domain(projected, recipe):
    """Raw rows may supply narrower kinematics and either split/weight recipe."""
    requested = population_definition(recipe)
    for key in ('features', 'sources'):
        require(projected[key] == requested[key], 'projection definition differs: ' + key)
    if recipe['system'] == 'AuAu':
        for key in ('minimum_bias_branch', 'minimum_bias_pass'):
            require(projected['training'][key] == requested['training'][key], 'projection MB selection differs')
    wider = projected['domain']; narrower = requested['domain']
    require(set(wider) == set(narrower) and all(wider[k] == narrower[k] for k in wider
            if k not in ('et_min', 'et_max_exclusive', 'eta_abs_max', 'vertex_abs_max_cm')) and
            wider['eta_abs_max'] >= narrower['eta_abs_max'] and
            wider['vertex_abs_max_cm'] >= narrower['vertex_abs_max_cm'] and
            wider['et_min'] <= narrower['et_min'] < narrower['et_max_exclusive'] <= wider['et_max_exclusive'],
            'projection does not contain the requested ET domain')


def cached_rows(item, recipe, projection, tracker, phase):
    """Filter memory-mapped cached rows without reading or hashing TG trees."""
    receipt = projection['receipt']; require_projection_domain(receipt['definition'], recipe)
    index_file = projection.get('index_file', 'projection.json')
    require(projection_witness(projection['folder'], index_file=index_file) == projection['witness'], 'verified projection changed')
    dtype = projection_dtype(recipe); n = receipt['rows']
    data = np.memmap(projection['folder'] / 'rows.bin', dtype=dtype, mode='r', shape=(n,),
                     offset=projection.get('byte_offset', 0)) if n else np.empty(0, dtype)
    tracker.write(phase, source_logical_id=item['logical_id'], source_processed=0,
                  source_total=n, source_remaining=n)
    domain = recipe.get('training_domain', recipe['domain'])
    for begin in range(0, n, 65536):
        end = min(n, begin + 65536); records = data[begin:end]
        keep = (records['et'] >= domain['et_min']) & (records['et'] < domain['et_max_exclusive'])
        keep &= np.abs(records['x'][:, recipe['features'].index('cluster_Eta')]) < domain['eta_abs_max']
        keep &= np.abs(records['x'][:, recipe['features'].index('vertexz')]) < domain['vertex_abs_max_cm']
        selected = records[keep]
        yield tuple(selected[key] for key in dtype.names)
        field = 'eligible_photons' if phase == 'prepare-count' else 'matrix_rows_written'
        tracker.write(phase, tracker.processed + end - begin, source_processed=end, source_remaining=n - end,
                      **{field: tracker.extra.get(field, 0) + len(selected)})
    require(projection_witness(projection['folder'], index_file=index_file) == projection['witness'], 'projection changed during preparation')


@dataclass
class Population:
    """One prepared photon-row population with aligned features and witnesses.

    Row i has the same meaning in x, labels, kinematics, weights and train.
    Preparation may back these arrays with read-only memory-mapped files;
    receipt binds their input population and event-grouped split.
    """
    x: np.ndarray
    y: np.ndarray
    et: np.ndarray
    centrality: np.ndarray
    weights: np.ndarray
    train: np.ndarray
    receipt: dict
    # Historical pp: 0=train, 1=validation, 2=test. Canonical: 0=train, 2=holdout.
    partition: np.ndarray | None = None
    split_rank: np.ndarray | None = None

    @property
    def holdout(self):
        return ~self.train if self.partition is None else self.partition == 2


def inverse_pdf(values, *, fixed_range=None, cap=None, progress=None):
    """Exact PPG12 twenty-bin spline weights, with bounded observable passes."""
    from scipy.interpolate import UnivariateSpline
    values = np.asarray(values, np.float64)
    n = len(values); require(n >= 20, 'insufficient finite rows for exact weights')
    def batches(phase):
        if progress: progress(phase, 0, n)
        for begin in range(0, n, 65536):
            end = min(n, begin + 65536)
            yield begin, end
            if progress: progress(phase, end, n)
    lo, hi = float('inf'), float('-inf')
    for begin, end in batches('scan'):
        part = values[begin:end]
        require(np.isfinite(part).all(), 'insufficient finite rows for exact weights')
        lo, hi = min(lo, float(part.min())), max(hi, float(part.max()))
    lo, hi = fixed_range or (lo, hi)
    require(hi > lo, 'zero-width weighting population')
    edges = np.linspace(lo, hi, 21); counts = np.zeros(20, np.int64)
    for begin, end in batches('histogram'):
        counts += np.histogram(values[begin:end], bins=edges)[0]
    # Preserve NumPy density=True's division order and the historical spline.
    histogram = counts / np.diff(edges) / counts.sum()
    spline = UnivariateSpline((edges[:-1] + edges[1:]) / 2, histogram * 20., s=0.)
    factors = np.empty(n, np.float64)
    for begin, end in batches('spline'):
        factors[begin:end] = 1. / np.clip(spline(values[begin:end]), 1e-3, None)
        if cap is not None: factors[begin:end] = np.minimum(factors[begin:end], cap)
    mean = factors.mean()
    for begin, end in batches('normalize'):
        factors[begin:end] /= mean
    return factors


def canonical_weights(x, y, features, *, progress=None):
    """Class balance times per-class eta/ET flattening, before the event split.

    These are training weights; generator normalization and analysis event
    weights do not enter the learner's objective.
    """
    weights = np.ones(len(y), np.float64)
    for klass in (0, 1):
        mask = y == klass; n = int(mask.sum())
        require(n >= 20, 'both training classes need at least twenty eligible rows')
        weights[mask] = len(y) / (2. * n)
        for name, settings in [('cluster_Eta', dict(fixed_range=(-.7, .7))), ('cluster_Et', dict(cap=800.))]:
            def observed(phase, processed, total):
                if progress: progress(f'prepare-weights-class{klass}-{name}-{phase}', processed, total)
            weights[mask] *= inverse_pdf(x[mask, features.index(name)], progress=observed, **settings)
    require(np.isfinite(weights).all() and (weights > 0).all(), 'invalid canonical training weights')
    return weights


def population_definition(recipe):
    """Settings that change eligible rows, their weights or event assignments."""
    definition = {key: recipe[key] for key in ('features', 'training', 'sources', 'weights')}
    definition['domain'] = recipe.get('training_domain', recipe['domain'])
    return definition


def population_identity(inputs, recipe, split_reference=None):
    """Learner and routing variations share preparation, weights and event split.

    A wider population also binds the nominal population identity whose event
    assignments it inherits. The reference does not change the model recipe.
    """
    identity = dict(inputs=inputs, **population_definition(recipe),
                    preparation_code=file_hash(__file__))
    if split_reference is not None:
        identity['split_reference'] = split_reference
    return digest(identity)


def reference_split(work, inputs, recipe):
    """Verify a narrower, completed population before borrowing its event split."""
    work = Path(work).resolve()
    receipt = read_json(work / 'population.json')
    require(receipt['sources'] == inputs, 'split reference has different logical sources')
    definition = population_definition(recipe)
    for key in ('features', 'training', 'sources', 'weights'):
        require(receipt['definition'][key] == definition[key], 'split reference differs: ' + key)
    narrower = receipt['definition']['domain']; wider = definition['domain']
    require(set(narrower) == set(wider) and
            all(narrower[k] == wider[k] for k in wider if k not in ('et_min', 'et_max_exclusive')) and
            wider['et_min'] <= narrower['et_min'] < narrower['et_max_exclusive'] <= wider['et_max_exclusive'],
            'split reference must cover the same selection with a narrower ET range')
    require(population_identity(inputs, dict(recipe, training_domain=narrower)) == receipt['input_recipe_sha256'],
            'split reference belongs to different preparation code or settings')
    require(file_hash(work / 'events.sqlite') == receipt['event_split_sha256'],
            'reference event split changed')
    return work, receipt


def assign_events(db, event_count, reference=None):
    """Freeze nominal rank-based 90/10 membership; extend it without moving events.

    New wider-range events use the nominal last-held-out (hash,event_key)
    cutoff. Consequently the wider event fraction is recorded, not forced to
    exactly 10 percent at the cost of changing common event assignments.
    """
    if reference is None:
        count = max(1, min(event_count - 1, round(event_count * .1)))
        db.execute('UPDATE events SET holdout=1 WHERE rowid IN '
                   '(SELECT rowid FROM events ORDER BY hash,event_key LIMIT ?)', (count,))
        db.commit()
        return count, None
    folder, receipt = reference
    db.execute('ATTACH DATABASE ? AS nominal', ((folder / 'events.sqlite').as_uri() + '?mode=ro',))
    shared = db.execute('SELECT COUNT(*) FROM nominal.events r JOIN events e '
                        'ON r.source=e.source AND r.entry=e.entry AND r.hash=e.hash AND r.event_key=e.event_key').fetchone()[0]
    require(shared == receipt['events'], 'wider population lost or changed a reference event')
    cutoff = db.execute('SELECT hash,event_key FROM nominal.events WHERE holdout=1 '
                        'ORDER BY hash DESC,event_key DESC LIMIT 1').fetchone()
    require(cutoff is not None, 'reference has no held-out events')
    db.execute('UPDATE events SET holdout = ((hash,event_key) <= (?,?))', cutoff)
    db.execute('UPDATE events SET holdout = (SELECT r.holdout FROM nominal.events r '
               'WHERE r.source=events.source AND r.entry=events.entry) '
               'WHERE EXISTS (SELECT 1 FROM nominal.events r '
               'WHERE r.source=events.source AND r.entry=events.entry)')
    differences = db.execute('SELECT COUNT(*) FROM nominal.events r JOIN events e '
                             'ON r.source=e.source AND r.entry=e.entry WHERE r.holdout<>e.holdout').fetchone()[0]
    require(differences == 0, 'common event split changed')
    count = db.execute('SELECT COUNT(*) FROM events WHERE holdout=1').fetchone()[0]
    db.commit(); db.execute('DETACH DATABASE nominal')
    return count, dict(input_recipe_sha256=receipt['input_recipe_sha256'],
                       event_split_sha256=receipt['event_split_sha256'], shared_events=shared,
                       additional_events=event_count - shared, common_assignment_differences=differences,
                       method='inherit_nominal_and_extend_hash_cutoff', cutoff=list(cutoff))


def ppg12_rows(data, inputs, recipe, tracker):
    """Historical per-sample low-ET flattening, with deterministic retained order.

    A published sample may contain many logical ROOT sources. Flatten the
    complete sample, as the old shapes_split_<sample>.txt loader did, rather
    than flattening each ROOT source independently. Occupied bins keep their
    minimum count; bins are right-closed as in pandas.cut(include_lowest=True).
    """
    settings = recipe['training']['background_flattening']; pieces = []; report = []
    for role in ('background', 'signal'):
        for sample in recipe['sources'][role]:
            source_ids = [i for i, item in enumerate(inputs) if str(item['source']['sample']).lower() == sample]
            rows = np.flatnonzero(np.isin(data['source_index'], source_ids))
            low = rows[data['et'][rows] < settings['et_threshold']]
            if role == 'background' and len(low):
                edges = np.linspace(float(data['et'][low].min()), settings['et_threshold'], settings['bins'] + 1)
                bins = np.searchsorted(edges, data['et'][low], side='left') - 1
                bins[data['et'][low] == edges[0]] = 0
                counts = np.bincount(bins, minlength=settings['bins'])
                minimum = int(counts[counts > 0].min()); rng = np.random.RandomState(settings['seed'])
                retained = [rng.choice(low[bins == b], minimum, replace=False)
                            for b in range(settings['bins']) if counts[b]]
                retained.append(rows[data['et'][rows] >= settings['et_threshold']])
                selected = np.concatenate(retained)
                report.append(dict(sample=sample, input_rows=len(rows), retained_rows=len(selected),
                                   low_et_input_rows=len(low), occupied_bin_minimum=minimum))
            else:
                selected = rows
            pieces.append(selected)
            tracker.stage('prepare-ppg12-flatten-' + sample, len(rows), len(rows), 'photons')
    indices = np.concatenate(pieces)
    require(len(np.unique(indices)) == len(indices), 'background flattening duplicated a photon')
    return indices, report


def row_partition(y, training):
    """The historical stratified 70/10/20 row split (test seed 42, val seed 43)."""
    from sklearn.model_selection import train_test_split
    rows = np.arange(len(y))
    trainval, test = train_test_split(rows, test_size=training['test_fraction'],
                                    random_state=training['seed'], stratify=y)
    train, validation = train_test_split(trainval,
        test_size=training['validation_fraction'] / (training['train_fraction'] + training['validation_fraction']),
        random_state=training['seed'] + 1, stratify=y[trainval])
    partition = np.full(len(y), 2, np.int8)
    partition[train] = 0; partition[validation] = 1
    order = np.concatenate((train, validation, test)); rank = np.empty(len(y), np.int64)
    rank[order] = np.arange(len(y))
    return partition, rank


def tracked_rows(item, recipe, tracker, phase):
    """Count every consumed event, including rejected and photon-free events."""
    seen = 0
    tracker.write(phase, source_logical_id=item['logical_id'], source_processed=0,
                  source_total=item['entries'], source_remaining=item['entries'])
    def observed(events, selected):
        nonlocal seen
        seen += events
        require(seen <= item['entries'], 'source event count exceeded its admitted total')
        field = 'eligible_photons' if phase == 'prepare-count' else 'matrix_rows_written'
        tracker.write(phase, tracker.processed + events, source_processed=seen,
                      source_remaining=item['entries'] - seen,
                      **{field: tracker.extra.get(field, 0) + selected})
    yield from training_rows(item, recipe, event_progress=observed)
    require(seen == item['entries'], 'source event count differs from its admitted total')


def prepare(inputs, recipe, work, *, split_reference=None, progress_identity=None, projections=None):
    """Create/reuse a disk-backed population with duplicate-free source keys.

    A wider-range preparation may inherit a verified nominal event split;
    class balance and ET/eta weights still use its own eligible photon rows.
    """
    reference = reference_split(split_reference, inputs, recipe) if split_reference is not None else None
    work = Path(work); work.mkdir(parents=True, exist_ok=True)
    pin = population_identity(inputs, recipe, reference[1]['input_recipe_sha256'] if reference else None)
    projection_pins = None
    if projections is not None:
        require(set(projections) == {item['logical_id'] for item in inputs}, 'projection coverage is incomplete or unexpected')
        for item in inputs:
            part = projections[item['logical_id']]
            require(part['receipt']['source'] == item, 'projection source population differs')
            require_projection_domain(part['receipt']['definition'], recipe)
        projection_pins = {key: value['receipt']['projection_sha256'] for key, value in projections.items()}
    done = work / 'population.json'
    if done.exists():
        receipt = read_json(done)
        require(receipt['input_recipe_sha256'] == pin, 'prepared population belongs to different inputs')
        require(receipt.get('projection_receipts') == projection_pins, 'prepared projection provenance differs')
        if projections is None:
            for item in inputs:
                require(file_hash(item['path']) == item['sha256'], 'training input changed')
        for filename, h in receipt['arrays'].items():
            require(file_hash(work / filename) == h, 'prepared matrix changed')
        require(file_hash(work / 'events.sqlite') == receipt['event_split_sha256'], 'prepared event split changed')
        return load_population(work, receipt)
    require(not (work / 'events.sqlite').exists(), 'incomplete population: retain evidence and use a fresh work directory')
    counted_inputs = [dict(item, entries=item.get('entries', None)) for item in inputs]
    for item in counted_inputs:
        if item['entries'] is None: item['entries'] = source_entries(metadata(item['path']))
        require(type(item['entries']) is int and item['entries'] >= 0, 'invalid admitted source event count')
    input_events = sum(item['entries'] for item in counted_inputs)
    visits = 2 * sum(part['receipt']['rows'] for part in projections.values()) if projections is not None else 2 * input_events
    tracker = WorkProgress(work / 'progress.json', visits,
                           'cached photon visits' if projections is not None else 'event visits', pin, progress_identity)
    tracker.write('prepare-start', input_events=input_events, input_files=len(inputs), eligible_photons=0)
    db = None
    total = 0; counts = [0, 0]
    try:
        db = sqlite3.connect(work / 'events.sqlite', uri=True)
        db.set_progress_handler(lambda: (tracker.beat(), 0)[1], 100000)
        db.execute('CREATE TABLE events(source TEXT,entry INTEGER,event_key TEXT NOT NULL UNIQUE,hash TEXT,holdout INTEGER DEFAULT 0,PRIMARY KEY(source,entry))')
        db.execute('CREATE TABLE photons(source TEXT,entry INTEGER,photon INTEGER,PRIMARY KEY(source,entry,photon))')
        def rows(item, phase):
            if projections is not None:
                return cached_rows(item, recipe, projections[item['logical_id']], tracker, phase)
            require(file_hash(item['path']) == item['sha256'], 'training input changed')
            return tracked_rows(item, recipe, tracker, phase)
        for item in counted_inputs:
            for x, y, et, cent, entries, photons, keys in rows(item, 'prepare-count'):
                for entry, photon, key in zip(entries, photons, keys):
                    key = bytes(key).decode('ascii')
                    require(len(key) == 32 and all(c in '0123456789abcdef' for c in key),
                            'invalid historical event identity')
                    hash_text = '|'.join(('event50', str(recipe['training']['seed']),
                                          recipe['training']['split_identity'], key))
                    h = hashlib.sha256(hash_text.encode()).hexdigest()[:8]
                    db.execute('INSERT INTO events(source,entry,event_key,hash) VALUES(?,?,?,?) '
                               'ON CONFLICT(source,entry) DO NOTHING', (item['logical_id'], int(entry), key, h))
                    require(db.execute('SELECT event_key FROM events WHERE source=? AND entry=?',
                                       (item['logical_id'], int(entry))).fetchone()[0] == key,
                            'one source entry has inconsistent event identities')
                    db.execute('INSERT INTO photons VALUES(?,?,?)', (item['logical_id'], int(entry), int(photon)))
                total += len(y)
                for klass in (0, 1): counts[klass] += int((y == klass).sum())
            db.commit()
        require(total > 0 and min(counts) >= 20, 'nominal training population is empty or missing a class')
        require(tracker.extra['eligible_photons'] == total, 'eligible photon telemetry differs from counted rows')
        tracker.write('prepare-event-split', phase_unit='eligible events')
        nevents = db.execute('SELECT COUNT(*) FROM events').fetchone()[0]
        historical = recipe['training'].get('split_method') == 'row_stratified'
        require(not historical or reference is None, 'row split cannot inherit an event split')
        nholdout, reference_receipt = (None, None) if historical else assign_events(db, nevents, reference)
        shapes = {'x': (total, len(recipe['features'])), 'y': (total,), 'et': (total,),
                  'centrality': (total,), 'weights': (total,), 'train': (total,), 'partition': (total,),
                  'source_index': (total,), 'entry': (total,), 'photon': (total,), 'event_key': (total,), 'split_rank': (total,)}
        types = {'x': np.float32, 'y': np.int8, 'et': np.float64, 'centrality': np.float64,
                 'weights': np.float64, 'train': np.bool_, 'partition': np.int8,
                 'source_index': np.int32, 'entry': np.int64, 'photon': np.int32, 'event_key': 'S32', 'split_rank': np.int64}
        data = {name: np.lib.format.open_memmap(work / (name + '.npy'), mode='w+',
                                              dtype=types[name], shape=shape) for name, shape in shapes.items()}
        offset = 0
        for source_index, item in enumerate(counted_inputs):
            for x, y, et, cent, entries, photons, keys in rows(item, 'prepare-fill'):
                n = len(y); sl = slice(offset, offset + n)
                data['x'][sl] = x; data['y'][sl] = y; data['et'][sl] = et; data['centrality'][sl] = cent
                data['source_index'][sl] = source_index; data['entry'][sl] = entries; data['photon'][sl] = photons
                data['event_key'][sl] = keys
                if not historical:
                    data['train'][sl] = [not db.execute('SELECT holdout FROM events WHERE source=? AND entry=?',
                                                      (item['logical_id'], int(e))).fetchone()[0] for e in entries]
                offset += n
            if projections is None:
                require(file_hash(item['path']) == item['sha256'], 'input changed during population preparation')
        require(offset == total, 'training population changed between passes')
        require(tracker.extra.get('matrix_rows_written', 0) == total, 'matrix photon telemetry differs from filled rows')
        flattening = None
        if historical:
            retained, flattening = ppg12_rows(data, counted_inputs, recipe, tracker)
            # Compact on writable analysis storage. Keep original source-entry/
            # photon identities; no tree re-read or scored-tree copy is needed.
            for name in data:
                source = data[name]; shape = (len(retained),) + source.shape[1:]
                target = np.lib.format.open_memmap(work / (name + '.selected.npy'), mode='w+', dtype=source.dtype, shape=shape)
                for begin in range(0, len(retained), 65536):
                    end = min(len(retained), begin + 65536)
                    target[begin:end] = source[retained[begin:end]]
                    tracker.stage('prepare-ppg12-compact-' + name, end, len(retained), 'photons')
                target.flush(); data[name] = target
            total = len(retained); counts = [int((data['y'] == k).sum()) for k in (0, 1)]
            partition, rank = row_partition(data['y'], recipe['training'])
            data['partition'][:] = partition; data['split_rank'][:] = rank
            data['train'][:] = data['partition'] == 0
        else:
            data['partition'][:] = np.where(data['train'], 0, 2)
            data['split_rank'][:] = np.arange(total)
        tracker.write('prepare-weights', phase_unit='eligible photons')
        data['weights'][:] = canonical_weights(data['x'], data['y'], recipe['features'],
            progress=lambda phase, done, count: tracker.stage(phase, done, count, 'photons'))
        if recipe['training'].get('weight_dtype') == 'float32':
            data['weights'][:] = data['weights'].astype(np.float32)
        tracker.write('prepare-seal', phase_unit='prepared arrays')
        for array in data.values(): array.flush()
        def sealed_hash(path):
            return file_hash(path, progress=lambda done, count:
                             tracker.stage('prepare-seal-' + Path(path).name, done, count, 'bytes'))
        files = {name: Path(array.filename).name for name, array in data.items()}
        array_hashes = {filename: sealed_hash(work / filename) for filename in files.values()}
        receipt = dict(schema='PhotonIDPopulationV1', input_recipe_sha256=pin, rows=total,
                       class_counts=counts, events=nevents, holdout_events=nholdout,
                       split_method='row_stratified_70_10_20' if historical else 'event_grouped_90_10',
                       partition_counts=[int((data['partition'] == k).sum()) for k in (0, 1, 2)],
                       array_files=files, background_flattening=flattening,
                       definition=population_definition(recipe), split_reference=reference_receipt,
                       event_split_sha256=sealed_hash(work / 'events.sqlite'),
                       event_key_convention='sha256_128(lane|sample|run|segment|producer_sequence)',
                       split_order='sha256_32(event50|seed|split_identity|event_key), then event_key',
                       arrays=array_hashes,
                       split_sha256=array_hashes[files['partition' if historical else 'train']], isolation_used=False,
                       event_weights_used=False, sources=inputs, input_events=input_events,
                       event_visits=0 if projections is not None else visits,
                       cached_photon_visits=visits if projections is not None else 0,
                       projection_receipts=projection_pins)
        save_json(done, receipt)
        tracker.write('complete', tracker.total, status='PASS')
    except BaseException as error:
        tracker.fail(error)
        raise
    finally:
        if db is not None: db.close()
    return load_population(work, receipt)


def load_population(work, receipt):
    """Open the prepared arrays without copying the complete photon matrix."""
    files = receipt.get('array_files', {})
    arrays = {n: np.load(Path(work) / files.get(n, n + '.npy'), mmap_mode='r')
              for n in ('x', 'y', 'et', 'centrality', 'weights', 'train')}
    partition = np.load(Path(work) / files['partition'], mmap_mode='r') if 'partition' in files else None
    rank = np.load(Path(work) / files['split_rank'], mmap_mode='r') if 'split_rank' in files else None
    return Population(**arrays, receipt=receipt, partition=partition, split_rank=rank)
