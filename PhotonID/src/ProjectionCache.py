"""Compact storage for the raw photon rows produced by Inputs.

Each bounded source shard becomes one binary row file and one immutable index.
Extraction uses temporary storage; the index retains every source receipt and
byte interval. Preparation verifies the bundle once, then reads its intervals
without reopening trees. This module changes storage, not photon selection.
"""
import fcntl
import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile

from . import Inputs
from .Models import WorkProgress, digest, file_hash, read_json, require, save_json

SCHEMA = 'PhotonIDProjectionBundleV1'


def projection_code(manifest=None):
    """Extraction provenance, independent of a later population consumer."""
    if manifest is None:
        return dict(projector_sha256=file_hash(Inputs.__file__), storage_sha256=file_hash(__file__))
    saved = read_json(manifest)
    require(saved.get('schema') == 'PhotonIDPublicCodeSnapshotV1', 'original extraction code manifest is invalid')
    result = {}
    for path, key in (('PhotonID/src/Inputs.py', 'projector_sha256'),
                      ('PhotonID/src/ProjectionCache.py', 'storage_sha256')):
        selected = [row['sha256'] for row in saved['members'] if row['path'] == path]
        require(len(selected) == 1 and re.fullmatch('[0-9a-f]{64}', selected[0]),
                'original extraction code manifest is incomplete')
        result[key] = selected[0]
    return result


def bundle_identity(sources, recipe, *, code=None):
    return digest(dict(sources=sources, definition=Inputs.population_definition(recipe),
                       **(code or projection_code())))


def load_bundle(folder, recipe, *, progress=None, code=None):
    """Verify all bytes and member receipts once; return source-aligned views."""
    folder = Path(folder)
    witness = Inputs.projection_witness(folder, index_file='bundle.json')
    receipt = read_json(folder / 'bundle.json')
    body = {k: v for k, v in receipt.items() if k != 'bundle_sha256'}
    require(receipt.get('schema') == SCHEMA and receipt.get('status') == 'PASS' and
            receipt.get('bundle_sha256') == digest(body), 'bundle receipt changed or incomplete')
    require(receipt['definition'] == Inputs.population_definition(recipe), 'bundle population differs')
    code = code or receipt.get('projection_code') or projection_code()
    require(set(code) == {'projector_sha256', 'storage_sha256'} and
            all(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) for value in code.values()) and
            receipt.get('projection_code', code) == code, 'bundle extraction code provenance differs')
    sources = [part['projection']['source'] for part in receipt['parts']]
    require(sources and sources == sorted(sources, key=lambda item: item['logical_id']) and
            len({s['logical_id'] for s in sources}) == len(sources), 'bundle sources must be sorted and unique')
    require(receipt['cache_identity'] == bundle_identity(sources, recipe, code=code), 'bundle source or code identity differs')
    dtype = Inputs.projection_dtype(recipe)
    require(digest(receipt['record_dtype']) == digest(dtype.descr) and
            receipt['source_events'] == sum(item['entries'] for item in sources) <= 100000,
            'bundle event count or layout differs')
    result = {}; offset = 0; rows = 0; aggregate = hashlib.sha256()
    with (folder / 'rows.bin').open('rb') as stream:
        for part in receipt['parts']:
            member = part['projection']; item = member['source']
            member_body = {k: v for k, v in member.items() if k != 'projection_sha256'}
            require(member['schema'] == 'PhotonIDSourceProjectionV1' and member['status'] == 'PASS' and
                    member['projection_sha256'] == digest(member_body) and
                    member['cache_identity'] == digest(dict(source=item, definition=Inputs.population_definition(recipe),
                                                           code_sha256=code['projector_sha256'])) and
                    member['definition'] == receipt['definition'], 'bundle member provenance differs')
            count = member['rows']
            require(type(count) is int and count >= 0 and member['source_events'] == item['entries'] and
                    digest(member['record_dtype']) == digest(dtype.descr), 'bundle member count or layout differs')
            length = count * dtype.itemsize
            require(part['byte_offset'] == offset and part['byte_length'] == length,
                    'bundle intervals overlap or have gaps')
            remaining = length; checksum = hashlib.sha256()
            while remaining:
                block = stream.read(min(4 * 1024**2, remaining))
                require(block, 'bundle rows truncated')
                checksum.update(block); aggregate.update(block); remaining -= len(block)
                if progress: progress()
            require(checksum.hexdigest() == member['rows_sha256'], 'bundle source rows changed')
            result[item['logical_id']] = dict(folder=folder, recipe=recipe, receipt=member,
                witness=witness, index_file='bundle.json', byte_offset=offset)
            offset += length; rows += count
        require(not stream.read(1) and aggregate.hexdigest() == receipt['rows_sha256'] and
                offset == receipt['bytes'] and rows == receipt['rows'], 'bundle total bytes or rows differ')
    require(Inputs.projection_witness(folder, index_file='bundle.json') == witness,
            'bundle changed during verification')
    return result


def load_campaign(inputs, recipe, root, *, progress_path=None, progress_identity=None, code_manifest=None):
    """Require complete coverage, admitting either packed or legacy source caches."""
    expected = {item['logical_id']: item for item in inputs}; result = {}
    code = projection_code(code_manifest) if code_manifest is not None else None
    tracker = WorkProgress(progress_path, len(expected), 'projection sources',
                           'projection-cache-' + recipe['system'].lower(), progress_identity) if progress_path else None
    verified_rows = 0
    if tracker: tracker.write('verify-projections', cached_photon_rows=0)
    try:
        for path in sorted((Path(root) / '.bundles').glob('*/*/bundle.json')):
            if tracker: tracker.beat()
            header = read_json(path)
            if header.get('definition') != Inputs.population_definition(recipe):
                continue
            ids = {part['projection']['source']['logical_id'] for part in header['parts']}
            if not ids.intersection(expected):
                continue
            for logical_id, part in load_bundle(path.parent, recipe, progress=tracker.beat if tracker else None, code=code).items():
                if logical_id not in expected: continue
                require(part['receipt']['source'] == expected[logical_id], 'bundle source population differs')
                require(logical_id not in result, 'source appears in multiple projection bundles')
                result[logical_id] = part
                verified_rows += part['receipt']['rows']
            if tracker:
                tracker.write('verify-projections', len(result), cached_photon_rows=verified_rows)
        for logical_id, item in expected.items():
            if logical_id not in result:
                try:
                    result[logical_id] = Inputs.load_projection(item, recipe, root)
                except FileNotFoundError as error:
                    raise ValueError('projection missing for ' + logical_id +
                                     '; project the exact source catalogue before preparation') from error
                verified_rows += result[logical_id]['receipt']['rows']
                if tracker:
                    if len(result) - tracker.processed >= 100:
                        tracker.write('verify-projections', len(result), cached_photon_rows=verified_rows)
                    else:
                        tracker.beat()
        if tracker: tracker.write('complete', len(result), cached_photon_rows=verified_rows, status='PASS')
        return result
    except Exception as error:
        if tracker: tracker.fail(error)
        raise


def project_shard(inputs, recipe, root, shard, *, scratch_root=None, progress_identity=None):
    """Project one whole-source shard; seal four persistent files on shared storage.

    Completed bundles are verified deterministic skips. Failed shared outputs
    and source staging are retained; a fresh cache root is needed to retry.
    """
    require(all(type(shard[k]) is int for k in ('index', 'begin', 'end', 'events')) and
            shard['index'] >= 0 and 0 <= shard['begin'] < shard['end'] <= len(inputs), 'invalid bundle shard interval')
    selected = inputs[shard['begin']:shard['end']]
    require(selected == sorted(selected, key=lambda item: item['logical_id']) and
            sum(item['entries'] for item in selected) == shard['events'] <= 100000, 'bundle event bound differs')
    staging = scratch_root or os.environ.get('_CONDOR_SCRATCH_DIR')
    require(staging is not None, 'packed projection requires --scratch-root or _CONDOR_SCRATCH_DIR; system /tmp is not selected automatically')
    staging = Path(staging).absolute()
    require(staging.is_dir() and not staging.is_symlink() and staging.resolve() == staging,
            'projection scratch must be an existing physical directory')
    pin = bundle_identity(selected, recipe)
    folder = Path(root) / '.bundles' / pin[:2] / pin
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.bundle.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        complete = (folder / 'bundle.json').exists()
        tracker = WorkProgress(folder / 'progress.json', shard['events'], 'source events',
                               'projection-bundle-' + pin, progress_identity,
                               skip_total=shard['events'] if complete else 0)
        tracker.write('bundle-verify' if complete else 'bundle-start', eligible_photons=0)
        if complete:
            parts = load_bundle(folder, recipe, progress=tracker.beat)
            tracker.write('complete', tracker.total, skip_processed=tracker.total, status='PASS',
                          eligible_photons=sum(p['receipt']['rows'] for p in parts.values()))
            return dict(status='REUSED', shard=shard, receipt=str(folder / 'bundle.json'),
                        progress=str(folder / 'progress.json'))
        require(not (folder / 'rows.bin').exists(), 'incomplete bundle: retain evidence and use a fresh cache root')
        scratch = Path(tempfile.mkdtemp(prefix='photonid-', dir=staging)); scratch.chmod(0o2755)
        tracker.write('bundle-extract', source_staging_directory=str(scratch))
        members = []; offset = 0; rows = 0; aggregate = hashlib.sha256()
        try:
            with (folder / 'rows.bin').open('xb') as output:
                for item in selected:
                    legacy = Inputs.projection_directory(item, recipe, root)
                    if legacy.exists():
                        part = Inputs.load_projection(item, recipe, root)
                        tracker.write('bundle-copy', tracker.processed + item['entries'],
                                      source_logical_id=item['logical_id'])
                    else:
                        def observed(events, photons):
                            tracker.write('bundle-extract', tracker.processed + events,
                                          source_logical_id=item['logical_id'],
                                          eligible_photons=tracker.extra['eligible_photons'] + photons)
                        part = Inputs.project_source(item, recipe, scratch,
                            progress_identity=progress_identity, event_progress=observed)
                    member = part['receipt']; count = member['rows']; length = count * Inputs.projection_dtype(recipe).itemsize
                    checksum = hashlib.sha256(); copied = 0
                    with (part['folder'] / 'rows.bin').open('rb') as source:
                        for block in iter(lambda: source.read(4 * 1024**2), b''):
                            output.write(block); aggregate.update(block); checksum.update(block); copied += len(block)
                            tracker.beat()
                    require(copied == length and checksum.hexdigest() == member['rows_sha256'] and
                            Inputs.projection_witness(part['folder']) == part['witness'], 'source changed while bundling')
                    members.append(dict(projection=member, byte_offset=offset, byte_length=length))
                    offset += length; rows += count
                output.flush(); os.fsync(output.fileno())
            require(tracker.processed == tracker.total, 'bundle source coverage incomplete')
            body = dict(schema=SCHEMA, status='PASS', cache_identity=pin,
                        projection_code=projection_code(),
                        definition=Inputs.population_definition(recipe), record_dtype=Inputs.projection_dtype(recipe).descr,
                        source_events=shard['events'], rows=rows, bytes=offset,
                        parts=members, rows_sha256=aggregate.hexdigest())
            save_json(folder / 'bundle.json', dict(body, bundle_sha256=digest(body)))
            load_bundle(folder, recipe, progress=tracker.beat)
            tracker.write('complete', tracker.total, status='PASS', eligible_photons=rows)
        except BaseException as error:
            tracker.fail(error)
            raise
        shutil.rmtree(scratch)
        return dict(status='PROJECTED', shard=shard, receipt=str(folder / 'bundle.json'),
                    progress=str(folder / 'progress.json'))
