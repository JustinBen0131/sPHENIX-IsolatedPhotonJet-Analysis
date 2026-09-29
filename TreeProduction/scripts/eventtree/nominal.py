"""Recompute nominal scores/classes from the saved tree and pinned packages."""
from types import SimpleNamespace

import awkward as ak
import numpy as np
import uproot

from .integrity import metadata, sha256, equal_arrays, copy_objects, prove_objects
from pathlib import Path
from .schema import Schema
from .tables import require


def check_arrays(arrays, meta, photon_freeze, isolation_freeze=None):
    from .photon_id import augment as id_augment
    from .isolation_class import augment as iso_augment, augment_pp
    schema=Schema(); values={}
    for name,b in schema.branches.items():
        a=arrays[name]
        if b.kind=='vector':
            counts=ak.to_numpy(ak.num(a,axis=1)).astype(np.int64)
            values[name]=(b.kind,ak.to_numpy(ak.flatten(a)),counts)
        else: values[name]=(b.kind,ak.to_numpy(a),None)
    ret=np.flatnonzero(ak.to_numpy(arrays['evt_record'])==0)
    verified=[]
    def vector(name, data, counts):
        _,saved,all_counts=values[name]
        require(np.array_equal(all_counts[ret],counts),'nominal collection alignment differs: '+name)
        require(saved.dtype==np.asarray(data).dtype,'nominal branch dtype differs: '+name)
        require(np.array_equal(saved,data,equal_nan=True),'saved nominal values differ: '+name)
        verified.append(name)
    cm=dict(meta,strides=meta['schema']['strides'])
    proxy=SimpleNamespace(out=SimpleNamespace(values=values,vector=vector),meta=cm,
                          ret_pos=ret,system=meta['schema']['system'])
    rid=id_augment(proxy,photon_freeze)
    require(meta['photon_id']['freeze_sha256']==rid['freeze_sha256'],'photon freeze differs')
    if proxy.system=='auau':
        require(isolation_freeze is not None,'AuAu isolation freeze required')
        iso=iso_augment(proxy,isolation_freeze)
    else: iso=augment_pp(proxy)
    require(meta['isolation']['freeze_sha256']==iso['freeze_sha256'],'isolation freeze differs')
    return dict(status='PASS',entries=len(arrays),photons=int(ak.sum(ak.num(arrays['pho_et']))),branches=verified)


def check(path, photon_freeze, isolation_freeze=None, step_size='16 MB'):
    entries=photons=0
    with uproot.open(path) as f:
        meta=metadata(f);t=f['eventTree']
        chunks=meta.get('consolidation',{}).get('chunks') or [dict(entry_begin=0,entry_end=t.num_entries,original_meta=meta)]
        for c in chunks:
            if c['entry_begin']==c['entry_end']:
                check_arrays(t.arrays(entry_start=0,entry_stop=0,library='ak'),c['original_meta'],photon_freeze,isolation_freeze)
            for arrays in t.iterate(entry_start=c['entry_begin'],entry_stop=c['entry_end'],step_size=step_size,library='ak'):
                r=check_arrays(arrays,c['original_meta'],photon_freeze,isolation_freeze)
                entries+=r['entries'];photons+=r['photons']
        require(entries==t.num_entries,'nominal regression coverage differs')
    return dict(status='PASS',entries=entries,photons=photons)


def finalize_arrays(arrays, meta, isolation_freeze=None):
    """Invoke the existing shared definitions once; check the resulting values.

    A complete input mapping is mandatory. Finalized reference trees are read
    by regression, not passed through this new-production fill a second time.
    """
    import awkward as ak
    import numpy as np
    from types import SimpleNamespace
    from eventtree.stream import Output
    from eventtree.schema import Schema
    from eventtree.photon_id import augment, FREEZE
    from eventtree.isolation_class import augment as augment_iso, augment_pp
    from eventtree.nominal import check_arrays
    nominal_names = ('pho_bdt_score', 'pho_bdt_state', 'pho_id_class', 'pho_iso_class')
    missing = set(meta['native_capture']['missing_branches']) - set(nominal_names)
    require(not missing, 'nominal finalization requires complete capture mapping: '+', '.join(sorted(missing)[:8]))
    ret = np.flatnonzero(ak.to_numpy(arrays.evt_record) == 0)
    output = Output(Schema(), ret, len(arrays))
    for name, branch in output.schema.branches.items():
        value = arrays[name]
        output.values[name] = (branch.kind, ak.to_numpy(ak.flatten(value)) if branch.kind == 'vector' else ak.to_numpy(value),
                              ak.to_numpy(ak.num(value)).astype(np.int64) if branch.kind == 'vector' else None)
    contract = dict(meta, strides=meta['schema']['strides'])
    proxy = SimpleNamespace(out=output, meta=contract, ret_pos=ret, system=meta['schema']['system'])
    meta['photon_id'] = augment(proxy, FREEZE)
    if proxy.system == 'auau':
        require(isolation_freeze is not None, 'AuAu requires the pinned isolation freeze')
        meta['isolation'] = augment_iso(proxy, isolation_freeze)
    else:
        meta['isolation'] = augment_pp(proxy)
    for name in nominal_names:
        _, values, counts = output.values[name]
        arrays = ak.with_field(arrays, ak.unflatten(values, counts), name)
    meta['schema']['filled_reserved'] = sorted(set(meta['schema']['filled_reserved']) | set(nominal_names))
    meta['schema']['absent_branches'] = sorted(set(meta['schema']['absent_branches']) - set(nominal_names))
    meta['native_capture']['missing_branches'] = []
    check_arrays(arrays, meta, FREEZE, isolation_freeze)
    return arrays, meta


def finalize_file(candidate, destination, isolation_freeze=None, batch_entries=128, progress=None):
    """Bounded shared finalization with saved-file conservation and nominal proof."""
    import awkward as ak
    import copy
    import uproot
    from eventtree.integrity import metadata
    from eventtree.integrity import equal_arrays
    from eventtree.schema import Schema
    from eventtree.stream import Writer
    from eventtree.tables import dumps
    from eventtree.check import check_arrays
    from eventtree.nominal import check as check_nominal
    from eventtree.photon_id import FREEZE
    require(1 <= batch_entries <= 1024, 'finalization batch size outside [1,1024]')
    require(Path(candidate).resolve() != Path(destination).resolve(), "finalization needs a separate candidate")
    require(not Path(destination).exists(), "finalization destination already exists")
    before = sha256(Path(candidate))
    schema = Schema()
    nominal = {'pho_bdt_score', 'pho_bdt_state', 'pho_id_class', 'pho_iso_class'}
    with uproot.open(candidate) as f:
        original = metadata(f)
        require(original["schema"]["schema_sha256"] == schema.sha256, "native schema digest differs")
        require(not (set(original["native_capture"]["missing_branches"]) - nominal), "native mapping is incomplete")
        writer = Writer(destination)
        combined = None
        processed = 0
        try:
            tree = f['eventTree']
            batches = tree.iterate(list(schema.branches), step_size=batch_entries, library='ak') if tree.num_entries else [tree.arrays(list(schema.branches), library='ak')]
            for arrays in batches:
                meta = copy.deepcopy(original)
                meta['schema'].update(entries=len(arrays), retained_events=int(ak.sum(arrays.evt_record == 0)))
                arrays, meta = finalize_arrays(arrays, meta, isolation_freeze)
                require(check_arrays(arrays, meta)['ok'], 'finalization batch has invalid structure')
                writer.append(arrays)
                processed += len(arrays)
                if progress:
                    progress('nominal_fill', processed, tree.num_entries)
                if combined is None:
                    combined = meta
                else:
                    for key in ('photon_id', 'isolation'):
                        current, total = meta[key], combined[key]
                        require(current['freeze_sha256'] == total['freeze_sha256'], 'nominal package changed within source')
                        total['photons'] += current['photons']
                        for tally in ('states', 'classes'):
                            if tally in current:
                                for label, number in current[tally].items():
                                    total[tally][label] = total[tally].get(label, 0) + number
        finally:
            writer.close()
        combined['schema'].update(entries=tree.num_entries,
            retained_events=original['schema']['retained_events'], counters=writer.counters)
        combined['native_capture']['status'] = 'FINALIZED_SAVED_CHECKS_PENDING'
        with uproot.update(destination) as out:
            copy_objects(f, out)
            for name, value in combined.items():
                out['meta/'+name] = dumps(value)
        checked = 0
        with uproot.open(destination) as final:
            require(metadata(final) == combined, "saved final metadata differs")
            prove_objects(f, final)
            require(final['eventTree'].num_entries == tree.num_entries, 'saved final entry count differs')
            require([name for name, cls in final.classnames(cycle=False).items() if cls == 'TTree'] == ['eventTree'],
                    'final output does not contain exactly one eventTree')
            for begin in range(0, tree.num_entries, batch_entries):
                stop = min(begin+batch_entries, tree.num_entries)
                a = tree.arrays(list(schema.branches), entry_start=begin, entry_stop=stop, library='ak')
                b = final['eventTree'].arrays(list(schema.branches), entry_start=begin, entry_stop=stop, library='ak')
                for name in set(schema.branches)-nominal:
                    require(equal_arrays(a[name], b[name]), 'saved final changed captured payload: '+name)
                cm = copy.deepcopy(combined)
                cm['schema'].update(entries=len(b), retained_events=int(ak.sum(b.evt_record == 0)))
                require(check_arrays(b, cm)['ok'], 'saved final eventTree structure failed')
                checked += len(b)
                if progress:
                    progress('saved_conservation', checked, tree.num_entries)
        require(checked == tree.num_entries, 'saved content conservation has incomplete coverage')
    require(sha256(Path(candidate)) == before, 'candidate changed during finalization')
    nominal_proof = check_nominal(destination, FREEZE, isolation_freeze, step_size=batch_entries)
    return dict(status='PASS_SAVED_FINALIZATION', input_sha256=before, output_sha256=sha256(Path(destination)),
                entries=checked, conserved_branches=len(schema.branches)-len(nominal), nominal=nominal_proof,
                replacement_certified=False)

