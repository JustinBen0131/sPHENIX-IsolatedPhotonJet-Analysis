"""Versioned reserved-branch augmentation with an exact event/jet join.

The first recipe uses the native calorimeter constituent energy divided by
the original raw jet energy. Negative SUB1 contributions are retained. Zero
or non-finite denominator/numerator produces NaN, never a fabricated zero.
"""
import copy
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from .integrity import equal_arrays
from .integrity import metadata
from .integrity import sha256
from .schema import Schema
from .stream import Writer
from .tables import dumps, require

RADII=('02','03','04')
FRACTIONS=('emcal_frac','ihcal_frac','ohcal_frac')
ALLOWED=tuple('jet'+r+'_'+f for r in RADII for f in FRACTIONS)
# Native raw key is the saved collection ordinal (converter inverse proof).
# eta/phi are calibrated values; comparing them to pre-JES float coordinates
# would incorrectly reject roundoff from the unchanged calibration operation.
FINGERPRINT=('rawpt','area')


def fractions(energy, components):
    require(len(components)==3,'three calorimeter components required')
    energy=np.float32(energy)
    components=np.asarray(components,dtype=np.float32)
    if not np.isfinite(energy) or energy==0:
        return [float('nan')]*3
    # Native Jet accessors, constituent accumulation and division are float32.
    # Promote only the finished values to the collaborator vector<double>.
    return [float(np.float32(v/energy)) if np.isfinite(v) else float('nan') for v in components]


def apply(input_path, output, replay_rows, *, recipe, qualification, input_version, step_size='16 MB', progress=None):
    """Replay records are strictly ordered by (source_entry, radius, jet_index).

    Each contains run, physical_event, radius, jet_index, the exact saved jet
    fingerprint and layer/jet energies. No nearest-jet matching or silent skips.
    The qualification binds the replay code/calibration/jet recipe digest.
    """
    require(set(recipe['allowed_branches'])==set(ALLOWED),'fraction recipe branch allowlist differs')
    require(qualification.get('status')=='PASS' and qualification.get('recipe_sha256')==recipe['sha256'],
            'DST replay is unqualified for this exact recipe')
    require(qualification.get('input_version')==input_version,'replay qualification input version differs')
    before=sha256(input_path); schema=Schema(); iterator=iter(replay_rows); pending=next(iterator,None)
    with uproot.open(input_path) as source:
        meta=metadata(source)
        require('consolidation' in meta,'augmentation requires a stage-2 logical source')
        if meta['source'].get('si_di_role')=='DI':
            require(qualification.get('di_deterministic_replay') is True,'DI replay has no deterministic qualification')
        original_meta=copy.deepcopy(meta)
        writer=Writer(output,schema); jets=0; entries=0
        for a in source['eventTree'].iterate(list(schema.branches),step_size=step_size,library='ak'):
            columns={name:[] for name in ALLOWED}
            for event in range(len(a)):
                event_entry=int(a['evt_src_entry'][event]); run=int(a['evt_run'][event])
                valid=int(a['evt_phys_seq_valid'][event]); physical=int(a['evt_phys_seq'][event])
                for radius in RADII:
                    values=[[] for _ in FRACTIONS]
                    n=len(a['jet'+radius+'_pt'][event])
                    for index in range(n):
                        require(valid==1,'physical event identity unavailable for replay join')
                        require(pending is not None,'replay missing a jet')
                        key=(int(pending['source_entry']),pending['radius'],int(pending['jet_index']))
                        require(key==(event_entry,radius,index),'replay jet order/key mismatch')
                        require(pending['run']==run and pending['physical_event']==physical,'replay event identity mismatch')
                        for field in FINGERPRINT:
                            expected=float(a['jet'+radius+'_'+field][event][index])
                            actual=float(pending[field])
                            require(np.float64(expected).tobytes()==np.float64(actual).tobytes(),
                                    'replay jet fingerprint differs: '+field)
                        result=fractions(float(pending['jet_energy']),pending['calo_energy'])
                        for out,value in zip(values,result): out.append(value)
                        jets+=1; pending=next(iterator,None)
                    for f,v in zip(FRACTIONS,values): columns['jet'+radius+'_'+f].append(v)
            for name,rows in columns.items():
                counts=np.array([len(v) for v in rows],dtype=np.int64)
                flat=np.array([v for row in rows for v in row],dtype=np.float64)
                a=ak.with_field(a,ak.unflatten(flat,counts),name)
            writer.append(a);entries+=len(a)
            if progress: progress(entries,source['eventTree'].num_entries)
        if not entries: writer.append(source['eventTree'].arrays(list(schema.branches),library='ak'))
        require(pending is None,'replay has extra or duplicate jets')
        writer.close()
        for m in [meta]+[c['original_meta'] for c in meta['consolidation']['chunks']]:
            m['schema']['absent_branches']=[b for b in m['schema']['absent_branches'] if b not in ALLOWED]
            m['schema']['filled_reserved']=sorted(set(m['schema']['filled_reserved'])|set(ALLOWED))
        meta['schema']['counters']=writer.counters
        meta['augmentation']=dict(schema='RJAugmentationV1',recipe=recipe,input_version=input_version,
            input_sha256=before,qualification=qualification,base_metadata=original_meta,
            undefined='NaN for zero/nonfinite denominator or nonfinite component; negative fractions retained')
        with uproot.update(output) as out:
            # Keep chunk objects and all other metadata; updated direct metadata is written once below.
            out.copy_from(source,filter_name=lambda n: not(n.startswith('meta/') and '/' not in n[5:]),
                filter_classname=lambda cls:cls not in ('TTree','TDirectory','TDirectoryFile'),require_matches=False)
            for key,value in meta.items(): out['meta/'+key]=dumps(value)
    require(sha256(input_path)==before,'stage-2 input changed')
    proof=prove(input_path,output,step_size=step_size)
    return dict(schema='RJAugmentationReceiptV1',status='PASS',recipe_sha256=recipe['sha256'],
        logical_source_id=meta['consolidation']['logical_source_id'],input_sha256=before,input_version=input_version,
        output=str(output),output_sha256=sha256(output),entries=entries,jets=jets,proof=proof)


def prove(original, augmented, *, step_size='16 MB'):
    from .check import check
    schema=Schema(); offset=0
    with uproot.open(original) as a,uproot.open(augmented) as b:
        before_meta=metadata(a);after_meta=metadata(b)
        require(after_meta['augmentation']['base_metadata']==before_meta,'augmentation lost original metadata')
        expected=copy.deepcopy(before_meta)
        for m in [expected]+[c['original_meta'] for c in expected['consolidation']['chunks']]:
            m['schema']['absent_branches']=[n for n in m['schema']['absent_branches'] if n not in ALLOWED]
            m['schema']['filled_reserved']=sorted(set(m['schema']['filled_reserved'])|set(ALLOWED))
        expected['schema']['counters']=after_meta['schema']['counters']
        require({k:v for k,v in after_meta.items() if k!='augmentation'}==expected,'augmentation changed protected metadata')
        for name,cls in a.classnames(recursive=True,cycle=False).items():
            if cls in ('TTree','TDirectory','TDirectoryFile') or (name.startswith('meta/') and '/' not in name[5:]):continue
            old,_=a.key(name).get_uncompressed_chunk_cursor();new,_=b.key(name).get_uncompressed_chunk_cursor()
            require(np.array_equal(old.raw_data,new.raw_data),'augmentation changed carried object: '+name)
        require(a['eventTree'].num_entries==b['eventTree'].num_entries,'augmentation changed event count')
        for before in a['eventTree'].iterate(list(schema.branches),step_size=step_size,library='ak'):
            after=b['eventTree'].arrays(list(schema.branches),entry_start=offset,entry_stop=offset+len(before),library='ak')
            for name in schema.branches:
                if name not in ALLOWED: require(equal_arrays(before[name],after[name]),'protected branch changed: '+name)
            for radius in RADII:
                counts=ak.num(after['jet'+radius+'_pt'])
                for field in FRACTIONS:
                    require(bool(ak.all(ak.num(after['jet'+radius+'_'+field])==counts)),'fraction radius alignment differs')
            offset+=len(before)
    structural=check(augmented)
    require(structural['ok'],'augmented structure failed: '+str(structural['problems'][:3]))
    return dict(status='PASS',entries=offset,protected_branches=len(schema.branches)-len(ALLOWED),allowed_branches=list(ALLOWED))
