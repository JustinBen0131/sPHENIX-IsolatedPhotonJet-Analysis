"""Metadata-aware, bounded-memory consolidation of one original logical source.

Never sums cumulative scaler snapshots or normalization metadata. Each original
chunk's metadata and serialized objects remain independently recoverable. Coded
strings have a union dictionary with reversible per-chunk index maps. All other
payloads, event order and event-local object references are bit-preserved.
"""
import copy
import hashlib
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from .integrity import equal_arrays
from .schema import Schema
from .check import INDEX_SUFFIXES
from .stream import Writer
from .tables import dumps, require
from .integrity import sha256


from .integrity import metadata


def definition(meta, schema=None):
    """Physics identity, independent of each chunk's statistics and ROOT version."""
    schema = schema or Schema()
    require(schema.accepts_digest(meta['schema']['schema_sha256']), 'source schema differs')
    p, i = meta.get('photon_id',{}), meta.get('isolation',{})
    return (schema.contract_sha256,p.get('freeze_sha256'),p.get('model_sha256'),i.get('freeze_sha256'))


def merge_shapes(metas):
    """Union observed shape constants without treating empty data as a recipe.

    Per-chunk metadata remains the inverse-proof authority. Zero strides carry
    no observation; a missing stride is allowed only for an absent branch.
    Collection identity is always required; optional constants are compared
    wherever they were observed. Conflicting nonempty observations refuse.
    """
    strides = {}
    collections = copy.deepcopy(metas[0]['collections'])
    for m in metas:
        for name, value in m['schema']['strides'].items():
            require(isinstance(value, int) and value >= 0, 'invalid stride: '+name)
            old = strides.get(name, 0)
            require(not old or not value or old == value, 'conflicting nonzero stride: '+name)
            strides[name] = old or value
        require(set(m['collections']) == set(collections), 'jet collection families differ')
        for family, rows in m['collections'].items():
            expected = collections[family]
            require(len(rows) == len(expected), 'jet collection roster differs')
            for a, b in zip(expected, rows):
                for key in ('name', 'input_identity', 'radius', 'radius_code'):
                    require(a.get(key) == b.get(key), 'jet collection identity differs: '+key)
                for key, value in b.items():
                    require(key not in a or a[key] == value, 'jet collection recipe differs: '+key)
                    a[key] = copy.deepcopy(value)
    for m in metas:
        for name in strides:
            require(name in m['schema']['strides'] or name in m['schema']['absent_branches'],
                    'non-absent branch lacks stride: '+name)
    return strides, collections


def remap(arrays, maps, reverse=False):
    """Map chunk-local codes to the logical-source union (or back). Scalar and
    vector coded fields are both remapped; an index field's -1 ("no entry")
    is kept as is."""
    for name, mapping in maps.items():
        vector = arrays[name].ndim > 1
        counts = ak.to_numpy(ak.num(arrays[name],axis=1)) if vector else None
        values = ak.to_numpy(ak.flatten(arrays[name]) if vector else arrays[name])
        lowest = -1 if name.endswith(INDEX_SUFFIXES) else 0
        present = values >= 0
        if reverse:
            invert = {int(v):k for k,v in enumerate(mapping)}
            require(all(int(v) in invert for v in np.unique(values[present])), 'unmapped consolidated code: '+name)
            require(not len(values) or values.min() >= lowest, 'invalid consolidated code: '+name)
            values = np.asarray([invert[int(v)] if v >= 0 else int(v) for v in values],dtype=values.dtype)
        else:
            require(not len(values) or (values.min() >= lowest and values.max() < len(mapping)), 'invalid source code: '+name)
            mapped = np.asarray(mapping,dtype=values.dtype)
            values = np.where(present, mapped[np.where(present, values, 0)] if len(mapped) else values, values).astype(values.dtype)
        arrays = ak.with_field(arrays,ak.unflatten(values,counts) if vector else values,name)
    return arrays


def expected_ranges(chunks, expected_end, exclusions=()):
    end = 0
    for c in sorted([*chunks,*exclusions],key=lambda c:(c['event_begin'],c['event_end'])):
        require(c['event_begin'] == end and c['event_end'] >= end,'source range gap, overlap or ordering error')
        end = c['event_end']
    require(end == expected_end,'source range does not reach expected end')


def consolidate(chunks, output, *, logical_source_id, expected_events, exclusions=(), step_size='16 MB', progress=None):
    """Each chunk binds path, sha256, original unit/piece and [event_begin,event_end).

    The caller must supply the frozen logical-source identity and exact expected
    interval. Hashes are rechecked before and after reading. All failures leave
    an unpublished candidate for diagnosis; existing files are never replaced.
    """
    require(chunks,'no source chunks')
    require(len({c['sha256'] for c in chunks}) == len(chunks),'same payload assigned to multiple chunks')
    require(len({c['piece_id'] for c in chunks}) == len(chunks),'duplicate piece identity')
    require(all(e.get('authority_json') for e in exclusions),'excluded ranges lack pinned authority')
    expected_ranges(chunks,expected_events,exclusions)
    expected_entries=sum(c['event_end']-c['event_begin'] for c in chunks)
    schema = Schema(); metas=[]; union={}; entries=0; retained=0; centrality_hist=[]
    for c in chunks:
        require(sha256(c['path']) == c['sha256'],'input changed: '+c['path'])
        with uproot.open(c['path']) as f:
            m=metadata(f)
            centrality_hist.append(('hist/SIM/h_centrality' in f,f['eventTree'].num_entries))
            require('consolidation' not in m,'source is already consolidated')
            require(schema.accepts_digest(m['schema']['schema_sha256']),'source schema differs')
            require(f['eventTree'].num_entries == m['schema']['entries'],'source entry metadata differs')
            require(m['schema']['entries'] == c['event_end']-c['event_begin'],
                    'source event accounting incomplete: retained plus rejected must cover interval')
            if metas:
                require(definition(m,schema) == definition(metas[0],schema),'nominal definition differs')
                for k in ('system','lane'):
                    require(m['schema'][k] == metas[0]['schema'][k], 'chunk schema property differs: '+k)
                for k in ('dataset','run','segment','sample','period','si_di_role'):
                    require(m['source'].get(k) == metas[0]['source'].get(k),'chunk source differs: '+k)
            for branch, labels in m.get('codes',{}).items():
                target=union.setdefault(branch,[])
                for label in labels:
                    if label not in target: target.append(label)
            metas.append(m); entries+=f['eventTree'].num_entries; retained+=m['schema']['retained_events']
    require(entries == expected_entries,'event count differs from admitted source coverage')
    if any(present for present,count in centrality_hist):
        require(all(present or count==0 for present,count in centrality_hist),
                'nonempty chunk is missing the source centrality histogram')
    combined=copy.deepcopy(metas[0]); combined['codes']=union
    combined['schema']['strides'], combined['collections'] = merge_shapes(metas)
    combined['schema'].update(schema_sha256=schema.sha256,entries=entries,retained_events=retained,
        absent_branches=sorted(set.intersection(*(set(m['schema']['absent_branches']) for m in metas))))
    trigger=[m.get('trigger_run_info') for m in metas]
    if not all(v==trigger[0] for v in trigger):combined.pop('trigger_run_info',None)
    if trigger[0] and 'name' in trigger[0]:
        selection={k:trigger[0][k] for k in ('name','initial_prescale')}
        require(all(v and {k:v[k] for k in selection}==selection for v in trigger),
                'trigger names/prescales vary within logical source')
        combined['trigger_selection']=selection
    combined['source']['logical_source_id']=logical_source_id
    combined['source']['consolidated']=True
    combined['consolidation']=dict(schema='RJLogicalSourceV1',logical_source_id=logical_source_id,
        expected_events=expected_entries,original_source_events=expected_events,exclusions=list(exclusions),chunks=[],
        metadata_policy='per-chunk originals; cumulative counters and normalization are not added',
        content_proof='payload bitwise after reversible dictionary mapping')
    writer=Writer(output,schema); offset=0; physical_seen=set()
    for index,(c,m) in enumerate(zip(chunks,metas)):
        maps={b:[union[b].index(v) for v in labels] for b,labels in m.get('codes',{}).items()}
        chunk=dict(c,entry_begin=offset,entry_end=offset+m['schema']['entries'],
                   original_meta=m, code_maps=maps,object_prefix=f'meta/chunks/{index:06d}/objects/')
        combined['consolidation']['chunks'].append(chunk)
        with uproot.open(c['path']) as f:
            t=f['eventTree']; local=0
            # Iterate still needs one empty typed batch to create all branches.
            batches=t.iterate(list(schema.branches),step_size=step_size,library='ak') if t.num_entries else [t.arrays(list(schema.branches),library='ak')]
            for a in batches:
                seq=ak.to_numpy(a['evt_src_entry'])
                require(np.array_equal(seq,np.arange(c['event_begin']+local,c['event_begin']+local+len(a))),
                        'missing, duplicated or wrongly associated source events')
                for run,physical,valid in zip(ak.to_numpy(a['evt_run']),ak.to_numpy(a['evt_phys_seq']),ak.to_numpy(a['evt_phys_seq_valid'])):
                    if valid == 1:
                        event_key=(int(run),int(physical))
                        require(event_key not in physical_seen,'duplicate physical event within logical source')
                        physical_seen.add(event_key)
                writer.append(remap(a,maps)); local+=len(a)
                if progress: progress(offset+local,entries)
        offset=chunk['entry_end']
        require(sha256(c['path']) == c['sha256'],'input changed during consolidation')
    writer.close(); combined['schema']['counters']=writer.counters
    combined['schema']['encoding']='counted'
    with uproot.update(output) as out:
        for key,value in combined.items(): out['meta/'+key]=dumps(value)
        for c in combined['consolidation']['chunks']:
            with uproot.open(c['path']) as f:
                prefix=c['object_prefix']
                out.copy_from(f,filter_classname=lambda cls:cls not in ('TTree','TDirectory','TDirectoryFile'),
                              rename=lambda name:prefix+name,require_matches=False)
    # This is the maintained SIM centrality reader's additive event histogram.
    # Original histograms remain byte-preserved below their chunk objects.
    from .rootio import root
    R=root(); destination=R.TFile.Open(str(output),'UPDATE'); aggregate=None
    key='hist/SIM/h_centrality'
    for c in combined['consolidation']['chunks']:
        h=destination.Get(c['object_prefix']+key)
        if not h:continue
        require(h.InheritsFrom('TH1'),'centrality object is not a histogram')
        if aggregate is None:aggregate=h.Clone('h_centrality');aggregate.SetDirectory(0)
        else:
            require(aggregate.GetNbinsX()==h.GetNbinsX() and all(aggregate.GetXaxis().GetBinLowEdge(i)==h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)),
                    'centrality histogram edges differ')
            require(bool(aggregate.Add(h)),'centrality histogram merge failed')
    if aggregate is not None:
        d=destination.mkdir('hist');d=d.mkdir('SIM');d.cd();aggregate.Write('h_centrality')
    destination.Close()
    proof=prove(output,step_size=step_size)
    require(proof['status']=='PASS','consolidation readback proof failed')
    return dict(schema='RJLogicalSourceReceiptV1',status='PASS',logical_source_id=logical_source_id,
                output=str(output),output_sha256=sha256(output),entries=entries,pieces=len(chunks),proof=proof)


def prove(path, *, step_size='16 MB'):
    """Read every saved branch and every serialized object, against each input."""
    schema=Schema(); checked=0; objects=0
    with uproot.open(path) as final:
        m=metadata(final); chunks=m['consolidation']['chunks']; t=final['eventTree']
        require(t.num_entries == m['consolidation']['expected_events'],'consolidated entry count differs')
        for c in chunks:
            require(sha256(c['path']) == c['sha256'],'source hash drift during proof')
            with uproot.open(c['path']) as f:
                require(metadata(f)==c['original_meta'],'chunk provenance changed')
                start=0
                for a in f['eventTree'].iterate(list(schema.branches),step_size=step_size,library='ak'):
                    b=t.arrays(list(schema.branches),entry_start=c['entry_begin']+start,
                               entry_stop=c['entry_begin']+start+len(a),library='ak')
                    b=remap(b,c['code_maps'],reverse=True)
                    for name in schema.branches:
                        require(equal_arrays(a[name],b[name]),'consolidated payload differs: '+name)
                    checked+=len(a); start+=len(a)
                for name,cls in f.classnames(recursive=True,cycle=False).items():
                    if cls in ('TTree','TDirectory','TDirectoryFile'): continue
                    a,_=f.key(name).get_uncompressed_chunk_cursor()
                    b,_=final.key(c['object_prefix']+name).get_uncompressed_chunk_cursor()
                    require(np.array_equal(a.raw_data,b.raw_data),'chunk object differs: '+name); objects+=1
            require(sha256(c['path']) == c['sha256'],'input changed during proof')
        require(checked==t.num_entries,'conservation coverage differs')
    structural=check_consolidated(path,step_size=step_size)
    require(structural['ok'],'consolidated structure: '+str(structural['problems'][:3]))
    return dict(status='PASS',entries=checked,objects=objects,chunks=len(chunks))


def check_consolidated(path, *, step_size='16 MB'):
    from .check import check_arrays, tree_keys
    problems=[]; entries=0
    with uproot.open(path) as f:
        m=metadata(f); t=f['eventTree']
        chunks=m['consolidation']['chunks']; end=0
        expected_ranges(chunks,m['consolidation'].get('original_source_events',m['consolidation']['expected_events']),
                        m['consolidation'].get('exclusions',()))
        strides,collections=merge_shapes([c['original_meta'] for c in chunks])
        require(strides==m['schema']['strides'] and collections==m['collections'],'consolidated shape union differs')
        for c in chunks:
            require(c['entry_begin']==end and c['entry_end']>=end,'invalid consolidated chunk coverage')
            require(c['entry_end']-end==c['original_meta']['schema']['entries'],'chunk accounting metadata differs')
            end=c['entry_end']
        require(end==t.num_entries==m['consolidation']['expected_events'],'consolidated manifest coverage differs')
        for c in chunks:
            for a in t.iterate(entry_start=c['entry_begin'],entry_stop=c['entry_end'],step_size=step_size,library='ak'):
                cm=copy.deepcopy(c['original_meta'])
                cm['codes']=m['codes']; cm['schema']['counters']=m['schema']['counters']
                cm['schema']['entries']=len(a); cm['schema']['retained_events']=int(ak.sum(a['evt_record']==0))
                report=check_arrays(a,cm); problems.extend(report['problems']); entries+=len(a)
        if entries != t.num_entries: problems.append('chunk entry coverage differs')
    if len(tree_keys(path)) != 1: problems.append('not exactly one TTree')
    return dict(ok=not problems,entries=entries,problems=problems)
