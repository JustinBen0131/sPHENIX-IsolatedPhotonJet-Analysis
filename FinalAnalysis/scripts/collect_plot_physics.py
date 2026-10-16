#!/usr/bin/env python3
"""Collect one completed native draw, never a glob of historical figures.

Scientific incompleteness is recorded, not guessed away. A corrupt or changed
receipt/input is an error and does not replace the current generation pointer.
This is metadata I/O only: it performs no histogramming, fitting or rendering.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import tempfile

from plot_physics_inventory import build, canonical, require, sha
from plot_physics_sources import source_inventory, source_roots
from scientific_origins import validate_origin, source_path, locate_source
from native_producer_ancestry import ancestry_errors


def atomic(path, value):
    path=Path(path)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False) as stream:
        stream.write(canonical(value)+'\n');stream.flush();os.fsync(stream.fileno());name=stream.name
    os.replace(name,path)


def native_metadata(path):
    """Read producer-written TObjString values, without evaluating input code."""
    import ROOT
    root=ROOT.TFile.Open(str(path),'READ')
    require(root and not root.IsZombie(),'cannot read native input metadata')
    try:
        values={}
        for key in ('native_parameter_exports','scientific_parameter_bindings','native_source_witness','native_source_ancestry'):
            saved=root.Get(key) or root.Get('meta/'+key)
            if saved:
                require(saved.InheritsFrom('TObjString'),'native metadata is not a JSON string')
                values[key]=json.loads(str(saved.GetString().Data()))
        return values
    finally:root.Close()


def roster(fig):
    rows=[]
    for line in (fig/'FIGURES.txt').read_text().splitlines():
        if not line.startswith('fig/'):continue
        name=line[4:]
        require(name and name not in rows,'empty or duplicate output in current FIGURES.txt')
        path=(fig/name).resolve()
        require(path.is_relative_to(fig) and path!=fig,'current figure path escapes generation')
        require(path.is_file() and Path(str(path)+'.json').is_file(),'current figure or receipt missing: '+name)
        rows.append(name)
    return rows


def validate_exports(exports,source_root):
    """Verify every native export; return valid, non-adoptable runtime sites.

    Generated run YAML is a real consumed input, but its absolute location is
    not authority to edit a canonical template. A changed source/span remains
    an integrity error for either kind of site.
    """
    import math,struct
    def numeric(value):
        if isinstance(value,str) and value.startswith('f64:'):
            require(len(value)==20,'malformed portable native value')
            try:decoded=struct.unpack('>d',bytes.fromhex(value[4:]))[0]
            except (ValueError,struct.error) as error:raise ValueError('malformed portable native value') from error
            require(math.isfinite(decoded),'nonfinite portable native value')
        elif isinstance(value,float):require(math.isfinite(value),'nonfinite native value')
        elif isinstance(value,list):
            for item in value:numeric(item)
        elif isinstance(value,dict):
            for item in value.values():numeric(item)
    external=[]
    for key,export in exports.items():
        require(isinstance(key,str) and key and isinstance(export,dict) and 'value' in export,'malformed native parameter export')
        numeric(export['value'])
        origin={k:export[k] for k in ('path','before_text','selector','canonical_key')}
        origin.update(kind='source_text',encoding='native_value',export_id=key)
        path=Path(export['path'])
        if path.is_absolute():
            path=path.resolve()
            require(path.is_file(),'native runtime origin is missing')
            require(sha(path)==export.get('source_sha256'),'native export source bytes differ')
            locate_source(path.read_text(),origin)
            external.append({'export_id':key,'path':str(path),'source_sha256':export['source_sha256'],
                'state':'observed_external_runtime','canonical_adoption':'UNAVAILABLE'})
        else:
            require(sha(source_path(source_root,export['path']))==export.get('source_sha256'),'native export source bytes differ')
            validate_origin(origin,export['value'],source_root,exports)
    return external


def collect(fig,source_root,registry_path,qualification_path=None,input_roots=(),*,reader=native_metadata):
    fig,source_root=Path(fig).resolve(),Path(source_root).resolve()
    names=roster(fig);roster_hash=sha(fig/'FIGURES.txt');sources=source_inventory(source_root)
    registry_path=Path(registry_path).resolve()
    registry_hash=sha(registry_path) if registry_path.is_file() else None
    registry=json.loads(registry_path.read_text()) if registry_path.is_file() else {'families':{}}
    qualification_path=Path(qualification_path).resolve() if qualification_path else source_root/'FinalAnalysis/plot_physics_qualification.json'
    qualification=json.loads(qualification_path.read_text()) if qualification_path.is_file() else None
    reasons=[];observed={};inputs=[];qualified_inputs=True
    qualified_inventories=[r.get('source_inventory') for f in registry.get('families',{}).values()
        for r in f.get('tested_native_releases',[]) if r.get('qualification_report_sha256')]
    if qualification and qualification.get('status')=='PASS':qualified_inventories.append(qualification.get('source_inventory'))
    if not registry_path.is_file():reasons.append('renderer-family registry absent')
    if not qualification:reasons.append('native qualification absent')
    elif (qualification.get('status')!='PASS' or qualification.get('source_inventory')!=sources or
          qualification.get('source_roots')!=source_roots() or qualification.get('source_root')!=str(source_root)):
        reasons.append('native qualification does not bind this exact producer source inventory');qualification=None
    for item in input_roots:
        path=Path(item).resolve();before=sha(path);metadata={}
        try:metadata=reader(path)
        except ImportError:reasons.append(str(path)+': ROOT metadata reader unavailable')
        require(sha(path)==before,'input changed during native metadata readback')
        exports=metadata.get('native_parameter_exports',{})
        require(isinstance(exports,dict),'native parameter exports must be an object')
        for key,value in exports.items():
            require(key not in observed or observed[key]==value,'contradictory native parameter exports')
            observed[key]=value
        inputs.append({'path':str(path),'sha256':before,'metadata':metadata})
        witness=metadata.get('native_source_witness',{})
        ancestry_problems=ancestry_errors(witness,metadata.get('native_source_ancestry'),qualified_inventories)
        if ancestry_problems:
            qualified_inputs=False
            reasons.extend(str(path)+': '+problem for problem in ancestry_problems)
        if not (qualification and witness.get('schema')=='NativeProducerSourceWitnessV1' and
                witness.get('state')=='recorded' and witness.get('source_roots')==source_roots() and
                witness.get('source_root')==str(source_root) and
                witness.get('source_inventory')==sources):
            qualified_inputs=False
            reasons.append(str(path)+': producing binary has no matching qualified compiled source witness')
    # The normal producing binary supplies its actual evaluated values and a
    # compile-time source witness. Matching the release fixture qualification
    # authorizes this deterministic seal, never values copied from fixtures.
    generation_qualification=None;external_exports=[]
    if qualification and qualified_inputs and inputs:
        external_exports=validate_exports(observed,source_root)
        for row in external_exports:
            reasons.append('observed runtime origin is outside canonical source adoption: '+row['path']+' ['+row['export_id']+']')
    if qualification and qualified_inputs and inputs:
        generation_qualification=dict(qualification)
        generation_qualification.update({'native_value_exports':observed,
            'release_qualification':{'path':str(qualification_path),'sha256':sha(qualification_path)},
            'native_inputs':[{'path':row['path'],'sha256':row['sha256']} for row in inputs],
            'qualification_kind':'exact_native_generation'})
    elif observed:reasons.append('observed native exports lack qualified producing-binary readback')
    if not inputs:reasons.append('exact native producer input readback is required for owner admission')
    receipts=[Path(str(fig/name)+'.json') for name in names]
    for path in receipts:
        receipt=json.loads(path.read_text());recorded=receipt.get('inputs',{})
        for item in inputs:require(item['sha256'] in recorded.values(),'producing receipt does not bind the exact native input: '+item['path'])
        parents=receipt.get('native_source_ancestry',{}).get('inputs',[])
        for item in inputs:
            matching=[row for row in parents if isinstance(row,dict) and row.get('source_path') and
                Path(row['source_path']).resolve()==Path(item['path'])]
            expected=item['metadata']
            if not any(row.get('witness')==expected.get('native_source_witness') and
                       row.get('ancestry')==expected.get('native_source_ancestry') for row in matching):
                reasons.append(str(path)+': renderer ancestry does not retain actual input metadata: '+item['path'])
        witness=receipt.get('native_source_witness',{})
        reasons.extend(str(path)+': '+problem for problem in ancestry_errors(
            witness,receipt.get('native_source_ancestry'),qualified_inventories))
        if not (qualification and witness.get('schema')=='NativeProducerSourceWitnessV1' and
                witness.get('state')=='recorded' and witness.get('source_roots')==source_roots() and
                witness.get('source_root')==str(source_root) and
                witness.get('source_inventory')==sources):
            reasons.append(str(path)+': renderer binary has no matching qualified compiled source witness')
    inventory=build(receipts,fig,registry,source_root,generation_qualification)
    inventory.update({'expected_outputs':names,'source_inventory':sources,'source_root':str(source_root),'source_roots':source_roots(),
        'input_metadata':inputs,'observed_native_value_exports':observed,'runtime_origin_exports':external_exports,
        'native_value_exports':observed if generation_qualification else {},
        'generation_roster':{'path':'FIGURES.txt','sha256':roster_hash},
        'qualification':{'path':str(qualification_path),'sha256':sha(qualification_path)} if qualification else None})
    inventory['unavailable_reasons']=sorted(set(inventory['unavailable_reasons']+reasons))
    inventory['owner_admissible']=inventory['owner_admissible'] and not reasons
    require(not registry_hash or sha(registry_path)==registry_hash,'renderer registry changed during collection')
    require(source_inventory(source_root)==sources,'producer source changed while collecting generation')
    require(sha(fig/'FIGURES.txt')==roster_hash and roster(fig)==names,'figure roster changed while collecting generation')
    for item in inputs:require(sha(item['path'])==item['sha256'],'native input changed during inventory collection')
    if generation_qualification:
        require(validate_exports(observed,source_root)==external_exports,'runtime origin changed during collection')
    for name,row in inventory['outputs'].items():
        require(sha(row['native_path'])==row['native_sha256'] and sha(row['receipt_path'])==row['receipt_sha256'],
                'current figure or receipt changed during inventory collection: '+name)
    data=(canonical(inventory)+'\n').encode();import hashlib
    digest=hashlib.sha256(data).hexdigest();directory=fig/'.physics';directory.mkdir(exist_ok=True)
    target=directory/(digest+'.json')
    try:
        with target.open('xb') as stream:stream.write(data);stream.flush();os.fsync(stream.fileno())
    except FileExistsError:require(target.read_bytes()==data,'immutable native inventory digest collision')
    sealed_qualification=None
    owner_qualification=None
    if generation_qualification:
        generation_qualification['status']='PASS' if inventory['owner_admissible'] else 'PARTIAL'
        generation_qualification['inventory_file_sha256']=digest
        generation_qualification['inventory_file']=str(target)
        generation_qualification['expected_outputs']=names
        qdata=(canonical(generation_qualification)+'\n').encode();qhash=hashlib.sha256(qdata).hexdigest()
        qpath=directory/(qhash+'.qualification.json')
        try:
            with qpath.open('xb') as stream:stream.write(qdata);stream.flush();os.fsync(stream.fileno())
        except FileExistsError:require(qpath.read_bytes()==qdata,'immutable qualification digest collision')
        sealed_qualification={'path':str(qpath.relative_to(fig)),'sha256':qhash}
        if inventory['owner_admissible']:
            require(bool(generation_qualification.get('mutation_cases')),'owner qualification requires executed dependency mutations')
            owner={'schema':'NativeDependencyQualificationV1','status':'PASS',
                'native_inventory_sha256':digest,'expected_outputs':names,
                'mutation_cases':generation_qualification['mutation_cases'],
                'source_root':str(source_root),'source_roots':source_roots(),'source_inventory':sources,
                'native_value_exports':observed,
                'renderer_coverage':{'path':str(registry_path),'sha256':registry_hash},
                'release_qualification':generation_qualification['release_qualification'],
                'generation_qualification':{'path':str(qpath),'sha256':qhash}}
            odata=(canonical(owner)+'\n').encode();ohash=hashlib.sha256(odata).hexdigest();opath=directory/(ohash+'.owner-qualification.json')
            try:
                with opath.open('xb') as stream:stream.write(odata);stream.flush();os.fsync(stream.fileno())
            except FileExistsError:require(opath.read_bytes()==odata,'immutable owner qualification collision')
            owner_qualification={'path':str(opath.relative_to(fig)),'sha256':ohash}
    pointer={'schema':'PlotPhysicsGenerationPointerV1','inventory':{'path':str(target.relative_to(fig)),'sha256':digest},
        'owner_admissible':inventory['owner_admissible'],'unavailable_reasons':inventory['unavailable_reasons'],
        'outputs':names,'generation_roster_sha256':roster_hash,'qualification':sealed_qualification,
        'owner_qualification':owner_qualification}
    atomic(fig/'PHYSICS_CURRENT.json',pointer)
    return pointer


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fig',type=Path,required=True)
    p.add_argument('--source-root',type=Path,required=True);p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--qualification',type=Path);p.add_argument('--input-root',type=Path,action='append',default=[])
    a=p.parse_args()
    try:print(json.dumps(collect(a.fig,a.source_root,a.registry,a.qualification,a.input_root),indent=2));return 0
    except (ValueError,OSError,KeyError,TypeError) as error:
        print('native generation collection failed: '+str(error),file=__import__('sys').stderr);return 2

if __name__=='__main__':raise SystemExit(main())
