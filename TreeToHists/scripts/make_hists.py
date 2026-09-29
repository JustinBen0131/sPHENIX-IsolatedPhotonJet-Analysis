"""Plan, reduce, resume and publish one immutable ROOT file per exact sample.

Only this command coordinates files. Physics definitions and event loops live in
the maintained C++ fillers; sample identities and normalizations come from the
selected input collection. A completed sample is reused only with identical pins.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import time
import uuid

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent))
import numpy as np
import uproot
import yaml
from TreeToHists.selection import DEFAULT, canonical, digest, file_sha, resolve, verify_resolved, AUAU_DATA_TRIGGERS
from TreeToHists.reduce import inspect_inputs, published_inputs
from TreeToHists.histogram_contract import inventory

COMMANDS = ('runPP', 'runAuAu', 'runAllData', 'runAllSim', 'runAll')


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save(path, value):
    """Atomic local state update. Immutable results use exclusive publication."""
    path = Path(path)
    temp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    with temp.open('x') as stream:
        stream.write(canonical(value) + '\n'); stream.flush(); os.fsync(stream.fileno())
    os.replace(temp, path)


def definitions(path):
    h = yaml.safe_load(Path(path).read_text())
    if h.get('version') != 1:
        raise ValueError('unsupported histogram contract version')
    fixed = {'native_measurement_et_edges_gev': list(range(15,36)),
             'native_qa_et_edges_gev': list(range(10,41)),
             'jet_thresholds_gev': [5,7,10,12], 'reported_photon_et_edges_gev': [17,23,33]}
    for key, value in fixed.items():
        if h.get(key) != value:
            raise ValueError(key + ' differs from the qualified filler contract')
    for key in ('native_xj_edges','native_reco_response_edges_gev','native_truth_response_edges_gev'):
        a = np.asarray(h[key], dtype=float)
        if len(a)<2 or not np.isfinite(a).all() or not (np.diff(a)>0).all():
            raise ValueError('invalid storage axis: '+key)
    if h['native_xj_edges'][0]!=0 or h['native_xj_edges'][-1]!=3:
        raise ValueError('xJ storage support must be [0,3]')
    for key, bounds in (('native_reco_response_edges_gev',(10,40)),('native_truth_response_edges_gev',(5,40))):
        if (h[key][0],h[key][-1]) != bounds or not set([15,17,23,33,35]).issubset(h[key]):
            raise ValueError(key+' must retain response support and all photon interval edges')
    if not h.get('extensions'):
        raise ValueError('this public contract requires its extension pass')
    return h


def source_info(path):
    with uproot.open(path) as f:
        tree=f['eventTree']; source=json.loads(str(f['meta/source']))
        meta={'source':source,'entries':tree.num_entries}
        for key in ('trigger_selection','trigger_run_info','normalization','weighting','event_accounting','accounting','consolidation'):
            if 'meta/'+key in f:
                value=json.loads(str(f['meta/'+key]))
                # Consolidation's original metadata stays in the immutable input;
                # retain source/chunk identity, avoiding repeated full dictionaries.
                if key=='consolidation':
                    value={k:v for k,v in value.items() if k!='chunks'}
                    value['input_metadata_sha256']=digest(json.loads(str(f['meta/'+key])))
                meta[key]=value
    ident=source.get('logical_source_id') or source.get('input_file_sha256')
    if not isinstance(ident,str) or not ident:
        raise ValueError('missing original source identity: '+str(path))
    return dict(path=str(path),sha256=file_sha(path),identity=ident,metadata=meta)


def plan(args):
    selection=resolve(args.selection,args.data_trigger); histograms=definitions(args.histograms)
    collection=Path(args.inputs or os.environ.get('PHOTONJET_TREES','')).expanduser()
    if not collection.is_file():
        raise ValueError('supply --inputs collection.json or PHOTONJET_TREES; canonical collection locations are selected by the user')
    catalogue=json.loads(collection.read_text())
    if catalogue.get('schema')!='PhotonJetInputCollectionV1' or not catalogue.get('samples'):
        raise ValueError('expected PhotonJetInputCollectionV1 with nonempty samples')
    code={str(p.relative_to(HERE)):file_sha(p) for p in sorted(HERE.rglob('*'))
          if p.is_file() and 'tests' not in p.parts and p.suffix in ('.py','.h','.cc','.C') and '__pycache__' not in p.parts}
    samples=[]; keys=set(); identities=set()
    for sample in catalogue['samples']:
        system,name,family=(sample[k] for k in ('system','name','family'))
        if system not in ('PP','AuAu') or family not in ('data','photonjet','inclusivejet') or not re.fullmatch('[A-Za-z0-9_-]+',name):
            raise ValueError('invalid system, sample name or family')
        if (name=='data') != (family=='data'):
            raise ValueError('each system has one data sample named data')
        selected=(args.command=='runAll' or args.command=='runAllSim' and family!='data' or
                  args.command=='runAllData' and family=='data' or
                  args.command=='runPP' and system=='PP' and family=='data' or
                  args.command=='runAuAu' and system=='AuAu' and family=='data')
        if not selected or args.sample and name not in args.sample:
            continue
        key=system+'/'+name
        if key in keys: raise ValueError('duplicate sample: '+key)
        keys.add(key)
        paths=[]
        for value in sample.get('files',[]):
            paths.append((collection.resolve().parent/value).resolve())
        for value in sample.get('groups',[]):
            paths.extend(published_inputs((collection.resolve().parent/value).resolve()))
        actual=inspect_inputs(paths,selection)
        if actual!=(system.lower(),{'data':'data','photonjet':'phosim','inclusivejet':'jetsim'}[family]):
            raise ValueError('sample identity disagrees with eventTree metadata: '+key)
        sources=[source_info(p) for p in paths]
        for s in sources:
            if s['identity'] in identities: raise ValueError('duplicate original source identity: '+s['identity'])
            identities.add(s['identity'])
        populations={s['metadata']['source'].get('sample') for s in sources}
        if family!='data' and (len(populations)!=1 or None in populations):
            raise ValueError('SIM package must contain exactly one named source population')
        if sample.get('source_sample') and populations!={sample['source_sample']}:
            raise ValueError('source_sample binding differs: '+key)
        weights=sample.get('centrality_weights',[])
        vertex=None
        if system=='AuAu' and family!='data':
            if len(weights)!=16 or not np.isfinite(weights).all() or not (np.asarray(weights)>0).all():
                raise ValueError(key+' needs its 16 production-derived centrality_weights')
            if not sample.get('centrality_evidence'):
                raise ValueError(key+' needs centrality_evidence identifying the production-specific derivation')
            expected_trigger=AUAU_DATA_TRIGGERS[selection['systems']['auau']['data_trigger']]
            if sample.get('centrality_data_trigger')!=expected_trigger:
                raise ValueError(key+' centrality weights must bind the selected DATA trigger: '+expected_trigger)
            vertex=sample.get('vertex_reweighting')
            if not isinstance(vertex,dict) or not vertex.get('evidence'):
                raise ValueError(key+' needs a qualified AuAu vertex_reweighting fit and closure evidence')
            vertex=dict(vertex)
            if (vertex.get('coordinate')!='truth_vtx_z' or
                vertex.get('selection_sha256')!=selection['selection_sha256'] or
                vertex.get('centrality_weights_sha256')!=digest(weights)):
                raise ValueError(key+' vertex fit must bind truth_vtx_z, the selection and centrality weights')
            fit=(collection.resolve().parent/vertex.get('path','')).resolve()
            if not fit.is_file() or file_sha(fit)!=vertex.get('sha256'):
                raise ValueError(key+' vertex fit is missing or its SHA256 differs')
            vertex['path']=str(fit)
            with uproot.open(fit) as f:
                hist=f.get(vertex.get('histogram',''))
                if hist is None or not hist.classname.startswith('TH1'):
                    raise ValueError(key+' vertex fit requires a one-dimensional histogram')
                values=hist.values()
                if len(values)<2 or not np.isfinite(values).all() or (values<0).any():
                    raise ValueError(key+' vertex fit contains invalid weights')
            if any(s['metadata'].get('weighting',{}).get('vertex_enabled') for s in sources):
                raise ValueError(key+' input already includes vertex weights')
        elif weights: raise ValueError('centrality weights apply only to AuAu SIM')
        elif sample.get('vertex_reweighting'):
            raise ValueError('additional vertex weights apply only to AuAu SIM; pp uses its saved factor')
        normalization=sample.get('normalization',{})
        if family!='data':
            generated=normalization.get('generated_events')
            if not isinstance(generated,(int,float)) or not np.isfinite(generated) or generated<=0:
                raise ValueError(key+' needs a positive generated-event count')
            if normalization.get('source_factor')=='family_reference':
                cross_section=catalogue.get('reference_cross_section_pb',{}).get(normalization.get('reference'))
            elif normalization.get('source_factor')=='cross_section_over_generated':
                cross_section=normalization.get('cross_section_pb')
            else: raise ValueError(key+' has an unsupported normalization rule')
            if not isinstance(cross_section,(int,float)) or not np.isfinite(cross_section) or cross_section<=0:
                raise ValueError(key+' needs a positive bound cross section')
            normalization=dict(normalization,numerator_pb=cross_section,applied_in_histograms=False)
        rows=inventory(system,'DATA' if family=='data' else family)
        for row in rows:
            row['kind']={'measurement':'spectra','response':'response','qa':'qa','extension':'extensions'}[row['group']]
            row['key']=row['old_key'] or row['path'].replace('/','__')
        samples.append(dict(key=key,system=system,name=name,family=family,sources=sources,weights=weights,vertex_reweighting=vertex,
                            normalization=normalization,centrality_evidence=sample.get('centrality_evidence'),
                            centrality_data_trigger=sample.get('centrality_data_trigger'),inventory=rows))
    if not samples: raise ValueError('no samples match the command/selection')
    library=os.environ.get('PHOTONJET_HISTS_LIBRARY')
    library=dict(path=str(Path(library).resolve()),sha256=file_sha(library)) if library else None
    runtime=dict(python=sys.version,numpy=np.__version__,uproot=uproot.__version__,library=library)
    value=dict(schema='PhotonJetReductionPlanV1',command=args.command,selection=selection,histograms=histograms,runtime=runtime,
               collection_sha256=file_sha(collection),code=code,samples=samples)
    value['plan_sha256']=digest(value)
    return value


def verify_code(code):
    for name,pin in code.items():
        if file_sha(HERE/name)!=pin: raise ValueError('code changed after preparation: '+name)


def root_job(path):
    """One isolated C++ worker, with an exact job, source and code receipt."""
    job=json.loads(Path(path).read_text()); start=time.monotonic()
    verify_code(job['code'])
    for source in job['pins']:
        if file_sha(source['path'])!=source['sha256']: raise ValueError('changed input: '+source['path'])
    import ROOT
    ROOT.gROOT.SetBatch(True)
    # Interpretation is the portable local fallback; deployments may point at the
    # same sources built with src/autogen.sh. No dependency install is performed.
    sources=['PhotonJetHists.cc','internal/Selection.cc','internal/Histograms.cc','internal/Package.cc','internal/Centrality.cc']
    library=job.get('library')
    if library:
        if file_sha(library['path'])!=library['sha256']: raise ValueError('compiled library changed')
        if ROOT.gSystem.Load(library['path'])<0: raise RuntimeError('cannot load '+library['path'])
        text='#include '+json.dumps(str(HERE/'src/PhotonJetHists.h'))
    else:
        text='\n'.join('#include '+json.dumps(str(HERE/'src'/p)) for p in sources)
    if not ROOT.gInterpreter.Declare(text): raise RuntimeError('C++ histogram package failed to load')
    if ROOT.photonjet.runHistogramJob(str(path)):
        raise RuntimeError('C++ histogram job refused completion; see the preceding error')
    verify_code(job['code'])
    for source in job['pins']:
        if file_sha(source['path'])!=source['sha256']: raise ValueError('input changed while processing')
    with uproot.open(job['output']) as f:
        if str(f['profile'])!=job['profile']: raise ValueError('saved selection profile differs')
        verify_resolved(json.loads(str(f['selection'])))
    save(job['receipt'],dict(status='PASS',job_sha256=file_sha(path),output=job['output'],
        sha256=file_sha(job['output']),elapsed_seconds=time.monotonic()-start,
        peak_rss_platform_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,root_version=ROOT.gROOT.GetVersion(),completed_at=now()))


def execute(job, directory, deadline):
    directory.mkdir(parents=True,exist_ok=True)
    recipe=digest(job); receipt=directory/'PASS.json'
    if receipt.exists():
        r=json.loads(receipt.read_text()); old=json.loads((directory/'job.json').read_text())
        if (old.get('recipe_sha256')!=recipe or old['output']!=r['output'] or
            file_sha(directory/'job.json')!=r['job_sha256'] or file_sha(r['output'])!=r['sha256']):
            raise ValueError('completed product has changed; select a new output directory: '+str(directory))
        return r
    # Preserve interrupted candidates for diagnosis. A resumed attempt is new.
    attempt=directory/('attempt_'+uuid.uuid4().hex[:12]);attempt.mkdir()
    expanded=dict(job,recipe_sha256=recipe,output=str(attempt/'product.root'),receipt=str(attempt/'PASS.json'),progress=str(attempt/'progress.json'))
    path=attempt/'job.json';save(path,expanded)
    with (attempt/'worker.log').open('x') as log:
        result=subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),'--worker',str(path)],
                              cwd=HERE.parent,stdout=log,stderr=subprocess.STDOUT,timeout=deadline)
    if result.returncode:
        save(attempt/'FAILED.json',dict(status='FAILED',returncode=result.returncode,job_sha256=file_sha(path),
             first_error='worker refused completion; see worker.log',failed_at=now()))
        raise RuntimeError('worker failed; inspect '+str(attempt/'worker.log'))
    r=json.loads((attempt/'PASS.json').read_text())
    save(directory/'job.json',expanded);save(receipt,r)
    return r


def check_package(path,sample,selection):
    """Saved-file validation: inventory, uncertainty cells and event references."""
    with uproot.open(path) as f:
        meta=json.loads(str(f['meta/package']))
        if meta['source_count']!=len(sample['sources']) or meta['profile']!=selection['profile']:
            raise ValueError('sample source count/profile changed')
        for row in sample['inventory']:
            if row['kind']=='extensions' and row['path'] not in f: continue
            if row['path'] not in f: raise ValueError('missing package object: '+row['path'])
            h=f[row['path']]
            if h.classname.startswith('TH') and not h.classname.startswith('THn'):
                if not np.isfinite(h.values(flow=True)).all() or h.variances(flow=True) is None or not np.isfinite(h.variances(flow=True)).all():
                    raise ValueError('invalid histogram payload: '+row['path'])
        roots=['records']
        if sample['system']=='AuAu' and sample['family']=='data':
            roots+=['records/triggers/'+k for k in AUAU_DATA_TRIGGERS]
        for root in roots:
            n=f[root+'/events'].num_entries
            for key in ('photons','jets'):
                for arrays in f[root+'/'+key].iterate(['event'],step_size='8 MB',library='np'):
                    if np.any(arrays['event']<0) or np.any(arrays['event']>=n): raise ValueError('broken event reference')
        if f['records/sources'].num_entries!=len(sample['sources']): raise ValueError('source records lost')
    return meta


@contextmanager
def lock(output):
    with (output/'.make_hists.lock').open('a') as f:
        try: fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('another coordinator owns this output directory')
        yield


def run(args, prepared):
    out=args.out.expanduser().resolve();out.mkdir(parents=True,exist_ok=True)
    with lock(out):
        work=out/'.work'; work.mkdir(exist_ok=True)
        plan_path=work/'plan.json'
        if plan_path.exists():
            previous=json.loads(plan_path.read_text())
            if previous!=prepared: raise ValueError('output is bound to a different plan; choose another --out')
            if not args.resume: raise ValueError('existing plan: use --resume to verify and reuse it')
        else: save(plan_path,prepared)
        sel=work/'selection.json';save(sel,prepared['selection'])
        manifest_path=out/'manifest.json'
        manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else dict(
            schema='PhotonJetHistogramCollectionV1',plan_sha256=prepared['plan_sha256'],samples={},created_at=now())
        products_done=0
        total=sum(len(s['sources'])*(3 if s['family']=='data' else 4)+1 for s in prepared['samples'])
        def progress(phase,active=None,error=None):
            save(out/'progress.json',dict(phase=phase,active=active,processed=products_done,total=total,
                remaining=total-products_done,updated_at=now(),first_error=error,
                samples_complete=len(manifest['samples']),samples_total=len(prepared['samples'])))
        try:
            progress('reducing')
            for sample in prepared['samples']:
                key=sample['key']; parts=[]
                if key in manifest['samples']:
                    if file_sha(out/manifest['samples'][key]['path'])!=manifest['samples'][key]['sha256']:
                        raise ValueError('published sample changed: '+key)
                    products_done+=len(sample['sources'])*(3 if sample['family']=='data' else 4)+1
                    progress('reducing');continue
                for i,source in enumerate(sample['sources']):
                    products=[]
                    for kind in ['spectra','qa']+(['response'] if sample['family']!='data' else [])+['extensions']:
                        progress('reducing',key+'/'+str(i)+'/'+kind)
                        job=dict(operation='reduce',input=source['path'],entries=source['metadata']['entries'],kind=kind,
                            selection=str(sel),profile=prepared['selection']['profile'],weights=sample['weights'],
                            histograms=prepared['histograms'],code=prepared['code'],
                            vertex_reweighting=sample['vertex_reweighting'],
                            pins=[source]+([sample['vertex_reweighting']] if sample['vertex_reweighting'] else []),library=prepared['runtime']['library'])
                        folder=work/key/('%03d'%(i//1000))/('%06d'%i)/kind
                        result=execute(job,folder,args.timeout)
                        products.append(dict(kind=kind,path=result['output'],sha256=result['sha256']))
                        products_done+=1
                    parts.append(dict(metadata=source,products=products))
                progress('packaging',key)
                metadata=dict(schema='PhotonJetSampleV1',system=sample['system'],name=sample['name'],family=sample['family'],
                    profile=prepared['selection']['profile'],selection=prepared['selection'],histograms=prepared['histograms'],
                    normalization=sample['normalization'],centrality_weights=sample['weights'],centrality_evidence=sample['centrality_evidence'],
                    centrality_data_trigger=sample['centrality_data_trigger'],vertex_reweighting=sample['vertex_reweighting'],
                    plan_sha256=prepared['plan_sha256'],runtime=prepared['runtime'],acceptance='LOCAL_STRUCTURE_ONLY_REQUIRES_PHYSICS_QUALIFICATION')
                job=dict(operation='package',sources=parts,metadata=metadata,inventory=sample['inventory'],
                    code=prepared['code'],pins=[p for part in parts for p in part['products']],profile=metadata['profile'],library=prepared['runtime']['library'])
                result=execute(job,work/key/'package',args.timeout)
                meta=check_package(result['output'],sample,prepared['selection'])
                target=out/(key+'.root');target.parent.mkdir(exist_ok=True)
                # Exclusive hard-link publication is atomic on this output volume.
                # A crash after link but before manifest is reconciled by exact hash.
                if target.exists():
                    if file_sha(target)!=result['sha256']: raise ValueError('publication collision: '+str(target))
                else: os.link(result['output'],target)
                manifest['samples'][key]=dict(path=key+'.root',sha256=result['sha256'],sources=len(sample['sources']),
                    histograms=meta['histograms_written'],availability=meta['availability'],published_at=now())
                manifest['updated_at']=now();save(manifest_path,manifest);products_done+=1;progress('reducing')
            progress('complete');return manifest
        except Exception as exc:
            progress('failed',error=str(exc));raise


def main():
    p=argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('command',choices=COMMANDS,nargs='?')
    p.add_argument('--inputs',help='published input collection JSON; default PHOTONJET_TREES')
    p.add_argument('--out',type=Path,default=Path('histograms'),help='user-owned output directory (default ./histograms)')
    p.add_argument('--sample',action='append',help='exact sample name; repeat to select several')
    p.add_argument('--selection',type=Path,default=DEFAULT)
    p.add_argument('--data-trigger',choices=tuple(AUAU_DATA_TRIGGERS),help='nominal AuAu scaled trigger (default photon12); both streams are retained')
    p.add_argument('--histograms',type=Path,default=HERE/'config/histograms.yaml')
    p.add_argument('--plan',action='store_true',help='verify inputs and show counts without running workers or writing output')
    p.add_argument('--resume',action='store_true',help='reuse only unchanged products with matching PASS receipts')
    p.add_argument('--timeout',type=int,default=3600,help='hard deadline in seconds per isolated worker')
    p.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    a=p.parse_args()
    try:
        if a.worker: root_job(a.worker);return 0
        if not a.command: p.error('a run command is required; see --help')
        if a.timeout<1: p.error('--timeout must be positive')
        prepared=plan(a)
        if a.plan:
            print(json.dumps(dict(plan_sha256=prepared['plan_sha256'],profile=prepared['selection']['profile'],
                samples=[dict(key=s['key'],files=len(s['sources']),events=sum(v['metadata']['entries'] for v in s['sources']),
                              target_histograms=len(s['inventory'])) for s in prepared['samples']]),indent=2))
        else: print(json.dumps(run(a,prepared),indent=2))
        return 0
    except Exception as exc:
        print(str(exc),file=sys.stderr);return 1


if __name__=='__main__':
    raise SystemExit(main())
