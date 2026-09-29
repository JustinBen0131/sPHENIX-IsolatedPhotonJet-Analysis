"""Reduce eventTrees using one resolved selection for spectra, QA and response.

    python -m TreeToHists.reduce --inputs files.txt --out results/nominal
    python -m TreeToHists.reduce --inputs files.txt --selection my_selection.yaml --check

Input files are never modified. Each output source has its own three products,
so the existing source-bootstrap and FinalAnalysis readers keep their identity.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import uproot
import yaml

from .selection import DEFAULT, canonical, file_sha, number, resolve, verify_resolved

HERE = Path(__file__).resolve().parent


def inspect_inputs(paths, selection):
    if not paths or len(paths)!=len(set(paths)):
        raise ValueError('input list must be nonempty with no duplicate paths')
    family = None
    for path in paths:
        if any(c.isspace() for c in str(path)):
            raise ValueError('the maintained ROOT file-list interface does not accept whitespace in paths: '+str(path))
        with uproot.open(path) as f:
            info = json.loads(str(f['meta/schema']))
            source = json.loads(str(f['meta/source']))
            system = info['system']; spec = selection['systems'][system]
            dataset = source['dataset']
            kind = 'data' if not dataset.startswith('isSim') else 'jetsim' if 'Inclusive' in dataset else 'phosim'
            if family is not None and family != (system,kind):
                raise ValueError('one reduction must contain a single system and DATA/SIM family')
            family = system,kind
            tree = f['eventTree']
            for name in (spec['score_branch'],*([spec['state_branch']] if spec['state_branch'] else []),'pho_et','isow_cone_sum','isow_threshold','isow_sideband_threshold'):
                if name not in tree:
                    raise ValueError(str(path)+': missing selection input '+name)
            receipt = json.loads(str(f[spec['score_metadata']]))
            optional = spec.get('optional_model')
            if optional:
                score = receipt.get(optional['model_name'], {})
                if any(score.get(k) != optional[k] for k in ('release_sha256', 'score_branch', 'domain')) or score.get('status') != 'PASS':
                    raise ValueError(str(path)+': optional score is not bound to the selected release/domain')
            else:
                if receipt['model_sha256'] != spec['model_sha256'] or receipt.get('freeze_sha256') != spec['photon_id_package_sha256']:
                    raise ValueError(str(path)+': saved score differs from the selected feature/model package')
                if spec['score_branch'] != 'pho_bdt_score' or spec['state_branch'] != 'pho_bdt_state':
                    if any(receipt.get(k) != spec[k] for k in ('score_branch', 'state_branch')):
                        raise ValueError(str(path)+': augmentation receipt does not bind score/state branches')
            if spec['verify_nominal_isolation'] and system=='auau':
                iso = json.loads(str(f['meta/isolation']))
                if iso['freeze_sha256'] != spec['isolation']['package_sha256']:
                    raise ValueError(str(path)+': saved nominal AuAu isolation package differs')
    return family


def resolve_group(group):
    """Consume the immutable publication manifest; no producer implementation import."""
    import hashlib
    group = Path(group).resolve(); current = group / 'CURRENT'
    if not current.is_symlink(): raise ValueError('group is not published: ' + str(group))
    directory = current.resolve(strict=True)
    if directory.parent != group / 'versions': raise ValueError('group pointer escapes its version directory')
    manifest = json.loads((directory / 'MANIFEST.json').read_text())
    body = {k: v for k, v in manifest.items() if k != 'version'}
    if hashlib.sha256(canonical(body).encode()).hexdigest() != manifest['version'] or directory.name != manifest['version']:
        raise ValueError('group publication manifest changed')
    names = [f['relative_path'] for f in manifest['files']]
    identities = [f['logical_id'] for f in manifest['files']]
    if len(set(names)) != len(names) or len(set(identities)) != len(identities):
        raise ValueError('duplicate file or source in group manifest')
    for name in names:
        if Path(name).is_absolute() or '..' in Path(name).parts or Path(name).suffix != '.root':
            raise ValueError('unsafe publication file name')
    return directory, manifest


def published_inputs(group):
    """Resolve one immutable group version once before preparing any reader jobs."""
    directory,manifest=resolve_group(group)
    paths=[]
    for item in manifest['files']:
        p=directory/item['relative_path']
        if file_sha(p)!=item['sha256']:raise ValueError('published input differs: '+str(p))
        paths.append(p)
    return paths


def write_json(path, value):
    with Path(path).open('x') as f:
        f.write(canonical(value)+'\n'); f.flush(); os.fsync(f.fileno())


def worker(path):
    job = json.loads(Path(path).read_text())
    selection = verify_resolved(json.loads(Path(job['selection']).read_text()))
    for p,pin in job['code'].items():
        if file_sha(p)!=pin:
            raise ValueError('reducer code changed after preparation: '+p)
    if file_sha(job['input']) != job['input_sha256']:
        raise ValueError('input changed after preparation')
    inspect_inputs([Path(job['input'])],selection)
    import ROOT
    ROOT.gROOT.SetBatch(True)
    macro = HERE/'src'/('FillPhoton'+{'spectra':'Jet','qa':'QA','response':'JetResponse'}[job['kind']]+'.C')
    if not ROOT.gInterpreter.Declare('#include '+json.dumps(str(macro))):
        raise RuntimeError('ROOT failed to load reducer')
    inp,out,sel = job['input'],job['output'],job['selection']
    weights = ' '.join(map(str,job['weights']))
    if job['kind']=='spectra':
        code = ROOT.FillPhotonJet(inp,out,'nominal',sel)
    elif job['kind']=='qa':
        code = ROOT.FillPhotonQA(inp,out,'nominal',weights,sel)
    else:
        edges = [' '.join(map(str,job['response'].get(k,[]))) for k in ('reco_ptgamma_edges','truth_ptgamma_edges','xj_edges')]
        code = ROOT.FillPhotonJetResponse(inp,out,'nominal',*edges,weights,'',sel)
    if code:
        raise RuntimeError('ROOT reducer returned '+str(code))
    with uproot.open(out) as f:
        if str(f['profile']) != selection['profile'] or verify_resolved(json.loads(str(f['selection']))) != selection:
            raise ValueError('saved output selection disagrees with requested definition')
    if file_sha(inp)!=job['input_sha256'] or any(file_sha(p)!=pin for p,pin in job['code'].items()):
        raise ValueError('input or reducer code changed during reduction')


def run(inputs, out, selection_path=DEFAULT, analysis_path=None, check=False):
    selection = resolve(selection_path)
    paths = [];pinned={}
    if Path(inputs).is_dir():paths=published_inputs(inputs)
    else:
        for line in Path(inputs).read_text().splitlines():
            value = line.strip()
            if value and not value.startswith('#'):
                path=Path(inputs).resolve().parent/value
                text=str(path)
                if '/CURRENT/' in text:
                    group,relative=text.split('/CURRENT/',1)
                    if group not in pinned:
                        pinned[group]=resolve_group(group)[0]
                    path=pinned[group]/relative
                paths.append(path.resolve())
    system,kind = inspect_inputs(paths,selection)
    response = {}
    analysis = None
    if analysis_path:
        analysis = yaml.safe_load(Path(analysis_path).read_text())
        response = dict(analysis.get('response') or {})
        # xJ belongs to measurement in the existing FinalAnalysis YAML. Leaving
        # it at the ROOT default can make its nonuniform edges impossible to read.
        if 'xj_edges' not in response and 'xj_edges' in analysis.get('measurement',{}):
            response['xj_edges'] = analysis['measurement']['xj_edges']
    for key in ('reco_ptgamma_edges','truth_ptgamma_edges','xj_edges'):
        if key not in response:
            continue
        edges = response[key]
        if not isinstance(edges,list) or len(edges)<2:
            raise ValueError(key+' needs at least two edges')
        edges = [number(v,key) for v in edges]
        if any(a>=b for a,b in zip(edges,edges[1:])):
            raise ValueError(key+' must increase strictly')
        response[key] = edges
    weights = response.get('auau_centrality_weights') if system=='auau' and kind!='data' else []
    if system=='auau' and kind!='data' and (not isinstance(weights,list) or len(weights)!=16 or
         any(isinstance(v,bool) or not isinstance(v,(float,int)) or not 0<v<float('inf') for v in weights)):
        raise ValueError('AuAu SIM requires its 16 finite positive centrality weights through --analysis-config')
    summary = dict(profile=selection['profile'],name=selection['name'],system=system,family=kind,files=len(paths),
                   products=['spectra','qa']+(['response'] if kind!='data' else []))
    if check:
        return summary
    out = Path(out).resolve(); out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic(); jobs=[]; done=[]
    def progress(error=None):
        value=dict(processed=len(done),total=len(jobs),remaining=len(jobs)-len(done),
                   phase='failed' if error else 'complete' if jobs and len(done)==len(jobs) else 'reducing',
                   updated_at=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),elapsed_seconds=time.monotonic()-started,
                   first_error=error)
        temp=out/'progress.next.json'; temp.write_text(canonical(value)+'\n'); os.replace(temp,out/'progress.json')
    try:
        sel=out/'selection.json';write_json(sel,selection)
        if analysis is not None:
            # FinalAnalysis already checks profile equality across DATA and SIM.
            # Supply the matching configuration without asking users to copy a hash.
            bound_analysis = dict(analysis,profile=selection['profile'])
            with (out/'analysis.yaml').open('x') as f:
                yaml.safe_dump(bound_analysis,f,sort_keys=False)
        code={str(p):file_sha(p) for p in sorted((HERE/'src').iterdir()) if p.suffix in ('.h','.C')}
        code.update({str(p):file_sha(p) for p in (HERE/'selection.py',HERE/'reduce.py')})
        write_json(out/'invocation.json',dict(**summary,inputs=[str(p) for p in paths],code=code))
        # Preserve one physical input per response output for source bootstrap.
        for i,path in enumerate(paths):
            folder=out/('source_%06d'%i);folder.mkdir()
            pin=file_sha(path)
            for product in summary['products']:
                job=dict(input=str(path),input_sha256=pin,output=str(folder/(product+'.root')),kind=product,
                         selection=str(sel),weights=weights,response=response,code=code)
                j=folder/(product+'.job.json');write_json(j,job);jobs.append(j)
        progress()
        for j in jobs:
            with j.with_suffix('.log').open('x') as log:
                result=subprocess.run([sys.executable,'-B','-m','TreeToHists.reduce','--worker',str(j)],
                    cwd=HERE.parent,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError('reduction failed; inspect '+str(j.with_suffix('.log')))
            job=json.loads(j.read_text()); done.append(dict(path=job['output'],sha256=file_sha(job['output']),
                input_sha256=job['input_sha256'],kind=job['kind']));progress()
        result=dict(**summary,status='PASS',outputs=done,selection_sha256=selection['selection_sha256'])
        write_json(out/'COMPLETE.json',result)
        return result
    except Exception as exc:
        progress(str(exc));write_json(out/'FAILED.json',dict(error=str(exc),completed_products=done))
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--inputs',type=Path,help='text file with one eventTree path per line; relative to the list')
    p.add_argument('--out',type=Path,help='new output directory; existing directories are refused')
    p.add_argument('--selection',type=Path,default=DEFAULT,help='readable selection YAML (default: config/selection.yaml)')
    p.add_argument('--analysis-config',type=Path,help='existing analysis YAML supplying response axes and AuAu SIM centrality weights')
    p.add_argument('--check',action='store_true',help='validate configuration and input metadata without writing outputs')
    p.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    a=p.parse_args()
    try:
        if a.worker:
            worker(a.worker);return 0
        if not a.inputs or (not a.check and not a.out):
            p.error('--inputs and --out are required (omit --out with --check)')
        result=run(a.inputs,a.out,a.selection,a.analysis_config,a.check)
        print(json.dumps(result,indent=2));return 0
    except Exception as exc:
        print(str(exc),file=sys.stderr);return 1


if __name__=='__main__':
    raise SystemExit(main())
