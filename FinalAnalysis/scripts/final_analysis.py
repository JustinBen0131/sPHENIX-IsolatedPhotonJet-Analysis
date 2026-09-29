#!/usr/bin/env python3
"""Public analysis commands over a saved histogram collection.

Preparation verifies the collection and resolves recorded normalization without
changing physics settings. Execution calls the maintained numerical kernels
once, saves their results and figure inputs, and records exact reuse evidence.
Rendering reads those saved inputs; it never calls the measurement calculation.
"""
from __future__ import annotations

import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys

import yaml

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent))
from FinalAnalysis.src import analysis
from FinalAnalysis.src.internal import inputs, systematics


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha(path):
    return analysis.sha256(Path(path))


def save(path, value):
    """Immutable small receipt; repeated identical writes are harmless."""
    path = Path(path); text = canonical(value) + '\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open('x') as stream:
            stream.write(text); stream.flush(); os.fsync(stream.fileno())
    except FileExistsError:
        if path.read_text() != text: raise ValueError('existing artifact differs: ' + str(path))


def progress(out, phase, done, total, **extra):
    path = Path(out) / 'progress.json'; pending = path.with_suffix('.pending')
    pending.write_text(canonical(dict(phase=phase, processed=done, total=total, remaining=total-done,
                                    updated_at=dt.datetime.now(dt.timezone.utc).isoformat(), **extra)) + '\n')
    os.replace(pending, path)


def prepare(args):
    collection = inputs.read_collection(args.manifest)
    config = yaml.safe_load(args.config.read_text()); resolved = copy.deepcopy(config)
    physical = {'runPP': {'pp'}, 'runAuAu': {'auau'}, 'runAll': {'pp', 'auau'}}[args.command]
    keys = [k for k in config['systems'] if analysis.physical_system(config, k) in physical]
    if not keys: raise ValueError('configuration has no requested systems')
    data, simulation, profiles = {}, {}, set()
    for item in collection.values():
        meta = item['metadata']; system = {'PP': 'pp', 'AuAu': 'auau'}[meta['system']]
        if system not in physical: continue
        profiles.add(meta['profile'])
        if config['profile'] == 'nominal':
            selection = meta['selection']['systems'][system]
            if not selection['verify_nominal_id'] or not selection['verify_nominal_isolation']:
                raise ValueError('nominal configuration received a varied photon selection')
        if meta['family'] == 'data':
            data.setdefault(system, []).append(str(item['path']))
        else:
            name = meta['name'] if meta['name'].startswith(system + '_') else system + '_' + meta['name']
            if name in simulation: raise ValueError('duplicate named simulation sample: ' + name)
            simulation[name] = [str(item['path'])]
            n = meta['normalization']
            if n.get('applied_in_histograms') is not False:
                raise ValueError('normalization must be recorded as unapplied: ' + name)
            target = resolved['samples'].get(name)
            if target is None: continue  # QA-only samples need no measurement weight.
            for field in ('generated_events', 'cross_section_pb', 'source_factor', 'reference'):
                if field in n:
                    if target.get(field) is not None and target[field] != n[field]:
                        raise ValueError('configured and recorded normalization differ: ' + name + '.' + field)
                    target[field] = n[field]
            if n.get('source_factor') == 'family_reference':
                numerator = n.get('numerator_pb')
                if numerator is None and n.get('weight_pb_per_event') is not None:
                    # Older packages record the evaluated factor. Keep their
                    # exact numerator in config; never reconstruct it by an
                    # inverse floating-point multiplication.
                    numerator = resolved.get('reference_cross_section_pb', {}).get(n['reference'])
                reference = resolved.setdefault('reference_cross_section_pb', {})
                if numerator is None: raise ValueError('missing exact reference cross section: ' + name)
                if reference.get(n['reference']) not in (None, numerator):
                    raise ValueError('reference cross section differs from the collection')
                reference[n['reference']] = numerator
            inputs.sample_weight(target, resolved.get('reference_cross_section_pb') or {})
            if system == 'auau':
                weights = meta['centrality_weights']
                if len(weights) != 16 or not meta.get('centrality_evidence'):
                    raise ValueError('AuAu response has no source-bound centrality weights')
                if resolved['response'].get('auau_centrality_weights') not in (None, weights):
                    raise ValueError('AuAu sample centrality weights disagree')
                resolved['response']['auau_centrality_weights'] = weights
    if len(profiles) != 1: raise ValueError('requested samples must share one selection profile')
    profile = profiles.pop()
    if config['profile'] not in ('nominal', profile):
        raise ValueError('configured profile differs from the saved collection')
    resolved['profile'] = profile
    jobs = []
    for key in keys:
        system = analysis.physical_system(resolved, key)
        if not data.get(system): raise ValueError('missing DATA sample for ' + system)
        missing = set(config['systems'][key]['samples']) - set(simulation)
        if missing: raise ValueError('missing response samples: ' + ', '.join(sorted(missing)))
        analysis.Settings.from_config(resolved, key)
        jobs.append(dict(system=key, data=data[system],
                         simulation={n: simulation[n] for n in config['systems'][key]['samples']}))
    pins = {str(item['path']): item['sha256'] for item in collection.values()}
    code = {str(p.relative_to(HERE)): sha(p) for folder in ('src', 'plots', 'scripts')
            for p in sorted((HERE / folder).rglob('*.py'))}
    plan = dict(schema='PhotonJetAnalysisPlanV1', command=args.command, config_sha256=sha(args.config),
                manifest_sha256=sha(args.manifest), resolved_config=resolved, jobs=jobs, inputs=pins, code=code)
    plan['sha256'] = hashlib.sha256(canonical(plan).encode()).hexdigest()
    return plan


def verify_pins(plan):
    for path, pin in plan['inputs'].items():
        if sha(path) != pin: raise ValueError('histogram input changed: ' + path)
    for name, pin in plan['code'].items():
        if sha(HERE / name) != pin: raise ValueError('analysis code changed: ' + name)


def execute(args, plan):
    from FinalAnalysis.plots import figure_package, native_response, render_qa
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    with (out / '.analysis.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (out / 'PLAN.json').exists() and not args.resume:
            raise ValueError('analysis already prepared; use --resume to verify and reuse it')
        save(out / 'PLAN.json', plan)
        config_path = out / '.work/resolved_config.json'; save(config_path, plan['resolved_config'])
        done = 0; completed = {}
        try:
            for job in plan['jobs']:
                key = job['system']; receipt_path = out / 'receipts' / (key + '.json')
                progress(out, 'analysis', done, len(plan['jobs']), active=key)
                if receipt_path.exists():
                    receipt = json.loads(receipt_path.read_text())
                    if receipt['plan_sha256'] != plan['sha256']:
                        raise ValueError('result receipt belongs to another plan')
                    for name, pin in receipt['files'].items():
                        if sha(out / name) != pin: raise ValueError('saved result changed: ' + name)
                    verify_pins(plan); completed[key] = receipt; done += 1; continue
                result_path = out / 'results' / (key + '.json')
                package_path = out / 'packages' / (key + '.npz')
                if any(p.exists() for p in (result_path, result_path.with_suffix('.csv'), package_path, package_path.with_suffix('.json'))):
                    raise ValueError('interrupted output lacks its receipt; preserve it and choose a fresh --out: ' + key)
                result_path.parent.mkdir(parents=True, exist_ok=True); package_path.parent.mkdir(parents=True, exist_ok=True)
                record = {}; verify_pins(plan)
                document, result, settings = analysis.compute(config_path, key, job['data'], job['simulation'], record=record)
                native = native_response.from_files(analysis, job['simulation'], plan['resolved_config'], key)
                package, checks = figure_package.build(result, record, key, analysis._prior.per_photon,
                                                       native=native, restrict=analysis._closure.restrict)
                save(result_path, analysis.jsonable(document))
                analysis.write_points(result, settings.xj_edges, result_path.with_suffix('.csv'))
                figure_package.write(package_path, package, checks, dict(schema='IANFigurePackageV1', label=key,
                                     result=str(result_path), result_sha256=sha(result_path), plan_sha256=plan['sha256']))
                verify_pins(plan)
                files = [result_path, result_path.with_suffix('.csv'), package_path, package_path.with_suffix('.json')]
                receipt = dict(schema='PhotonJetAnalysisReceiptV1', status='PASS', plan_sha256=plan['sha256'],
                               files={str(p.relative_to(out)): sha(p) for p in files}, figure_checks=checks)
                save(receipt_path, receipt); completed[key] = receipt; done += 1
            qa_path = out / 'packages/qa_inputs.npz'; qa_receipt = out / 'receipts/qa.json'
            if qa_receipt.exists():
                receipt = json.loads(qa_receipt.read_text())
                if receipt['plan_sha256'] != plan['sha256']: raise ValueError('QA receipt belongs to another plan')
                for name, pin in receipt['files'].items():
                    if sha(out / name) != pin: raise ValueError('saved QA input changed: ' + name)
            else:
                verify_pins(plan)
                render_qa.prepare(args.manifest, plan['resolved_config'], qa_path,
                                  {analysis.physical_system(plan['resolved_config'], j['system']) for j in plan['jobs']})
                verify_pins(plan)
                save(qa_receipt, dict(schema='PhotonJetQAReceiptV1', plan_sha256=plan['sha256'],
                     files={str(p.relative_to(out)): sha(p) for p in (qa_path, qa_path.with_suffix('.json'))}))
            save(out / 'RESULTS.json', dict(schema='PhotonJetAnalysisResultsV1', plan_sha256=plan['sha256'], systems=completed,
                                          qa=json.loads(qa_receipt.read_text())))
            progress(out, 'complete', done, len(plan['jobs']))
        except Exception as exc:
            progress(out, 'failed', done, len(plan['jobs']), first_error=str(exc)); raise
    return dict(status='PASS', output=str(out), systems=list(completed))


def render(args):
    from FinalAnalysis.plots import render as renderer, tables, render_qa
    results = json.loads((args.results / 'RESULTS.json').read_text())
    plan = json.loads((args.results / 'PLAN.json').read_text())
    body = {key: value for key, value in plan.items() if key != 'sha256'}
    if hashlib.sha256(canonical(body).encode()).hexdigest() != plan['sha256']:
        raise ValueError('saved plan digest differs')
    if results['plan_sha256'] != plan['sha256']: raise ValueError('results belong to another plan')
    if set(results['systems']) != {job['system'] for job in plan['jobs']}:
        raise ValueError('saved systems differ from the analysis plan')
    for receipt in [*results['systems'].values(), results['qa']]:
        if receipt['plan_sha256'] != plan['sha256']: raise ValueError('saved receipt belongs to another plan')
        for name, pin in receipt['files'].items():
            if sha(args.results / name) != pin: raise ValueError('saved rendering input changed: ' + name)
    if args.out.exists() and any(args.out.iterdir()): raise ValueError('render needs a new or empty --out')
    renderer.main(['--packages', str(args.results / 'packages'), '--output', str(args.out)])
    tables.main(['--results', str(args.results / 'results'), '--out', str(args.out / 'tables')])
    render_qa.main(['--inputs', str(args.results / 'packages/qa_inputs.npz'), '--output', str(args.out)])
    # Detector and model figures require their own recorded histogram/training
    # inputs. Absence stays visible in the figure inventory, never a blank PNG.
    inventory = yaml.safe_load(args.figures.read_text())
    made = {item['path'].split('fig/', 1)[-1] for item in json.loads((args.out / 'FIGURE_MANIFEST.json').read_text())['figures']}
    qa_manifest = json.loads((args.out / 'QA_FIGURE_MANIFEST.json').read_text())
    made.update(item['path'].removeprefix('fig/') for item in qa_manifest['figures'])
    rows = [dict(item, status='written' if item['path'].removeprefix('fig/') in made else 'missing_input')
            for item in inventory['figures']]
    save(args.out / 'FIGURE_AVAILABILITY.json', dict(schema='PhotonJetFigureAvailabilityV1', figures=rows,
                                                   qa_input_notes=qa_manifest['notes']))
    return dict(status='PASS', output=str(args.out), missing_figure_inputs=sum(r['status'] == 'missing_input' for r in rows))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['runPP', 'runAuAu', 'runAll', 'runSystematics', 'render'])
    parser.add_argument('--manifest', type=Path, help='TreeToHists manifest.json')
    parser.add_argument('--config', type=Path, default=HERE / 'config/final_analysis.yaml')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--check', action='store_true', help='validate and print the plan without calculating')
    parser.add_argument('--resume', action='store_true', help='verify and reuse completed results of the identical plan')
    parser.add_argument('--results', type=Path, help='saved analysis directory, for render')
    parser.add_argument('--figures', type=Path, default=HERE / 'config/figures.yaml')
    parser.add_argument('--system', help='system key for runSystematics')
    parser.add_argument('--ratio', nargs=2, help='numerator and denominator for runSystematics')
    args = parser.parse_args(argv)
    if args.command == 'render':
        if not args.results: parser.error('render requires --results')
        if args.check or args.resume: parser.error('render writes a new directory; --check/--resume do not apply')
        result = render(args)
    elif args.command == 'runSystematics':
        if bool(args.system) == bool(args.ratio): parser.error('runSystematics needs exactly one of --system or --ratio')
        if args.check or args.resume: parser.error('--check/--resume apply to measurement commands')
        systematics.main(['--config', str(args.config), '--output', str(args.out)] +
                         (['--system', args.system] if args.system else ['--ratio', *args.ratio]))
        result = dict(status='PASS', output=str(args.out))
    else:
        if not args.manifest: parser.error('measurement commands require --manifest')
        plan = prepare(args)
        result = plan if args.check else execute(args, plan)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as exc:
        print('FinalAnalysis: ' + str(exc), file=sys.stderr)
        raise SystemExit(2)
