#!/usr/bin/env python3
"""Run native fixture gates and pin only the renderer families they qualify.

A source edit changes qualification provenance, not scientific identity. Rerun
this command after such an edit; no family is admitted merely because the
generic serializer tests pass. The registry is replaced atomically only after
the complete fixture command succeeds against unchanged source bytes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import uuid
from plot_physics_inventory import sha, require, canonical
from plot_physics_sources import source_inventory, source_roots, families

PINS = {'IANFigures.C': 'plots/IANFigures.C', 'QAFigures.C': 'plots/QAFigures.C',
        'sPhenixStyle.C': 'plots/sPhenixStyle.C', 'sPhenixStyle.h': 'plots/sPhenixStyle.h',
        'RenderMetadata.h': 'src/RenderMetadata.h', 'PlotPhysics.h': 'src/PlotPhysics.h',
        'RootScientificTrace.h': 'src/RootScientificTrace.h',
        'RendererScientificRoster.h': 'src/RendererScientificRoster.h',
        'RenderScientificParameters.h': 'src/RenderScientificParameters.h',
        'ScientificTrace.h': '../TreeToHists/src/ScientificTrace.h',
        'HistogramScientificDefinitions.h': '../TreeToHists/src/HistogramScientificDefinitions.h'}
MUTATIONS = {'values', 'errors', 'covariance', 'scalars', 'selection', 'weights', 'models',
             'normalization', 'projection', 'series_swap', 'missing_measurement',
             'same_pixels_different_calculation', 'cosmetics', 'unrelated_input', 'unrelated_configuration',
             'per_output_reset', 'unbound_primitive', 'missing_series', 'upstream_definition'}


def pinset(root):
    return {name: sha(Path(root) / relative) for name, relative in PINS.items()}


def registry_from_report(report, pins, roster, previous=None, report_reference=None):
    require(report.get('schema') == 'NativeRendererFamilyQualificationV1' and report.get('status') == 'PASS',
            'native family fixture command did not pass')
    require(report.get('producing_pinset') == pins, 'fixture report does not bind exact tested source')
    require(set(report.get('families', {})) == set(roster), 'fixture family roster differs from current renderer roster')
    require(bool(report.get('source_inventory')), 'native qualification has no full producer source inventory')
    result = {'schema': 'RendererFamilyCoverageV1', 'families': {}}
    evidence = hashlib.sha256((canonical(report)+'\n').encode()).hexdigest()
    if report_reference:require(report_reference['sha256']==evidence,'retained qualification report digest differs')
    for family, row in report['families'].items():
        tested = set(row.get('passed_mutations', []))
        passed = (row.get('entrypoint_executed') is True and row.get('no_unbound_primitives') is True and
                  row.get('producer_trace_guards_passed') is True and not row.get('missing_series') and
                  tested >= MUTATIONS and row.get('status') == 'PASS')
        release={'producing_pinset':pins,'source_inventory':report['source_inventory'],
                 'qualification_report_sha256':evidence}
        if report_reference:release['qualification_report']=report_reference
        old=(previous or {}).get('families',{}).get(family,{})
        # Only independently admitted witness-backed releases survive. A
        # legacy/partial row is never upgraded by merging its old pinsets.
        retained=(old.get('tested_native_releases',[]) if old.get('state')=='admitted' and
                  old.get('contract_version')==1 and old.get('requires_native_source_witness') is True else [])
        releases=sorted({canonical(r):r for r in retained+[release]}.values(),key=canonical) if passed else []
        result['families'][family] = {
            'state': 'admitted' if passed else 'partial', 'contract_version': 1,
            'fixture_evidence_sha256': evidence if passed else None,
            'tested_producer_pinsets': sorted({canonical(r['producing_pinset']):r['producing_pinset']
                for r in releases}.values(),key=canonical),
            'requires_native_source_witness':True, 'requires_native_ancestry':True, 'tested_native_releases':releases,
            'missing_mutations': sorted(MUTATIONS - tested),
            'why': 'Native entrypoint, membership, numerical and producer-trace guard fixtures passed; each output still requires its own complete ancestry' if passed else
                   'Native family entrypoint/roster/trace guard fixture incomplete'}
    require(result['families'] and all(row['state']=='admitted' for row in result['families'].values()),
            'one or more native renderer families remain unqualified')
    return result


def retained_releases(previous):
    """Refuse to silently forget historical admission whose evidence was lost."""
    for row in previous.get('families',{}).values():
        if row.get('state')!='admitted' or not row.get('requires_native_source_witness'):continue
        for release in row.get('tested_native_releases',[]):
            reference=release.get('qualification_report',{})
            require(reference.get('path') and reference.get('sha256')==release.get('qualification_report_sha256'),
                    'historical native release lacks a retained qualification report')
            path=Path(reference['path'])
            require(path.is_file() and sha(path)==reference['sha256'],'historical native qualification evidence changed or is missing')
            report=json.loads(path.read_text())
            require(report.get('status')=='PASS' and report.get('producing_pinset')==release.get('producing_pinset') and
                    report.get('source_inventory')==release.get('source_inventory'),'historical native release qualification differs')


def retain_report(directory,report):
    data=(canonical(report)+'\n').encode();digest=hashlib.sha256(data).hexdigest()
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);path=directory/(digest+'.json')
    try:
        with path.open('xb') as stream:stream.write(data);stream.flush();os.fsync(stream.fileno())
    except FileExistsError:require(path.read_bytes()==data,'immutable native qualification collision')
    return {'path':str(path.resolve()),'sha256':digest}


def write_json_atomic(path, value):
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as temporary:
        temporary.write(json.dumps(value, indent=2) + '\n')
        temporary.flush()
        os.fsync(temporary.fileno())
        name = temporary.name
    os.replace(name, path)


def qualify(root, registry_path, *, report_path=None, root_config=None, compiler=None, evidence_dir=None):
    root, registry_path = Path(root).resolve(), Path(registry_path).resolve()
    report_path=Path(report_path).resolve() if report_path else root/'plot_physics_qualification.json'
    require(report_path is None or Path(report_path).resolve() != registry_path,
            'qualification report must not replace the admission registry')
    before = pinset(root)
    sources_before = source_inventory(root.parent)
    roster = json.loads((root / 'plots/physics_family_roster.json').read_text())['families']
    require(set(roster)==set(families(root)), 'declared family roster differs from maintained save entrypoints')
    # Fixed package-owned entry point. No commands supplied by receipts.
    evidence_dir=Path(evidence_dir).resolve() if evidence_dir else root/'.physics'/'fixture-runs'/str(uuid.uuid4())
    require(not evidence_dir.exists(),'qualification evidence destination must be a fresh directory')
    command = [sys.executable, str(root / 'tests/run_plot_physics_fixtures.py'), '--json','--work',str(evidence_dir)]
    if root_config is not None: command += ['--root-config', str(Path(root_config).resolve())]
    if compiler is not None: command += ['--compiler', str(Path(compiler).resolve())]
    env = dict(os.environ)
    env.pop('FINALANALYSIS_CACHE', None)
    completed = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True, check=False)
    evidence_dir.mkdir(parents=True,exist_ok=True)
    (evidence_dir/'qualification.command.json').write_text(json.dumps(command,indent=2)+'\n')
    (evidence_dir/'qualification.stdout.json').write_text(completed.stdout)
    (evidence_dir/'qualification.stderr.log').write_text(completed.stderr)
    try:
        report = json.loads(completed.stdout)
    except ValueError as error:
        raise ValueError('native fixtures returned no valid report; registry unchanged: ' + completed.stderr[-1500:]) from error
    # Preserve a PARTIAL/FAIL report for integration without weakening admission.
    if report_path is not None: write_json_atomic(report_path, report)
    require(completed.returncode == 0,
            'native qualification ' + str(report.get('status', 'FAIL')) + '; registry unchanged')
    after = pinset(root)
    require(before == after, 'source changed during qualification; registry unchanged')
    require(sources_before == source_inventory(root.parent), 'producer source inventory changed during qualification')
    require(report.get('source_inventory') == sources_before and report.get('source_roots') == source_roots() and
            report.get('source_root')==str(root.parent),
            'qualification does not bind the complete producer source inventory')
    previous=json.loads(registry_path.read_text()) if registry_path.is_file() else {}
    retained_releases(previous)
    reference=retain_report(registry_path.parent/'.physics'/'qualifications',report)
    result = registry_from_report(report, after, roster,previous,reference)
    result['tested_source_inventory'] = sources_before
    result['source_root']=str(root.parent)
    result['source_roots'] = source_roots()
    result['qualification_report_sha256'] = reference['sha256']
    result['qualification_report']=reference
    write_json_atomic(registry_path, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--registry', type=Path)
    p.add_argument('--report', type=Path, help='retain PASS/PARTIAL/FAIL native evidence separately from the registry')
    p.add_argument('--root-config', type=Path)
    p.add_argument('--compiler', type=Path)
    p.add_argument('--evidence-dir',type=Path,help='fresh retained fixture directory; defaults to .physics/fixture-runs/<uuid>')
    a = p.parse_args()
    try:
        result = qualify(a.root, a.registry or a.root / 'plot_physics_coverage.json',
                         report_path=a.report, root_config=a.root_config, compiler=a.compiler,evidence_dir=a.evidence_dir)
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps({'admitted': sum(row['state'] == 'admitted' for row in result['families'].values()),
                      'families': len(result['families'])}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
