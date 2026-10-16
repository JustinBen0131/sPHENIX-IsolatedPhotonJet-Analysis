#!/usr/bin/env python3
"""Bind a qualified numerical capture to an unchanged legacy native artifact.

The exact old renderer and exact retained inputs must have been replayed into
an isolated directory. Byte identity with the accepted native artifact is a
necessary binding check, never a physics-equality classifier. The capture must
also carry qualified complete definitions. This command runs no renderer,
changes no acceptance, and overwrites no old image or receipt.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from plot_physics_inventory import sha, physics, require, COVERAGE, qualified
from audit_plot_capture import audit as audit_capture


def create_binding(accepted_native, accepted_receipt, replay_manifest, registry):
    accepted_native, accepted_receipt = Path(accepted_native), Path(accepted_receipt)
    old = json.loads(accepted_receipt.read_text())
    replay = json.loads(Path(replay_manifest).read_text())
    require(replay.get('schema') == 'LegacyPlotReplayCaptureV1', 'wrong replay schema')
    native_sha = sha(accepted_native)
    require(native_sha == old['image_sha256'], 'accepted native bytes do not bind the old receipt')
    require(replay.get('baseline_native_sha256') == native_sha and
            replay.get('baseline_receipt_sha256') == sha(accepted_receipt), 'replay baseline differs')
    require(replay['output'] == old['plot_recipe']['output'], 'replay output differs')
    # Every original renderer/input pin must be verified, not merely named.
    for field, original in [('renderer_files', old.get('style_and_renderer', {})), ('inputs', old.get('inputs', {}))]:
        require(original and set(replay[field]) == set(original), 'replay ' + field + ' roster differs')
        for name, digest in original.items():
            row = replay[field][name]
            require(row['sha256'] == digest and sha(row['retained_path']) == digest, 'retained ' + field + ' differs: ' + name)
    require(sha(replay['replayed_native_path']) == native_sha,
            'replay is not byte-identical; no baseline binding is certified')
    capture_path = Path(replay['capture_receipt_path'])
    capture = json.loads(capture_path.read_text())
    body = physics(capture)
    require(body['output'] == replay['output'] and capture['image_sha256'] == native_sha, 'capture is not bound to the replayed native bytes')
    family = registry.get('families', {}).get(body['renderer_family'], {})
    require(family.get('state') == 'admitted' and family.get('contract_version') == body['contract_version'], 'capture family not admitted')
    require(qualified(capture,body,registry),'capture lacks the exact complete admitted release evidence')
    require(any(p and all(capture.get('style_and_renderer', {}).get(k) == v for k, v in p.items())
                for p in family.get('tested_producer_pinsets', [])), 'capture recorder is not qualified')
    require(set(body['coverage']) == COVERAGE and all(v.get('state') in ('complete', 'not_applicable') for v in body['coverage'].values()), 'capture remains scientifically incomplete')
    require(body['quantities'] and body['definitions'] and body['dependencies'], 'capture omits scientific content')
    require(replay.get('instrumentation_audit', {}).get('state') == 'PASS', 'instrument-only replay audit absent')
    audit = replay['instrumentation_audit']
    require(sha(audit['path']) == audit['sha256'], 'instrumentation audit changed')
    audit_body = json.loads(Path(audit['path']).read_text())
    audit_inputs = replay.get('instrumentation_sources', {})
    rebuilt = audit_capture(audit_inputs['original_files'], audit_inputs['capture_files'],
                            old['style_and_renderer'], capture['style_and_renderer'])
    require(rebuilt == audit_body, 'instrumentation audit is not reproducible from the retained exact sources')
    require(audit_body.get('schema') == 'PlotCaptureInstrumentationAuditV1' and
            audit_body.get('scientific_expressions_unchanged') is True and
            audit_body.get('original_renderer') == old['style_and_renderer'] and
            audit_body.get('capture_renderer') == capture['style_and_renderer'], 'instrumentation audit does not bind both exact renderers')
    require(sha(accepted_native) == native_sha, 'accepted native changed during the binding')
    return {'schema': 'LegacyPlotPhysicsBindingV1', 'output': replay['output'],
            'baseline_native_sha256': native_sha, 'baseline_receipt_sha256': sha(accepted_receipt),
            'capture_receipt_sha256': sha(capture_path), 'capture_physics': capture['physics'],
            'capture_renderer': capture['style_and_renderer'], 'replay_manifest_sha256': sha(replay_manifest),
            'instrumentation_audit_sha256': audit['sha256'], 'acceptance_effect': 'NONE', 'placement_effect': 'NONE'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--accepted-native', type=Path, required=True)
    p.add_argument('--accepted-receipt', type=Path, required=True)
    p.add_argument('--replay-manifest', type=Path, required=True)
    p.add_argument('--registry', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    result = create_binding(a.accepted_native, a.accepted_receipt, a.replay_manifest, json.loads(a.registry.read_text()))
    with a.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
