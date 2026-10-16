#!/usr/bin/env python3
"""Index native per-output scientific dependencies for the canonical owner.

Reads only receipts and their exact native products. It creates an inventory,
never a selection change, publication pointer or human acceptance. Origin
records are supplied by the producer; missing origins block owner admission.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
from scientific_origins import validate_origin
from native_producer_ancestry import ancestry_errors


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def physics(receipt):
    envelope = receipt['physics']
    text = envelope['canonical']
    require(hashlib.sha256(text.encode()).hexdigest() == envelope['sha256'], 'physics binding differs')
    body = json.loads(text)
    require(canonical(body) == text and body['schema'] == 'PlotPhysicsV1', 'physics format differs')
    require(body['output'] == receipt['plot_recipe']['output'], 'output binding differs')
    return body


COVERAGE = {'quantities','errors','covariance','scalar_annotations','selection',
            'weights','normalization','models','projection','dependencies'}
HASH=re.compile(r'[0-9a-f]{64}')


def qualified(receipt,body,registry):
    """The renderer admission/coverage gates shared with the review consumer.

    Return false for missing evidence. Scientific incompleteness is a recorded
    unknown state, never an exception used to discard the output inventory.
    """
    try:
        family=registry['families'][body['renderer_family']];pins=receipt.get('style_and_renderer',{})
        if (registry.get('schema')!='RendererFamilyCoverageV1' or family.get('state')!='admitted' or
            family.get('contract_version')!=body.get('contract_version') or
            not HASH.fullmatch(str(family.get('fixture_evidence_sha256','')))):return False
        if not any(p and all(pins.get(k)==v for k,v in p.items()) for p in family.get('tested_producer_pinsets',[])):return False
        if family.get('requires_native_ancestry') and ancestry_errors(receipt.get('native_source_witness'),
                receipt.get('native_source_ancestry'),[r.get('source_inventory') for r in family.get('tested_native_releases',[])]):return False
        if family.get('requires_native_source_witness'):
            witness=receipt.get('native_source_witness',{});sources=witness.get('source_inventory')
            if (witness.get('schema')!='NativeProducerSourceWitnessV1' or witness.get('state')!='recorded' or
                not isinstance(sources,dict) or not sources or
                not all(HASH.fullmatch(str(value)) for value in sources.values())):return False
            if not any(release.get('source_inventory')==sources and
                HASH.fullmatch(str(release.get('qualification_report_sha256',''))) and release.get('producing_pinset') and
                all(pins.get(k)==v for k,v in release['producing_pinset'].items())
                for release in family.get('tested_native_releases',[])):return False
        coverage=body['coverage']
        if set(coverage)!=COVERAGE or any(row.get('state') not in ('complete','not_applicable') or
                                         not row.get('evidence') for row in coverage.values()):return False
        if any(coverage[field]['state']!='complete' for field in ('quantities','dependencies')):return False
        return all(isinstance(body.get(field),dict) and body[field] for field in ('definitions','quantities','dependencies'))
    except (KeyError,TypeError,ValueError,AttributeError):return False

def unknown(value):
    if isinstance(value,dict):
        return value.get('state') in ('unknown','missing') or any(unknown(v) for v in value.values())
    return isinstance(value,list) and any(unknown(v) for v in value)


def build(receipts, native_root, registry, source_root=None, qualification=None):
    native_root = Path(native_root).resolve()
    definitions, outputs, reasons = {}, {}, []
    for receipt_path in sorted(map(Path, receipts)):
        receipt = json.loads(receipt_path.read_text())
        body = physics(receipt)
        output = body['output']
        require(output not in outputs, 'duplicate output; select one exact generation per inventory')
        native = (native_root / output).resolve()
        require(native_root in native.parents, 'native output escapes its declared generation')
        require(sha(native) == receipt['image_sha256'], 'native bytes differ: ' + output)
        family = registry.get('families', {}).get(body['renderer_family'], {})
        complete=bool(qualified(receipt,body,registry))
        capture_complete=complete
        if not complete:
            reasons.append(output + ': native scientific completeness is unqualified')
        origins = receipt.get('scientific_definition_origins', {})
        edges = set()
        for series in body['dependencies'].values():
            if not isinstance(series, dict):
                complete = False
                continue
            for root in series.values():
                for term_key, value in root.get('terms', {}).items():
                    # Producer term ids are canonical semantic addresses. The
                    # channel distinguishes physical populations, never runs.
                    category, local = term_key.split('/', 1)
                    identity = {'observable': root['identity']['observable'],
                                'channel': root['identity']['channel'], 'term': local}
                    if root['identity'].get('roles'):identity['roles']=root['identity']['roles']
                    definition_id = root.get('references', {}).get(term_key) or canonical(identity)
                    origin = origins.get(definition_id)
                    # A declared editable reference may not be downgraded to
                    # read-only merely by deleting its origin metadata.
                    readonly=(origin is None and term_key not in root.get('references',{}) and
                              capture_complete and not unknown(value))
                    if origin is None and not readonly:
                        complete = False
                        reasons.append(output + ': no qualified native definition witness for ' + definition_id)
                    origin_binding = []
                    if origin and source_root:
                        exports = {}
                        if qualification:
                            require(qualification.get('schema') == 'NativeRendererFamilyQualificationV1' and
                                    qualification.get('status') == 'PASS' and
                                    qualification.get('producing_pinset') in family.get('tested_producer_pinsets', []),
                                    'origin exports require the exact admitted executable qualification')
                            exports = qualification.get('native_value_exports', {})
                        try:
                            origin_binding = validate_origin(origin, value, source_root, exports)
                        except (ValueError, KeyError, IndexError, TypeError, OSError) as error:
                            complete = False
                            reasons.append(output + ': origin verification failed for ' + definition_id + ': ' + str(error))
                    elif not readonly:
                        complete = False
                        reasons.append(output + ': canonical source readback unavailable for ' + definition_id)
                    row = {'category': category, 'value': value, 'origin': origin,
                           'origin_binding': origin_binding,
                           'state': 'known' if origin_binding else 'native_readonly' if readonly else 'origin_unverified'}
                    if readonly:
                        row['witnesses']=[{'receipt_path':str(receipt_path.resolve()),'receipt_sha256':sha(receipt_path),
                            'physics_sha256':receipt['physics']['sha256'],'definition_id':definition_id,
                            'value_sha256':hashlib.sha256(canonical(value).encode()).hexdigest()}]
                    if definition_id in definitions:
                        previous=definitions[definition_id]
                        require({k:v for k,v in previous.items() if k!='witnesses'} ==
                                {k:v for k,v in row.items() if k!='witnesses'}, 'contradictory definitions in one generation')
                        if readonly:
                            row['witnesses']=sorted({canonical(w):w for w in previous['witnesses']+row['witnesses']}.values(),key=canonical)
                    definitions[definition_id] = row
                    edges.add(definition_id)
        if not edges:
            complete = False
            reasons.append(output + ': no typed dependency edges')
        outputs[output] = {'native_path': str(native), 'receipt_path': str(receipt_path.resolve()),
                           'renderer_family': body['renderer_family'],
                           'native_sha256': receipt['image_sha256'],
                           'receipt_sha256': sha(receipt_path),
                           'physics_sha256': receipt['physics']['sha256'],
                           'complete': complete, 'dependencies': sorted(edges)}
    reverse = {key: sorted(output for output, row in outputs.items() if key in row['dependencies'])
               for key in definitions}
    return {'schema': 'PhysicsNativeInventoryV1', 'definitions': definitions, 'outputs': outputs,
            'affected_outputs': reverse, 'owner_admissible': bool(outputs) and all(o['complete'] for o in outputs.values()),
            'unavailable_reasons': sorted(set(reasons))}


def impacted(before, changed_definitions):
    require(before.get('schema') == 'PhysicsNativeInventoryV1', 'wrong native inventory')
    require(before.get('owner_admissible') is True, 'inventory is incomplete; propagation is unavailable')
    require(set(changed_definitions) <= set(before['definitions']), 'undeclared changed definition')
    return sorted({output for key in changed_definitions for output in before['affected_outputs'][key]})


def proposal_roster(before, after, changed_definitions, declared_outputs):
    affected = impacted(before, changed_definitions)
    require(set(affected) == set(declared_outputs), 'declared change group differs from the dependency impact roster')
    require(after.get('owner_admissible') is True, 'candidate inventory is incomplete')
    require(set(affected) <= set(after['outputs']), 'affected candidate outputs have not all been produced')
    changed = {key for key in set(before['definitions']) | set(after['definitions'])
               if key not in before['definitions'] or key not in after['definitions'] or
               before['definitions'][key]['value'] != after['definitions'][key]['value']}
    require(changed <= set(changed_definitions), 'candidate changed an undeclared canonical definition')
    return {'schema': 'PhysicsChangeGroupRosterV1', 'affected_outputs': affected,
            'definition_ids': sorted(changed_definitions),
            'source_adoption_requires_new_qualification': sorted(key for key in changed_definitions
                if before['definitions'][key].get('state')=='native_readonly'),
            'baseline': {key: before['outputs'][key] for key in affected},
            'candidate': {key: after['outputs'][key] for key in affected},
            'authority': 'EXPLICIT_GROUP_APPROVAL_REQUIRED', 'canonical_effect': 'NOT_EXECUTED'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-root', type=Path, required=True)
    p.add_argument('--registry', type=Path, required=True)
    p.add_argument('--receipt', type=Path, action='append', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--qualification', type=Path)
    args = p.parse_args()
    result = build(args.receipt, args.native_root, json.loads(args.registry.read_text()), args.source_root,
                   json.loads(args.qualification.read_text()) if args.qualification else None)
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
