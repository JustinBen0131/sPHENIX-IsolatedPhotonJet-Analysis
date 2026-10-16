#!/usr/bin/env python3
"""Read-only validation of canonical scientific definition edit locations.

The location digest is provenance. The effective value is already recorded by
native instrumentation, and is never replaced by a digest of source text.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import re
import struct


def require(ok, message):
    if not ok:
        raise ValueError(message)


def portable(value):
    if isinstance(value, float):
        require(math.isfinite(value), 'nonfinite effective value')
        return 'f64:' + struct.pack('>d', value if value else 0.).hex()
    if isinstance(value, list):
        return [portable(v) for v in value]
    if isinstance(value, dict):
        return {k: portable(v) for k, v in value.items()}
    return value


def at_pointer(value, pointer):
    require(pointer == '' or pointer.startswith('/'), 'invalid JSON pointer')
    for token in pointer.split('/')[1:]:
        key = token.replace('~1', '/').replace('~0', '~')
        value = value[int(key)] if isinstance(value, list) else value[key]
    return value


def source_path(root, relative):
    root = Path(root).resolve()
    raw = Path(relative)
    require(not raw.is_absolute() and '..' not in raw.parts, 'origin must be repository relative')
    path = (root / raw).resolve()
    require(root in path.parents and path.is_file(), 'origin escapes source root or is missing')
    return path


def decoded_literal(text, encoding):
    require(text == text.strip(), 'literal cannot include surrounding whitespace')
    if encoding in ('binary64', 'binary32'):
        require(re.fullmatch(r'[+-]?(?:(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?|0[xX](?:[0-9a-fA-F]+(?:\.[0-9a-fA-F]*)?|\.[0-9a-fA-F]+)[pP][+-]?[0-9]+)[fFlL]?', text), 'numeric origin is not a C++ literal')
        raw = text[:-1] if text[-1:] in ('f','F','l','L') else text
        require(encoding != 'binary64' or text[-1:] not in ('f','F','l','L'), 'binary64 literal has a different declared precision')
        require(encoding != 'binary32' or text[-1:] in ('f','F'), 'binary32 origin requires its float suffix')
        value = float.fromhex(raw) if '0x' in raw.lower() else float(raw)
        if encoding == 'binary32':
            value = struct.unpack('>f', struct.pack('>f', value))[0]
        return portable(value)
    if encoding == 'integer':
        require(re.fullmatch(r'-?(?:0|[1-9][0-9]*)', text), 'integer origin is not a decimal literal')
        return int(text)
    if encoding == 'bool':
        require(text in ('true','false'), 'bool origin must be true or false')
        return text == 'true'
    if encoding == 'string':
        require(text.startswith('"') and text.endswith('"'), 'string origin is not one quoted string')
        value = json.loads(text)
        require(isinstance(value, str), 'string origin must decode to a string')
        return value
    raise ValueError('unavailable literal encoding: ' + str(encoding))


def locate_source(text, site):
    selector = site['selector']
    before, after, target = selector['before_context'], selector['after_context'], site['before_text']
    require(before and after and target, 'source selector requires bounded context and a nonempty span')
    exact = before + target + after
    require(text.count(exact) == 1, 'source selector is absent or ambiguous')
    start = text.index(exact) + len(before)
    return start, start + len(target)


def validate_origin(origin, effective_value, source_root, native_exports=None):
    """Verify all sites; returns bound readback provenance and performs no edit.

    Composite definitions can use value_pointer at each site. native_value
    expressions require an exact native export previously produced by the
    executable qualification, bound to source bytes and span. This validator
    cannot manufacture or infer that export from the expression.
    """
    require(isinstance(origin, dict), 'canonical origin missing')
    sites = origin.get('sites') if origin.get('kind') == 'multiple' else [origin]
    require(isinstance(sites, list) and sites, 'canonical origin has no sites')
    evidence = []
    for site in sites:
        path = source_path(source_root, site['path'])
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        expected = at_pointer(effective_value, site.get('value_pointer', ''))
        if site['kind'] == 'json_pointer':
            actual = at_pointer(json.loads(data), site['pointer'])
            require(actual == site['before'], 'JSON origin before value differs')
            require(portable(actual) == expected, 'JSON origin is not the recorded effective value')
            evidence.append({'kind':'json_pointer','path':site['path'],'pointer':site['pointer'], 'source_sha256':digest,'effective_value':expected})
        elif site['kind'] == 'source_text':
            text = data.decode('utf-8'); start, end = locate_source(text, site)
            encoding = site['encoding']
            if encoding == 'native_value':
                key = site.get('export_id')
                export = (native_exports or {}).get(key, {})
                require(export.get('source_sha256') == digest and export.get('path') == site['path'] and
                        export.get('before_text') == site['before_text'] and export.get('selector') == site['selector'] and
                        export.get('canonical_key') == site.get('canonical_key'),
                        'native expression lacks the exact qualified source/span export')
                require(export.get('value') == expected, 'native expression export differs from the effective value')
            else:
                require(decoded_literal(site['before_text'], encoding) == expected,
                        'source literal differs from the recorded effective value')
            evidence.append({'kind':'source_text','path':site['path'],'source_sha256':digest,'selector':site['selector'],
                             'before_text':site['before_text'],'start_utf8':len(text[:start].encode()),
                             'end_utf8':len(text[:end].encode()),'encoding':encoding,'effective_value':expected})
        else:
            raise ValueError('unsupported canonical origin kind: ' + str(site.get('kind')))
    return evidence
