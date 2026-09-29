"""Add declared branches to an immutable eventTree candidate.

The caller owns the calculation and its scientific validation. This interface
owns bounded writing, metadata/object preservation, saved-file comparison and
input identity checks. Publication remains a separate whole-group operation.
An existing candidate is checked, never overwritten or silently regenerated.
"""
from pathlib import Path

import numpy as np
import uproot

from .integrity import equal_arrays, metadata, sha256, copy_objects, prove_objects
from .stream import Writer
from .tables import dumps, require


def prove_preserved(original, candidate, *, added, expected_metadata):
    """Reopen and prove exact old payloads, ordering, types and carried objects."""
    entries = 0
    with uproot.open(original) as before, uproot.open(candidate) as after:
        require(metadata(after) == expected_metadata, 'updated metadata differs from the declared change')
        old, new = before['eventTree'], after['eventTree']
        require(old.num_entries == new.num_entries, 'update changed event count')
        require(set(old.keys()) <= set(new.keys()), 'update removed an existing branch')
        require(set(new.keys()) - set(old.keys()) <= set(added) | {'n_' + n for n in added},
                'update added an undeclared branch')
        for name in old.keys():
            require(old[name].typename == new[name].typename, 'update changed branch type: ' + name)
        # Counters may change names when a legacy shared counter is expanded;
        # all retained counters and every payload still have exact values.
        for arrays in old.iterate(step_size='16 MB', library='ak'):
            saved = new.arrays(list(arrays.fields), entry_start=entries,
                               entry_stop=entries + len(arrays), library='ak')
            for name in arrays.fields:
                require(equal_arrays(arrays[name], saved[name]), 'update changed protected branch: ' + name)
            entries += len(arrays)
        prove_objects(before, after)
        require(entries == old.num_entries, 'update proof did not cover every event')
        protected = len(old.keys())
    return dict(status='PASS', entries=entries, protected_branches=protected)


def add_branches(original, candidate, *, input_sha256, schema, added, chunks,
                 transform, new_metadata, progress=None):
    """Write one candidate and verify conservation; ``chunks`` is restartable.

    Each chunk is ``(arrays, context)``. ``transform`` returns all old payloads
    plus the declared additions. It is called once on typed empty input too.
    ``new_metadata`` is the complete intended metadata, copied by the caller.
    No published file is modified, and no output is declared scientifically
    valid merely because the generic conservation check succeeds.
    """
    original, candidate = Path(original), Path(candidate)
    require(original.resolve() != candidate.resolve(), 'update must use a separate candidate')
    require(sha256(original) == input_sha256, 'update input hash changed')
    added = set(added)
    with uproot.open(original) as source:
        payload = {n for n in source['eventTree'].keys() if not n.startswith('n_')}
        require(payload.isdisjoint(added), 'addition would overwrite an existing branch')
        require(payload | added == set(schema.branches), 'update schema has missing or undeclared payloads')
        # Derive counters before writing so restart validation has the same
        # exact expected metadata as the initial attempt.
        expected_counters = {n: 'n_' + n for n, b in schema.branches.items()
                             if b.kind == 'vector' and b.dtype is not object}
        for prefix in ('jet', 'jetunsub', 'truthjet'):
            for radius in ('02', '03', '04'):
                expected_counters[prefix + radius + '_pt'] = 'n_' + prefix + radius
        new_metadata['schema']['counters'] = expected_counters
        if not candidate.exists():
            candidate.parent.mkdir(parents=True, exist_ok=True)
            writer = Writer(candidate, schema)
            entries = 0
            try:
                for arrays, context in chunks():
                    output = transform(arrays, context)
                    require(len(output) == len(arrays), 'update calculation changed event count')
                    writer.append(output)
                    entries += len(arrays)
                    if progress:
                        progress(entries, source['eventTree'].num_entries)
                if not entries:
                    empty = source['eventTree'].arrays(entry_start=0, entry_stop=0, library='ak')
                    writer.append(transform(empty, metadata(source)))
            finally:
                writer.close()
            require(entries == source['eventTree'].num_entries, 'update omitted input events')
            require(writer.counters == expected_counters, 'writer counter contract changed')
            with uproot.update(candidate) as out:
                copy_objects(source, out)
                for name, value in new_metadata.items():
                    out['meta/' + name] = dumps(value)
    result = prove_preserved(original, candidate, added=added, expected_metadata=new_metadata)
    require(sha256(original) == input_sha256, 'update input changed during execution')
    return result
