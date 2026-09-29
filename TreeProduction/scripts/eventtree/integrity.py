"""File identity, metadata and exact payload comparison for eventTree products.

These primitives operate on the public file contract. They do not reconstruct
DSTs, convert historical tables or resolve a model's physics definition.
"""
import hashlib
import json

import awkward as ak
import numpy as np


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def metadata(file):
    """Read direct JSON metadata; nested carried objects keep their own types."""
    return {key[5:]: json.loads(str(file[key]))
            for key, kind in file.classnames(recursive=True, cycle=False).items()
            if key.startswith('meta/') and '/' not in key[5:] and kind == 'TObjString'}


def equal_arrays(a, b):
    """Compare types and packed payload bytes, including NaNs and signed zero."""
    if str(ak.type(a)) != str(ak.type(b)):
        return False
    fa, na, ba = ak.to_buffers(ak.to_packed(a))
    fb, nb, bb = ak.to_buffers(ak.to_packed(b))
    return (fa.to_json() == fb.to_json() and na == nb and ba.keys() == bb.keys()
            and all(ba[key].dtype == bb[key].dtype and
                    ba[key].tobytes() == bb[key].tobytes() for key in ba))


def copy_objects(source, destination):
    """Carry non-tree objects verbatim; direct JSON metadata is written by its owner."""
    destination.copy_from(source,
        filter_name=lambda n: not (n.startswith('meta/') and '/' not in n[5:]),
        filter_classname=lambda cls: cls not in ('TTree', 'TDirectory', 'TDirectoryFile'),
        require_matches=False)


def prove_objects(source, destination):
    """Verify preserved completion, trigger and other carried ROOT objects."""
    for name, kind in source.classnames(recursive=True, cycle=False).items():
        if kind in ('TTree', 'TDirectory', 'TDirectoryFile') or (name.startswith('meta/') and '/' not in name[5:]):
            continue
        a, _ = source.key(name).get_uncompressed_chunk_cursor()
        b, _ = destination.key(name).get_uncompressed_chunk_cursor()
        if not np.array_equal(a.raw_data, b.raw_data):
            raise ValueError('changed carried object: ' + name)
