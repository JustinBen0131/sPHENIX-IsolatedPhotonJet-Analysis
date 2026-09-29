"""Expand the maintained eventTree contract into concrete branches and types.

Production, augmentation and readers use config/eventTree_v1.yaml. Historical
conversion mappings belong to private regression fixtures, outside runtime.
File digests identify exact revisions; contract digests permit only explicitly
reviewed revisions with identical branch and output definitions.
"""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

SCHEMA_PATH = Path(__file__).resolve().parents[2] / 'config' / 'eventTree_v1.yaml'
RADII = ('02', '03', '04')
TEMPLATES = ('jet<RR>', 'jetunsub<RR>', 'truthjet<RR>')

# The earlier v1 file included a conversion-only sources table. Its removal
# leaves this exact output contract unchanged; arbitrary old digests refuse.
HISTORICAL_CONTRACTS = {
    '0ce3b08f37826282644ab6e3929582e920f53918288ad573801d2bb4cb8d4a44':
        '7129b8649b9273df2ebba079fd40a7f99a69912554bc5189881982be299f83dc',
}

DTYPES = {'int8': np.int8, 'uint8': np.uint8, 'int16': np.int16, 'uint16': np.uint16,
          'int32': np.int32, 'uint32': np.uint32, 'int64': np.int64, 'uint64': np.uint64,
          'float': np.float32, 'double': np.float64, 'string': object}


@dataclass(frozen=True)
class Branch:
    name: str
    group: str
    kind: str            # 'scalar', 'vector' or 'fixed'
    dtype: object        # numpy dtype, or object for strings
    length: int          # fixed-array length; 0 otherwise
    source: str
    title: str
    unit: str
    reserved: bool = False

    @property
    def leaf(self):
        return self.dtype is not object and np.dtype(self.dtype)


def parse_type(text):
    """'vector<int32>' -> ('vector', int32, 0); 'uint64[64]' -> ('fixed', uint64, 64); 'double' -> scalar."""
    if text.startswith('vector<') and text.endswith('>'):
        return 'vector', DTYPES[text[7:-1]], 0
    if text.endswith(']'):
        base, _, size = text[:-1].partition('[')
        return 'fixed', DTYPES[base], int(size)
    return 'scalar', DTYPES[text], 0


def expand(name):
    """Concrete names of a possibly templated branch or group name."""
    for template in TEMPLATES:
        if template in name:
            return [name.replace('<RR>', radius) for radius in RADII]
    return [name]


class Schema:
    def __init__(self, path=SCHEMA_PATH):
        self.path = Path(path)
        raw = self.path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.doc = yaml.safe_load(raw)
        contract = {key: value for key, value in self.doc.items() if key != 'sources'}
        self.contract_sha256 = hashlib.sha256(
            json.dumps(contract, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        self.version = self.doc['version']
        self.branches = {}
        for item in self.doc['branches']:
            self._add(item, reserved=False)
        for item in self.doc['reserved']:
            self._add(dict(item, unit='', source='reserved:' + item['filled_by']), reserved=True)

    def accepts_digest(self, digest):
        """Accept this revision or a pinned revision of the identical contract."""
        return digest == self.sha256 or HISTORICAL_CONTRACTS.get(digest) == self.contract_sha256

    def _add(self, item, reserved):
        kind, dtype, length = parse_type(item['type'])
        for name in expand(item['name']):
            group = name.split('_', 1)[0] if '<RR>' in item['group'] else item['group']
            if name in self.branches:
                raise ValueError('duplicate branch in schema: ' + name)
            self.branches[name] = Branch(name, group, kind, dtype, length, item['source'],
                                         item['title'], item.get('unit', ''), reserved)

    def group(self, group):
        return [b for b in self.branches.values() if b.group == group]
