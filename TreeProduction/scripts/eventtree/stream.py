"""Bounded-memory TTree writer used by consolidation and augmentation.

The 655 named payload branches are unchanged. Independent generated counters
avoid relying on a group that happened to be empty in the first input chunk.
"""
from pathlib import Path

import awkward as ak
import numpy as np

from . import rootio
from .schema import Schema
from .tables import offsets, require


class Writer:
    def __init__(self, path, schema=None, compression=404):
        self.schema = schema or Schema()
        self.path = Path(path)
        require(not self.path.exists(), 'output exists: '+str(path))
        self.root = rootio.root()
        self.writer = self.root.eventtree.Writer(str(path),'eventTree','RJEventTreeV1',compression,0)
        self.entries = 0
        self.counters = {n:'n_'+n for n,b in self.schema.branches.items()
                         if b.kind == 'vector' and b.dtype is not object}
        # Maintained ROOT response readers consume these public collection
        # counters directly. Anchor them to pT, including in empty batches.
        for prefix in ('jet','jetunsub','truthjet'):
            for radius in ('02','03','04'):
                self.counters[prefix+radius+'_pt']='n_'+prefix+radius

    def append(self, arrays):
        n = len(arrays)
        require(set(self.schema.branches) <= set(arrays.fields),'stream payload branches missing')
        self.writer.beginBatch(n)
        keep = []
        for name, b in self.schema.branches.items():
            a = arrays[name]
            if b.kind == 'vector':
                counts = ak.to_numpy(ak.num(a,axis=1)).astype(np.int64)
                off = offsets(counts).view(np.longlong); keep.append(off)
                if b.dtype is object:
                    data = self.root.std.vector['std::string']()
                    for v in ak.to_list(ak.flatten(a)): data.push_back(v)
                    starts = self.root.std.vector['long long']()
                    for v in off: starts.push_back(int(v))
                    self.writer.strings(name,data,starts,rootio.title_of(b)); continue
                self.writer.counter(self.counters[name],off)
                arr = rootio.buffer(ak.to_numpy(ak.flatten(a)),b.dtype); keep.append(arr)
                cxx,leaf,_ = rootio.CXX[np.dtype(b.dtype)]
                self.writer.counted[cxx](name,leaf,arr,off,self.counters[name],int(counts.max()) if n else 0,rootio.title_of(b))
            else:
                arr = rootio.buffer(ak.to_numpy(a).reshape(-1),b.dtype); keep.append(arr)
                cxx,leaf,_ = rootio.CXX[np.dtype(b.dtype)]
                if b.kind == 'scalar': self.writer.scalar[cxx](name,leaf,arr,rootio.title_of(b))
                else: self.writer.fixed[cxx](name,leaf,arr,b.length,rootio.title_of(b))
        require(self.writer.fill() == n,'stream write entry mismatch')
        self.entries += n

    def close(self):
        self.writer.close()


class Output:
    """Branch buffers over the entry axis (retained and upstream-rejected events in source order)."""

    def __init__(self, schema, retained_positions, n_entries):
        self.schema, self.pos, self.n = schema, retained_positions, n_entries
        self.values = {}

    def scalar(self, name, values):
        self._check(name, 'scalar', values)
        out = np.full(self.n, sentinel(self.schema.branches[name].dtype), dtype=self.schema.branches[name].dtype)
        out[self.pos] = values
        self.values[name] = ('scalar', out, None)

    def fixed(self, name, values):
        branch = self.schema.branches[name]
        self._check(name, 'fixed', values.reshape(-1))
        require(values.shape[1:] == (branch.length,), name + ': fixed length differs')
        out = np.zeros((self.n, branch.length), dtype=branch.dtype)
        out[self.pos] = values
        self.values[name] = ('fixed', out, None)

    def vector(self, name, values, counts_per_event):
        branch = self.schema.branches[name]
        if branch.dtype is not object:
            self._check(name, 'vector', values)
        require(len(values) == int(np.sum(counts_per_event)), name + ': value count differs from the group size')
        counts = np.zeros(self.n, dtype=np.int64)
        counts[self.pos] = counts_per_event
        self.values[name] = ('vector', values, counts)

    def _check(self, name, kind, values):
        branch = self.schema.branches.get(name)
        require(branch is not None, 'target branch missing from the schema: ' + name)
        require(branch.kind == kind, '%s: schema kind %s, conversion kind %s' % (name, branch.kind, kind))
        require(np.asarray(values).dtype == np.dtype(branch.dtype),
                '%s: source type %s differs from schema type %s' % (name, np.asarray(values).dtype, branch.dtype))


def sentinel(dtype):
    """Value of a scalar branch on an entry or product where its source does not apply."""
    return np.nan if np.dtype(dtype).kind == 'f' else 0

