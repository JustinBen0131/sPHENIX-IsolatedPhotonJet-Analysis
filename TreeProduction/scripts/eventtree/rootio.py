"""Write an eventTree file: the TTree through the C++ writer core, meta objects and carried objects with uproot."""
from pathlib import Path

import numpy as np
import uproot

from .tables import dumps, offsets, require

WRITER = Path(__file__).with_name('writer.hxx')
CXX = {np.dtype(np.int8): ('Char_t', 'B', np.int8), np.dtype(np.uint8): ('UChar_t', 'b', np.uint8),
       np.dtype(np.int16): ('Short_t', 'S', np.int16), np.dtype(np.uint16): ('UShort_t', 's', np.uint16),
       np.dtype(np.int32): ('int32_t', 'I', np.int32), np.dtype(np.uint32): ('uint32_t', 'i', np.uint32),
       np.dtype(np.int64): ('Long64_t', 'L', np.longlong), np.dtype(np.uint64): ('ULong64_t', 'l', np.ulonglong),
       np.dtype(np.float32): ('float', 'F', np.float32), np.dtype(np.float64): ('double', 'D', np.float64)}
_ROOT = None


def root():
    global _ROOT
    if _ROOT is None:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        if not hasattr(ROOT, 'eventtree'):
            ROOT.gInterpreter.Declare(WRITER.read_text())
        _ROOT = ROOT
    return _ROOT


def buffer(values, dtype):
    """Contiguous buffer typed for the cppyy converter (int64/uint64 as long long).

    cppyy cannot pass a zero-length array as a pointer; an empty column gets one unused element."""
    array = np.ascontiguousarray(values, dtype=dtype)
    if array.size == 0:
        array = np.zeros(1, dtype=dtype)
    return array.view(CXX[np.dtype(dtype)][2])


def title_of(branch):
    return branch.title + (' [%s]' % branch.unit if branch.unit and ('[' not in branch.title) else '')


