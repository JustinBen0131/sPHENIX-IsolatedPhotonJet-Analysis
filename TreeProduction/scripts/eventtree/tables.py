"""Shared helpers for the converter and the inverse: exact comparison, exact JSON values, and the
fixed facts of the production final layout (tree paths, identity domains, collection orders)."""
import json

import awkward as ak
import numpy as np

REPLAY = 'ReplayFoundationV1/'
# Production final-layout trees, in file order.
TREES = ('sources', 'events', 'photons', 'models', 'showerCells', 'showerViews', 'isolationConstituents',
         'isolationWitnesses', 'jets', 'jetConstituents', 'photonJets', 'truthPhotons', 'truthPhotonMisses',
         REPLAY + 'RJEmbeddedPhotonDiagnosticOccurrenceV1', REPLAY + 'RJPPG12DiagnosticOccurrenceV1',
         REPLAY + 'RJEventDisplaySnapshotV1', REPLAY + 'RJUpstreamRejectedEventV1', 'truthJets', 'recoTruthLinks',
         'weights', 'triggerScalers', 'triggerRunInfo', 'photonAssociations', 'photonResponseLinks',
         'eventAugmentations', 'photonAugmentations')
# Reco jet collections in meta/collections order; the analysis view first.
PP_JETS = (('jet02', 'towerinfo_calib', 2), ('jet03', 'towerinfo_calib', 3), ('jet04', 'towerinfo_calib', 4))
AUAU_JETS = (('jet02', 'towerinfo_sub1', 2), ('jet03', 'towerinfo_sub1', 3), ('jet04', 'towerinfo_sub1', 4),
             ('jetunsub02', 'towerinfo_retower_unsubtracted', 2), ('jetunsub03', 'towerinfo_retower_unsubtracted', 3),
             ('jetunsub04', 'towerinfo_retower_unsubtracted', 4))
TRUTH_JETS = (('truthjet02', 2), ('truthjet03', 3), ('truthjet04', 4))
LINK_JET, LINK_PHOTON, LINK_NONE = 2, 1, 0
NO_RESPONSE_ROW = -128
REPAIR_VERSION = 'CanonicalPhotonRepairTablesV1'


class ConversionError(ValueError):
    """The part does not satisfy an invariant the lossless restructure relies on; nothing is written."""


def require(ok, message):
    if not ok:
        raise ConversionError(message)


def bits(values):
    """Exact comparable view: floats by bit pattern, everything else as is."""
    values = np.asarray(values)
    if values.dtype == np.float64:
        return values.view(np.uint64)
    if values.dtype == np.float32:
        return values.view(np.uint32)
    return values


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind == 'f' or b.dtype.kind == 'f':
        return a.dtype == b.dtype and bool(np.array_equal(bits(a), bits(b)))
    if a.dtype.kind in 'OUS' or b.dtype.kind in 'OUS':
        return a.tolist() == b.tolist()
    return bool(np.array_equal(a, b))


def encode(value):
    """Exact JSON value: floats as their bit pattern, integers as integers, strings as strings."""
    if isinstance(value, (np.floating, float)):
        if isinstance(value, np.float32):
            return {'f32': int(np.float32(value).view(np.uint32))}
        return {'f64': int(np.float64(value).view(np.uint64))}
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode()
    return value


def decode(value, dtype=None):
    if isinstance(value, dict) and 'f64' in value:
        return np.uint64(value['f64']).view(np.float64)
    if isinstance(value, dict) and 'f32' in value:
        return np.uint32(value['f32']).view(np.float32)
    if dtype is not None and dtype is not object:
        return np.dtype(dtype).type(value)
    return value


def dumps(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def flat(column):
    """(flat numpy values, per-row counts) of a jagged column."""
    return ak.to_numpy(ak.flatten(column, axis=1)), ak.to_numpy(ak.num(column, axis=1)).astype(np.int64)


def offsets(counts):
    out = np.zeros(len(counts) + 1, dtype=np.int64)
    np.cumsum(counts, out=out[1:])
    return out


def group_counts(parent_index, n_parents, what):
    """Per-parent row counts; the rows must be grouped by parent in nondecreasing order."""
    parent_index = np.asarray(parent_index, dtype=np.int64)
    if len(parent_index):
        require(parent_index.min() >= 0 and parent_index.max() < n_parents, what + ': parent index out of range')
        require(bool(np.all(np.diff(parent_index) >= 0)), what + ': rows are not grouped by parent')
    return np.bincount(parent_index, minlength=n_parents).astype(np.int64)


def local_index(parent_counts):
    """Index of each child row within its parent, for rows grouped by parent."""
    starts = offsets(parent_counts)[:-1]
    return np.arange(int(parent_counts.sum()), dtype=np.int64) - np.repeat(starts, parent_counts)
