"""Fill the reserved isolation branch pho_iso_class from a frozen Au+Au isolation package.

The package (CanonicalIsolationFreezeV1, written by TreeToHists/isolation.py freeze) holds the witness it applies to
(R = 0.3, reconstructed, its subtraction method), the threshold line in centrality at the nominal efficiency, the
sideband gap and the centrality domain. Per photon, with t = intercept + slope * cent_value:

    2 isolated        cone sum < t
    1 gap             t <= cone sum <= t + gap (unclassified)
    0 non-isolated    cone sum > t + gap
   -1 invalid         witness not reconstructed or of another method, a non-finite cone sum or centrality, or a
                      centrality outside the package domain; never background

pp isolation is PPG12's R=0.4 definition, using the reconstructed producer
witness. Its cone already includes the producer's applicable SIM correction;
applying that correction a second time here would change the selection.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from .tables import require,decode

ISOLATED, GAP, NONISOLATED, INVALID = 2, 1, 0, -1
CODES = {'pho_iso_class': {ISOLATED: 'isolated', GAP: 'between isolated and non-isolated (excluded)',
                           NONISOLATED: 'non-isolated', INVALID: 'invalid input'}}
RECONSTRUCTED = 0  # isow_reconstructed_or_truth of a reconstructed witness

PP_DEFINITION = dict(schema='PPG12IsolationDefinitionV1', system='pp', radius=0.4,
                     stride_index=1, subtraction_method=1, threshold='stored reconstructed witness',
                     isolated='cone < threshold', nonisolated='cone > sideband_threshold',
                     boundaries='gap', sim_correction='already in stored witness; do not repeat')


def check_radius(converter,n,index,radius):
    if not n:return
    observed=converter.meta['constants'].get('isolationWitnesses.radius')
    require(observed is not None and len(observed)==2,'isolation witness radius constants missing')
    require(float(decode(observed[index]))==radius,'isolation witness radius differs from nominal definition')


def classify_pp(cone, threshold, sideband, reconstructed, method):
    cone, threshold, sideband = [np.asarray(v, dtype=float) for v in (cone, threshold, sideband)]
    valid = ((np.asarray(reconstructed) == RECONSTRUCTED) & (np.asarray(method) == 1)
             & np.isfinite(cone) & np.isfinite(threshold) & np.isfinite(sideband) & (sideband >= threshold))
    klass = np.where(cone < threshold, ISOLATED, np.where(cone > sideband, NONISOLATED, GAP))
    return np.where(valid, klass, INVALID).astype(np.int32)


def augment_pp(converter):
    require(converter.system == 'pp', 'PPG12 isolation requires pp')
    values = converter.out.values
    n = len(values['pho_et'][1])
    check_radius(converter,n,1,0.4)
    per_event = values['pho_et'][2][converter.ret_pos]
    def witness(name):
        require(len(values[name][1]) == 2*n, 'isolation witnesses are not 2 per photon: '+name)
        return values[name][1].reshape(n, 2)[:, 1]
    klass = classify_pp(*(witness(k) for k in ('isow_cone_sum', 'isow_threshold',
        'isow_sideband_threshold', 'isow_reconstructed_or_truth', 'isow_subtraction_method')))
    converter.out.vector('pho_iso_class', klass, per_event)
    pin = hashlib.sha256(json.dumps(PP_DEFINITION,sort_keys=True,separators=(',', ':')).encode()).hexdigest()
    return dict(schema='RJEventTreeIsolationReceiptV1', system='pp', photons=n,
                freeze_sha256=pin, derivation_sha256=pin, definition=PP_DEFINITION,
                codes={branch:{str(k):v for k,v in labels.items()} for branch,labels in CODES.items()},
                classes={str(k):int((klass == k).sum()) for k in (ISOLATED,GAP,NONISOLATED,INVALID)})


def load(path):
    path = Path(path)
    package = json.loads(path.read_text())
    require(package.get('schema') == 'CanonicalIsolationFreezeV1', 'not a canonical isolation freeze')
    require(package.get('system') == 'auau', 'the isolation freeze applies to Au+Au')
    package['_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return package


def classify(cone, centrality, reconstructed, method, package):
    """pho_iso_class of photons given their witness cone sum, event centrality, witness kind and method."""
    cone, centrality = np.asarray(cone, dtype=float), np.asarray(centrality, dtype=float)
    line, witness, domain = package['threshold'], package['witness'], package['domain']
    threshold = line['intercept_gev'] + line['slope_gev_per_percent'] * centrality
    valid = ((np.asarray(reconstructed) == RECONSTRUCTED) & (np.asarray(method) == witness['subtraction_method'])
             & np.isfinite(cone) & np.isfinite(centrality)
             & (centrality >= domain['centrality_min']) & (centrality < domain['centrality_max_exclusive']))
    klass = np.where(cone < threshold, ISOLATED, np.where(cone > threshold + package['sideband_gap_gev'], NONISOLATED, GAP))
    return np.where(valid, klass, INVALID).astype(np.int32)


def augment(converter, path):
    """Add pho_iso_class to a converter's outputs; returns the per-file receipt."""
    package = load(path)
    require(converter.system == 'auau', 'pp isolation is applied from the stored witness thresholds; no package')
    values = converter.out.values
    counts = values['pho_et'][2]
    per_event = counts[converter.ret_pos]
    n = len(values['pho_et'][1])
    index = package['witness']['stride_index']
    check_radius(converter,n,index,package['witness']['radius'])
    require(len(values['isow_cone_sum'][1]) == 2 * n, 'isolation witnesses are not 2 per photon')
    cone = values['isow_cone_sum'][1].reshape(n, 2)[:, index]
    kind = values['isow_reconstructed_or_truth'][1].reshape(n, 2)[:, index]
    method = values['isow_subtraction_method'][1].reshape(n, 2)[:, index]
    centrality = np.repeat(values['cent_value'][1][converter.ret_pos], per_event)
    klass = classify(cone, centrality, kind, method, package)
    converter.out.vector('pho_iso_class', klass, per_event)
    return dict(schema='RJEventTreeIsolationReceiptV1', freeze_sha256=package['_sha256'], system='auau', photons=int(n),
                derivation_sha256=package['derivation']['package_sha256'],
                codes={branch: {str(k): v for k, v in labels.items()} for branch, labels in CODES.items()},
                classes={str(k): int((klass == k).sum()) for k in (ISOLATED, GAP, NONISOLATED, INVALID)})
