"""Fill the reserved photon-ID branches (pho_bdt_score, pho_bdt_state, pho_id_class) from a frozen package.

The package is PhotonID/models/nominal.json: hash-bound models, the view, the domain and the
working points. Au+Au photons are scored with ROOT's RBDT on the H70 shower view, the same call the sealed
schema-10 package reducer makes. pp photons carry the PPG12 producer score, checked against a recompute.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from .tables import require

_PHOTON_ID = Path(__file__).resolve().parents[3] / 'PhotonID'
FREEZE = _PHOTON_ID / 'models' / 'nominal.json'
STATE_SCORED, STATE_OUTSIDE, STATE_INVALID = 1, 0, 2
TIGHT, GAP, NONTIGHT, FAIL, UNCLASSIFIED = 3, 2, 1, 0, -1
# Model input name -> eventTree branch. pho_ branches are per photon, shv_ branches per photon at the model's
# shower view, vtx_z and cent_value per event.
BRANCH = {'cluster_Et': 'pho_et', 'cluster_Eta': 'pho_eta', 'vertexz': 'vtx_z', 'centrality': 'cent_value',
          'cluster_weta_cogx': 'shv_weta_cogx', 'cluster_wphi_cogx': 'shv_wphi_cogx',
          'cluster_weta33_cogx': 'shv_weta33_cogx', 'cluster_wphi33_cogx': 'shv_wphi33_cogx',
          'e11_over_e33': 'shv_e11_over_e33', 'e32_over_e35': 'shv_e32_over_e35',
          'cluster_et1': 'shv_et1', 'cluster_et2': 'shv_et2', 'cluster_et3': 'shv_et3', 'cluster_et4': 'shv_et4'}
CODES = {'pho_bdt_state': {STATE_SCORED: 'scored', STATE_OUTSIDE: 'outside the model domain',
                           STATE_INVALID: 'inside the domain with an invalid input feature'},
         'pho_id_class': {TIGHT: 'tight', GAP: 'between non-tight and tight (excluded)', NONTIGHT: 'non-tight',
                          FAIL: 'below non-tight', UNCLASSIFIED: 'not scored'}}


def load(path=FREEZE):
    """Verify nominal manifests, or resolve an explicitly supplied legacy package.

    The original descriptor retains its historical payload names and hash;
    installed model paths come from the nominal layout without changing science.
    """
    path = Path(path)
    package = json.loads(path.read_text())
    require(package.get('schema') == 'CanonicalPhotonIDFreezeV1', 'not a canonical photon-ID freeze')
    package_hash=hashlib.sha256(path.read_bytes()).hexdigest();models={}
    relocated=path.name=='nominal.json' and path.parent.name=='models'
    for system in ('auau', 'pp'):
        model = path.parent / package[system]['model']['file']
        if relocated:
            folder=path.parent/({'auau':'AuAu','pp':'PP'}[system])/'nominal'
            manifest=json.loads((folder/'manifest.json').read_text())
            require(manifest.get('schema')=='PhotonIDFrozenNominalV1' and manifest.get('freeze_sha256')==package_hash
                    and manifest.get('definition')==package[system] and manifest.get('system')=={'auau':'AuAu','pp':'PP'}[system]
                    and manifest.get('artifacts',{}).get('model.tmva.root')==package[system]['model']['sha256'],
                    system+' nominal manifest differs from the frozen definition')
            model=folder/'model.tmva.root'
        require(hashlib.sha256(model.read_bytes()).hexdigest() == package[system]['model']['sha256'],
                system + ' model file differs from the freeze')
        models[system]=model
    package['_sha256'] = package_hash
    package['_models'] = models
    package['_dir'] = path.parent
    return package


def evaluator(package, system):
    from .rootio import root
    ROOT = root()
    ROOT.TMVA.Tools.Instance()
    model_path=package.get('_models',{}).get(system,package['_dir'] / package[system]['model']['file'])
    model = ROOT.TMVA.Experimental.RBDT('myBDT', str(model_path))
    features = ROOT.std.vector('float')()

    def score(values):
        features.clear()
        for value in values:
            features.push_back(float(value))
        return float(model.Compute(features)[0])
    return score, ROOT.gROOT.GetVersion()


def feature_parity(order, column, features):
    """Photons whose model inputs differ from the named branches they were built from.

    column(branch) returns one value per photon (the model's shower view for shv_ branches); features is the
    (photons x inputs) matrix handed to the model. Values compare as float32, NaN equal to NaN."""
    rebuilt = np.column_stack([np.asarray(column(BRANCH[name]), dtype=np.float32) for name in order])
    stored = np.asarray(features, dtype=np.float32)
    same = (stored == rebuilt) | (np.isnan(stored) & np.isnan(rebuilt))
    return int(np.count_nonzero(~same.all(axis=1)))


def view_column(values, position, per_event, n):
    """Per-photon column accessor over flat converter values (shv_ at one view, per-event branches repeated)."""
    def column(branch):
        data = values[branch][1]
        if branch.startswith('shv_'):
            return data.reshape(n, 7)[:, position]
        if branch in ('vtx_z', 'cent_value'):
            return np.repeat(data[per_event[0]], per_event[1])
        return data
    return column


def classify_auau(score, centrality, wp):
    t70 = wp['WP70']['intercept'] + wp['WP70']['slope'] * centrality
    t80 = wp['WP80']['intercept'] + wp['WP80']['slope'] * centrality
    return np.where(score > t70, TIGHT, np.where((0.1 < score) & (score < t80), NONTIGHT,
                                                 np.where(score >= t80, GAP, FAIL))).astype(np.int32)


def classify_pp(score, et):
    tight = 0.815625 - 0.0015625 * et
    low = 0.7333333333333333 - 0.01333333333333333 * et
    high = 0.684375 + 0.0015625 * et
    return np.where(score > tight, TIGHT, np.where((low < score) & (score < high), NONTIGHT,
                                                   np.where(score >= high, GAP, FAIL))).astype(np.int32)


def augment(converter, path=FREEZE):
    """Add the three reserved branches to a converter's outputs; returns the per-file receipt."""
    package = load(path)
    values, meta, system = converter.out.values, converter.meta, converter.system
    spec = package[system]
    et, counts = values['pho_et'][1], values['pho_et'][2]
    n = len(et)
    per_event = counts[converter.ret_pos]
    names = meta['constants'].get('showerViews.definition_name')
    receipt = dict(schema='RJEventTreePhotonIDReceiptV1', freeze_sha256=package['_sha256'], system=system,
                   model_sha256=spec['model']['sha256'], photons=int(n),
                   codes={branch: {str(k): v for k, v in labels.items()} for branch, labels in CODES.items()})
    score = np.full(n, np.nan)
    state = np.full(n, STATE_OUTSIDE, dtype=np.int32)
    klass = np.full(n, UNCLASSIFIED, dtype=np.int32)
    if n:
        require(names is not None and spec['features']['view'] in names, 'no H70 shower view in this file')
        position = names.index(spec['features']['view'])
        stride = meta['strides']['shv_features']
        require(stride == spec['features']['count'], 'H70 view carries %d features, the model %d' % (stride, spec['features']['count']))
        if 'view_semantic_sha256' in spec['features']:
            require(meta['constants']['showerViews.semantic_sha256'][position] == spec['features']['view_semantic_sha256'],
                    'H70 view semantic hash differs from the model contract')
        if 'floor_gev' in spec['features']:
            floor = np.uint64(meta['constants']['showerViews.floor_gev'][position]['f64']).view(np.float64)
            require(floor == spec['features']['floor_gev'], 'H70 view floor differs from the model contract')
        features = values['shv_features'][1].reshape(n, 7, stride)[:, position, :]
        column = view_column(values, position, (converter.ret_pos, per_event), n)
        mismatches = feature_parity(spec['features']['order'], column, features)
        require(mismatches == 0, '%d photons: model inputs differ from the named %s shower branches' % (mismatches, spec['features']['view']))
        receipt['feature_parity'] = dict(photons=int(n), mismatches=0)
        finite_state = values['shv_finite_feature_state'][1].reshape(n, 7)[:, position]
        valid = (finite_state == 1) & np.all(np.isfinite(features), axis=1)
        compute, version = evaluator(package, system)
        receipt['root_version'] = version
        if system == 'auau':
            centrality = np.repeat(values['cent_value'][1][converter.ret_pos], per_event)
            domain = spec['domain']
            inside = ((et >= domain['et_min']) & (et < domain['et_max_exclusive']) & np.isfinite(centrality)
                      & (centrality >= domain['centrality_min']) & (centrality < domain['centrality_max_exclusive']))
            state[inside & ~valid] = STATE_INVALID
            scored = inside & valid
            for i in np.flatnonzero(scored):
                score[i] = compute(features[i])
            require(bool(np.all(np.isfinite(score[scored]))), 'a scored photon has a non-finite score')
            state[scored] = STATE_SCORED
            klass[scored] = classify_auau(score[scored], centrality[scored], spec['working_points'])
        else:
            producer = values['pho_prodbdt_score'][1]
            finite = values['pho_prodbdt_finite'][1]
            inside = et >= spec['domain']['et_min']
            score_valid=(finite == 1) & np.isfinite(producer) & (producer != -1.0)
            scored = inside & valid & score_valid
            state[inside & ~(valid & score_valid)] = STATE_INVALID
            # The producer scored on the farm's ROOT build; a different build may round the last float32 bits
            # differently. Accept at most 2 ulp, and never a different region.
            index = np.flatnonzero(scored)
            recomputed = np.array([compute(features[i]) for i in index])
            ulp = np.abs(recomputed - producer[index]) / np.spacing(np.float32(producer[index])).astype(np.float64)
            worst = float(ulp.max()) if len(index) else 0.0
            flips = int(np.count_nonzero(classify_pp(recomputed, et[index]) != classify_pp(producer[index], et[index])))
            require(worst <= 2.0, 'pp producer score differs from the recompute by %.2f float32 ulp' % worst)
            require(flips == 0, 'pp recompute changes the photon-ID region of %d photons' % flips)
            receipt['pp_recompute_max_float32_ulp'] = worst
            receipt['pp_recompute_region_flips'] = flips
            score[scored] = producer[scored]
            state[scored] = STATE_SCORED
            klass[scored] = classify_pp(score[scored], et[scored])
    for name, array in (('pho_bdt_score', score), ('pho_bdt_state', state), ('pho_id_class', klass)):
        converter.out.vector(name, array, per_event)
    receipt['states'] = {str(k): int((state == k).sum()) for k in (STATE_SCORED, STATE_OUTSIDE, STATE_INVALID)}
    receipt['classes'] = {str(k): int((klass == k).sum()) for k in (TIGHT, GAP, NONTIGHT, FAIL, UNCLASSIFIED)}
    return receipt
