"""Resolve the readable selection YAML once for all ROOT reducers.

No training or tree rewriting. Alternative scores must already exist in the
eventTree with a matching model receipt. The resolved definition is identical
for DATA and SIM and its digest is the output profile used by FinalAnalysis.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import re

import yaml

HERE = Path(__file__).resolve().parent
DEFAULT = HERE / 'config/selection.yaml'
CANONICAL_ID = HERE.parent / 'PhotonID/models/nominal.json'
CANONICAL_ISO = HERE / 'config/auau_isolation_wp80.json'
SHA = re.compile(r'[0-9a-f]{64}')
BRANCH = re.compile(r'[A-Za-z][A-Za-z0-9_]*')
AUAU_DATA_TRIGGERS = {
    'photon12': 'Photon 12 GeV + MBD NS >= 2, vtx < 150 cm',
}


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',', ':'),allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()


class UniqueLoader(yaml.SafeLoader):
    pass


def _mapping(loader, node, deep=False):
    pairs = loader.construct_pairs(node,deep=deep)
    result = {}
    for k,v in pairs:
        if k in result:
            raise ValueError('duplicate YAML key: '+str(k))
        result[k] = v
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def keys(value, allowed, label):
    if not isinstance(value,dict):
        raise ValueError(label+' must be a mapping')
    extra = set(value)-set(allowed)
    if extra:
        raise ValueError(label+': unknown settings '+', '.join(sorted(extra)))


def number(value, label):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError(label+' must be a finite number')
    return float(value)


def line(value, label):
    keys(value,('intercept','slope'),label)
    if set(value) != {'intercept','slope'}:
        raise ValueError(label+' needs intercept and slope')
    return {k:number(value[k],label+'.'+k) for k in ('intercept','slope')}


def load_package(path, schema):
    raw = Path(path).read_bytes(); p = json.loads(raw)
    if p.get('schema') != schema:
        raise ValueError(str(path)+': expected '+schema)
    return p,hashlib.sha256(raw).hexdigest()


def optional_selection(path, system):
    """Read a sealed PhotonID release through its JSON contract, without training imports."""
    path = Path(path).resolve(); folder = path.parent
    spec, _ = load_package(path, 'PhotonIDSelectionV1')
    release, _ = load_package(folder / 'manifest.json', 'PhotonIDReleaseV1')
    body = {k: v for k, v in release.items() if k != 'release_sha256'}
    if release.get('status') != 'PASS' or digest(body) != release.get('release_sha256'):
        raise ValueError('optional model release is not validated or its manifest changed')
    for name, pin in release['artifacts'].items():
        artifact = (folder / name).resolve()
        if not artifact.is_relative_to(folder) or file_sha(artifact) != pin:
            raise ValueError('optional model artifact changed: ' + name)
    if path.name not in release['artifacts'] or spec['working_points_file'] not in release['artifacts']:
        raise ValueError('selection and working points must be sealed release artifacts')
    if spec['system'] != {'pp': 'PP', 'auau': 'AuAu'}[system] or spec['system'] != release['system']:
        raise ValueError('optional selection collision system differs')
    if spec['state_branch'] is not None or spec['score_branch'] != release['score_branch']:
        raise ValueError('optional selection must use the released single score')
    if not BRANCH.fullmatch(spec['score_branch']) or spec['score_branch'] == 'pho_bdt_score':
        raise ValueError('optional score collides with nominal or has an invalid name')
    if spec['cells'] != release['cells'] or spec['domain'] != release['recipe']['domain']:
        raise ValueError('selection routing/domain differs from its released model')
    if spec['non_tight_lower'] != .1 or spec['tight'] != 'score > WP70' or spec['non_tight'] != '0.1 < score < WP80':
        raise ValueError('unsupported optional photon-ID inequalities')
    wp_path = (folder / spec['working_points_file']).resolve()
    wp, _ = load_package(wp_path, 'PhotonIDWorkingPointsV1')
    d = spec['domain']; et_bounds = [number(d['et_min'], 'ET min'), number(d['et_max_exclusive'], 'ET max')]
    cent_bounds = [number(d.get('centrality_min', 0.), 'centrality min'),
                   number(d.get('centrality_max_exclusive', 80.), 'centrality max')]
    if et_bounds[0] >= et_bounds[1] or cent_bounds[0] >= cent_bounds[1]:
        raise ValueError('invalid optional model domain')
    cells = []
    for cell in spec['cells']:
        bounds = cell['bounds']; keys(bounds, ('et', 'centrality_coarse', 'centrality_fine'), 'model cell')
        et = bounds.get('et', et_bounds)
        cent = bounds.get('centrality_coarse', bounds.get('centrality_fine', cent_bounds))
        if 'centrality_coarse' in bounds and 'centrality_fine' in bounds:
            raise ValueError('cell has ambiguous centrality routing')
        for axis, limits in ((et, et_bounds), (cent, cent_bounds)):
            if len(axis) != 2 or not limits[0] <= axis[0] < axis[1] <= limits[1]:
                raise ValueError('model cell lies outside its domain')
        w = wp['cells'][cell['id']]
        if w.get('status') != 'PASS' or w['axis'] != ('pho_et' if system == 'pp' else 'cent_value'):
            raise ValueError('model cell has no applicable working-point validation')
        cuts = {target: {k: number(w['cuts'][source][k], source + '.' + k) for k in ('intercept', 'slope')}
                for target, source in (('tight', 'WP70'), ('non_tight_high', 'WP80'))}
        for x in (et if system == 'pp' else cent):
            tight, high = (cuts[k]['intercept'] + cuts[k]['slope'] * x for k in ('tight', 'non_tight_high'))
            if not .1 < high <= tight <= 1.:
                raise ValueError('optional working-point bands overlap or are invalid')
        cells.append(dict(id=cell['id'], et=et, centrality=cent, **cuts))
    if len({c['id'] for c in cells}) != len(cells) or set(wp['cells']) != {c['id'] for c in cells}:
        raise ValueError('duplicate or missing working-point cells')
    # Check each elementary rectangle: endpoints alone cannot detect gaps in a
    # two-dimensional grid. Half-open bounds use the same rule as scoring.
    ets = sorted(set(et_bounds + [x for c in cells for x in c['et']]))
    cents = sorted(set(cent_bounds + [x for c in cells for x in c['centrality']]))
    for et in [(a+b)/2 for a,b in zip(ets, ets[1:])]:
        for cent in [(a+b)/2 for a,b in zip(cents, cents[1:])]:
            if sum(c['et'][0] <= et < c['et'][1] and c['centrality'][0] <= cent < c['centrality'][1] for c in cells) != 1:
                raise ValueError('optional model cells overlap or leave a gap')
    return dict(score_branch=spec['score_branch'], state_branch=None, score_metadata='meta/photon_models',
                model_name=release['name'], release_sha256=release['release_sha256'], domain=d,
                cells=cells, selection_sha256=file_sha(path), working_points_sha256=file_sha(wp_path))


def resolve(path=DEFAULT, data_trigger=None):
    path = Path(path).resolve()
    config = yaml.load(path.read_text(),Loader=UniqueLoader)
    keys(config,('version','name','photon_id_package','auau_isolation_package','systems'),'selection')
    if config.get('version') != 1 or isinstance(config['version'],bool):
        raise ValueError('selection version must be 1')
    if not isinstance(config.get('name'),str) or not config['name'].strip():
        raise ValueError('selection name is required')
    def package_path(key, default):
        value = config.get(key)
        return (path.parent/value).resolve() if value is not None else default
    pid,pid_sha = load_package(package_path('photon_id_package',CANONICAL_ID),'CanonicalPhotonIDFreezeV1')
    iso,iso_sha = load_package(package_path('auau_isolation_package',CANONICAL_ISO),'CanonicalIsolationFreezeV1')
    if iso.get('system') != 'auau' or iso['witness']['reconstructed'] is not True:
        raise ValueError('AuAu isolation package must use a reconstructed witness')
    canonical_id,_ = load_package(CANONICAL_ID,'CanonicalPhotonIDFreezeV1')
    overrides = config.get('systems',{})
    keys(overrides,('pp','auau'),'systems')
    resolved = dict(schema='PhotonJetSelectionV1',version=1,systems={})
    for system in ('pp','auau'):
        spec = pid[system]
        if not SHA.fullmatch(spec['model']['sha256']):
            raise ValueError('invalid model sha256 for '+system)
        if system == 'auau':
            cuts = dict(tight=line(spec['working_points']['WP70'],'WP70'),
                        non_tight_low=dict(intercept=0.1,slope=0.0),
                        non_tight_high=line(spec['working_points']['WP80'],'WP80'))
            isolation = dict(radius=number(iso['witness']['radius'],'radius'),subtraction_method=2,
                threshold=dict(intercept=number(iso['threshold']['intercept_gev'],'isolation intercept'),
                               slope=number(iso['threshold']['slope_gev_per_percent'],'isolation slope')),
                sideband_gap_gev=number(iso['sideband_gap_gev'],'isolation gap'),
                centrality_min=number(iso['domain']['centrality_min'],'centrality min'),
                centrality_max=number(iso['domain']['centrality_max_exclusive'],'centrality max'),
                package_sha256=iso_sha)
            if iso['witness']['subtraction_method'] != 2 or iso['threshold']['variable'] != 'cent_value':
                raise ValueError('AuAu package requires SUB1 and a centrality threshold')
        else:
            # The frozen pp contract expresses these PPG12 lines as readable formulae.
            # A replacement package may instead supply explicit working_points.
            if spec.get('working_points'):
                cuts = {k:line(spec['working_points'][k],k) for k in ('tight','non_tight_low','non_tight_high')}
            else:
                if spec['regions'] != canonical_id['pp']['regions']:
                    raise ValueError('a different pp ID definition needs explicit working_points')
                cuts = dict(tight=dict(intercept=0.815625,slope=-0.0015625),
                    non_tight_low=dict(intercept=0.7333333333333333,slope=-0.01333333333333333),
                    non_tight_high=dict(intercept=0.684375,slope=0.0015625))
            isolation = dict(radius=0.4,subtraction_method=1,threshold='stored',sideband_gap_gev=None,
                             centrality_min=0.0,centrality_max=80.0,package_sha256=None)
        original_cuts,original_iso = copy.deepcopy(cuts),copy.deepcopy(isolation)
        change = dict(overrides.get(system,{}) or {})
        if system=='auau' and data_trigger is not None:
            change['data_trigger']=data_trigger
        keys(change,('score_branch','state_branch','score_metadata','tight','non_tight_low','non_tight_high','isolation','data_trigger','photon_id_selection'),system)
        trigger = change.get('data_trigger','photon12' if system=='auau' else 'photon4')
        if (system=='auau' and trigger not in AUAU_DATA_TRIGGERS) or (system=='pp' and trigger!='photon4'):
            raise ValueError('nominal DATA triggers are photon12 for AuAu and photon4 for pp')
        for k in cuts:
            if k in change:
                cuts[k] = line(change[k],system+'.'+k)
        iso_change = change.get('isolation',{})
        keys(iso_change,('radius','threshold','sideband_gap_gev'),system+'.isolation')
        for k,v in iso_change.items():
            isolation[k] = line(v,system+'.isolation.threshold') if k=='threshold' else number(v,system+'.isolation.'+k)
        if isolation['radius'] not in (0.3,0.4):
            raise ValueError('isolation radius must be a stored witness: 0.3 or 0.4')
        if system=='auau' and isolation['radius'] != original_iso['radius']:
            raise ValueError('a new AuAu radius requires a matching isolation package')
        if system=='pp' and isolation['radius'] != 0.4 and 'threshold' not in iso_change:
            raise ValueError('a changed pp radius needs its explicit threshold definition')
        if isolation['sideband_gap_gev'] is not None and isolation['sideband_gap_gev'] < 0:
            raise ValueError('isolation sideband gap cannot be negative')
        if not 0 <= isolation['centrality_min'] < isolation['centrality_max'] <= 80:
            raise ValueError('isolation centrality domain is outside the retained nominal domain')
        # Lines must not create overlapping tight/non-tight regions in the analysis domain.
        for x in ((15.0,35.0) if system=='pp' else (0.0,80.0)):
            vals = {k:v['intercept']+v['slope']*x for k,v in cuts.items()}
            if not vals['non_tight_low'] < vals['non_tight_high'] <= vals['tight']:
                raise ValueError(system+': tight and non-tight cuts overlap or have an empty band')
        branches = {k:change.get(k,v) for k,v in dict(score_branch='pho_bdt_score',state_branch='pho_bdt_state',score_metadata='meta/photon_id').items()}
        if not all(isinstance(branches[k],str) and BRANCH.fullmatch(branches[k]) for k in ('score_branch','state_branch')):
            raise ValueError('score/state must name eventTree branches')
        if not re.fullmatch(r'meta/[A-Za-z0-9_/]+',str(branches['score_metadata'])):
            raise ValueError('score_metadata must name a meta/ receipt')
        if (branches['score_branch'] != 'pho_bdt_score' or branches['state_branch'] != 'pho_bdt_state') and branches['score_metadata']=='meta/photon_id':
            raise ValueError('an alternative score needs its own augmentation metadata')
        nominal_score = (branches == dict(score_branch='pho_bdt_score',state_branch='pho_bdt_state',score_metadata='meta/photon_id')
                         and spec['model']['sha256']==canonical_id[system]['model']['sha256'])
        resolved['systems'][system] = dict(**branches,model_sha256=spec['model']['sha256'],
            photon_id_package_sha256=pid_sha,**cuts,isolation=isolation,data_trigger=trigger,
            event_selection=dict(vertex_abs_z_max_cm=10.0 if system=='auau' else 60.0,
                require_minimum_bias=system=='auau', trigger_word='scaled',
                trigger_name=AUAU_DATA_TRIGGERS[trigger] if system=='auau' else 'Photon 4 GeV + MBD NS >= 1',
                data_run_periods={} if system=='auau' else {'0mrad':[47289,51274],'1p5mrad':[51274,54000]}),
            verify_nominal_id=nominal_score and cuts==original_cuts and pid_sha==file_sha(CANONICAL_ID),
            verify_nominal_isolation=isolation==original_iso and (system=='pp' or iso_sha==file_sha(CANONICAL_ISO)))
        if change.get('photon_id_selection'):
            if set(change) & {'score_branch', 'state_branch', 'score_metadata', 'tight', 'non_tight_low', 'non_tight_high'}:
                raise ValueError('optional release owns its score and working points; do not override them')
            optional = optional_selection(path.parent / change['photon_id_selection'], system)
            resolved['systems'][system].update({k: optional[k] for k in ('score_branch', 'state_branch', 'score_metadata')})
            resolved['systems'][system].update(optional_model=optional, verify_nominal_id=False)
    resolved['selection_sha256'] = digest(resolved)
    # Exact identity is retained even if two configurations use the same human label.
    resolved['profile'] = 'selection:'+resolved['selection_sha256']
    resolved['name'] = config['name']
    return resolved


def verify_resolved(value):
    body = {k:v for k,v in value.items() if k not in ('selection_sha256','profile','name')}
    if value.get('selection_sha256') != digest(body) or value.get('profile') != 'selection:'+digest(body):
        raise ValueError('resolved selection digest/profile mismatch')
    return value


def compatible_outputs(paths):
    """Refuse a concatenation of results produced under different selections."""
    import uproot
    seen = set()
    for path in paths:
        with uproot.open(path) as f:
            if 'selection' in f:
                value = verify_resolved(json.loads(str(f['selection'])))
                if str(f['profile']) != value['profile']:
                    raise ValueError(str(path)+': selection/profile metadata differ')
                seen.add(value['profile'])
            else:
                profile = str(f['profile']) if 'profile' in f else 'legacy_unbound'
                if profile.startswith('selection:'):
                    raise ValueError(str(path)+': missing resolved selection metadata')
                seen.add(profile)
    if len(seen)>1:
        raise ValueError('cannot combine different selections: '+', '.join(sorted(seen)))
