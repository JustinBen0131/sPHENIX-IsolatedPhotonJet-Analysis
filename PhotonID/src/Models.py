"""Models: immutable photon-ID releases and their scientific definitions.

The command driver uses this module to resolve recipes, verify model bytes and
resume only work bound to the same inputs, code and environment. Generated
artifacts never replace the bundled nominal models. A release becomes usable
only after its numerical and working-point validation has passed.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys

import yaml

HERE = Path(__file__).resolve().parents[1]
SYSTEMS = {'AuAu': 'auau', 'PP': 'pp'}
NAMES = {
    'nominal_retrained_v1': ('bdt', 'pho_bdt_score_nominal_retrained_v1', ()),
    'centCoarse_v1': ('bdt', 'pho_bdt_score_cent_coarse_v1', ('centrality_coarse',)),
    'centFine_v1': ('bdt', 'pho_bdt_score_cent_fine_v1', ('centrality_fine',)),
    'ET_v1': ('bdt', 'pho_bdt_score_et_v1', ('et',)),
    'ETcentCoarse_v1': ('bdt', 'pho_bdt_score_et_cent_coarse_v1', ('et', 'centrality_coarse')),
    'ETcentFine_v1': ('bdt', 'pho_bdt_score_et_cent_fine_v1', ('et', 'centrality_fine')),
    'MLP_v1': ('mlp', 'pho_mlp_score_v1', ()),
    'logReg_v1': ('logreg', 'pho_logreg_score_v1', ()),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    """Create once: an existing different receipt is an error, not an overwrite."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = canonical(value) + '\n'
    try:
        with path.open('x') as out:
            out.write(raw)
            out.flush()
            os.fsync(out.fileno())
    except FileExistsError:
        require(path.read_text() == raw, 'immutable artifact differs: ' + str(path))


def read_json(path):
    return json.loads(Path(path).read_text())


def environment():
    versions = {}
    for name in ('numpy', 'scipy', 'uproot', 'awkward', 'xgboost', 'PyYAML', 'ROOT'):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return dict(python=sys.version, executable=sys.executable, packages=versions)


def load_config(path=None):
    path = Path(path or HERE / 'config/training.yaml').resolve()
    config = yaml.safe_load(path.read_text())
    require(config.get('version') == 1, 'unsupported training configuration')
    require(config['training']['train_fraction'] == .9, 'nominal training fraction must be 0.9')
    require(config['training']['seed'] == 13, 'nominal split seed must be 13')
    require(config['training']['minimum_bias_branch'] == 'cent_mb_decision' and
            config['training']['minimum_bias_pass'] == 1, 'corrected minimum-bias witness changed')
    require(len(config['features']['AuAu']) == 14 and len(config['features']['PP']) == 11,
            'nominal feature counts differ')
    expected = [n for n in config['features']['AuAu'] if n not in
                ('cluster_weta33_cogx', 'cluster_wphi33_cogx', 'centrality')]
    require(config['features']['PP'] == expected, 'pp feature order differs from nominal projection')
    require(config['weights'] == dict(bins=20, eta_range=[-.7, .7], et_cap=800),
            'canonical weighting settings changed')
    for edges in config['bins'].values():
        require(len(edges) >= 2 and all(a < b for a, b in zip(edges, edges[1:])), 'invalid bin edges')
    return config


def recipe(config, system, name):
    require(system in SYSTEMS and name in NAMES, 'unknown system/model')
    require(system == 'AuAu' or name == 'nominal_retrained_v1', 'this optional model is AuAu only')
    algorithm, branch, axes = NAMES[name]
    domain = dict(et_min=15., et_max_exclusive=35., eta_abs_max=.7, vertex_abs_max_cm=10.)
    if system == 'AuAu':
        domain.update(centrality_min=0., centrality_max_exclusive=80.)
    return dict(schema='PhotonIDRecipeV1', name=name, system=system,
                algorithm=algorithm, score_branch=branch, features=config['features'][system],
                domain=domain, axes={a: config['bins'][a] for a in axes},
                training=config['training'], weights=config['weights'],
                learner=config['learners'][algorithm], validation=config['validation'],
                working_points=config['working_points'],
                sources=config['sources'][system],
                feature_view='H70', isolation_in_training=False)


def request_identity(recipe, inputs):
    """Bind semantics and backend source bytes; code changes require a new release."""
    code = {p.name: file_hash(p) for p in sorted((HERE / 'src').glob('*.py'))}
    code['photon_id.py'] = file_hash(HERE / 'scripts/photon_id.py')
    return digest(dict(recipe=recipe, inputs=inputs, code=code, environment=environment()))


def release_dir(root, system, name):
    require(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', name) is not None, 'invalid release name')
    require(system in SYSTEMS, 'invalid collision system')
    return Path(root).resolve() / system / name


def load_release(folder, *, request=None):
    folder = Path(folder).resolve()
    manifest = read_json(folder / 'manifest.json')
    require(manifest.get('schema') == 'PhotonIDReleaseV1' and manifest.get('status') == 'PASS',
            'model release is not validated')
    body = {k: v for k, v in manifest.items() if k != 'release_sha256'}
    require(digest(body) == manifest['release_sha256'], 'model manifest changed')
    if request is not None:
        require(manifest['request_sha256'] == request,
                'existing model has different inputs/recipe/environment; choose a new release version')
    for relative, expected in manifest['artifacts'].items():
        p = folder / relative
        require(p.resolve().is_relative_to(folder), 'artifact path escapes model release')
        require(file_hash(p) == expected, 'model artifact changed: ' + relative)
    return manifest


def list_models(root, system=None):
    result = []
    for collision in ([system] if system else SYSTEMS):
        for p in sorted((Path(root) / collision).glob('*/manifest.json')):
            m = read_json(p)
            if m.get('schema') == 'PhotonIDReleaseV1':
                load_release(p.parent)
            elif m.get('schema') == 'PhotonIDFrozenNominalV1':
                for relative, h in m['artifacts'].items():
                    artifact = (p.parent / relative).resolve()
                    require(artifact.is_relative_to(p.parent.resolve()), 'nominal artifact path escapes package')
                    require(file_hash(artifact) == h, 'bundled nominal payload changed')
            else:
                raise ValueError('unknown model manifest: ' + str(p))
            result.append(dict(system=collision, name=p.parent.name,
                               score_branch=m['score_branch'], status=m['status'], path=str(p.parent)))
    return result


def progress(path, phase, processed, total, **extra):
    """Atomic in-process progress; terminal receipts remain the validation authority."""
    import datetime
    value = dict(phase=phase, processed=int(processed), total=int(total),
                 remaining=int(total - processed), updated_at=datetime.datetime.now(
                     datetime.timezone.utc).isoformat(), **extra)
    import resource
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    value['peak_memory_bytes'] = int(maximum if sys.platform == 'darwin' else maximum * 1024)
    require(0 <= processed <= total, 'invalid progress count')
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.pending')
    temp.write_text(canonical(value) + '\n')
    os.replace(temp, path)
