"""Run maintained DST populations through the direct eventTree producer.

The C++ producer owns reconstruction and captured branch filling. This driver
resolves source identities, launches Fun4All, checks its saved completion and
calls shared nominal finalization before collision-safe publication.
"""
import argparse
import datetime as dt
import fcntl
import time
from collections import deque
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

TREE_DIR = Path(__file__).resolve().parents[1]
BINDINGS = TREE_DIR / "config" / "production_populations.json"
MACRO = TREE_DIR / "macros" / "Fun4All_PhotonJetTree.C"
PIECE_EVENTS = {"auau_data": 35000, "default": 100000}

def format_modules():
    """Resolve the sibling public package, never private tests or agent paths."""
    path = TREE_DIR / 'scripts'
    require((path / 'eventtree/schema.py').is_file(), 'the TreeProduction eventTree module is required')
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


class FrontendError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise FrontendError(message)


def text(value, name, allow_empty=False):
    require(isinstance(value, str) and (allow_empty or value.strip()), f"missing {name}")
    require(not any(ord(c) < 32 for c in value), f"control character in {name}")
    return value


def integer(value, name, minimum=0, maximum=2147483647):
    require(type(value) is int and minimum <= value <= maximum,
            f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path):
    result = json.loads(path.read_text(), object_pairs_hook=unique_object)
    require(isinstance(result, dict), f"expected a JSON object: {path}")
    return result


def known_keys(value, allowed, context):
    unknown = set(value) - set(allowed.split())
    require(not unknown, f"unsupported {context} fields: {', '.join(sorted(unknown))}")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pinned_file(base, value, expected, name):
    path = (base / text(value, name)).resolve()
    require(path.is_file(), f"{name} is not a readable file: {path}")
    require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected),
            f"{name} requires a lowercase SHA256")
    require(sha256(path) == expected, f"{name} SHA256 mismatch: {path}")
    return path


def plan_population(lane, binding, settings, args):
    """Resolve job metadata only; never interpret DST columns or YAML policy."""


    manifest = pinned_file(BINDINGS.parent, binding.get("manifest"),
                           binding.get("manifest_sha256"), f"{lane} manifest")
    population = read_json(manifest)
    known_keys(population, "schema_version lane campaign sources", "manifest")
    require(population.get("schema_version") == 1, f"unsupported manifest schema: {manifest}")
    require(population.get("lane") == lane, f"manifest lane mismatch: {manifest}")
    require(population.get("campaign") == settings["campaign"], f"manifest campaign mismatch: {manifest}")
    sources = population.get("sources")
    require(isinstance(sources, list) and sources, f"empty/missing source population: {manifest}")


    profile = text(binding.get("profile"), f"lanes.{lane}.profile")
    config = (BINDINGS.parent / text(settings.get("production_config"), "production_config")).resolve()
    require(config.is_file(), f"production configuration is missing: {config}")


    jobs, ordinals, bundles = [], set(), set()
    for source in sources:
        require(isinstance(source, dict), f"{lane}: source must be an object")
        known_keys(source, "source_file_ordinal input_list input_list_sha256 run segment "
                          "first_entry event_count sample period si_di_role production_config "
                          "production_config_sha256 logical_source_id input_file_sha256", "source")
        require(("production_config" in source) == ("production_config_sha256" in source),
                "source production_config and production_config_sha256 must be supplied together")
        ordinal = integer(source.get("source_file_ordinal"), "source_file_ordinal", maximum=2**63 - 1)
        require(ordinal not in ordinals, f"{lane}: duplicate source ordinal {ordinal}")
        ordinals.add(ordinal)
        first = integer(source.get("first_entry"), "first_entry")
        events = integer(source.get("event_count"), "event_count", minimum=0)
        inputs = pinned_file(manifest.parent, source.get("input_list"),
                             source.get("input_list_sha256"), "input_list")


        require(source["input_list_sha256"] not in bundles, f"{lane}: duplicate input-list bundle")
        bundles.add(source["input_list_sha256"])
        job_config = config
        if "production_config" in source:
            job_config = pinned_file(manifest.parent, source["production_config"],
                                     source.get("production_config_sha256"), "source production_config")


        job = dict(profile=profile, input_list=str(inputs), output_file="",
                   run=integer(source.get("run"), "run"),
                   segment=integer(source.get("segment"), "segment"),
                   events=min(events, args.events) if args.events is not None else events,
                   first_entry=first, config=str(job_config), source_file_ordinal=ordinal,
                   input_file_sha256=source.get("input_file_sha256", ""), source_manifest_sha256=binding["manifest_sha256"],
                   sample=text(source.get("sample"), "sample"),
                   period=text(source.get("period"), "period", allow_empty=True),
                   si_di_role=text(source.get("si_di_role"), "si_di_role", allow_empty=True),
                   input_list_sha256=source["input_list_sha256"],
                   configuration_sha256=sha256(job_config), macro_sha256=sha256(MACRO), verbosity=0)
        job['logical_source_id'] = source.get('logical_source_id') or (
            lane + ':' + str(job['run']) + ':' + str(job['segment']) + ':' + str(ordinal))
        text(job['logical_source_id'], 'logical_source_id')
        require(not job['input_file_sha256'] or re.fullmatch(r'[0-9a-f]{64}', job['input_file_sha256']),
                'input_file_sha256 must be a SHA256 when supplied')
        job['source_events'] = events
        job['complete_source'] = first == 0 and job['events'] == events
        job['resume'] = getattr(args, 'resume', False)
        job['timeout'] = getattr(args, 'timeout', 43200)
        selected = ((args.source is None or ordinal == args.source)
                    and (getattr(args, 'run', None) is None or job['run'] == args.run)
                    and (getattr(args, 'sample', None) is None or job['sample'] == args.sample))
        if selected:
            jobs.append(job)

    require(jobs, f"{lane}: requested source/run/sample does not exist; nothing was selected")


    selected_output = getattr(args, 'output', None)
    output = (selected_output.resolve() if selected_output else
              (BINDINGS.parent / text(settings.get("output_root"), "output_root")).resolve())
    output /= settings["campaign"]
    development = args.events is not None or args.source is not None
    if development:
        output = output / "development" / lane / (
            f"source-{args.source if args.source is not None else 'all'}_"
            f"events-{args.events if args.events is not None else 'manifest'}")
    else:
        output /= lane
    if getattr(args, 'run', None) is not None:
        output /= f"run-{args.run}"
    if getattr(args, 'sample', None) is not None:
        require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.sample), "sample must be a safe exact identifier")
        output /= args.sample
    require(not output.exists() or getattr(args, "resume", False), f"output already exists; use --resume for exact verified work: {output}")

    for job in jobs:
        group = str(job['run']) if lane.endswith('_data') else job['sample']
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_. -]*', group), 'unsafe source group name')
        job['output_file'] = str(output / group / f"part-{job['source_file_ordinal']:08d}.root")

    return dict(lane=lane, profile=profile, manifest=str(manifest),
                manifest_sha256=binding["manifest_sha256"], output=str(output),
                development=development, population_sources=len(sources), jobs=jobs)


def macro_argument(job, check_only=False):
    # JSON string quoting is also valid for these C++ UTF-8 string literals.
    # subprocess receives an argument vector: no shell, eval or interpolation.
    values = [job[key] for key in (
        "profile", "input_list", "output_file", "run", "segment", "events",
        "first_entry", "config", "source_file_ordinal", "input_file_sha256",
        "source_manifest_sha256", "sample", "period", "si_di_role", "verbosity")]
    return str(MACRO) + "(" + ",".join(json.dumps(v, ensure_ascii=False) for v in values) + (",true" if check_only else "") + ")"


def key_values(obj):
    pairs = [line.split("=", 1) for line in str(obj).splitlines() if "=" in line]
    return unique_object(pairs)


def check_completion(job, uproot, library_hash):
    """Read the direct eventTree, source interval and producer terminal record."""
    import numpy as np
    format_modules()
    from eventtree.integrity import metadata as read_metadata
    from eventtree.schema import Schema
    with uproot.open(job["output_file"]) as product:
        completion = key_values(product["completion"])
        execution = key_values(product["metadata"])
        require(completion.get("completion_status") == "complete", "producer output is not complete")
        expected = {key: str(job[key]) for key in (
            "run", "segment", "source_file_ordinal", "first_entry", "source_manifest_sha256",
            "sample", "period", "si_di_role", "configuration_sha256", "macro_sha256")}
        expected.update(max_events=str(job["events"]), producer_library_sha256=library_hash)
        for key, value in expected.items():
            require(execution.get(key) == value, "output provenance mismatch: " + key)
        meta = read_metadata(product)
        schema = Schema()
        require(meta['schema']['schema_sha256'] == schema.sha256, 'native schema differs')
        nominal = {'pho_bdt_score', 'pho_bdt_state', 'pho_id_class', 'pho_iso_class'}
        require(not (set(meta['native_capture']['missing_branches']) - nominal), 'native mapping is incomplete')
        require(meta['source'].get('completed') is True, 'source is incomplete')
        for key in ('run', 'segment', 'source_file_ordinal', 'first_entry', 'sample', 'period', 'si_di_role', 'source_manifest_sha256'):
            require(meta['source'].get(key) == job[key], 'source metadata mismatch: ' + key)
        trees = [name for name, kind in product.classnames(recursive=True, cycle=False).items() if kind == 'TTree']
        require(trees == ['eventTree'], 'producer must write exactly one eventTree')
        tree = product['eventTree']
        require(set(schema.branches) <= set(tree.keys()), 'producer omitted schema branches')
        observed = int(completion['upstream_observed_events'])
        retained = int(completion['retained_events'])
        rejected = int(completion['upstream_rejected_events'])
        require(observed == job['events'], 'source ended before/after its declared event count')
        require(int(completion['encountered_events']) == retained and observed == retained + rejected,
                'producer completion accounting is inconsistent')
        require(tree.num_entries == observed == meta['schema']['entries'], 'saved event count differs')
        require(meta['schema']['retained_events'] == retained, 'retained metadata count differs')
        seen = kept = dropped = 0
        for batch in tree.iterate(['evt_record', 'evt_run', 'evt_src_file_ordinal', 'evt_src_entry'], step_size='8 MB', library='np'):
            n = len(batch['evt_record'])
            require(np.array_equal(batch['evt_src_entry'], np.arange(job['first_entry']+seen, job['first_entry']+seen+n)),
                    'source entries contain a gap, overlap or reordered event')
            require(bool(np.all(batch['evt_src_file_ordinal'] == job['source_file_ordinal'])), 'wrong source ordinal')
            require(bool(np.all(batch['evt_run'] == job['run'])), 'wrong event run')
            require(bool(np.isin(batch['evt_record'], [0, 1]).all()), 'unknown event record kind')
            kept += int(np.sum(batch['evt_record'] == 0)); dropped += int(np.sum(batch['evt_record'] == 1)); seen += n
        require((seen, kept, dropped) == (observed, retained, rejected), 'saved event accounting differs')
        return completion



def check_mapping(job, root, environment):
    """Query the actual loaded producer through the macro's read-only plan mode."""
    result = subprocess.run([root, '-l', '-b', '-q', macro_argument(job, check_only=True)],
                            cwd=TREE_DIR, env=environment, capture_output=True, text=True, timeout=120)
    prefix = 'TREEPRODUCTION_MAPPING_JSON='
    reports = [line[len(prefix):] for line in result.stdout.splitlines() if line.startswith(prefix)]
    require(len(reports) == 1, 'native mapping preflight did not return one report: ' + result.stderr[-1200:])
    report = json.loads(reports[0])
    require(report.get('ready') is True and result.returncode == 0,
            'native mapping is incomplete: ' + ', '.join(report.get('missing_branches', [])[:12]))
    require(report.get('schema_sha256') == sha256(TREE_DIR / 'config/eventTree_v1.yaml'), 'mapping schema differs')
    return report


def execute_piece(job, root, environment, uproot, library_hash):
    """The only per-source execution path for every command and restriction."""

    directory = Path(job["output_file"]).parent
    directory.mkdir()  # Exclusive reservation; partial outputs are never reused.
    work = directory / '.work'
    work.mkdir()
    log = directory / "production.log"


    require(sha256(Path(job["input_list"])) == job["input_list_sha256"], "input list changed after preflight")
    require(sha256(Path(job["config"])) == job["configuration_sha256"], "configuration changed after preflight")
    require(sha256(MACRO) == job["macro_sha256"], "macro changed after preflight")

    capture_job = dict(job, output_file=str(work / 'eventTree.native.root'))
    argv = [root, "-l", "-b", "-q", macro_argument(capture_job)]


    (directory / "invocation.json").write_text(json.dumps(dict(job=job, argv=argv), indent=2) + "\n")
    print(f"  source {job['source_file_ordinal']}: {job['events']} events; log {log}", flush=True)

    with log.open("x") as stream:
        result = subprocess.run(argv, cwd=TREE_DIR, env=environment, stdout=stream, stderr=subprocess.STDOUT, timeout=job.get("timeout", 43200))


    if result.returncode:
        print(f"ERROR: producer failed ({result.returncode}); see {log}", file=sys.stderr)
        with log.open(errors="replace") as stream:
            for line in deque(stream, maxlen=8):
                print("  " + line.rstrip(), file=sys.stderr)
        return result.returncode if result.returncode > 0 else 128 - result.returncode

    require(sha256(Path(job["input_list"])) == job["input_list_sha256"], "input list changed during execution")
    require(sha256(Path(job["config"])) == job["configuration_sha256"], "configuration changed during execution")
    require(sha256(MACRO) == job["macro_sha256"], "macro changed during execution")
    if job.get("producer_library"):
        require(sha256(Path(job["producer_library"])) == library_hash, "installed producer library changed during execution")

    try:
        completion = check_completion(capture_job, uproot, library_hash)
    except Exception as error:
        raise FrontendError(f"completion check failed: {error}; output retained for inspection: {directory}") from error

    format_modules()
    from eventtree.nominal import finalize_file
    native = Path(capture_job['output_file'])
    before = sha256(native)
    final = work / 'eventTree.final.root'
    report = finalize_file(native, final, job.get('isolation_freeze'))
    require(sha256(native) == before, 'native tree changed before publication')
    if job.get('isolation_freeze'):
        require(sha256(Path(job['isolation_freeze'])) == job['isolation_freeze_sha256'], 'isolation freeze changed before publication')
    os.link(final, job['output_file'])  # Atomic no-replace publication, same filesystem.
    (directory / 'finalization.json').write_text(json.dumps(dict(status='FINAL_EVENTTREE_CHECKED',
        output_sha256=sha256(Path(job['output_file'])), checks=report, replacement_certified=False), indent=2)+'\n')
    return 0


def execute_source(job, root, environment, uproot, library_hash):
    """Process bounded pieces, then publish exactly one file for this source."""
    format_modules()
    from eventtree.consolidate import consolidate
    output = Path(job['output_file'])
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent / '.work' / output.stem
    work.mkdir(parents=True, exist_ok=job.get('resume', False))
    with (work / 'source.lock').open('a') as lock:
        try: fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close()
            raise FrontendError('source is already owned by another invocation: ' + str(output))
        identity = {k:v for k,v in job.items() if k != 'resume'}
        identity['producer_library_sha256'] = library_hash
        identity['runtime'] = {str(p.relative_to(TREE_DIR)):sha256(p)
                               for p in sorted((TREE_DIR/'scripts').rglob('*'))
                               if p.is_file() and p.suffix in ('.py','.hxx')}
        identity_pin = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        identity_file = work / 'source_identity.json'
        if identity_file.exists():
            require(read_json(identity_file)['sha256'] == identity_pin, 'resume source or runtime identity differs')
        else:
            with identity_file.open('x') as f: json.dump(dict(sha256=identity_pin, job=identity), f, indent=2)
        receipt_path = work / 'source_receipt.json'
        if receipt_path.exists():
            receipt = read_json(receipt_path)
            require(receipt['identity_sha256'] == identity_pin, 'source receipt identity differs')
            require(output.is_file() and sha256(output) == receipt['output_sha256'], 'published source changed')
            lock.close()
            return 0
        limit = PIECE_EVENTS.get(job['profile'], PIECE_EVENTS['default'])
        chunks = []
        total = job['events']
        ranges = [(i, min(i + limit, total)) for i in range(0, total, limit)] or [(0, 0)]
        started = time.monotonic()
        def progress(phase, processed, **extra):
            record = dict(phase=phase, logical_source_id=job['logical_source_id'],
                          processed=processed, total=total, remaining=total-processed,
                          elapsed_seconds=time.monotonic()-started,
                          updated_at=dt.datetime.now(dt.timezone.utc).isoformat(), **extra)
            pending = work / 'progress.pending.json'
            pending.write_text(json.dumps(record, sort_keys=True) + '\n')
            os.replace(pending, work / 'progress.json')
        progress('preparing', 0)
        try:
            for index, (begin, end) in enumerate(ranges):
                piece = dict(job, output_file=str(work / f'piece-{index:06d}' / 'trees.root'),
                             first_entry=job['first_entry']+begin, events=end-begin)
                progress('reconstruction', begin, piece=index, pieces=len(ranges))
                piece_path = Path(piece['output_file'])
                if piece_path.parent.exists():
                    require(job.get('resume', False), 'existing piece requires explicit resume')
                    saved_receipt = piece_path.parent / 'finalization.json'
                    require(saved_receipt.is_file() and piece_path.is_file(),
                            'partial piece has no validated receipt; retained for inspection: ' + str(piece_path.parent))
                    saved = read_json(saved_receipt)
                    require(saved['status'] == 'FINAL_EVENTTREE_CHECKED' and saved['output_sha256'] == sha256(piece_path),
                            'piece receipt or saved output changed')
                    check_completion(piece, uproot, library_hash)
                else:
                    code = execute_piece(piece, root, environment, uproot, library_hash)
                    require(code == 0, f'piece {index} failed with exit {code}')
                chunks.append(dict(path=piece['output_file'], sha256=sha256(Path(piece['output_file'])),
                                   piece_id=str(index), original_unit_id=job['logical_source_id'],
                                   event_begin=piece['first_entry'], event_end=piece['first_entry']+piece['events']))
                progress('pieces_complete', end, piece=index, pieces=len(ranges))
            candidate = work / 'source.root'
            proof_path = work / 'source_validation.json'
            if proof_path.exists():
                saved = read_json(proof_path)
                require(saved['identity_sha256'] == identity_pin and sha256(candidate) == saved['sha256'], 'saved source candidate differs')
                proof = saved['proof']
            elif job['first_entry'] != 0:
                # A selected interval is a development product, never a full-source claim.
                require(len(chunks) == 1, 'partial-source development output must fit one piece')
                os.link(chunks[0]['path'], candidate)
                proof = dict(status='DEVELOPMENT_INTERVAL', first_entry=job['first_entry'], events=total)
            else:
                proof = consolidate(chunks, candidate, logical_source_id=job['logical_source_id'],
                                    expected_events=total, progress=lambda done,n: progress('source_validation', done))
            if not proof_path.exists():
                with proof_path.open('x') as f:
                    json.dump(dict(identity_sha256=identity_pin, sha256=sha256(candidate), proof=proof), f, indent=2)
            if output.exists():
                require(job.get('resume', False) and sha256(output) == sha256(candidate), 'source publication collision')
            else: os.link(candidate, output)
            receipt = dict(status='SOURCE_COMPLETE' if job['complete_source'] else 'DEVELOPMENT_INTERVAL',
                           identity_sha256=identity_pin, logical_source_id=job['logical_source_id'], output_sha256=sha256(output),
                           input_list_sha256=job['input_list_sha256'], source_manifest_sha256=job['source_manifest_sha256'],
                           pieces=len(chunks), entries=total, proof=proof, replacement_certified=False)
            with receipt_path.open('x') as f: json.dump(receipt, f, indent=2)
            progress('complete', total)
            lock.close()
            return 0
        except Exception as error:
            progress('failed', sum(c['event_end']-c['event_begin'] for c in chunks), error=str(error))
            lock.close()
            raise


def main():
    global BINDINGS
    parser = argparse.ArgumentParser(prog="produce_trees.sh", add_help=False)
    parser.add_argument("--events", type=int)
    parser.add_argument("--source", type=int)
    parser.add_argument("--run", type=int)
    parser.add_argument("--sample")
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--timeout", type=int, default=43200)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("lanes", nargs="+")
    args = parser.parse_args()
    integer(args.timeout, "--timeout", minimum=1, maximum=43200)
    if args.catalogue is not None:
        BINDINGS = args.catalogue.resolve()
    if args.run is not None:
        integer(args.run, "--run", minimum=1)
        require(all(lane in ('pp_data', 'auau_data') for lane in args.lanes), "--run requires a DATA command")
    if args.sample is not None:
        require(all(lane not in ('pp_data', 'auau_data') for lane in args.lanes), "--sample requires a SIM command")


    if args.events is not None:
        integer(args.events, "--events", minimum=1)
    if args.source is not None:
        integer(args.source, "--source", maximum=2**63 - 1)
        require(len(args.lanes) == 1, "--source requires a primitive command")


    settings = read_json(BINDINGS)
    known_keys(settings, "schema_version campaign output_root production_config producer_library lanes "
                        "isolation_freeze isolation_freeze_sha256", "binding")
    require(settings.get("schema_version") == 1, "unsupported population bindings schema")
    bindings = settings.get("lanes")
    require(isinstance(bindings, dict), "missing lanes in population bindings")

    missing = []
    for lane in args.lanes:
        require(lane in bindings and isinstance(bindings[lane], dict), f"unknown lane binding: {lane}")
        known_keys(bindings[lane], "profile manifest manifest_sha256", f"{lane} binding")
        if not bindings[lane].get("manifest"):
            missing.append(f"{lane}: canonical population manifest is not yet configured.\n"
                           f"  Expected binding: {BINDINGS} -> lanes.{lane}.manifest")
    require(not missing, "\n".join(missing))

    campaign = text(settings.get("campaign"), f"{BINDINGS} -> campaign")
    require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", campaign), "campaign must be a safe directory name")
    require(MACRO.is_file(), f"production macro is missing: {MACRO}")


    plans = [plan_population(lane, bindings[lane], settings, args) for lane in args.lanes]


    for plan in plans:
        print(f"Lane: {plan['lane']} | profile: {plan['profile']}")
        print(f"  Population: {plan['manifest']}")
        for config in sorted({job['config'] for job in plan['jobs']}):
            print(f"  Configuration: {config}")
        print(f"  Sources: {len(plan['jobs'])}/{plan['population_sources']} | "
              f"{'DEVELOPMENT RESTRICTION' if plan['development'] else 'declared full population'}")
        print(f"  Output: {plan['output']}/<run-or-sample>/part-<ordinal>.root")


    if args.dry_run:
        print("DRY RUN: bindings checked; no outputs created and no ROOT invoked.\n"
              "C++ configuration, conditions, environment and science readiness were not validated.")
        return 0

    if any(plan['lane'].startswith('auau') for plan in plans):
        isolation = pinned_file(BINDINGS.parent, settings.get('isolation_freeze'),
                               settings.get('isolation_freeze_sha256'), 'AuAu isolation freeze')
        for plan in plans:
            if plan['lane'].startswith('auau'):
                for job in plan['jobs']:
                    job.update(isolation_freeze=str(isolation), isolation_freeze_sha256=settings['isolation_freeze_sha256'])


    root = shutil.which("root")
    require(root, "ROOT is not on PATH; enter the configured sPHENIX environment")

    library_value = settings.get("producer_library")
    if library_value:
        library = (BINDINGS.parent / text(library_value, "producer_library")).absolute()
    else:
        prefix = os.environ.get("MYINSTALL")
        require(prefix, "set MYINSTALL to the TreeProduction installation prefix")
        library = Path(prefix).absolute() / 'lib/libPhotonJetTree.so'

    # Keep the installed .so symlink name: versioned targets are normal, while
    # R__LOAD_LIBRARY still resolves libPhotonJetTree.so in its parent directory.
    require(library.is_file() and library.name == "libPhotonJetTree.so",
            f"installed producer library is missing or misnamed: {library}")

    try:
        import uproot
    except ImportError as error:
        raise FrontendError("uproot is required to verify the producer's completion record") from error


    environment = os.environ.copy()
    include = library.parent.parent / 'include'
    require((include / 'TreeProduction/PhotonJetTree.h').is_file(), "installed TreeProduction headers are missing")
    environment['ROOT_INCLUDE_PATH'] = str(include) + (os.pathsep + environment['ROOT_INCLUDE_PATH']
                                                      if environment.get('ROOT_INCLUDE_PATH') else '')
    for variable in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        environment[variable] = str(library.parent) + (os.pathsep + environment[variable] if environment.get(variable) else "")

    library_hash = sha256(library)


    # The installed native implementation declares mapping readiness. Check
    # each resolved plan before creating an output directory or processing DSTs.
    for plan in plans:
        for job in plan['jobs']:
            limit = PIECE_EVENTS.get(job['profile'], PIECE_EVENTS['default'])
            check_mapping(dict(job, events=min(job['events'], limit)), root, environment)

    # Admit the whole aggregate before starting its first event loop.
    for plan in plans:
        directory = Path(plan["output"])
        directory.mkdir(parents=True, exist_ok=args.resume)
        manifest_copy = directory / 'population.json'
        if manifest_copy.exists():
            require(sha256(manifest_copy) == plan['manifest_sha256'], 'resume population differs')
        else:
            with manifest_copy.open('xb') as f: f.write(Path(plan['manifest']).read_bytes())


    for plan in plans:
        require(sha256(Path(plan["manifest"])) == plan["manifest_sha256"], "population changed after preflight")
        for job in plan["jobs"]:
            require(sha256(library) == library_hash, "installed producer library changed during this invocation")
            job["producer_library"] = str(library)
            status = execute_source(job, root, environment, uproot, library_hash)
            if status:
                return status

        label = "development selection" if plan["development"] else "declared population"
        print(f"FINAL TREE COMPLETE: {plan['lane']} {label}; {len(plan['jobs'])} saved outputs checked.\n"
              "Independent replacement certification is a separate release gate.")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FrontendError, OSError, ValueError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("ERROR: interrupted; partial output remains unaccepted.", file=sys.stderr)
        sys.exit(130)
