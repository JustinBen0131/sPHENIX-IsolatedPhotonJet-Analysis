"""Run maintained DST populations through the direct eventTree producer.

The C++ producer reconstructs, labels and writes every eventTree. This driver
resolves source identities, launches Fun4All in bounded pieces, and calls the
installed ROOT-only eventTree library to check each saved piece and assemble
one proven file per logical source before collision-safe publication.
"""
import argparse
import csv
from collections import deque
from contextlib import contextmanager
import datetime as dt
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

TREE_DIR = Path(__file__).resolve().parents[1]
BINDINGS = TREE_DIR / "config" / "production_populations.json"
MACRO = TREE_DIR / "macros" / "Fun4All_PhotonJetTree.C"
SCHEMA = TREE_DIR / "config" / "eventTree_v1.yaml"
PHOTON_ID = TREE_DIR.parent / "PhotonID" / "models" / "nominal.json"
PIECE_EVENTS = {"auau_data": 35000, "default": 100000}
SOURCE_INTEGERS = ('source_file_ordinal', 'source_global_entry_begin', 'run', 'segment', 'first_entry', 'event_count')
SOURCE_TEXT = ('sample', 'period', 'si_di_role', 'logical_source_id', 'input_list', 'input_list_sha256')
SOURCE_COLUMNS = SOURCE_INTEGERS + SOURCE_TEXT + tuple(f'rng_seed_{i}' for i in range(1, 6)) + (
    'rng_pedestal_sequence', 'cemc_map_path', 'cemc_map_sha256')
INDEX_COLUMNS = ('group', 'manifest', 'manifest_sha256', 'sources', 'events')
DATA_STREAMS = ('DST_JETCALO.list', 'DST_Jet.list')
SIM_STREAMS = ('DST_CALO_CLUSTER.list', 'G4Hits.list', 'DST_TRUTH_JET.list', 'DST_GLOBAL.list', 'DST_MBD_EPD.list')
EMBEDDED_STREAMS = ('DST_CALO.list', 'DST_TRUTH_G4HIT.list', 'DST_TRUTH_JET.list', 'DST_GLOBAL.list', 'DST_MBD_EPD.list')
DI_STREAMS = ('G4Hits.list', 'DST_TRUTH_JET.list')


def load_eventtree(library_dir):
    """Bind the installed ROOT-only eventTree library; no detector library is loaded."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    include = Path(library_dir).parent / "include"
    require((include / "TreeProduction/EventTree.h").is_file(),
            "installed eventTree header is missing: " + str(include))
    ROOT.gInterpreter.AddIncludePath(str(include))
    require(ROOT.gSystem.Load(str(Path(library_dir) / "libEventTree")) >= 0,
            "cannot load libEventTree from " + str(library_dir))
    require(ROOT.gInterpreter.Declare('#include "TreeProduction/EventTree.h"'),
            "cannot declare the eventTree interface")
    return ROOT.eventtree


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
    raw = gzip.decompress(path.read_bytes()).decode('utf-8') if path.suffix == '.gz' else path.read_text()
    result = json.loads(raw, object_pairs_hook=unique_object)
    require(isinstance(result, dict), f"expected a JSON object: {path}")
    return result


def tsv_rows(path, columns):
    """Read a fixed-column source/index table; reject malformed or extra fields."""
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t', strict=True)
        require(reader.fieldnames == list(columns), 'unexpected TSV columns: ' + str(path))
        rows = list(reader)
    require(rows and all(set(row) == set(columns) and all(v is not None for v in row.values()) for row in rows),
            'empty or malformed TSV: ' + str(path))
    return rows


def stream_lists(lane):
    """Declared DST stream names in the order used to align source rows."""
    if lane.endswith('_data'):
        return DATA_STREAMS
    if lane.endswith('_di'):
        return DI_STREAMS
    return EMBEDDED_STREAMS if lane.endswith('_embedded') else SIM_STREAMS


def read_population(path, lane, campaign):
    """Join aligned conventional DST lists with their pinned source metadata."""
    if path.suffix != '.tsv':
        return read_json(path)  # Existing custom JSON catalogues remain readable.
    rows = tsv_rows(path, SOURCE_COLUMNS)
    streams = [(path.parent / name).read_text().splitlines() for name in stream_lists(lane)]
    require(all(len(values) == len(rows) and all(values) for values in streams),
            'DST stream lists must have one aligned row per source: ' + str(path.parent))
    sources, bundles = [], {}
    for index, row in enumerate(rows):
        row = {key: '' if value == '-' else value for key, value in row.items()}
        source = {key: int(row[key]) for key in SOURCE_INTEGERS if row[key] != ''}
        source.update({key: row[key] for key in SOURCE_TEXT})
        rng_fields = [row[f'rng_seed_{i}'] for i in range(1, 6)] + [row['rng_pedestal_sequence']]
        require(all(rng_fields) or not any(rng_fields), 'incomplete archived DI seed record')
        if any(rng_fields):
            source['historical_rng'] = dict(seeds=list(map(int, rng_fields[:5])), pedestal_sequence=int(rng_fields[5]))
        require(bool(row['cemc_map_path']) == bool(row['cemc_map_sha256']), 'incomplete CEMC map binding')
        if row['cemc_map_path']:
            source['cemc_recovery_map'] = dict(path=row['cemc_map_path'], sha256=row['cemc_map_sha256'])
        name = source['input_list']
        require(name not in bundles, 'duplicate source bundle label: ' + name)
        bundles[name] = '\t'.join(values[index] for values in streams) + '\n'
        sources.append(source)
    return dict(schema_version=2, lane=lane, campaign=campaign, sources=sources, input_lists=bundles)


def selection_list(path, numeric=False):
    """A conventional run/sample list: one identifier per line, optional # comments."""
    values = []
    for line in path.read_text().splitlines():
        value = line.split('#', 1)[0].strip()
        if not value:
            continue
        if numeric:
            require(re.fullmatch(r'[0-9]+', value) is not None and int(value) > 0,
                    'run list requires one positive run number per line: ' + str(path))
            value = str(int(value))
        else:
            require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value) is not None,
                    'sample list requires one exact sample name per line: ' + str(path))
        values.append(value)
    require(values and len(values) == len(set(values)), 'selection list is empty or contains duplicates: ' + str(path))
    return set(values)


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


def source_group(lane, source):
    """The run (DATA) or exact sample (SIM) that groups a lane's sources."""
    group = str(source['run']) if lane.endswith('_data') else source['sample']
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_. -]*', group), 'unsafe source group name')
    return group


def plan_population(lane, binding, settings, args):
    """Resolve job metadata only; never interpret DST columns or YAML policy."""


    manifest = pinned_file(BINDINGS.parent, binding.get("manifest"),
                           binding.get("manifest_sha256"), f"{lane} manifest")
    population = read_population(manifest, lane, settings['campaign'])
    known_keys(population, "schema_version lane campaign sources input_lists", "manifest")
    require(population.get("schema_version") in (1, 2), f"unsupported manifest schema: {manifest}")
    require(population.get("lane") == lane, f"manifest lane mismatch: {manifest}")
    require(population.get("campaign") == settings["campaign"], f"manifest campaign mismatch: {manifest}")
    sources = population.get("sources")
    require(isinstance(sources, list) and sources, f"empty/missing source population: {manifest}")
    embedded = population.get('input_lists')
    require((population['schema_version'] == 2) == isinstance(embedded, dict),
            'schema 2 requires embedded input lists; schema 1 uses separate files')
    if embedded is not None:
        require(set(embedded) == {source['input_list'] for source in sources},
                'embedded input-list roster differs from source roster')


    profile = text(binding.get("profile"), f"lanes.{lane}.profile")
    config = (BINDINGS.parent / text(settings.get("production_config"), "production_config")).resolve()
    require(config.is_file(), f"production configuration is missing: {config}")


    jobs, ordinals, bundles, rosters = [], set(), set(), {}
    config_hash, macro_hash = sha256(config), sha256(MACRO)
    for source in sources:
        require(isinstance(source, dict), f"{lane}: source must be an object")
        known_keys(source, "source_file_ordinal source_global_entry_begin input_list input_list_sha256 "
                          "run segment first_entry event_count sample period si_di_role production_config "
                          "production_config_sha256 logical_source_id input_file_sha256 historical_rng "
                          "cemc_recovery_map", "source")
        require(("production_config" in source) == ("production_config_sha256" in source),
                "source production_config and production_config_sha256 must be supplied together")
        ordinal = integer(source.get("source_file_ordinal"), "source_file_ordinal", maximum=2**63 - 1)
        # Global entry of this file's entry 0 in the maintained source order;
        # -1 when the catalogue does not define one.
        global_begin = (integer(source["source_global_entry_begin"], "source_global_entry_begin", maximum=2**63 - 1)
                        if "source_global_entry_begin" in source else -1)
        first = integer(source.get("first_entry"), "first_entry")
        events = integer(source.get("event_count"), "event_count", minimum=0)
        if embedded is None:
            inputs = pinned_file(manifest.parent, source.get("input_list"),
                                 source.get("input_list_sha256"), "input_list")
        else:
            name = text(source.get('input_list'), 'input_list')
            require(Path(name).name == name, 'embedded input list must be a filename')
            content = embedded[name]
            require(isinstance(content, str) and
                    hashlib.sha256(content.encode('utf-8')).hexdigest() == source.get('input_list_sha256'),
                    'embedded input-list SHA256 mismatch: ' + name)
            inputs = str(manifest) + '#' + name


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
                   source_global_entry_begin=global_begin,
                   input_file_sha256=source.get("input_file_sha256", ""),
                   source_manifest_sha256=binding["manifest_sha256"],
                   sample=text(source.get("sample"), "sample"),
                   period=text(source.get("period"), "period", allow_empty=True),
                   si_di_role=text(source.get("si_di_role"), "si_di_role", allow_empty=True),
                   input_list_sha256=source["input_list_sha256"],
                   historical_rng=historical_rng(source.get("historical_rng")),
                   **cemc_recovery_map(source.get("cemc_recovery_map")),
                   configuration_sha256=(config_hash if job_config == config else sha256(job_config)),
                   macro_sha256=macro_hash, verbosity=0)
        if embedded is not None:
            job['input_list_content'] = content
            job['input_list_identity'] = str(inputs)
        job['logical_source_id'] = source.get('logical_source_id') or (
            lane + ':' + str(job['run']) + ':' + str(job['segment']) + ':' + str(ordinal))
        text(job['logical_source_id'], 'logical_source_id')
        # Ordinals number the files of one group (DATA run or SIM sample).
        group = source_group(lane, job)
        require((group, ordinal) not in ordinals, f"{lane}: duplicate source ordinal {ordinal} in {group}")
        ordinals.add((group, ordinal))
        rosters.setdefault(group, []).append(job['logical_source_id'])
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
    require(args.source is None or len({source_group(lane, job) for job in jobs}) == 1,
            f"{lane}: --source {args.source} exists in several groups; add --run or --sample")


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
    require(not output.exists() or getattr(args, "resume", False) or getattr(args, 'output_admitted', False),
            f"output already exists; use --resume for exact verified work: {output}")

    for job in jobs:
        job['output_file'] = str(output / source_group(lane, job) / f"part-{job['source_file_ordinal']:08d}.root")

    return dict(lane=lane, campaign=settings["campaign"], profile=profile, manifest=str(manifest),
                manifest_sha256=binding["manifest_sha256"], output=str(output),
                development=development, population_sources=len(sources), jobs=jobs,
                population_events=sum(source['event_count'] for source in sources),
                rosters={group: sorted(ids) for group, ids in rosters.items()})


def population_plans(settings, args):
    """Read one run/sample at a time, keeping full-campaign memory bounded."""
    routes = []
    for lane in args.lanes:
        binding = settings['lanes'][lane]
        if 'index' in binding:
            index = pinned_file(BINDINGS.parent, binding['index'], binding.get('index_sha256'), lane + ' index')
            groups = []
            for row in tsv_rows(index, INDEX_COLUMNS):
                groups.append(dict(group=row['group'], manifest=str((index.parent / row['manifest']).resolve()),
                                   manifest_sha256=row['manifest_sha256'], sources=int(row['sources']), events=int(row['events'])))
        else:
            groups = binding.get('groups')
        routes.append((lane, binding, groups))
    selected_groups = getattr(args, 'selected_groups', None)
    if selected_groups is not None:
        require(all(groups is not None for _, _, groups in routes), 'list selectors require indexed populations')
        known = {group['group'] for _, _, groups in routes for group in groups}
        require(not selected_groups - known, 'unknown run/sample in selection list: ' + ', '.join(sorted(selected_groups - known)))
    for lane, binding, groups in routes:
        if groups is None:
            yield plan_population(lane, binding, settings, args)
            continue
        require(isinstance(groups, list) and groups, lane + ': missing groups')
        names = set()
        selected = 0
        require(args.source is None or args.run is not None or args.sample is not None
                or (selected_groups is not None and len(selected_groups) == 1),
                '--source requires --run or --sample for a grouped catalogue')
        for group in groups:
            known_keys(group, 'group manifest manifest_sha256 sources events', 'population group')
            name = text(group.get('group'), 'population group')
            require(name not in names, lane + ': duplicate population group ' + name)
            names.add(name)
            if selected_groups is not None and name not in selected_groups:
                continue
            if args.run is not None and name != str(args.run):
                continue
            if args.sample is not None and name != args.sample:
                continue
            plan = plan_population(lane, dict(profile=binding['profile'],
                                             manifest=group['manifest'],
                                             manifest_sha256=group['manifest_sha256']), settings, args)
            require(set(plan['rosters']) == {name}, 'manifest group differs from index: ' + name)
            require(plan['population_sources'] == group['sources'], 'source count differs from index')
            # Restricted plans still bind the complete run/sample manifest.
            require(plan['population_events'] == group['events'],
                    'event count differs from index')
            plan['group'] = name
            selected += 1
            yield plan
        require(selected or selected_groups is not None, lane + ': requested run/sample does not exist')


@contextmanager
def input_bundle(job):
    """Materialize only the current source's exact DST list for Fun4All."""
    if 'input_list_content' not in job:
        yield job
        return
    with tempfile.TemporaryDirectory(prefix='treeproduction-input-') as directory:
        path = Path(directory) / 'source.list'
        path.write_bytes(job['input_list_content'].encode('utf-8'))
        require(sha256(path) == job['input_list_sha256'], 'materialized input list differs')
        yield dict(job, input_list=str(path))


def historical_rng(record):
    """An archived double-interaction source's replayed seeds as the producer reads
    them, "s1,s2,s3,s4,s5:pedestal"; empty for every other source."""
    if record is None:
        return ""
    require(isinstance(record, dict) and set(record) == {"seeds", "pedestal_sequence"},
            "historical_rng must hold exactly seeds and pedestal_sequence")
    seeds = record["seeds"]
    require(isinstance(seeds, list) and len(seeds) == 5, "historical_rng needs five seeds")
    values = [integer(seed, "historical_rng seed", maximum=2**32 - 1) for seed in seeds]
    pedestal = integer(record["pedestal_sequence"], "historical_rng pedestal_sequence", maximum=3259)
    return ",".join(map(str, values)) + ":" + str(pedestal)


def cemc_recovery_map(record):
    """An Au+Au DATA source's own CEMC bad-tower map, {"path": ..., "sha256": ...}; empty otherwise."""
    if record is None:
        return dict(cemc_recovery_map="", cemc_recovery_map_sha256="")
    require(isinstance(record, dict) and set(record) == {"path", "sha256"},
            "cemc_recovery_map must hold exactly path and sha256")
    path, digest = text(record["path"], "cemc_recovery_map path"), text(record["sha256"], "cemc_recovery_map sha256")
    require(path.startswith("/") and re.fullmatch(r"[0-9a-f]{64}", digest),
            "cemc_recovery_map needs an absolute path and a SHA256")
    return dict(cemc_recovery_map=path, cemc_recovery_map_sha256=digest)


def macro_argument(job, check_only=False):
    """Serialize one job in the steering macro's declared argument order."""
    # JSON string quoting is also valid for these C++ UTF-8 string literals.
    # subprocess receives an argument vector: no shell, eval or interpolation.
    values = [job[key] for key in (
        "profile", "input_list", "output_file", "run", "segment", "events",
        "first_entry", "config", "source_file_ordinal", "input_file_sha256",
        "source_manifest_sha256", "sample", "period", "si_di_role", "verbosity")]
    suffix = (",true" if check_only else ",false")
    suffix += "," + json.dumps(job.get('resolved_configuration_file', ''), ensure_ascii=False)
    suffix += "," + json.dumps(job.get('source_global_entry_begin', -1))
    suffix += "," + json.dumps(job.get('historical_rng', ''))
    suffix += "," + json.dumps(job.get('cemc_recovery_map', ''), ensure_ascii=False)
    suffix += "," + json.dumps(job.get('cemc_recovery_map_sha256', ''))
    return str(MACRO) + "(" + ",".join(json.dumps(v, ensure_ascii=False) for v in values) + suffix + ")"


def configuration_document(raw):
    """Parse configuration syntax once; all physics decisions remain native."""
    import yaml

    class ConfigurationLoader(yaml.SafeLoader):
        pass

    def mapping(loader, node):
        loader.flatten_mapping(node)
        values = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=True)
            require(type(key) in (str, int), 'configuration keys must be strings or integers')
            key = str(key)
            require(key not in values, 'duplicate configuration key: ' + key)
            values[key] = loader.construct_object(value_node, deep=True)
        return values

    ConfigurationLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, mapping)
    document = yaml.load(raw, Loader=ConfigurationLoader)
    require(isinstance(document, dict), 'configuration must be an object')
    # Refuse YAML-specific objects and non-finite numbers rather than changing
    # their meaning while crossing the JSON boundary.
    json.dumps(document, allow_nan=False)
    return document


@contextmanager
def resolved_configuration(job, directory=None):
    """Supply parsed configuration to the native plan; physics values remain native-owned."""
    config = Path(job['config'])
    raw = config.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == job['configuration_sha256'],
            'configuration changed before parsing')
    envelope = dict(schema='TreeProductionConfigurationV1',
                    source_sha256=job['configuration_sha256'], document=configuration_document(raw))
    content = json.dumps(envelope, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n'
    with tempfile.TemporaryDirectory(prefix='treeproduction-config-') as temporary:
        path = (Path(temporary) if directory is None else Path(directory)) / 'configuration.json'
        with path.open('x') as stream:
            stream.write(content)
        yield dict(job, resolved_configuration_file=str(path),
                   resolved_configuration_sha256=sha256(path))
        require(path.read_text() == content and sha256(config) == job['configuration_sha256'],
                'configuration changed during native execution')


def key_values(obj):
    pairs = [line.split("=", 1) for line in str(obj).splitlines() if "=" in line]
    return unique_object(pairs)


def check_completion(job, native, library_hash):
    """The producer's completion record, source interval and saved structure, read natively."""
    expected = {key: job[key] for key in (
        "run", "segment", "source_file_ordinal", "first_entry", "source_manifest_sha256",
        "sample", "period", "si_di_role", "configuration_sha256", "macro_sha256", "events")}
    expected["producer_library_sha256"] = library_hash
    completion = json.loads(str(native.checkProducerOutput(job["output_file"], str(SCHEMA), json.dumps(expected))))
    structure = json.loads(str(native.checkFile(job["output_file"], str(SCHEMA), str(PHOTON_ID))))
    require(structure["ok"], "saved eventTree structure: " + "; ".join(structure["problems"][:3]))
    return dict(completion=completion, structure=structure)


def consolidate_source(native, chunks, candidate, job, progress_path):
    """One proven file for the logical source: bit-preserved payloads, per-piece metadata and objects."""
    request = dict(chunks=chunks, output=str(candidate), logical_source_id=job["logical_source_id"],
                   expected_events=job["events"], schema_path=str(SCHEMA), photon_id_package=str(PHOTON_ID),
                   progress_path=str(progress_path))
    return json.loads(str(native.consolidate(json.dumps(request))))



def check_mapping(job, root, environment):
    """Query the actual loaded producer through the macro's read-only plan mode."""
    with resolved_configuration(job) as parsed_job:
        result = subprocess.run([root, '-l', '-b', '-q', macro_argument(parsed_job, check_only=True)],
                                cwd=TREE_DIR, env=environment, capture_output=True, text=True, timeout=120)
    prefix = 'TREEPRODUCTION_MAPPING_JSON='
    reports = [line[len(prefix):] for line in result.stdout.splitlines() if line.startswith(prefix)]
    require(len(reports) == 1, 'native mapping preflight did not return one report: ' + result.stderr[-1200:])
    report = json.loads(reports[0])
    require(report.get('ready') is True and result.returncode == 0,
            'native mapping is incomplete: ' + ', '.join(report.get('missing_branches', [])[:12]))
    require(report.get('schema_sha256') == sha256(TREE_DIR / 'config/eventTree_v1.yaml'), 'mapping schema differs')
    return report


def execute_piece(job, root, environment, native, library_hash):
    """Run one bounded reconstruction piece and check its native completion record."""

    directory = Path(job["output_file"]).parent
    directory.mkdir()  # Exclusive reservation; partial outputs are never reused.
    work = directory / '.work'
    work.mkdir()
    log = directory / "production.log"


    require(sha256(Path(job["input_list"])) == job["input_list_sha256"], "input list changed after preflight")
    require(sha256(Path(job["config"])) == job["configuration_sha256"], "configuration changed after preflight")
    require(sha256(MACRO) == job["macro_sha256"], "macro changed after preflight")

    with resolved_configuration(job, work) as parsed_job:
        capture_job = dict(parsed_job, output_file=str(work / 'eventTree.root'))
        argv = [root, "-l", "-b", "-q", macro_argument(capture_job)]
        (directory / "invocation.json").write_text(json.dumps(dict(job=capture_job, argv=argv), indent=2) + "\n")
        print(f"  source {job['source_file_ordinal']}: {job['events']} events; log {log}", flush=True)
        with log.open("x") as stream:
            result = subprocess.run(argv, cwd=TREE_DIR, env=environment, stdout=stream, stderr=subprocess.STDOUT,
                                    timeout=job.get("timeout", 43200))


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
        require(sha256(Path(job["producer_library"])) == library_hash,
                "installed producer library changed during execution")

    try:
        checks = check_completion(capture_job, native, library_hash)
    except Exception as error:
        raise FrontendError(f"completion check failed: {error}; output retained for inspection: {directory}") from error

    # The checked native file is the final piece: no rewrite, atomic no-replace
    # publication on the same filesystem.
    native_file = Path(capture_job['output_file'])
    os.link(native_file, job['output_file'])
    (directory / 'finalization.json').write_text(json.dumps(dict(status='FINAL_EVENTTREE_CHECKED',
        output_sha256=sha256(Path(job['output_file'])), checks=checks), indent=2)+'\n')
    return 0


def execute_source(job, root, environment, native, library_hash):
    """Process bounded pieces and prove one assembled logical-source file."""
    output = Path(job['output_file'])
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent / '.work' / output.stem
    work.mkdir(parents=True, exist_ok=job.get('resume', False))
    with (work / 'source.lock').open('a') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close()
            raise FrontendError('source is already owned by another invocation: ' + str(output))
        identity = {k:v for k,v in job.items() if k not in ('resume', 'input_list_content')}
        if 'input_list_identity' in job:
            # Temporary filenames are locators; the pinned catalogue and exact
            # list digest bind scientific identity and remain stable on resume.
            identity['input_list'] = job['input_list_identity']
        identity['producer_library_sha256'] = library_hash
        identity['runtime'] = {str(p.relative_to(TREE_DIR)):sha256(p)
                               for p in sorted((TREE_DIR/'scripts').rglob('*.py')) if p.is_file()}
        identity_pin = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        identity_file = work / 'source_identity.json'
        if identity_file.exists():
            require(read_json(identity_file)['sha256'] == identity_pin, 'resume source or runtime identity differs')
        else:
            with identity_file.open('x') as f:
                json.dump(dict(sha256=identity_pin, job=identity), f, indent=2)
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
                            'partial piece has no validated receipt; retained for inspection: '
                            + str(piece_path.parent))
                    saved = read_json(saved_receipt)
                    require(saved['status'] == 'FINAL_EVENTTREE_CHECKED'
                            and saved['output_sha256'] == sha256(piece_path),
                            'piece receipt or saved output changed')
                    check_completion(piece, native, library_hash)
                else:
                    code = execute_piece(piece, root, environment, native, library_hash)
                    require(code == 0, f'piece {index} failed with exit {code}')
                chunks.append(dict(path=piece['output_file'], sha256=sha256(Path(piece['output_file'])),
                                   piece_id=str(index), original_unit_id=job['logical_source_id'],
                                   event_begin=piece['first_entry'], event_end=piece['first_entry']+piece['events']))
                progress('pieces_complete', end, piece=index, pieces=len(ranges))
            candidate = work / 'source.root'
            proof_path = work / 'source_validation.json'
            if proof_path.exists():
                saved = read_json(proof_path)
                require(saved['identity_sha256'] == identity_pin and sha256(candidate) == saved['sha256'],
                        'saved source candidate differs')
                proof = saved['proof']
            elif job['first_entry'] != 0:
                # A selected interval is a development product, never a full-source claim.
                require(len(chunks) == 1, 'partial-source development output must fit one piece')
                os.link(chunks[0]['path'], candidate)
                proof = dict(status='DEVELOPMENT_INTERVAL', first_entry=job['first_entry'], events=total)
            else:
                progress('source_validation', 0)
                proof = consolidate_source(native, chunks, candidate, job, work / 'progress.json')
            if not proof_path.exists():
                with proof_path.open('x') as f:
                    json.dump(dict(identity_sha256=identity_pin, sha256=sha256(candidate), proof=proof), f, indent=2)
            if output.exists():
                require(job.get('resume', False) and sha256(output) == sha256(candidate),
                        'source publication collision')
            else:
                os.link(candidate, output)
            receipt = dict(status='SOURCE_COMPLETE' if job['complete_source'] else 'DEVELOPMENT_INTERVAL',
                           identity_sha256=identity_pin, logical_source_id=job['logical_source_id'],
                           output_sha256=sha256(output), input_list_sha256=job['input_list_sha256'],
                           source_manifest_sha256=job['source_manifest_sha256'],
                           pieces=len(chunks), entries=total, proof=proof)
            with receipt_path.open('x') as f:
                json.dump(receipt, f, indent=2)
            progress('complete', total)
            lock.close()
            return 0
        except Exception as error:
            progress('failed', sum(c['event_end']-c['event_begin'] for c in chunks), error=str(error))
            lock.close()
            raise


def publish_groups(native, plan, recipe):
    """Publish each complete run or sample group of a declared population (stage 2).

    A group is complete when every manifest source of its run or sample has its
    final file and a passing consolidation receipt. The receipt the version
    binds is that source's saved consolidation proof. A group that already has
    a published version, possibly advanced to stage 3 by PhotonID, is left as it
    is. Development selections are never published.
    """
    if plan['development']:
        return {}
    members = {}
    for job in plan['jobs']:
        members.setdefault(Path(job['output_file']).parent, []).append(job)
    versions = {}
    for directory, jobs in sorted(members.items()):
        roster = plan['rosters'][directory.name]
        require(sorted(job['logical_source_id'] for job in jobs) == roster,
                'a declared population must select whole groups: ' + str(directory))
        files = []
        for job in jobs:
            output = Path(job['output_file'])
            work = output.parent / '.work' / output.stem
            if not (output.is_file() and (work / 'source_receipt.json').is_file()):
                break
            receipt = read_json(work / 'source_receipt.json')
            require(receipt['status'] == 'SOURCE_COMPLETE' and receipt['output_sha256'] == sha256(output),
                    'source receipt does not bind its published file: ' + str(output))
            validation = work / 'group_receipt.json'
            content = json.dumps(receipt['proof'], indent=2, sort_keys=True) + '\n'
            if validation.exists():
                require(validation.read_text() == content,
                        'group receipt differs from the source proof: ' + str(validation))
            else:
                with validation.open('x') as f:
                    f.write(content)
            files.append(dict(logical_id=job['logical_source_id'], relative_path=output.name, path=str(output),
                              sha256=receipt['output_sha256'], receipt_path=str(validation),
                              receipt_sha256=sha256(validation)))
        if len(files) != len(jobs):
            continue
        if (directory / 'CURRENT').is_symlink():
            published = json.loads(str(native.resolveGroup(str(directory))))['manifest']
            require(sorted(f['logical_id'] for f in published['files']) == roster,
                    'published group roster differs: ' + str(directory))
            versions[directory.name] = published['version']
            continue
        request = dict(group=str(directory), stage=2, files=files, expected_sources=roster, recipe=recipe,
                       group_id=f"{plan['campaign']}/{plan['lane']}/{directory.name}")
        versions[directory.name] = json.loads(str(native.publishGroup(json.dumps(request))))['version']
    return versions


def main():
    """Validate the request, coordinate source jobs, then publish complete groups."""
    global BINDINGS
    parser = argparse.ArgumentParser(prog="produce_trees.sh", add_help=False)
    parser.add_argument("--events", type=int)
    parser.add_argument("--source", type=int)
    parser.add_argument("--run", type=int)
    parser.add_argument("--sample")
    parser.add_argument("--run-list", type=Path)
    parser.add_argument("--sample-list", type=Path)
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--timeout", type=int, default=43200)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("lanes", nargs="+")
    args = parser.parse_args()
    require(sum(value is not None for value in (args.run, args.sample, args.run_list, args.sample_list)) <= 1,
            'choose one of --run, --run-list, --sample or --sample-list')
    args.selected_groups = None
    if args.run_list is not None:
        require(all(lane.endswith('_data') for lane in args.lanes), '--run-list requires a DATA command')
        args.selected_groups = selection_list(args.run_list, numeric=True)
    if args.sample_list is not None:
        require(all(not lane.endswith('_data') for lane in args.lanes), '--sample-list requires a SIM command')
        args.selected_groups = selection_list(args.sample_list)
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
    known_keys(settings, "schema_version campaign output_root production_config producer_library lanes", "binding")
    require(settings.get("schema_version") == 1, "unsupported population bindings schema")
    bindings = settings.get("lanes")
    require(isinstance(bindings, dict), "missing lanes in population bindings")

    missing = []
    for lane in args.lanes:
        require(lane in bindings and isinstance(bindings[lane], dict), f"unknown lane binding: {lane}")
        known_keys(bindings[lane], "profile manifest manifest_sha256 groups index index_sha256", f"{lane} binding")
        require(sum(key in bindings[lane] for key in ('groups', 'index', 'manifest')) <= 1,
                'bind one population index, grouped manifest set or single manifest')
        if not any(bindings[lane].get(key) for key in ('manifest', 'groups', 'index')):
            missing.append(f"{lane}: no population manifest is bound.\n"
                           f"  Expected binding: {BINDINGS} -> lanes.{lane}.manifest")
    require(not missing, "\n".join(missing))

    campaign = text(settings.get("campaign"), f"{BINDINGS} -> campaign")
    require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", campaign), "campaign must be a safe directory name")
    require(MACRO.is_file(), f"production macro is missing: {MACRO}")


    totals = {}
    for plan in population_plans(settings, args):
        row = totals.setdefault(plan['lane'], [0, 0, 0, plan['output']])
        row[0] += 1
        row[1] += len(plan['jobs'])
        row[2] += plan['population_sources']
    for lane, (groups, selected, sources, output) in totals.items():
        print(f"Lane: {lane} | groups: {groups} | sources: {selected}/{sources}")
        print(f"  Output: {output}/<run-or-sample>/part-<ordinal>.root")


    if args.dry_run:
        print("DRY RUN: bindings checked; no outputs created and no ROOT invoked.\n"
              "C++ configuration, conditions, environment and science readiness were not validated.")
        return 0

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

    environment = os.environ.copy()
    include = library.parent.parent / 'include'
    require((include / 'TreeProduction/PhotonJetTree.h').is_file(), "installed TreeProduction headers are missing")
    environment['ROOT_INCLUDE_PATH'] = str(include) + (os.pathsep + environment['ROOT_INCLUDE_PATH']
                                                      if environment.get('ROOT_INCLUDE_PATH') else '')
    for variable in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        inherited = os.pathsep + environment[variable] if environment.get(variable) else ""
        environment[variable] = str(library.parent) + inherited

    library_hash = sha256(library)
    eventtree_library = library.parent / "libEventTree.so"
    require(eventtree_library.is_file(), f"installed eventTree library is missing: {eventtree_library}")
    eventtree_hash = sha256(eventtree_library)
    native = load_eventtree(library.parent)


    # The installed native implementation declares mapping readiness. Check
    # each resolved plan before creating an output directory or processing DSTs.
    for plan in population_plans(settings, args):
        for job in plan['jobs']:
            limit = PIECE_EVENTS.get(job['profile'], PIECE_EVENTS['default'])
            with input_bundle(job) as resolved:
                check_mapping(dict(resolved, events=min(job['events'], limit)), root, environment)

    # Admit the whole aggregate before starting its first event loop.
    args.output_admitted = True
    created_directories = set()
    for plan in population_plans(settings, args):
        directory = Path(plan["output"])
        if directory not in created_directories:
            directory.mkdir(parents=True, exist_ok=args.resume)
            created_directories.add(directory)
        population_dir = directory / 'populations' if 'group' in plan else directory
        population_dir.mkdir(exist_ok=True)
        copy_name = (plan['group'] + '.sources.tsv' if Path(plan['manifest']).suffix == '.tsv'
                     else Path(plan['manifest']).name) if 'group' in plan else 'population.json'
        manifest_copy = (population_dir / copy_name if 'group' in plan
                         else directory / 'population.json')
        if manifest_copy.exists():
            require(sha256(manifest_copy) == plan['manifest_sha256'], 'resume population differs')
        else:
            with manifest_copy.open('xb') as f:
                f.write(Path(plan['manifest']).read_bytes())


    for plan in population_plans(settings, args):
        require(sha256(Path(plan["manifest"])) == plan["manifest_sha256"], "population changed after preflight")
        for job in plan["jobs"]:
            require(sha256(library) == library_hash, "installed producer library changed during this invocation")
            require(sha256(eventtree_library) == eventtree_hash,
                    "installed eventTree library changed during this invocation")
            job["producer_library"] = str(library)
            job["eventtree_library_sha256"] = eventtree_hash
            with input_bundle(job) as resolved:
                status = execute_source(resolved, root, environment, native, library_hash)
            if status:
                return status

        label = "development selection" if plan["development"] else "declared population"
        print(f"FINAL TREE COMPLETE: {plan['lane']} {label}; {len(plan['jobs'])} saved outputs checked.")
        recipe = dict(stage="tree_production", population_manifest_sha256=plan["manifest_sha256"],
                      producer_library_sha256=library_hash, eventtree_library_sha256=eventtree_hash,
                      schema_sha256=sha256(SCHEMA))
        for group, version in publish_groups(native, plan, recipe).items():
            print(f"  Published group {group}: version {version}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FrontendError, OSError, ValueError, KeyError, TypeError, csv.Error, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("ERROR: interrupted; partial output remains unaccepted.", file=sys.stderr)
        sys.exit(130)
