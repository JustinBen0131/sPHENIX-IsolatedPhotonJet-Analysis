"""Private population/execution helper for produce_trees.sh; no event physics.

JSON keeps population planning dependency-free. YAML and physical input-list
interpretation belong exclusively to production::LoadPlan. ROOT and uproot are
needed only for execution and the terminal completion check, never for help or
dry-run. An exact, reviewed population must be supplied before anything runs.
"""

# /*
#  * produce.py
#  *
#  * Private execution layer behind produce_trees.sh.
#  *
#  * This module translates a user-level production request into a set of
#  * completely pinned per-source Fun4All jobs. It owns:
#  *
#  *   - canonical population manifest binding;
#  *   - source/job metadata validation;
#  *   - SHA256 pinning of mutable production inputs;
#  *   - deterministic output placement;
#  *   - invocation of the existing Fun4All steering macro;
#  *   - terminal provenance and accounting validation.
#  *
#  * It deliberately does NOT own:
#  *
#  *   - YAML reconstruction-policy interpretation;
#  *   - DST node semantics;
#  *   - detector calibration or reconstruction;
#  *   - photon, jet, truth or isolation algorithms;
#  *   - scientific event/object selection.
#  *
#  * Those remain downstream of this helper:
#  *
#  *   produce_trees.sh
#  *        |
#  *        v
#  *   produce.py                  [population + execution]
#  *        |
#  *        v
#  *   Fun4All_PhotonJetTree.C     [thin steering]
#  *        |
#  *        v
#  *   production::LoadPlan        [configuration resolution]
#  *   Production*.cc              [reconstruction registration]
#  *        |
#  *        v
#  *   PhotonJetTree               [scientific capture]
#  *
#  * All normal, development-restricted and aggregate commands eventually pass
#  * through execute_source(). There is one per-source production path.
#  */

import argparse
from collections import deque
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


# /*
#  * Canonical TreeProduction locations
#  *
#  * TREE_DIR is resolved from this file, never from the caller's cwd.
#  *
#  * BINDINGS owns the mapping from canonical lane names to reviewed population
#  * manifests.
#  *
#  * MACRO is the single Fun4All steering entrypoint used for every source.
#  */

TREE_DIR = Path(__file__).resolve().parents[1]
BINDINGS = TREE_DIR / "config" / "production_populations.json"
MACRO = TREE_DIR / "macros" / "Fun4All_PhotonJetTree.C"


# /*
#  * Frontend validation
#  *
#  * Expected configuration/input problems are represented as FrontendError so
#  * the command-line boundary can report them cleanly without exposing Python
#  * tracebacks to a production user.
#  */

class FrontendError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise FrontendError(message)


# /*
#  * Primitive value validation
#  *
#  * These helpers provide one rule for externally supplied strings and integer
#  * coordinates. Population manifests are treated as explicit contracts rather
#  * than permissive configuration files.
#  */

def text(value, name, allow_empty=False):
    require(isinstance(value, str) and (allow_empty or value.strip()), f"missing {name}")
    require(not any(ord(c) < 32 for c in value), f"control character in {name}")
    return value


def integer(value, name, minimum=0, maximum=2147483647):
    require(type(value) is int and minimum <= value <= maximum,
            f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


# /*
#  * Strict JSON loading
#  *
#  * Duplicate keys and unknown fields are rejected. A production manifest must
#  * therefore have one unambiguous interpretation; misspelled or obsolete
#  * fields cannot be silently ignored.
#  */

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


# /*
#  * Content-addressed file binding
#  *
#  * Every externally referenced production artifact is checked against its
#  * declared SHA256 before it is admitted into a job plan. Paths identify where
#  * a file can be read; hashes identify the bytes that were approved.
#  */

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


# /*
#  * Population planning
#  *
#  * Resolve one canonical production lane into immutable per-source Job inputs.
#  *
#  * This function interprets only the population contract. It does not inspect
#  * the physical DST input-list format and does not interpret tree_production
#  * YAML. Both remain the responsibility of production::LoadPlan.
#  *
#  * The result contains enough information to:
#  *
#  *   1. identify every selected physical source bundle;
#  *   2. verify every mutable control artifact before execution;
#  *   3. construct one canonical output location per source;
#  *   4. invoke the same steering macro for every production mode.
#  */

def plan_population(lane, binding, settings, args):
    """Resolve job metadata only; never interpret DST columns or YAML policy."""

    # /*
    #  * Bind and validate the lane's reviewed population manifest.
    #  */

    manifest = pinned_file(BINDINGS.parent, binding.get("manifest"),
                           binding.get("manifest_sha256"), f"{lane} manifest")
    population = read_json(manifest)
    known_keys(population, "schema_version lane campaign sources", "manifest")
    require(population.get("schema_version") == 1, f"unsupported manifest schema: {manifest}")
    require(population.get("lane") == lane, f"manifest lane mismatch: {manifest}")
    require(population.get("campaign") == settings["campaign"], f"manifest campaign mismatch: {manifest}")
    sources = population.get("sources")
    require(isinstance(sources, list) and sources, f"empty/missing source population: {manifest}")

    # /*
    #  * Resolve the shared production profile and default YAML configuration.
    #  */

    profile = text(binding.get("profile"), f"lanes.{lane}.profile")
    config = (BINDINGS.parent / text(settings.get("production_config"), "production_config")).resolve()
    require(config.is_file(), f"production configuration is missing: {config}")

    # /*
    #  * Convert each manifest source into one fully pinned macro invocation.
    #  *
    #  * Source ordinals and physical input-list bundles are unique within the
    #  * lane. Per-source production configuration overrides are allowed only
    #  * when both path and SHA256 are explicitly supplied.
    #  */

    jobs, ordinals, bundles = [], set(), set()
    for source in sources:
        require(isinstance(source, dict), f"{lane}: source must be an object")
        known_keys(source, "source_file_ordinal input_list input_list_sha256 run segment "
                          "first_entry event_count sample period si_di_role production_config "
                          "production_config_sha256", "source")
        require(("production_config" in source) == ("production_config_sha256" in source),
                "source production_config and production_config_sha256 must be supplied together")
        ordinal = integer(source.get("source_file_ordinal"), "source_file_ordinal", maximum=2**63 - 1)
        require(ordinal not in ordinals, f"{lane}: duplicate source ordinal {ordinal}")
        ordinals.add(ordinal)
        first = integer(source.get("first_entry"), "first_entry")
        events = integer(source.get("event_count"), "event_count", minimum=1)
        inputs = pinned_file(manifest.parent, source.get("input_list"),
                             source.get("input_list_sha256"), "input_list")

        # /*
        #  * One manifest row represents one physical input bundle.
        #  *
        #  * Sharding, retry and resume semantics are intentionally not inferred
        #  * from filenames. They require a future explicit operational contract.
        #  */

        require(source["input_list_sha256"] not in bundles, f"{lane}: duplicate input-list bundle")
        bundles.add(source["input_list_sha256"])
        job_config = config
        if "production_config" in source:
            job_config = pinned_file(manifest.parent, source["production_config"],
                                     source.get("production_config_sha256"), "source production_config")

        # /*
        #  * This dictionary mirrors the arguments of Fun4All_PhotonJetTree.C.
        #  * Scientific reconstruction choices are still resolved later by
        #  * production::LoadPlan from profile + YAML.
        #  */

        job = dict(profile=profile, input_list=str(inputs), output_file="",
                   run=integer(source.get("run"), "run"),
                   segment=integer(source.get("segment"), "segment"),
                   events=min(events, args.events) if args.events is not None else events,
                   first_entry=first, config=str(job_config), source_file_ordinal=ordinal,
                   input_file_sha256="", source_manifest_sha256=binding["manifest_sha256"],
                   sample=text(source.get("sample"), "sample"),
                   period=text(source.get("period"), "period", allow_empty=True),
                   si_di_role=text(source.get("si_di_role"), "si_di_role", allow_empty=True),
                   input_list_sha256=source["input_list_sha256"],
                   configuration_sha256=sha256(job_config), macro_sha256=sha256(MACRO), verbosity=0)
        if args.source is None or ordinal == args.source:
            jobs.append(job)

    require(jobs, f"{lane}: source ordinal {args.source} does not exist; nothing was selected")

    # /*
    #  * Output ownership
    #  *
    #  * Full populations live directly under <campaign>/<lane>.
    #  *
    #  * Any --events or --source restriction is segregated into a development
    #  * hierarchy so a canary output cannot be mistaken for a full population.
    #  *
    #  * Existing destinations are never overwritten or implicitly resumed.
    #  */

    output = (BINDINGS.parent / text(settings.get("output_root"), "output_root")).resolve()
    output /= settings["campaign"]
    development = args.events is not None or args.source is not None
    if development:
        output = output / "development" / lane / (
            f"source-{args.source if args.source is not None else 'all'}_"
            f"events-{args.events if args.events is not None else 'manifest'}")
    else:
        output /= lane
    require(not output.exists(), f"output already exists; refusing overwrite or automatic resume: {output}")

    for job in jobs:
        job["output_file"] = str(output / f"source-{job['source_file_ordinal']:08d}" / "trees.root")

    return dict(lane=lane, profile=profile, manifest=str(manifest),
                manifest_sha256=binding["manifest_sha256"], output=str(output),
                development=development, population_sources=len(sources), jobs=jobs)


# /*
#  * Fun4All macro invocation
#  *
#  * Construct exactly one ROOT macro expression from the resolved Job.
#  *
#  * json.dumps provides valid quoting for the C++ string arguments used here.
#  * subprocess later receives a direct argv vector: there is no shell eval,
#  * command interpolation or second parsing layer.
#  */

def macro_argument(job):
    # JSON string quoting is also valid for these C++ UTF-8 string literals.
    # subprocess receives an argument vector: no shell, eval or interpolation.
    values = [job[key] for key in (
        "profile", "input_list", "output_file", "run", "segment", "events",
        "first_entry", "config", "source_file_ordinal", "input_file_sha256",
        "source_manifest_sha256", "sample", "period", "si_di_role", "verbosity")]
    return str(MACRO) + "(" + ",".join(json.dumps(v, ensure_ascii=False) for v in values) + ")"


# /*
#  * Producer metadata reader
#  *
#  * PhotonJetTree writes compact key=value metadata objects. Parse them through
#  * the same duplicate-key guard used for JSON so terminal provenance remains
#  * unambiguous.
#  */

def key_values(obj):
    pairs = [line.split("=", 1) for line in str(obj).splitlines() if "=" in line]
    return unique_object(pairs)


# /*
#  * Terminal completion contract
#  *
#  * A successful ROOT exit and the existence of trees.root are not sufficient
#  * to accept a source.
#  *
#  * The completed product must independently demonstrate:
#  *
#  *   - completion_status == complete;
#  *   - exact source/run/configuration/macro/library provenance;
#  *   - exact declared input-event accounting;
#  *   - retained + upstream-rejected == observed;
#  *   - ROOT table counts consistent with terminal accounting;
#  *   - one completed source record.
#  *
#  * This is an operational completion check. It is not scientific acceptance.
#  */

def check_completion(job, uproot, library_hash):
    """Read the producer's terminal record, not a filename or success banner."""

    with uproot.open(job["output_file"]) as product:
        completion = key_values(product["completion"])
        metadata = key_values(product["metadata"])
        require(completion.get("completion_status") == "complete", "producer output is not complete")

        expected = {key: str(job[key]) for key in (
            "run", "segment", "source_file_ordinal", "first_entry", "source_manifest_sha256",
            "sample", "period", "si_di_role", "configuration_sha256", "macro_sha256")}
        expected.update(max_events=str(job["events"]), producer_library_sha256=library_hash)

        for key, value in expected.items():
            require(metadata.get(key) == value, f"output provenance mismatch: {key}")

        observed = int(completion["upstream_observed_events"])
        encountered = int(completion["encountered_events"])
        retained = int(completion["retained_events"])
        rejected = int(completion["upstream_rejected_events"])

        require(observed == job["events"], "source ended before/after its declared event count")
        require(encountered == retained and observed == retained + rejected,
                "producer completion accounting is inconsistent")
        require(product["Events"].num_entries == retained, "retained event row count mismatch")
        require(product["UpstreamRejectedEvents"].num_entries == rejected, "rejected event row count mismatch")
        require(product["Sources"].num_entries == 1, "expected exactly one source row")
        require(bool(product["Sources"]["completed"].array(library="np")[0]), "source row is incomplete")

        return completion


# /*
#  * Per-source executor
#  *
#  * This is the one execution path used by every primitive lane, aggregate
#  * command and development restriction.
#  *
#  * Each source receives its own exclusively created directory containing:
#  *
#  *   trees.root
#  *   production.log
#  *   invocation.json
#  *   completion.json
#  *
#  * All pinned control artifacts are rehashed immediately before ROOT starts.
#  * A failed process retains its directory and log for diagnosis; it is never
#  * promoted or reused automatically.
#  */

def execute_source(job, root, environment, uproot, library_hash):
    """The only per-source execution path for every command and restriction."""

    directory = Path(job["output_file"]).parent
    directory.mkdir()  # Exclusive reservation; partial outputs are never reused.
    log = directory / "production.log"

    # /*
    #  * Protect the interval between planning and execution. A job is not
    #  * allowed to run if its input list, production YAML or steering macro
    #  * changed after preflight.
    #  */

    require(sha256(Path(job["input_list"])) == job["input_list_sha256"], "input list changed after preflight")
    require(sha256(Path(job["config"])) == job["configuration_sha256"], "configuration changed after preflight")
    require(sha256(MACRO) == job["macro_sha256"], "macro changed after preflight")

    argv = [root, "-l", "-b", "-q", macro_argument(job)]

    # /*
    #  * Persist the exact resolved Job and argv before execution so the source
    #  * can be audited even if ROOT terminates before producing an output file.
    #  */

    (directory / "invocation.json").write_text(json.dumps(dict(job=job, argv=argv), indent=2) + "\n")
    print(f"  source {job['source_file_ordinal']}: {job['events']} events; log {log}", flush=True)

    with log.open("x") as stream:
        result = subprocess.run(argv, cwd=TREE_DIR, env=environment, stdout=stream, stderr=subprocess.STDOUT)

    # /*
    #  * Preserve normal process failure semantics and print only a small tail
    #  * of the log at the frontend. The complete log remains beside the source.
    #  */

    if result.returncode:
        print(f"ERROR: producer failed ({result.returncode}); see {log}", file=sys.stderr)
        with log.open(errors="replace") as stream:
            for line in deque(stream, maxlen=8):
                print("  " + line.rstrip(), file=sys.stderr)
        return result.returncode if result.returncode > 0 else 128 - result.returncode

    # /*
    #  * ROOT returned successfully; now independently verify the producer's own
    #  * terminal provenance and accounting before recording source completion.
    #  */

    try:
        completion = check_completion(job, uproot, library_hash)
    except Exception as error:
        raise FrontendError(f"completion check failed: {error}; output retained for inspection: {directory}") from error

    (directory / "completion.json").write_text(json.dumps(
        dict(status="producer_complete_checked", science_accepted=False, completion=completion), indent=2) + "\n")

    return 0


# /*
#  * Production orchestration
#  *
#  * main() performs the complete invocation in four phases:
#  *
#  *   1. validate CLI restrictions;
#  *   2. resolve and preflight every requested population;
#  *   3. establish the ROOT/runtime environment and reserve every lane output;
#  *   4. execute sources through execute_source() and verify completion.
#  *
#  * No event loop begins until every requested lane passes population preflight.
#  */

def main():
    parser = argparse.ArgumentParser(prog="produce_trees.sh", add_help=False)
    parser.add_argument("--events", type=int)
    parser.add_argument("--source", type=int)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("lanes", nargs="+")
    args = parser.parse_args()

    # /*
    #  * Development restriction validation.
    #  */

    if args.events is not None:
        integer(args.events, "--events", minimum=1)
    if args.source is not None:
        integer(args.source, "--source", maximum=2**63 - 1)
        require(len(args.lanes) == 1, "--source requires a primitive command")

    # /*
    #  * Read the repository-level population-binding document.
    #  *
    #  * An unbound lane is a hard preflight failure. There is no nearby-file
    #  * search, directory glob, fallback population or inferred calibration.
    #  */

    settings = read_json(BINDINGS)
    known_keys(settings, "schema_version campaign output_root production_config producer_library lanes", "binding")
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

    # /*
    #  * Resolve the complete request before starting any source.
    #  */

    plans = [plan_population(lane, bindings[lane], settings, args) for lane in args.lanes]

    # /*
    #  * Human-readable resolved plan.
    #  *
    #  * This is deliberately population-level information. Detailed detector
    #  * configuration remains visible in the subsequent C++ startup report
    #  * emitted by production::Print().
    #  */

    for plan in plans:
        print(f"Lane: {plan['lane']} | profile: {plan['profile']}")
        print(f"  Population: {plan['manifest']}")
        for config in sorted({job['config'] for job in plan['jobs']}):
            print(f"  Configuration: {config}")
        print(f"  Sources: {len(plan['jobs'])}/{plan['population_sources']} | "
              f"{'DEVELOPMENT RESTRICTION' if plan['development'] else 'declared full population'}")
        print(f"  Output/logs: {plan['output']}/source-<ordinal>/")

    # /*
    #  * Dry-run boundary
    #  *
    #  * A dry run proves only that the frontend bindings and population plan are
    #  * internally resolvable. ROOT, C++ configuration, conditions and scientific
    #  * readiness are intentionally not claimed.
    #  */

    if args.dry_run:
        print("DRY RUN: bindings checked; no outputs created and no ROOT invoked.\n"
              "C++ configuration, conditions, environment and science readiness were not validated.")
        return 0

    # /*
    #  * Runtime environment
    #  *
    #  * Full execution requires:
    #  *
    #  *   - ROOT on PATH;
    #  *   - the installed libPhotonJetTree.so named exactly as expected by
    #  *     R__LOAD_LIBRARY;
    #  *   - uproot for independent terminal-file validation.
    #  */

    root = shutil.which("root")
    require(root, "ROOT is not on PATH; enter the configured sPHENIX environment")

    library = (BINDINGS.parent / text(settings.get("producer_library"),
               f"{BINDINGS} -> producer_library (installed libPhotonJetTree.so)")).absolute()

    # Keep the installed .so symlink name: versioned targets are normal, while
    # R__LOAD_LIBRARY still resolves libPhotonJetTree.so in its parent directory.
    require(library.is_file() and library.name == "libPhotonJetTree.so",
            f"installed producer library is missing or misnamed: {library}")

    try:
        import uproot
    except ImportError as error:
        raise FrontendError("uproot is required to verify the producer's completion record") from error

    # /*
    #  * Make the explicitly bound producer library discoverable to the dynamic
    #  * loader without changing the caller's existing library path contents.
    #  */

    environment = os.environ.copy()
    for variable in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        environment[variable] = str(library.parent) + (os.pathsep + environment[variable] if environment.get(variable) else "")

    library_hash = sha256(library)

    # /*
    #  * Aggregate admission
    #  *
    #  * Reserve every requested lane directory and snapshot its population
    #  * manifest before the first event loop begins. This makes partial aggregate
    #  * admission impossible: every lane has passed frontend preflight first.
    #  */

    # Admit the whole aggregate before starting its first event loop.
    for plan in plans:
        directory = Path(plan["output"])
        directory.mkdir(parents=True, exist_ok=False)
        (directory / "population.json").write_bytes(Path(plan["manifest"]).read_bytes())

    # /*
    #  * Source execution
    #  *
    #  * Recheck manifest and installed-library identity while the invocation is
    #  * in progress, then execute every source through the single executor.
    #  */

    for plan in plans:
        require(sha256(Path(plan["manifest"])) == plan["manifest_sha256"], "population changed after preflight")
        for job in plan["jobs"]:
            require(sha256(library) == library_hash, "installed producer library changed during this invocation")
            status = execute_source(job, root, environment, uproot, library_hash)
            if status:
                return status

        label = "development selection" if plan["development"] else "declared population"
        print(f"SUCCESS: {plan['lane']} {label}; {len(plan['jobs'])} source completion records checked.")

    return 0


# /*
#  * Process boundary
#  *
#  * Expected operational/configuration failures produce one concise ERROR line.
#  * Ctrl-C receives the conventional exit status 130 and explicitly states that
#  * any partial output is not accepted.
#  */

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FrontendError, OSError, ValueError, KeyError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("ERROR: interrupted; partial output remains unaccepted.", file=sys.stderr)
        sys.exit(130)
