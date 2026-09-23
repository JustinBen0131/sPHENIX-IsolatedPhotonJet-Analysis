"""Operational dispatch only: configuration, completed inputs, output locations.

The shell entrypoints deliberately share this parser and execution route. This
file never defines features, labels, cuts, weights, or reconstruction constants.
"""

from pathlib import Path
import argparse
import sys

HERE = Path(__file__).resolve().parent
TRAIN_HELP = """sPHENIX Photon ID — candidate model training

Commands:
  trainAuAu              canonical_auau_v1
  trainPP                canonical_pp_v1
  trainPPG12Equivalent   ppg12_equivalent_v1 (historical mechanics)
  trainAll               the same three primitive paths, sequentially

Missing population/physics bindings fail closed. No automatic acceptance.
"""
AUGMENT_HELP = """sPHENIX Photon ID — immutable-base ROOT sidecars

Data:
  runPP  runAuAu  runAllData
Simulation:
  runPhotonJetSim  runInclusiveJetSim
  runEmbeddedPhotonJetSim  runEmbeddedInclusiveJetSim  runAllSim
Everything:
  runAll

Default: nominal_pp or nominal_auau selection for the lane's system.
--models requests raw scores without requiring a selection recipe.
"""


def parser(mode):
    p = argparse.ArgumentParser(
        prog="./"
        + ("train_photon_id.sh" if mode == "train" else "augment_photon_id.sh"),
        description=TRAIN_HELP if mode == "train" else AUGMENT_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("command", nargs="?", metavar="COMMAND")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="validate bindings and show plan; do not train/write",
    )
    p.add_argument(
        "--registry",
        type=Path,
        default=HERE / "model_registry.yaml",
        help="registry YAML",
    )
    p.add_argument("--output-root", type=Path, help="override output/photon_id root")
    p.add_argument(
        "--campaign",
        type=Path,
        help="development override: completed TreeProduction manifest",
    )
    p.add_argument("--campaign-sha256", help="required hash with --campaign override")
    if mode == "train":
        p.add_argument(
            "--profiles",
            type=Path,
            default=HERE / "training_profiles.yaml",
            help="training profile YAML",
        )
    else:
        p.add_argument(
            "--source",
            type=int,
            help="one-based source within a primitive lane; same execution path",
        )
        p.add_argument(
            "--models",
            help="comma-separated model names, or all-compatible; raw scoring",
        )
        p.add_argument(
            "--selections",
            "--selection",
            dest="selections",
            help="comma-separated recipe names, or all-compatible",
        )
    return p


def manifest_for(args, raw, base):
    from registry import campaign

    if bool(args.campaign) != bool(args.campaign_sha256):
        raise ValueError("--campaign and --campaign-sha256 must be supplied together")
    path = args.campaign or (
        base / raw["input_campaign"] if raw.get("input_campaign") else None
    )
    sha = args.campaign_sha256 or raw.get("input_campaign_sha256")
    return campaign(path, sha)


def model_choices(system, args, raw, specs):
    from registry import resolve_selection

    selection_names = []
    if args.selections == "all-compatible":
        selection_names = [
            n
            for n, r in raw["selection_recipes"].items()
            if specs[r["model"]].system == system
        ]
    elif args.selections:
        selection_names = args.selections.split(",")
    elif not args.models:
        selection_names = ["nominal_" + system]
    recipes = [
        resolve_selection(raw, n, specs, args.registry.parent) for n in selection_names
    ]
    if args.models == "all-compatible":
        names = [n for n, s in specs.items() if s.system == system and s.bound]
        skipped = [n for n, s in specs.items() if s.system == system and not s.bound]
        if skipped:
            print("Unbound models excluded from all-compatible: " + ", ".join(skipped))
    else:
        names = args.models.split(",") if args.models else []
    names = list(dict.fromkeys(names + [r["model"] for r in recipes]))
    if not names:
        raise ValueError(f"no bound compatible models for {system}")
    for name in names:
        if name not in specs:
            raise ValueError(f"unknown model {name}")
        if specs[name].system != system:
            raise ValueError(f"{name} is incompatible with {system}")
        if not specs[name].bound:
            raise ValueError(f"{name}: model artifact/receipt binding is unresolved")
    return [specs[n] for n in names], recipes


def main(mode, argv=None):
    p = parser(mode)
    args = p.parse_args(argv)
    if not args.command:
        p.print_help()
        return 0
    try:
        from registry import (
            read,
            TRAIN,
            LANES,
            expand,
            write_json,
            file_hash,
            digest,
            bound_file,
        )
        from features import load_registry

        raw = read(args.registry)
        specs = load_registry(args.registry)
        root = (args.output_root or args.registry.parent / raw["output_root"]).resolve()
        if mode == "train":
            from train import check_profile, train_package
            from inputs import load_training
            from augment import make_evaluator

            config = read(args.profiles)["profiles"]
            commands = list(TRAIN) if args.command == "trainAll" else [args.command]
            if any(c not in TRAIN for c in commands):
                p.print_help()
                raise ValueError("unknown training command")
            plans = []
            for c in commands:
                name = TRAIN[c]
                profile = config[name]
                spec = specs[name]
                print(
                    f'{c} → {name} | profile: {args.profiles} | output: {root/"training"/name}'
                )
                check_profile(spec, profile)
                camp = manifest_for(args, profile, args.profiles.parent)
                out = root / "training" / name
                if out.exists():
                    raise FileExistsError(
                        f"package already exists: {out}; reuse/review before a deliberate new version"
                    )
                # Input/identity validation belongs to preflight, including dry run.
                if args.dry_run:
                    load_training(camp, spec, profile)
                data = None  # trainAll never retains several full training populations
                import subprocess

                result = subprocess.run(
                    ["git", "rev-parse", "HEAD"],
                    cwd=HERE,
                    capture_output=True,
                    text=True,
                )
                provenance = dict(
                    campaign=camp["campaign"],
                    campaign_sha256=camp["_sha256"],
                    campaign_manifest=str(camp["_manifest"]),
                    source_manifests=[
                        {k: v for k, v in s.items() if not k.startswith("_")}
                        for s in camp["sources"]
                    ],
                    registry_sha256=file_hash(args.registry),
                    profile_file_sha256=file_hash(args.profiles),
                    model_version=raw["models"][name]["version"],
                    checked_inputs=[
                        dict(path=str(s["_path"]), sha256=s["sha256"])
                        for s in camp["sources"]
                    ]
                    + [
                        dict(
                            path=str(
                                bound_file(
                                    s["training_witness"], camp["_manifest"].parent
                                )
                            ),
                            sha256=s["training_witness"]["sha256"],
                        )
                        for s in camp["sources"]
                        if s.get("training_witness")
                    ],
                    software_git_commit=(
                        result.stdout.strip() if result.returncode == 0 else None
                    ),
                    software_sha256={q.name: file_hash(q) for q in HERE.glob("*.py")},
                    development_input_override=bool(args.campaign),
                )
                reference = None
                if name == "ppg12_equivalent_v1" and specs["ppg12_original"].bound:
                    ref = specs["ppg12_original"]
                    if (
                        ref.features != spec.features
                        or ref.shower_definition != spec.shower_definition
                        or ref.ratio_policy != spec.ratio_policy
                        or ref.domain != spec.domain
                    ):
                        raise ValueError(
                            "reference comparison requires the same bound feature adapter and evaluation domain"
                        )
                    reference = make_evaluator(ref)
                    provenance["reference_model_sha256"] = ref.model_sha256
                plans.append((spec, profile, data, out, provenance, reference, camp))
            if args.dry_run:
                print(
                    "DRY RUN: all requested profiles passed input/binding preflight; no model trained."
                )
                return 0
            for spec, profile, data, out, provenance, reference, camp in plans:
                data = load_training(camp, spec, profile)
                receipt = train_package(spec, profile, data, out, provenance, reference)
                del data
                print(
                    f"{spec.name}: BUILD/MECHANICAL PASS — candidate, SCIENTIFIC ACCEPTANCE NOT REVIEWED. {out}"
                )
        else:
            from augment import make_evaluator, augment_file

            try:
                commands = expand(args.command)
            except ValueError:
                p.print_help()
                raise
            if args.source is not None and (len(commands) != 1 or args.source < 1):
                raise ValueError(
                    "--source requires a primitive lane and a positive one-based index"
                )
            camp = manifest_for(args, raw, args.registry.parent)
            plans = []
            for command in commands:
                lane = LANES[command]
                system = "auau" if lane.startswith("auau_") else "pp"
                sources = [s for s in camp["sources"] if s["lane"] == lane]
                if not sources:
                    raise ValueError(
                        f"{lane}: completed campaign has no source binding"
                    )
                if args.source is not None:
                    if args.source > len(sources):
                        raise ValueError("--source exceeds bound lane population")
                    sources = [sources[args.source - 1]]
                models, recipes = model_choices(system, args, raw, specs)
                # Content/version-specific directory prevents ambiguous overwrites.
                configuration = digest(
                    dict(
                        models=[(s.name, s.model_sha256) for s in models],
                        recipes=recipes,
                    )
                )[:16]
                directory = (
                    root / "augmentation" / camp["campaign"] / lane / configuration
                )
                evaluators = {s.name: make_evaluator(s) for s in models}
                print(
                    f'{command} → {lane} | models: {", ".join(s.name for s in models)} | sources: {len(sources)} | output: {directory}'
                )
                for source in sources:
                    a, b = source["source_id"]
                    output = directory / f"{a:016x}{b:016x}.photon_id.root"
                    if output.exists():
                        raise FileExistsError(f"refusing existing sidecar: {output}")
                    plans.append((source, output, models, evaluators, recipes))
            if args.dry_run:
                print(
                    "DRY RUN: manifest/model/recipe preflight passed; no sidecar written."
                )
                return 0
            for source, out, models, evaluators, recipes in plans:
                result = augment_file(
                    source,
                    out,
                    models,
                    evaluators,
                    recipes,
                    raw["models"],
                    dict(name=camp["campaign"], manifest_sha256=camp["_sha256"]),
                )
                print(
                    f'Sidecar complete: {result["output"]} | SHA256 {result["sha256"]}'
                )
        return 0
    except (ValueError, KeyError, OSError, ImportError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1], sys.argv[2:]))
