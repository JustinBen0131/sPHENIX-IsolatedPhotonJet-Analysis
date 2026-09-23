#!/usr/bin/env bash

# /*
#  * produce_trees.sh
#  *
#  * Public command-line interface for canonical PhotonJetTree production.
#  *
#  * This file owns only the vocabulary exposed to a user:
#  *
#  *   ./produce_trees.sh runPP
#  *   ./produce_trees.sh runAuAu
#  *   ./produce_trees.sh runPhotonJetSim
#  *   ...
#  *   ./produce_trees.sh runAll
#  *
#  * It does not know how detector reconstruction is performed, which nodes are
#  * consumed, which calibration payloads are used, or how PhotonJetTree fills
#  * its output. Those responsibilities remain below this interface:
#  *
#  *   produce_trees.sh
#  *        |
#  *        v
#  *   scripts/produce.py
#  *        |
#  *        v
#  *   Fun4All_PhotonJetTree.C
#  *        |
#  *        v
#  *   Production.cc / Production*.cc
#  *        |
#  *        v
#  *   PhotonJetTree
#  *
#  * Aggregate commands are composed entirely from the same primitive lanes.
#  * There is therefore one production path for a single-source development
#  * test, a complete lane, and a full runAll production.
#  */

set -euo pipefail


# /*
#  * Repository location
#  *
#  * Resolve TreeProduction from the location of this script rather than from
#  * the caller's current working directory. All subsequent paths are therefore
#  * stable whether the command is launched from the repository root, this
#  * directory, or another working directory.
#  */

tree_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)


# /*
#  * Public interface
#  *
#  * These are the only production concepts a normal user needs to know.
#  * Population manifests, source ordinals, output paths and the Fun4All macro
#  * argument list are intentionally hidden behind this layer.
#  */

usage() {
    cat <<'EOF'
sPHENIX Isolated-Photon + Jet Tree Production

Usage:
  ./produce_trees.sh COMMAND [options]

Data:
  runPP
  runAuAu
  runAllData

Simulation:
  runPhotonJetSim
  runInclusiveJetSim
  runEmbeddedPhotonJetSim
  runEmbeddedInclusiveJetSim
  runAllSim

Everything:
  runAll

Options:
  --events N   Process at most N events per selected source (development).
  --source N   Select manifest source ordinal N (primitive commands only).
  --dry-run    Check bindings and show the plan; create nothing, run no ROOT.
  --help       Show this help.

Populations: config/production_populations.json
Unbound populations fail; no input population or calibration is guessed.
EOF
}


# /*
#  * Command -> production-lane mapping
#  *
#  * Primitive commands map directly onto the canonical lane names understood
#  * by production_populations.json.
#  *
#  * Aggregate commands recurse through those same primitive mappings. They do
#  * not have an independent implementation:
#  *
#  *   runAllData
#  *      -> runPP
#  *      -> runAuAu
#  *
#  *   runAllSim
#  *      -> runPhotonJetSim
#  *      -> runInclusiveJetSim
#  *      -> runEmbeddedPhotonJetSim
#  *      -> runEmbeddedInclusiveJetSim
#  *
#  *   runAll
#  *      -> runAllData
#  *      -> runAllSim
#  *
#  * This is what guarantees that development and full-production workflows do
#  * not drift into separate scientific code paths.
#  */

expand_command() {
    case "$1" in
        runPP)                       printf '%s\n' pp_data ;;
        runAuAu)                     printf '%s\n' auau_data ;;
        runPhotonJetSim)             printf '%s\n' pp_photon_sim ;;
        runInclusiveJetSim)          printf '%s\n' pp_inclusive_sim ;;
        runEmbeddedPhotonJetSim)     printf '%s\n' auau_photon_embedded ;;
        runEmbeddedInclusiveJetSim)  printf '%s\n' auau_inclusive_embedded ;;
        runAllData) expand_command runPP; expand_command runAuAu ;;
        runAllSim)
            expand_command runPhotonJetSim
            expand_command runInclusiveJetSim
            expand_command runEmbeddedPhotonJetSim
            expand_command runEmbeddedInclusiveJetSim
            ;;
        runAll) expand_command runAllData; expand_command runAllSim ;;
        *) return 2 ;;
    esac
}


# /*
#  * Command selection
#  *
#  * Resolve the public command into its complete primitive lane set before
#  * forwarding anything to the private helper.
#  */

if [[ $# -eq 0 || $1 == --help ]]; then
    usage
    exit 0
fi

command_name=$1
shift

if ! expanded=$(expand_command "$command_name"); then
    printf 'ERROR: unknown command: %s\n\n' "$command_name" >&2
    usage >&2
    exit 2
fi

lanes=()
while IFS= read -r lane; do lanes+=("$lane"); done <<< "$expanded"


# /*
#  * Development restrictions
#  *
#  * --events and --source narrow the canonical production path; they never
#  * select a different producer or macro.
#  *
#  * --source is meaningful only for one primitive population. Applying a source
#  * ordinal to an aggregate command would be ambiguous and is rejected here
#  * before population resolution begins.
#  *
#  * Keep helper_args nonempty: Bash 3.2 treats an empty array as unset under -u.
#  */

helper_args=("$tree_dir/scripts/produce.py")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help) usage; exit 0 ;;
        --dry-run) helper_args+=("$1"); shift ;;
        --events|--source)
            if [[ $# -lt 2 || ! $2 =~ ^[0-9]+$ ]]; then
                printf 'ERROR: %s requires an integer. See --help.\n' "$1" >&2
                exit 2
            fi
            if [[ $1 == --source && ${#lanes[@]} -ne 1 ]]; then
                printf 'ERROR: --source requires a primitive command. See --help.\n' >&2
                exit 2
            fi
            helper_args+=("$1" "$2")
            shift 2
            ;;
        *) printf 'ERROR: unknown option: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done


# /*
#  * Private orchestration dependency
#  *
#  * Python is used only for population binding, provenance validation and
#  * execution orchestration. Event reconstruction itself remains entirely in
#  * the existing C++/Fun4All TreeProduction path.
#  */

if ! command -v python3 >/dev/null 2>&1; then
    printf 'ERROR: python3 is required for population binding.\n' >&2
    exit 1
fi


# /*
#  * Handoff
#  *
#  * The complete expanded lane set is passed in one invocation so produce.py
#  * can preflight the entire request before the first source begins running.
#  * This prevents an aggregate production from partially starting when a later
#  * lane has an unresolved population or other binding problem.
#  */

exec python3 "${helper_args[@]}" "${lanes[@]}"
