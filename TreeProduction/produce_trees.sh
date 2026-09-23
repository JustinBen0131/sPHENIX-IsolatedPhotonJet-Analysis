#!/usr/bin/env bash
# Public command vocabulary only. The helper binds sources and invokes the
# existing macro; Production.cc continues to own all reconstruction policy.
set -euo pipefail

tree_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

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

# Aggregates expand recursively into the same primitive lanes, in this order.
expand_command() {
    case "$1" in
        runPP)                       printf '%s\n' pp_data ;;
        runAuAu)                     printf '%s\n' auau_data ;;
        runPhotonJetSim)              printf '%s\n' pp_photon_sim ;;
        runInclusiveJetSim)           printf '%s\n' pp_inclusive_sim ;;
        runEmbeddedPhotonJetSim)      printf '%s\n' auau_photon_embedded ;;
        runEmbeddedInclusiveJetSim)   printf '%s\n' auau_inclusive_embedded ;;
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
# Keep this array nonempty: Bash 3.2 treats an empty array as unset under -u.
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

if ! command -v python3 >/dev/null 2>&1; then
    printf 'ERROR: python3 is required for population binding.\n' >&2
    exit 1
fi
# One helper call preflights the entire aggregate before any source is run.
exec python3 "${helper_args[@]}" "${lanes[@]}"
