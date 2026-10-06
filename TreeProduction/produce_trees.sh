#!/usr/bin/env bash
# Public interface: select a maintained population and reconstruction preset.
# All commands use one production runner; detector algorithms remain in C++.


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
  runPhotonJetDISim
  runInclusiveJetDISim
  runEmbeddedPhotonJetSim
  runEmbeddedInclusiveJetSim
  runAllSim

Everything:
  runAll

Options:
  --run N      Select one DATA run from the maintained catalogue.
  --run-list FILE  Select DATA runs from a plain list (one run number per line).
  --sample ID  Select one exact SIM sample (for example run28_photonjet20).
  --sample-list FILE  Select SIM samples from a plain list (one name per line).
  --events N   Process at most N events per selected source (development).
  --source N   Select manifest source ordinal N (primitive commands only).
  --catalogue FILE  Use an explicitly prepared alternative catalogue.
  --output DIR Choose the output root; existing products are refused.
  --resume     Reuse only identity-matched, validated source/piece outputs.
  --timeout N  Maximum seconds per reconstruction piece (default 43200).
  --dry-run    Check bindings and show the plan; create nothing, run no ROOT.
  --help       Show this help.

Populations: config/production_populations.yaml (or --catalogue FILE).
Unbound populations fail; no input population or calibration is guessed.
EOF
}

expand_command() {
    case "$1" in
        runPP)                       printf '%s\n' pp_data ;;
        runAuAu)                     printf '%s\n' auau_data ;;
        runPhotonJetSim)             printf '%s\n' pp_photon_sim ;;
        runInclusiveJetSim)          printf '%s\n' pp_inclusive_sim ;;
        runPhotonJetDISim)           printf '%s\n' pp_photon_di ;;
        runInclusiveJetDISim)        printf '%s\n' pp_inclusive_di ;;
        runEmbeddedPhotonJetSim)     printf '%s\n' auau_photon_embedded ;;
        runEmbeddedInclusiveJetSim)  printf '%s\n' auau_inclusive_embedded ;;
        runAllData) expand_command runPP; expand_command runAuAu ;;
        runAllSim)
            expand_command runPhotonJetSim
            expand_command runInclusiveJetSim
            expand_command runPhotonJetDISim
            expand_command runInclusiveJetDISim
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

helper_args=("$tree_dir/scripts/produce.py")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help) usage; exit 0 ;;
        --dry-run|--resume) helper_args+=("$1"); shift ;;
        --events|--source|--run|--timeout)
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
        --sample|--run-list|--sample-list|--catalogue|--output)
            if [[ $# -lt 2 || -z $2 || $2 == --* ]]; then
                printf 'ERROR: %s requires a value. See --help.\n' "$1" >&2
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

exec python3 "${helper_args[@]}" "${lanes[@]}"
