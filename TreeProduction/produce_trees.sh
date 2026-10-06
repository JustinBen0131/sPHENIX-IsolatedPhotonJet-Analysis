#!/usr/bin/env bash
# Public interface: select a maintained population and reconstruction preset.
# All commands use one production runner; detector algorithms remain in C++.
#
# =============================================================================
# TreeProduction public launcher
# =============================================================================
#
# Purpose
# -------
# This script is the thin, human-facing entry point for TreeProduction.
#
# It deliberately does NOT:
#   - inspect DST contents,
#   - resolve physics configuration,
#   - choose calibration constants,
#   - run detector reconstruction itself,
#   - validate ROOT output,
#   - assemble production pieces, or
#   - publish completed tree groups.
#
# Instead, it has three narrowly defined responsibilities:
#
#   1. Translate readable public commands such as
#
#          runPP
#          runPhotonJetSim
#          runAll
#
#      into the internal production-lane names understood by produce.py.
#
#   2. Perform lightweight command-line validation so obvious invocation
#      errors fail immediately with a useful message.
#
#   3. Transfer control to:
#
#          scripts/produce.py
#
#      which owns population binding, source scheduling, reconstruction
#      orchestration, output validation, source assembly, and publication.
#
# The physics/reconstruction implementation itself remains in the C++ producer.
#
# Overall control flow
# --------------------
#
#   user
#     |
#     v
#   produce_trees.sh
#     |
#     |  public command -> internal lane(s)
#     |  validate/forward CLI options
#     v
#   scripts/produce.py
#     |
#     |  resolve maintained population
#     |  construct source jobs
#     |  invoke ROOT/Fun4All
#     |  validate and assemble outputs
#     v
#   C++ TreeProduction implementation
#
# Keeping this wrapper intentionally small prevents shell-level logic from
# becoming a second source of production or physics policy.


# =============================================================================
# Shell safety
# =============================================================================
#
# -e  Exit when an unhandled command fails.
# -u  Treat use of an unset variable as an error.
# -o pipefail
#     Propagate failures from any command in a pipeline.

set -euo pipefail


# =============================================================================
# Resolve the TreeProduction package directory
# =============================================================================
#
# BASH_SOURCE[0] identifies this launcher even when invoked through another
# working directory. Resolving its physical parent directory gives us a stable
# anchor for locating scripts/produce.py.
#
# CDPATH is cleared for the subshell so a user's shell configuration cannot
# alter the result of cd.

tree_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)


# =============================================================================
# Public command-line interface
# =============================================================================
#
# usage() is intentionally the user-facing vocabulary of TreeProduction.
#
# Public commands describe analysis populations in physicist-readable terms.
# They are translated later into internal lane identifiers.
#
# Options are mostly forwarded unchanged to produce.py. The shell performs only
# lightweight syntax checks; authoritative population and production validation
# belongs to the Python/C++ layers.

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


# =============================================================================
# Public command -> internal production lane translation
# =============================================================================
#
# produce.py works in terms of explicit internal lane names:
#
#   pp_data
#   auau_data
#   pp_photon_sim
#   pp_inclusive_sim
#   pp_photon_di
#   pp_inclusive_di
#   auau_photon_embedded
#   auau_inclusive_embedded
#
# Users should not need to remember that internal vocabulary, so this function
# provides the stable public interface.
#
# Aggregate commands are expanded recursively:
#
#   runAllData -> runPP + runAuAu
#   runAllSim  -> every maintained simulation lane
#   runAll     -> runAllData + runAllSim
#
# The function prints one internal lane per line. The caller captures those
# lines and converts them into the Bash `lanes` array.
#
# Returning status 2 for an unknown command allows the caller to distinguish
# invalid public commands from successful expansion.

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


# =============================================================================
# Handle the trivial help/no-command cases
# =============================================================================
#
# With no arguments there is no production request to resolve, so display the
# interface and exit successfully.
#
# A top-level --help behaves the same way and does not invoke Python.

if [[ $# -eq 0 || $1 == --help ]]; then
    usage
    exit 0
fi


# =============================================================================
# Consume the required public COMMAND
# =============================================================================
#
# The first positional argument is always the public TreeProduction command.
# Remove it from "$@" before processing the remaining options.

command_name=$1
shift


# =============================================================================
# Expand the public command into one or more internal lanes
# =============================================================================
#
# Command substitution captures the newline-separated output of expand_command.
# Failure means the user supplied no recognized public production command.

if ! expanded=$(expand_command "$command_name"); then
    printf 'ERROR: unknown command: %s\n\n' "$command_name" >&2
    usage >&2
    exit 2
fi


# =============================================================================
# Convert the expanded lane list into a Bash array
# =============================================================================
#
# One lane is read per line. Using an array keeps each lane as a distinct
# argument when control is eventually transferred to produce.py.

lanes=()
while IFS= read -r lane; do lanes+=("$lane"); done <<< "$expanded"


# =============================================================================
# Initialize the Python-driver argument vector
# =============================================================================
#
# helper_args begins with the absolute package-relative location of produce.py.
# Every recognized user option below is appended to this array unchanged.
#
# The internal lane names are intentionally kept separate until the final exec.

helper_args=("$tree_dir/scripts/produce.py")


# =============================================================================
# Parse and forward public options
# =============================================================================
#
# This shell layer performs only validation that is cheap and unambiguous:
#
#   - flag-only options are forwarded directly;
#   - integer-valued options must have a non-negative integer token;
#   - value options must actually have a following value;
#   - --source is restricted here to a primitive command;
#   - unknown options fail immediately.
#
# Semantic checks such as whether a run/sample exists, whether a DATA/SIM
# selector is appropriate, whether a catalogue is valid, and whether production
# identities match are performed by produce.py.

while [[ $# -gt 0 ]]; do
    case "$1" in

        # ---------------------------------------------------------------------
        # Help requested after COMMAND.
        # ---------------------------------------------------------------------
        --help) usage; exit 0 ;;

        # ---------------------------------------------------------------------
        # Boolean flags: no following value.
        #
        # --dry-run
        #   Ask produce.py to resolve/check the requested population without
        #   creating output or invoking ROOT.
        #
        # --resume
        #   Permit produce.py to reuse only work whose saved identities and
        #   validation receipts match the current request.
        # ---------------------------------------------------------------------
        --dry-run|--resume) helper_args+=("$1"); shift ;;

        # ---------------------------------------------------------------------
        # Integer-valued options.
        #
        # The shell verifies only lexical integer form here. More restrictive
        # range and contextual checks remain the responsibility of produce.py.
        # ---------------------------------------------------------------------
        --events|--source|--run|--timeout)
            if [[ $# -lt 2 || ! $2 =~ ^[0-9]+$ ]]; then
                printf 'ERROR: %s requires an integer. See --help.\n' "$1" >&2
                exit 2
            fi

            # --source addresses a source ordinal inside one concrete lane.
            # Aggregate commands expand to several lanes and therefore cannot
            # identify an unambiguous source ordinal at this wrapper layer.
            if [[ $1 == --source && ${#lanes[@]} -ne 1 ]]; then
                printf 'ERROR: --source requires a primitive command. See --help.\n' >&2
                exit 2
            fi

            # Preserve the option/value pair exactly for the Python driver.
            helper_args+=("$1" "$2")
            shift 2
            ;;

        # ---------------------------------------------------------------------
        # String/path-valued options.
        #
        # Refuse a missing value and refuse another "--..." token in the value
        # position. Detailed interpretation is deferred to produce.py.
        # ---------------------------------------------------------------------
        --sample|--run-list|--sample-list|--catalogue|--output)
            if [[ $# -lt 2 || -z $2 || $2 == --* ]]; then
                printf 'ERROR: %s requires a value. See --help.\n' "$1" >&2
                exit 2
            fi

            # Preserve the option/value pair exactly for the Python driver.
            helper_args+=("$1" "$2")
            shift 2
            ;;

        # ---------------------------------------------------------------------
        # Anything not explicitly part of the public interface is rejected.
        # ---------------------------------------------------------------------
        *) printf 'ERROR: unknown option: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done


# =============================================================================
# Verify the only launcher-level runtime dependency
# =============================================================================
#
# This wrapper needs Python 3 to enter the production controller. ROOT,
# libPhotonJetTree, libEventTree, the sPHENIX environment, configurations,
# manifests, calibrations, and production identities are validated downstream
# by produce.py and the native production implementation.

if ! command -v python3 >/dev/null 2>&1; then
    printf 'ERROR: python3 is required for population binding.\n' >&2
    exit 1
fi


# =============================================================================
# Transfer control to the production controller
# =============================================================================
#
# `exec` replaces this Bash process with Python rather than creating an extra
# long-lived wrapper process.
#
# The resulting conceptual invocation is:
#
#   python3 scripts/produce.py [forwarded options] <lane> [<lane> ...]
#
# For example:
#
#   ./produce_trees.sh runPhotonJetSim --sample run28_photonjet20
#
# becomes conceptually:
#
#   python3 scripts/produce.py --sample run28_photonjet20 pp_photon_sim
#
# From this point forward, produce.py owns the production lifecycle.

exec python3 "${helper_args[@]}" "${lanes[@]}"
