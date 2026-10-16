#!/bin/sh
# FinalAnalysis entry point. Needs ROOT, libyaml, nlohmann/json.hpp and Python 3 in the selected environment.
#
#   final_analysis.sh build                                     compile the analysis only; print the executable path
#   final_analysis.sh all         --config PIPELINE              everything below in order, into one pinned run folder
#   final_analysis.sh compact     --collection COLLECTION.yaml --out DIRECTORY
#                                                                a delivered collection without the DATA packages'
#                                                                per-source records: the collection to analyse
#   final_analysis.sh qa          --config SETTINGS              the photon and event QA inputs: the QA file named in SETTINGS
#   final_analysis.sh qa          --models MODELS --out FILE     saved training QA, independent of histogram packages
#   final_analysis.sh run         --config SETTINGS              the measurement: the results file named in SETTINGS
#   final_analysis.sh run         --corrections-only SETTINGS    diagnostic correction objects only; no fit or unfolding
#   final_analysis.sh systematics --config SYSTEMATICS           the systematic uncertainties of results of one change each
#   final_analysis.sh draw        [--results FILE] [--qa FILE] [--systematics FILE] --fig DIR
#                                                                the figures of the given files, sPHENIX style, into DIR
#                                                                with the IAN's names
#                               --only combinatoric [--comparisons FILE]
#                                                                Fig27 panels, optionally with pinned comparison recipes
#   final_analysis.sh tables      --results FILE [--systematics FILE] --out DIR
#                                                                the IAN result tables (LaTeX and CSV) into DIR
#
# Everything is compiled on first use into build/ (FINALANALYSIS_BUILD overrides it) and again only where the code
# changed, by content and not by time. build/ also holds saved work that is made again when missing (the simulated
# source rows, file hashes), so deleting build/ is always safe; FINALANALYSIS_CACHE= (empty) saves nothing.
set -eu
package=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build="${FINALANALYSIS_BUILD:-$package/build}"
export FINALANALYSIS_CACHE="${FINALANALYSIS_CACHE-$build/cache}"
usage_draw="usage: final_analysis.sh draw [--results FILE] [--qa FILE] [--systematics FILE] --fig DIR"
usage_tables="usage: final_analysis.sh tables --results FILE [--systematics FILE] --out DIR"

# sha: a short SHA-256 of the standard input
sha() { { sha256sum 2>/dev/null || shasum -a 256; } | cut -c1-20; }

# program NAME OPTIMISATION INCLUDED LIBRARIES SOURCE...: sets $program to build/NAME.<hash>, linked from one object per
# source file, compiled side by side. INCLUDED names a function that prints the project files the sources include.
# Every file in build/ carries in its name a hash of what it is made from (the compiler and ROOT by version, the flags,
# this script, the included files and the source by path and content; the program its objects and libraries) and gets
# that name by one rename when it is complete. So a file that is there is the right one, whatever the times on the
# files say, whoever else builds into the same folder and however an earlier build ended.
# (Only Inputs.cc reads FINALANALYSIS_INPUT_CODE. Its complete local include closure is already covered by that
# object's key, so the saved source-row identity changes with every source or included-header byte change.)
program() {
  name=$1; optimisation=$2; included=$3; libraries=$4; shift 4
  mkdir -p "$build/objects"
  # The saved source rows (Inputs.cc readSources) carry a hash of the code that makes them.
  native_witness=$(python3 "$package/scripts/write_native_build_witness.py" \
    --source-root "$package/.." --output-dir "$build/native-source")
  cxx=$(root-config --cxx); flags="$optimisation -std=c++17 $(root-config --cflags) -include $native_witness"
  input_pins=$(dependencies "$package/src/Inputs.cc")
  input_code=$(printf '%s\n' "$input_pins" | sha | cut -c1-16)
  made_with=$({ "$cxx" --version; root-config --version; printf '%s\n' "$flags"; cat "$0" "$package/dependencies.py"; } | sha)
  key() {
    # Keep failure outside the hash pipeline: a missing include must never acquire a valid cache key.
    pins=$("$included" "$1") || return 1
    { printf '%s\n' "$made_with" "$1" "$pins"; cat "$1"; } | sha
  }
  compilers=""; count=$#
  for source in "$@"; do
    made_from=$(key "$source")
    object="$build/objects/$(basename "$source").$made_from.o"
    if [ ! -f "$object" ]; then
      # An object keeps its name only if the code is the same after the compilation as before it.
      # shellcheck disable=SC2086
      ( $cxx $flags -I"$package/src" -I"$package/plots" -DFINALANALYSIS_INPUT_CODE="\"$input_code\"" -x c++ -c "$source" -o "$object.$$" &&
        [ "$(key "$source")" = "$made_from" ] && mv "$object.$$" "$object" || { rm -f "$object.$$"; exit 1; } ) &
      compilers="$compilers $!"
    fi
    set -- "$@" "$object"
  done
  shift "$count"  # the sources' objects
  failed=""
  for compiler in $compilers; do wait "$compiler" || failed=1; done  # all of them: none is left running
  [ -z "$failed" ] || { echo "FinalAnalysis: compiling $name failed, or its code changed meanwhile" >&2; exit 1; }
  program="$build/$name.$(printf '%s\n' "$made_with" "$libraries" "$@" | sha)"
  if [ ! -x "$program" ]; then
    # shellcheck disable=SC2046,SC2086
    $cxx "$@" -o "$program.$$" $libraries -lyaml
    mv "$program.$$" "$program"
  fi
}

# The compiler owns include resolution, including local headers in the sibling packages.
dependencies() { python3 "$package/dependencies.py" "$cxx" "$flags" "$package" "$@"; }
analysis_headers() { dependencies "$1"; }
compile() {
  program final_analysis -O2 analysis_headers "$(root-config --libs)" "$package"/src/*.cc
}

# The figures and tables: the ROOT macros of plots/, built as ordinary programs. `root 'IANFigures.C+(...)'` draws the
# same files, but its dictionary step takes most of a minute after every change; this takes seconds. Drawing spends
# its time inside ROOT, so the macros are compiled without optimisation.
# Use the compiler include closure, including scientific tracing headers in sibling packages.
plot_files() { dependencies "$1"; }
compile_plots() {
  program figures -O0 plot_files "$(root-config --glibs)" "$package/plots/IANFigures.C"; figures=$program
  program tables -O0 plot_files "$(root-config --glibs)" "$package/plots/IANTables.C"; tables=$program
}

case "${1:-}" in
  build)
    [ "$#" -eq 1 ] || { echo 'usage: final_analysis.sh build' >&2; exit 2; }
    compile
    printf '%s\n' "$program" ;;
  all)
    [ "$#" -eq 3 ] && [ "$2" = "--config" ] || { echo 'usage: final_analysis.sh all --config PIPELINE' >&2; exit 2; }
    compile
    exec "$program" --all "$3" "$package" ;;
  compact)
    [ "$#" -eq 5 ] && [ "$2" = "--collection" ] && [ "$4" = "--out" ] || { echo 'usage: final_analysis.sh compact --collection COLLECTION.yaml --out DIRECTORY' >&2; exit 2; }
    compile
    exec "$program" --compact "$3" "$5" ;;
  qa)
    if [ "$#" -eq 5 ] && [ "$2" = "--models" ] && [ "$4" = "--out" ]; then
      compile
      exec "$program" --model-qa "$3" "$5"
    fi
    [ "$#" -eq 3 ] && [ "$2" = "--config" ] || { echo 'usage: final_analysis.sh qa --config FILE | --models MODELS --out FILE' >&2; exit 2; }
    compile
    exec "$program" --qa "$3" ;;
  run)
    shift; compile
    exec "$program" "$@" ;;
  systematics)
    [ "$#" -eq 3 ] && [ "$2" = "--config" ] || { echo "usage: final_analysis.sh $1 --config FILE" 1>&2; exit 2; }
    compile
    exec "$program" "--$1" "$3" ;;
  draw)
    shift
    results=""; fig=""; qa=""; systematics=""; only=""; comparisons=""
    while [ "$#" -gt 1 ]; do
      case "$1" in
        --results) results=$2 ;;
        --fig) fig=$2 ;;
        --qa) qa=$2 ;;
        --systematics) systematics=$2 ;;
        --only) only=$2 ;;
        --comparisons) comparisons=$2 ;;
        *) echo "$usage_draw" 1>&2; exit 2 ;;
      esac
      shift 2
    done
    [ "$#" -eq 0 ] || { echo "$usage_draw" >&2; exit 2; }
    [ -n "$fig" ] && { [ -n "$results" ] || [ -n "$qa" ]; } || { echo "$usage_draw" 1>&2; exit 2; }
    compile_plots
    "$figures" "$results" "$fig" "$systematics" "$qa" "$only" "$comparisons"
    # Inventory only this draw's explicit roster. Missing upstream evidence is
    # recorded as incomplete; it never grants acceptance or canonical adoption.
    set -- --fig "$fig" --source-root "$package/.." --registry "$package/plot_physics_coverage.json"
    if [ -f "$package/plot_physics_qualification.json" ]; then
      set -- "$@" --qualification "$package/plot_physics_qualification.json"
    fi
    for input in "$results" "$qa" "$systematics"; do
      [ -z "$input" ] || set -- "$@" --input-root "$input"
    done
    python3 "$package/scripts/collect_plot_physics.py" "$@" ;;
  tables)
    shift
    results=""; out=""; systematics=""
    while [ "$#" -gt 1 ]; do
      case "$1" in
        --results) results=$2 ;;
        --out) out=$2 ;;
        --systematics) systematics=$2 ;;
        *) echo "$usage_tables" 1>&2; exit 2 ;;
      esac
      shift 2
    done
    [ "$#" -eq 0 ] || { echo "$usage_tables" >&2; exit 2; }
    [ -n "$results" ] && [ -n "$out" ] || { echo "$usage_tables" 1>&2; exit 2; }
    compile_plots
    exec "$tables" "$results" "$out" "$systematics" ;;
  *)
    sed -n '2,20p' "$0" 1>&2; exit 2 ;;
esac
