#!/bin/sh
# FinalAnalysis entry point. Needs ROOT, libyaml and nlohmann/json.hpp in the selected environment.
#
#   final_analysis.sh run         --config SETTINGS              the measurement: the results file named in SETTINGS
#   final_analysis.sh qa          --config SETTINGS              the photon and event QA inputs: the QA file named in SETTINGS
#   final_analysis.sh systematics --config SYSTEMATICS           the systematic uncertainties of results of one change each
#   final_analysis.sh draw        [--results FILE] [--qa FILE] [--systematics FILE] --fig DIR
#                                                                the figures of the given files, sPHENIX style, into DIR
#                                                                with the IAN's names
#   final_analysis.sh tables      --results FILE [--systematics FILE] --out DIR
#                                                                the IAN result tables (LaTeX and CSV) into DIR
#   final_analysis.sh all         --config PIPELINE               reduction through figures, tables and a pinned run record
#
# The C++ is compiled on first use into build/ (FINALANALYSIS_BUILD overrides it); the figures and tables are ROOT
# macros compiled there too.
set -eu
package=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build="${FINALANALYSIS_BUILD:-$package/build}"
usage_draw="usage: final_analysis.sh draw [--results FILE] [--qa FILE] [--systematics FILE] --fig DIR"
usage_tables="usage: final_analysis.sh tables --results FILE [--systematics FILE] --out DIR"
compile() {
  mkdir -p "$build"
  if [ ! -x "$build/final_analysis" ] || [ "$package/../TreeToHists/src/Configuration.h" -nt "$build/final_analysis" ] || \
     [ -n "$(find "$package/src" -type f -newer "$build/final_analysis")" ]; then
    # shellcheck disable=SC2046
    $(root-config --cxx) -O2 -std=c++17 $(root-config --cflags) -I"$package/src" "$package"/src/*.cc \
      -o "$build/final_analysis" $(root-config --libs) -lyaml
  fi
}
cpp_string() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
case "${1:-}" in
  all)
    [ "$#" -eq 3 ] && [ "$2" = "--config" ] || { echo 'usage: final_analysis.sh all --config PIPELINE' >&2; exit 2; }
    compile
    exec "$build/final_analysis" --all "$3" "$package" ;;
  run)
    shift; compile
    exec "$build/final_analysis" "$@" ;;
  systematics|qa)
    [ "$#" -eq 3 ] && [ "$2" = "--config" ] || { echo "usage: final_analysis.sh $1 --config FILE" 1>&2; exit 2; }
    compile
    exec "$build/final_analysis" "--$1" "$3" ;;
  draw)
    shift
    results=""; fig=""; qa=""; systematics=""
    while [ "$#" -gt 1 ]; do
      case "$1" in
        --results) results=$2 ;;
        --fig) fig=$2 ;;
        --qa) qa=$2 ;;
        --systematics) systematics=$2 ;;
        *) echo "$usage_draw" 1>&2; exit 2 ;;
      esac
      shift 2
    done
    [ "$#" -eq 0 ] || { echo "$usage_draw" >&2; exit 2; }
    [ -n "$fig" ] && { [ -n "$results" ] || [ -n "$qa" ]; } || { echo "$usage_draw" 1>&2; exit 2; }
    mkdir -p "$build"
    exec root -l -b -q -e "gSystem->SetBuildDir(\"$(cpp_string "$build")\", true); gSystem->AddIncludePath(\"-I$(cpp_string "$package")/plots\");" \
      "$package/plots/IANFigures.C+(\"$(cpp_string "$results")\", \"$(cpp_string "$fig")\", \"$(cpp_string "$systematics")\", \"$(cpp_string "$qa")\")" ;;
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
    mkdir -p "$build"
    exec root -l -b -q -e "gSystem->SetBuildDir(\"$(cpp_string "$build")\", true);" \
      "$package/plots/IANTables.C+(\"$(cpp_string "$results")\", \"$(cpp_string "$out")\", \"$(cpp_string "$systematics")\")" ;;
  *)
    sed -n '2,14p' "$0" 1>&2; exit 2 ;;
esac
