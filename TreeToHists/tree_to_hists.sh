#!/bin/bash
# TreeToHists entry point: eventTree v1 files -> one histogram package per sample.
# Needs ROOT (root-config), libyaml and nlohmann/json.hpp in the selected environment.
#
#   tree_to_hists.sh build
#   tree_to_hists.sh plan   --settings YAML --inputs INPUTS.yaml --work DIR [--sample SYSTEM/NAME]...
#   tree_to_hists.sh fill   --work DIR --job N
#   tree_to_hists.sh fill-data-jes --work NOMINAL --up-work UP --down-work DOWN --job N
#   tree_to_hists.sh fill-data-jes-group --work NOMINAL --up-work UP --down-work DOWN --group G0 --scratch DIR
#   tree_to_hists.sh selection-groups --work REF
#   tree_to_hists.sh selection-plan --work REF --out VIEW --variant NAME --prescriptions JSON
#   tree_to_hists.sh fill-selection-group --work REF --view-work VIEW --group G0 --scratch DIR
#   tree_to_hists.sh validate-selection-group --work REF --view-work VIEW --group G0
#   tree_to_hists.sh merge  --work DIR --sample SYSTEM/NAME --out DIR
#   tree_to_hists.sh centrality --settings YAML --inputs INPUTS.yaml --data AuAu/data.root --out weights.json
#   tree_to_hists.sh vertex population --config YAML --out population.root
#   tree_to_hists.sh vertex derive --config YAML --population population.root --out DIR
#   tree_to_hists.sh local  --settings YAML --inputs INPUTS.yaml --work DIR --out DIR [--parallel N] [--sample ...]
#   tree_to_hists.sh campaign --config CAMPAIGN.yaml [--prepare-only|--submit]
#   tree_to_hists.sh campaign-finish --work DIR --out DIR [--parallel N]   (the last job of a flat submission)
#   tree_to_hists.sh collection --work DIR --out DIR
#   tree_to_hists.sh collection --from COLLECTION.yaml [--from COLLECTION.yaml]... --out DIR
#   tree_to_hists.sh isolation population --out population.root EVENTTREE...
#   tree_to_hists.sh isolation derive --out package.json POPULATION...
#   tree_to_hists.sh isolation freeze --out freeze.json PACKAGE [--efficiency 0.8] [--gap 0.8]
#   tree_to_hists.sh isolation check  --out check.json FREEZE POPULATION...
#
# Build once before farm launch. Workers source the pinned setup and use the
# prebuilt executable. The campaign DAG ends with the validated collection.
set -euo pipefail
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
binary="$package/build/tree_to_hists"

build() {
  mkdir -p "$package/build"
  if [ ! -x "$package/src/configure" ] || [ "$package/src/configure.ac" -nt "$package/src/configure" ] || [ "$package/src/Makefile.am" -nt "$package/src/Makefile.in" ]; then
    "$package/src/autogen.sh"
  fi
  (cd "$package/build"; "$package/src/configure" CXX="$(root-config --cxx)"; make -j "${TREETOHISTS_BUILD_JOBS:-2}")
}
current() {  # build when a source is newer than the binary (never inside a fill job)
  if [ ! -x "$binary" ] || [ "$package/../TreeProduction/src/internal/Configuration.h" -nt "$binary" ] || [ -n "$(find "$package/src" -type f -newer "$binary")" ]; then build; fi
}
value() {  # value NAME ARGS...: the argument after NAME
  local name=$1; shift
  while [ $# -gt 1 ]; do [ "$1" = "$name" ] && { echo "$2"; return; }; shift; done
  echo "tree_to_hists.sh: missing $name" >&2; exit 2
}

command=${1:-}; [ $# -gt 0 ] && shift
case "$command" in
  selection-groups|selection-plan|fill-selection-group|validate-selection-group|fill|fill-data-jes|fill-data-jes-group|fill-group|merge-group|campaign-inspect|campaign-plan|campaign-finish|collection)
    if [ -n "${TREETOHISTS_SETUP:-}" ]; then set +u; source ${TREETOHISTS_SETUP}; set -u; fi ;;
esac
case "$command" in
  build) build ;;
  plan) current; exec "$binary" plan --home "$package" "$@" ;;
  fill|fill-data-jes|fill-data-jes-group)
    [ -x "$binary" ] || { echo "tree_to_hists.sh: no build/tree_to_hists; run tree_to_hists.sh build first" >&2; exit 1; }
    exec "$binary" "$command" "$@" ;;
  merge) current; exec "$binary" merge "$@" ;;
  centrality) current; exec "$binary" centrality --home "$package" "$@" ;;
  vertex) current; exec "$binary" vertex "$@" --home "$package" ;;
  isolation) current; exec "$binary" isolation "$@" ;;
  campaign)
    current
    launch=0; args=()
    for item in "$@"; do
      case "$item" in --submit) launch=1 ;; --prepare-only) ;; *) args+=("$item") ;; esac
    done
    "$binary" campaign --home "$package" "${args[@]}"
    cfg=$(value --config "${args[@]}")
    # The native command reports the frozen directory. Read the YAML with the
    # same driver rather than introducing a second configuration parser.
    if [ "$launch" -eq 1 ]; then
      launch_work=$("$binary" campaign-submit-check --config "$cfg" --home "$package")
      (cd "$launch_work/condor"; condor_submit_dag campaign.dag)
    fi ;;
  selection-groups|selection-plan|fill-selection-group|validate-selection-group|campaign-inspect|campaign-plan|merge-group|collection)
    [ -x "$binary" ] || { echo 'TreeToHists must be built before launching farm workers' >&2; exit 1; }
    exec "$binary" "$command" "$@" ;;
  campaign-finish)
    # Where the campaign graph cannot go to DAGMan: the fill groups are queued as plain jobs together with this one,
    # which waits for them, then merges and writes the validated collection (campaign_finish.py).
    [ -x "$binary" ] || { echo 'TreeToHists must be built before launching farm workers' >&2; exit 1; }
    exec python3 "$package/campaign_finish.py" "$binary" "$@" ;;
  local)
    current
    work=$(value --work "$@"); out=$(value --out "$@"); parallel=1
    args=(); while [ $# -gt 0 ]; do
      case "$1" in --out) shift 2 ;; --parallel) parallel=$2; shift 2 ;; *) args+=("$1"); shift ;; esac
    done
    "$binary" plan --home "$package" "${args[@]}"
    if [ -s "$work/jobs.txt" ]; then
      cut -d' ' -f1 "$work/jobs.txt" | xargs -P "$parallel" -I{} "$binary" fill --work "$work" --job {}
    fi
    while IFS= read -r sample; do "$binary" merge --work "$work" --sample "$sample" --out "$out"; done < "$work/samples.txt"
    "$binary" collection --work "$work" --out "$out" ;;
  *) sed -n '2,28p' "$0" >&2; exit 2 ;;
esac
