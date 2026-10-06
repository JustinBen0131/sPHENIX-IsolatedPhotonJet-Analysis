#!/bin/bash
# TreeToHists entry point: eventTree v1 files -> one histogram package per sample.
# Needs ROOT (root-config), libyaml and nlohmann/json.hpp in the selected environment.
#
#   tree_to_hists.sh build
#   tree_to_hists.sh plan   --settings YAML --inputs INPUTS.yaml --work DIR [--sample SYSTEM/NAME]...
#   tree_to_hists.sh fill   --work DIR --job N
#   tree_to_hists.sh merge  --work DIR --sample SYSTEM/NAME --out DIR
#   tree_to_hists.sh centrality --settings YAML --inputs INPUTS.yaml --data AuAu/data.root --out weights.json
#   tree_to_hists.sh local  --settings YAML --inputs INPUTS.yaml --work DIR --out DIR [--parallel N] [--sample ...]
#   tree_to_hists.sh condor --work DIR --logs DIR --events-log /tmp/$USER/<campaign>/events.log --setup "SCRIPT [ARGS]"
#   tree_to_hists.sh isolation population --out population.root EVENTTREE...
#   tree_to_hists.sh isolation derive --out package.json POPULATION...
#   tree_to_hists.sh isolation freeze --out freeze.json PACKAGE [--efficiency 0.8] [--gap 0.8]
#   tree_to_hists.sh isolation check  --out check.json FREEZE POPULATION...
#
# build compiles the driver once into build/; a fill job never compiles. local plans, fills and merges on this machine.
# condor writes WORK/condor/fill.sub for the planned fill jobs and prints the command that submits it, following the
# SDCC rules: no inherited environment (the job sources SETUP), the Condor event log in /tmp/$USER/<campaign>, stdout
# and stderr persistent outside TG and sharded below 1000 files per directory, no file transfer, no runtime limit.
# After the jobs, merge each sample (the printed commands). centrality derives the Au+Au simulation centrality weights
# from the DATA package, before the Au+Au simulation is planned. isolation derives, freezes and checks the Au+Au
# isolation package from the truth-signal photons of embedded simulation (src/Isolation.cc).
set -euo pipefail
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
binary="$package/build/tree_to_hists"
sources=("$package/src/TreeToHists.cc" "$package/src/Isolation.cc")

build() {
  mkdir -p "$package/build"
  # shellcheck disable=SC2046
  "$(root-config --cxx)" -O2 -std=c++17 $(root-config --cflags) -I"$package/src" "${sources[@]}" \
    -o "$binary" $(root-config --libs) -lTreePlayer -lyaml
}
current() {  # build when a source is newer than the binary (never inside a fill job)
  if [ ! -x "$binary" ] || [ -n "$(find "$package/src" -type f -newer "$binary")" ]; then build; fi
}
value() {  # value NAME ARGS...: the argument after NAME
  local name=$1; shift
  while [ $# -gt 1 ]; do [ "$1" = "$name" ] && { echo "$2"; return; }; shift; done
  echo "tree_to_hists.sh: missing $name" >&2; exit 2
}
submit_description() {  # the Condor description of the planned fill jobs
  cat <<SUB
# TreeToHists fill jobs of $work: one per eventTree file (jobs.txt).
universe              = vanilla
executable            = $package/tree_to_hists.sh
arguments             = fill --work $work --job \$(Process)
getenv                = False
environment           = "TREETOHISTS_SETUP='$setup'"
shard                 = \$INT(Process / 1000, %03d)
log                   = $events
output                = $logs/\$(shard)/fill_\$(Process).out
error                 = $logs/\$(shard)/fill_\$(Process).err
should_transfer_files = NO
request_memory        = 4GB
queue $jobs
SUB
}

command=${1:-}; [ $# -gt 0 ] && shift
case "$command" in
  build) build ;;
  plan) current; exec "$binary" plan --home "$package" "$@" ;;
  fill)
    if [ -n "${TREETOHISTS_SETUP:-}" ]; then set +u; source ${TREETOHISTS_SETUP}; set -u; fi
    [ -x "$binary" ] || { echo "tree_to_hists.sh: no build/tree_to_hists; run tree_to_hists.sh build first" >&2; exit 1; }
    exec "$binary" fill "$@" ;;
  merge) current; exec "$binary" merge "$@" ;;
  centrality) current; exec "$binary" centrality --home "$package" "$@" ;;
  isolation) current; exec "$binary" isolation "$@" ;;
  local)
    current
    work=$(value --work "$@"); out=$(value --out "$@"); parallel=1
    args=(); while [ $# -gt 0 ]; do
      case "$1" in --out) shift 2 ;; --parallel) parallel=$2; shift 2 ;; *) args+=("$1"); shift ;; esac
    done
    "$binary" plan --home "$package" "${args[@]}"
    cut -d' ' -f1 "$work/jobs.txt" | xargs -P "$parallel" -I{} "$binary" fill --work "$work" --job {}
    for sample in $(cut -d' ' -f2 "$work/jobs.txt" | sort -u); do "$binary" merge --work "$work" --sample "$sample" --out "$out"; done ;;
  condor)
    current
    work=$(cd -- "$(value --work "$@")" && pwd); logs=$(value --logs "$@"); events=$(value --events-log "$@"); setup=$(value --setup "$@")
    case "$logs" in /sphenix/tg*|/sphenix/lustre*) echo "tree_to_hists.sh: stdout/stderr stay outside TG/Lustre (SDCC storage tiers)" >&2; exit 2 ;; esac
    case "$events" in "/tmp/${USER}/"*) ;; *) echo "tree_to_hists.sh: the Condor event log goes in /tmp/\$USER/<campaign>/ (SDCC rule)" >&2; exit 2 ;; esac
    jobs=$(wc -l < "$work/jobs.txt"); shards=$(( (jobs + 999) / 1000 ))
    for ((s = 0; s < shards; ++s)); do mkdir -p "$logs/$(printf %03d $s)"; done
    mkdir -p "$work/condor"
    submit_description > "$work/condor/fill.sub"
    echo "Wrote $work/condor/fill.sub ($jobs fill jobs). Submit it from the submit host, after mkdir -p $(dirname "$events"):"
    echo "  condor_submit $work/condor/fill.sub"
    echo "When every job has finished (a done.json in each job directory), merge each sample:"
    for sample in $(cut -d' ' -f2 "$work/jobs.txt" | sort -u); do echo "  $package/tree_to_hists.sh merge --work $work --sample $sample --out OUT"; done ;;
  *) sed -n '2,23p' "$0" >&2; exit 2 ;;
esac
