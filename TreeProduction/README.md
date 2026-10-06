# Photon–jet tree production

TreeProduction turns sPHENIX DSTs into analysis trees for the isolated-photon +
jet measurement. One Fun4All producer, `PhotonJetTree`, serves pp and Au+Au,
DATA and simulation, including archived pp double-interaction (DI) samples and
Au+Au embedded simulation. Every logical DST source becomes exactly one file
holding one schema-bound `eventTree`
([branch reference](SCHEMA.md), [contract](config/eventTree_v1.yaml)).

The reconstruction follows the established sPHENIX RecoilJets definitions
and writes the retained analysis quantities directly to the eventTree layout.

## Build

Use the sPHENIX release named in `config/tree_production.yaml` (`ana.568`),
with ROOT, FastJet, HepMC, CLHEP, Autotools and Python 3 with PyYAML (the build
expands the branch schema). `OFFLINE_MAIN` identifies the installed release;
`MYINSTALL` is your installation prefix. From `TreeProduction`:

```bash
mkdir -p build
cd build
../src/autogen.sh --prefix="$MYINSTALL"
make -j4
make install
cd ..
./produce_trees.sh --help
```

This builds and installs `libEventTree.so` (the saved-file library, ROOT only)
and `libPhotonJetTree.so` (the producer, linked to it), with their headers
under `include/TreeProduction/`. Reading, checking or extending saved trees on a
machine without the sPHENIX release needs only the first:

```bash
../src/autogen.sh --prefix="$MYINSTALL" --disable-producer
make -j4 && make install
```

The runner needs Python 3 with ROOT's Python bindings and PyYAML. The producer
also reads the frozen photon-ID package in `../PhotonID` and the Au+Au isolation
freeze in `../TreeToHists/config`, both pinned by digest in the configuration.

## Running

| Command | Collection | Configuration profile |
|---|---|---|
| `runPP` | pp DATA | `pp_data` |
| `runAuAu` | Au+Au DATA | `auau_data` |
| `runPhotonJetSim` | pp photon+jet simulation | `pp_photon_sim` |
| `runInclusiveJetSim` | pp inclusive-jet simulation | `pp_inclusive_sim` |
| `runPhotonJetDISim` | pp double-interaction photon+jet simulation | `pp_photon_di` |
| `runInclusiveJetDISim` | pp double-interaction inclusive-jet simulation | `pp_inclusive_di` |
| `runEmbeddedPhotonJetSim` | Au+Au embedded photon+jet simulation | `auau_photon_embedded` |
| `runEmbeddedInclusiveJetSim` | Au+Au embedded inclusive-jet simulation | `auau_inclusive_embedded` |
| `runAllData`, `runAllSim`, `runAll` | the DATA, the simulation, or every collection | |

```bash
./produce_trees.sh runAll --dry-run                          # check every binding; run nothing
./produce_trees.sh runPP --run 51109                         # one DATA run of the catalogue
./produce_trees.sh runPP --run 51109 --resume                # reuse what is complete and verified
./produce_trees.sh runPhotonJetDISim --sample run28_photonjet20_double --source 2 --events 10
```

| Option | Meaning |
|---|---|
| `--run N` | one DATA run |
| `--run-list FILE` | DATA runs from a plain list, one run number per line |
| `--sample ID` | one exact simulation sample |
| `--sample-list FILE` | SIM samples from a plain list, one exact name per line |
| `--source N` | one source by its ordinal (with `--run` or `--sample` when several groups use N) |
| `--events N` | at most N events per selected source |
| `--catalogue FILE` | a prepared catalogue instead of `config/production_populations.json` |
| `--output DIR` | output root; existing products are refused |
| `--resume` | reuse completed pieces and sources whose identities match; refuse anything else |
| `--timeout N` | seconds per reconstruction piece (at most 43,200) |
| `--dry-run` | check bindings and print the plan; writes nothing and runs no ROOT |

A selection with `--source` or `--events` is a development run: it writes below
`development/` and is never published as a finished group. An aggregate command
checks every selected population before the first event loop.

## Input lists

Ordinary users select a collection with `produce_trees.sh`; input lists and
reconstruction settings are chosen automatically. The distributed inputs cover
all 1,568 pp and 2,459 Au+Au good-run-list runs, 14 ordinary simulation samples
and eight archived double-interaction samples: 891,012 complete DST sources.

```text
lists/
  pp_data/
    runs.list
    <run>/DST_JETCALO.list
    <run>/DST_Jet.list
    <run>/sources.tsv
  auau_data/
    runs.list
    <run>/DST_JETCALO.list
    <run>/DST_Jet.list
    <run>/sources.tsv
  simulation/
    samples.list
    <sample>/<DST-stream>.list
    <sample>/sources.tsv
```

A DST `.list` contains one filename per line. The lists of a run or sample
have the same row order: row N in every stream belongs to the same source.
`NONE` marks a stream absent from an embedded source. Keep paired streams
aligned; the runner checks their exact contents before invoking Fun4All.

To choose runs, create a plain file containing one run number per line:

```bash
./produce_trees.sh runPP --run-list my_pp_runs.list
./produce_trees.sh runAuAu --run-list my_auau_runs.list
./produce_trees.sh runAllSim --sample-list my_sim_samples.list
```

A sample list contains exact sample names, one per line, such as
`run28_photonjet20` or `run28_jet20_double`. Blank lines and `#` comments are
allowed. Unknown or duplicate identifiers refuse the command before execution.
Without a selector, the command processes its full maintained population.

The `sources.tsv` sidecars carry source identities, event counts, archived DI
seeds and Au+Au calibration bindings. The runner loads these automatically;
users do not assemble metadata or edit JSON to select their inputs.
`config/populations/*.tsv` pins the source tables, and the small
`config/production_populations.json` binds the lane profiles and output root.
Existing custom JSON catalogues remain supported through `--catalogue`.
Input coverage does not imply completed or scientifically accepted output trees.

## Configuration

`config/tree_production.yaml` holds the reconstruction presets; ordinary users
select one through the command. Its sections, in order: `contract` and
`runtime` (the release, and the PPG12 relative event weights: slices, samples
and periods); `capture` (the event capture policy); `nodes`; `calorimeter` and
`photons` (photon builder, named models, preselection, the seven shower views,
isolation); `pp` and `auau` (event classification, isolation witnesses and
nominal isolation, photon identification, centrality, CEMC status,
underlying-event subtraction, jet views, the established producer's
diagnostics); `jets` (calibration and matching); `profiles` (the eight
production profiles); and `output` (the optional per-object tables). Every
external payload is pinned by SHA-256 and verified before the event loop; a
missing payload is an error, never a default.

## Output

Each source is processed in bounded pieces (35,000 events for Au+Au DATA,
100,000 otherwise). A piece that starts inside its source opens every input
stream at its first entry, and each event is checked to be at the source entry
it is recorded as in every stream, so a misplaced input stops the piece rather
than relabelling events. Every piece is checked where it was written
(completion record, provenance, source interval, event accounting, structure
and the photon-ID inputs) and the pieces are assembled into

    <output>/<campaign>/<lane>/<run-or-sample>/part-<ordinal>.root

Assembly preserves event order, event-local references, per-piece metadata and
objects, and reversible code dictionaries; cumulative run trigger tables stay
per piece, and only the additive simulation centrality denominator is summed.
The candidate is reread against every piece before it is published. Pieces,
logs and receipts stay under the group's `.work/` directory.

When a declared population finishes (a lane, or one `--run` or `--sample`), each
group whose sources all completed is published as a version:
`versions/<version>/` hard-links the files beside a `MANIFEST.json` that binds
each file's digest and receipt, and `CURRENT` points to it. A published group is
never replaced. Receipts bind the input, configuration, macro, library and
output digests; the producer also records the loaded library, the resolved
configuration and the compiled source digest in every file.

## What a tree holds

One entry per event handed to the producer, in source order, including events
without a photon (the exposure denominator); see [SCHEMA.md](SCHEMA.md) for
every branch.

- **Event.** Identity and source coordinates, GL1 trigger words and the run
  trigger table (`meta/trigger_run_info`), vertex, MBD, centrality, calorimeter
  energy sums and the producer-owned weights. pp simulation stores the relative
  PPG12 weight `slice x vertex x SI/DI fraction x period share` in `w_event`
  with its factors (`w_xsec` and `w_exposure` are aliases, not extra factors);
  absolute normalization and the Au+Au vertex and centrality weights are applied
  downstream.
- **Event classification.** `evt_terminal_status` records the established
  producer's event classification and never removes anything. pp: 0 accepted,
  1 trigger, 2 |vz|, 3 period, 4 sample ownership. Au+Au: 0 accepted,
  1 trigger, 2 |vz|, 3 minimum bias, 4 generator-slice stitch, 5 centrality.
  The reader predicate is
  `evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && evt_terminal_status == 0`.
- **Photons.** Every candidate of the bundled `photonjet::PhotonClusterBuilder`
  inside the capture domain (5 <= ET < 40 GeV, |eta| < 0.7; Au+Au DATA from
  12 GeV, as its established production), with its kinematics, seven shower-shape views
  (`H70`, `H0`, `G70`, `G0`, `O70`, `O0`, `R70`), isolation cones (pp
  topological clusters, Au+Au underlying-event subtracted towers) with their
  historical working points (`isow_*`), cells and cone constituents, the named
  model scores (PPG12 NPB, and the PPG12 classifier in pp) and the
  preselection. Every candidate carries the frozen nominal labels
  (`pho_bdt_score`, `pho_bdt_state`, `pho_id_class`, `pho_iso_class`) of
  `../PhotonID/models/nominal.json` and the pinned Au+Au isolation freeze; no
  candidate is removed by them.
- **Jets.** Anti-kT R = 0.2, 0.3 and 0.4, raw and calibrated (the JES applied
  once), and their EMCal, inner-HCal and outer-HCal energy fractions. Au+Au
  keeps the underlying-event subtracted view (`jet*`) and the unsubtracted view
  (`jetunsub*`). Photon-jet pairs carry `xJ` and the azimuthal separation.
- **Truth.** In simulation: truth photons with their ancestry and isolation,
  truth jets, truth vertices, and the photon and jet reconstruction-to-truth
  links.
- **The established producer's diagnostics.** pp: its PPG12 photon-yield
  bookkeeping, the audit occurrences of the fills it checked (`pdiag_*`) and, in
  simulation, its weighted and raw fill counters (`ppg12_wfill*`,
  `ppg12_rfill*`). Au+Au embedded: the generator-slice stitch decision and the
  truth-isolation audit (`embdiag_*`), and the photon-count unfolding's truth
  misses (`truthmiss_*`): a truth signal photon with no selected isolated tight
  photon, or whose selected photon has another dominant truth primary.

Reconstruction follows the established producer exactly. In particular:

- **Archived DI.** The producer replays, first in first out, the five random
  seeds the archived production consumed (catalogue `historical_rng`): the MBD
  digitization, the pedestal file draw, and the CEMC, inner and outer HCal
  waveform simulations over that pedestal file, followed by the ordinary
  calibration chain under the run-28 conditions. A source whose seeds choose a
  different pedestal file than recorded is refused.
- **Au+Au DATA.** Centrality is reconstructed from the run's own divisions, run
  scale and vertex scale in the frozen calibration set, each verified against
  the set's pinned `MANIFEST.json`; an event has one only when the minimum-bias
  classifier passed (`cent_value` is NaN otherwise). The CEMC status of the paired input is
  checked against the bad-tower map its DST was reconstructed with; a tower the
  map rejects but whose saved status lacks the hot flag is restored only from
  the run's own map bound by the catalogue (`cemc_recovery_map`). Underlying-
  event subtraction runs in two passes before isolation and the subtracted jets.
- **Au+Au embedded.** Conditions come from the data run the simulation is
  embedded in, which the DST name carries (`...-data-<run>-...`); the source
  keeps its simulation identity (run 28).

## Reading the trees

`examples/quick_start.py` (uproot) and `examples/quick_start.C` (TTreeReader)
apply the reader predicate, read photons, jets and pairs, and histogram the
back-to-back photon + R = 0.4 jet `xJ`:

```bash
python3 examples/quick_start.py part-00070661.root
root -l -b -q 'examples/quick_start.C("part-00070661.root")'
```

`EventTree.h` documents the saved-file library: after
`ROOT.gSystem.Load("libEventTree")` and including `TreeProduction/EventTree.h`,
`ROOT.eventtree.checkFile(path, schema, "")` returns a JSON report of one file.
The library refuses a schema other than the one it was built from.

## Compatibility

Jet energy fractions are stored for R = 0.2, 0.3 and 0.4; their numerical
equivalence to the established downstream calculation remains unqualified.
FastJet area ghosts depend on the event history within a reconstruction job.
Splitting an input source can therefore change jet areas and the final digits
of jet and photon–jet kinematics relative to a differently split production.

## Source map

| File | Responsibility |
|---|---|
| `produce_trees.sh` | Public commands and option forwarding |
| `scripts/produce.py` | Catalogue checks, execution, receipts, resume and group publication; calls the eventTree library |
| `scripts/check_environment.py` | Compare the active Python environment with the repository's `environment.yml` |
| `macros/Fun4All_PhotonJetTree.C` | Resolve the plan, configure conditions and reconstruction, register the producer and inputs, run |
| `config/tree_production.yaml` | Reconstruction presets and production profiles |
| `config/eventTree_v1.yaml` | The eventTree branch contract ([readable reference](SCHEMA.md)) |
| `config/production_populations.json`, `config/populations/`, `lists/` | Full GRL and simulation source populations |
| `src/configure.ac`, `src/Makefile.am`, `src/autogen.sh` | Autotools build of the two libraries and their headers |
| `src/PhotonJetTree.h/.cc` | Producer lifecycle and the event-processing order |
| `src/PhotonClusterBuilder.h/.cc` | Photon reconstruction: kinematics, shower shapes, isolation |
| `src/EventTree.h/.cc` | Saved-file operations (ROOT only): schema, writing, checks, assembly, branch additions, publication |
| `src/internal/Production.h/.cc` | Job, plan and conditions: the configuration resolved and validated |
| `src/internal/Setup.cc` | Conditions, reconstruction modules, the input streams and their per-event entry check |
| `src/internal/Types.h` | Configuration and record definitions, units and missing-value states |
| `src/internal/Configuration.h` | Read-only access to the parsed configuration |
| `src/internal/Event.cc` | Event context, triggers, weights and the event classification |
| `src/internal/Objects.cc` | Photons, jets, pairs, nominal labels and the PPG12 source arbitration |
| `src/internal/Truth.cc` | Simulation truth, associations and the Au+Au embedded truth diagnostics |
| `src/internal/PPG12.h/.cc`, `PPG12Simulation.cc` | The established producer's PPG12 bookkeeping |
| `src/internal/Output.cc` | The eventTree writer: file metadata, one entry per event, completion |
| `examples/` | Python and C++ readers |
