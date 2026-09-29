# Photon–jet tree production

One Fun4All producer for pp and AuAu, DATA and simulation. The public command
selects an input collection and reconstruction preset; users do not assemble
standard input lists or edit detector settings.

**Current qualification:** the package layout, frontend, bundled photon builder
and native fraction arithmetic have local checks. A full sPHENIX build and
source-matched reproduction of the final single `eventTree` remain required.
The C++ producer writes one schema-bound `eventTree` directly. The runner
checks its source coverage and calls the shared nominal finalizer. Several required capture fields
are still unresolved; execution refuses known mapping gaps before reading DSTs.
Inspection candidates are not substitutes for published, qualified products.

## Build

Use the configured sPHENIX release, ROOT, FastJet, yaml-cpp, HepMC, CLHEP,
Autotools and Python 3. `OFFLINE_MAIN` identifies the installed release;
`MYINSTALL` is your local installation prefix. Python execution additionally
requires numpy, awkward, uproot and the bundled eventTree module. No other
package's build is replaced.

From `TreeProduction`:

```bash
mkdir -p build
cd build
../src/autogen.sh --prefix="$MYINSTALL"
make -j4
make install
cd ..
./produce_trees.sh --help
```

The three maintained build files live in `src/`. They build one
`libPhotonJetTree.so` and install headers under `include/TreeProduction/`.
The runner uses `$MYINSTALL/lib/libPhotonJetTree.so` unless the catalogue pins
another installation. It supplies that library and its headers to the ROOT
child. Missing dependencies stop configuration with an error.

`photonjet::PhotonClusterBuilder` is compiled into this library from the one
bundled source pair in `src/`. Its namespace prevents collisions with the
installed `::PhotonClusterBuilder`; all other reconstruction modules use the
normal coresoftware installation. No `LD_PRELOAD`, CaloReco replacement or
coresoftware checkout is required. Generic builder defaults and optional
model interfaces remain separate from the selected production preset.

## Commands

| Command | Input collection |
|---|---|
| `runPP` | pp DATA |
| `runAuAu` | AuAu DATA |
| `runPhotonJetSim` | pp photon–jet simulation |
| `runInclusiveJetSim` | pp inclusive-jet simulation |
| `runPhotonJetDISim` | pp double-interaction photon–jet simulation |
| `runInclusiveJetDISim` | pp double-interaction inclusive-jet simulation |
| `runEmbeddedPhotonJetSim` | AuAu embedded photon–jet simulation |
| `runEmbeddedInclusiveJetSim` | AuAu embedded inclusive-jet simulation |
| `runAllData` | Both DATA collections |
| `runAllSim` | All six SIM collections, including DI |
| `runAll` | Every DATA and SIM collection |

DI commands are exposed explicitly; the archived-input reconstruction route
still refuses execution until its required replay implementation is qualified.
An aggregate validates every selected population before starting any event loop.

```bash
./produce_trees.sh runPP --dry-run
./produce_trees.sh runAuAu --run 53285 --dry-run
./produce_trees.sh runPhotonJetSim --sample photon5 --dry-run
./produce_trees.sh runAuAu --source 0 --events 10 --dry-run
```

`--run` selects an exact DATA run; `--sample` selects an exact SIM sample.
`--source` is the stable manifest ordinal and requires a primitive command.
`--events` caps the count per selected source; zero and negative counts fail.
Source/event restrictions receive a separate development output path.
`--resume` reuses only source- and runtime-matched validated work.
`--timeout N` bounds each reconstruction process to at most 43,200 seconds.
`--output DIR` changes the output root. An advanced user may supply an already
prepared catalogue with `--catalogue FILE`. `--dry-run` writes nothing and does
not run ROOT; it checks operational bindings, not detector reconstruction.

pp SIM saves the established relative PPG12 weight
`slice × vertex × SI/DI fraction × period share` in `w_event`, with the
individual factors retained. The period selects its pinned truth-vertex curve;
DI uses both independent truth vertices. `w_xsec` aliases `w_slice` and
`w_exposure` aliases `w_period`: these are not extra factors. Absolute sample
normalization remains downstream. AuAu trees retain the vertices and unit
producer vertex factor; TreeToHists applies its separately qualified AuAu
vertex and centrality fits. AuAu SIM also carries the unweighted, pre-analysis
vertex-cut centrality denominator at `hist/SIM/h_centrality`.

## Maintained inputs and configuration

- `config/production_populations.json` binds each command to a profile, a
  source manifest and its SHA256. The current distributed bindings are
  explicitly unbound pending the qualified collection release; commands refuse
  them rather than discovering arbitrary files.
- `config/tree_production.yaml` owns reconstruction, nodes, calibrations,
  jet radii and capture settings. Ordinary users select a preset through the
  command. Missing required payloads remain errors.

A version-1 population manifest contains `lane`, `campaign` and `sources`.
Each source supplies `source_file_ordinal`, `input_list`, `input_list_sha256`,
`run`, `segment`, `first_entry`, `event_count`, `sample`, `period` and `si_di_role`.
`logical_source_id` preserves the original logical source identity;
`input_file_sha256` may additionally pin its input payload.
An optional per-source `production_config` requires its SHA256. Relative paths
resolve against the containing catalogue/manifest. The one-bundle input-list
layout is interpreted only by `production::LoadPlan` and `Setup.cc`.
Duplicate identities, changed pins, missing metadata and existing outputs fail.

The compiled source digest, actually loaded library digest, steering macro and
configuration digest are recorded automatically. Optional non-null runtime
pins in YAML must match those observed values. Editing installed sources does
not relabel an older library as rebuilt.

## Output contract and qualification

The delivery target is exactly one `eventTree` per logical DST source, grouped
by DATA run or exact SIM sample. Its authoritative branch definitions live in
[eventTree schema](SCHEMA.md), shared writer and maintained readers.
Nominal labels use the frozen AuAu and PPG12 photon ID; AuAu isolation uses the
frozen R=0.3 80% fit and pp uses PPG12. Candidate retention is independent of
those labels. Later model changes are versioned augmentations.

The new jet reconstruction captures EMCal, inner-HCal and outer-HCal energy
fractions for the analysis jets at R=0.2, 0.3 and 0.4. Fractions use signed
constituent energies and the original raw jet energy, with float arithmetic
before widening to double. Undefined fractions remain NaN. Capture reads the
exact current-event raw jet; it never matches to a nearby calibrated jet.
Farm replay and final saved-file regression are still required.

Each source is processed in bounded pieces (35,000 events for AuAu DATA;
100,000 for other profiles), finalized and reread, then assembled into
`<campaign>/<lane>/<run-or-sample>/part-<ordinal>.root`. Piece files and logs
remain below the group's ignored `.work/` directory. Consolidation preserves
event order, per-piece metadata and reversible code dictionaries. It proves
source coverage before publishing one final file.

Receipts bind input, configuration, macro, library, runtime and output hashes.
An explicit `--resume` can reuse matching completed pieces or a completed source;
partial or changed work is retained and refused. Publication never replaces a
different existing file. Source-matched reconstruction regression is still
required before replacing the production baseline.

The maintained catalogue pins the AuAu isolation freeze by path and SHA256;
the shared PhotonID package owns the frozen BDT models. Reconstruction captures
all seven shower views through one configurable builder calculation; each
view's definition and semantic digest live in one table (`internal/Types.h`).
View arithmetic retains the reference float accumulation order and disables
fused multiply-add contraction; local arithmetic tests do not certify DST replay.

The builder evaluates the configured named models on every candidate, each
bound to its file digest: the PPG12 NPB preselection score and, for pp, the
frozen PPG12 classifier from `PhotonID/models/PP/nominal` (-1 below 7 GeV),
which nominal finalization copies into `pho_bdt_score` after a recompute check.
Models add scores only; no candidate is selected. `pho_presel_active` records
PPG12's newPPG12 preselection and `isow_*` the historical isolation working
points; the nominal classes are defined by PhotonID and the isolation freeze.

The producer applies no event gate. Every event keeps its GL1 words, the run
trigger table (`meta/trigger_run_info`), vertex, centrality and minimum-bias
decisions; `evt_terminal_status` is 0. Triggers are selected by exact name and
events are selected downstream.

Photon truth association retains the maximum-energy primary identity separately
from the producer's signal subset. Canonical response records cover all retained
objects; missing capture evidence remains unknown. Several clusters may match
one truth photon. Capture census counts and event-local indices are checked
before serialization, including empty populations and out-of-domain events.
The legacy in-producer diagnostic records (`pdiag_*`, `embdiag_*`,
`truthmiss_*`, PPG12 fill counters) and the PPG12 photon-yield tower-mask and
source-arbitration witnesses are not captured natively yet, so the mapping is
incomplete and the frontend refuses completion until they are resolved.

## Source map

| File | Responsibility |
|---|---|
| `produce_trees.sh` | Public commands and option forwarding |
| `scripts/produce.py` | Population preflight, execution, source assembly and shared finalization |
| `macros/Fun4All_PhotonJetTree.C` | Resolve, configure, register, run and finish Fun4All |
| `src/PhotonJetTree.h/.cc` | Producer lifecycle and visible event-processing order |
| `src/PhotonClusterBuilder.h/.cc` | Authoritative photon reconstruction |
| `src/internal/Production.h/.cc` | Shared job contract and validated plan |
| `src/internal/Setup.cc` | Conditions, reconstruction and synchronized inputs |
| `src/internal/Types.h` | Configuration, identities and record definitions |
| `src/internal/Event.cc` | Event context, exposure, triggers and weights |
| `src/internal/Objects.cc` | Photon, jet and pairing capture |
| `src/internal/Truth.cc` | Simulation truth and associations |
| `src/internal/Output.cc` | Persistence, provenance and completion |
