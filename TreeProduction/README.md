# TreeProduction

Reconstruct sPHENIX DSTs into analysis trees for pp and AuAu DATA and simulation.
One Fun4All producer writes one schema-bound `eventTree` per logical DST source.
Read [EVENT_TREE.md](EVENT_TREE.md) for the tree layout and indexing;
[SCHEMA.md](SCHEMA.md) lists every branch, type and unit.

## Build

Use the configured sPHENIX release (`ana.568`), ROOT, FastJet, HepMC, CLHEP,
Autotools and Python 3 with PyYAML. Set `OFFLINE_MAIN` to the release and
`MYINSTALL` to your installation prefix. From the repository root:

```bash
mkdir -p TreeProduction/build
cd TreeProduction/build
../src/autogen.sh --prefix="$MYINSTALL"
make -j4
make install
cd ../..
```

This installs `libPhotonJetTree` (the producer) and `libEventTree` (saved-file
operations). For a ROOT-only reader build, add `--disable-producer` to the
configure command. Running the producer also needs Python's ROOT bindings.

## Run

```bash
./TreeProduction/produce_trees.sh --help
./TreeProduction/produce_trees.sh runPP --run 51109
./TreeProduction/produce_trees.sh runAuAu --run-list my_auau_runs.list
./TreeProduction/produce_trees.sh runPhotonJetSim --sample run28_photonjet20
./TreeProduction/produce_trees.sh runAll --dry-run
```

| Command | Input |
|---|---|
| `runPP`, `runAuAu` | pp or AuAu DATA |
| `runPhotonJetSim`, `runInclusiveJetSim` | pp single-interaction simulation |
| `runPhotonJetDISim`, `runInclusiveJetDISim` | Archived pp double-interaction simulation |
| `runEmbeddedPhotonJetSim`, `runEmbeddedInclusiveJetSim` | AuAu embedded simulation |
| `runAllData`, `runAllSim`, `runAll` | All maintained collections in the chosen category |

The command chooses the input lists and reconstruction profile automatically.
Use `--run-list` or `--sample-list` for plain files with one identifier per line;
blank lines and `#` comments are allowed. `--output` selects the output root,
`--resume` reuses verified completed work, and `--dry-run` checks the plan
without reconstruction. `--source` and `--events` make an inspection run under
`development/`. Existing products are not overwritten.

## Inputs and settings

The maintained input lists cover the full pp/AuAu good-run lists and simulation samples:

```text
lists/pp_data/<run>/         paired DST .list files and sources.tsv
lists/auau_data/<run>/       paired DST .list files and sources.tsv
lists/simulation/<sample>/  DST stream .list files and sources.tsv
```

Each DST list contains one filename per line; corresponding rows across
streams belong to the same source. Keep paired lists aligned. The source
sidecars are read automatically.

| Configuration | Purpose |
|---|---|
| [tree_production.yaml](config/tree_production.yaml) | Reconstruction profiles, conditions and pinned definitions |
| [production_populations.yaml](config/production_populations.yaml) | Maintained populations and default output location |
| [eventTree_v1.yaml](config/eventTree_v1.yaml) | Output branch contract |

`lists/pp_data/runs.list`, `lists/auau_data/runs.list` and
`lists/simulation/samples.list` contain the available identifiers.
Custom catalogues are supported with `--catalogue`.
The nominal photon-ID and AuAu isolation definitions come from the sibling
PhotonID and TreeToHists packages.

## Output and reading

```text
<output>/<campaign>/<lane>/<run-or-sample>/part-<ordinal>.root
```

Trees retain event context, photons, jets, simulation truth and nominal
selection labels. Empty events remain in the exposure accounting.
Completed run/sample groups have immutable versions selected by `CURRENT`;
downstream tools can read the group directory.

The [Python](examples/quick_start.py) and [C++](examples/quick_start.C) examples
show the event predicate and photon/jet indexing. They make unweighted
illustrative distributions; TreeToHists supplies the measurement selections
and weights.

## Source map

Follow the command driver, steering macro and `PhotonJetTree::process_event`.
The latter shows the capture order; the owning helpers are below.

| File | Responsibility |
|---|---|
| [produce.py](scripts/produce.py) | Input selection, source scheduling and saved-file checks |
| [Fun4All_PhotonJetTree.C](macros/Fun4All_PhotonJetTree.C) | Native job setup and execution |
| [PhotonJetTree.cc](src/PhotonJetTree.cc) | Producer lifecycle and event order |
| [PhotonClusterBuilder.cc](src/PhotonClusterBuilder.cc) | Photon reconstruction, shapes and isolation |
| [EventTree.h](src/EventTree.h), [EventTree.cc](src/EventTree.cc) | ROOT-only writing, checking, assembly and publication |
| [Production.h](src/internal/Production.h), [Production.cc](src/internal/Production.cc) | Job/plan definitions and configuration resolution |
| [Configuration.h](src/internal/Configuration.h), [Types.h](src/internal/Types.h) | Configuration access, units, identities and records |
| [Setup.cc](src/internal/Setup.cc) | Conditions, reconstruction and synchronized inputs |
| [Event.cc](src/internal/Event.cc) | Event context, classification, weights and bookkeeping |
| [Objects.cc](src/internal/Objects.cc) | Photons, jets, pairs and nominal labels |
| [Truth.cc](src/internal/Truth.cc) | Truth and reconstructed-to-truth associations |
| [Output.cc](src/internal/Output.cc) | Record-to-branch mapping and completion metadata |

Input coverage and a successful build do not establish production acceptance.
Calorimeter-fraction equivalence remains unqualified; reconstruction piece
boundaries can affect FastJet areas and final kinematic digits.
