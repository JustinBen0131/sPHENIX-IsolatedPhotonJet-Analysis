# TreeToHists

Reduce published eventTree files to one ROOT histogram package per sample.
Packages preserve spectra, responses, QA, Sumw2, event records and original
source contributions for FinalAnalysis.

## Build

Load the ROOT environment on your computer, or the selected sPHENIX release
on SDCC. The build requires a C++20 compiler, ROOT 6, Autoconf, Automake,
libyaml and `nlohmann/json.hpp`.

```bash
./TreeToHists/tree_to_hists.sh build
```

The conventional build files are in [src/](src/): `autogen.sh`, `configure.ac`
and `Makefile.am`. Build once before launching farm jobs.

## Inputs and selection

Copy [inputs.yaml](config/inputs.yaml) and [treetohists.yaml](config/treetohists.yaml).
Name the published tree groups or ROOT files in the input YAML. Paths are
relative to that file. DST `.list` files belong to TreeProduction.

```yaml
samples:
  - system: PP
    name: data
    family: data
    groups: [/path/to/published_pp_group]
```

The selection defaults to the new AuAu 15–35 GeV nominal model, the original
PPG12 model for pp, the frozen isolation definitions, pp Photon4, and AuAu
Photon12 with minimum-bias selection. Simulation families
are `photonjet` and `inclusivejet`; keep each named SI or DI population separate.

DATA produces `<system>/run_<run>.root` for every input run, followed by
`<system>/data.root` combined from the validated run files. Keep the run files
for offline QA. `exclude_runs` omits a run only from the combined analysis;
its run file is retained. Optional `run_population_sources` and `run_luminosity`
maps record expected source counts and qualified luminosity by run.
Missing values stay unknown. Trigger-study spectra use each run's GL1
live/scaled correction before combination; nominal physics weights are unchanged.

To add newly published inputs, set `reuse_collection` in the input YAML to a
previous `collection.yaml` and use new work/output directories. Unchanged run
and simulation packages are reused only when their inputs, code, selection and
corrections match. Changed runs are rebuilt; combined DATA is rebuilt from the
validated run files.

Choose the model and its working points in `treetohists.yaml`:

```yaml
Selection:
  PhotonIDModels: /path/to/models  # defaults to PhotonID/models
pp:
  Model: ppg12_original
auau:
  Model: nominal_retrained_v1     # wideET_v1 or any other released AuAu model
  TightWP: WP70
  NonTightWP: WP80
  NonTightMinimum: 0.1
```

Normal reduction reads the selected stored score. Add missing scores with
PhotonID's `scoreCampaign` first; an unavailable score or mismatched release
receipt fails clearly. ID, candidate ranking and regions are derived from
that score and its own cuts, independently of historical tight/non-tight labels.
The pp original model retains its established cuts; released models provide
WP70, WP80 and WP90. Overlapping working-point bands are rejected.

`ScoreMode: native` explicitly evaluates a released BDT from saved H70 features
during reduction. It uses the same evaluator, routing and cuts as stored
scoring. Existing `PhotonIDSelection: /path/to/selection.json` recipes retain
this native route; choose a named `Model` or an explicit path, once per system.

p+p simulation weights sum the 0 and 1.5 mrad periods using their luminosity
shares, interaction fractions and truth-vertex fits from TreeProduction's
`ppg12_weights`. Absolute normalization uses that family's reference cross
section divided by all entries of the selected trees, including inactive
entries. FinalAnalysis applies it once. AuAu simulation needs the qualified
centrality weights and truth-vertex fit specified in the input template.

## Produce on SDCC

Copy [campaign.yaml](config/campaign.yaml), set your GPFS directories, setup
script and measured timing file, then prepare and inspect the frozen campaign:

```bash
./TreeToHists/tree_to_hists.sh campaign --config campaign.yaml --prepare-only
./TreeToHists/tree_to_hists.sh campaign --config campaign.yaml --submit
```

Launch after site preflight and production approval. The DAG checks input
headers in parallel (`InspectSourcesPerJob`, `InspectConcurrency`), assembles
the plan on a compute node, fills bounded source groups, merges disjoint sources and validates the
completed collection. `MaxEventsPerJob` defaults to 50,000 (at most 100,000);
`MaxSourcesPerJob` limits the small inputs read serially by each worker. Samples
and DATA runs stay separate. Per-source phase files use private worker scratch;
GPFS receives one validated ROOT chunk per group, in sharded directories.
`MergeFanIn` defaults to 200: independent jobs merge completed chunks into
partials, then merge those partials through further levels. Merges start as
their inputs finish. They open one ROOT file at a time and stream source/event
records; each new package passes the full check before publication.
`FillConcurrency` and `MergeConcurrency` limit each stage. Set the separate
`FillMemoryMB`, `MergeMemoryMB` and `CollectionMemoryMB` requests from your
measurements. The final collection job checks small validation receipts and
file hashes rather than loading the complete input inventory again.
A source larger than the event cap uses deterministic ranges with its original
identity retained. Timing includes all fill
passes and must bind the current code, settings and ROOT build; ranges target eight
hours with a 1.5 safety factor. Workers use the prebuilt executable.

Completed fill phases and packages are reused when their hashes and bindings
match. Input reads have three bounded attempts. Other failures stop the DAG
and leave checkpoint and rescue evidence for an explicit restart. Changed
inputs or settings require a new work directory.

## Transfer and analyze

The finished GPFS folder contains `collection.yaml`, `selection.yaml`,
`histograms.yaml`, and `<system>/<sample>.root` with its JSON sidecar. Transfer
that folder through SDCC SFTP, then run [FinalAnalysis](../FinalAnalysis/README.md)
locally with [pipeline_histograms.yaml](../FinalAnalysis/config/pipeline_histograms.yaml).
The sidecars and ROOT metadata are read automatically. Partial or unknown
source coverage stays labeled; nominal p+p analysis requires both SI and DI
packages for each included family.

For a local reduction, use the same implementation:

```bash
./TreeToHists/tree_to_hists.sh local --settings treetohists.yaml \
    --inputs inputs.yaml --work work --out histograms --parallel 4
```

`centrality` derives AuAu weights from a DATA package. The `isolation`
commands derive and check an isolation definition. Run the entry point
without arguments for their syntax.

For the AuAu truth-vertex fit, copy [auau_vertex.yaml](config/auau_vertex.yaml)
and set the pinned inputs and centrality-weight file. `vertex population`
retains the needed tree sums; `vertex derive` writes `vertex.root` and
`DERIVATION.json`. The fit is a candidate until event-level vertex and
centrality closure checks qualify the production weight receipt.

```bash
./TreeToHists/tree_to_hists.sh vertex population --config auau_vertex.yaml --out population.root
./TreeToHists/tree_to_hists.sh vertex derive --config auau_vertex.yaml \
    --population population.root --out vertex_fit
```

## Source map

| File | Responsibility |
|---|---|
| [TreeToHists.cc](src/TreeToHists.cc) | Selection, input planning, filling and raw merge |
| [Campaign.cc](src/Campaign.cc) | Frozen campaign, Condor dependencies and collection publication |
| [Package.cc](src/Package.cc) | Source-aware hierarchical merge and package checks |
| [Corrections.cc](src/Corrections.cc) | AuAu vertex population and PPG12 fit candidate |
| [PhotonJet.h](src/PhotonJet.h) | Shared cuts, weights, truth relations and bounded readers |
| [FillPhotonJet.C](src/FillPhotonJet.C) | Measurement histograms and event records |
| [FillPhotonJetResponse.C](src/FillPhotonJetResponse.C) | Response, misses and fakes |
| [FillPhotonQA.C](src/FillPhotonQA.C) | Established QA and retained-observable extensions |
| [Isolation.cc](src/Isolation.cc) | Isolation derivation and checks |

Package validation checks structure, hashes and source/event accounting.
Scientific acceptance also requires adequate statistics and analysis closure.
