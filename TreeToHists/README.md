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
completed collection. `MaxEventsPerJob` defaults to 50,000 (at most 1,000,000);
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
If samples were produced in separate campaigns, combine their validated collections
with `tree_to_hists.sh collection --from first/collection.yaml --from second/collection.yaml --out final`.
This copies final sample packages and their metadata; retained DATA run packages stay
in the original campaign folders, with their run lineage preserved in `final/collections/`.
The sidecars and ROOT metadata are read automatically. Partial or unknown
source coverage stays labeled; nominal p+p analysis requires both SI and DI
packages for each included family.

An explicit collection recipe may pair nominal DATA with one exact SIM response:
`collection --from DATA/collection.yaml --from SIM/collection.yaml --out combined
--simulation-response-variation 'RECIPE'`. Without a recipe, histogram settings must
agree exactly. The legacy `PhotonJetSimulationResponseCollectionV1` remains supported.

For a native EMS recipe containing the exact token
`photon_selection=ems_selection_response`, also pass
`--ems-selection-contract ems.json`. This is the same complete, source-bound JSON
contract supplied to the native fill. Only SIM may carry that exact contract and
the fixed `ems_qa_extensions=nominal_archival_only` and
`ems_candidate_archive=original_nominal_inputs_replay_unqualified` labels.
The declaration preserves the mechanical, unadopted cone status.

To refit a nuisance with varied DATA and immutable SIM, add
`--data-jet-scale-recipe data-scale.json`, for example
`{"system":"pp","nominal":1.0,"varied":1.02}`. Use `auau` for an AuAu DATA
scale. Both scales must be distinct, positive and finite. Supply DATA of the
declared system filled at the varied value; other-system DATA and every SIM
package retain the nominal value in that field. All other settings stay equal.
The output histogram settings come from the declared system's varied DATA, even
when p+p DATA supplies the collection's selection. The two options compose in one
`PhotonJetRecipeCollectionV2` declaration under `simulation_response`, and travel
unchanged through nested assembled collections. Omitting the SIM recipe permits
the DATA scale recipe with nominal SIM. Package files and metadata are copied
unchanged; these options never rescale histograms or grant scientific acceptance.

The metadata-only regression suite is `TreeToHists/tests/test_portable_collection.py`.
Set `TREETOHISTS_BINARY` to the built native assembler and
`FINALANALYSIS_INGEST_BINARY` to the compiled fixture
`FinalAnalysis/tests/collection_ingest.cc`. The latter includes the real
`Pipeline.cc` ingestion function; link with ROOT/libyaml and unused-function
elimination (`-Wl,-dead_strip` on macOS, or function sections and `--gc-sections`
on Linux). It validates packages without running analysis fits or event loops.

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

## Candidate records and systematic replay

The supplied nominal configuration performs a full fill with candidate records enabled.
It evaluates the selected AuAu H70 model natively and keeps the selected PPG12 saved score
for pp. The selected isolation definition and working-point lines are in
`config/treetohists.yaml`; DATA, both SIM families, spectra, QA and responses use that
same configuration. Set `Fill.JetRadius` to 0.4, 0.3 or 0.2 for each collection.
There is one full-fill path; no partial/add-on package or output-mode setting.

`Fill.PhotonJetRecords: 1`, `Fill.RecordMinEt: 0` and `Fill.RecordMaxRows: 10000000`
retain finite nonnegative-energy photons and the complete selected-radius jet and truth
collections of accepted events. A row ceiling fails explicitly; it never truncates a
package. Keep the zero floor for resolution variations so low-energy photons remain
available for migration tests.

Packages store `PhotonJetRecordsV1` under `records/candidates/`, separately from the
measurement tables in `records/`. Its tables are `events`, `photons`, `isolation`,
`jets`, `truth_photons`, `truth_jets`, `photon_links`, `jet_links` and `unsub_jets`.
Jet fields have radius-independent names such as `jet_pt` and `truthjet_pt`;
`jet_radius` in every event and in `meta/package.candidate_records` binds the view.
Events retain source/file, original entry, run, sequence, centrality, reconstructed
and truth vertex, MBD time and individual weight factors. Candidate `input_index`
retains its original per-event index. Sequence alone is insufficient to join files:
use source/part identity, original entry and original jet index for timing side tables.

SIM link indices and collection ranges are validated with the native matcher before
saving direct photon/jet truth matches. Its match witness is one for SIM and zero for
DATA. Redundant raw link and unsubtracted-jet tables are empty after SIM validation;
DATA retains its original representation. Both package merge levels preserve the same
record schema, radius, model-input order, view and semantic hash.

Photon records preserve raw producer width/NPB fields and the exact 14 selected H70
inputs separately. They are different definitions: raw `pho_npb_cluster_et2-4` and
named widths must not be substituted for the H70 model vector. Availability flags
preserve missing optional fields; NaN is not converted to zero. Raw AuAu isolation
is archival only and is excluded from downstream selection, tuning, response and
uncertainty calculations.

`CandidateMeasurementReplay.h` and `CandidateResponseReplay.h` replay records through
the same measurement and response accumulators used by a full fill. Replay requires
the original source identities, selected-radius records, selection and complete sample
normalization; it applies the saved event-weight factors once. Missing required fields,
wrong radii, discarded migrating photons and invalid match witnesses fail explicitly.
The uniform `leading_global` histograms describe the global leading jet followed by
fiducial/recoil vetoes. `scoreiso_fine` provides 100 score bins alongside the unchanged
20-bin `scoreiso` histograms and is included in the package inventory.

For local pp simulation, `pp.VertexFitDirectory` may point to a directory containing
`0mrad/reweight.root` and `1p5mrad/reweight.root`. This changes only the file locations;
the expected hashes and histogram names remain those of `TreeProduction`'s frozen
PPG12 configuration and are checked before use. Omit it on SDCC to use the original
paths. Store the same fit bytes, never a replacement fit.

The compiled `producer_code` records per-file hashes of the scientific source.
Compare these with the local source used for the build and retain the configuration,
input hashes and validated collection sidecars. Structural validation and software
fixture checks are separate from full-statistics physics acceptance.
# Optional native response-component grid

`Fill.ResponseGrid: /absolute/path/grid.json` opts the existing `plan` / campaign
planner, `fill` / grouped fill, and package assemblers into a bounded SIM response
grid. The JSON record schema is `PhotonJetNativeResponseGridSettingsV1`, with
`schema`, `variants`, `limits`, and an optional `ems_selection_contract`.
Each variant has exactly `name` and `variation`.
Include `{ "name": "nominal", "variation": "" }`; the other names are distinct
safe labels, and their resolved recipes must differ. Supported recipes are the
existing named photon resolution with `photon_selection=response_only` (optionally
the existing named jet resolution), or `centrality_weight=off` for AuAu SIM,
alone or with the nominal response-only resolutions. EMS components require
`photon_selection=ems_selection_response`, the nominal named resolutions, and
exactly one explicit source-bound `EmsNativeSelectionContractV1` under
`ems_selection_contract`. Other grids must omit that field. Leave
`Fill.ResponseVariation` and `Fill.EMSContract` empty when the grid is selected.

`config/response_grid_variations.json` records the R04 prescriptions: JER
0.100 +/- 0.018, PPG18 section 6.4 EMR low (no additional smear) and high
(0.13, 0.08, 0.08), and section 6.3 EMS factors 0.9852 and 1.0148.
Complete an existing six-component record without changing its limits or
existing recipe strings:

```bash
python3 TreeToHists/response_grid.py --base base_grid.json --system AuAu \
  --ems-contract ems.json --out complete_grid.json
```

The generator adds `smems_low`, `smems_high`, and, for AuAu, `centralityoff`;
it refuses an altered nominal/seed, duplicate names, or an existing output.
The historical empty `nominal` and corrected `smnom` retain their meanings.
Both EMS endpoints share seed 480 and the nominal source/entry/object draws.
EMS rescales original photon ET for the existing NPB, ID-domain, BDT and
isolation evaluator, and rescales the same nominal response realization.
Only the EMCal fraction of each jet receives the EMS factor. EMR uses original
ET for classification. The retained analysis-cone contract remains explicit;
raw AuAu isolation never replaces the analysis witness.

`fill-selection-group` retains the existing DATA view names for SIM selection
products and borrows their response/measurement sinks into the same response
reader. QA, extensions, input-header inspection and any stitching denominator
retain their own traversals. This is one response source pass per row, not one
total producer traversal. Named EMS components include responses, misses,
fakes, photon response/truth/reco, leakage and simulated Delta-phi transfer
ingredients. The grid does not claim complete varied QA/spectra packages.

Declare all five positive finite integer limits: `max_histogram_bytes`,
`max_resident_bytes`, `max_output_bytes`, `max_sparse_bins_per_histogram`, and
`max_seconds`. They remain inside the native envelope and are checked against
the actual axes, variant count, resident memory and observed sparse/output sizes.
`max_seconds` is the planned duration: progress reports `planned_wall_seconds`
and `over_plan`, and production continues when that estimate is exceeded.
Bounded foreground and canary invocations enforce their own process deadline.
The source reader retains the original logical source, entry range and object
indices, including consolidated sources whose logical ID differs from their
historical input hash. Every grid requires its actual compiled/runtime source
roster and a fully resolved selection, period/vertex fits, ownership, centrality
map and axes. There is no inferred resolution, calibration or farm admission.

The response callback visits each original entry once for the declared variants.
An additional file checksum pass is disclosed separately. Existing measurement,
QA, extensions, planner inspections and stitching-denominator reads retain their
own paths; this is **response-component fanout**, not one total producer read.
DATA uses its unchanged fill path. Native progress separates original entries
from event/variation work and publishes `FIRST_BAD` without a PASS receipt on
failure. A failed create-once grid remains evidence; choose a new work directory
when the incomplete grid cannot be resumed through its completed checkpoints.

The nominal grid product supplies the ordinary `response.root`. Named products
are bound in `done.json` and retained at
`response_components/<name>/response/...` in both raw and hierarchical sample
packages, with their own source-contribution cells and exact configuration in
`meta/package.response_components`. The package inventory and immutable binding
include these optional components; common nominal histogram settings remain
unchanged. These extra paths carry `RESPONSE_COMPONENT_ONLY` scope and do not
select a FinalAnalysis systematic recipe automatically.

Full EMS measurement/QA fanout remains refused. The separately source-bound EMS
single-variation path requires `Fill.EMSContract`; retained EMS replay, scientific
cone adoption, production inputs, other radii and site resources remain separate
qualification gates. Synthetic pp SI/DI and two-period fixtures qualify the
current production axes; their fit values are fixture data, never calibrations.



### Shared pp DATA JES fill

`fill-data-jes --work NOMINAL --up-work UP --down-work DOWN --job N` is an
opt-in fill of one row from three independently planned ordinary pp DATA
work directories. Prepare each plan with the same original catalogue, selection,
axes and ranges; only `Fill.DataJetScale.pp` differs: nominal `1/0.916`,
energy-up `1/0.907`, energy-down `1/0.923`. Set AuAu DATA scale to 1 and leave
the response recipe empty. Supply numeric factors in YAML. Existing products,
foreign/reused samples, duplicate identities and other settings differences
are refused before filling. The row command accepts ordinary plans. The grouped command below consumes
the existing native campaign group descriptors.

The same original spectra reader fills three complete measurement and timing
products before the threshold and xJ cuts. Candidate archives, when requested,
retain identical original inputs. QA and extensions use the nominal inputs and
are explicitly labelled `nominal_reused` in the two varied package metadata
records; they do not represent varied QA. Each directory then uses ordinary
`merge` and `collection`. Compose a varied pp DATA collection with unchanged
nominal SIM and other-system DATA through the existing
`--data-jet-scale-recipe` declaration. A partial failure requires exact row
reconciliation; this command does not infer rerun authority or resume an
incomplete three-view transaction.


`fill-data-jes-group --work NOMINAL --up-work UP --down-work DOWN --group G0
--scratch DIR` binds three already sealed native fill-group descriptors. Their
source order, entry ranges, settings and code must agree except for the exact
pp DATA scales above. It reuses the existing group scratch, row fill, native
merge and source-aware chunk publication lifecycle. Default `fill-group` and
its generated campaign recipes keep their ordinary behavior.

One invocation publishes three ordinary validated chunks, each bound to its
own plan/descriptor. Only nominal owns the physical progress total: one shared
science pass plus nominal QA and extensions, three event-passes per original
entry. Varied progress reports zero additional physical passes and links to
the shared owner; varied completion timing is null because it is shared.
The controller must bind the three descriptor nodes to this single invocation;
three unchanged campaign DAGs must not independently dispatch all three fills.
Merge nodes and typed collection composition retain their existing interfaces.

The command refuses existing result tuples, work/output aliases and mismatched
descriptor hashes, settings, selections, source identities/hashes or ranges
before opening an original. Source-file checking follows the existing native
metadata/stamp and sealed-publication contract; no extra original hash/read
pass is added. An interrupted publication leaves exact artifacts/scratch for
reconciliation and never authorizes an implicit retry. Successful completion
removes only its private, newly created group scratch after all three chunks
validate. Resource configuration and real-source validation remain required before
production.
