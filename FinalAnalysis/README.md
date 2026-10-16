# FinalAnalysis

Calculate the per-photon xJγ measurement from TreeToHists packages, including
background corrections, unfolding, uncertainties, QA, figures and tables.
Requires ROOT 6 with `root` and `root-config` on PATH, a POSIX shell and
`sha256sum` or `shasum`, libyaml, `nlohmann/json.hpp` and Python 3.
The supplied conda environment includes these dependencies.

## Analyze a transferred collection locally

Run TreeToHists on SDCC, download its finished collection with SFTP, then run
FinalAnalysis on your computer.

From the repository root:

```bash
./FinalAnalysis/final_analysis.sh all --config FinalAnalysis/config/pipeline_histograms.yaml
```

In [pipeline_histograms.yaml](config/pipeline_histograms.yaml), point
`Pipeline.HistogramCollection` to the downloaded `collection.yaml`.
The default output is `FinalAnalysis/outputs/analysis_v1/`, with plots in `fig/`.
Use a new output version for each iteration; [analysis.yaml](config/analysis.yaml)
holds the analysis settings.
The manifest's files, hashes, selection, sample identities and coverage are
checked before QA, analysis, figures and tables run. Paths are relative to
the configuration file.

`Pipeline.UseCollectionSamples: 1` selects the actual photon and inclusive-jet
sample names in that collection. Nominal p+p families require both SI and DI
packages. Partial coverage is preserved, with each package's exact exposure.
Analysis windows, corrections, toys and closure settings stay in `analysis.yaml`.

### Make a compact collection first

A delivered DATA package is almost all per-source bookkeeping (`records/source_contributions`, 20 to 30 GB) that the
measurement never reads; it reads the event records, a few MB. Make a compact twin once per delivery, on a fast
local disk, and point `Pipeline.HistogramCollection` at its `collection.yaml`:

```bash
./FinalAnalysis/final_analysis.sh compact --collection /delivered/collection.yaml --out /fast/disk/collection_compact
```

`compact` checks every delivered file against the manifest (the one full read of the delivery), then copies the
collection with each DATA package written without `records/source_contributions` and `records/sources`. The two
studies that need those records, the photon-trigger turn-on and the Au+Au centrality target, are computed once and
stored in the twin (`derived/`), by the same function the QA step uses on a full package.
`<package>.compact.json` beside each twin records the delivered hash it was made from. A collection whose file
matches neither its manifest nor such a record is refused, so a compact collection is as bound to the delivery as
the delivery itself. Simulated packages are copied unchanged.

### What a pass costs

Measured on a 12-core laptop with the four results (p+p and three Au+Au classes), a compact collection on the local
disk and the nominal toy counts (2000 data, 1000 response, 200 closure, and the iteration scan's own 2000 and 1000):

| Step | Time |
|---|---|
| `compact`, once per delivery | about 45 minutes |
| Building after a change to `src/` | 10 s for everything, a few seconds for one file |
| `qa` | 4 s |
| `run`, the nominal measurement | 12 s |
| `run`, one variation | 5 s |
| `draw`, 140 figures | 8 s, 15 s after a change to `plots/` |
| `all` with eight variations | 1 minute |

The first `run` on new packages or new bins takes about a minute longer, once: it reads every simulated source out of
the packages and saves what it needs (below). Each `run` prints every phase that took a second or more, so its log
shows where the time went.

Why it is fast, and what to keep true when changing the code:

- **Nothing is read twice.** The source bootstrap needs every simulated source's own contribution to the response.
  Reading those out of the packages is the slowest thing an analysis does, so `readSources` (Inputs.cc) keeps each
  source as one row of numbers under `build/cache/`, named by the packages, the bins and the code that made it, and
  reads the packages again only when one of those changes. Each row set must add up to the package's response in every
  number before it is used. A file's SHA-256 is kept the same way, under its path, size and time of writing.
- **A response toy is one pass over numbers.** A resampled response is the sum of the sources' rows with integer
  multiplicities (`combine`), not a sum of thousands of small matrices.
- **One process per result, one per variation.** `run` forks per entry of `Systems`; each result has its own random
  stream (`Uncertainty.Seed` and its position in `Systems`), so its numbers do not depend on what runs beside it.
  `all` runs `Pipeline.Parallel` variations at a time.
- **A variation computes what the systematics step uses.** That is its central value: its response toys and closure
  are switched off (set `Uncertainty.ResponseToys` or `Closure.Toys` in the variation to keep them), where the nominal
  takes its iterations from the iteration scan the variation uses the k the nominal selected, and a variation that changes only analysis
  settings (prior, unfolding, toys, combinatoric estimate, iterations, purity inputs) uses the nominal QA file. Its
  data toys stay, since the prior fit uses their covariance, so the central values are those of a full run.
  A variation may set `<variation>.HistogramCollection` to its own validated `collection.yaml`; its selection and
  input directory come from that collection, and it receives its own QA. The pipeline pins the manifest and payloads.
  An unvalidated directory or simultaneous tree reduction and collection override is rejected.
- **Only changed code is rebuilt.** One object per source file, compiled side by side; the figures and tables are
  ordinary programs built from the macros in `plots/`, without ROOT's dictionary step.

`build/` holds only what can be made again: delete it at any time. `FINALANALYSIS_CACHE=` (empty) saves nothing.
`./FinalAnalysis/final_analysis.sh build` compiles without running an analysis and prints the executable path.
The launcher uses [dependencies.py](dependencies.py) to ask the selected C++ compiler for every local included
header, including sibling-package headers; their bytes also identify the saved `readSources` rows.

### Iterating

To try a setting, copy [analysis.yaml](config/analysis.yaml), change it, give it a new `Output.Results` and run it;
then draw:

```bash
./FinalAnalysis/final_analysis.sh run --config my_analysis.yaml
./FinalAnalysis/final_analysis.sh draw --results my_results.root --fig my_figures
```

To change a figure of a finished run, edit the macro in `plots/` and draw again from the saved results; nothing is
recomputed:

```bash
./FinalAnalysis/final_analysis.sh draw --results outputs/analysis_v1/results/nominal.root \
    --qa outputs/analysis_v1/results/nominal_qa.root --fig my_figures
```

`all` pins the analysis code (`src/`) when it starts and the plot code (`plots/`) only when drawing starts, and fails
at the end if a pinned file changed meanwhile. So `plots/` can be edited while results are being computed, and `src/`
between runs.

Declared analysis variations and systematics run offline from these packages.
A variation that changes tree selection or response filling needs a new
TreeToHists reduction. [pipeline.yaml](config/pipeline.yaml) retains the local
tree-to-results interface; the SDCC campaign finishes at histogram production.

## Run individual steps

```bash
./FinalAnalysis/final_analysis.sh compact --collection collection.yaml --out collection_compact
./FinalAnalysis/final_analysis.sh run --config FinalAnalysis/config/analysis.yaml
./FinalAnalysis/final_analysis.sh qa --config FinalAnalysis/config/analysis.yaml
./FinalAnalysis/final_analysis.sh systematics --config FinalAnalysis/config/systematics.yaml
./FinalAnalysis/final_analysis.sh draw --results results.root --qa qa.root \
    --systematics systematics.root --fig figures
./FinalAnalysis/final_analysis.sh tables --results results.root \
    --systematics systematics.root --out tables
```

Existing results and QA files are not overwritten. Failed steps stop the full chain and leave their logs.

Plot comparisons use the same C++ drawing function as their baseline, with a
different declared population or projection. The renderer writes stable PNG
names and records the input hashes, source function and comparison recipes in
each sidecar. For example, the AuAu BDT score view draws truth-tagged background
and all inclusive-jet candidates through `savedScoreReview` in
[QAFigures.C](plots/QAFigures.C). Running `draw` regenerates both together;
changes to their shared style apply to both.

## Inputs and settings

Inputs are TreeToHists' DATA and simulation ROOT packages. Source records,
event records, selection definitions and normalization are checked on reading.
The portable manifest uses relative paths, so no SDCC tree paths are needed.
Generated JSON sidecars and ROOT metadata are read automatically.

For photon-ID training figures, set `Pipeline.PhotonIDModels` to the model
campaign's `models.yaml` (`QA.PhotonIDModels` for the individual QA step).
Copy the models directory with its relative layout intact. FinalAnalysis checks
the releases, imports their saved diagnostics and draws them without refitting
or recomputing training statistics. Routing and depth comparisons require
matching sources, event splits and training weights.

To draw saved training diagnostics before the histogram packages are available:

```bash
./FinalAnalysis/final_analysis.sh qa --models /path/to/models/models.yaml --out model_qa.root
./FinalAnalysis/final_analysis.sh draw --qa model_qa.root --fig figures
```

The AuAu collection supplies the training-weight, feature, routing and depth
figures. The pp collection supplies the pp/AuAu held-out comparison and the
original-versus-retrained PPG12 comparison. Bind the canonical AuAu reference
with PhotonID's `collectQA PP --auau-reference` before using the pp collection;
FinalAnalysis checks that reference and the shared learner settings. Each QA
output imports one model collection.

The BDT working-point and score-distribution figures use reduced histogram QA,
with the chosen release's `selection.json` supplied to TreeToHists. Saved
training diagnostics supply the training figures; they do not replace those
model-dependent analysis histograms.

| Configuration | Purpose |
|---|---|
| [analysis.yaml](config/analysis.yaml) | Packages, bins, corrections, prior, iterations, toys and outputs |
| [systematics.yaml](config/systematics.yaml) | Nominal and variation results, nuisance treatment and correlations |
| [pipeline.yaml](config/pipeline.yaml) | End-to-end inputs, commands, output directory and required figure list |
| [pipeline_histograms.yaml](config/pipeline_histograms.yaml) | Local analysis of a transferred, verified histogram collection |
| [pipeline_systematics.yaml](config/pipeline_systematics.yaml) | Variations used by the full chain |

These YAML files use nested settings and ordinary lists. Change settings for a study;
the same implementation performs the calculation. Unknown keys are errors.
Analysis bins must be exact unions of the package bins. The reported xJ lower
edge is derived from the jet threshold and photon window; its upper edge is configured.

`Systematics.Coverage: coverage.json` optionally requires an explicit uncertainty
roster. The JSON has schema `PhotonJetSystematicCoverageV1`, the exact
`nominal_sha256`, a unique `sources` list, and a `results` object keyed by every
nominal result. Each result declares every listed source with `status` and a
nonempty `reason`. Status is `included`, `not_yet_included`, or `not_applicable`.
Included entries require the configured variation products for that result;
excluded entries refuse supplied products. Missing declarations fail before
publication. The checked table is stored as `systematic_coverage` in the output,
and its file hash joins the input provenance. These labels describe coverage;
they do not establish scientific acceptance.

A source may explicitly name `Reference` and a nonempty `ReferenceReason` together.
Its signed variation shifts are then measured against that saved reference and
added to the configured nominal for bands, integrated quantities and ratios.
This is an additive delta transfer, which can be used for a declared first-pass
approximation when a source was evaluated with an earlier response. It does not
recompute that source with the new response. Every included result and its bins
must be present and compatible in the reference; missing inputs fail. The output
stores `systematic_source_references`, its treatment and reason, and hashes every
reference file. Omit both settings for ordinary variations of the nominal.

QA projection bins use `QA.AbcdEdges`. Drawing windows use `QA.Windows.ShowerShapes`, `IsolationOverlay`, and
`NcbScore` (two increasing bounds in GeV), and `TruthIsolation` (consecutive low/high pairs, for example
`[15, 20, 25, 30]`). Omitted values retain the defaults in `analysis.yaml`. The QA output records the resolved values
in each package's render selection, under `qa.abcd_edges_gev` and `qa.windows`, alongside the input selection.

## Output

### Selecting a saved named response component

For a saved typed histogram collection containing
`response_components/<label>/response/...`, a result may opt in with
`ResponseComponent`, `ResponseRecipe`, and `ResponseGridSHA256` together.
Copy the exact label/recipe and the SHA256 of the canonical JSON configuration
from `meta/package.response_components.configuration`; these fields select
persisted counts and never apply smearing again. Omit all three for nominal
behavior. Direct analysis also requires `Input.Collection: collection.yaml`;
the existing histogram Pipeline binds it automatically. A variation sets the
fields under `<variation>.Analysis.<result>` and can reuse nominal QA.

Typed collection/package hashes, grid/recipe/source compatibility and response
inventory/axes are checked before consumption. The selected response, leakage,
dphi transfer and source rows use the same namespace and retain once-only
sample normalization. DATA measurements and nominal spectra/QA keep their
existing paths. Selection enters response-cache identity and saved provenance.
This consumer currently requires `Combinatoric.Method: data_dphi` and
`Combinatoric.Shape: simulation`. Pure centrality-off is an exclusive AuAu
alternative. EventTree reductions, replay/other response overlays and full EMS
selection fanout are refused. A component's availability does not establish
scientific acceptance, closure or a systematic prescription.

- `results.root`: corrected and unfolded measurements, toys, scans and closures.
- `qa.root`: photon/event QA and saved study quantities.
- `systematics.root`: systematic covariances and bands.
- Figures and LaTeX/CSV tables in the requested directories.
- `PIPELINE.json` records a full run; `FIGURES.txt` lists written and unavailable figures.

Rendering uses the shared sPHENIX style and saved results.
[fig/README.md](fig/README.md) maps output names to IAN figures.
Outputs are local files; the commands do not edit a manuscript.

## Source map

Start at `Analysis.cc::main`; `runSystem` shows one measurement's step order.
The shared interfaces, array layouts and units are in `FinalAnalysis.h`.

| File | Responsibility |
|---|---|
| [FinalAnalysis.h](src/FinalAnalysis.h) | Settings, data types and step interfaces |
| [Inputs.cc](src/Inputs.cc) | Checked inputs, exact-bin projections and the saved simulated source rows |
| [Purity.cc](src/Purity.cc) | ABCD purity and combinatoric subtraction |
| [Prior.cc](src/Prior.cc) | Data-fitted joint prior |
| [Unfolding.cc](src/Unfolding.cc) | Bayes unfolding, normalization and the iteration scan |
| [Analysis.cc](src/Analysis.cc) | Ordered chain, toys, closures and saved results |
| [Systematics.cc](src/Systematics.cc) | Variation covariances, bands and ratios |
| [QA.cc](src/QA.cc) | Saved QA and study quantities; `sourceStudies` holds the two that read per-source records |
| [Pipeline.cc](src/Pipeline.cc) | Full-chain orchestration, collection checks and the compact twin |
| [RenderMetadata.cc](src/RenderMetadata.cc) | Saved plot semantics and provenance |
| [Phases.h](src/Phases.h) | The time each phase of an analysis took |
| [plots/](plots/README.md) | Figures, tables and shared style |

Full scientific acceptance requires adequate input statistics, qualified
normalization and weighting, and independent purity/combinatoric closure.
The estimator closure covers the prior fit and unfolding after corrections.
