# FinalAnalysis

Calculate the per-photon xJγ measurement from TreeToHists packages, including
background corrections, unfolding, uncertainties, QA, figures and tables.
Requires ROOT 6 with `root` and `root-config` on PATH, a POSIX shell and
`sha256sum` or `shasum`, libyaml and `nlohmann/json.hpp`.
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

Declared analysis variations and systematics run offline from these packages.
A variation that changes tree selection or response filling needs a new
TreeToHists reduction. [pipeline.yaml](config/pipeline.yaml) retains the local
tree-to-results interface; the SDCC campaign finishes at histogram production.

## Run individual steps

```bash
./FinalAnalysis/final_analysis.sh run --config FinalAnalysis/config/analysis.yaml
./FinalAnalysis/final_analysis.sh qa --config FinalAnalysis/config/analysis.yaml
./FinalAnalysis/final_analysis.sh systematics --config FinalAnalysis/config/systematics.yaml
./FinalAnalysis/final_analysis.sh draw --results results.root --qa qa.root \
    --systematics systematics.root --fig figures
./FinalAnalysis/final_analysis.sh tables --results results.root \
    --systematics systematics.root --out tables
```

The command compiles the C++ locally on first use. Existing results and QA files are
not overwritten. Failed steps stop the full chain and leave their logs.
For live toy counts during analysis, set `FINALANALYSIS_PROGRESS_PATH` to an
absolute JSON file path.

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

## Output

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
| [Inputs.cc](src/Inputs.cc) | Checked inputs and exact-bin projections |
| [Purity.cc](src/Purity.cc) | ABCD purity and combinatoric subtraction |
| [Prior.cc](src/Prior.cc) | Data-fitted joint prior |
| [Unfolding.cc](src/Unfolding.cc) | Bayes unfolding, normalization and iteration rule |
| [Analysis.cc](src/Analysis.cc) | Ordered chain, toys, closures and saved results |
| [Systematics.cc](src/Systematics.cc) | Variation covariances, bands and ratios |
| [QA.cc](src/QA.cc) | Saved QA and study quantities |
| [Pipeline.cc](src/Pipeline.cc) | Full-chain orchestration |
| [RenderMetadata.cc](src/RenderMetadata.cc) | Saved plot semantics and provenance |
| [plots/](plots/README.md) | Figures, tables and shared style |

Full scientific acceptance requires adequate input statistics, qualified
normalization and weighting, and independent purity/combinatoric closure.
The estimator closure covers the prior fit and unfolding after corrections.
