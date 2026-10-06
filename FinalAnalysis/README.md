# FinalAnalysis

Calculate the per-photon xJγ measurement from TreeToHists packages, including
background corrections, unfolding, uncertainties, QA, figures and tables.
Requires ROOT 6 with `root` and `root-config` on PATH, a POSIX shell and
`sha256sum` or `shasum`, libyaml and `nlohmann/json.hpp`.
The supplied conda environment includes these dependencies.

## Run the full chain

From the repository root:

```bash
./FinalAnalysis/final_analysis.sh all --config FinalAnalysis/config/pipeline.yaml
```

Copy the supplied configuration files into a run directory and set your
inputs, outputs and requested figures. Paths are relative to the configuration.
`all` runs reduction, QA, nominal analysis, declared variations, systematics,
figures and tables. To start from existing histograms, set
`Pipeline.HistogramCollection` to the histogram-package directory instead of
setting `Pipeline.InputCollection` and reduction settings.

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

The command compiles the C++ on first use. Existing results and QA files are
not overwritten. Failed steps stop the full chain and leave their logs.

## Inputs and settings

Inputs are TreeToHists' DATA and simulation ROOT packages. Source records,
event records, selection definitions and normalization are checked on reading.
For the full chain, `Pipeline.InputCollection` points to TreeToHists' `inputs.yaml`.

| Configuration | Purpose |
|---|---|
| [analysis.yaml](config/analysis.yaml) | Packages, bins, corrections, prior, iterations, toys and outputs |
| [systematics.yaml](config/systematics.yaml) | Nominal and variation results, nuisance treatment and correlations |
| [pipeline.yaml](config/pipeline.yaml) | End-to-end inputs, commands, output directory and required figure list |
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
