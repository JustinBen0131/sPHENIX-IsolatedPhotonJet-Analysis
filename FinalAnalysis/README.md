# FinalAnalysis

Calculate the per-photon xJγ measurement from TreeToHists packages, including
background corrections, unfolding, uncertainties, QA, figures and tables.
Requires ROOT 6 with `root` and `root-config` on PATH, a POSIX shell and
`sha256sum` or `shasum`. The JSON headers are bundled.

## Run the full chain

From the repository root:

```bash
./FinalAnalysis/final_analysis.sh all --config FinalAnalysis/config/pipeline.cfg
```

Copy the supplied configuration files into a run directory and set your
inputs, outputs and requested figures. Paths are relative to the configuration.
`all` runs reduction, QA, nominal analysis, declared variations, systematics,
figures and tables. To start from existing histograms, set
`Pipeline.HistogramCollection` instead of a tree catalogue and reduction settings.

## Run individual steps

```bash
./FinalAnalysis/final_analysis.sh run --config FinalAnalysis/config/analysis.cfg
./FinalAnalysis/final_analysis.sh qa --config FinalAnalysis/config/analysis.cfg
./FinalAnalysis/final_analysis.sh systematics --config FinalAnalysis/config/systematics.cfg
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

| Configuration | Purpose |
|---|---|
| [analysis.cfg](config/analysis.cfg) | Packages, bins, corrections, prior, iterations, toys and outputs |
| [systematics.cfg](config/systematics.cfg) | Nominal and variation results, nuisance treatment and correlations |
| [pipeline.cfg](config/pipeline.cfg) | End-to-end inputs, commands, output directory and required figure list |
| [pipeline_systematics.cfg](config/pipeline_systematics.cfg) | Variations used by the full chain |

These are ROOT TEnv-style `Key: value` files. Change settings for a study;
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
