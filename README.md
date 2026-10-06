# sPHENIX photon–jet analysis

Four packages for isolated-photon and jet analysis in pp and AuAu DATA and simulation:

**TreeProduction → TreeToHists → FinalAnalysis**

TreeProduction includes the nominal photon-ID scores and isolation labels.
Use PhotonID when you want to train another model or add its scores.

| Package | Purpose | Guide |
|---|---|---|
| TreeProduction | Reconstruct DSTs and write eventTrees | [Build and run](TreeProduction/README.md) |
| PhotonID | Train models and add named photon-ID scores | [Train and score](PhotonID/README.md) |
| TreeToHists | Apply selections and fill histogram packages | [Make histograms](TreeToHists/README.md) |
| FinalAnalysis | Calculate measurements, uncertainties, figures and tables | [Analyze results](FinalAnalysis/README.md) |

Start with TreeToHists if you already have eventTrees, or FinalAnalysis if you
already have histogram packages. Keep the four-package directory layout intact.

## Setup

```bash
conda env create -f environment.yml
conda activate photonjet
python TreeProduction/scripts/check_environment.py
```

DST reconstruction also requires the sPHENIX release and the
[TreeProduction libraries](TreeProduction/README.md#build).
PhotonID's additional training dependencies are listed in its guide.

## Run

Examples below run from the repository root. Set the input and output paths
in the package configuration files before running.

```bash
./TreeProduction/produce_trees.sh --help
./PhotonID/photon_id.sh --help
./TreeToHists/tree_to_hists.sh local --settings TreeToHists/config/treetohists.cfg \
    --inputs catalogue.json --work work --out histograms
./FinalAnalysis/final_analysis.sh all --config FinalAnalysis/config/pipeline.cfg
```

Each package has its own configuration and a short source map in its README.
Calculations live in C++/ROOT; PhotonID uses Python for model training and scoring.
Figures read saved analysis results.

## Products

- **Trees:** one ROOT `eventTree` per logical DST source, with source metadata.
- **Optional scores:** additional photon-aligned branches in a new published tree version.
- **Histograms:** one ROOT package and JSON sidecar per DATA collection or SIM sample.
- **Analysis:** results, QA and systematic ROOT files, then figures and tables.

TreeProduction's `.list` files select DST inputs. TreeToHists' catalogue selects
the resulting ROOT files or published groups. The package guides explain the
required settings and current scientific qualification limits.
