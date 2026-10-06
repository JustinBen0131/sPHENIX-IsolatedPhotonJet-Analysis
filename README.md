# sPHENIX photon–jet analysis

Produce and analyze isolated photons and jets in pp and AuAu collisions,
using DATA and simulation through four packages.

**TreeProduction → TreeToHists → FinalAnalysis**

Trees include the nominal photon-identification scores and isolation labels.
Use **PhotonID** to train additional models or add their scores before making
histograms.

## Packages

| Package | Purpose | User guide |
|---|---|---|
| **TreeProduction** | Reconstruct DST inputs and write final eventTrees. | [Build and produce trees](TreeProduction/README.md) |
| **PhotonID** | Train photon-ID models and add named scores to trees. | [Train and score models](PhotonID/README.md) |
| **TreeToHists** | Apply selections and fill spectra, QA and simulation response histograms. | [Make histogram collections](TreeToHists/README.md) |
| **FinalAnalysis** | Calculate measurements and uncertainties, then produce figures and tables. | [Analyze and render results](FinalAnalysis/README.md) |

Start with TreeToHists if you already have eventTrees, or FinalAnalysis if you
already have a histogram collection. Each package's README explains its
commands, inputs, outputs and configuration.

## Setup

Create the analysis environment from the repository root:

```bash
conda env create -f environment.yml
conda activate photonjet
python TreeProduction/scripts/check_environment.py
```

DST reconstruction uses the sPHENIX software environment. Follow the
[TreeProduction build instructions](TreeProduction/README.md) to compile its
C++ libraries using the Autotools files in `TreeProduction/src/`; a ROOT-only
build (`--disable-producer`) provides the eventTree library that PhotonID
uses to read and extend saved trees. Model-training dependencies are
described in the [PhotonID guide](PhotonID/README.md).

## Commands

Run these from the repository root to see the available datasets, commands
and options:

```bash
./TreeProduction/produce_trees.sh --help
./PhotonID/photon_id.sh --help
./TreeToHists/tree_to_hists.sh
./FinalAnalysis/final_analysis.sh
```

For example, fill the histogram packages of a tree catalogue, measure, and draw the IAN's figures and tables:

```bash
./TreeToHists/tree_to_hists.sh local --settings TreeToHists/config/treetohists.cfg --inputs catalogue.json --work work --out histograms
./FinalAnalysis/final_analysis.sh run --config FinalAnalysis/config/analysis.cfg
./FinalAnalysis/final_analysis.sh qa  --config FinalAnalysis/config/analysis.cfg
./FinalAnalysis/final_analysis.sh draw --results results.root --qa qa.root --fig FinalAnalysis/fig
./FinalAnalysis/final_analysis.sh tables --results results.root --out FinalAnalysis/tables
```

## Data products

- **Trees:** one ROOT file per logical DST source, containing one `eventTree`
  with event, photon, jet and simulation-truth information, nominal scores and
  classes, and calorimeter fractions for R = 0.2, 0.3 and 0.4. ROOT metadata
  records the source and reconstruction settings.
- **Additional scores:** named photon-aligned branches added by PhotonID,
  preserving the nominal definitions and dataset layout.
- **Histograms:** one ROOT package per DATA collection or SIM sample, with a
  JSON sidecar and event/source records for uncertainty calculations.
- **Results:** `results.root` (the measurement), `qa.root` (photon and event
  QA) and `systematics.root`; figures and tables are drawn from them.

Configuration lives inside its owning package. Standard production commands
select maintained input collections and reconstruction presets; custom inputs,
analysis selections and output locations use the options in each user guide.
