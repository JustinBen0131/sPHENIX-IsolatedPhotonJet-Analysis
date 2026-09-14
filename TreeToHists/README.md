# TreeToHists: nominal histograms

`make_histograms.py` reads scored collaborator trees and writes, per photon-pT
bin (15-20, 20-25, 25-35 GeV):

* event-leading photon counts in the ABCD regions, with Sumw2, and
* the recoil xJgamma spectrum in regions A and C, with Sumw2,

as a JSON payload (the input format of `FinalAnalysis/run_corrections.py`) and
as ROOT histograms.  All cuts come from `config/nominal.yaml`; the thresholds
come from the configuration by default (`--thresholds-from config`) so that
editing a threshold and rerunning this script is enough to regenerate the
histograms.  `--thresholds-from tree` uses the per-candidate threshold
branches instead and exists for mechanism checks.  `--truth-signal-only`
restricts simulation to candidates linked to a nominal truth-signal photon,
which gives the prompt-photon leakage inputs of the purity correction.
