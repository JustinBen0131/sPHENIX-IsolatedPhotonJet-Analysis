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
branches instead and exists for mechanism checks.

`sample_weights.py` is the one implementation of the complete per-event
analysis weight, `event_weight x source factor x centrality factor`, from
`config/samples.yaml`: the family cross-section scale for p+p simulation, the
cross section over generated events for embedded photon+jet samples, the
ownership stitching of embedded inclusive-jet slices by the leading R=0.4
truth-jet pT, and the 5 percent centrality maps for Au+Au simulation.  Pass
`--sample <name>` (or a second column in `--input-list`) to apply it; the
response builds use the same module.

`--truth-signal-only` restricts simulation to the best-matched candidate of
each truth photon satisfying the nominal contract, one per truth photon.  The
resulting ABCD counts are the prompt-photon leakage of the purity correction
and their Region-A row normalises the combinatoric template.
