# sPHENIX isolated photon + jet analysis

Code to go from sPHENIX DSTs to the unfolded photon-jet momentum balance
`xJgamma = pT(jet) / ET(photon)` in p+p and Au+Au collisions, organised as four
stages that hand off through documented files:

```text
DST                      TreeProduction/   Fun4All producer + tree builder
  -> collaborator trees  PhotonID/         train, score, derive working points
  -> scored trees        TreeToHists/      nominal selection, ABCD and xJ histograms
  -> histograms          FinalAnalysis/    purity correction, unfolding, overlay
  -> final xJgamma comparison plot
```

One configuration file, `config/nominal.yaml`, holds the physics choices shared
by every stage: the truth-photon definition, the R = 0.3 isolation cone, the
photon and jet kinematic windows, the per-system model inputs and the
working-point thresholds.  Changing a threshold regenerates flags and
histograms without retraining a model or re-reading a DST.

Status: **first-pass candidate.**  Every stage has a real command and the
Python stages are exercised end to end on the small fixture trees in
`tests/fixtures`.  Models and working points for the nominal selection are
not yet derived (see *What is still needed*), so the configuration ships with
`null` thresholds and the stages stop with an explicit message until they are
filled in.

## Physics definitions

* **Signal photon (truth).**  A Geant photon with a valid generator
  association, PPG12 photon class in {-1, 0, 1, 2}, truth isolation energy in
  a cone of R = 0.3 below 4 GeV and |eta| < 0.7.  Missing associations are
  never admitted as signal.  This one definition is the BDT training target and
  the selected truth photon in simulation.
* **Reconstructed isolation.**  Cone R = 0.3 (`iso_r03` in the trees), with the
  established estimator and underlying-event subtraction from the producer.
  R = 0.4 values are kept in the trees as a control.
* **Photon identification.**  One gradient-boosted decision tree per collision
  system on the ordered shower-shape inputs stored per candidate: 11 inputs in
  p+p, 14 in Au+Au (the p+p list plus the 3x3 widths and centrality).
* **Working points.**  Tight: `score > tight_threshold`.  Non-tight (bounded):
  `nontight_low < score <= nontight_high`, which may leave a gap below the tight
  threshold.  Isolated: `iso_r03 < isolated_max_gev`.  Non-isolated:
  `iso_r03 > nonisolated_min_gev`.  Candidates between two thresholds belong to
  neither class; candidates with a missing score or threshold are invalid, not
  background.
* **ABCD regions.**  A = tight & isolated (signal region), B = tight &
  non-isolated, C = non-tight & isolated, D = non-tight & non-isolated.  The
  event-leading photon is chosen independently in every region among photons
  with 15 <= ET < 35 GeV and |eta| < 0.7.
* **Recoil jets.**  Anti-kT R = 0.4, pT > 5 GeV, |eta| < 0.7, back to back
  (delta phi > 7pi/8).  xJgamma is histogrammed in three photon-pT bins
  (15-20, 20-25, 25-35 GeV).
* **Weights.**  The stored `event_weight` is applied exactly once per event
  (photon counts) or pair (recoil spectra); Sumw2 is carried everywhere.

## Requirements

* sPHENIX software environment with ROOT >= 6.26 and TMVA (PyROOT is used for
  scoring, working points and ROOT histogram output).
* Python 3.10+ with the packages in `requirements.txt`
  (`pip install -r requirements.txt`).  `xgboost` and `scikit-learn` are only
  needed for training.
* Building the producer needs the sPHENIX `coresoftware` install
  (`Fun4All`, `calobase`, `caloreco`, `jetbase`, `calotrigger`, `centrality`).

Every command below is run from the repository root.

## 1. TreeProduction: DST to collaborator trees

Build the two analysis modules once in an sPHENIX shell:

```bash
cd TreeProduction/src && ./autogen.sh --prefix=$MYINSTALL && make -j4 install && cd ../..
cd TreeProduction/src_AuAu && ./autogen.sh --prefix=$MYINSTALL && make -j4 install && cd ../..
```

Run the producer on one DST list (the macro reads `TreeProduction/macros/analysis_config.yaml`):

```bash
root -l -b -q 'TreeProduction/macros/Fun4All_recoilJets.C(0, "pp_dsts.list", "pp_producer.root")'
root -l -b -q 'TreeProduction/macros/Fun4All_recoilJets_AuAu.C(0, "auau_dsts.list", "auau_producer.root")'
```

Convert producer output into the collaborator trees (no model required; the
score branches are left unscored and are filled by PhotonID):

```bash
python TreeProduction/build_photonjet_collaboration_tree.py --system pp \
    --input pp_producer.root --output trees/pp_data.root
python TreeProduction/validate_photonjet_collaboration_tree.py --system pp trees/pp_data.root
```

The output contains eight trees (`events`, `eventTree`, `photons`, `jets`,
`photonJets`, `truthPhotons`, `truthJets`, `recoTruthLinks`); the branch
contract is in `contracts/photonjet_trees_v1_branches.json`.  See
`TreeProduction/README.md` for the producer options and the site paths that
must be edited before building elsewhere.

## 2. PhotonID: train, score, derive working points

```bash
# Train one model per system from photon+jet (signal) and inclusive-jet
# (background) simulation trees.  Needs xgboost and scikit-learn.
python PhotonID/train_photon_bdt.py --system pp --config config/nominal.yaml \
    --signal trees/pp_photonjet_sim.root --background trees/pp_inclusive_sim.root \
    --output-dir models/pp

# Score data and simulation with the same per-system model; flags are
# regenerated from config/nominal.yaml and every row is joined by candidate id.
python PhotonID/score_trees.py --system pp --config config/nominal.yaml \
    --model models/pp/pp_photon_bdt.root --input trees/pp_data.root \
    --input trees/pp_photonjet_sim.root --output-dir scored/

# Derive the tight threshold at the configured signal efficiency and paste the
# printed fragment into config/nominal.yaml (systems.pp.working_points).
python PhotonID/derive_working_points.py --system pp --config config/nominal.yaml \
    --input scored/pp_photonjet_sim.root --target-efficiency 0.80 \
    --output models/pp/working_points.yaml

# After editing thresholds, rerun score_trees.py without --model to refresh
# the flags from the stored score.
python PhotonID/score_trees.py --system pp --config config/nominal.yaml \
    --input trees/pp_data.root --output-dir scored/
```

`train_photon_bdt.py` writes a manifest with the ordered features, the label
mapping, the sample hashes and the held-out AUC; pin the model with
`model_file` and `model_sha256` in the configuration.

## 3. TreeToHists: nominal histograms

```bash
python TreeToHists/make_histograms.py --system pp --config config/nominal.yaml \
    --input scored/pp_data.root --output-json hists/pp_data.json --output-root hists/pp_data.root

# Prompt-photon leakage from the photon+jet simulation (truth-signal only)
python TreeToHists/make_histograms.py --system pp --config config/nominal.yaml \
    --input scored/pp_photonjet_sim.root --truth-signal-only \
    --output-json hists/pp_leakage.json
```

Output: ABCD event-leading photon counts per photon-pT bin and the xJgamma
spectra in regions A and C, with Sumw2, as JSON and ROOT histograms.

## 4. FinalAnalysis: purity correction, unfolding, overlay

```bash
# Response from the photon+jet simulation trees (2D: photon pT x xJgamma)
python FinalAnalysis/photonjet/cli.py response build --system pp --dimension 2D \
    --input scored/pp_photonjet_sim.root --output-stem response/pp

# Leakage-aware ABCD purity, corrected spectrum, iterative-Bayes unfolding
python FinalAnalysis/run_corrections.py --config config/nominal.yaml \
    --data hists/pp_data.json --leakage hists/pp_leakage.json \
    --response response/pp --iterations 4 --output results/pp.json

# Final comparison (one photon-pT bin per figure)
python FinalAnalysis/plot_overlay.py --result "p+p=results/pp.json" \
    --result "Au+Au 0-20%=results/auau.json" --pt-bin 0 --output results/xjgamma_pt15_20.png
```

`run_corrections.py` reports, per photon-pT bin, the ABCD solution, the
purity-corrected spectrum and, when a response is given, the unfolded spectrum
with data-statistics uncertainties and the refolding chi2/ndf.  The iteration
count is a configuration choice; set `unfolding.iterations` once the
refold/toy scan (`photonjet.analysis.unfolding.scan_problem`) has been run on
the nominal samples.

## Tests

```bash
python tests/run_tests.py                      # everything
python tests/run_tests.py test_chain_on_fixtures.py
```

The suite checks the selection boundaries, that TreeToHists reproduces the
independent reference reducers bin by bin on the fixture trees, that scoring
preserves the tree contract and joins by candidate identity, that the
correction and unfolding driver runs on a real response bundle, and the reused
kernel tests for the ABCD, response and unfolding arithmetic.

## What is still needed before nominal results

1. Collaborator trees produced with the current producer for data and
   simulation.  The fixture trees in `tests/` predate the R = 0.3 truth
   isolation branches and cannot define the nominal training signal.
2. Trained models: run `train_photon_bdt.py` once those trees exist; the
   configuration ships without a model.
3. Working points: derive the tight threshold with
   `derive_working_points.py`; choose and record the non-tight band and the
   R = 0.3 isolation thresholds in the configuration.
4. Unfolding iteration count from the refold/toy scan.
5. Optional combinatoric (unmatched recoil) subtraction and photon efficiency
   correction, both supported by the kernels but not applied in this pass.

## Layout

```text
config/nominal.yaml       one nominal configuration (pp and auau sections)
contracts/                collaborator tree branch contract
TreeProduction/           Fun4All producer (src, src_AuAu, macros) and tree builder
PhotonID/                 photon_selection.py, train_photon_bdt.py, score_trees.py, derive_working_points.py
TreeToHists/              make_histograms.py
FinalAnalysis/            run_corrections.py, plot_overlay.py, photonjet/ (analysis kernels and cli)
tests/                    unit and connected-chain tests, small fixture trees
tools/                    fixture generator used by the tests
```
