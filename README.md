# sPHENIX isolated photon + jet analysis

Code to go from sPHENIX DSTs to the unfolded photon-jet momentum balance
`xJgamma = pT(jet) / ET(photon)` in p+p and Au+Au collisions, organised as four
stages that hand off through documented files:

```text
DST                      TreeProduction/   Fun4All producer + tree builder
  -> collaborator trees  PhotonID/         train, score, derive working points
  -> scored trees        TreeToHists/      nominal selection, weights, ABCD and xJ histograms
  -> histograms          FinalAnalysis/    responses, purity correction, unfolding, overlay
  -> final xJgamma comparison plot
```

Two files hold the physics choices shared by every stage: `config/nominal.yaml`
(truth-photon definition, R = 0.3 isolation, kinematic windows, model inputs,
working-point thresholds, training recipe, correction settings) and
`config/samples.yaml` (cross sections, stitching windows, centrality maps).
Changing a threshold regenerates flags and histograms without retraining a
model or re-reading a DST.  `docs/ANALYSIS.md` explains the physics behind
every choice; `docs/STATUS.md` says what is validated and what is pending.

Status: **first-pass candidate, downstream stages complete, inputs pending.**
Every stage has a real command and the Python stages are exercised end to end
on synthetic inputs and the small fixture trees in `tests/`.  Models and
working points for the nominal selection are not yet derived, and the
producer in `TreeProduction/` will be refreshed from the frozen production
source when the current Au+Au production completes.

## Physics definitions

* **Signal photon (truth).**  A Geant photon with a valid generator
  association, photon class in {-1, 0, 1, 2} (direct, fragmentation, and
  unclassified; decay photons are class 3 and rejected), truth isolation
  energy in a cone of R = 0.3 below 4 GeV and |eta| < 0.7.  Missing
  associations are never admitted as signal.  This one definition is the BDT
  training target and the selected truth photon in simulation.
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
* **Weights.**  One complete weight per event, producer weight times source
  factor times centrality factor, applied exactly once; Sumw2 is its square.

## Requirements

* sPHENIX software environment with ROOT >= 6.26 and TMVA (PyROOT is used for
  scoring, working points and ROOT histogram output).
* Python 3.10+ with the packages in `requirements.txt`
  (`pip install -r requirements.txt`).  `xgboost`, `scikit-learn` and `scipy`
  are only needed for training.
* Building the producer needs the sPHENIX `coresoftware` install
  (`Fun4All`, `calobase`, `caloreco`, `jetbase`, `calotrigger`, `centrality`).

Every command below is run from the repository root.

## 1. TreeProduction: DST to collaborator trees

Build the two analysis modules once in an sPHENIX shell (see
`TreeProduction/README.md` for the environment and the external inputs):

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
contract is in `contracts/photonjet_trees_v1_branches.json`.

## 2. PhotonID: train, score, derive working points

```bash
# Train one model per system from photon+jet (signal) and inclusive-jet
# (background) simulation trees with the established recipe.
python PhotonID/train_photon_bdt.py --system pp --config config/nominal.yaml \
    --signal trees/pp_photonjet_sim.root --background trees/pp_inclusive_sim.root \
    --output-dir models/pp

# Score data and simulation with the same per-system model; flags and leader
# indices are regenerated from config/nominal.yaml, joined by candidate id.
python PhotonID/score_trees.py --system pp --config config/nominal.yaml \
    --model models/pp/pp_photon_bdt.root --input trees/pp_data.root \
    --input trees/pp_photonjet_sim.root --input trees/pp_inclusive_sim.root --output-dir scored/

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

## 3. TreeToHists: nominal histograms with the complete weights

```bash
python TreeToHists/make_histograms.py --system pp --config config/nominal.yaml \
    --sample pp_data --input scored/pp_data.root \
    --output-json hists/pp_data.json --output-root hists/pp_data.root

# Prompt-photon leakage from the photon+jet simulation (nominal truth signal only)
python TreeToHists/make_histograms.py --system pp --config config/nominal.yaml \
    --sample pp_photonjet --input scored/pp_photonjet_sim.root --truth-signal-only \
    --output-json hists/pp_leakage.json
```

For Au+Au embedded inclusive-jet slices give one line per file in an input
list, `path sample`, so each slice gets its ownership stitching and centrality
factor (`--sample` names one sample for every input).

## 4. FinalAnalysis: responses, corrections, unfolding, overlay

```bash
# Pair response (photon pT x xJ) and per-event photon response from the
# photon+jet simulation, with the nominal truth definition and sample weights
python FinalAnalysis/photonjet/cli.py response build --config config/nominal.yaml --system pp \
    --sample pp_photonjet --input scored/pp_photonjet_sim.root --output-stem response/pp_pairs
python FinalAnalysis/photonjet/cli.py photon-response build --config config/nominal.yaml --system pp \
    --sample pp_photonjet --input scored/pp_photonjet_sim.root --output-stem response/pp_photons

# Purity correction, combinatoric subtraction (Au+Au), unfolding of pairs and
# photons, toy uncertainties, iteration selection; final points as JSON + CSV
python FinalAnalysis/run_corrections.py --config config/nominal.yaml --system pp \
    --data hists/pp_data.json --leakage hists/pp_leakage.json \
    --pair-response response/pp_pairs --photon-response response/pp_photons \
    --output results/pp.json

# Final comparison, optionally with published reference points
python FinalAnalysis/plot_overlay.py --result "p+p=results/pp.json" \
    --result "Au+Au 0-20%=results/auau.json" \
    --reference "ATLAS p+p, 63-80 GeV=FinalAnalysis/reference/atlas_plb789_167_table1_xjgamma.csv:pp" \
    --output results/xjgamma_overlay.png
```

`FinalAnalysis/README.md` describes the chain step by step.

## Tests

```bash
python tests/run_tests.py                      # everything
python tests/run_tests.py test_chain.py        # one module
```

The suite covers the selection boundaries and truth contract, the complete
weight contract (source factors, ownership windows, centrality support), the
per-event photon response, the correction and unfolding chain with its gate
and iteration selection on synthetic self-consistent responses, the training
recipe (weighting, flattening, split, label mapping), scoring against the tree
contract, agreement of TreeToHists with the independent reference reducers on
the fixture trees, the producer scrub rule, and the tree-builder interface.

## What is still needed before nominal results

1. Collaborator trees produced with the current producer for data and
   simulation.  The fixture trees in `tests/` predate the R = 0.3 truth
   isolation branches; the nominal training signal needs the new trees.
2. Trained models: run `train_photon_bdt.py` once those trees exist; the
   configuration ships without a model.
3. Working points: derive the tight threshold with `derive_working_points.py`;
   choose and record the non-tight band and the R = 0.3 isolation thresholds.
4. `generated_events` for the p+p simulation samples in `config/samples.yaml`,
   and re-derived Au+Au centrality maps for the new production.
5. Systematic uncertainties, which this repository does not evaluate.

## Layout

```text
config/nominal.yaml        one nominal configuration (pp and auau sections)
config/samples.yaml        sample normalisation: cross sections, stitching, centrality maps
contracts/                 collaborator tree branch contract
docs/ANALYSIS.md           the physics of every choice; docs/STATUS.md what is validated
TreeProduction/            Fun4All producer (src, src_AuAu, macros) and tree builder
PhotonID/                  photon_selection.py, train_photon_bdt.py, score_trees.py, derive_working_points.py
TreeToHists/               make_histograms.py, sample_weights.py
FinalAnalysis/             run_corrections.py, plot_overlay.py, reference/, photonjet/ (kernels and cli)
tests/                     unit and connected-chain tests, fixture trees, synthetic tree helper
tools/                     producer scrub script, fixture generator
```
