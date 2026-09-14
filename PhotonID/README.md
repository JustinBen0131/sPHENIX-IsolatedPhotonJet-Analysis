# PhotonID: training, scoring, working points

* `photon_selection.py` is the one implementation of the truth-signal
  contract, the working-point thresholds and the tight / non-tight /
  isolated / non-isolated flags.  Every other stage imports it.
* `train_photon_bdt.py` trains one XGBoost model per collision system from
  photon+jet (signal) and inclusive-jet (background) trees with the
  established recipe (class balance, spline-based eta and ET flattening with
  the 800 cap computed before the training window, row-level stratified
  split with the per-system seeds, per-system hyperparameters) and exports it
  as a TMVA RBDT file plus a manifest with the label mapping, weighting steps,
  sample hashes and AUC.  Needs `xgboost`, `scikit-learn`, `scipy`.
* `score_trees.py` evaluates the model for every candidate with complete
  inputs, writes the configured thresholds into the threshold branches,
  regenerates the flags and the stored leader indices, and propagates
  everything to `photonJets` and `eventTree` by candidate identity.  Without
  `--model` it keeps the stored score and only refreshes thresholds and flags.
* `derive_working_points.py` finds the tight threshold that gives the
  configured signal efficiency per photon-pT bin (and centrality bin in
  Au+Au) from scored signal simulation and prints a YAML fragment for
  `config/nominal.yaml`.

Rules kept by these scripts: data and simulation are scored with the same
per-system model; every join uses the candidate identity; thresholds live in
the configuration, never in the code; a candidate with a missing score or
threshold is invalid rather than background; the producer event weight is not
used in training.
