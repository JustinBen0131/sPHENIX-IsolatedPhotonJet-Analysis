# PhotonID

Consumes completed base trees and produces frozen model artifacts or candidate-keyed score sidecars. It never rewrites base trees or assigns tight/non-tight/ABCD regions.

`features.py` owns the ordered feature vocabulary, shower view, ratio policy, domain and finite-input rules used by both `train.py` and `augment.py`. The feature identity hashes that definition and the implementation bytes. Inputs are narrowed to the actual float32 inference representation before validity checks.

## Registry

Five distinct **unbound** slots exist: `ppg12_npb_reference`, `ppg12_tight_pp_reference`, `ppg12_tight_auau_reference`, `canonical_pp_v1`, `canonical_auau_v1`. Reference and canonical models must not be conflated. Feature lists for unbound entries are proposals requiring model-specific parity review. Canonical training targets 15–35 GeV; application-domain coverage of wider response support remains unresolved.

A bound registry entry needs an approved status, binding evidence, exact model file and SHA256. A CDB key records a separately resolved import; Python inference does not silently download one. No model weights were trained, imported or published in this pass.

## Training

Training is blocked until the registry supplies an explicitly approved complete recipe: label rule, ET window, class/kinematic weighting, optional thinning parameters, split fractions/seed/rule and XGBoost hyperparameters. The currently supported recipe uses known truth relations and source sample roles, fits weights/thinning on training observations only and evaluates held-out AUC with unit weights. This supported implementation is **not** an accepted pp/AuAu training recipe by itself.

Splits hash source/event identity rather than candidate row order. Candidates from a shared event stay together. Alternative representations with different source/event keys require a reviewed common physical-event grouping before mixing; the code cannot infer those aliases. Incomplete truth censuses and duplicate candidate identities are refused.

## Sidecars

`augment.py` writes one `PhotonScores` TTree per base file/model. Rows carry source/event/candidate identity, model name, raw score and states for finite evaluation, nonfinite evaluation, missing inputs, out of domain and unavailable model. Metadata binds the exact base-file SHA256, model SHA256, feature-definition SHA256 and input recipe. No join may rely on row position. An unbound model produces only unavailable-model diagnostics, never a usable score product.

Outputs must be new. Base mutation during scoring is rejected. Training/inference parity, TMVA export, numerical score behavior, model bindings and working-point derivation remain runtime/science gates. The scripts currently load a file's candidate features in memory; only the histogram loop has an event-batched implementation.
