# PhotonID

Train reproducible photon-identification models from completed single
`eventTree` datasets, validate them numerically, derive their working points
and add named scores to DATA and SIM. The original AuAu schema10 and pp PPG12
scores remain the default: `pho_bdt_score`. Optional training never overwrites
that branch or changes the nominal isolation definition.

## Commands

Use an already configured environment with Python, NumPy, SciPy, PyYAML,
uproot, awkward and ROOT, plus TreeProduction's installed eventTree library:
set `MYINSTALL` to its installation prefix (a `--disable-producer` build
needs only ROOT). BDT fitting additionally requires XGBoost. This package does
not install dependencies or activate a different environment.
`PHOTONID_PYTHON` may select the Python executable in your chosen environment.

```bash
./photon_id.sh --help
./photon_id.sh listModels AuAu
./photon_id.sh reTrainNominal AuAu
./photon_id.sh reTrainNominal PP
./photon_id.sh reTrainNominal AuAuAndPP
./photon_id.sh trainNew centAndETdepModels AuAu
./photon_id.sh trainNew MLP AuAu
./photon_id.sh trainNew logReg AuAu
./photon_id.sh score MLP_v1 AuAu
```

Training commands perform preparation, fitting, saved-model validation,
working-point derivation and complete-group score augmentation. `score` applies
an existing validated release to additional groups. Matching completed work is
reused. No plotting or scheduler submission occurs here.

Standard commands resolve the maintained published group paths in
`config/training.yaml`; ordinary users should not construct file lists. For a
custom completed dataset, repeat `--group /path/to/published/group`. Each group
must have the shared eventTree publication `CURRENT` pointer and immutable `MANIFEST.json`.
The exact source roster, sample role, system and file hashes are checked.

Useful options:

- `--check`: read-only input/recipe/dependency preview.
- `--models /path/to/models`: use an external writable model release root.
- `--release-version v2`: explicitly start a new release for changed inputs,
  recipe or environment; existing versions are immutable.
- `--no-score`: finish the validated model package without publishing scores.
- `--config /path/to/training.yaml`: explicit maintained recipe/deployment file.

The paired command uses the separate AuAu/PP collections in the configuration.
It reuses an AuAu release only when its input, recipe, code and environment
identity matches a completed, hash-verified release. A bundled frozen model is
not evidence that eventTree retraining has already happened.

## Model names and branches

| Selection name | Score branch | Definition |
|---|---|---|
| `nominal` | `pho_bdt_score` | Original schema10 in AuAu; original PPG12 in pp |
| `nominal_retrained_v1` | `pho_bdt_score_nominal_retrained_v1` | Shared nominal BDT recipe |
| `centCoarse_v1` | `pho_bdt_score_cent_coarse_v1` | Three centrality cell models |
| `centFine_v1` | `pho_bdt_score_cent_fine_v1` | Seven centrality cell models |
| `ET_v1` | `pho_bdt_score_et_v1` | Eight ET cell models |
| `ETcentCoarse_v1` | `pho_bdt_score_et_cent_coarse_v1` | 24 ET/centrality cell models |
| `ETcentFine_v1` | `pho_bdt_score_et_cent_fine_v1` | 56 ET/centrality cell models |
| `MLP_v1` | `pho_mlp_score_v1` | Global nominal-input neural network |
| `logReg_v1` | `pho_logreg_score_v1` | Global nominal-input logistic regression |

Every optional release adds **one** photon-aligned score branch. Values outside
its applicability domain are NaN. Missing required inputs or a nonfinite
prediction **inside** the domain block publication; there is no new state
branch or silent fallback model. The original nominal state/class branches
remain unchanged.

Model selection must bind the score, domain, model hashes and its own working
points together. All AuAu optional models use strict ID comparisons:
non-tight `0.1 < score < WP80`, tight `score > WP70`. The band between WP80 and
WP70 is excluded; values at or below 0.1 fail ID. Original nominal models keep
their frozen definitions. Newly trained pp working points are fitted against
ET; the package does not copy original PPG12 numerical cuts onto new scores.

## Shared training definition

- Embedded AuAu Photon12/20 supply signal with photon class 1 or 2. Inclusive
  Jet12/20/30/40 supply valid nonprompt background. Cross-role and unresolved
  classifications are discarded. Labels use the dominant-primary barcode and
  truth photon inventory, independently of analysis-isolation matching.
- AuAu embedded events must pass the corrected-replay minimum-bias witness,
  `cent_mb_decision == 1`, as in the maintained eventTree training projection.
  Missing/unknown decisions and failures are rejected. The retired
  `evt_embedded_mb_decision` is retained as provenance, never used as a fallback.
- Require the retained, ready, valid-vertex event predicate, `15 <= ET < 35`
  GeV, `|eta| < 0.7`, `|z| < 10` cm; AuAu centrality is `[0,80)` percent.
- No truth or reconstructed isolation requirement enters training.
- Use the fourteen ordered H70 AuAu features. pp omits weta33, wphi33 and
  centrality, retaining the eleven PPG12 base features. Named branches must
  agree with the stored H70 feature vector and its semantic hash.
- Prepare class balancing times the canonical spline ET/eta flattening once,
  before splitting or routing. Do not add event, cross-section, centrality or
  vertex weights. No per-cell rebalance, extra ET weights or row caps.
- The source-qualified event-grouped split is 90/10, seed 13. Every learner
  and cell uses the same saved event assignments; different candidates of one
  event cannot cross partitions. This reproduces the method on current
  eventTrees, not the exact fitted weights or old training population.

The MLP starts with 48/24 tanh hidden nodes and a sigmoid output, Adam learning
rate 0.001, L2 0.0001, batches of 4096 and 80 fixed epochs. Logistic regression
uses a linear sigmoid, Adam learning rate 0.018, L2 0.002, batches of 262144 and
60 fixed epochs. Scaling is fitted on training rows only and saved with the
model. Logistic inputs retain standardized clipping at 8; MLP clipping,
jitter, distillation, searches and temperature fitting are disabled. Neither
learner selects epochs using the outer holdout. These starting recipes need
qualification on the intended production populations.

Routed BDTs use half-open bins:

- ET: 15, 17, 19, 21, 23, 25, 28, 32, 35 GeV.
- Coarse centrality: 0, 20, 50, 80 percent.
- Fine centrality: 0, 10, 20, 30, 40, 50, 60, 80 percent.

The five families contain 98 constituent models but produce only five score
branches. Missing classes, failed cells or incomplete routing block the
affected family; cells are never silently dropped.

## Validation and publication

Numerical checks cover population identity, split/class coverage, weighted
train/holdout discrimination and saved-artifact score agreement. The configured
AUC/gap gates are engineering starting criteria; a PASS does not establish
physics closure or promote a model to the analysis baseline.

WP70/WP80/WP90 use **all current eligible signal**, including train and holdout,
with canonical training weights. The historical 4,665,004 candidates are not a
required new count. AuAu lines are fitted versus centrality. Each routed cell
uses its own calibration; a single centrality-bin cell has a constant local
threshold. Efficiency, line residuals and cut ordering are checked. There is
no separate untouched evaluation partition after full-signal WP calibration.

Input ROOT reads and writes are bounded. The population matrix is disk-backed;
SQLite checks source/event/photon identities. The eventTree library writes
each score candidate, then reopens it and compares every original payload
branch and carried object; every score is recomputed independently.
Publication uses the same library's complete-group roster validation, DATA
event uniqueness and atomic compare-and-swap pointer, after rehashing every
file against its proof receipt. A partial augmentation leaves the published
group unchanged, with the same folder/file count. Failed candidates are
retained for inspection. Never delete or overwrite them to force progress.

## Files and generated releases

| File | Responsibility |
|---|---|
| `photon_id.sh` | Public command entrance |
| `scripts/photon_id.py` | Readable phase coordinator and reuse |
| `config/training.yaml` | Shared recipes and deployment group bindings |
| `src/Inputs.py` | H70 projection, labels, weights and event split |
| `src/Training.py` | BDT, NumPy MLP and logistic fitting |
| `src/Validation.py` | Numerical and saved-model checks |
| `src/WorkingPoints.py` | Model/cell efficiency thresholds |
| `src/Scoring.py` | Aligned score evaluation and shared publication adapter |
| `src/Models.py` | Immutable releases, hashes and model names |

`models/AuAu/nominal/` and `models/PP/nominal/` describe and carry the bundled
original models. Generated releases are under `models/<system>/<name>/` or
`--models`; they are gitignored. Each contains a manifest, model artifacts,
working points, validation and a downstream selection definition. Preparation
arrays and preparation progress stay in the ignored `models/.populations/`
cache, shared by learners with identical input preparation and event splits.

`models/nominal.json` pins the shared original definitions used by TreeProduction and
TreeToHists. Its bytes, model hashes and working points are preserved from the
original package; actual payload paths come from each system's nominal
manifest. For pp the bundled snapshot is hash-pinned; CDB identity, conditions
and correspondence must be recorded from actual retrieval evidence before
claiming CDB qualification. No query of a live “latest” CDB model may silently
change the nominal payload.

## Current qualification boundary

Published dataset groups are deployment bindings and are not guessed while
storage is unavailable. The pp training sample population remains explicitly
unbound pending its approved manifest. BDT fitting refuses an environment
without XGBoost. Local tests do not establish farm or production qualification.

The emitted `selection.json` supplies one-score, domain and routed-WP bindings.
TreeToHists consumes that sealed file with `systems.auau.photon_id_selection`
(or `systems.pp.photon_id_selection` for a retrained pp nominal model) in its
selection YAML. Omit score/state/working-point overrides for that system;
the release supplies them. The reader checks the saved score's release/domain
receipt and routes working points without a separate state branch. Local
adapter qualification does not establish production model qualification.
PhotonID supplies the definition and scores; TreeToHists owns histogram reduction.
