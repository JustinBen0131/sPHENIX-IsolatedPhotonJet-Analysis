# Photon identification

Train candidate models, then attach scores and selection decisions as keyed ROOT
sidecars. **Base trees remain immutable and model-independent.** New models do
not require DST production or changes to the base-tree schema.

From PhotonID:

```sh
./train_photon_id.sh trainAuAu
./train_photon_id.sh trainPP
./train_photon_id.sh trainPPG12Equivalent
./train_photon_id.sh trainAll

./augment_photon_id.sh runPP
./augment_photon_id.sh runAuAu
./augment_photon_id.sh runAll
```

Both entrypoints print help without arguments or with --help. Fish users invoke
these executable bash scripts normally. Aggregates reuse primitive executors;
trainAll currently runs sequentially to bound local memory.

**Readiness: local infrastructure and synthetic fixtures only.** All supplied
production profiles remain unbound. No production model was trained or accepted.
Commands fail with a nonzero error until their campaign, label and model bindings
are supplied. Canonical names express intended identities, not acceptance.

## Installation and options

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
./train_photon_id.sh --help
./augment_photon_id.sh --help
```

The wrappers find this local .venv automatically; PHOTONID_PYTHON can override it.
Help needs only the standard library. PyROOT is needed only for bound historical
TMVA models; the new XGBoost JSON path uses uproot for ROOT I/O.

```sh
./train_photon_id.sh trainAll --dry-run
./augment_photon_id.sh runPhotonJetSim --models canonical_pp_v1 --source 1 --dry-run
./augment_photon_id.sh runPP --models canonical_pp_v1,ppg12_equivalent_v1
./augment_photon_id.sh runAuAu --selection nominal_auau
./augment_photon_id.sh runAll --models all-compatible
```

--models requests raw scoring independent of recipes. --selection and --selections
are synonyms and accept comma-separated names. Both model and recipe options may
be supplied together. Without either, the system's nominal recipe is required.
Model all-compatible selects bound compatible models and prints unbound exclusions.
Explicit incompatible/unbound names fail. Recipe all-compatible requires all
selected compatible recipes to be complete.

--source is one-based within a primitive lane and only narrows input. --dry-run
checks bindings and displays the plan without writing. Training dry-run also
reads schemas/features/witnesses. Augmentation dry-run checks manifest hashes,
models and recipes; full base-schema checks occur during execution. No batch
submission mode exists.

--campaign PATH --campaign-sha256 SHA is an explicit development override.
--registry, training --profiles, and --output-root allow portable bindings.
Document paths resolve relative to their configuration/manifest. Existing outputs
are never overwritten. Failed training retains INCOMPLETE.json for inspection
before a deliberate retry/version change.

## Lane and profile ownership

| Augmentation command | Completed TreeProduction lane |
|---|---|
| runPP | pp_data |
| runAuAu | auau_data |
| runPhotonJetSim | pp_photon_sim |
| runInclusiveJetSim | pp_inclusive_sim |
| runEmbeddedPhotonJetSim | auau_photon_embedded |
| runEmbeddedInclusiveJetSim | auau_inclusive_embedded |
| runAllData | runPP + runAuAu |
| runAllSim | four simulation primitives |
| runAll | runAllData + runAllSim |

training_profiles.yaml owns training policy and expected sample categories.
model_registry.yaml owns model, ID-WP, isolation and selection bindings.
Actual populations belong to TreeProduction. Its production_populations.json
describes DST inputs, **not completed outputs**. The expected completed-tree
interface is in [CONTRACTS.md](CONTRACTS.md); no population is fabricated here.

| Identity | Features | Split | Weight fit population | XGBoost | Binding |
|---|---:|---|---|---|---|
| canonical_auau_v1 | 14 H70: BaseV3E + widths + centrality | physical-event 80/10/10, seed 13 | TRAIN | 450 trees, depth 4, rate .035 | campaign/labels unbound |
| canonical_pp_v1 | 13 H70: BaseV3E + widths | physical-event 80/10/10, seed 13 | TRAIN | same 450-tree recipe | campaign/labels unbound |
| ppg12_equivalent_v1 | 11 BaseV3E | row-stratified 70/10/20, seeds 42/43 | all before split, historical exception | 750 trees, depth 5, rate .1 | domain/label/campaign adapter unbound |
| ppg12_original | 11 BaseV3E | historical reference, no retraining | historical | frozen TMVA | expected SHA only; artifact/import receipt unbound |
| historical_auau_h70_the134 | 14 H70 | historical 90/10 | historical | frozen 450-tree model | expected SHA only; artifact/import receipt unbound |

Weights are class balancing times per-class ET/eta inverse-density splines:
PDF floor 1e-3, ET-weight cap 800, no event/cross-section/vertex/centrality weight.
Canonical splines and normalizations are fitted on TRAIN only and frozen.
Parity intentionally retains historical full-population weighting and per-source
background flattening below 15 GeV; it is not event-grouped or leakage-free validation.

Canonical parameters additionally use subsample=.85, colsample_bytree=.85,
reg_alpha=5, reg_lambda=.3, tree_method=hist, grow_policy=lossguide, max_bin=256,
seed 13 and two threads. Parity uses subsample=.5, colsample_bytree=.6,
colsample_bylevel=1, seed 42 and four threads. No early stopping, scaler or tuning.

Canonical domains: 15 <= ET < 35 GeV, |eta| < .7, and for AuAu 0 <= c < 80.
Parity's training domain remains unbound: the registry's 15–35 scoring window must
not silently truncate historical training. Bind the approved full domain
consistently in profile and registry (and original-model adapter for comparison).

Connected physical-event identities prevent candidates sharing either the hard
event or an embedded underlying event from crossing canonical partitions. Missing
identity witnesses fail closed. Source/event output keys alone do not prove
physical independence. Assignments, groups, row order, seed and counts are saved.
A changed population is new split provenance; cross-population membership
stability is not promised.

The named dominant-primary prompt adapter defines BDT signal independently of
reconstructed isolation; it never uses TruthPhotons.analysis_signal for BDT labels.
Unknown ancestry is excluded. Exact historical-label equivalence remains an
explicit binding. The analysis-signal flag needs a separate reviewed binding for
isolation calibration.

## Calibration and recipes

TRAIN fits the BDT/weights. VALIDATION alone derives WP70/80/90 and, when its truth
denominator is bound, ISO70/80/90. FROZEN_CHOICES.json precedes test prediction.
TEST evaluates frozen choices; ALL is labeled **full-sample / not held-out**.

AuAu BDT calibration saves six pT bins [15,17,19,21,23,26,35] in each 5% centrality
slice, constant pT fits, pooled slice quantiles and their linear centrality fit.
pp uses a declared linear ET fit. Historical weighted quantiles/errors are reused;
old coefficients are not. Crossing/out-of-range curves fail instead of clipping.

AuAu isolation consumes retained SUB1 R=.4 with candidate removal recorded.
Its independently selected truth-isolated analysis-signal population is unit
weighted, with no BDT completeness/score cut. The recovered estimator uses a
flow-aware 700-bin histogram on [-20,50] GeV, constant pT fits, then linear
centrality fits. Quantiles in under/overflow fail.

| Recipe | Model | ID package | Isolation | Unresolved |
|---|---|---|---|---|
| nominal_pp | ppg12_equivalent_v1 | its own validation WPs | ppg12_topocluster_r04_v1 | model and model-specific WPs |
| nominal_auau | canonical_auau_v1 | its own validation WPs | canonical_auau_sub1_r04_v1 | model/WPs, denominator, nominal efficiency and sideband |

pp isolation is versioned once: isolated I < .490 + .037 ET, non-isolated
I > 1.290 + .037 ET; equality/intermediate values are GAP.
Historical AuAu executable evidence uses strict I > sideband_threshold, and a
completion report describes zero gap. This does not establish the new ISO90 alias.
After review, same_as_isolated supports that zero-gap binding. No numeric AuAu
sideband has been invented.

ID: score > WP70 is TIGHT; .1 < score < WP80 is NONTIGHT;
WP80 <= score <= WP70 is EXCLUDED; score <= .1 is BELOW_FLOOR.
Invalid inputs remain INVALID. WP90 is retained. TreeToHists will combine these
regions with isolation to form ABCD; ABCD is not a model result.

An unbound analysis denominator omits isolation calibration/its two plots and
writes isolation/UNBOUND.json; otherwise-bound model/ID training still works.
An unbound nominal sideband likewise does not block scoring or isolated calibration.

## Output

```text
output/photon_id/
  training/<model identity>/
    training_profile.json
    model/model.json
    split/{assignments.npz,weighting.json}
    working_points/id.json
    isolation/calibration.json       # AuAu, denominator bound
    isolation/selection.json         # nominal/sideband also bound
    isolation/UNBOUND.json           # instead, if denominator unbound
    FROZEN_CHOICES.json
    diagnostics/{validation,test,all}/
      payload.npz
      numerical.json
      01_score_distributions_and_roc.png
      02_bdt_wp_pt_centrality_grid.png  # pp: 02_bdt_wp_pt_grid.png
      03_bdt_wp_centrality_fits.png     # pp: 03_bdt_wp_pt_fits.png
      04_bdt_score_regions.png
      05_isolation_wp_pt_centrality_grid.png  # AuAu only
      06_isolation_wp_centrality_fits.png     # AuAu only
      07_ppg12_reference_comparison.png       # bound comparison only
    TRAINING_RECEIPT.json
  augmentation/<campaign>/<lane>/<configuration hash>/
    <source identity>.photon_id.root
```

Validation grids show quantiles/errors. Test/all grids show achieved efficiencies
at frozen functions; fit panels explicitly display the validation reference,
never a refit. numerical.json contains fit tables and efficiency/error reports.
Render again without a model:

```sh
.venv/bin/python diagnostics.py /path/to/package/diagnostics/validation
```

Receipts bind model/features, configuration, inputs/witness hashes, splits,
weighting, hyperparameters, metrics, software commit and source hashes, library
versions, calibrations and diagnostics. BUILD/MECHANICAL PASS is separate from
SCIENTIFIC/CANONICAL ACCEPTANCE: NOT_REVIEWED. Nothing auto-promotes or changes
aliases. Registry imports require exact artifact/feature/receipt hashes and
model-specific WPs; historical references require a reviewed import receipt.

## Files and tests

Originally this directory contained README, features.py, train.py, augment.py and
model_registry.yaml. Those are retained/refined, with focused modules:

- features.py: shared ordered training/inference adapter.
- inputs.py: source joins, prompt labels and retained isolation.
- train.py: weights, splits and fit/freeze/evaluate lifecycle.
- calibration.py + numerics.py: calibration, region codes, recovered pure functions.
- diagnostics.py: saved numerical payloads and rendering.
- augment.py: evaluators and long-form ROOT serialization.
- registry.py: bindings, compatibility, composition.
- frontend.py: operational CLI for the two public shell scripts.

No alternate framework, manager hierarchy or historical runtime imports.
[PROVENANCE.md](PROVENANCE.md) and [HISTORICAL_SOURCES.json](HISTORICAL_SOURCES.json)
record predecessors. [CONTRACTS.md](CONTRACTS.md) defines interfaces. The old
single-model CLI/schema is superseded by PhotonIDSidecarV2. Supporting it in the
old TreeToHists reader is future work; no downstream compatibility is claimed.

```sh
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python -m compileall -q .
bash -n train_photon_id.sh augment_photon_id.sh
```

Synthetic trees, calibration fixtures and tiny XGBoost fits do not establish
SDCC, PyROOT/TMVA, Fun4All, population, physics or downstream closure.

Training and augmentation currently hold one selected population/file in memory.
Large-population memory use and streaming are not certified by these tiny tests.
