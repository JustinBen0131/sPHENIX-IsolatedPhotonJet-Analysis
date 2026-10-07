# PhotonID

Train the AuAu and pp photon-identification models used by the analysis.
PhotonID produces portable model releases, working points and saved training
diagnostics. TreeToHists evaluates a selected release during its normal tree
read; FinalAnalysis draws the saved diagnostics and analysis QA.

## Setup

Use Python with NumPy, SciPy, PyYAML, uproot, awkward and ROOT. BDT fitting
also needs XGBoost and scikit-learn, plus ROOT's `SaveXGBoost`/`RBDT` interface.
The repository's [environment.yml](../environment.yml) supplies the local
analysis environment. From the repository root:

```bash
conda env create -f environment.yml
conda activate photonjet
python TreeProduction/scripts/check_environment.py
```

An existing compatible environment can be used instead; the commands with
`--check` below report missing training dependencies before opening input trees.
Set `MYINSTALL` to TreeProduction's installed eventTree
library prefix; a `--disable-producer` build is sufficient.
`PHOTONID_PYTHON` can select your Python executable.

For reference training on SDCC, use the established Python 3.13 / ROOT 6.32.06
release environment and the numerical versions in
[config/training-requirements.txt](config/training-requirements.txt).
Reuse an existing environment with these versions. If it needs restoring,
create the prefix on writable analysis storage:

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.562
mkdir -p /path/to/analysis/pip-staging
export TMPDIR=/path/to/analysis/pip-staging
export PIP_NO_CACHE_DIR=1
python3 -m venv --system-site-packages /path/to/analysis/training-env
/path/to/analysis/training-env/bin/python -m pip install --no-compile -r PhotonID/config/training-requirements.txt
export PHOTONID_PYTHON=/path/to/analysis/training-env/bin/python
export PYTHONPATH="/path/to/analysis/training-env/lib/python3.13/site-packages${PYTHONPATH:+:$PYTHONPATH}"
```

The path ordering keeps the pinned numerical packages ahead of the release's
packages while retaining its ROOT bindings. On Linux the requirements use the
CPU distribution of XGBoost 3.2.0; the learner settings are the same.

Run fitting on an execute node. The public command never installs packages
or changes the selected environment itself.

## Reference training structure

The nominal AuAu model retains the fourteen ordered H70 inputs and source-role
labels. Its fixed learner uses 450 trees, depth 4, seed 13 and two CPU threads;
all settings live in [config/training.yaml](config/training.yaml).
Class balance and twenty-bin ET/eta spline weights are computed once over the
complete selected range, before any routing-cell slice. Every routed learner
reuses those global weights and the nominal event split.

Whole events use the original lane, sample, run, segment and producer sequence
to form the historical 128-bit identity. A seeded SHA256 ordering selects the
nominal 10% holdout, with the event identity resolving hash ties. Publishing a
new file or reordering the source list does not change that identity.
Training saves weighted train/holdout AUC and log-loss at every round and uses
all 450 rounds. The established XGBoost 3 export adapter supplies ROOT's scalar
configuration view without altering the saved booster or fitted predictions.
Saved ROOT scores are checked against XGBoost before calibrating working points.
Weighted signal quantiles and the accepted centrality-line fit produce
WP70/80/90; tight uses WP70 and non-tight uses `0.1 < score < WP80`.

## Run

Copy [config/training.yaml](config/training.yaml) to an analysis configuration
and fill `datasets.AuAu` and `datasets.PP` with the completed group paths.
The supplied configuration leaves these locations unbound. From the repository
root, inspect the selected groups once:

```bash
./PhotonID/photon_id.sh --help
./PhotonID/photon_id.sh listModels AuAu
./PhotonID/photon_id.sh planCampaign AuAu
./PhotonID/photon_id.sh inspectCampaign AuAu --config /path/to/training.yaml --models /path/to/analysis/models
```

Use the generated `inputs.yaml` in the remaining stages. These commands assume
an existing physical data-staging directory on the worker; use that worker's
`_CONDOR_SCRATCH_DIR` instead of `--scratch-root` when available.

```bash
./PhotonID/photon_id.sh projectCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models --check
./PhotonID/photon_id.sh projectCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models --shard-index 0 --scratch-root /path/to/worker/data-staging
# After every shard in the projection plan has finished:
./PhotonID/photon_id.sh verifyProjections AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models
./PhotonID/photon_id.sh prepareCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models
./PhotonID/photon_id.sh trainCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models --check
./PhotonID/photon_id.sh trainCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models
```

Training consumes verified preparation, fits models, validates the saved artifacts,
derives working points and saves diagnostics. Matching completed work is reused.
For bounded batch jobs, select one recipe per invocation after preparation:

```bash
./PhotonID/photon_id.sh trainCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models --model nominal_retrained_v1
./PhotonID/photon_id.sh trainCampaign AuAu --config /path/to/training.yaml --inputs /path/to/inputs.yaml --models /path/to/analysis/models --model wideET_v1
# Run the remaining recipes listed by planCampaign, then collect their saved QA.
./PhotonID/photon_id.sh collectQA AuAu --models /path/to/analysis/models
```

`--model` can be repeated and preserves the campaign's recipe order. Each selected
recipe reuses complete preparation and keeps all of its routing cells together;
it does not change the learner, weights, event split or working points.
Independent `--model` invocations can run in parallel after preparation. They
reuse the verified arrays on disk and each fit the complete selected recipe.
Selected invocations seal their own releases; `collectQA` requires the complete
campaign before writing its comparison index.
`--check` checks the selected runtime before opening input groups and previews
source and recipe bindings. Training commands leave the tree groups unchanged.
The separate `score MODEL SYSTEM` command is available when a published tree
copy with an additional score branch is explicitly needed.

`inspectCampaign` verifies the stage-2 source roster and reports a sealed
`inputs.yaml` under `<models>/.inputs/<identity>/`. Keep its generated
`inputs.json.gz` alongside it; YAML indexes the compressed source/provenance
record. Compression preserves the source identities and content digest, and
the reader also accepts existing uncompressed `inputs.json` catalogues.
`--inputs` reuses that exact roster without inspecting TG again and implies
`--no-score`.
When the producer's saved `population.json` or `populations/<sample>.sources.tsv`
is available, inspection checks its hash against the publication, matches every
logical source and the event total, and uses its identities and counts without
opening tree headers. Add `--require-roster` to require this metadata-only path;
missing or incompatible producer rosters then stop before tree access. Without
that option, older publications use tree metadata. Projection still
checks each tree's payload hash, source identity, collision system, contiguous
event coverage and exact entry count before accepting its cached rows.

`projectCampaign --check` plans contiguous whole-source shards, each with at
most 100,000 events. `--max-events` can lower that bound. Run each planned
`--shard-index` once through your execution environment; the command executes
one shard and performs no fanout. Source extraction still needs measured
runtime, memory and storage qualification. A source larger than the bound is
refused pending qualified source-range partitioning.

Projection reads each source's training fields once over the campaign support, preserving
original source-entry/photon IDs and the historical event identity, and excluding
isolation and training weights.
It saves hash-bound, memory-mapped raw records under `<models>/.projections/`.
The default packed layout creates one binary row file, one immutable source
index, one lock and one terminal progress record per shard. The index retains
each source's receipt and byte interval. Source extraction uses the declared
data-staging directory; completed staging is removed after bundle verification,
and interrupted staging is retained. No scripts, executables or builds are
created there. System `/tmp` is not chosen automatically. `--cache-layout source`
retains the older per-source layout, and preparation can read either layout.
Each source payload is hashed before projection; inspection does not repeat
that full-file read.
Use `--projection-root DIR` to choose another shared cache. Verified parts are
reused; incomplete or changed parts are refused.

`verifyProjections` runs the same complete-cache verification used by
preparation and writes `<models>/<system>/projection_census.json`. The report
binds the input catalogue, source and photon counts, record bytes, and verified
source receipts. Missing sources, duplicate bundles or changed bytes fail
before writing the report. This command requires `--inputs` and reads only the
cached rows and catalogue. It does not reopen trees or fit a model. Its `--check`
preview lists the expected source population without verifying cache contents.

`prepareCampaign` requires complete projection coverage and creates the
two campaign populations and saves `<models>/<system>/preparation_census.json`.
This report records the prepared array sizes and the actual background/signal
counts in each model's training, validation and holdout cells. It supports
resource sizing before batch fitting and performs no fits itself.
Preparation and fitting then read GPFS products without reading
or hashing TG trees. A later `trainCampaign` with the same sources, code and
`--models` root verifies and reuses the populations. Inspection, projection
and preparation need the eventTree library and ROOT, but not XGBoost.

[config/training.yaml](config/training.yaml) holds the recipes and standard
group bindings. For custom inputs, repeat `--group /path/to/published/group`.
Use `--models DIR` for a writable release root and `--release-version v2`
for an explicitly new release.

For pp, use the same inspect → project → prepare → train commands with `PP`
instead of `AuAu`. Both pp models use the eleven ordered PPG12 features and
the source-role labels; truth or reconstructed isolation does not enter training.
The standard source pool is photonjet5/10/20 and jet8/12/20/30/40. The supplied
Jet8 sample replaces the absent historical Jet5 input, as in the corrected
PPG12 comparison. Updated trees and source supply are recorded separately
from the historical training mechanics in each release.

| pp model | Support | Learner and preparation | Use |
|---|---|---|---|
| `nominal_retrained_v1` | 15–35 GeV | 450 trees, depth 4; event-grouped 90/10; canonical ET/eta weights | Comparison with the canonical AuAu model |
| `ppg12_retrained_v1` | 6–40 GeV, historical closed acceptance | 750 trees, depth 5; per-sample low-ET flattening; stratified 70/10/20 rows, seeds 42/43 | Comparison with the bundled frozen PPG12 model |

The PPG12 validation and test rows remain distinct. Test rows alone supply its
held-out plots. Background flattening groups all logical files of a sample
before subsampling. Its weights use the same twenty-bin ET/eta splines, with
the historical float32 training-weight conversion.

After both systems finish, bind the pp/AuAu comparison without fitting again:

```bash
./PhotonID/photon_id.sh collectQA PP --models /path/to/analysis/models --auau-reference /path/to/analysis/models/AuAu/nominal_retrained_v1
```

This emits a hash-bound `models.yaml` for FinalAnalysis. `trainCampaign PP`
also accepts `--auau-reference` when that release is already available.

`planCampaign AuAu` describes 14 configurations: the global 15–35 GeV
depth-4 baseline, a global 12–40 GeV variant, four centrality/ET routing
families, and depths 2–10. Depth 4 reuses the nominal release. The plan opens
no trees and qualifies no inputs. `routingStudy` selects C3, C7, ET×C3 and
ET×C7; `centAndETdepModels` also includes the ET-only family.

Use `--no-score` to fit and validate models without changing the input groups.
`trainCampaign` always leaves the input groups unchanged. It resolves the
stage-2 sources once and reuses verified preparation across routing and depth.
Keep `--models` on writable analysis storage. The 12–40 GeV preparation retains
every nominal event's training/holdout assignment; additional events use the
nominal hash cutoff. Its receipt records the wider holdout fraction and verifies
zero differences for common events. Each model trains, calibrates its working
points and applies over its own full range, with its own class-balance and
flattening weights:

| Model | Training, working-point calibration and application | Analysis results | Underflow / overflow support |
|---|---|---|---|
| `nominal_retrained_v1` | 15–35 GeV | 17–33 GeV | 15–17 / 33–35 GeV |
| `wideET_v1` | 12–40 GeV | 14–38 GeV | 12–14 / 38–40 GeV |

Input inspection reports checked files in `<models>/.progress/inputs.json`.
Before preparation or training, `.progress/verify-projections-auau.json` (or
`verify-projections-pp.json`) counts verified source caches and reports their
total photon rows. Verification emits progress while hashing packed row files.
Each projection shard's `progress.json` counts source events, including events
with no selected photons; verified cache reuse is a deterministic skip with
separate counters. Each population's `progress.json` counts cached photon
visits across its two GPFS passes. `eligible_photons` and
`matrix_rows_written` remain separate counts. The record reports its phase,
source, remaining work, rate, ETA, peak memory and terminal errors; it is written
by the processing invocation. Event visits reaching their total does not finish
the event split, weighting or artifact sealing phases.

For a controlled job, supply `--progress-context job.yaml` containing
`task_id`, `workstream_id`, `campaign_tag` and `executor` as nonempty strings.
Fitting reports completed boosting rounds for each model cell, with the same
job identity carried through native-score checks, calibration and saved QA.
Saved QA uses the counter name `qa-<system>-<model>` and retains the exact
training request hash in its progress record.
Photon checks and calibration bins have separate phase counters.
`collectQA` counts verified model releases in `.progress/collect-auau.json`
or `.progress/collect-pp.json`; a bound AuAu comparison adds one release to the
pp total. Its `--check` preview remains read-only.
These labels affect telemetry only. Keep preparation, models and logs on writable
analysis storage; input trees are read without modification. Work totals and
measured runtime must fit the execution environment before a full campaign.

On shared storage, set `PHOTONID_PROGRESS_MODE=stdout`. The processing invocation
emits `PHOTONID_PROGRESS` JSON records for its controller to mirror; it writes
one terminal record and preserves an existing record on verified reuse. This
avoids periodic progress-file rewrites on GPFS. Local interactive work uses
atomic progress files by default. Keep compiled code, the installed reader and
compiler temporary files in the declared persistent build workspace. Data
staging is separate from `TMPDIR` and never changes the compiler environment.

AuAu and canonical pp ranges include the lower endpoint and exclude the upper
endpoint. The PPG12 comparison retains its historical closed acceptance and
does not prescribe measurement or guard bins. The recipe,
campaign plan and `selection.json` record these ranges as `et_ranges_gev`;
working-point metadata also records the calibration domain. No interior
analysis bin edges are prescribed by this range declaration. A comparison on
the common 15–35 GeV interval is an additional diagnostic and does not restrict
the wider model's calibration or application.

## Models

| Model | Score branch |
|---|---|
| `nominal` | `pho_bdt_score` |
| `nominal_retrained_v1` | `pho_bdt_score_nominal_retrained_v1` |
| `ppg12_retrained_v1` (pp) | `pho_bdt_score_ppg12_retrained_v1` |
| `wideET_v1` | `pho_bdt_score_et12to40_v1` |
| `depth2_v1`, `depth3_v1`, `depth5_v1` … `depth10_v1` | `pho_bdt_score_depth_2_v1` … `pho_bdt_score_depth_10_v1` |
| `centCoarse_v1`, `centFine_v1` | `pho_bdt_score_cent_coarse_v1`, `pho_bdt_score_cent_fine_v1` |
| `ET_v1` | `pho_bdt_score_et_v1` |
| `ETcentCoarse_v1`, `ETcentFine_v1` | `pho_bdt_score_et_cent_coarse_v1`, `pho_bdt_score_et_cent_fine_v1` |
| `MLP_v1`, `logReg_v1` | `pho_mlp_score_v1`, `pho_logreg_score_v1` |

Optional scores are photon-aligned and NaN outside their model's domain.
The recipe fixes the features, source roles, training weights and split.
AuAu variants share event assignments; the historical pp recipe keeps its
separate row split. Working points belong to their model; they are not
interchangeable between releases.

## Outputs and use in TreeToHists

Releases live in `models/<system>/<name>/` or the `--models` directory.
They contain the manifest, model artifacts, validation, working points and
`selection.json`. Keep that directory intact.

Each new release also saves `diagnostics.root` and `diagnostics.json` for
FinalAnalysis: raw and training-weighted ET/eta distributions, exact ROOT
feature moments and standardized mean differences, XGBoost total gain, and
weighted train/holdout AUC, log loss and Brier score. The ROOT file contains
ordinary histograms and scalar tables. It reuses validated model scores and
does not read or rewrite trees. Full distributions and metrics use the model's
own support. The ET/centrality QA grid extends to the wider model's support
and guard boundaries; these diagnostic edges do not prescribe analysis bins. Training weights are
separate from downstream generator and analysis weights.

Saved class score distributions and ROC tables support the pp/AuAu comparison
at 15–35 GeV, using AuAu centrality 50–80%. The PPG12 release also saves a paired
test-photon table with the original score, retrained score and their difference;
the frozen model is evaluated on the same H70 vector. Original source, entry
and photon IDs remain attached to those rows.

`trainCampaign` also writes an immutable `models.yaml` under
`<models>/.campaigns/<identity>/` and reports its path. It indexes the sealed
releases and QA hashes for FinalAnalysis's model selector and routing/depth
comparisons. Paths resolve relative to the YAML file; keep the models directory
layout intact when copying it. Routing/depth comparisons require matching
sources, prepared rows, event split and training weights. The wider model keeps
its own population and support.

Set `Pipeline.PhotonIDModels` to this `models.yaml` in the FinalAnalysis pipeline
configuration, or `QA.PhotonIDModels` for its standalone QA step. See
[FinalAnalysis's run guide](../FinalAnalysis/README.md) for the histogram inputs
and command. The model index supplies the nominal, range, routing and depth
choices; no training or tree projection is repeated by the plotting step.

BDT cells contain the saved XGBoost reference and a ROOT RBDT export. All
applicable predictions must agree within 1e-6 in absolute score. The ROOT
evaluator supplies release validation, working points and optional scores;
reloading it must pass the separate `saved_score_atol` gate. Each cell's JSON
records the ROOT file, object key and input count in `native_evaluator`.
Histogram readers can evaluate that same model on the H70 float32 feature
vector during their normal tree read. Reader integration and runtime
compatibility must be qualified before production use.

Select an optional model in a copy of TreeToHists' `treetohists.yaml`:

```yaml
auau:
  PhotonIDSelection: /path/to/models/AuAu/ETcentFine_v1/selection.json
```

Use `pp.PhotonIDSelection` for a retrained pp model. Omit score, state and
working-point overrides for that system; the release supplies them.
Scoring creates a checked group version and advances `CURRENT` only when
every source is ready. Partial work leaves the published group unchanged.

## Source map

The driver follows preparation → fitting → validation → working points →
scoring. Training labels and weights are defined in `Inputs.py`; histogram
selection belongs to TreeToHists.

| File | Responsibility |
|---|---|
| [photon_id.py](scripts/photon_id.py) | Commands and phase order |
| [Inputs.py](src/Inputs.py) | Features, labels, weights and event split |
| [ProjectionCache.py](src/ProjectionCache.py) | Compact source-row storage and verification |
| [Training.py](src/Training.py) | BDT, MLP and logistic fitting |
| [Validation.py](src/Validation.py) | Numerical and saved-model checks |
| [WorkingPoints.py](src/WorkingPoints.py) | Efficiency thresholds |
| [Diagnostics.py](src/Diagnostics.py) | Saved training QA for FinalAnalysis |
| [Scoring.py](src/Scoring.py) | Photon-aligned scores and group publication |
| [Models.py](src/Models.py) | Recipes, release names and artifact checks |

Unbound dataset locations are refused; set the published group paths in the
configuration or supply `--group`. Local model checks do not establish production or
downstream physics qualification.
