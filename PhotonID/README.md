# PhotonID

PhotonID trains the AuAu and pp photon-identification models used by the analysis.
It turns finalized simulation TTrees into **a fitted model, its ID cuts, and saved
validation/QA**. TreeToHists applies that release; FinalAnalysis draws the figures.

## The scientific flow

```text
TreeProduction: finalized simulation TTrees
    ↓
PhotonID: extract → prepare → fit → validate → calibrate
    ├─ model + ID cuts → TreeToHists → analysis histograms ─┐
    └─ saved training diagnostics ─────────────────────────┤
                                                          ↓
                                               FinalAnalysis: figures
```

Source extraction is shared across a campaign's variants. Preparation builds the
population, weights and split once per configured range; fitting reuses those
arrays. Preparation is database-free, with bounded passes and worker scratch for
temporaries. Finalized TG TTrees are valid inputs and may be read again when needed.

## What is fixed

The canonical recipe is defined in `config/training.yaml` and the named model recipe.

| Choice | Canonical definition |
|---|---|
| Inputs | 14 ordered AuAu features; 11 ordered pp features |
| Labels | Prompt photons from signal sources; non-prompt clusters from background sources, using the dominant-primary truth association. No isolation requirement |
| Fit | XGBoost, 450 trees, depth 4, seed 13; all 450 boosting rounds retained |
| Split | Whole events assigned to 90% training / 10% holdout using stable source/event identity |
| Weights | Class balance and 20-bin ET/eta spline weights over the full selected range, before routing |
| ID regions | Tight: `score > WP70`; non-tight: `0.1 < score < WP80`; WP90 is diagnostic |

WP70/80/90 come from weighted signal-efficiency quantiles and the configured
threshold fit. Each release carries its own working points and applicability.

## What varies

**One AuAu baseline** anchors the comparisons. A routing or depth scan changes only
its named knob, reusing the nominal features, weights and event assignments.

| AuAu recipe | Changed choice | Train / calibrate / apply | Result range; guard support |
|---|---|---|---|
| `nominal_retrained_v1` | Baseline | 15–35 GeV | 17–33; 15–17 / 33–35 GeV |
| `wideET_v1` | ET support | 12–40 GeV | 14–38; 12–14 / 38–40 GeV |
| `centCoarse_v1`, `centFine_v1` | Centrality routing: C3 or C7 | 15–35 GeV | As baseline |
| `ETcentCoarse_v1`, `ETcentFine_v1` | ET × C3 or ET × C7 routing | 15–35 GeV | As baseline |
| `depth2_v1` … `depth10_v1` | Maximum tree depth 2–10; depth 4 reuses baseline | 15–35 GeV | As baseline |

C3/C7 mean three/seven centrality intervals; ET routing uses the configured bins.
These ranges include the lower endpoint and exclude the upper. The wide recipe
retains nominal assignments for common events and has its own full-range weights
and calibration. Guard support does not prescribe the interior analysis binning.

Canonical pp uses the 11-feature, 15–35 GeV recipe. The separate
`ppg12_retrained_v1` comparison reproduces the configured historical mechanics:
6–40 GeV closed acceptance, 750 trees, depth 5, and stratified 70/10/20 rows with
seeds 42/43. Its validation and test rows stay separate; it defines no result/guard bins.

## What a release contains

Releases live under `<models>/<system>/<name>/`. Keep each directory intact.

| Product | What it makes inspectable or usable |
|---|---|
| Model files + `manifest.json` | Fitted cells and hashes binding sources, recipe, code, runtime and artifacts |
| `working_points.json` + `selection.json` | Model-specific thresholds, ID regions and application range for TreeToHists |
| `validation.json`, `diagnostics.root`, `diagnostics.json` | Saved score checks, train/holdout metrics, feature/weight distributions and calibration quality |
| `models.yaml` (campaign collection) | Verified releases and saved QA for FinalAnalysis comparisons |

Stable event identities, fixed seeds/splits, saved weights and artifact hashes make
results reproducible and comparisons traceable. Saved ROOT predictions must agree
with XGBoost within `1e-6`; reloading the saved evaluator has its own numerical gate.
Validation and calibration are mandatory. Invalid scores/artifacts or insufficient
coverage stop release. Measured performance or threshold-line quality deviations
are retained as `REVIEW_REQUIRED` for review; physics acceptance remains separate.

## Using PhotonID

The single entrance is `./PhotonID/photon_id.sh`, run from the repository root.
Copy `config/training.yaml`, set the completed simulation groups and output directory,
then follow these stages in a compatible Python/ROOT environment:

| Stage | Command | Result |
|---|---|---|
| Select inputs | `inspectCampaign` | Exact source catalogue |
| Extract | `projectCampaign`, then `verifyProjections` | Shared, verified training inputs |
| Prepare | `prepareCampaign` | Weights, event split and arrays for each configured range |
| Train | `trainCampaign` | Fit, validate, calibrate and publish each selected model |
| Collect QA | `collectQA` | Saved model collection for FinalAnalysis |

Choose `AuAu` or `PP` and a named recipe from the tables above.
`planCampaign AuAu --summary` previews the full variant matrix; `--check` previews
extraction or training before execution. Matching completed work is verified and reused.
Training leaves the input trees unchanged.

[Follow the run guide](docs/run-guide.md#run) for the exact commands, shard execution
and worker scratch setup.

Before the usual histogram pass, add all released scores together:

```bash
./PhotonID/photon_id.sh scoreCampaign AuAu --group /path/to/group
./PhotonID/photon_id.sh scoreCampaign PP --group /path/to/group
```

Each source is copied once for all missing scores. Earlier scores, fractions
and original fields are preserved; matching completed groups are reused.
For parallel production, `--candidate-root /path/to/work` stages the checked
copies and receipts without changing the published group. Its publisher commits
the complete group against the same input version.
TreeToHists selects one `Model` and its working points in its YAML. The default
is the new AuAu nominal model and the original PPG12 model for pp.

For saved analysis QA, give FinalAnalysis the collected `models.yaml`. Drawing saved model diagnostics does not refit.
Working-point/score histogram QA requires applying the model through TreeToHists
and rendering those histograms through FinalAnalysis.

[Full run guide: setup, campaign commands, model details and downstream use](docs/run-guide.md).
