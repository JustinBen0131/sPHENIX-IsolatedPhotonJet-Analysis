# PhotonID first-pass implementation report

2026-09-22. Local implementation only; no production, external publication or scientific acceptance.

## A. Existing code

Originally five files: README.md, features.py, train.py, augment.py,
model_registry.yaml. No PhotonID-local test suite or public shell dispatcher.
The existing feature adapter and evaluator were retained; training and sidecar
orchestration were refined in place. The old training label depended on the
isolation-bearing analysis_signal flag; it is now an independent named adapter.

## B. Historical reuse

[PROVENANCE.md](PROVENANCE.md) records exact historical paths, functional reuse,
model identities and remaining limits. [HISTORICAL_SOURCES.json](HISTORICAL_SOURCES.json)
contains freshly computed source/receipt SHA256 values. Pure weighted quantile,
uncertainty, constant/linear fit and isolation-histogram functions were copied;
weighting, flattening and split mechanics were adapted to the shared pipeline.
No runtime imports from the historical repository exist.

## C. Files changed

Modified: README.md, features.py, train.py, augment.py, model_registry.yaml.

Created: train_photon_id.sh, augment_photon_id.sh, frontend.py, registry.py,
inputs.py, calibration.py, numerics.py, diagnostics.py, training_profiles.yaml,
requirements.txt, .gitignore, tests/test_photon_id.py, CONTRACTS.md,
PROVENANCE.md, HISTORICAL_SOURCES.json, FIRST_PASS_REPORT.md.

No TreeProduction, TreeToHists or FinalAnalysis file was changed by this pass.
Pre-existing TreeProduction changes remain untouched.

## D. Actual frontend help

### ./train_photon_id.sh --help

```text
usage: ./train_photon_id.sh [-h] [--dry-run] [--registry REGISTRY]
                            [--output-root OUTPUT_ROOT] [--campaign CAMPAIGN]
                            [--campaign-sha256 CAMPAIGN_SHA256]
                            [--profiles PROFILES]
                            [COMMAND]

sPHENIX Photon ID — candidate model training

Commands:
  trainAuAu              canonical_auau_v1
  trainPP                canonical_pp_v1
  trainPPG12Equivalent   ppg12_equivalent_v1 (historical mechanics)
  trainAll               the same three primitive paths, sequentially

Missing population/physics bindings fail closed. No automatic acceptance.

positional arguments:
  COMMAND

options:
  -h, --help            show this help message and exit
  --dry-run             validate bindings and show plan; do not train/write
  --registry REGISTRY   registry YAML
  --output-root OUTPUT_ROOT
                        override output/photon_id root
  --campaign CAMPAIGN   development override: completed TreeProduction
                        manifest
  --campaign-sha256 CAMPAIGN_SHA256
                        required hash with --campaign override
  --profiles PROFILES   training profile YAML
```

### ./augment_photon_id.sh --help

```text
usage: ./augment_photon_id.sh [-h] [--dry-run] [--registry REGISTRY]
                              [--output-root OUTPUT_ROOT]
                              [--campaign CAMPAIGN]
                              [--campaign-sha256 CAMPAIGN_SHA256]
                              [--source SOURCE] [--models MODELS]
                              [--selections SELECTIONS]
                              [COMMAND]

sPHENIX Photon ID — immutable-base ROOT sidecars

Data:
  runPP  runAuAu  runAllData
Simulation:
  runPhotonJetSim  runInclusiveJetSim
  runEmbeddedPhotonJetSim  runEmbeddedInclusiveJetSim  runAllSim
Everything:
  runAll

Default: nominal_pp or nominal_auau selection for the lane's system.
--models requests raw scores without requiring a selection recipe.

positional arguments:
  COMMAND

options:
  -h, --help            show this help message and exit
  --dry-run             validate bindings and show plan; do not train/write
  --registry REGISTRY   registry YAML
  --output-root OUTPUT_ROOT
                        override output/photon_id root
  --campaign CAMPAIGN   development override: completed TreeProduction
                        manifest
  --campaign-sha256 CAMPAIGN_SHA256
                        required hash with --campaign override
  --source SOURCE       one-based source within a primitive lane; same
                        execution path
  --models MODELS       comma-separated model names, or all-compatible; raw
                        scoring
  --selections, --selection SELECTIONS
                        comma-separated recipe names, or all-compatible
```

## E–F. Profiles and selection recipes

The complete profile table, weighting/hyperparameter details, and nominal recipe
table are in [README.md](README.md). There are three trainable candidate identities:
14-feature canonical AuAu and 13-feature canonical pp use physical-event-grouped
80/10/10; 11-feature PPG12-equivalent retains row-stratified 70/10/20 and historical
weighting/flattening. Reference model artifacts remain unbound.

nominal_pp names ppg12_equivalent_v1 and its own model-specific ID package, with
the versioned PPG12 topocluster R=.4 isolation formula. nominal_auau names
canonical_auau_v1 and its own ID/isolation packages. Its nominal isolation alias
and sideband remain explicit unresolved bindings. Historical zero-gap evidence
is recorded; a new ISO90 identity was not inferred from it.

## G–H. Exact schema and output structure

[CONTRACTS.md](CONTRACTS.md) lists every ROOT field and type, enums, key semantics,
manifest interface and registry bindings. Tables: Models, PhotonScores,
SelectionRecipes, PhotonSelections. All candidate relations carry source/event/
photon hi/lo words. No ABCD field or per-model base-tree branch is introduced.

[README.md](README.md) gives the full deterministic training/augmentation output
tree, diagnostic filenames, saved numerical payloads, freeze and receipt files.

## I. Runnable versus fail-closed

Help, configuration parsing, numerical calibration/region functions, saved-payload
rendering and synthetic training/sidecar paths are locally exercised. Actual
population training and augmentation remain intentionally blocked by missing
reviewed bindings. Raw scoring does not require a selection recipe. Missing AuAu
analysis-isolation denominator disables only its calibration/plots; missing
nominal sideband does not block an otherwise complete model/ID package.

Actual default-command observations (no output generated):

- trainAuAu: exit 2; ERROR: canonical_auau_v1: unresolved feature_binding_evidence in training_profiles.yaml
- trainPP: exit 2; ERROR: canonical_pp_v1: unresolved feature_binding_evidence in training_profiles.yaml
- trainPPG12Equivalent: exit 2; ERROR: ppg12_equivalent_v1: unresolved feature_binding_evidence in training_profiles.yaml
- trainAll: exit 2; ERROR: canonical_auau_v1: unresolved feature_binding_evidence in training_profiles.yaml
- runPP: exit 2; ERROR: TreeProduction completed campaign manifest is unbound; configure input_campaign and input_campaign_sha256
- runAuAu: exit 2; ERROR: TreeProduction completed campaign manifest is unbound; configure input_campaign and input_campaign_sha256
- runAll: exit 2; ERROR: TreeProduction completed campaign manifest is unbound; configure input_campaign and input_campaign_sha256

## J. Validation

31 tests passed in 4.400 seconds before semantics-preserving Black formatting.
Black validated AST equivalence. Final Python compile, bash syntax and git diff
whitespace checks passed after formatting. Tests cover:

- grouped splits, shared embedding components and partition disjointness;
- exact historical row-split indices, per-source thinning order and weighting;
- train-only weight fits and validation-only calibration;
- strict ID/isolation boundaries, ties, curve crossing and histogram overflow;
- shared feature joins and uint64 identities above 2^63;
- independent prompt labels versus analysis-isolation labels;
- incompatible system, model/WP/receipt mismatch and incomplete bindings;
- tiny real XGBoost fit/serialization/reload;
- a 3,000-row synthetic training package, freeze before test prediction, saved
  diagnostics and file hashes (25 trees solely for the fixture);
- multi-model keyed ROOT sidecar round-trip, raw-score-only sidecar,
  invalid states, source consistency, overwrite refusal and unchanged base hash;
- aggregate expansion, no-argument/help behavior from another directory,
  unknown commands and fail-closed defaults.

No production-sized training, SDCC, Fun4All, PyROOT/TMVA or downstream runtime
compatibility was tested. Large-population memory use remains unqualified.

Local dependency versions: numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, scikit-learn 1.9.1, xgboost 3.4.1, uproot 5.7.6, matplotlib 3.11.2, PyYAML 6.0.3.

## K. Remaining human scientific decisions

Accept/revise the prospective model-specific WP fits after real validation;
confirm the new nominal AuAu isolation efficiency and its supported zero-gap
sideband binding; approve the final historical-to-canonical feature/label and
analysis-denominator equivalence. Missing manifests/artifact paths are technical
bindings, not reasons to invent a physics choice. None blocks this completed
first-pass infrastructure deliverable.
