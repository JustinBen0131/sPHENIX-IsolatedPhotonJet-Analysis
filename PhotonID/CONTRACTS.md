# PhotonID interfaces

These are local implementation contracts, not claims that a production campaign
already implements them. TreeProduction scientific behavior is unchanged.

## Completed campaign input

Registry/training profile fields input_campaign and input_campaign_sha256 bind one
JSON/YAML PhotonJetTreeCampaignV1 manifest. It must have:

- schema: PhotonJetTreeCampaignV1
- status: complete
- campaign: a portable identifier
- sources: nonempty ordered list of completed base ROOT files
- each source: lane, path, sha256, source_id (two exact uint64 words)
- for training sources: sample, training_witness (path and sha256).

Lane names are the six names in README. sample names are the explicit training
recipe's Photon/Jet categories, not directory-name guesses. Paths resolve from
the manifest directory. Source identity and ROOT completion/accounting must
match. Duplicate files or source identities fail. Each input base file currently
contains exactly one completed Source. The manifest must represent the complete
lane; directory globbing and silent skipping are not supported.

This is a completed-output interface. Existing TreeProduction DST manifests
cannot be substituted. Emitting and validating this manifest is a later runtime
binding task; no real population or generator has been invented in this pass.

### Training witness boundary

The base currently does not certify all historical training/preselection and
cross-production physical-event mappings. Training therefore additionally needs
one hash-bound witness document per source:

- base_sha256: exact base-file hash.
- label_adapter: dominant_prompt_v1.
- events: records with event_id [event_hi,event_lo], physical_components
  (nonempty globally meaningful strings), and for embedded AuAu
  embedding_components_complete: true.
- photons: records with photon_id [event_hi,event_lo,photon_hi,photon_lo],
  training_eligible and isolation_eligible as explicit booleans.

Source identity is bound by the manifest/base hash; it is not omitted from final
candidate keys. Every candidate must have a policy witness; canonical profiles
also require its event witness. Parity uses the historical row split and does
not require physical grouping witnesses.

physical_components must name both hard-scatter and underlying embedded event
identities, including repeated use across files. Connected components determine
canonical group membership. A hash of the output row or arbitrary unique string
is not scientific evidence of distinct physical events. Changing the population
can change connected components and therefore constitutes a new split.

This is an explicit provisional adapter boundary, **not a new repair production
stage**. Nothing fabricates witnesses or labels. Bind reviewed predicates and
identities from the final source authority before running real training; replace
this boundary with direct retained-record bindings when established. The profile
separately requires feature_binding_evidence, label_binding_evidence,
preselection_binding_evidence, and
domain_binding_evidence. A truth-isolated analysis denominator additionally needs
analysis_signal_binding_evidence; missing it disables only isolation calibration.

The BDT adapter reads dominant_truth_state/pid/track_id/embedding_id from Photons,
joins TruthPhotons by event + native track + embedding, and requires a complete
truth census. Valid photon direct/fragmentation classes (1/2) are signal, known
hadronic class 3 is nonprompt. A measured NoPrimary or valid non-photon primary
is nonprompt. Unknown/unclassified/missing associations stay unknown. Signal
sources keep prompt labels, background sources keep nonprompt labels. No
reconstructed-isolation condition enters this adapter.

This mapping is implemented and unit-tested, but its historical equivalence is
unbound. In particular, it does not pretend that historical class<3, classes{1,2},
and the current analysis_signal flag are interchangeable.

## ROOT sidecar: PhotonIDSidecarV2

All four products are explicit ROOT TTrees, including zero-row tables. Metadata
is a JSON TObjString. There is one PhotonScores row per base photon × requested
model, and one PhotonSelections row per photon × requested recipe. No base
candidate is dropped due to invalid inputs. No join relies on row order.

Types below are uproot branch declarations: uint64 is ULong64_t, int32 is Int_t,
float64 is Double_t; string is a scalar variable-length ROOT string branch.

### Shared candidate key

Both candidate tables carry six uint64 fields, copied exactly from TreeProduction:

source_hi, source_lo, event_hi, event_lo, photon_hi, photon_lo.

### Models

| Field | Type |
|---|---|
| model_id | uint64 |
| model_name, model_version | string |
| collision_system | int32 (1 pp, 2 AuAu) |
| model_sha256, feature_schema_sha256, training_receipt_sha256 | string |

model_id is a deterministic hash-derived uint64 of name, artifact hash and
feature identity. Collisions/duplicates within a product fail. Cross-product
consumers must verify its metadata binding, not assume uint64 equality alone
proves model equivalence.

### PhotonScores

| Field | Type |
|---|---|
| six shared key fields | uint64 |
| model_id | uint64 |
| score | float64 |
| score_valid, state, complete, in_domain, shower_valid | int32 |

score_valid is 1 only for finite evaluated scores. States: 1 finite,
2 evaluated_nonfinite, 3 missing_inputs, 4 out_of_domain, 5 model_unavailable.
Unbound models fail preflight, so state 5 exists for internal diagnostics and
compatibility but is not emitted by a successful production command.
Invalid scores are NaN. Applicability and completeness are retained separately.

### SelectionRecipes

| Field | Type |
|---|---|
| selection_id, model_id | uint64 |
| selection_name | string |
| id_package_sha256, isolation_package_sha256 | string |

selection_id binds name, model_id and both package hashes. Model-specific ID
packages must match both model SHA and feature SHA. External packages use the
actual file SHA; inline isolation configuration uses canonical JSON SHA.

### PhotonSelections

| Field | Type |
|---|---|
| six shared key fields | uint64 |
| selection_id, model_id | uint64 |
| bdt_score | float64 |
| tight_threshold, nontight_upper_threshold, nontight_lower_threshold | float64 |
| id_region | int32 |
| isolation_value, isolation_radius | float64 |
| isolation_method | int32 (2 SUB1, 3 topocluster) |
| isolated_threshold, nonisolated_threshold | float64 |
| isolation_region, valid | int32 |

ID codes: INVALID=0, BELOW_FLOOR=1, NONTIGHT=2, EXCLUDED=3, TIGHT=4.
Isolation codes: INVALID=0, ISOLATED=1, GAP=2, NONISOLATED=3.
valid=1 means both classifications are valid, including an excluded/gap result.
It does not mean selected for physics. Strict inequalities put equal isolation
boundaries in GAP, even for zero gap. An invalid score forces id_region=INVALID;
valid isolation can remain independently recorded.

Metadata binds base SHA/source identity, campaign manifest SHA, all model
identities, full recipe/calibration payloads, state counts and enum meanings.
No model-dependent branches or ABCD field are added to the base.

Sidecar publication is exclusive and atomic on a local filesystem. Before
publication, uproot readback checks table schemas/counts and every identity word.
The base hash is checked before/after writing. Existing outputs cause failure.

## Packages and registry bindings

A model takes ordered features to scores. ID and isolation packages each map
stored values to regions. A selection recipe names a combination; a nominal
alias change is a registry edit, not distributed conditional code.

New ID packages store system, axis (et/centrality), half-open domain, linear
curves WP70/WP80/WP90, non_tight_lower, validation fit population, model SHA,
feature SHA, tables and candidate status. Ordering is checked over the complete
linear domain. Historical formulas are registered separately and remain unbound;
the new retrained model never silently borrows original PPG12 thresholds.

Isolation packages store system, method, radius, axis, domain, isolated and
nonisolated linear coefficients. AuAu calibration additionally preserves
ISO70/ISO80/ISO90 and fit tables, independently of the nominal alias.

An imported model requires model_file, model_sha256, feature_schema_sha256,
training_manifest, training_receipt_sha256, candidate/approved status and binding
evidence. The receipt's model/feature hashes must match. Frozen references also
must match expected_artifact_sha256. No importer downloads from personal paths.

## Downstream boundary

TreeToHists must join sidecars on all six candidate-key words plus model_id or
selection_id, reject base/model/package hash mismatches, and select the requested
recipe. It will combine TIGHT/NONTIGHT with ISOLATED/NONISOLATED into ABCD.
Its existing single-model reader has not been migrated in this task.

Purity, leakage, final recoil cuts, normalization, responses, unfolding and plots
of final measurements remain downstream. Model diagnostics do not perform those
analysis operations. Score-domain extension for response boundary populations
remains a later explicit contract decision, not silent extrapolation.
