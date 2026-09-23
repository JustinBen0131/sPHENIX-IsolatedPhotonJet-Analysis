# Historical predecessors and first-pass limits

Evidence was read locally from /Users/patsfan753/Desktop/ThesisAnalysis on
2026-09-22. That repository was not modified. Paths below are relative to that
read-only library. HISTORICAL_SOURCES.json records the actual SHA256 of every
listed source/receipt. No historical runtime dependency or source tree was copied.

## Reused implementations

| Historical path | Reuse in this directory |
|---|---|
| scripts/ml/parity/train_corrected_ppg12_parity.py | Explicit corrected pp recipe and source substitution |
| ppg12codeGit/FunWithxgboost/data_loader.py | Per-source ET flattening, seeded sampling and row order in train.py |
| ppg12codeGit/FunWithxgboost/reweighting.py | Class/eta/ET spline formula; split into fit/apply for canonical TRAIN-only use |
| ppg12codeGit/FunWithxgboost/model_builder.py | Historical two-stage stratified split and XGBoost settings |
| ppg12codeGit/FunWithxgboost/variant_configs/config_base_v3E_split.yaml | Ordered BaseV3E features and recipe configuration |
| scripts/ml/parity/validation/compare_heldout_ppg12_models.py | Same-population weighted ROC, paired-score correlation and delta-score families |
| scripts/plotting/efficiency/derive_the327_auau_h70_wp_pt_centrality_grid.py | Verbatim weighted_threshold, weighted_efficiency, threshold_uncertainty, constant_fit, linear_fit in numerics.py |
| scripts/plotting/efficiency/make_the96_auau_isolation_efficiency_fits.py | Verbatim HistogramData and flow-aware quantile_from_histogram in numerics.py |
| scripts/plotting/efficiency/derive_schema14_auau_isolation_thresholds.py | Retained-isolation histogram/denominator reference; exact new truth mapping remains unbound |
| scripts/plotting/photon_performance/jstg_bdt_recap.py | Score/ROC diagnostic family |
| scripts/slides/auau_isolation/make_the327_auau_wp_fit_slides.py | Centrality/pT grid and fitted-curve diagnostic families |
| scripts/plotting/truth_purity/make_auau_bdt_score_separation_centrality_grid.py | Score-region family and strict boundary interpretation |

Only the small pure numerical functions were copied verbatim. Training weighting
and flattening were minimally adapted to explicit arrays and frozen fit records.
Rendering was adapted into one saved-payload route; no claim of pixel-identical
historical PNG reproduction is made. Slides 15–18 mix model/package generations;
the new package deliberately reports one model and one frozen validation package.

## AuAu H70 reference

Campaign:
the134_h70_factorial_14model_candidate_20260805_v16_auau_source_closure_system_fix.

Model: centAsFeatBase3x3_pt15to35.
Expected TMVA SHA256:
8328af75235c8bcbf63659ba791ceff4e1165018de269c889472c6a80a26a41f.

The model identity is recorded in:
dataOutput/plots/jstg_frozen_h70_cuts_20260909/nominal_correction_receipt.json:6-10.
Its model_metadata.json (nominal_source/ in the same directory) records the
campaign's TMVA/XGBoost artifact paths at lines 109–110 and full training metadata.

Recovered executable:
agent_context/local/the327_auau_bdt_retrain_20260915/sealed_the134_code_120d401d/training/train_auau_photon_bdt.py.

Actual trainer SHA256:
fd1da822af2778b2402a25d38d6eaf18fdacf85782799367c68bfc136c53b2c8.

This same hash is pinned in the neighboring stage2/science_admission_receipt.json:50
and remote/v16_runtime_readiness.json:100-101. These corroborate the sealed
predecessor used for the recovered recipe; this pass did not independently rerun
the historical artifact or claim bitwise model reproduction.

Relevant executable functions are ppg12_exact_inverse_pdf_weights (line 1598),
compute_ppg12_exact_global_weights (1661), event_level_train_test_split (2501),
and the model constructor near 2716. The new 80/10/10 split and validation-only
calibration intentionally replace the old 90/10 workflow. The expected historical
model hash is registered only as a reference, never assigned to a new model.

## pp parity reference

Original expected TMVA SHA256:
7679e634260402fb3815b2733767182690eec7587f9e09bffc307a05d00d59df.

Historical independently retrained XGBoost SHA256:
df2bbf795fd24f89340811911daf9ddf45e370dcd022b975a090c8d770a74fb2.

Recovered training receipt/config:
agent_context/local/the269_corrected_ppg12_parity_v1/final_local_chain_v1/postmerge_run_v4/training/TRAINING_RECEIPT.json
and DERIVED_CONFIG.yaml.

The actual parity trainer hash is:
c86d180409b3f1ef78d288280e4f404490c7548c85f239a513988aeead5b6ea9.

The new recipe preserves historical row splitting (seeds 42/43), full-population
weight fitting before the split, per-source low-ET flattening and source order,
and the 750-tree parameters. Jet8 substitutes for unavailable historical Jet5.
The loader labels original pid values {1,2} as signal; exact correspondence to
the new dominant-primary records is still a reviewed binding, not an assumption.

The executable eta weighting uses 20 bins on [-.7,.7], despite the old YAML
listing a wider eta range. The implementation is the evidence used here.
Parity ordering preserves the actual train_test_split returned training order,
not merely split membership. Source concatenation must be bound by the manifest.

The new model-specific WP70/80/90 calibration is separately versioned and does
not make its thresholds identical to original PPG12 formulas. The original
historical ID package remains unbound in the registry. The old retrained hash
is evidence, not an expected hash for a new training population.

## AuAu non-isolated boundary evidence

Executable:
agent_context/local/the249_chatgpt_pp_workbench/unified_response_currenttree_workbench_v1/eventtree_reduction/purity_factorial_candidate_packet_v4_cent0to80/auau_background_reader.py:270-271

It applies isolated = cone_sum < threshold and
nonisolated = cone_sum > sideband_threshold.

The completion-report generator:
agent_context/local/the249_chatgpt_pp_workbench/unified_response_currenttree_workbench_v1/eventtree_reduction/data_working_sample_46467_ppg12_dual_view_packet_v4/finalize_packet_v4.py:108

describes “AuAu method-2 centrality threshold 7.57-0.0658*c and zero gap.”
Thus zero gap has explicit historical evidence. This is not a fresh validation
of all accepted AuAu purity outputs, nor proof that a newly derived ISO90 curve
is the nominal alias. The new registry leaves that binding unresolved.
The historical numeric curve is not copied into the prospective model package.

## Deliberately open

1. Completed TreeProduction output manifest and exact source hashes/identities.
2. Feature/ratio/preselection equivalence and dominant-primary label mapping.
3. Cross-file hard/underlying-event identity witnesses for grouped splits.
4. Exact parity training domain and ordering of retained source populations.
5. AuAu truth-isolated analysis denominator binding; no invented isolation rule.
6. Human approval of the new nominal AuAu isolation target/sideband binding.
7. Model artifacts/import receipts and model-specific WP package bindings.
8. Local TMVA/PyROOT and real SDCC execution remain untested.
9. New sidecar integration into TreeToHists is intentionally deferred.

Items 1–5 and 7 are evidence/binding work, not reasons to alter scientific
definitions. Human scientific review is needed to accept new calibrations,
confirm the AuAu nominal alias/sideband, and approve the final field/label map.
No production or scientific closure is claimed.
