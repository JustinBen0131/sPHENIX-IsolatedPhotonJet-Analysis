# Figures

Every file here is written by `plots/IANFigures.C` (the measurement, from `results.root` of `final_analysis.sh run`,
with `systematics.root` when given) and `plots/QAFigures.C` (photon and event QA, from `qa.root` of
`final_analysis.sh qa`), into the same folders and file names as the IAN's `fig/` folder, so the folder can be copied
into the manuscript as it is. `FIGURES.txt` lists what one `draw` wrote and every figure it could not draw, with the
reason.

Every figure follows the IAN's layout: the sPHENIX mark and its physics lines (collision system and energy, photon
window, recoil-jet definition, all compiled from the saved metadata) at the top left inside the frame, the legend at the
top right, a y range opened until no point reaches either block, and the IAN's palette and markers. Maps add a
data-free band above the data for the header instead of covering cells.

Every xJ figure shows only the reported range, where the distribution is fully populated: from the first xJ edge at or
above the recoil-jet threshold over the lowest photon pT (derived by the analysis, so it follows any rebinning) to
`Measurement.ReportedXjMax`. The axis always starts at 0.

The iteration-stability figures draw the frozen iteration rule's saved table (data and response statistics, the change
from k - 1 and their quadrature Q_k, with the selected k); a result unfolded with configured iterations has no table and
lists them as not drawn.

The tables below map generated output filenames to the measurement and QA figures.

| File | IAN figure | Note |
|---|---|---|
| `purityAndCombSub/pp_xj_purity_correction.png` | fig:xjg_purity_correction |  |
| `purityAndCombSub/pp_xj_combinatoric_crosscheck.png` | fig:xjg_combinatorial_subtraction |  |
| `Unfolding/pp_photon_response.png` | fig:photon_response |  |
| `Unfolding/pp_joint_response.png` | fig:joint_response |  |
| `Unfolding/pp_photon_yield_iteration_stability_publication.png` | fig:photon_iteration_stability |  |
| `Unfolding/pp_joint_iteration_stability.png` | fig:joint_iteration_stability |  |
| `Unfolding/pp_photon_closure_overlay.png` | fig:photon_closure |  |
| `Unfolding/pp_photon_yield_unfolding_nominal_2bin_k2.png` | fig:photon_yield_unfolding | IAN slot pp_photon_yield_unfolding_nominal_2bin_k3.png: the file name carries the iteration count |
| `Unfolding/pp_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `Unfolding/pp_xj_before_after_unfolding.png` | fig:xjg_before_after_unfolding |  |
| `purityAndCombSub/auau020_xj_purity_correction.png` | fig:xjg_purity_correction |  |
| `purityAndCombSub/auau020_xj_combinatoric_subtraction.png` | fig:xjg_combinatorial_subtraction |  |
| `Unfolding/auau_0to20_photon_response.png` | fig:photon_response |  |
| `Unfolding/auau_0to20_joint_response.png` | fig:joint_response |  |
| `Unfolding/auau020_photon_yield_iteration_stability.png` | fig:photon_iteration_stability |  |
| `Unfolding/auau020_joint_iteration_stability.png` | fig:joint_iteration_stability |  |
| `Unfolding/auau020_photon_closure_overlay.png` | fig:photon_closure |  |
| `Unfolding/auau020_photon_yield_unfolding_nominal_2bin_k2.png` | fig:photon_yield_unfolding |  |
| `Unfolding/auau020_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `Unfolding/auau020_xj_before_after_unfolding.png` | fig:xjg_before_after_unfolding |  |
| `purityAndCombSub/auau2050_xj_purity_correction.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau2050_xj_combinatoric_subtraction.png` | new: no IAN slot yet |  |
| `Unfolding/auau_20to50_photon_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau_20to50_joint_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_photon_yield_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_joint_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_photon_closure_overlay.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_photon_yield_unfolding_nominal_2bin_k2.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_xj_before_after_unfolding.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau5080_xj_purity_correction.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau5080_xj_combinatoric_subtraction.png` | new: no IAN slot yet |  |
| `Unfolding/auau_50to80_photon_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau_50to80_joint_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_photon_yield_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_joint_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_photon_closure_overlay.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_photon_yield_unfolding_nominal_2bin_k2.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_xj_before_after_unfolding.png` | new: no IAN slot yet |  |
| `Unfolding/pp_vs_auau_xj_unfolded.png` | fig:final_xjg_pp_auau | p+p and the most central Au+Au class |
| `FinalResults/pp_vs_auau_centrality_xj_unfolded.png` | new: no IAN slot yet | every class |
| `LHCComparison/pp_xj_atlas_overlay.png` | new: no IAN slot yet | ATLAS p+p 5.02 TeV, PLB 789 (2019) 167 |
| `LHCComparison/auau020_xj_atlas_overlay.png` | new: no IAN slot yet | ATLAS 0-10% Pb+Pb 5.02 TeV, PLB 789 (2019) 167 |

## Photon and event QA

Drawn by `plots/QAFigures.C` from `qa.root` once packages with QA fills exist (the new production; the legacy July
histograms carry none). Every IAN QA slot has its own file name below; the shower-shape figures are also written for the
other Au+Au classes (`cent1`, `cent2`). The photon-ID training figures for IAN
40–44 are drawn from an AuAu collection supplied through `PhotonIDModels`:
ET/eta reweighting, feature separation, split gain, routing and depth scans.
The pp collection draws figures 45 and 46 from its saved held-out score/ROC
records and paired original-versus-retrained PPG12 test scores. Figure 45 uses
15–35 GeV and AuAu 50–80%, with a separately sealed canonical AuAu reference;
figure 46 uses the PPG12 model's test rows and excludes its validation rows.
Missing variants or undefined statistics are recorded in `FIGURES.txt`.
The cluster-width panels compare DATA with simulation after the same tight and
isolated photon selection; the package retains the fine width distributions
for each centrality slice, alongside the preselection distributions.
BDT working-point and per-centrality score figures (13 and 14) use model-selected
histogram QA from TreeToHists. Their source is separate from training diagnostics.

| Model QA file under `photonOnly/` | Input |
|---|---|
| `auau_training_et_reweighting.png`, `auau_training_eta_reweighting.png` | Saved AuAu training distributions |
| `auau_bdt_single_feature_separation.png`, `auau_bdt_feature_split_gain.png` | Saved feature moments and split gain |
| `auau_bdt_routing_auc_by_et_centrality.png` | Matched routing releases |
| `auau_bdt_tree_depth_train_holdout_diagnostic.png` | Matched depth releases |
| `pp_auau_bdt_heldout_score_separation.png`, `pp_auau_bdt_heldout_roc.png` | Canonical pp and peripheral AuAu held-out records |
| `pp_bdt_original_vs_retrained_scores.png`, `pp_bdt_retraining_score_difference.png` | Paired PPG12 test-photon scores |

Alternative models' weight and feature figures are under `models/<name>/`.
The nominal model uses the stable top-level names. PNG sidecars bind the model
collection and imported references; drawing a figure does not establish
scientific acceptance.

| File | Kind |
|---|---|
| `event_obs/hMaxClusterEnergy_groupOverlay.png` | event and sample QA |
| `event_obs/hMaxClusterEnergy_groupTurnOnOverlay.png` | event and sample QA |
| `event_obs/pp_max_cluster_energy_overlay.png` | event and sample QA |
| `event_obs/pp_max_cluster_energy_turnon.png` | event and sample QA |
| `simulation_stitching/pp_photonjet_generator_stitching.png` | event and sample QA |
| `simulation_stitching/pp_inclusivejet_generator_stitching.png` | event and sample QA |
| `simulation_stitching/auau_photonjet_generator_stitching.png` | event and sample QA |
| `simulation_stitching/auau_inclusivejet_generator_stitching.png` | event and sample QA |
| `simulation_stitching/auau_centrality_data_embeddedMC.png` | event and sample QA |
| `simulation_stitching/pp_photonjet_reco_tag_stitching.png` | event and sample QA |
| `simulation_stitching/pp_inclusivejet_reco_tag_stitching.png` | event and sample QA |
| `simulation_stitching/auau_photonjet_reco_tag_stitching.png` | event and sample QA |
| `simulation_stitching/auau_inclusivejet_reco_tag_stitching.png` | event and sample QA |
| `event_obs/centrality.pdf` | event and sample QA |
| `event_obs/vertexz.pdf` | event and sample QA |
| `event_obs/cluster_pt.pdf` | event and sample QA |
| `event_obs/cluster_eta_phi_data.pdf` | event and sample QA |
| `event_obs/cluster_eta_phi_photon.pdf` | event and sample QA |
| `event_obs/jet_pt.pdf` | event and sample QA |
| `event_obs/jet_eta_phi_data.pdf` | event and sample QA |
| `event_obs/jet_eta_phi_photon.pdf` | event and sample QA |
| `photonOnly/auau_minbias_pretraining_before_after.png` | photon QA |
| `photonOnly/pp_ncb_score_data_mc.png` | photon QA |
| `photonOnly/auau_0to80_ncb_score_data_mc.png` | photon QA |
| `photonOnly/pp_ncb_eta_time_before.png` | photon QA |
| `shower_shapes/dis_e11_to_e33_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_e32_to_e35_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_bdt_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_weta_cogx_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_wphi_cogx_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_et1_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_et2_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_et3_cent0_pt0_cut1.pdf` | photon QA |
| `shower_shapes/dis_et4_cent0_pt0_cut1.pdf` | photon QA |
| `photonOnly/auau_signal_full_cluster_width.png` | photon QA |
| `photonOnly/auau_signal_3x3_width.png` | photon QA |
| `photonOnly/auau_background_full_cluster_width.png` | photon QA |
| `photonOnly/auau_background_3x3_width.png` | photon QA |
| `photonOnly/auau_h70_bdt_pt_centrality_grid.png` | photon QA |
| `photonOnly/auau_h70_bdt_centrality_fits.png` | photon QA |
| `photonOnly/auau_bdt_score_perCent.png` | photon QA |
| `shower_shapes/e11_e33_selection_flow.png` | photon QA |
| `photonOnly/pp_auau_truth_isolation_direct_fragmentation_2x2.png` | photon QA |
| `photonOnly/pp_truth_photon_pt_after_isolation.png` | photon QA |
| `photonOnly/auau_50to80_truth_photon_pt_after_isolation.png` | photon QA |
| `photonOnly/auau_20to50_truth_photon_pt_after_isolation.png` | photon QA |
| `photonOnly/auau_0to20_truth_photon_pt_after_isolation.png` | photon QA |
| `isolation/isoCentOverlay_noUE_pT_14_16.png` | photon QA |
| `isolation/isoCentOverlay_withUE_pT_14_16.png` | photon QA |
| `isolation/auau_r03_isolation_pt_centrality_grid.png` | photon QA |
| `isolation/auau_r03_isolation_centrality_fits.png` | photon QA |
| `isolation/pp_nominal_topo_isoStack.png` | photon QA |
| `isolation/auau_50to80_isoStack.png` | photon QA |
| `isolation/auau_20to50_isoStack.png` | photon QA |
| `isolation/auau_0to20_isoStack.png` | photon QA |
| `photonOnly/pp_auau_abcd_region_counts_1x4.png` | photon QA |
| `photonOnly/data_purity_and_prompt_leakage_2x4_style.png` | photon QA |
| `photonOnly/purity_qa_pp_auau_simCLOSURE.png` | photon QA |
| `photonOnly/pp_si_di_auau_photon_efficiency_1x4.png` | photon QA |
