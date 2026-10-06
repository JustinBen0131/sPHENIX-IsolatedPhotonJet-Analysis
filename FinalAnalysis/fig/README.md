# Figures

Every file here is written by `plots/IANFigures.C` (the measurement, from `results.root` of `final_analysis.sh run`,
with `systematics.root` when given) and `plots/QAFigures.C` (photon and event QA, from `qa.root` of
`final_analysis.sh qa`), into the same folders and file names as the IAN's `fig/` folder, so the folder can be copied
into the manuscript as it is. `FIGURES.txt` lists what one `draw` wrote and every figure it could not draw, with the
reason.

Every xJ figure shows only the reported range, where the distribution is fully populated: from the first xJ edge at or
above the recoil-jet threshold over the lowest photon pT (derived by the analysis, so it follows any rebinning) to
`Measurement.ReportedXjMax`. The axis always starts at 0.

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
| `FinalResults/pp_xj_before_after_unfolding.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau020_xj_purity_correction.png` | fig:xjg_purity_correction |  |
| `purityAndCombSub/auau020_xj_combinatoric_subtraction.png` | fig:xjg_combinatorial_subtraction |  |
| `Unfolding/auau_0to20_photon_response.png` | fig:photon_response |  |
| `Unfolding/auau_0to20_joint_response.png` | fig:joint_response |  |
| `Unfolding/auau020_photon_yield_iteration_stability.png` | fig:photon_iteration_stability |  |
| `Unfolding/auau020_joint_iteration_stability.png` | fig:joint_iteration_stability |  |
| `Unfolding/auau020_photon_closure_overlay.png` | fig:photon_closure |  |
| `Unfolding/auau020_photon_yield_unfolding_nominal_2bin_k2.png` | fig:photon_yield_unfolding |  |
| `Unfolding/auau020_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `FinalResults/auau020_xj_before_after_unfolding.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau2050_xj_purity_correction.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau2050_xj_combinatoric_subtraction.png` | new: no IAN slot yet |  |
| `Unfolding/auau_20to50_photon_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau_20to50_joint_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_photon_yield_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_joint_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_photon_closure_overlay.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_photon_yield_unfolding_nominal_2bin_k2.png` | new: no IAN slot yet |  |
| `Unfolding/auau2050_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `FinalResults/auau2050_xj_before_after_unfolding.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau5080_xj_purity_correction.png` | new: no IAN slot yet |  |
| `purityAndCombSub/auau5080_xj_combinatoric_subtraction.png` | new: no IAN slot yet |  |
| `Unfolding/auau_50to80_photon_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau_50to80_joint_response.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_photon_yield_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_joint_iteration_stability.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_photon_closure_overlay.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_photon_yield_unfolding_nominal_2bin_k2.png` | new: no IAN slot yet |  |
| `Unfolding/auau5080_xj_prior_reweighting.png` | new: no IAN slot yet |  |
| `FinalResults/auau5080_xj_before_after_unfolding.png` | new: no IAN slot yet |  |
| `FinalResults/pp_vs_auau_centrality_xj_unfolded.png` | new: no IAN slot yet |  |
| `LHCComparison/pp_xj_atlas_overlay.png` | new: no IAN slot yet | ATLAS p+p 5.02 TeV, PLB 789 (2019) 167 |
| `LHCComparison/auau020_xj_atlas_overlay.png` | new: no IAN slot yet | ATLAS 0-10% Pb+Pb 5.02 TeV, PLB 789 (2019) 167 |

## Photon and event QA

Drawn by `plots/QAFigures.C` from `qa.root` once packages with QA fills exist (the new production; the legacy July
histograms carry none). Every IAN QA slot has its own file name below; the shower-shape figures are also written for the
other Au+Au classes (`cent1`, `cent2`). The frozen photon-ID model's training and validation figures
(`photonOnly/auau_training_*`, `auau_bdt_*`, `pp_auau_bdt_heldout_*`, `pp_bdt_*`) belong to the frozen model and are not
drawn here.

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
