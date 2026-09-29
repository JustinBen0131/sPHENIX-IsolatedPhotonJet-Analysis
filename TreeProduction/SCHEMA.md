# eventTree v1

*Generated from `TreeProduction/config/eventTree_v1.yaml` by `TreeProduction/tests/make_schema_doc.py`; do not edit.*

One TTree per file. DATA files are one run segment; SIM files are one product part. One entry per event occurrence, in source order (evt_src_entry). Per-event quantities are scalars; each collection is a group of equal-length vectors named <group>_<field>; sub-objects of photons are flat vectors with a parent index or a fixed stride; links are event-local indices (-1 = none).

## Rules

- **Restructure.** Pure restructure of the production final layout. Dropped: 128-bit row identities (rebuilt from row order and indices), per-row strings (stored as integer codes with a dictionary in meta/codes), per-file and per-collection constants (meta), and exact duplicates of another stored branch. Everything else, computed values included, is stored at its native type. Every drop is verified per file by the converter, which fails closed, and the inverse rebuilds every original table bit for bit.
- **Storage.** Collections (the vector<T> types in this schema) are C arrays sized by an Int_t counter: n_<group> for the group, or n_<branch> for a member whose length differs from the rest of its group (a strided or empty member). meta/schema records the branch-to-counter map. Hash strings are std::vector<std::string>. Files are written with LZMA level 4.
- **Entries.** evt_record = 0 for retained events and 1 for upstream-rejected events (DATA); rejected entries carry only evt_ fields.
- **Reserved.** Reserved branches are declared by the schema before they are filled. A file that has not been augmented carries them empty; meta/schema lists which reserved groups a file fills.
- **Versioning.** Minor versions only append branches or groups. Readers ignore branches they do not know. Nothing is ever renamed, retyped or removed within a major version.
- **Reader predicate** for analysis events: `evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && evt_terminal_status == 0`.
- **Units.** energies, pT, ET: GeV; z: cm; angles: rad; times: ns; raw times: ADC samples (x17.6 = ns); centrality: percent; MBD charge: a.u.

## Collections

| System | Reco jets | Truth jets |
|---|---|---|
| pp | jet02 (towerinfo_calib, R=0.2), jet03 (R=0.3), jet04 (R=0.4) | truthjet02, truthjet03, truthjet04 |
| Au+Au | jet02/jet03/jet04 (towerinfo_sub1, UE-subtracted), jetunsub02/03/04 (towerinfo_retower_unsubtracted) | truthjet02, truthjet03, truthjet04 |

Analysis jets carry the same names in pp and Au+Au so one analysis chain reads both; the view is recorded in meta/collections.

Shower views: `H70`, `H0`, `G70`, `G0`, `O70`, `O0`, `R70` (stride 7); shv_ branches have 7 entries per photon in this order; the canonical photon-ID features are the H70 view.

Isolation witnesses: radii [0.3, 0.4] (stride 2).

`evt_terminal_status` codes are system specific:

- pp: 0 none, 1 trigger, 2 vz, 3 period, 4 ownership.
- auau: 0 none, 1 trigger, 2 vz, 3 minimum_bias, 4 stitch, 5 centrality.

System-specific codes; 0 also marks an early return (evt_nodes_ready == 0).

## Branch counts

Every file carries all 642 branches (the superset schema; a group that does not apply to a product is empty) plus 13 reserved branches. Branches with entries in current production, per product:

| Product | Populated branches |
|---|---:|
| pp_data | 578 |
| auau_data | 549 |
| pp_si_photon | 534 |
| pp_si_jet | 534 |
| pp_di | 614 |
| auau_emb_photon | 516 |
| auau_emb_jet | 516 |

## Branches

### Event identity and bookkeeping (`evt_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `evt_aug_donor_evidence` | int32 |  | donor evidence index | `eventAugmentations.donor_evidence_index` |
| `evt_run` | int32 |  | run | `events.run` |
| `evt_seq` | int64 |  | event sequence | `events.event_sequence` |
| `evt_phys_seq` | int64 |  | physical event sequence | `events.physical_event_sequence` |
| `evt_collab_iface_version` | int32 |  | collaborator interface version | `events.collaborator_interface_version` |
| `evt_embedded_mb_decision` | int32 |  | embedded minbias decision | `events.embedded_minbias_decision` |
| `evt_nodes_ready` | int32 |  | nodes ready | `events.nodes_ready` |
| `evt_terminal_status` | int32 |  | terminal status | `events.terminal_status` |
| `evt_n_candidates` | int32 |  | candidate count | `events.candidate_count` |
| `evt_n_tags` | int32 |  | tag count | `events.tag_count` |
| `evt_n_recoils` | int32 |  | recoil count | `events.recoil_count` |
| `evt_src_file_ordinal` | int64 |  | source file ordinal | `events.source_file_ordinal` |
| `evt_src_entry` | int64 |  | source entry ordinal | `events.source_entry_ordinal` |
| `evt_src_global_entry` | int64 |  | source global entry ordinal | `events.source_global_entry_ordinal` |
| `evt_phys_seq_valid` | int32 |  | physical event sequence valid | `events.physical_event_sequence_valid` |
| `evt_record` | int32 |  | 0 = retained event; 1 = upstream-rejected event (DATA) | `new` |

### Trigger (`trig_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `trig_live` | uint64 |  | live trigger bits | `events.live_trigger_bits` |
| `trig_scaled` | uint64 |  | scaled trigger bits | `events.scaled_trigger_bits` |
| `trig_packet_status` | uint64 |  | trigger packet status | `events.trigger_packet_status` |
| `trig_capture_version` | int32 |  | trigger capture version | `events.trigger_capture_version` |
| `trig_decision_state` | int32 |  | trigger decision state | `events.trigger_decision_state` |
| `trig_decision_available` | uint32 |  | trigger decision available | `events.trigger_decision_available` |
| `trig_packet_version` | int32 |  | trigger packet version | `events.trigger_packet_version` |

### GL1 scalers (`gl1_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `gl1_bco` | uint64 | GL1 beam clock count | first bco [GL1 beam clock count] | `triggerScalers.first_bco` |
| `gl1_state` | int32 |  | state | `triggerScalers.state` |
| `gl1_discontinuity` | uint32 |  | discontinuity flags | `triggerScalers.discontinuity_flags` |
| `gl1_raw` | uint64[64] | counts | raw [counts] | `triggerScalers.raw` |
| `gl1_live` | uint64[64] | counts | live [counts] | `triggerScalers.live` |
| `gl1_scaled` | uint64[64] | counts | scaled [counts] | `triggerScalers.scaled` |

### Vertex (`vtx_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `vtx_z` | double | cm | vertex z [cm] | `events.vertex_z` |
| `vtx_valid` | int32 |  | reco vertex valid | `events.reco_vertex_valid` |

### Centrality (`cent_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `cent_value` | double | percent (0-100) | centrality [percent (0-100)] | `events.centrality` |
| `cent_replay_version` | int32 |  | centrality replay version | `events.centrality_replay_version` |
| `cent_pmt_available` | int32 |  | centrality pmt available | `events.centrality_pmt_available` |
| `cent_inputs_valid` | int32 |  | centrality inputs valid | `events.centrality_inputs_valid` |
| `cent_mb_decision` | int32 |  | centrality mb decision | `events.centrality_mb_decision` |
| `cent_mbd_z_valid` | int32 | cm | centrality mbd z valid [cm] | `events.centrality_mbd_z_valid` |
| `cent_mbd_event` | int32 |  | centrality mbd event | `events.centrality_mbd_event` |
| `cent_mbd_clock` | int32 |  | centrality mbd clock | `events.centrality_mbd_clock` |
| `cent_mbd_femclock` | int32 |  | centrality mbd femclock | `events.centrality_mbd_femclock` |
| `cent_native_bin` | int32 |  | centrality native bin | `events.centrality_native_bin` |
| `cent_native_valid` | int32 |  | centrality native valid | `events.centrality_native_valid` |
| `cent_mbd_z` | double | cm | centrality mbd z [cm] | `events.centrality_mbd_z` |
| `cent_selected_charge` | double | MBD calibrated charge (a.u.) | centrality selected charge [MBD calibrated charge (a.u.)] | `events.centrality_selected_charge` |
| `cent_native_centile` | double | fraction (0-1) | centrality native centile [fraction (0-1)] | `events.centrality_native_centile` |
| `cent_pmt_time` | vector<double> | ns (MBD PMT time) | centrality pmt time [ns (MBD PMT time)] | `events.centrality_pmt_time` |
| `cent_pmt_valid` | vector<int32> |  | centrality pmt valid | `events.centrality_pmt_valid` |
| `cent_pmt_selected` | vector<int32> |  | centrality pmt selected | `events.centrality_pmt_selected` |

### MBD (`mbd_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `mbd_native_valid` | int32 |  | mbd native is valid | `eventAugmentations.mbd_native_is_valid` |
| `mbd_npmt_north` | int32 |  | mbd npmt north | `eventAugmentations.mbd_npmt_north` |
| `mbd_npmt_south` | int32 |  | mbd npmt south | `eventAugmentations.mbd_npmt_south` |
| `mbd_t0` | double | ns | mbd t0 [ns] | `eventAugmentations.mbd_t0` |
| `mbd_t0_reason` | int32 |  | mbd t0 reason index | `eventAugmentations.mbd_t0_reason_index` |
| `mbd_t0_state` | int32 |  | mbd t0 state | `eventAugmentations.mbd_t0_state` |
| `mbd_time_north` | double | ns | mbd time north [ns] | `eventAugmentations.mbd_time_north` |
| `mbd_time_north_reason` | int32 |  | mbd time north reason index | `eventAugmentations.mbd_time_north_reason_index` |
| `mbd_time_north_state` | int32 |  | mbd time north state | `eventAugmentations.mbd_time_north_state` |
| `mbd_time_south` | double | ns | mbd time south [ns] | `eventAugmentations.mbd_time_south` |
| `mbd_time_south_reason` | int32 |  | mbd time south reason index | `eventAugmentations.mbd_time_south_reason_index` |
| `mbd_time_south_state` | int32 |  | mbd time south state | `eventAugmentations.mbd_time_south_state` |
| `mbd_charge_sum` | double | MBD calibrated charge (a.u.) | mbd total charge [MBD calibrated charge (a.u.)] | `eventAugmentations.mbd_total_charge` |
| `mbd_charge_sum_origin` | int32 |  | mbd total charge origin | `eventAugmentations.mbd_total_charge_origin` |
| `mbd_charge_sum_state` | int32 |  | mbd total charge state | `eventAugmentations.mbd_total_charge_state` |
| `mbd_total_charge_native` | double | MBD calibrated charge (a.u.) | mbd total charge [MBD calibrated charge (a.u.)] | `events.mbd_total_charge` |
| `mbd_pmt_available` | int32 |  | mbd pmt available | `events.mbd_pmt_available` |
| `mbd_pmt_q` | vector<double> | MBD calibrated charge (a.u.) | mbd pmt charge [MBD calibrated charge (a.u.)] | `events.mbd_pmt_charge` |
| `mbd_native_capture_version` | int32 |  | mbd native capture version | `events.mbd_native_capture_version` |
| `mbd_native_available` | int32 |  | mbd native available | `events.mbd_native_available` |
| `mbd_native_south_time_ns` | double | ns | mbd native south time ns [ns] | `events.mbd_native_south_time_ns` |
| `mbd_native_north_time_ns` | double | ns | mbd native north time ns [ns] | `events.mbd_native_north_time_ns` |
| `mbd_native_t0_finite` | int32 |  | mbd native t0 finite | `events.mbd_native_t0_finite` |
| `mbd_native_south_time_finite` | int32 |  | mbd native south time finite | `events.mbd_native_south_time_finite` |
| `mbd_native_north_time_finite` | int32 |  | mbd native north time finite | `events.mbd_native_north_time_finite` |
| `mbd_native_north_npmt` | int32 |  | mbd native north npmt | `events.mbd_native_north_npmt` |

### Calorimeter sums (`calo_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `calo_emcal` | double | GeV | emcal total energy [GeV] | `eventAugmentations.emcal_total_energy` |
| `calo_emcal_reason` | int32 |  | emcal total energy reason index | `eventAugmentations.emcal_total_energy_reason_index` |
| `calo_emcal_state` | int32 |  | emcal total energy state | `eventAugmentations.emcal_total_energy_state` |
| `calo_available_mask` | uint32 |  | event calo available mask | `eventAugmentations.event_calo_available_mask` |
| `calo_capture_version` | int32 |  | event calo capture version | `eventAugmentations.event_calo_capture_version` |
| `calo_require_isgood` | int32 |  | event calo require isgood | `eventAugmentations.event_calo_require_isgood` |
| `calo_valid_mask` | uint32 |  | event calo valid mask | `eventAugmentations.event_calo_valid_mask` |
| `calo_ihcal` | double | GeV | ihcal total energy [GeV] | `eventAugmentations.ihcal_total_energy` |
| `calo_ihcal_reason` | int32 |  | ihcal total energy reason index | `eventAugmentations.ihcal_total_energy_reason_index` |
| `calo_ihcal_state` | int32 |  | ihcal total energy state | `eventAugmentations.ihcal_total_energy_state` |
| `calo_ohcal` | double | GeV | ohcal total energy [GeV] | `eventAugmentations.ohcal_total_energy` |
| `calo_ohcal_reason` | int32 |  | ohcal total energy reason index | `eventAugmentations.ohcal_total_energy_reason_index` |
| `calo_ohcal_state` | int32 |  | ohcal total energy state | `eventAugmentations.ohcal_total_energy_state` |
| `calo_total` | double | GeV | total calo energy [GeV] | `eventAugmentations.total_calo_energy` |
| `calo_total_reason` | int32 |  | total calo energy reason index | `eventAugmentations.total_calo_energy_reason_index` |
| `calo_total_state` | int32 |  | total calo energy state | `eventAugmentations.total_calo_energy_state` |
| `calo_emcal_native` | double | GeV | emcal total energy [GeV] | `events.emcal_total_energy` |
| `calo_ihcal_native` | double | GeV | ihcal total energy [GeV] | `events.ihcal_total_energy` |
| `calo_ohcal_native` | double | GeV | ohcal total energy [GeV] | `events.ohcal_total_energy` |
| `calo_total_native` | double | GeV | total calo energy [GeV] | `events.total_calo_energy` |
| `calo_max_cluster_e` | double | GeV | max cluster energy [GeV] | `events.max_cluster_energy` |
| `calo_require_isgood_native` | int32 |  | event calo require isgood | `events.event_calo_require_isgood` |
| `calo_native_mask` | int32 |  | bit i set when the native calorimeter sum i (emcal, ihcal, ohcal, total) was captured natively | `new` |

### Underlying event (`ue_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `ue_contract_state` | int32 |  | jet background contract state | `events.jet_background_contract_state` |
| `ue_flow_mode` | int32 |  | jet background flow mode | `events.jet_background_flow_mode` |
| `ue_eta_bins` | int32 | GeV | jet background eta bins [GeV] | `events.jet_background_eta_bins` |
| `ue_v2` | double |  | jet background v2 | `events.jet_background_v2` |
| `ue_psi2` | double | rad | jet background psi2 [rad] | `events.jet_background_psi2` |
| `ue_emcal` | vector<float> |  | jet background ue emcal | `events.jet_background_ue_emcal` |
| `ue_ihcal` | vector<float> |  | jet background ue ihcal | `events.jet_background_ue_ihcal` |
| `ue_ohcal` | vector<float> |  | jet background ue ohcal | `events.jet_background_ue_ohcal` |
| `ue_n_strips` | int32 |  | jet background n strips | `events.jet_background_n_strips` |
| `ue_n_towers` | int32 |  | jet background n towers | `events.jet_background_n_towers` |
| `ue_flow_failure_state` | int32 |  | jet background flow failure state | `events.jet_background_flow_failure_state` |

### PPG12 event quantities (`ppg12_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `ppg12_weight_active_state` | int32 |  | ppg12 weight active state | `events.ppg12_weight_active_state` |
| `ppg12_weight_lane_code` | int32 |  | ppg12 weight lane component code | `events.ppg12_weight_lane_component_code` |
| `ppg12_wfill` | int64 |  | ppg12 weighted fill count | `events.ppg12_weighted_fill_count` |
| `ppg12_rfill` | int64 |  | ppg12 raw fill count | `events.ppg12_raw_fill_count` |
| `ppg12_wfill_code1` | int64 |  | ppg12 weighted fill count code1 | `events.ppg12_weighted_fill_count_code1` |
| `ppg12_rfill_code1` | int64 |  | ppg12 raw fill count code1 | `events.ppg12_raw_fill_count_code1` |
| `ppg12_wfill_code2` | int64 |  | ppg12 weighted fill count code2 | `events.ppg12_weighted_fill_count_code2` |
| `ppg12_rfill_code2` | int64 |  | ppg12 raw fill count code2 | `events.ppg12_raw_fill_count_code2` |
| `ppg12_wfill_code3` | int64 |  | ppg12 weighted fill count code3 | `events.ppg12_weighted_fill_count_code3` |
| `ppg12_rfill_code3` | int64 |  | ppg12 raw fill count code3 | `events.ppg12_raw_fill_count_code3` |
| `ppg12_wfill_code4` | int64 |  | ppg12 weighted fill count code4 | `events.ppg12_weighted_fill_count_code4` |
| `ppg12_rfill_code4` | int64 |  | ppg12 raw fill count code4 | `events.ppg12_raw_fill_count_code4` |
| `ppg12_arb_state` | int32 |  | ppg12 source arbitration state | `events.ppg12_source_arbitration_state` |
| `ppg12_lead_key` | uint32 |  | ppg12 leading source cluster key | `events.ppg12_leading_source_cluster_key` |
| `ppg12_lead_encounter` | int32 |  | ppg12 leading source encounter ordinal | `events.ppg12_leading_source_encounter_ordinal` |
| `ppg12_lead_et` | double | GeV | ppg12 leading source et [GeV] | `events.ppg12_leading_source_et` |
| `ppg12_lead_eta` | double | dimensionless (pseudorapidity) | ppg12 leading source eta [dimensionless (pseudorapidity)] | `events.ppg12_leading_source_eta` |
| `ppg12_lead_phi` | double | rad | ppg12 leading source phi [rad] | `events.ppg12_leading_source_phi` |

### Weights (`w_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `w_event` | double | dimensionless weight | event weight [dimensionless weight] | `events.event_weight` |
| `w_sample_valid` | int32 |  | sample weight valid | `events.sample_weight_valid` |
| `w_vertex_valid` | int32 |  | vertex weight valid | `events.vertex_weight_valid` |
| `w_legacy` | double | dimensionless weight | legacy event weight [dimensionless weight] | `events.legacy_event_weight` |
| `w_slice` | double | dimensionless weight | slice weight [dimensionless weight] | `weights.slice_weight` |
| `w_xsec` | double | dimensionless weight | cross section weight [dimensionless weight] | `weights.cross_section_weight` |
| `w_vertex` | double | dimensionless weight | vertex weight [dimensionless weight] | `weights.vertex_weight` |
| `w_sidi` | double | dimensionless weight | si di weight [dimensionless weight] | `weights.si_di_weight` |
| `w_period` | double | dimensionless weight | period weight [dimensionless weight] | `weights.period_weight` |
| `w_exposure` | double | dimensionless weight | exposure weight [dimensionless weight] | `weights.exposure_weight` |
| `w_application_count` | int32 |  | application count | `weights.application_count` |

### Photon capture bookkeeping (`recopho_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `recopho_capture_version` | int32 |  | reco photon capture version | `events.reco_photon_capture_version` |
| `recopho_capture_domain` | int32 |  | reco photon capture domain | `events.reco_photon_capture_domain` |
| `recopho_capture_state` | int32 |  | reco photon capture state | `events.reco_photon_capture_state` |
| `recopho_inspected` | int32 |  | reco photon inspected | `events.reco_photon_inspected` |
| `recopho_eligible` | int32 |  | reco photon eligible | `events.reco_photon_eligible` |
| `recopho_written` | int32 |  | reco photon written | `events.reco_photon_written` |
| `recopho_unclassifiable` | int32 |  | reco photon unclassifiable | `events.reco_photon_unclassifiable` |

### Jet views (`jetview_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `jetview_available` | vector<int32> |  | reco jet view available | `events.reco_jet_view_available` |

### Photons (`pho_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `pho_prodbdt_score` | vector<double> | dimensionless (classifier output) | raw score [dimensionless (classifier output)] | `models.raw_score` |
| `pho_prodbdt_finite` | vector<int32> | dimensionless (classifier output) | finite score [dimensionless (classifier output)] | `models.finite_score` |
| `pho_prodbdt_applicability` | vector<int32> |  | applicability state | `models.applicability_state` |
| `pho_prodbdt_wp70` | vector<double> | dimensionless (classifier output) | wp70 [dimensionless (classifier output)] | `models.wp70` |
| `pho_prodbdt_wp80` | vector<double> | dimensionless (classifier output) | wp80 [dimensionless (classifier output)] | `models.wp80` |
| `pho_prodbdt_wp90` | vector<double> | dimensionless (classifier output) | wp90 [dimensionless (classifier output)] | `models.wp90` |
| `pho_prodbdt_delta_wp70` | vector<double> | dimensionless (classifier output) | delta wp70 [dimensionless (classifier output)] | `models.delta_wp70` |
| `pho_prodbdt_delta_wp80` | vector<double> | dimensionless (classifier output) | delta wp80 [dimensionless (classifier output)] | `models.delta_wp80` |
| `pho_prodbdt_delta_wp90` | vector<double> | dimensionless (classifier output) | delta wp90 [dimensionless (classifier output)] | `models.delta_wp90` |
| `pho_assoc_truthpho` | vector<int32> |  | truth photon id lo | `photonAssociations.truth_photon_id_lo` |
| `pho_assoc_state` | vector<int32> |  | association state | `photonAssociations.association_state` |
| `pho_assoc_dr` | vector<double> | dimensionless (DeltaR) | diagnostic delta r [dimensionless (DeltaR)] | `photonAssociations.diagnostic_delta_r` |
| `pho_time_ntowers` | vector<int32> |  | cluster time contributing towers | `photonAugmentations.cluster_time_contributing_towers` |
| `pho_time_den_e` | vector<double> | GeV | cluster time denominator energy [GeV] | `photonAugmentations.cluster_time_denominator_energy` |
| `pho_time_den_e_reason` | vector<int32> |  | cluster time denominator energy reason index | `photonAugmentations.cluster_time_denominator_energy_reason_index` |
| `pho_time_den_e_state` | vector<int32> |  | cluster time denominator energy state | `photonAugmentations.cluster_time_denominator_energy_state` |
| `pho_time_ns` | vector<double> | ns | cluster time ns [ns] | `photonAugmentations.cluster_time_ns` |
| `pho_time_raw` | vector<double> | ADC time samples (x17.6 = ns) | cluster time raw [ADC time samples (x17.6 = ns)] | `photonAugmentations.cluster_time_raw` |
| `pho_time_raw_reason` | vector<int32> |  | cluster time raw reason index | `photonAugmentations.cluster_time_raw_reason_index` |
| `pho_time_raw_state` | vector<int32> |  | cluster time raw state | `photonAugmentations.cluster_time_raw_state` |
| `pho_aug_donor_evidence` | vector<int32> |  | donor evidence index | `photonAugmentations.donor_evidence_index` |
| `pho_npb_applicability` | vector<int32> |  | npb applicability state | `photonAugmentations.npb_applicability_state` |
| `pho_npb_cut` | vector<double> | dimensionless (classifier output) | npb configured cut [dimensionless (classifier output)] | `photonAugmentations.npb_configured_cut` |
| `pho_npb_evaluation` | vector<int32> |  | npb evaluation state | `photonAugmentations.npb_evaluation_state` |
| `pho_npb_cluster_et` | vector<double> | GeV | npb input cluster Et [GeV] | `photonAugmentations.npb_input_cluster_Et` |
| `pho_npb_cluster_et_evidence_index` | vector<int32> |  | npb input cluster Et evidence index | `photonAugmentations.npb_input_cluster_Et_evidence_index` |
| `pho_npb_cluster_et_origin_index` | vector<int32> |  | npb input cluster Et origin index | `photonAugmentations.npb_input_cluster_Et_origin_index` |
| `pho_npb_cluster_et_reason_index` | vector<int32> |  | npb input cluster Et reason index | `photonAugmentations.npb_input_cluster_Et_reason_index` |
| `pho_npb_cluster_et_state` | vector<int32> |  | npb input cluster Et state | `photonAugmentations.npb_input_cluster_Et_state` |
| `pho_npb_cluster_eta` | vector<double> | dimensionless (pseudorapidity) | npb input cluster Eta [dimensionless (pseudorapidity)] | `photonAugmentations.npb_input_cluster_Eta` |
| `pho_npb_cluster_eta_evidence_index` | vector<int32> |  | npb input cluster Eta evidence index | `photonAugmentations.npb_input_cluster_Eta_evidence_index` |
| `pho_npb_cluster_eta_origin_index` | vector<int32> |  | npb input cluster Eta origin index | `photonAugmentations.npb_input_cluster_Eta_origin_index` |
| `pho_npb_cluster_eta_reason_index` | vector<int32> |  | npb input cluster Eta reason index | `photonAugmentations.npb_input_cluster_Eta_reason_index` |
| `pho_npb_cluster_eta_state` | vector<int32> | dimensionless (pseudorapidity) | npb input cluster Eta state [dimensionless (pseudorapidity)] | `photonAugmentations.npb_input_cluster_Eta_state` |
| `pho_npb_cluster_et1` | vector<double> | GeV | npb input cluster et1 [GeV] | `photonAugmentations.npb_input_cluster_et1` |
| `pho_npb_cluster_et1_evidence_index` | vector<int32> |  | npb input cluster et1 evidence index | `photonAugmentations.npb_input_cluster_et1_evidence_index` |
| `pho_npb_cluster_et1_origin_index` | vector<int32> |  | npb input cluster et1 origin index | `photonAugmentations.npb_input_cluster_et1_origin_index` |
| `pho_npb_cluster_et1_reason_index` | vector<int32> |  | npb input cluster et1 reason index | `photonAugmentations.npb_input_cluster_et1_reason_index` |
| `pho_npb_cluster_et1_state` | vector<int32> |  | npb input cluster et1 state | `photonAugmentations.npb_input_cluster_et1_state` |
| `pho_npb_cluster_et2` | vector<double> | GeV | npb input cluster et2 [GeV] | `photonAugmentations.npb_input_cluster_et2` |
| `pho_npb_cluster_et2_evidence_index` | vector<int32> |  | npb input cluster et2 evidence index | `photonAugmentations.npb_input_cluster_et2_evidence_index` |
| `pho_npb_cluster_et2_origin_index` | vector<int32> |  | npb input cluster et2 origin index | `photonAugmentations.npb_input_cluster_et2_origin_index` |
| `pho_npb_cluster_et2_reason_index` | vector<int32> |  | npb input cluster et2 reason index | `photonAugmentations.npb_input_cluster_et2_reason_index` |
| `pho_npb_cluster_et2_state` | vector<int32> |  | npb input cluster et2 state | `photonAugmentations.npb_input_cluster_et2_state` |
| `pho_npb_cluster_et3` | vector<double> | GeV | npb input cluster et3 [GeV] | `photonAugmentations.npb_input_cluster_et3` |
| `pho_npb_cluster_et3_evidence_index` | vector<int32> |  | npb input cluster et3 evidence index | `photonAugmentations.npb_input_cluster_et3_evidence_index` |
| `pho_npb_cluster_et3_origin_index` | vector<int32> |  | npb input cluster et3 origin index | `photonAugmentations.npb_input_cluster_et3_origin_index` |
| `pho_npb_cluster_et3_reason_index` | vector<int32> |  | npb input cluster et3 reason index | `photonAugmentations.npb_input_cluster_et3_reason_index` |
| `pho_npb_cluster_et3_state` | vector<int32> |  | npb input cluster et3 state | `photonAugmentations.npb_input_cluster_et3_state` |
| `pho_npb_cluster_et4` | vector<double> | GeV | npb input cluster et4 [GeV] | `photonAugmentations.npb_input_cluster_et4` |
| `pho_npb_cluster_et4_evidence_index` | vector<int32> |  | npb input cluster et4 evidence index | `photonAugmentations.npb_input_cluster_et4_evidence_index` |
| `pho_npb_cluster_et4_origin_index` | vector<int32> |  | npb input cluster et4 origin index | `photonAugmentations.npb_input_cluster_et4_origin_index` |
| `pho_npb_cluster_et4_reason_index` | vector<int32> |  | npb input cluster et4 reason index | `photonAugmentations.npb_input_cluster_et4_reason_index` |
| `pho_npb_cluster_et4_state` | vector<int32> |  | npb input cluster et4 state | `photonAugmentations.npb_input_cluster_et4_state` |
| `pho_npb_cluster_w32` | vector<double> |  | npb input cluster w32 | `photonAugmentations.npb_input_cluster_w32` |
| `pho_npb_cluster_w32_evidence_index` | vector<int32> |  | npb input cluster w32 evidence index | `photonAugmentations.npb_input_cluster_w32_evidence_index` |
| `pho_npb_cluster_w32_origin_index` | vector<int32> |  | npb input cluster w32 origin index | `photonAugmentations.npb_input_cluster_w32_origin_index` |
| `pho_npb_cluster_w32_reason_index` | vector<int32> |  | npb input cluster w32 reason index | `photonAugmentations.npb_input_cluster_w32_reason_index` |
| `pho_npb_cluster_w32_state` | vector<int32> |  | npb input cluster w32 state | `photonAugmentations.npb_input_cluster_w32_state` |
| `pho_npb_cluster_w52` | vector<double> |  | npb input cluster w52 | `photonAugmentations.npb_input_cluster_w52` |
| `pho_npb_cluster_w52_evidence_index` | vector<int32> |  | npb input cluster w52 evidence index | `photonAugmentations.npb_input_cluster_w52_evidence_index` |
| `pho_npb_cluster_w52_origin_index` | vector<int32> |  | npb input cluster w52 origin index | `photonAugmentations.npb_input_cluster_w52_origin_index` |
| `pho_npb_cluster_w52_reason_index` | vector<int32> |  | npb input cluster w52 reason index | `photonAugmentations.npb_input_cluster_w52_reason_index` |
| `pho_npb_cluster_w52_state` | vector<int32> |  | npb input cluster w52 state | `photonAugmentations.npb_input_cluster_w52_state` |
| `pho_npb_cluster_w72` | vector<double> |  | npb input cluster w72 | `photonAugmentations.npb_input_cluster_w72` |
| `pho_npb_cluster_w72_evidence_index` | vector<int32> |  | npb input cluster w72 evidence index | `photonAugmentations.npb_input_cluster_w72_evidence_index` |
| `pho_npb_cluster_w72_origin_index` | vector<int32> |  | npb input cluster w72 origin index | `photonAugmentations.npb_input_cluster_w72_origin_index` |
| `pho_npb_cluster_w72_reason_index` | vector<int32> |  | npb input cluster w72 reason index | `photonAugmentations.npb_input_cluster_w72_reason_index` |
| `pho_npb_cluster_w72_state` | vector<int32> |  | npb input cluster w72 state | `photonAugmentations.npb_input_cluster_w72_state` |
| `pho_npb_cluster_weta_cogx` | vector<double> |  | npb input cluster weta cogx | `photonAugmentations.npb_input_cluster_weta_cogx` |
| `pho_npb_cluster_weta_cogx_evidence_index` | vector<int32> |  | npb input cluster weta cogx evidence index | `photonAugmentations.npb_input_cluster_weta_cogx_evidence_index` |
| `pho_npb_cluster_weta_cogx_origin_index` | vector<int32> |  | npb input cluster weta cogx origin index | `photonAugmentations.npb_input_cluster_weta_cogx_origin_index` |
| `pho_npb_cluster_weta_cogx_reason_index` | vector<int32> |  | npb input cluster weta cogx reason index | `photonAugmentations.npb_input_cluster_weta_cogx_reason_index` |
| `pho_npb_cluster_weta_cogx_state` | vector<int32> |  | npb input cluster weta cogx state | `photonAugmentations.npb_input_cluster_weta_cogx_state` |
| `pho_npb_cluster_wphi_cogx` | vector<double> |  | npb input cluster wphi cogx | `photonAugmentations.npb_input_cluster_wphi_cogx` |
| `pho_npb_cluster_wphi_cogx_evidence_index` | vector<int32> |  | npb input cluster wphi cogx evidence index | `photonAugmentations.npb_input_cluster_wphi_cogx_evidence_index` |
| `pho_npb_cluster_wphi_cogx_origin_index` | vector<int32> |  | npb input cluster wphi cogx origin index | `photonAugmentations.npb_input_cluster_wphi_cogx_origin_index` |
| `pho_npb_cluster_wphi_cogx_reason_index` | vector<int32> |  | npb input cluster wphi cogx reason index | `photonAugmentations.npb_input_cluster_wphi_cogx_reason_index` |
| `pho_npb_cluster_wphi_cogx_state` | vector<int32> |  | npb input cluster wphi cogx state | `photonAugmentations.npb_input_cluster_wphi_cogx_state` |
| `pho_npb_e11_over_e13` | vector<double> |  | npb input e11 over e13 | `photonAugmentations.npb_input_e11_over_e13` |
| `pho_npb_e11_over_e13_evidence_index` | vector<int32> |  | npb input e11 over e13 evidence index | `photonAugmentations.npb_input_e11_over_e13_evidence_index` |
| `pho_npb_e11_over_e13_origin_index` | vector<int32> |  | npb input e11 over e13 origin index | `photonAugmentations.npb_input_e11_over_e13_origin_index` |
| `pho_npb_e11_over_e13_reason_index` | vector<int32> |  | npb input e11 over e13 reason index | `photonAugmentations.npb_input_e11_over_e13_reason_index` |
| `pho_npb_e11_over_e13_state` | vector<int32> |  | npb input e11 over e13 state | `photonAugmentations.npb_input_e11_over_e13_state` |
| `pho_npb_e11_over_e15` | vector<double> |  | npb input e11 over e15 | `photonAugmentations.npb_input_e11_over_e15` |
| `pho_npb_e11_over_e15_evidence_index` | vector<int32> |  | npb input e11 over e15 evidence index | `photonAugmentations.npb_input_e11_over_e15_evidence_index` |
| `pho_npb_e11_over_e15_origin_index` | vector<int32> |  | npb input e11 over e15 origin index | `photonAugmentations.npb_input_e11_over_e15_origin_index` |
| `pho_npb_e11_over_e15_reason_index` | vector<int32> |  | npb input e11 over e15 reason index | `photonAugmentations.npb_input_e11_over_e15_reason_index` |
| `pho_npb_e11_over_e15_state` | vector<int32> |  | npb input e11 over e15 state | `photonAugmentations.npb_input_e11_over_e15_state` |
| `pho_npb_e11_over_e17` | vector<double> |  | npb input e11 over e17 | `photonAugmentations.npb_input_e11_over_e17` |
| `pho_npb_e11_over_e17_evidence_index` | vector<int32> |  | npb input e11 over e17 evidence index | `photonAugmentations.npb_input_e11_over_e17_evidence_index` |
| `pho_npb_e11_over_e17_origin_index` | vector<int32> |  | npb input e11 over e17 origin index | `photonAugmentations.npb_input_e11_over_e17_origin_index` |
| `pho_npb_e11_over_e17_reason_index` | vector<int32> |  | npb input e11 over e17 reason index | `photonAugmentations.npb_input_e11_over_e17_reason_index` |
| `pho_npb_e11_over_e17_state` | vector<int32> |  | npb input e11 over e17 state | `photonAugmentations.npb_input_e11_over_e17_state` |
| `pho_npb_e11_over_e22` | vector<double> |  | npb input e11 over e22 | `photonAugmentations.npb_input_e11_over_e22` |
| `pho_npb_e11_over_e22_evidence_index` | vector<int32> |  | npb input e11 over e22 evidence index | `photonAugmentations.npb_input_e11_over_e22_evidence_index` |
| `pho_npb_e11_over_e22_origin_index` | vector<int32> |  | npb input e11 over e22 origin index | `photonAugmentations.npb_input_e11_over_e22_origin_index` |
| `pho_npb_e11_over_e22_reason_index` | vector<int32> |  | npb input e11 over e22 reason index | `photonAugmentations.npb_input_e11_over_e22_reason_index` |
| `pho_npb_e11_over_e22_state` | vector<int32> |  | npb input e11 over e22 state | `photonAugmentations.npb_input_e11_over_e22_state` |
| `pho_npb_e11_over_e31` | vector<double> |  | npb input e11 over e31 | `photonAugmentations.npb_input_e11_over_e31` |
| `pho_npb_e11_over_e31_evidence_index` | vector<int32> |  | npb input e11 over e31 evidence index | `photonAugmentations.npb_input_e11_over_e31_evidence_index` |
| `pho_npb_e11_over_e31_origin_index` | vector<int32> |  | npb input e11 over e31 origin index | `photonAugmentations.npb_input_e11_over_e31_origin_index` |
| `pho_npb_e11_over_e31_reason_index` | vector<int32> |  | npb input e11 over e31 reason index | `photonAugmentations.npb_input_e11_over_e31_reason_index` |
| `pho_npb_e11_over_e31_state` | vector<int32> |  | npb input e11 over e31 state | `photonAugmentations.npb_input_e11_over_e31_state` |
| `pho_npb_e11_over_e33` | vector<double> |  | npb input e11 over e33 | `photonAugmentations.npb_input_e11_over_e33` |
| `pho_npb_e11_over_e33_evidence_index` | vector<int32> |  | npb input e11 over e33 evidence index | `photonAugmentations.npb_input_e11_over_e33_evidence_index` |
| `pho_npb_e11_over_e33_origin_index` | vector<int32> |  | npb input e11 over e33 origin index | `photonAugmentations.npb_input_e11_over_e33_origin_index` |
| `pho_npb_e11_over_e33_reason_index` | vector<int32> |  | npb input e11 over e33 reason index | `photonAugmentations.npb_input_e11_over_e33_reason_index` |
| `pho_npb_e11_over_e33_state` | vector<int32> |  | npb input e11 over e33 state | `photonAugmentations.npb_input_e11_over_e33_state` |
| `pho_npb_e11_over_e51` | vector<double> |  | npb input e11 over e51 | `photonAugmentations.npb_input_e11_over_e51` |
| `pho_npb_e11_over_e51_evidence_index` | vector<int32> |  | npb input e11 over e51 evidence index | `photonAugmentations.npb_input_e11_over_e51_evidence_index` |
| `pho_npb_e11_over_e51_origin_index` | vector<int32> |  | npb input e11 over e51 origin index | `photonAugmentations.npb_input_e11_over_e51_origin_index` |
| `pho_npb_e11_over_e51_reason_index` | vector<int32> |  | npb input e11 over e51 reason index | `photonAugmentations.npb_input_e11_over_e51_reason_index` |
| `pho_npb_e11_over_e51_state` | vector<int32> |  | npb input e11 over e51 state | `photonAugmentations.npb_input_e11_over_e51_state` |
| `pho_npb_e11_over_e71` | vector<double> |  | npb input e11 over e71 | `photonAugmentations.npb_input_e11_over_e71` |
| `pho_npb_e11_over_e71_evidence_index` | vector<int32> |  | npb input e11 over e71 evidence index | `photonAugmentations.npb_input_e11_over_e71_evidence_index` |
| `pho_npb_e11_over_e71_origin_index` | vector<int32> |  | npb input e11 over e71 origin index | `photonAugmentations.npb_input_e11_over_e71_origin_index` |
| `pho_npb_e11_over_e71_reason_index` | vector<int32> |  | npb input e11 over e71 reason index | `photonAugmentations.npb_input_e11_over_e71_reason_index` |
| `pho_npb_e11_over_e71_state` | vector<int32> |  | npb input e11 over e71 state | `photonAugmentations.npb_input_e11_over_e71_state` |
| `pho_npb_e22_over_e33` | vector<double> |  | npb input e22 over e33 | `photonAugmentations.npb_input_e22_over_e33` |
| `pho_npb_e22_over_e33_evidence_index` | vector<int32> |  | npb input e22 over e33 evidence index | `photonAugmentations.npb_input_e22_over_e33_evidence_index` |
| `pho_npb_e22_over_e33_origin_index` | vector<int32> |  | npb input e22 over e33 origin index | `photonAugmentations.npb_input_e22_over_e33_origin_index` |
| `pho_npb_e22_over_e33_reason_index` | vector<int32> |  | npb input e22 over e33 reason index | `photonAugmentations.npb_input_e22_over_e33_reason_index` |
| `pho_npb_e22_over_e33_state` | vector<int32> |  | npb input e22 over e33 state | `photonAugmentations.npb_input_e22_over_e33_state` |
| `pho_npb_e22_over_e35` | vector<double> |  | npb input e22 over e35 | `photonAugmentations.npb_input_e22_over_e35` |
| `pho_npb_e22_over_e35_evidence_index` | vector<int32> |  | npb input e22 over e35 evidence index | `photonAugmentations.npb_input_e22_over_e35_evidence_index` |
| `pho_npb_e22_over_e35_origin_index` | vector<int32> |  | npb input e22 over e35 origin index | `photonAugmentations.npb_input_e22_over_e35_origin_index` |
| `pho_npb_e22_over_e35_reason_index` | vector<int32> |  | npb input e22 over e35 reason index | `photonAugmentations.npb_input_e22_over_e35_reason_index` |
| `pho_npb_e22_over_e35_state` | vector<int32> |  | npb input e22 over e35 state | `photonAugmentations.npb_input_e22_over_e35_state` |
| `pho_npb_e22_over_e37` | vector<double> |  | npb input e22 over e37 | `photonAugmentations.npb_input_e22_over_e37` |
| `pho_npb_e22_over_e37_evidence_index` | vector<int32> |  | npb input e22 over e37 evidence index | `photonAugmentations.npb_input_e22_over_e37_evidence_index` |
| `pho_npb_e22_over_e37_origin_index` | vector<int32> |  | npb input e22 over e37 origin index | `photonAugmentations.npb_input_e22_over_e37_origin_index` |
| `pho_npb_e22_over_e37_reason_index` | vector<int32> |  | npb input e22 over e37 reason index | `photonAugmentations.npb_input_e22_over_e37_reason_index` |
| `pho_npb_e22_over_e37_state` | vector<int32> |  | npb input e22 over e37 state | `photonAugmentations.npb_input_e22_over_e37_state` |
| `pho_npb_e22_over_e53` | vector<double> |  | npb input e22 over e53 | `photonAugmentations.npb_input_e22_over_e53` |
| `pho_npb_e22_over_e53_evidence_index` | vector<int32> |  | npb input e22 over e53 evidence index | `photonAugmentations.npb_input_e22_over_e53_evidence_index` |
| `pho_npb_e22_over_e53_origin_index` | vector<int32> |  | npb input e22 over e53 origin index | `photonAugmentations.npb_input_e22_over_e53_origin_index` |
| `pho_npb_e22_over_e53_reason_index` | vector<int32> |  | npb input e22 over e53 reason index | `photonAugmentations.npb_input_e22_over_e53_reason_index` |
| `pho_npb_e22_over_e53_state` | vector<int32> |  | npb input e22 over e53 state | `photonAugmentations.npb_input_e22_over_e53_state` |
| `pho_npb_e32_over_e35` | vector<double> |  | npb input e32 over e35 | `photonAugmentations.npb_input_e32_over_e35` |
| `pho_npb_e32_over_e35_evidence_index` | vector<int32> |  | npb input e32 over e35 evidence index | `photonAugmentations.npb_input_e32_over_e35_evidence_index` |
| `pho_npb_e32_over_e35_origin_index` | vector<int32> |  | npb input e32 over e35 origin index | `photonAugmentations.npb_input_e32_over_e35_origin_index` |
| `pho_npb_e32_over_e35_reason_index` | vector<int32> |  | npb input e32 over e35 reason index | `photonAugmentations.npb_input_e32_over_e35_reason_index` |
| `pho_npb_e32_over_e35_state` | vector<int32> |  | npb input e32 over e35 state | `photonAugmentations.npb_input_e32_over_e35_state` |
| `pho_npb_vertexz` | vector<double> | cm | npb input vertexz [cm] | `photonAugmentations.npb_input_vertexz` |
| `pho_npb_vertexz_evidence_index` | vector<int32> |  | npb input vertexz evidence index | `photonAugmentations.npb_input_vertexz_evidence_index` |
| `pho_npb_vertexz_origin_index` | vector<int32> |  | npb input vertexz origin index | `photonAugmentations.npb_input_vertexz_origin_index` |
| `pho_npb_vertexz_reason_index` | vector<int32> |  | npb input vertexz reason index | `photonAugmentations.npb_input_vertexz_reason_index` |
| `pho_npb_vertexz_state` | vector<int32> |  | npb input vertexz state | `photonAugmentations.npb_input_vertexz_state` |
| `pho_npb_pass` | vector<int32> |  | npb pass | `photonAugmentations.npb_pass` |
| `pho_npb_pass_reason` | vector<int32> |  | npb pass reason index | `photonAugmentations.npb_pass_reason_index` |
| `pho_npb_pass_state` | vector<int32> |  | npb pass state | `photonAugmentations.npb_pass_state` |
| `pho_npb_score` | vector<double> | dimensionless (classifier output) | npb score [dimensionless (classifier output)] | `photonAugmentations.npb_score` |
| `pho_npb_score_reason` | vector<int32> | dimensionless (classifier output) | npb score reason index [dimensionless (classifier output)] | `photonAugmentations.npb_score_reason_index` |
| `pho_npb_score_state` | vector<int32> | dimensionless (classifier output) | npb score state [dimensionless (classifier output)] | `photonAugmentations.npb_score_state` |
| `pho_npb_scoring_policy` | vector<int32> |  | npb scoring binding policy | `photonAugmentations.npb_scoring_binding_policy` |
| `pho_npb_scoring_pt` | vector<double> | GeV | npb scoring cluster pt [GeV] | `photonAugmentations.npb_scoring_cluster_pt` |
| `pho_npb_scoring_encounter` | vector<int64> |  | npb scoring encounter ordinal | `photonAugmentations.npb_scoring_encounter_ordinal` |
| `pho_npb_scoring_eta` | vector<double> | dimensionless (pseudorapidity) | npb scoring eta [dimensionless (pseudorapidity)] | `photonAugmentations.npb_scoring_eta` |
| `pho_npb_scoring_evidence` | vector<int32> |  | npb scoring evidence index | `photonAugmentations.npb_scoring_evidence_index` |
| `pho_npb_scoring_et_bits` | vector<uint32> |  | npb scoring input et bits | `photonAugmentations.npb_scoring_input_et_bits` |
| `pho_npb_scoring_eta_bits` | vector<uint32> | dimensionless (pseudorapidity) | npb scoring input eta bits [dimensionless (pseudorapidity)] | `photonAugmentations.npb_scoring_input_eta_bits` |
| `pho_npb_scoring_key` | vector<uint32> |  | npb scoring native cluster key | `photonAugmentations.npb_scoring_native_cluster_key` |
| `pho_npb_scoring_node` | vector<int32> |  | npb scoring node index | `photonAugmentations.npb_scoring_node_index` |
| `pho_npb_scoring_phi` | vector<double> | rad | npb scoring phi [rad] | `photonAugmentations.npb_scoring_phi` |
| `pho_resp_class` | vector<int32> |  | link class | `photonResponseLinks.link_class` |
| `pho_encounter` | vector<int32> |  | encounter ordinal | `photons.encounter_ordinal` |
| `pho_rank_keys` | vector<double> |  | rank keys | `photons.rank_keys` |
| `pho_et` | vector<double> | GeV | cluster et [GeV] | `photons.cluster_et` |
| `pho_eta` | vector<double> | dimensionless (pseudorapidity) | eta [dimensionless (pseudorapidity)] | `photons.eta` |
| `pho_phi` | vector<double> | rad | phi [rad] | `photons.phi` |
| `pho_weta_cogx` | vector<double> |  | native weta cogx | `photons.native_weta_cogx` |
| `pho_wphi_cogx` | vector<double> |  | native wphi cogx | `photons.native_wphi_cogx` |
| `pho_weta33_cogx` | vector<double> |  | native weta33 cogx | `photons.native_weta33_cogx` |
| `pho_wphi33_cogx` | vector<double> |  | native wphi33 cogx | `photons.native_wphi33_cogx` |
| `pho_weta35_cogx` | vector<double> |  | native weta35 cogx | `photons.native_weta35_cogx` |
| `pho_wphi53_cogx` | vector<double> |  | native wphi53 cogx | `photons.native_wphi53_cogx` |
| `pho_et1` | vector<double> | GeV | native et1 [GeV] | `photons.native_et1` |
| `pho_e11_over_e33` | vector<double> |  | native e11 over e33 | `photons.native_e11_over_e33` |
| `pho_e32_over_e35` | vector<double> |  | native e32 over e35 | `photons.native_e32_over_e35` |
| `pho_npb_score_native` | vector<double> | dimensionless (classifier output) | npb score [dimensionless (classifier output)] | `photons.npb_score` |
| `pho_auau_npb_score_native` | vector<double> | dimensionless (classifier output) | auau npb score [dimensionless (classifier output)] | `photons.auau_npb_score` |
| `pho_prodbdt_inputs` | vector<float> |  | ordered features | `photons.ordered_features` |
| `pho_presel_bits` | vector<uint64> |  | preselection bitmask | `photons.preselection_bitmask` |
| `pho_finite_feature_state` | vector<int32> |  | finite feature state | `photons.finite_feature_state` |
| `pho_below15_state` | vector<int32> |  | below15 retention state | `photons.below15_retention_state` |
| `pho_presel_reference` | vector<int32> |  | reference preselection state | `photons.reference_preselection_state` |
| `pho_presel_active` | vector<int32> |  | active preselection state | `photons.active_preselection_state` |
| `pho_truth_match_state` | vector<int32> |  | truth signal match state | `photons.truth_signal_match_state` |
| `pho_truth_match_barcode` | vector<int32> |  | truth signal match barcode | `photons.truth_signal_match_barcode` |
| `pho_ppg12_tower_mask_state` | vector<int32> |  | ppg12 tower mask state | `photons.ppg12_tower_mask_state` |
| `pho_ppg12_source_eligible` | vector<int32> |  | ppg12 source eligible state | `photons.ppg12_source_eligible_state` |
| `pho_key` | vector<uint32> |  | native cluster key | `photons.native_cluster_key` |
| `pho_dom_state` | vector<int32> |  | dominant truth state | `photons.dominant_truth_state` |
| `pho_dom_evaluator_mode` | vector<int32> |  | dominant truth evaluator mode | `photons.dominant_truth_evaluator_mode` |
| `pho_dom_track_id` | vector<int32> |  | dominant truth track id | `photons.dominant_truth_track_id` |
| `pho_dom_pid` | vector<int32> |  | dominant truth pid | `photons.dominant_truth_pid` |
| `pho_dom_barcode` | vector<int32> |  | dominant truth barcode | `photons.dominant_truth_barcode` |
| `pho_dom_embedding_id` | vector<int32> |  | dominant truth embedding id | `photons.dominant_truth_embedding_id` |
| `pho_dom_energy` | vector<double> |  | dominant truth energy contribution | `photons.dominant_truth_energy_contribution` |
| `pho_dom_vertex_id` | vector<int32> |  | dominant truth vertex id | `photons.dominant_truth_vertex_id` |
| `pho_ntime_version` | vector<int32> |  | native timing capture version | `photons.native_timing_capture_version` |
| `pho_ntime_den` | vector<double> | GeV | native time energy denominator [GeV] | `photons.native_time_energy_denominator` |
| `pho_ntime_num` | vector<double> | GeV | native time energy numerator [GeV] | `photons.native_time_energy_numerator` |
| `pho_ntime_ntowers` | vector<int32> |  | native time contributing towers | `photons.native_time_contributing_towers` |
| `pho_ntime_finite` | vector<int32> |  | native time finite | `photons.native_time_finite` |
| `pho_producer_vtx_z` | vector<double> | cm | producer vertex z [cm] | `photons.producer_vertex_z` |
| `pho_npbcap_association` | vector<int32> |  | npb capture association | `photons.npb_capture_association` |
| `pho_npbcap_configured_cut` | vector<double> | dimensionless (classifier output) | npb configured cut [dimensionless (classifier output)] | `photons.npb_configured_cut` |
| `pho_aunpbcap_configured_cut` | vector<double> | dimensionless (classifier output) | auau npb configured cut [dimensionless (classifier output)] | `photons.auau_npb_configured_cut` |
| `pho_npbcap_capture_version` | vector<int32> |  | npb capture capture version | `photons.npb_capture_capture_version` |
| `pho_npbcap_evaluation_state` | vector<int32> |  | npb capture evaluation state | `photons.npb_capture_evaluation_state` |
| `pho_npbcap_domain_mode` | vector<int32> |  | npb capture domain mode | `photons.npb_capture_domain_mode` |
| `pho_npbcap_in_producer_domain` | vector<int32> |  | npb capture in producer domain | `photons.npb_capture_in_producer_domain` |
| `pho_npbcap_model_available` | vector<int32> |  | npb capture model available | `photons.npb_capture_model_available` |
| `pho_npbcap_input_order_is_npb25` | vector<int32> |  | npb capture input order is npb25 | `photons.npb_capture_input_order_is_npb25` |
| `pho_npbcap_inputs_captured` | vector<int32> |  | npb capture inputs captured | `photons.npb_capture_inputs_captured` |
| `pho_npbcap_inputs_finite` | vector<int32> |  | npb capture inputs finite | `photons.npb_capture_inputs_finite` |
| `pho_npbcap_min_et` | vector<float> | GeV | npb capture min et [GeV] | `photons.npb_capture_min_et` |
| `pho_npbcap_max_et` | vector<float> | GeV | npb capture max et [GeV] | `photons.npb_capture_max_et` |
| `pho_npbcap_max_abs_eta` | vector<float> | dimensionless (pseudorapidity) | npb capture max abs eta [dimensionless (pseudorapidity)] | `photons.npb_capture_max_abs_eta` |
| `pho_npbcap_raw_score` | vector<float> | dimensionless (classifier output) | npb capture raw score [dimensionless (classifier output)] | `photons.npb_capture_raw_score` |
| `pho_npbcap_source_photon_id` | vector<uint32> |  | npb capture source photon id | `photons.npb_capture_source_photon_id` |
| `pho_npbcap_source_cluster_et` | vector<float> | GeV | npb capture source cluster et [GeV] | `photons.npb_capture_source_cluster_et` |
| `pho_npbcap_source_eta` | vector<float> | dimensionless (pseudorapidity) | npb capture source eta [dimensionless (pseudorapidity)] | `photons.npb_capture_source_eta` |
| `pho_npbcap_source_phi` | vector<float> | rad | npb capture source phi [rad] | `photons.npb_capture_source_phi` |
| `pho_npbcap_binding_policy` | vector<int32> |  | npb capture binding policy | `photons.npb_capture_binding_policy` |
| `pho_npbcap_source_identity_available` | vector<int32> |  | npb capture source identity available | `photons.npb_capture_source_identity_available` |
| `pho_npbcap_source_native_cluster_key` | vector<uint32> |  | npb capture source native cluster key | `photons.npb_capture_source_native_cluster_key` |
| `pho_npbcap_source_encounter_ordinal` | vector<int64> |  | npb capture source encounter ordinal | `photons.npb_capture_source_encounter_ordinal` |
| `pho_npbcap_kinematic_input_bits_available` | vector<int32> |  | npb capture kinematic input bits available | `photons.npb_capture_kinematic_input_bits_available` |
| `pho_npbcap_input_et_bits` | vector<uint32> |  | npb capture input et bits | `photons.npb_capture_input_et_bits` |
| `pho_npbcap_input_eta_bits` | vector<uint32> | dimensionless (pseudorapidity) | npb capture input eta bits [dimensionless (pseudorapidity)] | `photons.npb_capture_input_eta_bits` |
| `pho_npbcap_inputs` | vector<float> |  | npb capture ordered inputs | `photons.npb_capture_ordered_inputs` |
| `pho_npbcap_input_finite` | vector<int32> |  | npb capture ordered input finite | `photons.npb_capture_ordered_input_finite` |
| `pho_aunpbcap_capture_version` | vector<int32> |  | auau npb capture capture version | `photons.auau_npb_capture_capture_version` |
| `pho_aunpbcap_evaluation_state` | vector<int32> |  | auau npb capture evaluation state | `photons.auau_npb_capture_evaluation_state` |
| `pho_aunpbcap_domain_mode` | vector<int32> |  | auau npb capture domain mode | `photons.auau_npb_capture_domain_mode` |
| `pho_aunpbcap_in_producer_domain` | vector<int32> |  | auau npb capture in producer domain | `photons.auau_npb_capture_in_producer_domain` |
| `pho_aunpbcap_model_available` | vector<int32> |  | auau npb capture model available | `photons.auau_npb_capture_model_available` |
| `pho_aunpbcap_input_order_is_npb25` | vector<int32> |  | auau npb capture input order is npb25 | `photons.auau_npb_capture_input_order_is_npb25` |
| `pho_aunpbcap_inputs_captured` | vector<int32> |  | auau npb capture inputs captured | `photons.auau_npb_capture_inputs_captured` |
| `pho_aunpbcap_inputs_finite` | vector<int32> |  | auau npb capture inputs finite | `photons.auau_npb_capture_inputs_finite` |
| `pho_aunpbcap_min_et` | vector<float> | GeV | auau npb capture min et [GeV] | `photons.auau_npb_capture_min_et` |
| `pho_aunpbcap_max_et` | vector<float> | GeV | auau npb capture max et [GeV] | `photons.auau_npb_capture_max_et` |
| `pho_aunpbcap_max_abs_eta` | vector<float> | dimensionless (pseudorapidity) | auau npb capture max abs eta [dimensionless (pseudorapidity)] | `photons.auau_npb_capture_max_abs_eta` |
| `pho_aunpbcap_raw_score` | vector<float> | dimensionless (classifier output) | auau npb capture raw score [dimensionless (classifier output)] | `photons.auau_npb_capture_raw_score` |
| `pho_aunpbcap_source_photon_id` | vector<uint32> |  | auau npb capture source photon id | `photons.auau_npb_capture_source_photon_id` |
| `pho_aunpbcap_source_cluster_et` | vector<float> | GeV | auau npb capture source cluster et [GeV] | `photons.auau_npb_capture_source_cluster_et` |
| `pho_aunpbcap_source_eta` | vector<float> | dimensionless (pseudorapidity) | auau npb capture source eta [dimensionless (pseudorapidity)] | `photons.auau_npb_capture_source_eta` |
| `pho_aunpbcap_source_phi` | vector<float> | rad | auau npb capture source phi [rad] | `photons.auau_npb_capture_source_phi` |
| `pho_aunpbcap_binding_policy` | vector<int32> |  | auau npb capture binding policy | `photons.auau_npb_capture_binding_policy` |
| `pho_aunpbcap_source_identity_available` | vector<int32> |  | auau npb capture source identity available | `photons.auau_npb_capture_source_identity_available` |
| `pho_aunpbcap_source_native_cluster_key` | vector<uint32> |  | auau npb capture source native cluster key | `photons.auau_npb_capture_source_native_cluster_key` |
| `pho_aunpbcap_source_encounter_ordinal` | vector<int64> |  | auau npb capture source encounter ordinal | `photons.auau_npb_capture_source_encounter_ordinal` |
| `pho_aunpbcap_kinematic_input_bits_available` | vector<int32> |  | auau npb capture kinematic input bits available | `photons.auau_npb_capture_kinematic_input_bits_available` |
| `pho_aunpbcap_input_et_bits` | vector<uint32> |  | auau npb capture input et bits | `photons.auau_npb_capture_input_et_bits` |
| `pho_aunpbcap_input_eta_bits` | vector<uint32> | dimensionless (pseudorapidity) | auau npb capture input eta bits [dimensionless (pseudorapidity)] | `photons.auau_npb_capture_input_eta_bits` |
| `pho_aunpbcap_inputs` | vector<float> |  | auau npb capture ordered inputs | `photons.auau_npb_capture_ordered_inputs` |
| `pho_aunpbcap_input_finite` | vector<int32> |  | auau npb capture ordered input finite | `photons.auau_npb_capture_ordered_input_finite` |

### Photon shower views (stride 7) (`shv_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `shv_features` | vector<float> |  | ordered features | `showerViews.ordered_features` |
| `shv_cog_eta` | vector<double> | dimensionless (pseudorapidity) | cog eta [dimensionless (pseudorapidity)] | `showerViews.cog_eta` |
| `shv_cog_phi` | vector<double> | rad | cog phi [rad] | `showerViews.cog_phi` |
| `shv_raw_center_eta` | vector<double> | dimensionless (pseudorapidity) | raw center eta [dimensionless (pseudorapidity)] | `showerViews.raw_center_eta` |
| `shv_raw_center_phi` | vector<double> | rad | raw center phi [rad] | `showerViews.raw_center_phi` |
| `shv_center_eta_index` | vector<int32> |  | center eta index | `showerViews.center_eta_index` |
| `shv_center_phi_index` | vector<int32> | rad | center phi index [rad] | `showerViews.center_phi_index` |
| `shv_e11` | vector<double> | GeV | e11 [GeV] | `showerViews.e11` |
| `shv_e33` | vector<double> | GeV | e33 [GeV] | `showerViews.e33` |
| `shv_e32` | vector<double> | GeV | e32 [GeV] | `showerViews.e32` |
| `shv_e35` | vector<double> | GeV | e35 [GeV] | `showerViews.e35` |
| `shv_e11_over_e33` | vector<double> |  | e11 over e33 | `showerViews.e11_over_e33` |
| `shv_e32_over_e35` | vector<double> |  | e32 over e35 | `showerViews.e32_over_e35` |
| `shv_weta_cogx` | vector<double> |  | weta cogx | `showerViews.weta_cogx` |
| `shv_wphi_cogx` | vector<double> | dimensionless (classifier output) | wphi cogx [dimensionless (classifier output)] | `showerViews.wphi_cogx` |
| `shv_weta33_cogx` | vector<double> |  | weta33 cogx | `showerViews.weta33_cogx` |
| `shv_wphi33_cogx` | vector<double> | dimensionless (classifier output) | wphi33 cogx [dimensionless (classifier output)] | `showerViews.wphi33_cogx` |
| `shv_et1` | vector<double> | GeV | native et1 [GeV] | `showerViews.native_et1` |
| `shv_et2` | vector<double> | GeV | native et2 [GeV] | `showerViews.native_et2` |
| `shv_et3` | vector<double> | GeV | native et3 [GeV] | `showerViews.native_et3` |
| `shv_et4` | vector<double> | GeV | native et4 [GeV] | `showerViews.native_et4` |
| `shv_moment_eta_numerator` | vector<double> | dimensionless (pseudorapidity) | moment eta numerator [dimensionless (pseudorapidity)] | `showerViews.moment_eta_numerator` |
| `shv_moment_phi_numerator` | vector<double> | rad | moment phi numerator [rad] | `showerViews.moment_phi_numerator` |
| `shv_moment_denominator` | vector<double> | GeV | moment denominator [GeV] | `showerViews.moment_denominator` |
| `shv_moment33_eta_numerator` | vector<double> | dimensionless (pseudorapidity) | moment33 eta numerator [dimensionless (pseudorapidity)] | `showerViews.moment33_eta_numerator` |
| `shv_moment33_phi_numerator` | vector<double> | rad | moment33 phi numerator [rad] | `showerViews.moment33_phi_numerator` |
| `shv_moment33_denominator` | vector<double> | GeV | moment33 denominator [GeV] | `showerViews.moment33_denominator` |
| `shv_energy_source` | vector<int32> | GeV | energy source [GeV] | `showerViews.energy_source` |
| `shv_rectangular_membership` | vector<int32> |  | rectangular membership | `showerViews.rectangular_membership` |
| `shv_moment_membership` | vector<int32> | GeV | moment membership [GeV] | `showerViews.moment_membership` |
| `shv_finite_feature_state` | vector<int32> |  | finite feature state | `showerViews.finite_feature_state` |
| `shv_good_cell_count` | vector<int32> |  | good cell count | `showerViews.good_cell_count` |
| `shv_owned_cell_count` | vector<int32> |  | owned cell count | `showerViews.owned_cell_count` |
| `shv_active_sum_cell_count` | vector<int32> |  | active sum cell count | `showerViews.active_sum_cell_count` |
| `shv_active_moment_cell_count` | vector<int32> |  | active moment cell count | `showerViews.active_moment_cell_count` |
| `shv_exact_zero_count` | vector<int32> |  | exact zero count | `showerViews.exact_zero_count` |
| `shv_negative_count` | vector<int32> |  | negative count | `showerViews.negative_count` |
| `shv_nonfinite_count` | vector<int32> |  | nonfinite count | `showerViews.nonfinite_count` |

### Photon isolation witnesses (stride 2) (`isow_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `isow_subtraction_method` | vector<int32> |  | subtraction method | `isolationWitnesses.subtraction_method` |
| `isow_reconstructed_or_truth` | vector<int32> |  | reconstructed or truth | `isolationWitnesses.reconstructed_or_truth` |
| `isow_cone_sum` | vector<double> | GeV | cone sum [GeV] | `isolationWitnesses.cone_sum` |
| `isow_threshold` | vector<double> | GeV | threshold [GeV] | `isolationWitnesses.threshold` |
| `isow_sideband_threshold` | vector<double> | GeV | sideband threshold [GeV] | `isolationWitnesses.sideband_threshold` |
| `isow_pass_state` | vector<int32> |  | pass state | `isolationWitnesses.pass_state` |
| `isow_constituent_count` | vector<int32> |  | constituent count | `isolationWitnesses.constituent_count` |

### Isolation constituents (parent index isoc_pho) (`isoc_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `isoc_pho` | vector<int32> |  | candidate id lo | `isolationConstituents.candidate_id_lo` |
| `isoc_delta_eta` | vector<double> | dimensionless (pseudorapidity) | delta eta [dimensionless (pseudorapidity)] | `isolationConstituents.delta_eta` |
| `isoc_delta_phi` | vector<double> | rad | delta phi [rad] | `isolationConstituents.delta_phi` |
| `isoc_delta_r` | vector<double> | dimensionless (DeltaR) | delta r [dimensionless (DeltaR)] | `isolationConstituents.delta_r` |
| `isoc_subsystem` | vector<int32> |  | subsystem | `isolationConstituents.subsystem` |
| `isoc_raw_energy` | vector<double> | GeV | raw energy [GeV] | `isolationConstituents.raw_energy` |
| `isoc_calibrated_energy` | vector<double> | GeV | calibrated energy [GeV] | `isolationConstituents.calibrated_energy` |
| `isoc_sub1_energy` | vector<double> | GeV | sub1 energy [GeV] | `isolationConstituents.sub1_energy` |
| `isoc_phosub_residual` | vector<double> | GeV | phosub residual [GeV] | `isolationConstituents.phosub_residual` |
| `isoc_quality_state` | vector<int32> |  | quality state | `isolationConstituents.quality_state` |
| `isoc_mask_state` | vector<int32> |  | mask state | `isolationConstituents.mask_state` |
| `isoc_candidate_removal_state` | vector<int32> |  | candidate removal state | `isolationConstituents.candidate_removal_state` |
| `isoc_native_source` | vector<int32> |  | native source | `isolationConstituents.native_source` |
| `isoc_native_key` | vector<uint32> |  | native key | `isolationConstituents.native_key` |

### Shower cells (parent index cell_pho) (`cell_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `cell_pho` | vector<int32> |  | candidate id lo | `showerCells.candidate_id_lo` |
| `cell_local_ieta` | vector<int32> |  | local eta index | `showerCells.local_eta_index` |
| `cell_local_iphi` | vector<int32> | rad | local phi index [rad] | `showerCells.local_phi_index` |
| `cell_ieta` | vector<int32> |  | tower eta index | `showerCells.tower_eta_index` |
| `cell_iphi` | vector<int32> | rad | tower phi index [rad] | `showerCells.tower_phi_index` |
| `cell_e` | vector<double> | GeV | calibrated energy [GeV] | `showerCells.calibrated_energy` |
| `cell_rawcluster_value` | vector<double> |  | rawcluster map value | `showerCells.rawcluster_map_value` |
| `cell_is_good` | vector<int32> |  | is good | `showerCells.is_good` |
| `cell_is_zero` | vector<int32> |  | is zero | `showerCells.is_zero` |
| `cell_is_negative` | vector<int32> |  | is negative | `showerCells.is_negative` |
| `cell_is_nonfinite` | vector<int32> |  | is nonfinite | `showerCells.is_nonfinite` |
| `cell_seed_state` | vector<int32> |  | seed state | `showerCells.seed_state` |
| `cell_denominator_member` | vector<int32> |  | denominator membership | `showerCells.denominator_membership` |
| `cell_rawcluster_owned` | vector<int32> |  | rawcluster owned | `showerCells.rawcluster_owned` |
| `cell_rawcluster_value_present` | vector<int32> |  | rawcluster value present | `showerCells.rawcluster_value_present` |
| `cell_floor0_member` | vector<int32> |  | floor0 membership | `showerCells.floor0_membership` |
| `cell_floor70_member` | vector<int32> |  | floor70 membership | `showerCells.floor70_membership` |
| `cell_grid_bits` | vector<int32> |  | grid membership bitmask | `showerCells.grid_membership_bitmask` |
| `cell_time_version` | vector<int32> |  | tower time status capture version | `showerCells.tower_time_status_capture_version` |
| `cell_time_available` | vector<int32> |  | tower time status available | `showerCells.tower_time_status_available` |
| `cell_time_samples` | vector<double> | ADC time samples (x17.6 = ns) | tower time samples [ADC time samples (x17.6 = ns)] | `showerCells.tower_time_samples` |
| `cell_time_finite` | vector<int32> |  | tower time finite | `showerCells.tower_time_finite` |
| `cell_tower_status` | vector<uint32> |  | tower status | `showerCells.tower_status` |

### Reco jets, analysis view (`jet<RR>_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `jet<RR>_rawpt` | vector<double> | GeV | raw pt [GeV] | `jets.raw_pt` |
| `jet<RR>_eta` | vector<double> | dimensionless (pseudorapidity) | eta [dimensionless (pseudorapidity)] | `jets.eta` |
| `jet<RR>_phi` | vector<double> | rad | phi [rad] | `jets.phi` |
| `jet<RR>_area` | vector<double> | dimensionless (eta-phi area) | area [dimensionless (eta-phi area)] | `jets.area` |
| `jet<RR>_pt` | vector<double> | GeV | corrected pt [GeV] | `jets.corrected_pt` |
| `jet<RR>_jes_valid` | vector<int32> |  | 1 when the jet energy scale applied to this jet is valid; pair_xjgamma is NaN for invalid-JES jets | `new` |

### Reco jets, Au+Au unsubtracted view (`jetunsub<RR>_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `jetunsub<RR>_rawpt` | vector<double> | GeV | raw pt [GeV] | `jets.raw_pt` |
| `jetunsub<RR>_eta` | vector<double> | dimensionless (pseudorapidity) | eta [dimensionless (pseudorapidity)] | `jets.eta` |
| `jetunsub<RR>_phi` | vector<double> | rad | phi [rad] | `jets.phi` |
| `jetunsub<RR>_area` | vector<double> | dimensionless (eta-phi area) | area [dimensionless (eta-phi area)] | `jets.area` |
| `jetunsub<RR>_pt` | vector<double> | GeV | corrected pt [GeV] | `jets.corrected_pt` |

### Jet constituents (empty in current production) (`jcon_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `jcon_constituent_ordinal` | vector<int32> |  | constituent ordinal | `jetConstituents.constituent_ordinal` |
| `jcon_subsystem` | vector<int32> |  | subsystem | `jetConstituents.subsystem` |
| `jcon_energy` | vector<double> | GeV | energy [GeV] | `jetConstituents.energy` |
| `jcon_eta` | vector<double> | dimensionless (pseudorapidity) | eta [dimensionless (pseudorapidity)] | `jetConstituents.eta` |
| `jcon_phi` | vector<double> | rad | phi [rad] | `jetConstituents.phi` |
| `jcon_quality_state` | vector<int32> |  | quality state | `jetConstituents.quality_state` |
| `jcon_native_source` | vector<int32> |  | native source | `jetConstituents.native_source` |
| `jcon_native_key` | vector<uint32> |  | native key | `jetConstituents.native_key` |

### Photon-jet pairs (`pair_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `pair_delta_phi` | vector<double> | rad | delta phi [rad] | `photonJets.delta_phi` |
| `pair_recoil_state` | vector<int32> |  | recoil state | `photonJets.recoil_state` |
| `pair_photon_rank` | vector<int32> |  | photon rank | `photonJets.photon_rank` |
| `pair_jet_rank` | vector<int32> |  | jet rank | `photonJets.jet_rank` |
| `pair_wrong_photon_class` | vector<int32> |  | wrong photon class | `photonJets.wrong_photon_class` |
| `pair_wrong_recoil_class` | vector<int32> |  | wrong recoil class | `photonJets.wrong_recoil_class` |
| `pair_xjgamma` | vector<double> |  | xjgamma | `photonJets.xjgamma` |
| `pair_pho` | vector<int32> |  | photon index of the photon-jet pair | `new` |
| `pair_jet` | vector<int32> |  | jet index of the pair within the concatenated reco jet collections (meta/collections order) | `new` |

### Truth event (`truth_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truth_vtx_z` | double | cm | truth vertex z [cm] | `events.truth_vertex_z` |
| `truth_mbvtx_z` | double | cm | truth mb vertex z [cm] | `events.truth_mb_vertex_z` |
| `truth_denominator_complete` | int32 |  | truth denominator complete | `events.truth_denominator_complete` |

### Truth vertices (`truthvtx_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truthvtx_capture_state` | int32 |  | truth vertex capture state | `events.truth_vertex_capture_state` |
| `truthvtx_primary_id` | int32 |  | truth primary vertex id | `events.truth_primary_vertex_id` |
| `truthvtx_id` | vector<int32> |  | truth vertex id | `events.truth_vertex_id` |
| `truthvtx_embedding_id` | vector<int32> |  | truth vertex embedding id | `events.truth_vertex_embedding_id` |
| `truthvtx_embedding_valid` | vector<int32> |  | truth vertex embedding valid | `events.truth_vertex_embedding_valid` |
| `truthvtx_z` | vector<double> | cm | truth vertex z values [cm] | `events.truth_vertex_z_values` |
| `truthvtx_z_valid` | vector<int32> | cm | truth vertex z valid [cm] | `events.truth_vertex_z_valid` |

### Truth centrality (`truthcen_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truthcen_capture_state` | int32 |  | truth photon capture state | `events.truth_photon_capture_state` |
| `truthcen_native_embedded_primary_pid22_count` | int32 |  | truth photon native embedded primary pid22 count | `events.truth_photon_native_embedded_primary_pid22_count` |
| `truthcen_serialized_raw_count` | int32 |  | truth photon serialized raw count | `events.truth_photon_serialized_raw_count` |
| `truthcen_rejected_count` | int32 |  | truth photon rejected count | `events.truth_photon_rejected_count` |
| `truthcen_rejected_invalid_kinematics_count` | int32 |  | truth photon rejected invalid kinematics count | `events.truth_photon_rejected_invalid_kinematics_count` |
| `truthcen_rejected_nonfinite_isolation_count` | int32 |  | truth photon rejected nonfinite isolation count | `events.truth_photon_rejected_nonfinite_isolation_count` |
| `truthcen_duplicate_track_count` | int32 |  | truth photon duplicate track count | `events.truth_photon_duplicate_track_count` |
| `truthcen_isolation_input_incomplete_count` | int32 |  | truth photon isolation input incomplete count | `events.truth_photon_isolation_input_incomplete_count` |
| `truthcen_analysis_signal_count` | int32 |  | truth photon analysis signal count | `events.truth_photon_analysis_signal_count` |

### Truth-jet capture (`truthjetcap_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truthjetcap_state` | int32 |  | truth jet capture state | `events.truth_jet_capture_state` |
| `truthjetcap_container_valid` | vector<int32> |  | truth jet container valid | `events.truth_jet_container_valid` |

### Truth photons (`truthpho_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truthpho_resp_class` | vector<int32> |  | link class | `photonResponseLinks.link_class` |
| `truthpho_pt` | vector<double> | GeV | pt [GeV] | `truthPhotons.pt` |
| `truthpho_eta` | vector<double> | dimensionless (pseudorapidity) | eta [dimensionless (pseudorapidity)] | `truthPhotons.eta` |
| `truthpho_phi` | vector<double> | rad | phi [rad] | `truthPhotons.phi` |
| `truthpho_prompt_class` | vector<int32> | GeV | prompt class [GeV] | `truthPhotons.prompt_class` |
| `truthpho_source_role` | vector<int32> |  | source role | `truthPhotons.source_role` |
| `truthpho_g4_photon_valid` | vector<int32> |  | g4 photon valid | `truthPhotons.g4_photon_valid` |
| `truthpho_hepmc_association_valid` | vector<int32> |  | hepmc association valid | `truthPhotons.hepmc_association_valid` |
| `truthpho_truth_isolation_valid` | vector<int32> |  | truth isolation valid | `truthPhotons.truth_isolation_valid` |
| `truthpho_analysis_signal_r03` | vector<int32> |  | analysis signal r03 | `truthPhotons.analysis_signal_r03` |
| `truthpho_generator_barcode` | vector<int32> |  | generator barcode | `truthPhotons.generator_barcode` |
| `truthpho_truth_isolation_witness` | vector<double> |  | truth isolation witness | `truthPhotons.truth_isolation_witness` |
| `truthpho_truth_isolation_r03` | vector<double> |  | truth isolation r03 | `truthPhotons.truth_isolation_r03` |
| `truthpho_truth_isolation_r04` | vector<double> |  | truth isolation r04 | `truthPhotons.truth_isolation_r04` |
| `truthpho_reporting_guard_state` | vector<int32> |  | reporting guard state | `truthPhotons.reporting_guard_state` |
| `truthpho_gen_embedding_id` | vector<int32> |  | generator occurrence embedding id | `truthPhotons.generator_occurrence_embedding_id` |
| `truthpho_sample_source_role` | vector<int32> |  | sample source role | `truthPhotons.sample_source_role` |
| `truthpho_native_track_id` | vector<int32> |  | native track id | `truthPhotons.native_track_id` |
| `truthpho_embedding_id` | vector<int32> |  | embedding id | `truthPhotons.embedding_id` |
| `truthpho_native_vertex_id` | vector<int32> |  | native vertex id | `truthPhotons.native_vertex_id` |

### Truth jets (`truthjet<RR>_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truthjet<RR>_pt` | vector<double> | GeV | pt [GeV] | `truthJets.pt` |
| `truthjet<RR>_eta` | vector<double> | dimensionless (pseudorapidity) | eta [dimensionless (pseudorapidity)] | `truthJets.eta` |
| `truthjet<RR>_phi` | vector<double> | rad | phi [rad] | `truthJets.phi` |
| `truthjet<RR>_guard_state` | vector<int32> |  | reporting guard state | `truthJets.reporting_guard_state` |

### Truth photons without a reco candidate (`truthmiss_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `truthmiss_centrality_bin` | vector<int32> |  | centrality bin index | `truthPhotonMisses.centrality_bin_index` |
| `truthmiss_truth_pt` | vector<double> | GeV | truth pt [GeV] | `truthPhotonMisses.truth_pt` |

### Jet reco-truth links (`jlink_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `jlink_reco` | vector<int32> |  | reco id lo | `recoTruthLinks.reco_id_lo` |
| `jlink_truth` | vector<int32> |  | truth id lo | `recoTruthLinks.truth_id_lo` |
| `jlink_metric` | vector<double> | dimensionless (DeltaR) | match metric [dimensionless (DeltaR)] | `recoTruthLinks.match_metric` |
| `jlink_class` | vector<int32> |  | link class | `recoTruthLinks.link_class` |

### Photon reco-truth links (`plink_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `plink_reco` | vector<int32> |  | reco id lo | `recoTruthLinks.reco_id_lo` |
| `plink_truth` | vector<int32> |  | truth id lo | `recoTruthLinks.truth_id_lo` |
| `plink_metric` | vector<double> | dimensionless (DeltaR) | match metric [dimensionless (DeltaR)] | `recoTruthLinks.match_metric` |
| `plink_class` | vector<int32> |  | link class | `recoTruthLinks.link_class` |
| `plink_policy` | vector<int32> |  | photon association policy | `recoTruthLinks.photon_association_policy` |

### PPG12 diagnostic occurrences (`pdiag_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `pdiag_trigger` | vector<int32> |  | trigger | `RJPPG12DiagnosticOccurrenceV1.trigger` |
| `pdiag_kind` | vector<int32> |  | occurrence kind | `RJPPG12DiagnosticOccurrenceV1.occurrence_kind` |
| `pdiag_stage_code` | vector<int32> |  | stage code | `RJPPG12DiagnosticOccurrenceV1.stage_code` |
| `pdiag_source_code` | vector<int32> |  | source code | `RJPPG12DiagnosticOccurrenceV1.source_code` |
| `pdiag_sample_code` | vector<int32> |  | sample code | `RJPPG12DiagnosticOccurrenceV1.sample_code` |
| `pdiag_decision` | vector<int32> |  | decision | `RJPPG12DiagnosticOccurrenceV1.decision` |
| `pdiag_cluster_ordinal` | vector<int32> |  | cluster ordinal | `RJPPG12DiagnosticOccurrenceV1.cluster_ordinal` |
| `pdiag_truth_track_id` | vector<int32> |  | truth track id | `RJPPG12DiagnosticOccurrenceV1.truth_track_id` |
| `pdiag_generator_barcode` | vector<int32> |  | generator barcode | `RJPPG12DiagnosticOccurrenceV1.generator_barcode` |
| `pdiag_source_photon_pt` | vector<double> | GeV | source photon pt [GeV] | `RJPPG12DiagnosticOccurrenceV1.source_photon_pt` |
| `pdiag_flow_value` | vector<double> |  | flow value | `RJPPG12DiagnosticOccurrenceV1.flow_value` |
| `pdiag_reco_photon_pt` | vector<double> | GeV | reco photon pt [GeV] | `RJPPG12DiagnosticOccurrenceV1.reco_photon_pt` |
| `pdiag_response_photon_pt` | vector<double> | GeV | response photon pt [GeV] | `RJPPG12DiagnosticOccurrenceV1.response_photon_pt` |
| `pdiag_truth_photon_pt` | vector<double> | GeV | truth photon pt [GeV] | `RJPPG12DiagnosticOccurrenceV1.truth_photon_pt` |
| `pdiag_truth_prior_weight` | vector<double> | dimensionless weight | truth prior weight [dimensionless weight] | `RJPPG12DiagnosticOccurrenceV1.truth_prior_weight` |
| `pdiag_reco_vertex_z` | vector<double> | cm | reco vertex z [cm] | `RJPPG12DiagnosticOccurrenceV1.reco_vertex_z` |
| `pdiag_hard_truth_vertex_z` | vector<double> | cm | hard truth vertex z [cm] | `RJPPG12DiagnosticOccurrenceV1.hard_truth_vertex_z` |
| `pdiag_mb_truth_vertex_z` | vector<double> | cm | mb truth vertex z [cm] | `RJPPG12DiagnosticOccurrenceV1.mb_truth_vertex_z` |
| `pdiag_vertex_weight` | vector<double> | dimensionless weight | vertex weight [dimensionless weight] | `RJPPG12DiagnosticOccurrenceV1.vertex_weight` |
| `pdiag_period_event_weight` | vector<double> | dimensionless weight | period event weight [dimensionless weight] | `RJPPG12DiagnosticOccurrenceV1.period_event_weight` |
| `pdiag_cluster_energy` | vector<double> | GeV | cluster energy [GeV] | `RJPPG12DiagnosticOccurrenceV1.cluster_energy` |
| `pdiag_cluster_eta` | vector<double> | dimensionless (pseudorapidity) | cluster eta [dimensionless (pseudorapidity)] | `RJPPG12DiagnosticOccurrenceV1.cluster_eta` |
| `pdiag_cluster_et` | vector<double> | GeV | cluster et [GeV] | `RJPPG12DiagnosticOccurrenceV1.cluster_et` |
| `pdiag_truth_energy` | vector<double> | GeV | truth energy [GeV] | `RJPPG12DiagnosticOccurrenceV1.truth_energy` |
| `pdiag_truth_eta` | vector<double> | dimensionless (pseudorapidity) | truth eta [dimensionless (pseudorapidity)] | `RJPPG12DiagnosticOccurrenceV1.truth_eta` |
| `pdiag_truth_et` | vector<double> | GeV | truth et [GeV] | `RJPPG12DiagnosticOccurrenceV1.truth_et` |
| `pdiag_energy_contribution` | vector<double> |  | energy contribution | `RJPPG12DiagnosticOccurrenceV1.energy_contribution` |
| `pdiag_occurrence_weight` | vector<double> | dimensionless weight | occurrence weight [dimensionless weight] | `RJPPG12DiagnosticOccurrenceV1.occurrence_weight` |
| `pdiag_selection_bits` | vector<uint64> |  | selection bitmask | `RJPPG12DiagnosticOccurrenceV1.selection_bitmask` |

### Embedded-photon diagnostic occurrences (`embdiag_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `embdiag_kind` | vector<int32> |  | occurrence kind | `RJEmbeddedPhotonDiagnosticOccurrenceV1.occurrence_kind` |
| `embdiag_audit_bin` | vector<int32> |  | audit bin | `RJEmbeddedPhotonDiagnosticOccurrenceV1.audit_bin` |
| `embdiag_sample_code` | vector<int32> |  | sample code | `RJEmbeddedPhotonDiagnosticOccurrenceV1.sample_code` |
| `embdiag_centrality_bin_index` | vector<int32> |  | centrality bin index | `RJEmbeddedPhotonDiagnosticOccurrenceV1.centrality_bin_index` |
| `embdiag_photon_class` | vector<int32> |  | photon class | `RJEmbeddedPhotonDiagnosticOccurrenceV1.photon_class` |
| `embdiag_truth_pt` | vector<double> | GeV | truth pt [GeV] | `RJEmbeddedPhotonDiagnosticOccurrenceV1.truth_pt` |
| `embdiag_truth_isolation_witness` | vector<double> |  | truth isolation witness | `RJEmbeddedPhotonDiagnosticOccurrenceV1.truth_isolation_witness` |
| `embdiag_stitch_decision` | vector<int32> |  | stitch decision | `RJEmbeddedPhotonDiagnosticOccurrenceV1.stitch_decision` |
| `embdiag_stitch_photon_pt` | vector<double> | GeV | stitch photon pt [GeV] | `RJEmbeddedPhotonDiagnosticOccurrenceV1.stitch_photon_pt` |

### Event-display snapshots (empty in current production) (`snap_`)

| Branch | Type | Unit | Meaning | Source |
|---|---|---|---|---|
| `snap_selection_reason` | vector<int32> |  | selection reason | `RJEventDisplaySnapshotV1.selection_reason` |
| `snap_quota_class` | vector<int32> |  | quota class | `RJEventDisplaySnapshotV1.quota_class` |
| `snap_serialized_payload_hash` | vector<string> |  | serialized payload hash | `RJEventDisplaySnapshotV1.serialized_payload_hash` |

## Reserved branches

| Branch | Type | Filled by | Meaning |
|---|---|---|---|
| `pho_bdt_score` | vector<double> | canonical photon-ID augmentation | canonical BDT score (Au+Au: H70 model; pp: PPG12 model) |
| `pho_bdt_state` | vector<int32> | canonical photon-ID augmentation | score evaluation state (1 scored; other values explain why not) |
| `pho_id_class` | vector<int32> | canonical photon-ID augmentation | photon-ID region from the hash-bound working-point package: 3 tight, 2 excluded gap, 1 non-tight, 0 fail, -1 invalid |
| `pho_iso_class` | vector<int32> | canonical isolation package | isolation region from the hash-bound isolation package: 2 isolated, 1 gap, 0 non-isolated, -1 invalid |
| `jet<RR>_emcal_frac` | vector<double> | per-jet calorimeter fractions (next pass) | fraction of the jet energy in the EMCal |
| `jet<RR>_ihcal_frac` | vector<double> | per-jet calorimeter fractions (next pass) | fraction of the jet energy in the inner HCal |
| `jet<RR>_ohcal_frac` | vector<double> | per-jet calorimeter fractions (next pass) | fraction of the jet energy in the outer HCal |

## File-level objects

| Key | Content |
|---|---|
| `meta/schema` | schema name, version and sha256 of this file; product and layout profile; filled reserved groups |
| `meta/source` | the single sources row of the part, verbatim JSON |
| `meta/trigger_run_info` | the 64-row GL1 trigger table of the run segment (DATA), JSON with exact floats |
| `meta/centrality_binding` | centrality calibration binding (Au+Au DATA; required) |
| `meta/collections` | jet and truth-jet collections: name, input identity, algorithm, radius, quality bitmask |
| `meta/constants` | per-file constants hoisted from rows, verified constant by the converter |
| `meta/codes` | dictionaries for integer-coded strings |
| `meta/policies` | mbd_pmt_arm policy, mbd_pmt_charge_valid constant, invalid-JES policy string |
| `meta/source_tree_titles` | titles of the source trees, verbatim |
| `meta/replay/*` | replay-foundation named objects, verbatim |
| `meta/manifests/*` | augmentation, local-qualification, photon-repair and storage manifests, verbatim |
| `meta/converter_receipt` | converter version, input and output sha256, inverse-proof result |
| `hist/*` | producer histograms, byte copies of the source objects |

## Source coverage

All 809 columns of the production final layout are accounted for:

| Disposition | Columns |
|---|---:|
| DROP_DUPLICATE | 22 |
| DROP_IDENTITY | 120 |
| ENTRY | 3 |
| HOIST_META | 64 |
| OMITTED | 3 |
| STORE | 587 |
| STORE_CODE | 5 |
| STORE_INDEX | 5 |

The per-column list, with the rebuild rule of every column that is not stored as a branch, is the `sources` section of the schema file.
