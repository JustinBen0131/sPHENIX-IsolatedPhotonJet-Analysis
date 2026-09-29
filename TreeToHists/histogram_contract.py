"""Public sample histogram inventory; the C++ fillers own calculation and binning."""

# Public sample object inventory. Physics filling remains in the C++ reducers.
SHAPES = ['e11_over_e33', 'e32_over_e35', 'bdt', 'weta_cogx', 'wphi_cogx',
          'weta33_cogx', 'wphi33_cogx', 'et1', 'et2', 'et3', 'et4']
WIDTHS = ['weta_cogx', 'wphi_cogx', 'weta33_cogx', 'wphi33_cogx']
THRESHOLDS = [5, 7, 10, 12]
RESPONSE = {
    'photon_reco': 'reco_et', 'photon_truth': 'truth_et',
    'photon_response_truth_x_reco': 'truth_et,reco_et',
    'photon_fakes': 'reco_et', 'photon_misses': 'truth_et',
    'xj_reco': 'reco_et,xj', 'xj_truth': 'truth_et,xj',
    'xj_response_truth_global_x_reco_global': 'truth_cell,reco_cell',
    'xj_fakes': 'reco_et,xj', 'xj_fake_photon_reco': 'reco_et,xj',
    'xj_detector_fakes_reco': 'reco_et,xj', 'xj_misses': 'truth_et,xj',
    'xj_combinatoric_reco': 'reco_et,xj',
    'signal_leakage_truth_matched_abcd': 'region_abcd,reco_et',
    **{f'dphi_{kind}_{side}': 'reco_et,xj,dphi' for kind in ['truthjet', 'combinatoric']
       for side in ['near', 'far']}}

def inventory(system, family):
    auau, sim = system == 'AuAu', family != 'DATA'
    classes = ['cent0_20', 'cent20_50', 'cent50_80'] if auau else ['all']
    labels = ['data'] if not sim else (['sig', 'phoother'] if family == 'photonjet' else ['bkg', 'jetprompt'])
    F = 'data' if not sim else ('phosim' if family == 'photonjet' else 'jetsim')
    rows = []
    def add(path, axes, group, old=None, status='existing_booking'):
        rows.append(dict(path=path, axes=axes.split(','), group=group,
                         old_key=old, implementation=status, sumw2=True))
    for c in classes:
        for n, axes in [('events', 'count'), ('photons', 'et'), ('abcd', 'region_abcd,et')]:
            add(f'measurement/{c}/{n}', axes, 'measurement', f'{n}_{c}')
        for t in THRESHOLDS:
            add(f'measurement/{c}/jet{t}/recoil', 'region_ac,et,xj', 'measurement', f'recoil_{c}_jet{t}')
            for region in ['A', 'C']:
                for side in ['near', 'far']:
                    n = f'dphi{region}_{side}'
                    add(f'measurement/{c}/jet{t}/{n}', 'et,xj,dphi', 'measurement', f'dphi{region}_{c}_jet{t}_{side}')
            if sim:
                for n, axes in RESPONSE.items():
                    add(f'response/{c}/jet{t}/{n}', axes, 'response', f'{n}_{c}_jet{t}')
        for label in labels:
            dest = f'qa/photons/{c}/{label}'
            tag = f'{label}_{c}'
            for v in SHAPES:
                add(f'{dest}/shape_{v}', f'{v},et,stage', 'qa', f'ss_{v}_{tag}')
            for n, axes in [('ncb', 'npb,et'), ('iso', 'iso,et,id'),
                            ('abcd', 'region_abcd,et'), ('recoil', 'region_ac,et,xj'),
                            ('scoreiso', 'score,iso,et'), ('clu_et', 'wide_et,stage'),
                            ('clu_etaphi', 'eta,phi,stage')]:
                add(f'{dest}/{n}', axes, 'qa', f'{n}_{tag}')
            if auau:
                add(f'{dest}/isonoue', 'iso,et,id', 'qa', f'isonoue_{tag}')
            if sim:
                add(f'{dest}/abcdown', 'region_abcd,et,generator_window', 'qa', f'abcdown_{tag}')
                add(f'{dest}/recoilown', 'region_ownership,et,xj', 'qa', f'recoilown_{tag}')
        for n, axes in [('evt_vz', 'vz'), ('jet_pt', 'jet_pt,event_selection'),
                        ('jet_etaphi', 'jet_eta,jet_phi,event_selection')]:
            add(f'qa/events/{c}/{n}', axes, 'qa', f'{n}_{F}_{c}')
        if not sim:
            add(f'qa/photons/{c}/data/e11e33_dt', 'e11_over_e33,et,time_minus_mbd', 'qa', f'e11e33_dt_data_{c}')
            add(f'qa/photons/{c}/data/etatime', 'eta,time_minus_mbd', 'qa', f'etatime_data_{c}')
        if family == 'photonjet':
            add(f'qa/truth/{c}/truth_iso', 'truth_iso,truth_et,prompt_class', 'qa', f'truth_iso_{c}')
            add(f'qa/truth/{c}/eff', 'truth_et,efficiency_stage,reco_et', 'qa', f'eff_{c}')
        # Three fraction distributions cover all nine stage-3 branches by component axis.
        for radius in ['02', '03', '04']:
            add(f'qa/jets/{c}/R{radius}/calo_fractions', 'fraction,jet_pt,calo_component',
                'extension', status='proposed_requires_stage3')
        if sim:
            for obj in ['photon', 'jet']:
                add(f'qa/matching/{c}/{obj}', 'delta_r,reco_pt,match_class',
                    'extension', status='proposed')
        else:
            add(f'qa/photons/{c}/data/e11e33_npb_timing', 'e11_over_e33,et,time_minus_mbd',
                'extension', status='proposed_requires_explicit_tag_predicate')
    for label in labels:
        dest = f'qa/slices/{label}'
        for n, axes in [('bdt', 'score,et,centrality_slice'), ('isoslice', 'iso,et,centrality_slice')]:
            add(f'{dest}/{n}', axes, 'qa', f'{n}_{label}')
        if auau:
            add(f'{dest}/isonoueslice', 'iso,et,centrality_slice', 'qa', f'isonoueslice_{label}')
        for v in WIDTHS:
            add(f'{dest}/width_{v}', f'{v},et,centrality_slice', 'qa', f'width_{v}_{label}')
    add('qa/events/calo_total', 'calo_total,centrality,mb_state', 'qa', f'evt_calo_{F}')
    if auau:
        add('qa/events/centrality', 'centrality,weight_mode', 'qa', f'evt_cent_{F}')
        if not sim:
            add('qa/events/centrality_scaled_bits', 'centrality,scaled_bit', 'qa', 'cent_bits_data')
    if sim:
        add('qa/generator/leading_pt', 'generator_pt,generator_kind', 'qa', f'gen_{F}')
        add('qa/photon_id/score_roles', 'score,et,centrality_slice,truth_role',
            'extension', status='proposed_no_truth_isolation_requirement')
    else:
        add('qa/events/trigger_bits_raw', 'max_cluster_energy,scaled_bit,live_bit', 'qa', 'trig_data')
        add('qa/events/trigger_max_cluster_energy', 'max_cluster_energy,semantic_trigger,trigger_state',
            'extension', status='proposed_requires_per_run_trigger_blocks')
    for n, axes in [('cutflow', 'gate,centrality,weight_mode'), ('weights', 'event_weight,weight_state'),
                    ('mbd_charge', 'mbd_charge,centrality')]:
        add(f'qa/events/{n}', axes, 'extension', status='proposed')
    for detector in ['emcal', 'ihcal', 'ohcal', 'total']:
        add(f'qa/events/{detector}_vs_mbd_charge', f'calo_{detector},mbd_charge,mb_state',
            'extension', status='proposed')
        if detector != 'total':
            add(f'qa/events/{detector}_vs_centrality', f'calo_{detector},centrality,mb_state',
                'extension', status='proposed')
    if auau and not sim:
        # Retain the nominal Photon12 stream and its top-level aliases.
        measurement = [r for r in rows if r['group']=='measurement']
        for trigger in ('photon12',):
            rows.extend(dict(r,path='triggers/'+trigger+'/'+r['path'],
                             old_key='triggers/'+trigger+'/'+r['old_key']) for r in measurement)
    assert len({r['path'] for r in rows}) == len(rows)
    return rows
