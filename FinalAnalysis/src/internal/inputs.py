"""The bridge from TreeToHists sample packages to FinalAnalysis/src/analysis.py:

    event_table_arrays   the event tables of FillPhotonJet outputs in the layout of the whole-event bootstrap
                         (FinalAnalysis event_bootstrap.EventTable), on any photon-pT and xJ binning: they keep exact
                         photon E_T and xJ
    simulation_sums      one FillPhotonJetResponse output summed onto the measurement's response axes, named as
                         FinalAnalysis closure.package_window takes them
    simulation_dphi      the response's Delta-phi fills on those axes (the simulation closure of the data-driven
                         combinatoric estimate)
    sample_weight        the per-event normalisation of a simulated sample (FinalAnalysis/config/samples.yaml)
"""
from __future__ import annotations

import json
from pathlib import Path
import hashlib

import numpy as np
import uproot

THRESHOLDS = (5.0, 7.0, 10.0, 12.0)
NATIVE_PT_EDGES = np.arange(15.0, 36.0)


CENTRALITY_SLICE_LABELS = tuple("%d_%d" % (5 * i, 5 * i + 5) for i in range(16))


def event_table_arrays(paths, xj_edges, all_jets: bool = False, *, data_trigger=None, runs=None) -> dict[str, np.ndarray]:
    """The event tables of FillPhotonJet outputs (one path or several, concatenated in order) in the event-contribution
    layout FinalAnalysis/event_bootstrap.py reads: photon E_T bins on the native 1 GeV grid, recoil jets (window 1 only)
    binned on the measurement xJ edges, ROOT bin numbers, and jet__thresholds_passed as the number of thresholds a jet
    passes (they are nested). With ``all_jets`` the jets at every Delta phi are kept, with jet__window (1 for a recoil
    jet), jet__dphi and jet__far (1 for |Delta eta| > 0.8), for the data-driven combinatoric estimate.
    AuAu DATA defaults to its saved nominal trigger. Set data_trigger='photon10'
    or 'photon12' to read a retained stream, optionally restricted to runs.
    Event references are remapped after filtering; an event firing both triggers
    appears once in either choice. Legacy files need no new trigger argument."""
    paths = [paths] if isinstance(paths, (str, Path)) else list(paths)
    compatible_outputs(paths)
    xj_edges = np.asarray(xj_edges, dtype=float)
    parts = []
    for path in paths:
        with uproot.open(path) as f:
            keys = ("records/events","records/photons","records/jets") if 'meta/package' in f else ("table_events","table_photons","table_jets")
            if 'meta/package' in f:
                meta=json.loads(str(f['meta/package']))
                if meta['system']=='AuAu' and meta['family']=='data' and (data_trigger or meta.get('data_trigger')):
                    chosen=data_trigger or meta['data_trigger']
                    selected,_=trigger_sources(f,chosen,runs)
                    keys=tuple('records/triggers/'+chosen+'/'+k for k in ('events','photons','jets'))
                    if any(k not in f for k in keys): raise ValueError('package lacks retained trigger event tables')
                elif data_trigger is not None:
                    raise ValueError('named trigger switching requires a retained AuAu DATA package')
                else: selected=None
            else:
                if data_trigger is not None: raise ValueError('legacy file has no retained trigger streams')
                selected=None
            ev,ph,jt=(f[k].arrays(library='np') for k in keys)
            keep=np.ones(len(ev['run']),dtype=bool)
            if selected is not None: keep &= np.isin(ev['source'],selected)
            if runs is not None: keep &= np.isin(ev['run'],tuple(runs))
            remap=np.full(len(keep),-1,dtype=np.int64);remap[keep]=np.arange(keep.sum())
            ev={k:v[keep] for k,v in ev.items()}
            def objects(table):
                mask=keep[table['event']]
                out={k:v[mask] for k,v in table.items()}
                out['event']=remap[out['event']];return out
            parts.append((ev,objects(ph),objects(jt)))
    offsets = np.cumsum([0] + [len(ev["run"]) for ev, _, _ in parts])
    ev = {k: np.concatenate([p[0][k] for p in parts]) for k in ("run", "seq", "centrality", "weight", "slice")}
    ph = {k: np.concatenate([p[1][k] for p in parts]) for k in ("region", "et_bin")}
    ph["event"] = np.concatenate([p[1]["event"] + o for p, o in zip(parts, offsets)])
    if all_jets and not all("far" in p[2] for p in parts):
        raise ValueError("the jets at every Delta phi need FillPhotonJet outputs with the far column")
    recoil = [np.ones(len(p[2]["window"]), dtype=bool) if all_jets else p[2]["window"] == 1 for p in parts]
    keys = ("region", "et_bin", "xj", "pt", "thresholds") + (("window", "dphi", "far") if all_jets else ())
    jt = {k: np.concatenate([p[2][k][r] for p, r in zip(parts, recoil)]) for k in keys}
    jt["event"] = np.concatenate([p[2]["event"][r] + o for p, r, o in zip(parts, recoil, offsets)])
    xj_bin = np.searchsorted(xj_edges, jt["xj"], side="right")  # ROOT bin number; 0 underflow, n+1 overflow
    extra = {"jet__window": jt["window"].astype(np.int8), "jet__dphi": jt["dphi"], "jet__far": jt["far"].astype(np.int8)} if all_jets else {}
    return {**extra, "pt_edges_gev": NATIVE_PT_EDGES.copy(), "xj_edges": xj_edges, "jet_thresholds_gev": np.asarray(THRESHOLDS),
            "cent_slice_labels": np.asarray(CENTRALITY_SLICE_LABELS),
            "event__run": ev["run"], "event__seq": ev["seq"], "event__centrality": ev["centrality"], "event__weight": ev["weight"],
            "event__cent_slice": ev["slice"].astype(np.int8),
            "photon__event_index": ph["event"].astype(np.int64), "photon__region_ABCD": ph["region"].astype(np.int8),
            "photon__pt_bin": ph["et_bin"].astype(np.int8),
            "jet__event_index": jt["event"].astype(np.int64), "jet__region_A0_C1": jt["region"].astype(np.int8),
            "jet__pt_bin": jt["et_bin"].astype(np.int8), "jet__xj_bin": xj_bin.astype(np.int16),
            "jet__thresholds_passed": np.array([int(m).bit_length() for m in jt["thresholds"]], dtype=np.int8),
            "jet__pt": jt["pt"], "jet__xj": jt["xj"]}


def trigger_sources(f, data_trigger, runs=None):
    """Enabled original sources and exposure records for one exact AuAu trigger.

    Run membership is separate from an event firing the scaled trigger. Disabled
    and absent triggers contribute zero even if an inconsistent event bit is set.
    """
    if data_trigger not in AUAU_DATA_TRIGGERS: raise ValueError('unknown AuAu trigger')
    meta=json.loads(str(f['meta/package']))
    if meta['system']!='AuAu' or meta['family']!='data': raise ValueError('expected AuAu DATA package')
    records=[];selected=[]
    requested=None if runs is None else set(runs)
    raw_rows=(raw for chunk in f['records/sources'].iterate(['json'],step_size='8 MB',library='np') for raw in chunk['json'])
    for i,raw in enumerate(raw_rows):
        source=json.loads(raw);m=source['metadata'];table=m.get('trigger_selection') or m.get('trigger_run_info')
        if not table or len(table['name'])!=64 or len(table['initial_prescale'])!=64:
            raise ValueError('missing or malformed per-source trigger metadata')
        bits=[b for b,n in enumerate(table['name']) if n==AUAU_DATA_TRIGGERS[data_trigger]]
        if len(bits)>1: raise ValueError('ambiguous exact trigger name')
        run=m.get('source',{}).get('run')
        if requested is not None and run is None: raise ValueError('missing source run identity')
        enabled=bool(bits) and int(table['initial_prescale'][bits[0]])>=0
        included=enabled and (requested is None or run in requested)
        records.append(dict(source=i,run=run,bit=bits[0] if bits else None,enabled=enabled,included=included,
                            initial_prescale=table['initial_prescale'][bits[0]] if bits else None,input_sha256=source['sha256']))
        if included:selected.append(i)
    return selected,records


def measurement_histograms(path, *, data_trigger=None, runs=None):
    """Saved AuAu trigger measurements restricted by original run/source.

    Returns canonical measurement paths mapped to (values, variances), plus the
    selected/excluded source records. No eventTree or ROOT filling is required.
    """
    with uproot.open(path) as f:
        meta=json.loads(str(f['meta/package']))
        chosen=data_trigger or meta.get('data_trigger')
        selected,records=trigger_sources(f,chosen,runs)
        prefix='triggers/'+chosen+'/'
        keys=[r['path'] for r in json.loads(str(f['meta/inventory'])) if r['path'].startswith(prefix+'measurement/')]
        if not keys or meta.get('trigger_retention')!='named_scaled_trigger_measurements_and_source_contributions':
            raise ValueError('package lacks retained trigger measurement contributions')
        sums={k:(np.zeros(f[k].values().shape),np.zeros(f[k].values().shape)) for k in keys}
        for source in selected:
            for k,pair in source_histograms(f,keys,source).items():
                for a,b in zip(sums[k],pair): a+=b
        return {k[len(prefix):]:v for k,v in sums.items()},records


def sample_weight(sample: dict, reference_cross_section_pb: dict) -> float:
    """Per-event normalisation of a simulated sample, in pb per generated event, from its FinalAnalysis/config/samples.yaml entry:
    cross_section_over_generated  cross_section_pb / generated_events;
    family_reference              reference_cross_section_pb[reference] / generated_events (the producer event weights
                                  carry the samples' relative cross sections)."""
    rule, generated = sample.get("source_factor"), sample.get("generated_events")
    if not generated or generated <= 0:
        raise ValueError("a simulated sample needs its generated event count")
    if rule == "cross_section_over_generated":
        numerator = sample.get("cross_section_pb")
    elif rule == "family_reference":
        numerator = reference_cross_section_pb.get(sample.get("reference"))
    else:
        raise ValueError("unsupported sample normalisation: %s" % rule)
    if not numerator or numerator <= 0:
        raise ValueError("a simulated sample needs a positive cross section")
    return float(numerator) / float(generated)


# Response histogram -> numerical input (see internal/closure.py package_window).
RESPONSE_NAMES = {
    "photon_response_truth_x_reco": "photon_response", "photon_misses": "photon_misses", "photon_truth": "truth_photons",
    "photon_reco": "photon_reco", "photon_fakes": "photon_fakes",
    "xj_response_truth_global_x_reco_global": "pair_response", "xj_misses": "pair_misses", "xj_truth": "pair_truth",
    "xj_reco": "pair_reco", "xj_fake_photon_reco": "pair_fakes_UNMATCHED_RECO", "xj_combinatoric_reco": "pair_fakes_COMBINATORIC",
    "xj_detector_fakes_reco": "pair_boundary_reco_ASSIGNED_HARD_NONFIDUCIAL", "signal_leakage_truth_matched_abcd": "leakage_abcd_reco"}
SUMW2_NAMES = {"xj_combinatoric_reco": "pair_fakes_COMBINATORIC_sumw2", "signal_leakage_truth_matched_abcd": "leakage_abcd_reco_sumw2"}


def aggregation(fine: np.ndarray, coarse: np.ndarray) -> np.ndarray:
    """Coarse bin of every fine bin. The coarse edges must be fine edges and span the same range: the response
    categories depend on the axis range, so only then is the sum exact."""
    fine, coarse = np.asarray(fine, float), np.asarray(coarse, float)
    if fine[0] != coarse[0] or fine[-1] != coarse[-1] or not np.all(np.isin(coarse, fine)):
        raise ValueError("coarse edges must be fine edges with the same range")
    return np.searchsorted(coarse, fine[:-1], side="right") - 1


def name_suffix(class_name: str, threshold: float) -> str:
    """The name suffix of the per-class, per-threshold FillPhotonJetResponse histograms: _<class>_jet<threshold>."""
    return "_%s_jet%d" % (class_name, int(threshold))


def source_views(path):
    """Stable logical-source views for source bootstrap, without unpacking files."""
    with uproot.open(path) as f:
        if 'meta/package' not in f: return [path]
        count = f['records/sources'].num_entries
    return [str(path)+'::source='+str(i) for i in range(count)]


def _source_path(path):
    text=str(path)
    if '::source=' not in text: return path,None
    base,index=text.rsplit('::source=',1)
    return base,int(index)


def _response_histograms(f,names,class_name,threshold,source=None):
    """Read aggregate histograms or reconstruct one source's exact sparse cells."""
    packaged='meta/package' in f
    keys={n:'response/%s/jet%d/%s'%(class_name,int(threshold),n) if packaged else n+name_suffix(class_name,threshold) for n in names}
    values=source_histograms(f,list(keys.values()),source)
    return {n:values[k] for n,k in keys.items()}


def source_histograms(f,keys,source=None):
    """Read named aggregate histograms or exact saved source contributions.

    Source slices are indexed in the source table, so reading all sources is
    linear in the stored contributions rather than scanning the file per source.
    """
    if source is None:
        return {k:(f[k].values(),f[k].variances()) for k in keys}
    packaged='meta/package' in f
    if not packaged or not 0<=source<f['records/sources'].num_entries:
        raise ValueError('source view is outside this package')
    arrays={k:(np.zeros(f[k].values(flow=True).shape),np.zeros(f[k].values(flow=True).shape)) for k in keys}
    records=f['records/sources']
    limits={}
    prefix='contribution' if 'contribution_begin' in records else 'response'
    if prefix+'_begin' in records:
        row=records.arrays([prefix+'_begin',prefix+'_end'],entry_start=source,entry_stop=source+1,library='np')
        limits={'entry_start':int(row[prefix+'_begin'][0]),'entry_stop':int(row[prefix+'_end'][0])}
    table='records/source_contributions' if 'records/source_contributions' in f else 'records/response_contributions'
    for chunk in f[table].iterate(step_size='16 MB',library='np',**limits):
        choose=chunk['source']==source
        for k,(v,w) in arrays.items():
            mask=choose&(chunk['histogram']==k)
            cells=chunk['bin'][mask]
            if not len(cells): continue
            coords=np.unravel_index(cells,v.shape,order='F')
            np.add.at(v,coords,chunk['value'][mask]);np.add.at(w,coords,chunk['variance'][mask])
    return {k:tuple(a[tuple(slice(1,-1) for _ in a.shape)] for a in pair) for k,pair in arrays.items()}


def simulation_sums(path: Path, class_name: str, threshold: float, reco_pt_edges, truth_pt_edges, xj_edges) -> dict[str, np.ndarray]:
    """The response inputs of one FillPhotonJetResponse output for one class and threshold, summed onto the given
    axes and named as FinalAnalysis takes them (pair arrays flattened as pT bin * xJ bins + xJ bin)."""
    path,source=_source_path(path)
    compatible_outputs([path])
    suffix = name_suffix(class_name, threshold)
    with uproot.open(path) as f:
        fine = json.loads(str(f["grids"]))
        hists = _response_histograms(f,RESPONSE_NAMES,class_name,threshold,source)
    reco, truth, xj = aggregation(fine["reco_pt"], reco_pt_edges), aggregation(fine["truth_pt"], truth_pt_edges), aggregation(fine["xj"], xj_edges)
    nr, nt, nx = len(reco_pt_edges) - 1, len(truth_pt_edges) - 1, len(xj_edges) - 1
    cell_reco = (reco[:, None] * nx + xj[None, :]).ravel()
    cell_truth = (truth[:, None] * nx + xj[None, :]).ravel()
    out = {}
    for name, (values, variances) in hists.items():
        for target, array in ((RESPONSE_NAMES[name], values), (SUMW2_NAMES.get(name), variances)):
            if target is None:
                continue
            if name in ("photon_reco", "photon_fakes"):
                summed = np.zeros(nr); np.add.at(summed, reco, array)
            elif name in ("photon_truth", "photon_misses"):
                summed = np.zeros(nt); np.add.at(summed, truth, array)
            elif name == "photon_response_truth_x_reco":
                summed = np.zeros((nt, nr)); np.add.at(summed, (truth[:, None], reco[None, :]), array)
            elif name in ("xj_reco", "xj_fake_photon_reco", "xj_combinatoric_reco", "xj_detector_fakes_reco"):
                summed = np.zeros(nr * nx); np.add.at(summed, cell_reco, array.ravel())
            elif name in ("xj_truth", "xj_misses"):
                summed = np.zeros(nt * nx); np.add.at(summed, cell_truth, array.ravel())
            elif name == "xj_response_truth_global_x_reco_global":
                summed = np.zeros((nt * nx, nr * nx)); np.add.at(summed, (cell_truth[:, None], cell_reco[None, :]), array)
            else:  # leakage: region x reco pT
                summed = np.zeros((4, nr)); np.add.at(summed, (slice(None), reco), array)
            out[target] = summed
    return out


DPHI_FILLS = ("dphi_truthjet", "dphi_combinatoric")


def simulation_dphi(path: Path, class_name: str, threshold: float, reco_pt_edges, xj_edges) -> dict[str, np.ndarray] | None:
    """The Delta-phi fills of one FillPhotonJetResponse output (jets of matched photons with and without a truth jet
    nearby) for one class and threshold, summed onto the given axes: (reco pT, xJ, Delta phi, side 0 near / 1 far), with
    their sumw2. None for an output written before these fills existed."""
    path,source=_source_path(path)
    compatible_outputs([path])
    suffix = name_suffix(class_name, threshold)
    with uproot.open(path) as f:
        keys = {name.split(";")[0] for name in f.keys()}
        if 'meta/package' not in f and not all("%s_%s%s" % (name, side, suffix) in keys for name in DPHI_FILLS for side in ("near", "far")):
            return None
        fine = json.loads(str(f["grids"]))
        names=[name+'_'+side for name in DPHI_FILLS for side in ('near','far')]
        values=_response_histograms(f,names,class_name,threshold,source)
        hists = {(name,s):values[name+'_'+side] for name in DPHI_FILLS for s,side in enumerate(('near','far'))}
    reco, xj = aggregation(fine["reco_pt"], reco_pt_edges), aggregation(fine["xj"], xj_edges)
    n_dphi = next(iter(hists.values()))[0].shape[2]
    out = {}
    for name in DPHI_FILLS:
        for key, index in ((name, 0), (name + "_sumw2", 1)):
            summed = np.zeros((len(reco_pt_edges) - 1, len(xj_edges) - 1, n_dphi, 2))
            for s in (0, 1):
                np.add.at(summed[..., s], (reco[:, None], xj[None, :]), hists[(name, s)][index])
            out[key] = summed
    return out


AUAU_DATA_TRIGGERS = {
    'photon10': 'Photon 10 GeV + MBD NS >= 2, vtx < 150 cm',
    'photon12': 'Photon 12 GeV + MBD NS >= 2, vtx < 150 cm',
}


def compatible_outputs(paths):
    """Verify each saved selection identity and require a common profile."""
    profiles = set()
    for path in paths:
        path, _ = _source_path(path)
        with uproot.open(path) as f:
            profile = str(f['profile']) if 'profile' in f else 'legacy_unbound'
            if 'selection' in f:
                selection = json.loads(str(f['selection']))
                body = {k: v for k, v in selection.items() if k not in ('selection_sha256', 'profile', 'name')}
                pin = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
                if selection['selection_sha256'] != pin or selection['profile'] != 'selection:' + pin or profile != selection['profile']:
                    raise ValueError('saved selection digest/profile differs: ' + str(path))
            elif profile.startswith('selection:'):
                raise ValueError('missing resolved selection metadata: ' + str(path))
            profiles.add(profile)
    if len(profiles) > 1: raise ValueError('cannot combine different selections')


def source_inputs(paths):
    """One response view per original source, independent of physical file count."""
    return [view for path in paths for view in source_views(path)]


def response_identity(view, *, profile, centrality_weights, normalization):
    """Check a source's saved recipe and return its stable bootstrap identity."""
    path, index = _source_path(view)
    with uproot.open(path) as f:
        if str(f['profile']) != profile: raise ValueError('response selection differs: ' + str(path))
        if 'meta/package' in f:
            meta = json.loads(str(f['meta/package']))
            weights = meta['centrality_weights']; n = meta['normalization']
            if n.get('applied_in_histograms') is not False:
                raise ValueError('response normalization is already applied or unspecified')
            factor = sample_weight(normalization['sample'], normalization['reference'])
            if 'weight_pb_per_event' in n and n['weight_pb_per_event'] != factor:
                raise ValueError('saved normalization differs from analysis configuration')
            if 'numerator_pb' in n and n['numerator_pb'] / n['generated_events'] != factor:
                raise ValueError('saved normalization ingredients differ from analysis configuration')
            if index is None: raise ValueError('packaged response needs one source view per bootstrap row')
            row = json.loads(f['records/sources'].arrays(['json'], entry_start=index, entry_stop=index+1,
                                                        library='np')['json'][0])
            source = row['metadata']['source']
            identity = (source.get('input_file_sha256') or row.get('sha256'),
                        str(source.get('logical_source_id') or row.get('identity') or row.get('path')))
        else:
            weights = json.loads(str(f['centrality_weights']))['per_slice']
            sources = json.loads(str(f['sources']))
            if len(sources) != 1: raise ValueError('legacy response must contain exactly one source')
            identity = sources[0]['input_file_sha256'], sources[0]['path']
        if list(map(float, weights)) != list(map(float, centrality_weights)):
            raise ValueError('response centrality weights differ from configuration')
    return identity


def read_collection(manifest_path):
    """Read the public histogram manifest, verify all files and pin their metadata."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema') != 'PhotonJetHistogramCollectionV1' or not manifest.get('samples'):
        raise ValueError('expected a nonempty TreeToHists manifest')
    output = {}; paths = set()
    for key, item in manifest['samples'].items():
        path = (manifest_path.parent / item['path']).resolve()
        if not path.is_relative_to(manifest_path.parent) or path in paths:
            raise ValueError('histogram manifest contains an escaped or duplicate path')
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(8 << 20), b''): h.update(chunk)
        if h.hexdigest() != item['sha256']: raise ValueError('histogram file changed: ' + str(path))
        with uproot.open(path) as f:
            m = json.loads(str(f['meta/package']))
            if m.get('schema') != 'PhotonJetSampleV1' or m['plan_sha256'] != manifest['plan_sha256']:
                raise ValueError('sample is not bound to this histogram collection')
            if m.get('system') not in ('PP', 'AuAu') or m.get('family') not in ('data', 'photonjet', 'inclusivejet'):
                raise ValueError('invalid sample system or family')
            if key != m['system'] + '/' + m['name'] or (m['name'] == 'data') != (m['family'] == 'data'):
                raise ValueError('manifest key and sample identity differ: ' + key)
            if f['records/sources'].num_entries != item['sources']:
                raise ValueError('sample source count differs from its manifest')
            if m.get('source_count') != item['sources']:
                raise ValueError('sample metadata source count differs from its manifest')
            compatible_outputs([path])
        paths.add(path); output[key] = dict(path=path, sha256=item['sha256'], metadata=m)
    if len({item['metadata']['profile'] for item in output.values()}) != 1:
        raise ValueError('histogram collection mixes selection profiles')
    return output
