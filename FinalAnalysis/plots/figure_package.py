"""IAN figure package from the canonical chain: one npz per result, built from the objects FinalAnalysis
measurement.measure(record=...) records while it computes the result. No correction is reimplemented here; the package
only arranges recorded arrays and their toy covariances, and it refuses to write unless the recorded unfoldings give the
result's per-photon values exactly.

Layout: photon-pT-major, (n_pt, n_xj) arrays; flattened vectors are pT-major (pT bin outer, xJ inner), as FinalAnalysis
flattens them. Covariances come from the accepted data toys of the result (the statistical spread of each quantity);
the headline per-photon uncertainty is the result's own sigma_total (data toys plus response toys). The package holds
what FinalAnalysis/plots/render.py draws; the result document holds everything else.
"""
import json

import numpy as np


def _cov(rows):
    rows = np.asarray(rows, dtype=float)
    return np.cov(rows, rowvar=False, ddof=1)


def _cross(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return (a - a.mean(axis=0)).T @ (b - b.mean(axis=0)) / (len(a) - 1)


def per_interval(pairs, photons, xj_edges):
    """Per-photon xJ density in each photon-pT interval: pairs of the bin / photons of the bin / xJ width."""
    width = np.diff(np.asarray(xj_edges, dtype=float))
    n_pt = len(photons)
    return np.asarray(pairs, dtype=float).reshape(n_pt, -1) / np.asarray(photons, dtype=float)[:, None] / width[None, :]


def _same(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return bool(a.shape == b.shape and np.allclose(a, b, rtol=1e-10, atol=1e-12 * max(1.0, float(np.max(np.abs(b))) if b.size else 1.0)))


def xj_projection(native_edges, edges):
    """The analysis xJ bin of every native xJ bin as an (n, n_native) 0/1 matrix; the analysis edges must be native
    edges over the same range, so that the sum is exact."""
    native_edges, edges = np.asarray(native_edges, dtype=float), np.asarray(edges, dtype=float)
    if not (np.isclose(native_edges[0], edges[0]) and np.isclose(native_edges[-1], edges[-1])
            and all(np.any(np.isclose(native_edges, e, rtol=0, atol=1e-9)) for e in edges)):
        raise ValueError('the analysis xJ edges are not native xJ edges over the same range')
    index = np.searchsorted(edges, native_edges[:-1] + 1e-9, side='right') - 1
    return (index[None, :] == np.arange(len(edges) - 1)[:, None]).astype(float)


def build(result, record, label, per_photon, native=None, restrict=None):
    """The package of one result. ``per_photon`` is FinalAnalysis prior.per_photon (the result's own definition).
    ``native`` (FinalAnalysis/plots/native_response.py) adds the response on its native grid for the response figures, after checking
    with ``restrict`` (FinalAnalysis unfolding.restrict) that it sums exactly onto the measurement window."""
    s, w, c, data = record['settings'], record['window'], record['central'], record['data']
    xj_edges, pt_edges = np.asarray(s.xj_edges, dtype=float), np.asarray(s.pt_edges, dtype=float)
    n_pt, n_xj = len(pt_edges) - 1, len(xj_edges) - 1
    entry = result['entry']
    toy_pairs, toy_photons = np.asarray(record['toy_corrected_pairs']), np.asarray(record['toy_corrected_photons'])
    toy_upairs, toy_uphotons = np.asarray(record['toy_unfolded_pairs']), np.asarray(record['toy_unfolded_photons'])
    checks = {
        'per_photon_equals_result': bool(np.array_equal(per_photon(record['unfolded_pairs'], record['unfolded_photons'], xj_edges),
                                                        np.asarray(entry['central']))),
        'corrected_equals_result': bool(np.array_equal(c['corrected'], np.asarray(result['corrected_ptgamma_xj']))),
        'photons_equal_result': bool(np.array_equal(c['photons'], np.asarray(result['purity_corrected_photons_ptgamma']))),
        'toys': int(len(toy_upairs)),
    }
    if not (checks['per_photon_equals_result'] and checks['corrected_equals_result'] and checks['photons_equal_result']):
        raise ValueError('the recorded chain does not reproduce the result: %s' % checks)
    interval = per_interval(record['unfolded_pairs'], record['unfolded_photons'], xj_edges)
    toy_interval = np.array([per_interval(p, g, xj_edges).ravel() for p, g in zip(toy_upairs, toy_uphotons)])
    reco_interval = per_interval(c['corrected'].ravel(), c['photons'], xj_edges)
    toy_reco_interval = np.array([per_interval(p, g, xj_edges).ravel() for p, g in zip(toy_pairs, toy_photons)])
    package = dict(
        label=np.array(label), pt_edges=pt_edges, xj_edges=xj_edges,
        reported_mask=np.asarray(result['reported_xj_mask'], dtype=bool), jet_pt_min_gev=np.array(s.jet_pt_min_gev),
        pair_iterations=np.array(record['pair_iterations']), photon_iterations=np.array(record['photon_iterations']),
        # region A and the corrections (the purity and combinatoric figures)
        region_a=np.asarray(data['recoil_spectra'][0], dtype=float), region_a_sumw2=np.asarray(data['recoil_spectra_sumw2'][0], dtype=float),
        scaled_region_c=c['scaled_region_c'], purity_corrected=c['purity_corrected'],
        purity_corrected_variance=c['purity_corrected_variance'], template_combinatoric=c['combinatoric_subtracted'],
        corrected=c['corrected'], corrected_cov=_cov(toy_pairs),
        # photon normalisation
        photons_corrected=c['photons'], photons_corrected_cov=_cov(toy_photons),
        photons_unfolded=record['unfolded_photons'], photons_unfolded_cov=_cov(toy_uphotons),
        photons_unfolded_corrected_cross=_cross(toy_uphotons, toy_photons),
        # the joint unfolding refolded (the refold figure)
        pairs_refolded=record['refolded_pairs'],
        # per photon: the result itself, and per photon-pT interval from the same unfolding
        per_photon=np.asarray(entry['central'], dtype=float), sigma_total=np.asarray(entry['sigma_total'], dtype=float),
        per_photon_interval=interval, per_photon_interval_cov=_cov(toy_interval),
        reco_per_photon_interval=reco_interval, reco_per_photon_interval_cov=_cov(toy_reco_interval),
        # responses (the measurement window)
        xj_response=w.xj_response, xj_misses=w.xj_misses, photon_response=w.photon_response, photon_misses=w.photon_misses,
    )
    if native is not None:
        view = restrict(native['package'])
        checks['native_photon_restricts_to_window'] = _same(view.photon_response, w.photon_response) and _same(view.photon_misses, w.photon_misses)
        # restrict aggregates photon pT; the native xJ grid is summed onto the analysis xJ edges here
        k = np.kron(np.eye(n_pt), xj_projection(native['package'].grids.xj_edges, w.xj_edges))
        checks['native_pairs_restrict_to_window'] = (_same(k @ view.xj_response @ k.T, w.xj_response)
                                                     and _same(k @ np.asarray(view.xj_misses).ravel(), np.asarray(w.xj_misses).ravel()))
        if not (checks['native_photon_restricts_to_window'] and checks['native_pairs_restrict_to_window']):
            raise ValueError('the native response does not sum onto the measurement window: %s' % checks)
        g, a = native['package'].grids, native['package'].arrays
        package.update(native_source=np.array(native['source']), native_truth_pt_edges=np.asarray(g.truth_pt_edges, dtype=float),
                       native_reco_pt_edges=np.asarray(g.reco_pt_edges, dtype=float), native_xj_edges=np.asarray(g.xj_edges, dtype=float),
                       native_photon_response=np.asarray(a['photon_response'], dtype=float),
                       native_photon_misses=np.asarray(a['photon_misses'], dtype=float),
                       native_pair_response=np.asarray(a['pair_response'], dtype=float),
                       native_pair_misses=np.asarray(native['package'].all_pair_misses(), dtype=float))
        for name in ('photon_neff', 'pair_neff'):
            if native.get(name) is not None:
                package['native_' + name] = np.asarray(native[name], dtype=float)
    for k, values in (record.get('unfolded_photons_by_k') or {}).items():
        package['photons_unfolded_k_%s' % k] = np.asarray(values, dtype=float)
    # the photon-yield closure of the unfolding alone at the result's photon iterations (full: an identity by
    # construction; the halves carry the information)
    unfolding = (result.get('unfolding_closure') or {}).get('results') or {}
    for mode in ('full', 'A_to_B', 'B_to_A'):
        view = unfolding.get('%s|1D|response_only|k%d' % (mode, int(record['photon_iterations'])))
        if view:
            package['photon_closure_%s_truth' % mode] = np.asarray(view['test_truth'], dtype=float)
            package['photon_closure_%s_unfolded' % mode] = np.asarray(view['unfolded'], dtype=float)
            package['photon_closure_%s_sigma' % mode] = np.asarray(view['bootstrap_sigma_of_difference'], dtype=float)
    # the diagnostic iteration scan of the result's own prior form
    for k, values in entry['per_photon_by_k'][s.prior_form].items():
        package['per_photon_k_%s' % k] = np.asarray(values, dtype=float)
    combinatoric = result.get('combinatoric') or {}
    if combinatoric.get('method') == 'data_dphi':
        # the data-driven estimate: what it removed from the purity-corrected spectrum
        package['combinatoric_purity_corrected_background'] = np.asarray(combinatoric['purity_corrected_background'], dtype=float)
    closure = result.get('closure') or {}
    for mode, view in (closure.get('modes') or {}).items():
        for key in ('truth', 'unfolded', 'sigma_of_difference'):   # photon-pT bins, then their integral
            package['closure_%s_%s' % (mode, key)] = np.asarray(view[key], dtype=float)
    return package, checks


def write(path, package, checks, meta):
    np.savez_compressed(path, **package)
    path.with_suffix('.json').write_text(json.dumps(dict(meta, checks=checks, keys=sorted(package)), indent=1, default=str) + '\n')
