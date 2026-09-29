"""The simulation response on the finest grid its source carries, for the IAN response figures: every response is drawn
on its native grid with the analysis edges overlaid, never rebinned onto the unequal analysis bins, and each truth row is
divided by its row sum plus misses (the unfolder's denominator) so that it sums to that bin's efficiency.

    from_files(M, files, config, key)     the FillPhotonJetResponse outputs of the entry's samples, streamed on the
                                          production's own axes (each file's `grids`), summed per sample and weighted as
                                          the measurement weights them (package.sample_weight, closure.weighted_package),
                                          with the Kish effective entries of the photon and pair responses from the
                                          files' own sumw2
    from_samples(M, samples, settings)    the measurement's samples on their response grids, for sealed inputs that carry
                                          no finer grid (no sumw2, so no effective entries)

Each returns {package, source, photon_neff, pair_neff}; figure_package.build(native=...) checks that restricting the
package reproduces the measurement window before it stores anything.
"""
import json

import numpy as np
import uproot


def _grids(M, settings, fine):
    return M._contract.Grids(pt_edges=tuple(settings.pt_edges), xj_edges=tuple(float(x) for x in fine['xj']),
                             reco_pt_edges=tuple(float(x) for x in fine['reco_pt']),
                             truth_pt_edges=tuple(float(x) for x in fine['truth_pt']))


def _neff(total, sumw2):
    return np.divide(total ** 2, sumw2, out=np.zeros_like(total, dtype=float), where=sumw2 > 0)


def from_files(M, files, config, key):
    """``files``: {sample: [FillPhotonJetResponse outputs]} of one systems entry (the files compute() measured)."""
    block = config['systems'][key]
    settings = M.Settings.from_config(config, key)
    suffix = M._package.name_suffix(block['class'], settings.jet_pt_min_gev)
    fine = None
    sums, weights = {}, {}
    photon_w, photon_w2, pair_w, pair_w2 = 0.0, 0.0, 0.0, 0.0
    for name in block['samples']:
        weight = M._package.sample_weight(config['samples'][name], config.get('reference_cross_section_pb') or {})
        weights[name] = weight
        total = None
        for path in files[name]:
            with uproot.open(path) as f:
                grids_here = json.loads(str(f['grids']))
                prefix = 'response/%s/jet%d/' % (block['class'], settings.jet_pt_min_gev) if 'meta/package' in f else ''
                tail = '' if prefix else suffix
                photon = f[prefix + 'photon_response_truth_x_reco' + tail]
                pairs = f[prefix + 'xj_response_truth_global_x_reco_global' + tail]
                photon_values, photon_var = photon.values(), photon.variances()
                pair_values, pair_var = pairs.values(), pairs.variances()
            if fine is None:
                fine = grids_here
            elif grids_here != fine:
                raise ValueError('%s: native grids differ from the first file of the entry' % path)
            row = M._package.simulation_sums(path, block['class'], settings.jet_pt_min_gev, fine['reco_pt'], fine['truth_pt'], fine['xj'])
            total = row if total is None else {k: total[k] + row[k] for k in total}
            photon_w, photon_w2 = photon_w + weight * photon_values, photon_w2 + weight * weight * photon_var
            pair_w, pair_w2 = pair_w + weight * pair_values, pair_w2 + weight * weight * pair_var
        if total is None:
            raise ValueError('no simulation files for %s' % name)
        sums[name] = total
    package = M._closure.weighted_package(_grids(M, settings, fine), weights, sums)
    return dict(package=package, source='production grid (FillPhotonJetResponse axes)',
                photon_neff=_neff(photon_w, photon_w2), pair_neff=_neff(pair_w, pair_w2))


def from_samples(M, samples, settings):
    """The measurement's own samples on its response grids (the finest grid of sealed reference inputs)."""
    weights = {name: sample.weight for name, sample in samples.items()}
    sums = M._resampling.occurrence_sums(M._closure.sample_arrays(samples))
    return dict(package=M._closure.weighted_package(settings.grids(), weights, sums),
                source='measurement response grid (the inputs carry no finer grid)', photon_neff=None, pair_neff=None)
