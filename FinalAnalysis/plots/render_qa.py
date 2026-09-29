#!/usr/bin/env python3
"""Prepare and draw photon-selection and event QA from saved sample packages.

Preparation routes histograms through meta/inventory, applies source-recorded
sample normalization, and saves ordinary histogram and numerical figure inputs.
Drawing reads those saved inputs; the ABCD solutions, quantiles and descriptive
fits are computed only during preparation. No result changes the measurement.

    python -m FinalAnalysis.plots.render_qa --inputs RESULTS/packages/qa_inputs.npz --output DIR
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import uproot
import yaml

from FinalAnalysis.src.internal import qa as QA
from FinalAnalysis.src.internal import centrality as CENTRALITY        # the data target of the centrality weights
from FinalAnalysis.src.internal import inputs as PACKAGE              # the samples' normalisation

from . import render as R  # house style: brand, colours, palette, save

plt = R.plt
CLASSES = {'pp': [('all', 'pp')], 'auau': [('cent0_20', 'Au+Au 0–20%'), ('cent20_50', 'Au+Au 20–50%'), ('cent50_80', 'Au+Au 50–80%')]}
FILE_TAG = {'cent0_20': '0to20', 'cent20_50': '20to50', 'cent50_80': '50to80'}
LABELS = {'data': ('data',), 'phosim': ('sig', 'phoother'), 'jetsim': ('bkg', 'jetprompt')}
PHOTON_SAMPLE_PREFIXES = ('pp_photon', 'auau_photon')
NOTES = {}
NUMERICS = {}
PREPARING = False


def numerical(kind, function, *args):
    """Preparation computes once; rendering requires the exact saved numerical input."""
    key = QA.identity(kind, *args)
    if key not in NUMERICS:
        if not PREPARING: raise ValueError('missing prepared QA numerical input: ' + kind)
        NUMERICS[key] = function(*args)
    return NUMERICS[key]



class Hist:
    """Summed histogram: values and variances without flow, and the axis edges."""

    def __init__(self, values, variances, edges):
        self.values, self.variances, self.edges = values, variances, edges

    def __add__(self, other):
        if other is None:
            return self
        if len(self.edges) != len(other.edges) or any(not np.array_equal(a, b) for a, b in zip(self.edges, other.edges)):
            raise ValueError('cannot sum QA histograms with different axes')
        return Hist(self.values + other.values, self.variances + other.variances, self.edges)

    __radd__ = __add__

    def scaled(self, w):
        return Hist(self.values * w, self.variances * w * w, self.edges)


def read(f, key):
    h = f[key]
    values = np.asarray(h.values(flow=False), dtype=float)
    variances = h.variances(flow=False)
    variances = values.copy() if variances is None else np.asarray(variances, dtype=float)
    return Hist(values, variances, [np.asarray(ax.edges(), dtype=float) for ax in h.axes])


class Store:
    """One verified TreeToHists collection, routed by sample metadata and inventory."""

    def __init__(self, root, config, sample_weights, unit_weights):
        self.root, self.config, self.unit = Path(root), config, unit_weights
        self.override = sample_weights or {}
        self.cache = {}
        self.saved = False
        manifest = self.root if self.root.is_file() else self.root / 'manifest.json'
        self.collection = PACKAGE.read_collection(manifest)
        self.trigger_records = {}
        self.centrality_target = None

    def files(self, kind, system=None, product=None, sample=None):
        if kind not in ('qa_data', 'data', 'qa_sim', 'response'): raise ValueError(kind)
        out = []
        for item in self.collection.values():
            meta = item['metadata']; physical = {'PP': 'pp', 'AuAu': 'auau'}[meta['system']]
            name = meta['name'] if meta['name'].startswith(physical + '_') else physical + '_' + meta['name']
            if kind in ('qa_data', 'data'):
                if physical == system and meta['family'] == 'data': out.append(item['path'])
            elif name == sample and product == physical + '_' + meta['family']:
                out.append(item['path'])
        return sorted(out)

    def trigger_rows(self, system):
        if system in self.trigger_records: return self.trigger_records[system]
        rows = []
        for path in self.files('qa_data', system=system):
            with uproot.open(path) as f:
                key = 'qa/events/trigger_bits_raw'
                if key not in f: continue
                edges = f[key].axes[0].edges()
                for index, raw in enumerate(f['records/sources']['json'].array(library='np')):
                    source = json.loads(raw); meta = source['metadata']
                    info = meta.get('trigger_selection') or meta.get('trigger_run_info')
                    run = meta['source'].get('run')
                    if info is None or run is None: raise ValueError('trigger QA lacks source-bound run/table')
                    values, variance = PACKAGE.source_histograms(f, [key], index)[key]
                    rows.append((run, info, Hist(values, variance, [edges])))
        self.trigger_records[system] = rows
        return rows

    def target(self):
        if self.centrality_target is None:
            values, edges, records = CENTRALITY.trigger_histogram(self.files('qa_data', system='auau'), CENTRALITY.DATA_TRIGGER)
            self.centrality_target = (values, edges, records)
        return self.centrality_target

    def samples(self, system):
        """{(product, sample): family} of the system's simulation."""
        out = {}
        for item in self.collection.values():
            meta = item['metadata']; physical = {'PP': 'pp', 'AuAu': 'auau'}[meta['system']]
            if physical != system or meta['family'] == 'data': continue
            name = meta['name'] if meta['name'].startswith(system + '_') else system + '_' + meta['name']
            out[(system + '_' + meta['family'], name)] = 'phosim' if meta['family'] == 'photonjet' else 'jetsim'
        return out

    def weight(self, sample):
        if self.unit:
            return 1.0
        meta = next(v['metadata'] for v in self.collection.values()
                    if v['metadata']['name'] == sample or
                    {'PP': 'pp', 'AuAu': 'auau'}[v['metadata']['system']] + '_' + v['metadata']['name'] == sample)
        recorded = meta['normalization']
        if recorded.get('applied_in_histograms') is not False:
            raise ValueError('QA sample normalization must be recorded as unapplied: ' + sample)
        if sample in self.override:
            weight = float(self.override[sample])
            if not np.isfinite(weight) or weight <= 0: raise ValueError('invalid sample weight')
            if recorded.get('weight_pb_per_event') not in (None, weight):
                raise ValueError('QA override differs from recorded sample weight: ' + sample)
            return weight
        spec = (self.config.get('samples') or {}).get(sample)
        if spec is None:
            raise KeyError('no normalisation for sample %s (configuration or --sample-weights)' % sample)
        for field in ('source_factor', 'generated_events', 'cross_section_pb', 'reference'):
            if recorded.get(field) is not None and recorded[field] != spec.get(field):
                raise ValueError('QA normalization ingredient differs: ' + sample + '.' + field)
        weight = PACKAGE.sample_weight(spec, self.config.get('reference_cross_section_pb') or {})
        if recorded.get('weight_pb_per_event') not in (None, weight):
            raise ValueError('QA normalization differs from recorded sample weight: ' + sample)
        if recorded.get('numerator_pb') is not None and recorded['numerator_pb'] / recorded['generated_events'] != weight:
            raise ValueError('QA normalization ingredients differ: ' + sample)
        return weight

    def summed(self, kind, key, names=None, **where):
        """Every histogram of a file group (or only `names`), summed (cached per group)."""
        key = key + (tuple(sorted(names)),) if names else key
        key = json.dumps(key)
        if key not in self.cache:
            if self.saved: raise ValueError('missing prepared QA histogram group: ' + key)
            sums = {}
            for path in self.files(kind, **where):
                with uproot.open(path) as f:
                    rows = json.loads(str(f['meta/inventory']))
                    wanted = {'qa_data': 'qa', 'qa_sim': 'qa', 'data': 'spectra', 'response': 'response'}[kind]
                    mapped = {row['key']: row['path'] for row in rows
                              if row['kind'] == wanted and row['path'] in f}
                    for name in (set(mapped) & set(names)) if names else mapped:
                        obj = f[mapped[name]]
                        if not hasattr(obj, 'values') or not hasattr(obj, 'axes'):
                            continue
                        h = read(f, mapped[name])
                        if np.any(h.values) or np.any(h.variances): sums[name] = h + sums.get(name)
            self.cache[key] = sums
        return self.cache[key]

    def data(self, system, name, kind='qa_data'):
        return self.summed(kind, (kind, system), system=system).get(name)

    def sim(self, system, family, name, kind='qa_sim', per_sample=False, products=None):
        """The weighted sum over the system's samples of one family (or {sample: hist} with per_sample)."""
        out = {}
        for (product, sample), fam in self.samples(system).items():
            if fam != family or (products and product not in products):
                continue
            h = self.summed(kind, (kind, product, sample), product=product, sample=sample).get(name)
            if h is not None:
                out[sample] = h.scaled(self.weight(sample))
        if per_sample:
            return out
        total = None
        for h in out.values():
            total = h + total
        return total

    def labels(self, system, family, stem, c, **kw):
        """Sum over the candidate labels of a family: <stem>_<label>_<c>."""
        total = None
        for label in LABELS[family]:
            name = '%s_%s_%s' % (stem, label, c) if c else '%s_%s' % (stem, label)
            h = self.data(system, name) if family == 'data' else self.sim(system, family, name, **kw)
            total = h + total if h is not None else total
        return total


# ----------------------------------------------------------------------------------------------------------------
# helpers

def window(edges, lo, hi):
    return (edges[:-1] >= lo - 1e-9) & (edges[1:] <= hi + 1e-9)


def project(h, axis, ranges=None):
    """Sum a histogram over every axis but `axis`, restricting the axes named in ranges {index: (lo, hi)}."""
    ranges = ranges or {}
    v, var = h.values, h.variances
    for ax in sorted(set(range(v.ndim)) - {axis}, reverse=True):
        sel = np.ones(v.shape[ax], dtype=bool)
        if ax in ranges:
            sel &= window(h.edges[ax], *ranges[ax])
        v, var = np.compress(sel, v, axis=ax).sum(axis=ax), np.compress(sel, var, axis=ax).sum(axis=ax)
    return v, var, h.edges[axis]


def unit(v, e=None):
    s = v.sum()
    if s <= 0:
        return v, e
    return v / s, (None if e is None else e / s)


def watermark(fig, store):
    if store.unit:
        fig.text(.5, .5, 'UNIT SAMPLE WEIGHTS (mechanical test)', ha='center', va='center', rotation=30, fontsize=26,
                 color='#C62828', alpha=.25)


def save(fig, store, name):
    if PREPARING:
        plt.close(fig)
        return
    watermark(fig, store)
    if name.endswith('.pdf'):
        p = R.OUT / 'fig' / name
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p)
        R.FILES.append({'path': str(p.relative_to(R.OUT)), 'sha256': R.sha(p)})
        plt.close(fig)
        print('PDF', name, flush=True)
    else:
        R.save(fig, name)


def missing(name, why):
    NOTES[name] = 'not drawn: ' + why
    print('SKIP', name, why, flush=True)


def text_block(a, lines, x=.04, y=.95, size=10.5):
    a.text(x, y, '\n'.join(lines), transform=a.transAxes, va='top', fontsize=size, linespacing=1.35)


def energy(system):
    return r'$\sqrt{s}=200$ GeV' if system == 'pp' else r'$\sqrt{s_{NN}}=200$ GeV'


def classes_of(system):
    return CLASSES[system]


def summed_classes(store, system, fetch):
    total = None
    for c, _ in classes_of(system):
        h = fetch(c)
        total = h + total if h is not None else total
    return total


# ----------------------------------------------------------------------------------------------------------------
# photon selection

def ncb_score(store, system):
    name = 'photonOnly/%s_ncb_score_data_mc.png' % ('pp' if system == 'pp' else 'auau_0to80')
    got = {f: summed_classes(store, system, lambda c: store.labels(system, f, 'ncb', c)) for f in ('data', 'phosim', 'jetsim')}
    if got['data'] is None:
        return missing(name, 'no DATA QA')
    f, a = R.axes()
    for fam, color, text in (('phosim', R.RED, 'Photon+jet simulation'), ('jetsim', R.BLUE, 'Inclusive-jet simulation')):
        if got[fam] is None:
            continue
        v, var, edges = project(got[fam], 0, {1: (16, 35)})
        keep = edges[:-1] >= 0
        v, e = unit(v[keep], np.sqrt(var[keep]))
        a.stairs(v, edges[np.append(keep, True)], color=color, lw=1.5, label=text)
    v, var, edges = project(got['data'], 0, {1: (16, 35)})
    keep = edges[:-1] >= 0
    v, e = unit(v[keep], np.sqrt(var[keep]))
    x = (edges[:-1] + edges[1:])[keep] / 2
    R.points(a, x, v, e, R.BLACK, 'o', 'DATA')
    a.axvline(.5, ls='--', color=R.BLACK, lw=.9)
    a.set_yscale('log'); a.set_xlim(0, 1.0); a.set_xlabel('NCB classifier score'); a.set_ylabel('Normalized entries')
    R.brand(a); text_block(a, [('pp, ' if system == 'pp' else 'Au+Au 0–80%, ') + energy(system), r'$16<p_T^\gamma<35$ GeV'], y=.88)
    a.legend(loc='upper right', fontsize=9.5)
    save(f, store, name)


def ncb_eta_time(store):
    name = 'photonOnly/pp_ncb_eta_time_before.png'
    h = store.data('pp', 'etatime_data_all')
    if h is None:
        return missing(name, 'no pp DATA QA')
    from matplotlib.colors import LogNorm
    f, a = R.axes(matrix=True)
    v = h.values
    image = a.pcolormesh(h.edges[0], h.edges[1], np.ma.masked_less_equal(v.T, 0), norm=LogNorm(max(v[v > 0].min(), 1e-12), v.max()),
                         cmap='viridis', shading='flat', rasterized=True)
    a.set_xlabel(r'Cluster $\eta$'); a.set_ylabel('Cluster time − MBD time [ns]')
    R.brand(a, box=True); text_block(a, ['pp, ' + energy('pp'), r'$E_T^{\mathrm{cluster}}>10$ GeV', 'before the NCB requirement'], y=.88)
    f.colorbar(image, ax=a, pad=.025, fraction=.045).set_label('Clusters')
    save(f, store, name)


SHAPES = [('e11_to_e33', 'e11_over_e33', r'$E_{11}/E_{33}$'), ('e32_to_e35', 'e32_over_e35', r'$E_{32}/E_{35}$'),
          ('bdt', 'bdt', r'photon-ID score $s_{\mathrm{BDT}}$'), ('weta_cogx', 'weta_cogx', r'$w_\eta^{\mathrm{COG}}$'),
          ('wphi_cogx', 'wphi_cogx', r'$w_\phi^{\mathrm{COG}}$'), ('et1', 'et1', r'$E_{T,1}$ fraction'), ('et2', 'et2', r'$E_{T,2}$ fraction'),
          ('et3', 'et3', r'$E_{T,3}$ fraction'), ('et4', 'et4', r'$E_{T,4}$ fraction')]


def shower_shapes(store):
    """dis_<var>_cent<k>_pt0_cut1.pdf: after preselection (stage 1), 15-20 GeV, the three Au+Au classes (cent0..2)."""
    for k, (c, text) in enumerate(classes_of('auau')):
        for slot, var, xlabel in SHAPES:
            name = 'shower_shapes/dis_%s_cent%d_pt0_cut1.pdf' % (slot, k)
            data = store.labels('auau', 'data', 'ss_' + var, c)
            if data is None:
                missing(name, 'no Au+Au DATA QA')
                continue
            f, a = R.axes()
            for label, fam, color, legend in (('sig', 'phosim', R.RED, 'Photon+jet simulation (prompt)'),
                                              ('bkg', 'jetsim', R.BLUE, 'Inclusive-jet simulation')):
                h = store.sim('auau', fam, 'ss_%s_%s_%s' % (var, label, c))
                if h is None:
                    continue
                v, var_, edges = project(h, 0, {1: (15, 20), 2: (1, 2)})
                a.stairs(unit(v)[0], edges, color=color, lw=1.5, label=legend)
            v, var_, edges = project(data, 0, {1: (15, 20), 2: (1, 2)})
            v, e = unit(v, np.sqrt(var_))
            R.points(a, (edges[:-1] + edges[1:]) / 2, v, e, R.BLACK, 'o', 'DATA')
            a.set_xlabel(xlabel); a.set_ylabel('Normalized counts')
            R.brand(a); text_block(a, [text + ', ' + energy('auau'), r'$15<p_T^\gamma<20$ GeV, $|\eta^\gamma|<0.7$', 'after preselection'], y=.88)
            a.legend(loc='upper right', fontsize=9.5)
            save(f, store, name)


def widths(store):
    for label, fam, who in (('sig', 'phosim', 'signal'), ('bkg', 'jetsim', 'background')):
        for var, slot, xlabel in (('weta_cogx', 'full_cluster', r'$w_\eta^{\mathrm{COG}}$'), ('weta33_cogx', '3x3', r'$w_\eta^{3\times3,\mathrm{COG}}$')):
            name = 'photonOnly/auau_%s_%s_width.png' % (who, slot)
            h = store.sim('auau', fam, 'width_%s_%s' % (var, label))
            if h is None:
                missing(name, 'no Au+Au %s simulation QA' % fam)
                continue
            f, a = R.axes()
            for (lo, hi), color, text in (((0, 4), R.BLUE, '0–20%'), ((4, 10), R.ORANGE, '20–50%'), ((10, 16), R.GREEN, '50–80%')):
                v, _, edges = project(h, 0, {1: (15, 35), 2: (lo, hi)})
                a.stairs(unit(v)[0], edges, color=color, lw=1.6, label=text)
            a.set_xlim(0, .7); a.set_xlabel(xlabel); a.set_ylabel('Unit-normalized candidates / 0.01')
            R.brand(a); text_block(a, ['Au+Au, ' + energy('auau'), ('Truth-matched prompt photons' if label == 'sig' else 'Inclusive-jet candidates')
                                        + ', preselected', r'$15<p_T^\gamma<35$ GeV'], y=.88)
            a.legend(loc='upper right', title='Centrality', fontsize=9.5)
            save(f, store, name)


def bdt_per_cent(store):
    name = 'photonOnly/auau_bdt_score_perCent.png'
    sig, bkg = store.sim('auau', 'phosim', 'bdt_sig'), store.sim('auau', 'jetsim', 'bdt_bkg')
    if sig is None or bkg is None:
        return missing(name, 'needs Au+Au photon and inclusive simulation QA')
    f, axs = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
    f.subplots_adjust(left=.06, right=.99, bottom=.14, top=.93, wspace=.04)
    for a, ((lo, hi), text) in zip(axs, (((0, 4), '0–20%'), ((4, 10), '20–50%'), ((10, 16), '50–80%'))):
        R.style(a)
        for h, color, legend in ((sig, R.BLUE, 'Prompt photons (photon+jet)'), (bkg, R.ORANGE, 'Inclusive-jet candidates')):
            v, _, edges = project(h, 0, {1: (15, 35), 2: (lo, hi)})
            a.stairs(unit(v)[0], edges, color=color, lw=1.4, label=legend)
        a.set_xlabel(r'BDT score, $s_{\mathrm{BDT}}$'); a.set_xlim(0, 1)
        a.text(.04, .95, text, transform=a.transAxes, va='top', fontsize=13, weight='bold')
    axs[0].set_ylabel('Fraction of candidates'); axs[0].set_yscale('log')
    R.brand(axs[1]); axs[1].legend(loc='upper left', bbox_to_anchor=(0, .86), fontsize=9.5)
    axs[2].text(.04, .86, 'Embedded PYTHIA 8\n' + r'$15\leq p_T^\gamma<35$ GeV, no isolation', transform=axs[2].transAxes, va='top', fontsize=9.5)
    save(f, store, name)


def e11_flow(store, timing_window=(-40.0, -5.0)):
    """Rows 0-20, 20-50, 50-80 Au+Au and pp; columns before preselection, after preselection, after tight ID (stages
    0, 1, 2); 15-35 GeV. Green: DATA with cluster time - MBD time in the timing window (NCB-enriched), before
    preselection, scaled to the DATA of the panel."""
    name = 'shower_shapes/e11_e33_selection_flow.png'
    rows = [('auau', c, t) for c, t in classes_of('auau')] + [('pp', 'all', 'pp reference')]
    f, axs = plt.subplots(len(rows), 3, figsize=(12, 3.0 * len(rows)), sharex=True)
    f.subplots_adjust(left=.07, right=.99, bottom=.06, top=.96, hspace=.08, wspace=.12)
    drawn = False
    for r, (system, c, text) in enumerate(rows):
        data = store.labels(system, 'data', 'ss_e11_over_e33', c)
        for s, title in enumerate(('Before preselection', 'After preselection', 'After tight ID')):
            a = axs[r][s]
            R.style(a)
            if r == 0:
                a.set_title(title, fontsize=11)
            if data is None:
                continue
            drawn = True
            for label, fam, color in (('sig', 'phosim', R.RED), ('bkg', 'jetsim', R.BLUE)):
                h = store.sim(system, fam, 'ss_e11_over_e33_%s_%s' % (label, c))
                if h is not None:
                    v, _, edges = project(h, 0, {1: (15, 35), 2: (s, s + 1)})
                    a.stairs(unit(v)[0], edges, color=color, lw=1.2)
            v, var, edges = project(data, 0, {1: (15, 35), 2: (s, s + 1)})
            total = v.sum()
            v, e = unit(v, np.sqrt(var))
            R.points(a, (edges[:-1] + edges[1:]) / 2, v, e, R.BLACK, 'o')
            if s == 0:
                a.text(.03, .93, text, transform=a.transAxes, va='top', fontsize=10, weight='bold')
                tagged = store.data(system, 'e11e33_dt_data_%s' % c)
                if tagged is not None and total > 0:
                    t, _, t_edges = project(tagged, 0, {1: (15, 35), 2: timing_window})
                    a.stairs(t / total, t_edges, color=R.GREEN, lw=1.2)
                    a.text(.03, .80, 'NCB timing tag: %d' % round(t.sum()), transform=a.transAxes, va='top', fontsize=8.5, color=R.GREEN)
    if not drawn:
        plt.close(f)
        return missing(name, 'no DATA QA')
    for a in axs[-1]:
        a.set_xlabel(r'$E_{11}/E_{33}$')
    axs[len(rows) // 2][0].set_ylabel('Fraction of candidates')
    R.brand(axs[0][2], size=12)
    handles = [plt.Line2D([], [], color=R.BLACK, marker='o', ls='none', label='DATA'),
               plt.Line2D([], [], color=R.RED, label='Photon+jet simulation'), plt.Line2D([], [], color=R.BLUE, label='Inclusive-jet simulation'),
               plt.Line2D([], [], color=R.GREEN, label='NCB timing-tagged DATA (scaled)')]
    axs[0][2].legend(handles=handles, loc='upper left', bbox_to_anchor=(0, .82), fontsize=8)
    NOTES[name] = 'NCB timing tag: cluster - MBD time in [%g, %g) ns' % timing_window
    save(f, store, name)


def truth_isolation(store):
    name = 'photonOnly/pp_auau_truth_isolation_direct_fragmentation_2x2.png'
    panels = [('pp', 'all', 'pp, ' + energy('pp'))] + [('auau', c, t + ', ' + energy('auau')) for c, t in reversed(classes_of('auau'))]
    got = [(s, c, t, store.sim(s, 'phosim', 'truth_iso_%s' % c)) for s, c, t in panels]
    if all(h is None for *_, h in got):
        return missing(name, 'no photon-family simulation QA')
    f, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    f.subplots_adjust(left=.09, right=.98, bottom=.09, top=.92, hspace=.18, wspace=.08)
    cuts = np.arange(1, 11)
    for a, (system, c, text, h) in zip(axs.ravel(), got):
        R.style(a)
        a.set_title(text, fontsize=11, loc='left')
        if h is None:
            continue
        for cls, color, who in ((1, R.BLUE, 'Direct'), (2, R.ORANGE, 'Fragmentation')):
            for (lo, hi), marker in (((15, 20), 'o'), ((25, 30), 's')):
                v, var, edges = project(h, 0, {1: (lo, hi), 2: (cls, cls + 1)})
                total = v.sum()
                if total <= 0:
                    continue
                frac = np.array([v[edges[1:] <= cut + 1e-9].sum() / total for cut in cuts])
                err = np.sqrt(frac * (1 - frac) / max(total, 1))
                R.points(a, cuts, frac, err, color, marker, '%s, %g–%g GeV' % (who, lo, hi))
        a.axvline(4, ls='--', color=R.BLACK, lw=.9)
    axs[0][0].legend(loc='lower right', fontsize=8.5)
    for a in axs[1]:
        a.set_xlabel(r'Truth isolation threshold $E_{T,\mathrm{cut}}^{\mathrm{iso}}$ [GeV]')
    for a in axs[:, 0]:
        a.set_ylabel(r'Fraction passing $E_T^{\mathrm{iso,truth}}<E_{T,\mathrm{cut}}^{\mathrm{iso}}$')
    f.text(.09, .955, r'$\mathbf{\mathit{sPHENIX}}$ Internal      Photon+jet simulation, truth $R=0.3$, $|\eta^\gamma|<0.7$', fontsize=12)
    save(f, store, name)


def truth_pt_after_isolation(store):
    for system, c, text in [('pp', 'all', 'pp')] + [('auau', c, t) for c, t in classes_of('auau')]:
        name = 'photonOnly/%s_truth_photon_pt_after_isolation.png' % ('pp' if system == 'pp' else 'auau_' + FILE_TAG[c])
        h = store.sim(system, 'phosim', 'truth_iso_%s' % c)
        if h is None:
            missing(name, 'no photon-family simulation QA')
            continue
        f, a, b = R.axes(ratio=True)
        parts = {}
        for cls in (1, 2):
            v, var, edges = project(h, 1, {0: (0, 4), 2: (cls, cls + 1)})
            parts[cls] = (v, var)
        edges = h.edges[1]
        width = np.diff(edges)
        x = (edges[:-1] + edges[1:]) / 2
        allv, allvar = parts[1][0] + parts[2][0], parts[1][1] + parts[2][1]
        m = window(edges, 15, 35)
        R.points(a, x[m], (allv / width)[m], (np.sqrt(allvar) / width)[m], R.BLACK, 'o', 'All prompt photons')
        R.points(a, x[m], (parts[1][0] / width)[m], (np.sqrt(parts[1][1]) / width)[m], R.RED, 's', 'Direct')
        R.points(a, x[m], (parts[2][0] / width)[m], (np.sqrt(parts[2][1]) / width)[m], R.BLUE, '^', 'Fragmentation')
        a.set_yscale('log'); a.set_ylabel('Weighted photons / GeV'); a.legend(loc='upper right')
        R.brand(a); text_block(a, [text + ', ' + energy(system), r'$E_{T,\mathrm{truth}}^{\mathrm{iso}}<4$ GeV, $R=0.3$, $|\eta^\gamma_{\mathrm{truth}}|<0.7$'], y=.88)
        frac = np.divide(parts[1][0], allv, out=np.zeros_like(allv), where=allv > 0)
        R.points(b, x[m], frac[m], None, R.RED, 's')
        b.set_ylabel('Direct / total', fontsize=11); b.set_xlabel(r'Truth photon $p_T$ [GeV]'); b.set_xlim(15, 35)
        save(f, store, name)


def isolation_overlays(store):
    for kind, stem in (('noUE', 'isonoueslice_data'), ('withUE', 'isoslice_data')):
        name = 'isolation/isoCentOverlay_%s_pT_14_16.png' % kind
        h = store.data('auau', stem)
        if h is None:
            missing(name, 'no Au+Au DATA QA')
            continue
        f, a = R.axes()
        for (lo, hi), color, text in (((0, 2), R.RED, '0–10%'), ((2, 4), R.BLUE, '10–20%'), ((4, 8), R.GREEN, '20–40%'),
                                      ((8, 12), R.PURPLE, '40–60%'), ((12, 16), R.ORANGE, '60–80%')):
            v, var, edges = project(h, 0, {1: (14, 16), 2: (lo, hi)})
            v, var = v.reshape(-1, 5).sum(axis=1), var.reshape(-1, 5).sum(axis=1)   # 0.5 GeV bins
            e5 = edges[::5]
            v, err = unit(v, np.sqrt(var))
            shown = v > 0
            R.points(a, ((e5[:-1] + e5[1:]) / 2)[shown], v[shown], err[shown], color, 'o', text)
        a.set_xlim(-10, 40); a.set_xlabel(r'$E_T^{\mathrm{iso}}$ [GeV]'); a.set_ylabel('Normalized to unit area')
        R.brand(a); text_block(a, ['Au+Au, ' + energy('auau'), r'$14<p_T^\gamma<16$ GeV, $\Delta R_{\mathrm{cone}}<0.3$',
                                   'with SUB1 UE subtraction' if kind == 'withUE' else 'without UE subtraction'], y=.88)
        a.legend(loc='upper right', ncol=1, fontsize=9.5)
        save(f, store, name)


def isolation_stacks(store, sideband=6.0):
    for system, c, text in [('pp', 'all', 'pp')] + [('auau', c, t) for c, t in classes_of('auau')]:
        name = 'isolation/%s_isoStack.png' % ('pp_nominal_topo' if system == 'pp' else 'auau_' + FILE_TAG[c])
        data = store.data(system, 'iso_data_%s' % c)
        sig = store.sim(system, 'phosim', 'iso_sig_%s' % c)
        if data is None or sig is None:
            missing(name, 'needs DATA and photon-family simulation QA')
            continue
        tight, tvar, edges = project(data, 0, {1: (15, 35), 2: (0, 1)})
        nontight, _, _ = project(data, 0, {1: (15, 35), 2: (1, 2)})
        signal, _, _ = project(sig, 0, {1: (15, 35), 2: (0, 1)})
        side = edges[:-1] >= sideband
        scale_b = tight[side].sum() / nontight[side].sum() if nontight[side].sum() > 0 else 0.0
        background = nontight * scale_b
        rest = tight.sum() - background.sum()
        signal = signal * (rest / signal.sum() if signal.sum() > 0 else 0.0)
        # 1 GeV bins
        k = int(round(1.0 / np.diff(edges)[0]))
        n = len(tight) // k * k
        re = edges[:n + 1:k]
        rebin = lambda y: y[:n].reshape(-1, k).sum(axis=1)
        f, a = R.axes()
        a.stairs(rebin(background), re, fill=True, color='#F4B6B6', label='Non-tight data control')
        a.stairs(rebin(background) + rebin(signal), re, baseline=rebin(background), fill=True, color='#B8C0EE', label='Signal MC')
        a.stairs(rebin(background) + rebin(signal), re, color=R.BLUE, lw=1.2)
        R.points(a, (re[:-1] + re[1:]) / 2, rebin(tight), np.sqrt(rebin(tvar)), R.BLACK, 'o', 'Tight-photon data')
        a.set_xlim(-20 if system == 'auau' else -5, 20); a.set_xlabel(r'$E_T^{\mathrm{iso,reco}}$ [GeV]'); a.set_ylabel('Counts / 1 GeV')
        R.brand(a); text_block(a, [text + ', ' + energy(system), r'$15<E_T^\gamma<35$ GeV, $|\eta^\gamma|<0.7$',
                                   'non-tight normalized at $E_T^{iso}\\geq%g$ GeV' % sideband], y=.88)
        a.legend(loc='upper right', fontsize=9.5)
        NOTES[name] = dict(nontight_scale=scale_b, signal_scale_total=float(rest))
        save(f, store, name)


ABCD_EDGES = np.array([15, 17, 19, 21, 23, 25, 35], dtype=float)
FINE_EDGES = np.array([15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35], dtype=float)   # the purity closure studies' 2 GeV grid


def regrouped(values, edges, target):
    """Sum native 1 GeV bins into target bins (target edges must be native edges)."""
    idx = np.searchsorted(target, edges[:-1] + 1e-9, side='right') - 1
    out = np.zeros(values.shape[:-1] + (len(target) - 1,))
    ok = (idx >= 0) & (idx < len(target) - 1)
    np.add.at(out, (..., idx[ok]), values[..., ok]) if values.ndim > 1 else np.add.at(out, idx[ok], values[ok])
    return out


def abcd_counts(store):
    name = 'photonOnly/pp_auau_abcd_region_counts_1x4.png'
    panels = [('pp', 'all', 'pp')] + [('auau', c, t) for c, t in reversed(classes_of('auau'))]
    f, axs = plt.subplots(1, 4, figsize=(16, 4.4), sharey=True)
    f.subplots_adjust(left=.05, right=.99, bottom=.16, top=.9, wspace=.05)
    drawn = False
    x = (ABCD_EDGES[:-1] + ABCD_EDGES[1:]) / 2
    for a, (system, c, text) in zip(axs, panels):
        R.style(a); a.set_title(text, fontsize=11, loc='left')
        h = store.data(system, 'abcd_%s' % c, kind='data')
        if h is None:
            continue
        drawn = True
        counts = regrouped(h.values, h.edges[1], ABCD_EDGES) / np.diff(ABCD_EDGES)[None, :]
        errs = np.sqrt(regrouped(h.variances, h.edges[1], ABCD_EDGES)) / np.diff(ABCD_EDGES)[None, :]
        for r, (color, label) in enumerate(((R.BLACK, 'A: tight, isolated'), (R.RED, 'B: tight, non-isolated'),
                                            (R.BLUE, 'C: non-tight, isolated'), (R.GREEN, 'D: non-tight, non-isolated'))):
            R.points(a, x, counts[r], errs[r], color, 'o', label, np.diff(ABCD_EDGES) / 2)
        a.set_yscale('log'); a.set_xlabel(r'$E_T^{\gamma,\mathrm{reco}}$ [GeV]'); a.set_xlim(15, 35)
    if not drawn:
        plt.close(f)
        return missing(name, 'no DATA analysis fills')
    axs[0].set_ylabel('DATA event-leading candidates / GeV'); axs[-1].legend(loc='upper right', fontsize=8)
    R.brand(axs[0], size=12)
    save(f, store, name)


def leakage_of(store, system, c):
    """Signal leakage per region (A..D) x native reco pT, as the measurement takes it (the responses)."""
    total = None
    key = 'signal_leakage_truth_matched_abcd' + PACKAGE.name_suffix(c, store.config['measurement']['jet_pt_min_gev'])
    for (product, sample), fam in store.samples(system).items():
        if fam != 'phosim':
            continue
        h = store.summed('response', ('response', product, sample), names=[key], product=product, sample=sample).get(key)
        if h is not None:
            total = h.scaled(store.weight(sample)) + total
    return total


def purity_rows(counts, leakage):
    return numerical('purity', QA.purity_rows, counts, leakage)


def purity(store):
    name = 'photonOnly/data_purity_and_prompt_leakage_2x4_style.png'
    panels = [('pp', 'all', 'pp')] + [('auau', c, t) for c, t in reversed(classes_of('auau'))]
    f, axs = plt.subplots(2, 4, figsize=(16, 7.5), sharex=True)
    f.subplots_adjust(left=.05, right=.99, bottom=.09, top=.94, wspace=.18, hspace=.06)
    x = (ABCD_EDGES[:-1] + ABCD_EDGES[1:]) / 2
    drawn = False
    for k, (system, c, text) in enumerate(panels):
        a, b = axs[0][k], axs[1][k]
        R.style(a); R.style(b); a.set_title(text, fontsize=11, loc='left')
        h, leak = store.data(system, 'abcd_%s' % c, kind='data'), leakage_of(store, system, c)
        if h is None or leak is None:
            continue
        drawn = True
        counts = regrouped(h.values, h.edges[1], ABCD_EDGES)
        leakage = regrouped(leak.values, leak.edges[1], ABCD_EDGES)
        raw, corrected = purity_rows(counts, leakage)
        R.points(a, x, raw, None, R.BLACK, 'D', 'Raw ABCD', np.diff(ABCD_EDGES) / 2)
        R.points(a, x, corrected, None, R.BLUE, 'D', 'Leakage-corrected ABCD', np.diff(ABCD_EDGES) / 2, open_=True)
        a.set_ylim(0, 1.2)
        for r, color, label in ((1, R.BLACK, r'$S_B/S_A$'), (2, R.RED, r'$S_C/S_A$'), (3, R.BLUE, r'$S_D/S_A$')):
            frac = np.divide(leakage[r], leakage[0], out=np.full(len(x), np.nan), where=leakage[0] > 0)
            R.points(b, x, frac, None, color, 'o', label, np.diff(ABCD_EDGES) / 2)
        b.set_xlabel(r'$E_T^{\gamma,\mathrm{reco}}$ [GeV]'); b.set_xlim(15, 35); b.set_ylim(0, .6)
    if not drawn:
        plt.close(f)
        return missing(name, 'needs DATA analysis fills and the simulation responses')
    axs[0][0].set_ylabel('DATA photon purity'); axs[1][0].set_ylabel('Prompt-signal leakage fraction')
    axs[0][0].legend(loc='lower right', fontsize=9); axs[1][0].legend(loc='upper left', fontsize=9)
    R.brand(axs[0][0], size=12)
    save(f, store, name)


def simulation_purity(store, system, c, edges):
    """The simulation closure of the purity on ``edges`` (native edges): pseudo-data = every weighted simulation
    candidate (photon and inclusive families), its raw and leakage-corrected ABCD purity, and the truth purity =
    prompt-labelled Region-A leaders / all Region-A leaders; with the Region-A leaders on ``edges`` and on the native
    bins. None without both families and the responses."""
    parts = {label: store.sim(system, fam, 'abcd_%s_%s' % (label, c)) for fam in ('phosim', 'jetsim') for label in LABELS[fam]}
    leak = leakage_of(store, system, c)
    if any(v is None for v in parts.values()) or leak is None:
        return None
    total = sum(regrouped(h.values, h.edges[1], edges) for h in parts.values())
    prompt = sum(regrouped(parts[l].values, parts[l].edges[1], edges) for l in ('sig', 'jetprompt'))
    raw, corrected = purity_rows(total, regrouped(leak.values, leak.edges[1], edges))
    truth = np.divide(prompt[0], total[0], out=np.full(len(edges) - 1, np.nan), where=total[0] > 0)
    return dict(raw=raw, corrected=corrected, truth=truth, leaders=total[0], native_leaders=sum(h.values[0] for h in parts.values()),
                native_edges=next(iter(parts.values())).edges[1])


def purity_closure(store):
    """Simulation closure of the purity (simulation_purity) per class, on the ABCD bins."""
    name = 'photonOnly/purity_qa_pp_auau_simCLOSURE.png'
    panels = [('pp', 'all', 'pp')] + [('auau', c, t) for c, t in reversed(classes_of('auau'))]
    f, axs = plt.subplots(2, 4, figsize=(16, 7.5), sharex=True)
    f.subplots_adjust(left=.05, right=.99, bottom=.09, top=.94, wspace=.18, hspace=.06)
    x = (ABCD_EDGES[:-1] + ABCD_EDGES[1:]) / 2
    drawn = False
    for k, (system, c, text) in enumerate(panels):
        a, b = axs[0][k], axs[1][k]
        R.style(a); R.style(b); a.set_title(text, fontsize=11, loc='left')
        closure = simulation_purity(store, system, c, ABCD_EDGES)
        if closure is None:
            continue
        drawn = True
        raw, corrected, truth = closure['raw'], closure['corrected'], closure['truth']
        R.points(a, x, raw, None, R.BLACK, 'D', 'Raw ABCD', np.diff(ABCD_EDGES) / 2)
        R.points(a, x, corrected, None, R.BLUE, 'D', 'Leakage-corrected ABCD', np.diff(ABCD_EDGES) / 2, open_=True)
        R.points(a, x, truth, None, R.RED, 'o', 'Truth purity', np.diff(ABCD_EDGES) / 2)
        R.points(b, x, corrected / truth, None, R.BLUE, 'D', xerr=np.diff(ABCD_EDGES) / 2)
        b.axhline(1, ls='--', color=R.BLACK, lw=.8); b.set_ylim(.6, 1.6)
        a.set_ylim(0, 1.3); b.set_xlabel(r'$E_T^{\gamma,\mathrm{reco}}$ [GeV]'); b.set_xlim(15, 35)
    if not drawn:
        plt.close(f)
        return missing(name, 'needs both simulation families and the responses')
    axs[0][0].set_ylabel('Photon purity (simulation)'); axs[1][0].set_ylabel('Corrected ABCD / truth')
    axs[0][0].legend(loc='lower right', fontsize=9); R.brand(axs[0][0], size=12)
    save(f, store, name)


def efficiency(store):
    name = 'photonOnly/pp_si_di_auau_photon_efficiency_1x4.png'
    panels = [('pp', 'all', 'pp')] + [('auau', c, t) for c, t in reversed(classes_of('auau'))]
    f, axs = plt.subplots(1, 4, figsize=(16, 4.4), sharey=True)
    f.subplots_adjust(left=.05, right=.99, bottom=.16, top=.9, wspace=.05)
    drawn = False
    for a, (system, c, text) in zip(axs, panels):
        R.style(a); a.set_title(text, fontsize=11, loc='left')
        h = store.sim(system, 'phosim', 'eff_%s' % c)
        if h is None:
            continue
        drawn = True
        s = [project(h, 0, {1: (k, k + 1)})[0] for k in range(6)]
        x = (h.edges[0][:-1] + h.edges[0][1:]) / 2
        m = window(h.edges[0], 15, 35)
        ratio = lambda num, den: np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        eps_reco = ratio(s[1], s[0])
        eps_id_iso = ratio(s[5], s[3])
        eps_sel = ratio(s[5], s[0])
        R.points(a, x[m], eps_reco[m], None, R.BLACK, 's', r'$\varepsilon_{\mathrm{reco}}$')
        R.points(a, x[m], (eps_reco * eps_id_iso)[m], None, '#C2185B', 'o', r'$\varepsilon_{\mathrm{reco}}\,\varepsilon_{\mathrm{ID|iso}}$')
        R.points(a, x[m], eps_sel[m], None, R.GREEN, 'P', r'$\varepsilon_{\mathrm{sel}}$')
        a.set_xlabel(r'Truth photon $p_T$ [GeV]'); a.set_xlim(15, 35); a.set_ylim(0, 1.1)
    if not drawn:
        plt.close(f)
        return missing(name, 'no photon-family simulation QA')
    axs[0].set_ylabel('Photon efficiency'); axs[0].legend(loc='lower left', fontsize=9); R.brand(axs[0], size=12)
    save(f, store, name)


def minbias_pretraining(store):
    name = 'photonOnly/auau_minbias_pretraining_before_after.png'
    h = store.sim('auau', 'phosim', 'evt_calo_phosim')
    if h is None:
        return missing(name, 'no Au+Au photon-family simulation QA')
    from matplotlib.colors import LogNorm
    before = h.values.sum(axis=2)
    after = h.values[:, :, 1] + h.values[:, :, 3]          # cent_mb_decision == 1
    f, axs = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    f.subplots_adjust(left=.07, right=.92, bottom=.13, top=.95, wspace=.05)
    vmax = before.max()
    norm = LogNorm(max(vmax * 1e-5, 1e-12), vmax)
    for a, v, title in zip(axs, (before, after), ('Before: no classifier requirement', 'After: minimum-bias pass')):
        image = a.pcolormesh(h.edges[1], h.edges[0], np.ma.masked_less_equal(v, 0), norm=norm, cmap=R.BIRD, shading='flat', rasterized=True)
        a.set_xlim(0, 80); a.set_xlabel('Centrality percentile [%]')
        a.text(.03, .96, title, transform=a.transAxes, va='top', fontsize=11, weight='bold', bbox={'facecolor': 'white', 'edgecolor': 'none'})
    axs[0].set_ylabel(r'$E^{\mathrm{total}}_{\mathrm{calo}}$ [GeV]')
    cax = f.add_axes([.93, .13, .015, .82]); f.colorbar(image, cax=cax).set_label('Weighted events / bin')
    axs[1].text(.97, .86, 'sPHENIX Internal\nAu+Au embedded simulation\n' + energy('auau'), transform=axs[1].transAxes, ha='right', va='top', fontsize=10)
    save(f, store, name)


# ----------------------------------------------------------------------------------------------------------------
# event and dataset figures

def centrality_data_vs_mc(store):
    name = 'simulation_stitching/auau_centrality_data_embeddedMC.png'
    # the DATA target of the centrality weights (TreeToHists/centrality.py): the nominal trigger's events from the QA fills
    values, edges, _ = store.target()
    data = Hist(values, values.copy(), [edges]) if values is not None and values.sum() > 0 else None
    NOTES[name] = dict(data_trigger=CENTRALITY.DATA_TRIGGER, data_events=float(values.sum()) if values is not None else 0.0)
    sims = {fam: store.sim('auau', fam, 'evt_cent_%s' % fam) for fam in ('phosim', 'jetsim')}
    if data is None and all(v is None for v in sims.values()):
        return missing(name, 'no DATA target events and no simulation QA')
    f, a = R.axes()
    def two_percent(v, edges):
        k = max(1, int(round(2.0 / np.diff(edges)[0])))
        n = len(v) // k * k
        return v[:n].reshape(-1, k).sum(axis=1), edges[:n + 1:k]
    if data is not None:
        v, e = two_percent(data.values, data.edges[0])
        m = e[:-1] < 80
        R.points(a, ((e[:-1] + e[1:]) / 2)[m], unit(v[m])[0], None, R.BLACK, 'D', 'DATA (weight-contract reference)')
    for (fam, h), color, marker, text in zip(sims.items(), (R.BLUE, R.ORANGE), ('s', '^'), ('Embedded photon+jet', 'Embedded inclusive jet')):
        if h is None:
            continue
        v, _, edges = project(h, 0, {1: (0, 1)})
        v, e = two_percent(v, edges)
        m = e[:-1] < 80
        R.points(a, ((e[:-1] + e[1:]) / 2)[m] + (.3 if fam == 'jetsim' else 0), unit(v[m])[0], None, color, marker, text)
    a.set_xlim(0, 80); a.set_xlabel('Centrality percentile [%]'); a.set_ylabel('Fraction per 2% bin')
    R.brand(a); text_block(a, ['Au+Au, ' + energy('auau'), 'simulation before centrality reweighting'], y=.88)
    a.legend(loc='upper right', fontsize=9.5)
    save(f, store, name)


def event_basics(store):
    """event_obs/*.pdf: Au+Au 0-80%, DATA vs photon+jet and inclusive-jet simulation."""
    def all_classes(fetch):
        return summed_classes(store, 'auau', fetch)
    # centrality
    name = 'event_obs/centrality.pdf'
    data = store.data('auau', 'evt_cent_data')
    if data is None:
        missing(name, 'no Au+Au DATA QA')
    else:
        f, a = R.axes()
        v, _, e = project(data, 0, {1: (0, 1)})
        R.points(a, (e[:-1] + e[1:]) / 2, unit(v)[0], None, R.BLACK, 'o', 'DATA')
        for fam, color, text in (('phosim', R.RED, 'Photon+jet simulation'), ('jetsim', R.BLUE, 'Inclusive-jet simulation')):
            h = store.sim('auau', fam, 'evt_cent_%s' % fam)
            if h is not None:
                v, _, e = project(h, 0, {1: (1, 2)})
                a.stairs(unit(v)[0], e, color=color, lw=1.4, label=text + ' (reweighted)')
        a.set_xlim(0, 100); a.set_xlabel('Centrality [%]'); a.set_ylabel('Normalized events')
        R.brand(a); text_block(a, ['Au+Au, ' + energy('auau')], y=.88); a.legend(loc='upper right', fontsize=9)
        save(f, store, name)
    # vertex
    name = 'event_obs/vertexz.pdf'
    data = all_classes(lambda c: store.data('auau', 'evt_vz_data_%s' % c))
    if data is None:
        missing(name, 'no Au+Au DATA QA')
    else:
        f, a = R.axes()
        R.points(a, (data.edges[0][:-1] + data.edges[0][1:]) / 2, unit(data.values)[0], None, R.BLACK, 'o', 'DATA')
        for fam, color, text in (('phosim', R.RED, 'Photon+jet simulation'), ('jetsim', R.BLUE, 'Inclusive-jet simulation')):
            h = all_classes(lambda c: store.sim('auau', fam, 'evt_vz_%s_%s' % (fam, c)))
            if h is not None:
                a.stairs(unit(h.values)[0], h.edges[0], color=color, lw=1.4, label=text)
        a.set_xlabel(r'$v_z$ [cm]'); a.set_ylabel('Normalized events'); a.set_xlim(-30, 30)
        R.brand(a); text_block(a, ['Au+Au 0–80%, ' + energy('auau')], y=.88); a.legend(loc='upper right', fontsize=9)
        save(f, store, name)
    # cluster pT (DATA counts, simulation normalised to the DATA integral), after preselection
    name = 'event_obs/cluster_pt.pdf'
    data = all_classes(lambda c: store.labels('auau', 'data', 'clu_et', c))
    if data is None:
        missing(name, 'no Au+Au DATA QA')
    else:
        f, a = R.axes()
        dv, dvar, e = project(data, 0, {1: (1, 2)})
        x = (e[:-1] + e[1:]) / 2
        R.points(a, x, dv, np.sqrt(dvar), R.BLACK, 'o', 'DATA')
        for fam, color, text in (('phosim', R.RED, 'Photon+jet simulation'), ('jetsim', R.BLUE, 'Inclusive-jet simulation')):
            h = all_classes(lambda c: store.labels('auau', fam, 'clu_et', c))
            if h is not None:
                v, _, _ = project(h, 0, {1: (1, 2)})
                a.stairs(v * (dv.sum() / v.sum() if v.sum() > 0 else 0), e, color=color, lw=1.4, label=text)
        a.set_yscale('log'); a.set_xlabel(r'Cluster $E_T$ [GeV]'); a.set_ylabel('Clusters'); a.set_xlim(5, 50)
        R.brand(a); text_block(a, ['Au+Au 0–80%, ' + energy('auau'), 'preselected candidates'], y=.88); a.legend(loc='upper right', fontsize=9)
        save(f, store, name)
    # eta-phi maps
    for tag, fam, text in (('data', 'data', 'DATA'), ('photon', 'phosim', 'Photon+jet simulation')):
        name = 'event_obs/cluster_eta_phi_%s.pdf' % tag
        h = all_classes(lambda c: store.labels('auau', fam, 'clu_etaphi', c))
        if h is None:
            missing(name, 'no Au+Au %s QA' % fam)
            continue
        m = h.values[:, :, 1]
        f, a = R.axes(matrix=True)
        image = a.pcolormesh(h.edges[0], h.edges[1], np.ma.masked_less_equal(m.T, 0), cmap=R.BIRD, shading='flat', rasterized=True)
        a.set_xlabel(r'Cluster $\eta$'); a.set_ylabel(r'Cluster $\phi$ [rad]')
        R.brand(a, box=True); text_block(a, [text + ', Au+Au 0–80%', r'$E_T^{\mathrm{cluster}}>10$ GeV, preselected'], y=.88)
        f.colorbar(image, ax=a, pad=.025, fraction=.045).set_label('Clusters')
        save(f, store, name)
    # jets
    name = 'event_obs/jet_pt.pdf'
    data = all_classes(lambda c: store.data('auau', 'jet_pt_data_%s' % c))
    if data is None:
        missing(name, 'no Au+Au DATA QA')
    else:
        f, a = R.axes()
        dv, dvar, e = project(data, 0, {1: (1, 2)})
        R.points(a, (e[:-1] + e[1:]) / 2, dv, np.sqrt(dvar), R.BLACK, 'o', 'DATA')
        for fam, color, text in (('phosim', R.RED, 'Photon+jet simulation'), ('jetsim', R.BLUE, 'Inclusive-jet simulation')):
            h = all_classes(lambda c: store.sim('auau', fam, 'jet_pt_%s_%s' % (fam, c)))
            if h is not None:
                v, _, _ = project(h, 0, {1: (1, 2)})
                a.stairs(v * (dv.sum() / v.sum() if v.sum() > 0 else 0), e, color=color, lw=1.4, label=text)
        a.set_yscale('log'); a.set_xlabel(r'Jet $p_T$ [GeV]'); a.set_ylabel('Jets'); a.set_xlim(0, 60)
        R.brand(a); text_block(a, ['Au+Au 0–80%, ' + energy('auau'), r'anti-$k_T$ $R=0.4$, events with a Region-A photon'], y=.88)
        a.legend(loc='upper right', fontsize=9)
        save(f, store, name)
    for tag, fam, text in (('data', 'data', 'DATA'), ('photon', 'phosim', 'Photon+jet simulation')):
        name = 'event_obs/jet_eta_phi_%s.pdf' % tag
        h = all_classes(lambda c: store.data('auau', 'jet_etaphi_data_%s' % c) if fam == 'data'
                        else store.sim('auau', fam, 'jet_etaphi_%s_%s' % (fam, c)))
        if h is None:
            missing(name, 'no Au+Au %s QA' % fam)
            continue
        m = h.values[:, :, 1]
        f, a = R.axes(matrix=True)
        image = a.pcolormesh(h.edges[0], h.edges[1], np.ma.masked_less_equal(m.T, 0), cmap=R.BIRD, shading='flat', rasterized=True)
        a.set_xlabel(r'Jet $\eta$'); a.set_ylabel(r'Jet $\phi$ [rad]')
        R.brand(a, box=True); text_block(a, [text + ', Au+Au 0–80%', r'$p_T^{\mathrm{jet}}\geq5$ GeV, events with a Region-A photon'], y=.88)
        f.colorbar(image, ax=a, pad=.025, fraction=.045).set_label('Jets')
        save(f, store, name)


def smooth_fit(x, y, e):
    """Descriptive fit exp[a - (b + c ln x + d x) ln x] (x in GeV), weighted least squares on ln y; returns the model."""
    coef = numerical('stitching_fit', QA.stitching_fit, x, y, e)
    if not np.isfinite(coef).all(): return None
    return lambda z: np.exp(coef[0] - (coef[1] + coef[2] * np.log(z) + coef[3] * z) * np.log(z))


def stitching(store, system, family, kind):
    """kind 'generator': gen_<F> (leading prompt photon for the photon family, leading truth jet for the inclusive
    family); 'reco': Region-A leaders of 0-20% (Au+Au) or all (pp), every label."""
    fam_tag = 'photonjet' if family == 'phosim' else 'inclusivejet'
    name = 'simulation_stitching/%s_%s_%s_stitching.png' % (system, fam_tag, 'generator' if kind == 'generator' else 'reco_tag')
    c = 'all' if system == 'pp' else 'cent0_20'
    if kind == 'generator':
        parts = store.sim(system, family, 'gen_%s' % family, per_sample=True)
        parts = {s: Hist(h.values[:, 0 if family == 'phosim' else 1], h.variances[:, 0 if family == 'phosim' else 1], [h.edges[0]])
                 for s, h in parts.items()}
        xlabel = 'Maximum truth photon $p_T$ [GeV]' if family == 'phosim' else 'Maximum truth jet $p_T$ [GeV]'
        lo, hi = (15, 40) if family == 'phosim' else (10, 60)
    else:
        parts = {}
        for label in LABELS[family]:
            for s, h in store.sim(system, family, 'abcd_%s_%s' % (label, c), per_sample=True).items():
                parts[s] = Hist(h.values[0], h.variances[0], [h.edges[1]]) + parts.get(s)
        xlabel, (lo, hi) = r'Reconstructed photon $E_T$ [GeV]', (15, 35)
    if not parts:
        return missing(name, 'no %s %s simulation QA' % (system, family))
    edges = next(iter(parts.values())).edges[0]
    width = np.diff(edges)
    x = (edges[:-1] + edges[1:]) / 2
    m = window(edges, lo, hi)
    total = sum(parts.values(), None)
    f, a, b = R.axes(ratio=True)
    palette = [R.ORANGE, R.PURPLE, R.GREEN, R.RED, R.BLUE, '#8D6E63', '#00838F', '#7CB342']
    for (sample, h), color in zip(sorted(parts.items()), palette):
        R.points(a, x[m], (h.values / width)[m], (np.sqrt(h.variances) / width)[m], color, 'o', sample, open_=True)
    y, e = total.values / width, np.sqrt(total.variances) / width
    R.points(a, x[m], y[m], e[m], R.BLACK, 'D', 'Stitched sum')
    model = smooth_fit(x[m], y[m], e[m])
    if model is not None:
        xs = np.linspace(lo, hi, 200)
        a.plot(xs, model(xs), color=R.BLACK, lw=1.0, ls=':')
        ratio = y[m] / model(x[m])
        good = np.isfinite(ratio) & (y[m] > 0)
        R.points(b, x[m][good], ratio[good], (e[m] / model(x[m]))[good], R.BLACK, 'o')
        rms = float(np.sqrt(np.mean((ratio[good] - 1) ** 2))) if good.any() else float('nan')
        b.text(.97, .85, 'RMS %.1f%%' % (100 * rms), transform=b.transAxes, ha='right', fontsize=9.5)
        NOTES[name] = dict(rms=rms)
    b.axhline(1, ls='--', color=R.BLACK, lw=.8); b.set_ylim(.5, 1.5); b.set_ylabel('Spectrum / fit', fontsize=11)
    b.set_xlabel(xlabel); b.set_xlim(lo, hi)
    a.set_yscale('log'); a.set_ylabel('Weighted spectrum / GeV')
    R.brand(a)
    text_block(a, [('pp PYTHIA 8, ' if system == 'pp' else 'Au+Au embedded PYTHIA 8, ') + energy(system),
                   'Photon+jet samples' if family == 'phosim' else 'Inclusive-jet samples',
                   'Generator truth' if kind == 'generator' else ('Leading tight + isolated photon tag' + ('' if system == 'pp' else ', 0–20%'))], y=.88)
    a.legend(loc='upper right', fontsize=9)
    save(f, store, name)


# ----------------------------------------------------------------------------------------------------------------
# working-point derivations (checks of the frozen packages)

def thresholds(h, efficiencies, pt_groups, slice_range, above=True):
    """Signal-efficiency thresholds per photon-pT group: the score (above=True) or isolation (above=False) value at
    which the requested fraction of the signal passes, by linear interpolation of the cumulative distribution."""
    out = {}
    for lo, hi in pt_groups:
        v, _, edges = project(h, 0, {1: (lo, hi), 2: slice_range})
        total = v.sum()
        if total <= 0:
            continue
        levels = numerical('quantiles', QA.quantiles, v, edges, efficiencies, above)
        for eff, level in zip(efficiencies, levels): out[(lo, hi, eff)] = float(level)
    return out


def working_points(store, what):
    """what 'bdt': auau_h70_bdt_pt_centrality_grid / _centrality_fits; 'iso': auau_r03_isolation_*."""
    stem = 'bdt_sig' if what == 'bdt' else 'isoslice_sig'
    names = (('photonOnly/auau_h70_bdt_pt_centrality_grid.png', 'photonOnly/auau_h70_bdt_centrality_fits.png') if what == 'bdt'
             else ('isolation/auau_r03_isolation_pt_centrality_grid.png', 'isolation/auau_r03_isolation_centrality_fits.png'))
    h = store.sim('auau', 'phosim', stem)
    if h is None:
        for n in names:
            missing(n, 'no Au+Au photon-family simulation QA')
        return
    groups = [(15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 30), (30, 35)]
    effs = (.7, .8, .9)
    colors = {.7: R.BLUE, .8: R.GREEN, .9: R.ORANGE}
    per_slice = {s: thresholds(h, effs, groups, (s, s + 1), above=(what == 'bdt')) for s in range(16)}
    f, axs = plt.subplots(4, 4, figsize=(14, 11), sharex=True, sharey='row')
    f.subplots_adjust(left=.07, right=.99, bottom=.07, top=.95, hspace=.28, wspace=.08)
    summary = {e: [] for e in effs}
    for s in range(16):
        a = axs[s // 4][s % 4]
        R.style(a); a.set_title('%d–%d%%' % (5 * s, 5 * s + 5), fontsize=10)
        for e in effs:
            pts = [(0.5 * (lo + hi), per_slice[s].get((lo, hi, e))) for lo, hi in groups]
            pts = [(x, y) for x, y in pts if y is not None]
            if not pts:
                summary[e].append(np.nan)
                continue
            xs, ys = zip(*pts)
            a.plot(xs, ys, 'o', color=colors[e], ms=3.5)
            level = float(np.median(ys))
            a.axhline(level, ls='--', color=colors[e], lw=.9)
            summary[e].append(level)
    for a in axs[-1]:
        a.set_xlabel(r'Photon $p_T$ [GeV]')
    f.text(.01, .5, 'BDT score threshold' if what == 'bdt' else r'$R=0.3$ isolation threshold [GeV]', rotation=90, va='center', fontsize=13)
    axs[0][3].legend(handles=[plt.Line2D([], [], color=colors[e], marker='o', ls='--', label='%d%%' % round(100 * e)) for e in effs], fontsize=8)
    save(f, store, names[0])
    f, a = R.axes()
    centers = 5 * np.arange(16) + 2.5
    lines = []
    for e in effs:
        y = np.array(summary[e])
        good = np.isfinite(y)
        if good.sum() < 2:
            continue
        slope, intercept = numerical('centrality_fit', QA.centrality_fit, centers, y)
        a.plot(centers[good], y[good], 'o', color=colors[e])
        a.plot([0, 80], [intercept, intercept + 80 * slope], color=colors[e], lw=1.4,
               label='%d%%: %.4f %+.5f C' % (round(100 * e), intercept, slope))
        lines.append(dict(efficiency=e, intercept=float(intercept), slope=float(slope)))
    a.set_xlim(0, 80); a.set_xlabel('Centrality percentile [%]')
    a.set_ylabel('BDT score threshold' if what == 'bdt' else r'$R=0.3$ isolation threshold [GeV]')
    R.brand(a); text_block(a, ['Au+Au, ' + energy('auau') + r', $15\leq p_T^\gamma<35$ GeV', 'check of the frozen package: descriptive lines'], y=.88)
    a.legend(loc='lower right' if what == 'bdt' else 'upper right', fontsize=9)
    NOTES[names[1]] = dict(lines=lines)
    save(f, store, names[1])


TRIGGER_STUDIES = {
    # system: (reference trigger name, photon triggers (label, name prefix), first run, output names)
    'auau': ('MBD N&S >= 2, vtx < 150 cm', (('Photon 10', 'Photon 10 GeV + MBD NS >= 2'), ('Photon 12', 'Photon 12 GeV + MBD NS >= 2')), 0,
             ('event_obs/hMaxClusterEnergy_groupOverlay.png', 'event_obs/hMaxClusterEnergy_groupTurnOnOverlay.png')),
    'pp': ('MBD N&S >= 1', (('Photon 3', 'Photon 3 GeV + MBD NS >= 1'), ('Photon 4', 'Photon 4 GeV + MBD NS >= 1'),
                            ('Photon 5', 'Photon 5 GeV + MBD NS >= 1')), 47289,
           ('event_obs/pp_max_cluster_energy_overlay.png', 'event_obs/pp_max_cluster_energy_turnon.png')),
}


def bit_of(names, prefix, prefer='vtx < 150 cm'):
    """The bit of a trigger by name: the exact name, else the names starting with the prefix (preferring `prefer`)."""
    if prefix in names:
        return names.index(prefix)
    hits = [i for i, n in enumerate(names) if n.startswith(prefix)]
    if not hits:
        return None
    preferred = [i for i in hits if prefer in names[i]]
    return (preferred or hits)[0]


def trigger_study(store, system):
    """Maximum-cluster-energy spectra of the reference and photon-trigger samples (scaled bits) and the
    reference-based turn-ons (reference scaled and photon trigger live), over the runs in which every trigger of the
    study is enabled (initial prescale >= 0) and the reference was recorded, from the run's own trigger table."""
    reference, photons, first_run, names = TRIGGER_STUDIES[system]
    spectra, turnon_num, turnon_den, runs, edges = defaultdict(lambda: 0), defaultdict(lambda: 0), 0, [], None
    skipped = defaultdict(int)
    for run, info, h in store.trigger_rows(system):
        if run < first_run:
            skipped['before run %d' % first_run] += 1
            continue
        bits = [bit_of(info['name'], reference)] + [bit_of(info['name'], p) for _, p in photons]
        if any(b is None for b in bits) or any(info['initial_prescale'][b] < 0 for b in bits):
            skipped['a trigger of the study disabled'] += 1
            continue
        ref, pho = bits[0], bits[1:]
        if edges is not None and not np.array_equal(edges, h.edges[0]):
            raise ValueError('trigger QA source axes differ')
        edges = h.edges[0]
        if h.values[:, ref, 64].sum() <= 0:
            skipped['no reference events'] += 1
            continue
        runs.append(run)
        spectra['ref'] = spectra['ref'] + h.values[:, ref, 64]
        turnon_den = turnon_den + h.values[:, ref, 64]
        for (label, _), b in zip(photons, pho):
            spectra[label] = spectra[label] + h.values[:, b, 64]
            turnon_num[label] = turnon_num[label] + h.values[:, ref, b]
    NOTES[names[0]] = NOTES[names[1]] = dict(runs=len(runs), first=min(runs) if runs else None, last=max(runs) if runs else None,
                                             skipped=dict(skipped), reference=reference, photons=[p for _, p in photons])
    if not runs:
        for n in names:
            missing(n, 'no run with every trigger of the study enabled')
        return
    x = (edges[:-1] + edges[1:]) / 2
    colors = [R.RED, R.BLUE, R.GREEN]
    f, a = R.axes()
    a.stairs(spectra['ref'], edges, color=R.BLACK, lw=1.5, label=reference)
    for (label, _), color in zip(photons, colors):
        a.stairs(spectra[label], edges, color=color, lw=1.5, label=label)
    a.set_yscale('log'); a.set_xlabel('Maximum cluster energy [GeV]'); a.set_ylabel('Scaled-trigger events')
    R.brand(a); text_block(a, [('pp, ' if system == 'pp' else 'Au+Au, ') + energy(system), 'scaled-trigger study, %d runs' % len(runs)], y=.88)
    a.legend(loc='upper right', fontsize=9.5)
    save(f, store, names[0])
    f, a = R.axes()
    for (label, _), color in zip(photons, colors):
        eff = np.divide(turnon_num[label], turnon_den, out=np.full(len(x), np.nan), where=turnon_den > 0)
        err = np.sqrt(np.clip(eff * (1 - eff), 0, None) / np.where(turnon_den > 0, turnon_den, 1))
        crossing = next((float(xx) for xx, ee in zip(x, eff) if np.isfinite(ee) and ee >= .9), None)
        R.points(a, x, eff, err, color, 'o', label + ' / MBD' + (', 90%% at %.1f GeV' % crossing if crossing else ''))
    a.axhline(1, ls='--', color=R.BLACK, lw=.8); a.set_ylim(0, 1.15); a.set_xlim(0, 30 if system == 'auau' else 20)
    a.set_xlabel('Maximum cluster energy [GeV]'); a.set_ylabel('Trigger / MBD efficiency')
    R.brand(a); text_block(a, [('pp, ' if system == 'pp' else 'Au+Au, ') + energy(system), 'scaled-trigger study, %d runs' % len(runs)], y=.88)
    a.legend(loc='lower right', fontsize=9.5)
    save(f, store, names[1])


FIGURES = {
    'triggers': lambda s: (trigger_study(s, 'auau'), trigger_study(s, 'pp')),
    'ncb': lambda s: (ncb_score(s, 'pp'), ncb_score(s, 'auau'), ncb_eta_time(s)),
    'shapes': shower_shapes, 'widths': widths, 'bdt_per_cent': bdt_per_cent, 'e11_flow': e11_flow,
    'truth_isolation': truth_isolation, 'truth_pt': truth_pt_after_isolation, 'iso_overlays': isolation_overlays,
    'iso_stacks': isolation_stacks, 'abcd': abcd_counts, 'purity': purity, 'purity_closure': purity_closure,
    'efficiency': efficiency, 'minbias': minbias_pretraining, 'centrality': centrality_data_vs_mc, 'basics': event_basics,
    'stitching': lambda s: [stitching(s, sys_, fam, kind) for sys_ in ('pp', 'auau') for fam in ('phosim', 'jetsim')
                            for kind in ('generator', 'reco')],
    'bdt_wp': lambda s: working_points(s, 'bdt'), 'iso_wp': lambda s: working_points(s, 'iso'),
}


def prepare(manifest, config, output, systems=('pp', 'auau')):
    """Save ordinary histogram arrays and numerical QA results once.

    The existing figure orchestration requests exactly the projections used by
    the plots. During preparation its save operation only closes the temporary
    canvas; no figure is published. Rendering reads this package without ROOT.
    """
    global PREPARING
    output = Path(output)
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('QA input package already exists')
    NUMERICS.clear(); NOTES.clear()
    store = Store(manifest, config, {}, False)
    store.collection = {key: item for key, item in store.collection.items()
                        if {'PP': 'pp', 'AuAu': 'auau'}[item['metadata']['system']] in systems}
    PREPARING = True
    try:
        for key, fn in FIGURES.items():
            try: fn(store)
            except KeyError as exc: NOTES[key] = 'not drawn: ' + str(exc)
    finally:
        PREPARING = False
        plt.close('all')
    arrays = {}; groups = {}; triggers = {}
    def pack(h):
        token = 'h' + str(len(groups)) + '_' + str(len(arrays))
        arrays[token + '_value'] = h.values; arrays[token + '_variance'] = h.variances
        for i, edges in enumerate(h.edges): arrays[token + '_edge' + str(i)] = edges
        return dict(token=token, dimensions=len(h.edges))
    for group, histograms in store.cache.items():
        groups[group] = {name: pack(h) for name, h in histograms.items()}
    for system in ('pp', 'auau'):
        triggers[system] = [dict(run=run, info=info, histogram=pack(h)) for run, info, h in store.trigger_rows(system)]
    values, edges, records = store.target()
    target = None if values is None else pack(Hist(values, values.copy(), [edges]))
    for key, value in NUMERICS.items(): arrays[key] = value
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as stream: np.savez_compressed(stream, **arrays)
    metadata = dict(schema='PhotonJetQAInputsV1', sha256=R.sha(output), groups=groups, triggers=triggers,
                    centrality_target=target, centrality_records=records, config=config, notes=NOTES,
                    numerical_keys=sorted(NUMERICS),
                    collection={key: dict(v, path=str(v['path'])) for key, v in store.collection.items()},
                    inputs={str(v['path']): v['sha256'] for v in store.collection.values()})
    with output.with_suffix('.json').open('x') as stream:
        json.dump(metadata, stream, allow_nan=False, sort_keys=True)
    return metadata


def saved_store(path):
    """Restore only saved arrays; source histogram files are never reopened."""
    path = Path(path); metadata = json.loads(path.with_suffix('.json').read_text())
    if metadata.get('schema') != 'PhotonJetQAInputsV1' or R.sha(path) != metadata['sha256']:
        raise ValueError('saved QA figure inputs changed')
    with np.load(path, allow_pickle=False) as f: arrays = {k: f[k] for k in f.files}
    def unpack(row):
        token = row['token']
        return Hist(arrays[token + '_value'], arrays[token + '_variance'],
                    [arrays[token + '_edge' + str(i)] for i in range(row['dimensions'])])
    store = object.__new__(Store); store.root = path; store.config = metadata['config']; store.unit = False; store.saved = True
    store.override = {}; store.collection = {}
    store.cache = {group: {name: unpack(h) for name, h in histograms.items()} for group, histograms in metadata['groups'].items()}
    store.trigger_records = {system: [(row['run'], row['info'], unpack(row['histogram'])) for row in rows]
                             for system, rows in metadata['triggers'].items()}
    h = unpack(metadata['centrality_target']) if metadata['centrality_target'] else None
    store.centrality_target = (None if h is None else h.values, None if h is None else h.edges[0], metadata['centrality_records'])
    # Saved sample metadata supports the same single application of the weight.
    store.collection = metadata.get('collection', {})
    NUMERICS.clear(); NUMERICS.update({key: arrays[key] for key in metadata['numerical_keys']})
    return store, metadata


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', type=Path, required=True, help='saved QA figure inputs from analysis preparation')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args(argv)
    R.OUT = args.output.resolve()
    R.FILES.clear()      # this manifest lists this renderer's figures only
    NOTES.clear()
    store, metadata = saved_store(args.inputs)
    for key, fn in FIGURES.items():
        if args.only and key not in args.only:
            continue
        try:
            fn(store)
        except KeyError as exc:
            NOTES[key] = 'not drawn: %s' % exc
            print('SKIP', key, exc, flush=True)
    R.OUT.mkdir(parents=True, exist_ok=True)
    (R.OUT / 'QA_FIGURE_MANIFEST.json').write_text(json.dumps({
        'renderer': str(Path(__file__).resolve()), 'renderer_sha256': R.sha(__file__), 'inputs_sha256': R.sha(args.inputs),
        'figures': R.FILES, 'notes': NOTES}, indent=1, default=str) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
