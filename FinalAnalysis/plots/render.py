#!/usr/bin/env python3
"""Render the IAN γ+jet figures from canonical figure packages (FinalAnalysis/plots/figure_package.py), one per result.

Port of the live drop-in renderer (IAN_GammaJet_DropIn_20260925/reproducibility/analysis_source/render_package.py):
the presentation (style, branding, axes, colours, file names) is kept; every figure reads the canonical package, so
binning, reported range, iterations and classes come from the result, not from constants. Nothing is computed here
beyond sums over photon-pT bins, bin-width division and linear error propagation of recorded covariances.

    python -m FinalAnalysis.plots.render --packages DIR --output DIR [--atlas FinalAnalysis/plots/inputs/atlas_reference.json]
"""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm  # noqa: E402
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, MultipleLocator  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, GREEN, RED, BLACK, PURPLE = '#0072B2', '#D55E00', '#009E73', '#C62828', '#111111', '#6A3D9A'
BIRD = LinearSegmentedColormap.from_list('root_bird_like',
    ['#352A87', '#1552AD', '#1688BA', '#25A5AD', '#56B99A', '#A9D17A', '#F2E955', '#FFF99A'])
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11.5, 'mathtext.fontset': 'dejavusans',
    'axes.labelsize': 15, 'axes.labelpad': 7, 'axes.linewidth': 1.05,
    'legend.fontsize': 10.3, 'legend.frameon': False, 'figure.facecolor': 'white',
    'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'xtick.direction': 'in',
    'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 6, 'ytick.major.size': 6, 'xtick.minor.size': 3,
    'ytick.minor.size': 3, 'axes.formatter.use_mathtext': True,
})
FILES = []
OUT = None
# result label -> (file stem of the drop-in names, response stem, legend text)
TAGS = {'pp': ('pp', 'pp', 'pp'), 'auau_cent0_20': ('auau020', 'auau_0to20', 'Au+Au 0–20%'),
        'auau_cent20_50': ('auau2050', 'auau_20to50', 'Au+Au 20–50%'), 'auau_cent50_80': ('auau5080', 'auau_50to80', 'Au+Au 50–80%')}
CLASS_COLORS = {'auau_cent0_20': RED, 'auau_cent20_50': PURPLE, 'auau_cent50_80': GREEN}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(fig, name):
    p = OUT / 'fig' / name
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=300)
    FILES.append({'path': str(p.relative_to(OUT)), 'sha256': sha(p),
                  'width': round(fig.get_figwidth() * 300), 'height': round(fig.get_figheight() * 300)})
    plt.close(fig)
    print('PNG', name, flush=True)


def style(a):
    a.xaxis.set_minor_locator(AutoMinorLocator())
    a.yaxis.set_minor_locator(AutoMinorLocator())
    a.tick_params(which='both', top=True, right=True)


def axes(ratio=False, matrix=False):
    if ratio:
        f, (a, b) = plt.subplots(2, 1, figsize=(7.1, 6.1), sharex=True, gridspec_kw={'height_ratios': [3.15, 1], 'hspace': .035})
        f.subplots_adjust(left=.145, right=.975, bottom=.12, top=.975)
        style(a); style(b)
        b.yaxis.set_major_locator(MaxNLocator(4))
        return f, a, b
    f, a = plt.subplots(figsize=(7.2, 6.2) if matrix else (7.1, 5.4))
    f.subplots_adjust(left=.145, right=.85 if matrix else .975, bottom=.15, top=.975)
    style(a)
    return f, a


def brand(a, size=16, box=False, x=.035):
    pack = HPacker(children=[TextArea('sPHENIX', textprops={'size': size, 'weight': 'bold', 'style': 'italic'}),
                             TextArea('Internal', textprops={'size': size - 1})], align='baseline', pad=0, sep=5)
    anchor = AnchoredOffsetbox(loc='upper left', child=pack, frameon=box, bbox_to_anchor=(x, .978),
                               bbox_transform=a.transAxes, borderpad=0, pad=.16 if box else 0)
    if box:
        anchor.patch.set_facecolor('white'); anchor.patch.set_edgecolor('none')
    a.add_artist(anchor)


def system_text(label):
    return r'$p+p,\ \sqrt{s}=200$ GeV' if label == 'pp' else TAGS[label][2] + r', $\sqrt{s_{NN}}=200$ GeV'


def label(a, p, groups=None, sim=False, extra=None, box=False):
    brand(a, box=box)
    lines = [system_text(str(p['label']))]
    if groups is not None:
        lines.append(rf'${groups[0]:g}<p_T^\gamma<{groups[-1]:g}$ GeV')
    if sim:
        lines.append('PYTHIA8 simulation' if str(p['label']) == 'pp' else 'Embedded simulation')
    if extra:
        lines += extra
    a.text(.04, .882, '\n'.join(lines), transform=a.transAxes, va='top', fontsize=10.9, linespacing=1.34,
           bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 2} if box else None)


def jet_text(p):
    return r'$R=0.4,\ p_T^{\mathrm{jet}}\geq%g$ GeV' % float(p['jet_pt_min_gev'])


def points(a, x, y, e, color=BLACK, marker='o', text=None, xerr=None, open_=False):
    return a.errorbar(x, y, yerr=e, xerr=xerr, fmt=marker, linestyle='none', ms=4.8, mew=1.1, mfc='white' if open_ else color,
                      mec=color, color=color, elinewidth=.95, capsize=0, label=text, zorder=4)


def centers(p):
    e = p['xj_edges']
    return (e[:-1] + e[1:]) / 2, np.diff(e)


def mask(p):
    return p['reported_mask'].astype(bool)


def summed(p, key, variance_key=None, cov_key=None):
    """An (n_pt, n_xj) array summed over photon pT, divided by the xJ width, with its error."""
    _, w = centers(p)
    y = p[key].sum(axis=0) / w
    if cov_key is not None:
        n_pt, n_xj = p[key].shape
        proj = np.kron(np.ones((1, n_pt)), np.eye(n_xj))
        e = np.sqrt(np.maximum(np.diag(proj @ p[cov_key] @ proj.T), 0)) / w
    elif variance_key is not None:
        e = np.sqrt(np.maximum(p[variance_key].sum(axis=0), 0)) / w
    else:
        e = np.zeros_like(y)
    return y, e


def plot_spec(a, p, y, e, color, marker, text, open_=False):
    x, w = centers(p)
    m = mask(p)
    points(a, x[m], y[m], e[m], color, marker, text, w[m] / 2, open_)


def step_spec(a, p, y, color, text, fill=False):
    m = mask(p); i = np.flatnonzero(m)
    edges = p['xj_edges'][i[0]:i[-1] + 2]
    if fill:
        a.stairs(y[m], edges, fill=True, color=color, alpha=.18, lw=0, zorder=1)
    a.stairs(y[m], edges, color=color, lw=1.5, label=text, zorder=2)


def drawx(a):
    a.set_xlim(0, 1.5)
    a.set_xticks(np.arange(0, 1.51, .3))
    a.set_xlabel(r'$x_{J\gamma}$')


def room(a, ys, errors=None, fraction=.69):
    hi = max(float(np.max(y)) for y in ys)
    lo = min(0., min(float(np.min(y)) for y in ys))
    if errors is not None:
        hi = max(float(np.max(y + e)) for y, e in zip(ys, errors))
        lo = min(lo, min(float(np.min(y - e)) for y, e in zip(ys, errors)))
    a.set_ylim(min(0, lo * 1.15), max(hi, 1e-10) / fraction)


def groups(p):
    return p['pt_edges']


def correction(p):
    stem = TAGS[str(p['label'])][0]
    m = mask(p)
    A, ae = summed(p, 'region_a', 'region_a_sumw2')
    C, _ = summed(p, 'scaled_region_c')
    P, pe = summed(p, 'purity_corrected', 'purity_corrected_variance')
    Y, ye = summed(p, 'corrected', cov_key='corrected_cov')
    data_driven = 'combinatoric_purity_corrected_background' in p
    if data_driven:
        K, _ = summed(p, 'combinatoric_purity_corrected_background')
        before = P + K
        k_text = 'Data-driven combinatorial'
    else:
        K, _ = summed(p, 'template_combinatoric')
        before = P
        k_text = 'Scaled combinatorial (simulated template)'
    f, a = axes()
    plot_spec(a, p, A, ae, BLACK, 'o', 'Region A')
    step_spec(a, p, C, ORANGE, r'$\alpha C$')
    plot_spec(a, p, before, pe, BLUE, 'D', 'Leakage-corrected')
    label(a, p, groups(p), extra=[jet_text(p)])
    a.legend(loc='upper right', handlelength=1.6)
    drawx(a); a.set_ylabel(r'$dN_{\mathrm{pair}}/dx_{J\gamma}$')
    room(a, [A[m], C[m], before[m]], [ae[m], np.zeros(m.sum()), pe[m]], .65)
    save(f, 'purityAndCombSub/%s_xj_purity_correction.png' % stem)
    f, a = axes()
    plot_spec(a, p, before, pe, BLUE, 'D', 'Purity-corrected', True)
    step_spec(a, p, K, ORANGE, k_text, True)
    plot_spec(a, p, Y, ye, BLACK, 'o', 'After subtraction')
    label(a, p, groups(p), extra=[jet_text(p)])
    a.legend(loc='upper right', handlelength=1.6)
    drawx(a); a.set_ylabel(r'$dN_{\mathrm{pair}}/dx_{J\gamma}$')
    room(a, [before[m], Y[m], K[m]], [pe[m], ye[m], np.zeros(m.sum())], .65)
    save(f, 'purityAndCombSub/%s_xj_combinatoric_subtraction.png' % stem)


GREY = '#C4C4C4'
NOTES = {}


def conditional(matrix, misses):
    """P(reco | truth) with the unfolder's denominator, the row sum plus misses, so that each truth row sums to its
    efficiency; and that efficiency (NaN for a row with no truth)."""
    matrix, misses = np.asarray(matrix, dtype=float), np.asarray(misses, dtype=float)
    matched = matrix.sum(axis=1)
    denominator = matched + misses
    prob = np.divide(matrix, denominator[:, None], out=np.zeros_like(matrix), where=denominator[:, None] > 0)
    efficiency = np.divide(matched, denominator, out=np.full_like(matched, np.nan), where=denominator > 0)
    return prob, efficiency


def shared_norm(*arrays):
    """One log norm for every panel of a figure: the top from the data, the floor vmax / 10**decades with the decades
    set by the 5% quantile of the positive cells (2 to 4)."""
    positive = np.concatenate([np.asarray(x)[np.asarray(x) > 0].ravel() for x in arrays])
    vmax = float(positive.max())
    decades = int(np.clip(np.ceil(np.log10(vmax / np.quantile(positive, .05))), 2, 4))
    return LogNorm(vmax / 10 ** decades, vmax)


def mesh(a, x_edges, y_edges, prob, norm):
    """Three states: colour at or above the floor, grey where the cell holds entries below it, white where it holds
    none."""
    below = (prob > 0) & (prob < norm.vmin)
    if below.any():
        a.pcolormesh(x_edges, y_edges, np.ma.masked_where(~below, np.ones_like(prob)), cmap=ListedColormap([GREY]),
                     vmin=0, vmax=1, shading='flat', rasterized=True)
    return a.pcolormesh(x_edges, y_edges, np.ma.masked_less(prob, norm.vmin), norm=norm, cmap=BIRD, shading='flat',
                        rasterized=True)


def states_text(a, norm):
    """The three colour states, in the figure's bottom-right corner (clear of the axis labels)."""
    a.figure.text(.985, .012, 'grey: below %.0e;  white: empty' % norm.vmin, ha='right', va='bottom', fontsize=8.6,
                  color='#555555')


def merge_edges(edges, anchors):
    """Pairs of neighbouring bins between consecutive anchors (every anchor stays an edge), and the merged bin of every
    native bin."""
    edges = np.asarray(edges, dtype=float)
    keep = set(np.flatnonzero(np.isin(edges, np.asarray(anchors, dtype=float))).tolist()) | {0, len(edges) - 1}
    marks = sorted(keep)
    merged = [0]
    for lo, hi in zip(marks[:-1], marks[1:]):
        merged += list(range(lo + 2, hi, 2)) + [hi]
    merged = sorted(set(merged))
    return edges[merged], np.searchsorted(np.asarray(merged), np.arange(len(edges) - 1), side='right') - 1


def aggregate(matrix, rows, cols):
    out = np.zeros((rows.max() + 1, cols.max() + 1))
    np.add.at(out, (rows[:, None], cols[None, :]), matrix)
    return out


def native_photon(p):
    """The native photon response, merged by two (analysis edges kept) when the median effective entries of its core
    cells are below 30."""
    te, rex = p['native_truth_pt_edges'], p['native_reco_pt_edges']
    matrix, misses = p['native_photon_response'], p['native_photon_misses']
    note = dict(source=str(p['native_source']), merged=False, median_core_neff=None)
    if 'native_photon_neff' in p:
        prob, _ = conditional(matrix, misses)
        inside = np.zeros_like(prob, dtype=bool)
        lo, hi = p['pt_edges'][0], p['pt_edges'][-1]
        inside[(te[:-1] >= lo) & (te[1:] <= hi)] = True
        core = inside & (prob >= prob.max() * 1e-2)
        note['median_core_neff'] = float(np.median(p['native_photon_neff'][core])) if core.any() else 0.0
        if note['median_core_neff'] < 30:
            te, t_map = merge_edges(te, p['pt_edges'])
            rex, r_map = merge_edges(rex, p['pt_edges'])
            matrix = aggregate(matrix, t_map, r_map)
            misses = np.bincount(t_map, weights=misses, minlength=len(te) - 1)
            note['merged'] = True
    return te, rex, matrix, misses, note


def header(f, p, x):
    """Brand and sample line above the axes, so that no cell of a matrix is covered."""
    energy = r'$\sqrt{s}=200$ GeV' if str(p['label']) == 'pp' else r'$\sqrt{s_{NN}}=200$ GeV'
    f.text(x, .95, 'sPHENIX', fontsize=14, weight='bold', style='italic', va='bottom')
    f.text(x + .135, .95, 'Internal', fontsize=13, va='bottom')
    f.text(x, .905, TAGS[str(p['label'])][2] + ', ' + energy + ', ' + ('PYTHIA8' if str(p['label']) == 'pp' else 'embedded')
           + ' simulation', fontsize=10.2, va='bottom')


def photon_response(p):
    """Native-grid photon response (the IAN slot) and its analysis-bin twin on the same norm and colour bar."""
    if 'native_photon_response' not in p:
        return
    rstem = TAGS[str(p['label'])][1]
    te, rex, matrix, misses, note = native_photon(p)
    prob, eff = conditional(matrix, misses)
    window, window_eff = conditional(p['photon_response'], p['photon_misses'])
    norm = shared_norm(prob, window)
    NOTES['%s_photon_response' % rstem] = dict(note, floor=norm.vmin, vmax=norm.vmax)
    colour_label = r'$P(p_T^{\gamma,\mathrm{reco}}\mid p_T^{\gamma,\mathrm{truth}})$, rows sum to $\varepsilon$'
    f = plt.figure(figsize=(7.8, 6.6))
    a = f.add_axes([.115, .125, .585, .755]); strip = f.add_axes([.712, .125, .10, .755]); cax = f.add_axes([.86, .125, .024, .755])
    style(a)
    image = mesh(a, rex, te, prob, norm)
    rows = np.flatnonzero(np.isfinite(eff) & (prob.sum(axis=1) > 0)); cols = np.flatnonzero(prob.sum(axis=0) > 0)
    y0, y1 = te[rows[0]], te[rows[-1] + 1]
    x0, x1 = rex[cols[0]], rex[cols[-1] + 1]
    a.set_xlim(x0, x1); a.set_ylim(y0, y1)
    a.xaxis.set_major_locator(MultipleLocator(5)); a.yaxis.set_major_locator(MultipleLocator(5))
    for edge in p['pt_edges']:
        a.axvline(edge, color='white', ls='--', lw=1.2); a.axhline(edge, color='white', ls='--', lw=1.2)
    lo, hi = max(x0, y0), min(x1, y1)
    a.plot([lo, hi], [lo, hi], color='#EEEEEE', ls=':', lw=.9)
    strip.stairs(eff[rows], te[rows[0]:rows[-1] + 2], orientation='horizontal', color=BLACK, lw=1.4, baseline=None)
    strip.set_ylim(y0, y1); strip.set_xlim(0, 1.0); strip.set_xticks([0, .5, 1]); strip.set_xticklabels(['', '0.5', '1'])
    strip.tick_params(labelleft=False, labelsize=9); strip.set_xlabel(r'$\varepsilon$')
    for edge in p['pt_edges']:
        strip.axhline(edge, color=BLACK, ls='--', lw=.7, alpha=.5)
    header(f, p, .115)
    a.set_xlabel(r'$p_T^{\gamma,\mathrm{reco}}$ [GeV]'); a.set_ylabel(r'$p_T^{\gamma,\mathrm{truth}}$ [GeV]')
    cb = f.colorbar(image, cax=cax); cb.set_label(colour_label, fontsize=10.5)
    states_text(a, norm)
    save(f, 'Unfolding/%s_photon_response.png' % rstem)
    # the analysis-bin twin: same norm, the probability printed in every cell, the efficiency beside every row
    f = plt.figure(figsize=(7.4, 6.6))
    a = f.add_axes([.15, .125, .56, .755]); cax = f.add_axes([.86, .125, .024, .755])
    n = window.shape[0]
    grid = np.arange(n + 1) - .5
    image = mesh(a, grid, grid, window, norm)
    ed = groups(p)
    ticks = [f'{lo:g}–{hi:g}' for lo, hi in zip(ed[:-1], ed[1:])]
    a.set_xticks(range(n), ticks, fontsize=10); a.set_yticks(range(n), ticks, fontsize=10)
    a.tick_params(which='both', top=False, right=False)
    a.set_xlim(-.5, n - .5); a.set_ylim(-.5, n - .5)
    for i, j in np.ndindex(window.shape):
        if window[i, j] > 0:
            a.text(j, i, f'{window[i, j]:.3f}', ha='center', va='center', fontsize=11,
                   color='white' if window[i, j] < norm.vmax * .3 else BLACK)
    right = a.secondary_yaxis('right')
    right.set_yticks(range(n), [rf'$\varepsilon$ = {e:.2f}' for e in window_eff], fontsize=10)
    right.tick_params(length=0)
    header(f, p, .15)
    a.set_xlabel(r'$p_T^{\gamma,\mathrm{reco}}$ [GeV]'); a.set_ylabel(r'$p_T^{\gamma,\mathrm{truth}}$ [GeV]')
    cb = f.colorbar(image, cax=cax); cb.set_label(colour_label, fontsize=10.5)
    save(f, 'Unfolding/%s_photon_response_analysis_bins.png' % rstem)


def joint_response(p):
    """The joint (pT, xJ) response in global-bin index on its native grid (the IAN slot): P(reco | truth) with the
    unfolder's denominator, the photon-pT blocks and the analysis window marked; and the measurement window's matrix as
    its twin on the same norm."""
    if 'native_pair_response' not in p:
        return
    rstem = TAGS[str(p['label'])][1]
    nx = len(p['native_xj_edges']) - 1
    prob, _ = conditional(p['native_pair_response'], p['native_pair_misses'])
    window, _ = conditional(p['xj_response'], p['xj_misses'])
    norm = shared_norm(prob, window)
    NOTES['%s_joint_response' % rstem] = dict(source=str(p['native_source']), floor=norm.vmin, vmax=norm.vmax)
    rows = np.flatnonzero(prob.sum(axis=1) > 0); cols = np.flatnonzero(prob.sum(axis=0) > 0)
    for name, matrix, truth_pt, reco_pt, n_xj in (
            ('joint_response', prob, p['native_truth_pt_edges'], p['native_reco_pt_edges'], nx),
            ('joint_response_analysis_bins', window, p['pt_edges'], p['pt_edges'], len(p['xj_edges']) - 1)):
        f, a = axes(matrix=True)
        grid_y, grid_x = np.arange(matrix.shape[0] + 1) - .5, np.arange(matrix.shape[1] + 1) - .5
        image = mesh(a, grid_x, grid_y, matrix, norm)
        for k in range(1, len(reco_pt) - 1):
            a.axvline(k * n_xj - .5, color=BLACK, lw=.45, alpha=.5)
        for k in range(1, len(truth_pt) - 1):
            a.axhline(k * n_xj - .5, color=BLACK, lw=.45, alpha=.5)
        if name == 'joint_response':
            for edge in (p['pt_edges'][0], p['pt_edges'][-1]):
                if np.any(np.isclose(reco_pt, edge)):
                    a.axvline(int(np.flatnonzero(np.isclose(reco_pt, edge))[0]) * n_xj - .5, color='white', ls='--', lw=1.2)
                if np.any(np.isclose(truth_pt, edge)):
                    a.axhline(int(np.flatnonzero(np.isclose(truth_pt, edge))[0]) * n_xj - .5, color='white', ls='--', lw=1.2)
            a.set_xlim(cols[0] - .5, cols[-1] + .5); a.set_ylim(rows[0] - .5, rows[-1] + .5)
        a.set_xlabel(r'Reconstructed $(p_T^\gamma,x_{J\gamma})$ global bin'); a.set_ylabel(r'Truth $(p_T^\gamma,x_{J\gamma})$ global bin')
        label(a, p, sim=True, extra=[jet_text(p)], box=True)
        cb = f.colorbar(image, ax=a, pad=.025, fraction=.045)
        cb.set_label(r'$P(\mathrm{reco}\mid\mathrm{truth})$, rows sum to $\varepsilon$', fontsize=11)
        states_text(a, norm)
        save(f, 'Unfolding/%s_%s.png' % (rstem, name))


def xj_block_responses(p):
    """Per analysis photon-pT interval, the xJ response of pairs reconstructed in that interval on the native xJ grid:
    P(reco xJ, reco pT in the interval | truth xJ, truth pT in the interval) with the unfolder's denominator (rows sum
    to the fraction reconstructed in the interval), the analysis xJ edges, the row-mean profile and y = x. Never summed
    over photon-pT intervals: the jet threshold puts the kinematic wall at xJ = pT_min / ET, which moves between them."""
    if 'native_pair_response' not in p:
        return
    rstem = TAGS[str(p['label'])][1]
    xe = p['native_xj_edges']; nx = len(xe) - 1; centers_x = (xe[:-1] + xe[1:]) / 2
    te, rex = p['native_truth_pt_edges'], p['native_reco_pt_edges']
    matrix, misses = p['native_pair_response'], p['native_pair_misses']
    neff = p.get('native_pair_neff')
    blocks = []
    for lo, hi in zip(p['pt_edges'][:-1], p['pt_edges'][1:]):
        t = np.flatnonzero((te[:-1] >= lo - 1e-9) & (te[1:] <= hi + 1e-9)); r = np.flatnonzero((rex[:-1] >= lo - 1e-9) & (rex[1:] <= hi + 1e-9))
        t_rows = (t[:, None] * nx + np.arange(nx)[None, :]); r_cols = (r[:, None] * nx + np.arange(nx)[None, :])
        denominator = (matrix[t_rows.ravel()].sum(axis=1) + misses[t_rows.ravel()]).reshape(len(t), nx).sum(axis=0)
        block = matrix[np.ix_(t_rows.ravel(), r_cols.ravel())].reshape(len(t), nx, len(r), nx).sum(axis=(0, 2))
        raw = block.sum(axis=1)
        prob = np.divide(block, denominator[:, None], out=np.zeros_like(block), where=denominator[:, None] > 0)
        if neff is not None:
            sparse = neff[np.ix_(t_rows.ravel(), r_cols.ravel())].reshape(len(t), nx, len(r), nx).sum(axis=(0, 2)).sum(axis=1) < 30
        else:
            sparse = raw < raw.max() * 1e-3
        blocks.append((lo, hi, prob, sparse))
    norm = shared_norm(*[b[2] for b in blocks])
    for lo, hi, prob, sparse in blocks:
        f, a = axes(matrix=True)
        image = mesh(a, xe, xe, prob, norm)
        for edge in p['xj_edges']:
            if xe[0] < edge < xe[-1]:
                a.axvline(edge, color='white', ls='--', lw=.9, alpha=.9); a.axhline(edge, color='white', ls='--', lw=.9, alpha=.9)
        keep = (~sparse) & (prob.sum(axis=1) > 0)
        mean = np.divide(prob @ centers_x, prob.sum(axis=1), out=np.full(nx, np.nan), where=prob.sum(axis=1) > 0)
        a.plot(mean[keep], centers_x[keep], 'o', color='white', mec=BLACK, ms=4.2, mew=.8, label='Row mean', zorder=5)
        top = min(2.0, xe[-1])
        a.plot([0, top], [0, top], color='#EEEEEE', ls=':', lw=1.0, label='$y=x$')
        a.set_xlim(xe[0], top); a.set_ylim(xe[0], top)
        a.set_xlabel(r'$x_{J\gamma}^{\mathrm{reco}}$'); a.set_ylabel(r'$x_{J\gamma}^{\mathrm{truth}}$')
        label(a, p, [lo, hi], sim=True, extra=[jet_text(p)], box=True)
        a.legend(loc='lower right', fontsize=9, frameon=True, facecolor='white', edgecolor='none')
        cb = f.colorbar(image, ax=a, pad=.025, fraction=.045)
        cb.set_label(r'$P(x_{J\gamma}^{\mathrm{reco}}\mid x_{J\gamma}^{\mathrm{truth}})$, same $p_T^\gamma$ interval', fontsize=10.5)
        states_text(a, norm)
        save(f, 'Unfolding/%s_xj_response_pt%gto%g.png' % (rstem, lo, hi))


def iteration_stability(p):
    """Relative change of the per-photon result between successive iterations of the diagnostic scan, the result's
    iterations marked (the statistics of the frozen rule, when it ran, are in the result's iteration_rule)."""
    stem = TAGS[str(p['label'])][0]
    m = mask(p)
    ks = sorted(int(k.rsplit('_', 1)[1]) for k in p if k.startswith('per_photon_k_'))
    ys = {k: p['per_photon_k_%d' % k] for k in ks}
    change = [np.sqrt(np.mean(((ys[b] - ys[a]) / np.where(np.abs(ys[b]) > 0, ys[b], 1))[m] ** 2)) for a, b in zip(ks[:-1], ks[1:])]
    rel_stat = float(np.sqrt(np.mean((p['sigma_total'] / np.where(np.abs(p['per_photon']) > 0, p['per_photon'], 1))[m] ** 2)))
    f, a = axes()
    a.plot(ks[1:], change, marker='s', color=ORANGE, ms=5.2, lw=1.5, label=r'$D_k$ (change from the previous $k$)')
    a.axhline(rel_stat, color=BLUE, ls=':', lw=1.3, label=r'$S_{\mathrm{stat}}$ at the result')
    chosen = int(p['pair_iterations'])
    a.axvline(chosen, ymin=0, ymax=.67, ls='--', color=BLACK, lw=1)
    label(a, p, groups(p))
    a.legend(loc='upper right', handlelength=2)
    a.text(.965, .60, rf'$k_{{\mathrm{{Bayes}}}}={chosen}$', transform=a.transAxes, ha='right', fontsize=12)
    a.set_xticks(ks); a.set_xlabel('Bayesian iterations, $k$'); a.set_ylabel('Relative change (reported bins, rms)')
    a.set_yscale('log')
    a.set_ylim(min(change) * .5, max(rel_stat, max(change)) * 4)
    save(f, 'Unfolding/%s_joint_iteration_stability.png' % stem)


def photon_iteration_stability(p):
    """The unfolded photon yield per photon-pT bin: relative change between successive iterations of the diagnostic
    scan and its statistical precision at the result's iterations (data toys), the result's k marked."""
    stem = TAGS[str(p['label'])][0]
    ks = sorted(int(k.rsplit('_', 1)[1]) for k in p if k.startswith('photons_unfolded_k_'))
    if len(ks) < 2:
        return
    ys = {k: p['photons_unfolded_k_%d' % k] for k in ks}
    change = [float(np.sqrt(np.mean(((ys[b] - ys[a]) / ys[b]) ** 2))) for a, b in zip(ks[:-1], ks[1:])]
    rel_stat = float(np.sqrt(np.mean(np.diag(p['photons_unfolded_cov']) / p['photons_unfolded'] ** 2)))
    chosen = int(p['photon_iterations'])
    f, a = axes()
    a.plot(ks[1:], change, marker='s', color=ORANGE, ms=5.2, lw=1.5, mfc='white', label=r'$D_k$ (change from the previous $k$)')
    a.axhline(rel_stat, color=BLUE, ls=':', lw=1.3, label=r'$S_{\mathrm{stat}}$ at the result')
    a.axvline(chosen, ymin=0, ymax=.67, ls='--', color=BLACK, lw=1)
    label(a, p, groups(p))
    a.legend(loc='upper right', handlelength=2)
    a.text(.965, .60, rf'$k_{{\mathrm{{Bayes}}}}={chosen}$', transform=a.transAxes, ha='right', fontsize=12)
    a.set_xticks(ks); a.set_xlabel('Bayesian iterations, $k$'); a.set_ylabel('Relative change of the photon yield (rms)')
    a.set_yscale('log')
    a.set_ylim(max(min(change) * .5, 1e-6), max(rel_stat, max(change)) * 4)
    save(f, 'Unfolding/%s_photon_yield_iteration_stability%s.png' % (stem, '_publication' if stem == 'pp' else ''))


def photon_closure(p, mode):
    """The photon-yield closure of the unfolding alone: full (an identity by construction) or a half."""
    if 'photon_closure_%s_truth' % mode not in p:
        return
    stem = TAGS[str(p['label'])][0]
    f, a, b = axes(ratio=True)
    ed = groups(p); x = (ed[:-1] + ed[1:]) / 2; ex = np.diff(ed) / 2
    truth, u, err = p['photon_closure_%s_truth' % mode], p['photon_closure_%s_unfolded' % mode], p['photon_closure_%s_sigma' % mode]
    points(a, x, truth, None, BLUE, 'o', 'Truth' if mode == 'full' else 'Truth of the other half', ex, True)
    points(a, x, u, err, RED, 'o', 'Unfolded' if mode == 'full' else 'Unfolded half', ex)
    a.set_yscale('log'); a.set_ylim(min(truth.min(), u.min()) * .6, max(truth.max(), u.max()) * 7)
    a.set_ylabel('Weighted photons')
    points(b, x, u / truth, err / truth, RED, 'o', xerr=ex)
    b.set_xlim(ed[0], ed[-1]); b.set_xticks(ed)
    b.set_xlabel(r'$p_T^{\gamma,\mathrm{truth}}$ [GeV]')
    text = 'Same-population closure' if mode == 'full' else 'Half closure: ' + mode.replace('_to_', ' unfolds, truth of ')
    label(a, p, sim=True, extra=[text, rf'$k_{{\mathrm{{Bayes}}}}={int(p["photon_iterations"])}$'])
    a.legend(loc='upper right'); b.axhline(1, color=BLACK, ls='--', lw=.8)
    b.set_ylim(.7, 1.3); b.set_ylabel('Unfolded / truth', fontsize=11)
    save(f, ('Unfolding/%s_photon_closure_overlay.png' % stem) if mode == 'full'
         else ('Closure/%s_photon_half_closure_%s.png' % (stem, mode)))


def photon_yield(p):
    stem = TAGS[str(p['label'])][0]
    f, a, b = axes(ratio=True)
    ed = groups(p); x = (ed[:-1] + ed[1:]) / 2; ex = np.diff(ed) / 2
    r, u = p['photons_corrected'], p['photons_unfolded']
    vr, vu = np.diag(p['photons_corrected_cov']), np.diag(p['photons_unfolded_cov'])
    points(a, x, r, np.sqrt(vr), RED, 'o', 'Corrected reconstructed', ex)
    points(a, x, u, np.sqrt(vu), BLACK, 'o', 'Unfolded', ex, True)
    label(a, p, extra=[rf'$k_{{\mathrm{{Bayes}}}}={int(p["photon_iterations"])}$'])
    a.legend(loc='upper right'); a.set_ylabel('Photon yield')
    a.set_yscale('log'); a.set_ylim(min(r.min(), u.min()) * .55, max(r.max(), u.max()) * 6)
    cross = np.diag(p['photons_unfolded_corrected_cross'])
    ratio_var = vu / r ** 2 + u ** 2 * vr / r ** 4 - 2 * u * cross / r ** 3
    points(b, x, u / r, np.sqrt(np.maximum(ratio_var, 0)), BLUE, 'o', xerr=ex, open_=True)
    b.axhline(1, color=BLACK, lw=.8, ls='--')
    b.set_ylabel('Unfolded / reco', fontsize=11); b.set_xlabel(r'$p_T^\gamma$ [GeV]')
    b.set_xlim(ed[0], ed[-1]); b.set_xticks(ed)
    b.set_ylim(.75, max(2.5, (u / r).max() * 1.22))
    save(f, 'Unfolding/%s_photon_yield_unfolding_nominal_%dbin_k%d.png' % (stem, len(r), int(p['photon_iterations'])))


def closure(p, mode):
    """The estimator's simulation closure: full (same population) or a half (one half unfolds, the other is truth)."""
    if 'closure_%s_truth' % mode not in p:
        return
    stem = TAGS[str(p['label'])][0]
    f, a, b = axes(ratio=True)
    m = mask(p); x, w = centers(p)
    truth, u, sig = p['closure_%s_truth' % mode][-1], p['closure_%s_unfolded' % mode][-1], p['closure_%s_sigma_of_difference' % mode][-1]
    plot_spec(a, p, truth, np.zeros_like(truth), BLUE, 'o', 'Truth', True)
    plot_spec(a, p, u, sig, RED, 'o', 'Unfolded')
    room(a, [truth[m], u[m]], [np.zeros(m.sum()), sig[m]], .60)
    a.set_ylabel(r'$(1/N_\gamma)\,dN_{\mathrm{pair}}/dx_{J\gamma}$')
    good = m & (np.abs(truth) > 0)
    points(b, x[good], u[good] / truth[good], sig[good] / np.abs(truth[good]), RED, 'o', xerr=w[good] / 2)
    drawx(b)
    text = 'Same-population closure' if mode == 'full' else 'Half closure: ' + mode.replace('_to_', ' unfolds, truth of ')
    label(a, p, groups(p), sim=True, extra=[text])
    a.legend(loc='upper right'); b.axhline(1, color=BLACK, ls='--', lw=.8)
    b.set_ylim(.6, 1.4); b.set_ylabel('Unfolded / truth', fontsize=11)
    save(f, ('Unfolding/%s_xj_closure_overlay.png' if mode == 'full' else 'Closure/%s_xj_half_closure_' + mode + '.png') % stem)


def refold(p):
    stem = TAGS[str(p['label'])][0]
    f, a, b = axes(ratio=True)
    m = mask(p); x, w = centers(p)
    n_pt, n_xj = p['corrected'].shape
    y, ye = summed(p, 'corrected', cov_key='corrected_cov')
    ref = p['pairs_refolded'][:n_pt * n_xj].reshape(n_pt, n_xj).sum(axis=0) / w
    plot_spec(a, p, y, ye, BLACK, 'o', 'Corrected measured')
    step_spec(a, p, ref, BLUE, 'Refolded')
    label(a, p, groups(p)); a.legend(loc='upper right')
    a.set_ylabel(r'$dN_{\mathrm{pair}}/dx_{J\gamma}$')
    room(a, [y[m], ref[m]], [ye[m], np.zeros(m.sum())], .65)
    good = m & (ye > 0)
    points(b, x[good], (ref[good] - y[good]) / ye[good], None, BLUE, 'o', xerr=w[good] / 2)
    b.axhline(0, color=BLACK, lw=.8); b.axhspan(-2, 2, color=GREEN, alpha=.07)
    b.set_ylim(-4, 4); b.set_ylabel('Residual / stat.', fontsize=11); drawx(b)
    save(f, 'Closure/%s_joint_refold.png' % stem)


def reco_per_photon(p):
    """The corrected reconstructed per-photon density integrated over photon pT, with its data-toy error."""
    _, w = centers(p)
    n_pt, n_xj = p['corrected'].shape
    y = p['corrected'].sum(axis=0) / p['photons_corrected'].sum() / w
    proj = np.kron(np.ones((1, n_pt)), np.eye(n_xj)) / p['photons_corrected'].sum()
    e = np.sqrt(np.maximum(np.diag(proj @ p['corrected_cov'] @ proj.T), 0)) / w
    return y, e


def before_after(p, pt=None):
    stem = TAGS[str(p['label'])][0]
    f, a, b = axes(ratio=True); m = mask(p); x, w = centers(p)
    n_xj = len(w)
    if pt is None:
        r, re = reco_per_photon(p)
        u, ue = p['per_photon'], p['sigma_total']
        g = groups(p)
    else:
        r = p['reco_per_photon_interval'][pt]
        re = np.sqrt(np.maximum(np.diag(p['reco_per_photon_interval_cov'])[pt * n_xj:(pt + 1) * n_xj], 0))
        u = p['per_photon_interval'][pt]
        ue = np.sqrt(np.maximum(np.diag(p['per_photon_interval_cov'])[pt * n_xj:(pt + 1) * n_xj], 0))
        g = groups(p)[pt:pt + 2]
    plot_spec(a, p, r, re, BLACK, 'o', 'Corrected reconstructed', True)
    plot_spec(a, p, u, ue, BLUE, 'D', 'Unfolded')
    label(a, p, g, extra=[jet_text(p)])
    a.legend(loc='upper right'); a.set_ylabel(r'$(1/N_\gamma)\,dN_{\mathrm{pair}}/dx_{J\gamma}$')
    room(a, [r[m], u[m]], [re[m], ue[m]], .62)
    good = m & (r > 0)
    points(b, x[good], u[good] / r[good], None, BLUE, 'o', xerr=w[good] / 2)
    b.axhline(1, color=BLACK, ls='--', lw=.8); b.set_ylim(0, max(2.5, np.max(u[good] / r[good]) * 1.18))
    b.set_ylabel('Unfolded / reco', fontsize=11); drawx(b)
    save(f, 'FinalResults/%s_xj_before_after_unfolding%s.png' % (stem, '' if pt is None else f'_pt{g[0]:g}to{g[-1]:g}'))


def money(packages, pt=None):
    f, a = axes(); ys = []; errs = []
    entries = [('pp', BLUE, 'o')] + [(k, CLASS_COLORS[k], 's') for k in ('auau_cent0_20', 'auau_cent20_50', 'auau_cent50_80') if k in packages]
    for key, color, marker in entries:
        p = packages[key]
        if pt is None:
            y, e, g = p['per_photon'], p['sigma_total'], groups(p)
        else:
            n_xj = len(p['xj_edges']) - 1
            y = p['per_photon_interval'][pt]
            e = np.sqrt(np.maximum(np.diag(p['per_photon_interval_cov'])[pt * n_xj:(pt + 1) * n_xj], 0))
            g = groups(p)[pt:pt + 2]
        plot_spec(a, p, y, e, color, marker, f'{TAGS[key][2]}, {g[0]:g}–{g[-1]:g} GeV photons')
        ys.append(y[mask(p)]); errs.append(e[mask(p)])
    thresholds = {float(packages[key]['jet_pt_min_gev']) for key, _, _ in entries}
    if len(thresholds) != 1:
        raise ValueError('the overlaid results use different jet thresholds: %s' % sorted(thresholds))
    brand(a)
    a.text(.04, .882, r'$\sqrt{s_{NN}}=200$ GeV' + '\n' + r'anti-$k_T$, $R=0.4$, $p_T^{\mathrm{jet}}\geq%g$ GeV' % thresholds.pop() + '\n'
           + r'$|\eta^{\gamma,\mathrm{jet}}|<0.7$, $\Delta\phi_{\gamma\mathrm{jet}}>7\pi/8$', transform=a.transAxes, va='top',
           fontsize=10.5, linespacing=1.32)
    a.legend(loc='upper right', bbox_to_anchor=(1, .73), fontsize=9.8)
    a.set_ylabel(r'$(1/N_\gamma)\,dN_{\mathrm{pair}}/dx_{J\gamma}$'); drawx(a)
    room(a, ys, errs, .64)
    save(f, 'FinalResults/pp_vs_auau_xj_unfolded' + ('' if pt is None else '_upper_photon_bin') + '.png')


def lhc(p, ref):
    key = 'pp' if str(p['label']) == 'pp' else ('pbpb_0_10' if str(p['label']) == 'auau_cent0_20' else None)
    if key is None:
        return
    stem = TAGS[str(p['label'])][0]
    f, a = axes(); m = mask(p)
    y, e, g = p['per_photon'], p['sigma_total'], groups(p)
    plot_spec(a, p, y, e, BLUE, 'o', f'{TAGS[str(p["label"])][2]}, 200 GeV, photons {g[0]:g}–{g[-1]:g} GeV')
    rows = ref['bins']
    ed = np.array([rows[0]['low']] + [r['high'] for r in rows]); x = (ed[:-1] + ed[1:]) / 2
    yy = np.array([r[key]['value'] for r in rows]); er = np.array([r[key]['stat'] for r in rows]); sy = np.array([r[key]['sys'] for r in rows])
    for lo, hi, val, ee in zip(ed[:-1], ed[1:], yy, sy):
        a.add_patch(Rectangle((lo, val - ee), hi - lo, 2 * ee, facecolor=ORANGE, edgecolor='none', alpha=.14))
    sel = ref['selection']
    system = 'pp' if key == 'pp' else 'Pb+Pb %g–%g%%' % tuple(sel['pbpb_centrality_percent'])
    points(a, x, yy, er, ORANGE, 's', 'ATLAS %s, %g TeV' % (system, sel['sqrt_s_nn_gev'] / 1000), np.diff(ed) / 2, True)
    brand(a)
    a.text(.04, .882, r'anti-$k_T$, $R=0.4$, $\Delta\phi>7\pi/8$' + '\n' + r'sPHENIX: $p_T^{\mathrm{jet}}\geq%g$ GeV' % float(p['jet_pt_min_gev'])
           + '\n' + r'ATLAS: $%g<p_T^\gamma<%g$ GeV, $p_T^{\mathrm{jet}}\geq%g$ GeV' % (*sel['photon_pt_gev'], sel['jet_pt_min_gev']),
           transform=a.transAxes, va='top', fontsize=9.7, linespacing=1.35)
    a.legend(loc='upper right', bbox_to_anchor=(1, .735), fontsize=9.2)
    a.text(.97, .56, 'ATLAS shading: systematic uncertainty', transform=a.transAxes, ha='right', fontsize=8.8, color=ORANGE)
    a.set_ylabel(r'$(1/N_\gamma)\,dN_{\mathrm{pair}}/dx_{J\gamma}$'); drawx(a)
    room(a, [y[m], yy], [e[m], er + sy], .54)
    save(f, 'LHCComparison/%s_xj_atlas_overlay.png' % stem)


def main(argv=None):
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument('--packages', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--atlas', type=Path, default=Path(__file__).with_name('inputs') / 'atlas_reference.json')
    args = ap.parse_args(argv)
    OUT = args.output.resolve()
    FILES.clear(); NOTES.clear()      # the manifest lists this call's figures only
    packages = {path.stem: dict(np.load(path, allow_pickle=False)) for path in sorted(args.packages.glob('*.npz')) if path.stem != 'qa_inputs'}
    if not packages or set(packages) - set(TAGS):
        raise SystemExit('figure packages %s in %s: each must be one of %s' % (sorted(packages), args.packages, sorted(TAGS)))
    ref = json.loads(args.atlas.read_text())
    for key, p in packages.items():
        correction(p); photon_response(p); joint_response(p); xj_block_responses(p); iteration_stability(p); photon_iteration_stability(p); photon_yield(p)
        for mode in ('full', 'A_to_B', 'B_to_A'):
            closure(p, mode)
            photon_closure(p, mode)
        refold(p)
        for pt in (None,) + tuple(range(len(p['pt_edges']) - 1)):
            before_after(p, pt)
        lhc(p, ref)
    if 'pp' in packages and any(k.startswith('auau') for k in packages):
        money(packages); money(packages, len(packages['pp']['pt_edges']) - 2)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'FIGURE_MANIFEST.json').write_text(json.dumps({
        'renderer': str(Path(__file__).resolve()), 'renderer_sha256': sha(__file__),
        'packages': {k: sha(args.packages / (k + '.npz')) for k in packages}, 'atlas_reference_sha256': sha(args.atlas),
        'figures': FILES, 'response_notes': NOTES, 'no_curve_smoothing': True, 'no_shape_rescaling': True}, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
