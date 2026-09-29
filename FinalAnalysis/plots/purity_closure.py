"""Purity of the photon sample: the background correlation R_bkg and the closure of the ABCD method, from the QA fills
of one production (TreeToHists FillPhotonQA.C: event-leading region photons abcd_<label>_<c> and their recoil jets
recoil_<label>_<c>, labelled by the analysis truth signal), with every simulation sample weighted as the measurement
weights it.

Method (ATLAS Phys. Rev. C 93 (2016) 034914 Eq. 4; PPG12 IAN Sec. 4.3; PPG18 closure at the xJ level):
  pseudo-data     every simulated candidate, photon and inclusive families (data-like mixture)
  truth signal    Region-A leaders matched to an analysis-signal truth photon (source role 1), the signal the leakage
                  is normalised to and the measurement unfolds
  leakage         signal leaders per region / signal leaders in A
  R_bkg           bA bD / (bB bC) of the non-signal leaders
  solution        the measurement's own code: FinalAnalysis solve_leakage_abcd(counts, leakage, R) for the photons and
                  correct_measured (per_pt, R, no combinatoric subtraction) for the purity-corrected recoil pairs
Closures, per class and photon-pT bin: (0) R = 1 (the method of ATLAS and PPG12, nominal until decided); (1) R and leakage
from the same simulation (an identity: must close to rounding); (2) R and leakage from one half of the output groups
applied to the pseudo-data of the other half, both ways (the real test); photon-yield purity and the per-photon xJ
distribution (A - k C)/(1 - k l_C) / s against the truth-signal photons' own per-photon xJ, as PPG18. R's statistical
uncertainty is a bootstrap over output groups. With DATA QA, the factorisation is also tested where every region is
background: non-tight score sub-bands x non-isolated sub-bands in DATA and in the simulated background.

    python -m FinalAnalysis.plots.purity_closure --hists ROOT/PhotonJetHists_v1 --config FinalAnalysis/config/final_analysis.yaml --out DIR
        [--sample-weights W.json] [--unit-weights] [--bootstrap 200]
Writes PURITY_CLOSURE.json (and background_correlation per systems entry, the FinalAnalysis config values) and figures.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import uproot
import yaml

from FinalAnalysis.src.internal import corrections as CORRECTIONS

from . import render as R
from . import render_qa as Q

HERE = Path(__file__).resolve().parent

plt = R.plt
SIGNAL = {'phosim': 'sig', 'jetsim': 'jetprompt'}
BACKGROUND = {'phosim': 'phoother', 'jetsim': 'bkg'}
FINE = Q.FINE_EDGES


def owned(values, edges, window, axis):
    """Sum the ownership classes inside [lo, hi) (hi None: open) along `axis` of an ownership-class histogram."""
    lo, hi = window
    hi = np.inf if hi is None else hi
    keep = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= hi + 1e-9) if np.isfinite(hi) else (edges[:-1] >= lo - 1e-9)
    return np.compress(keep, values, axis=axis).sum(axis=axis)


def groups(store, system, ownership):
    """Per output group of every simulated sample: weighted {label: abcd (4, 20)} and {label: recoil (2, 20, n_xj)} per
    class, with the sample's family; the unit of the bootstrap and of the half split. A sample with an ownership window
    keeps only events whose generator variable lies in it (the abcdown/recoilown fills)."""
    out = []
    for (product, sample), family in sorted(store.samples(system).items()):
        w = store.weight(sample)
        window = ownership.get(sample)
        for path in store.files('qa_sim', product=product, sample=sample):
            with uproot.open(path) as f:
                keys = {k.split(';')[0] for k in f.keys()}
                row = dict(sample=sample, family=family, path=str(path), window=window, abcd={}, recoil={})
                for c, _ in Q.CLASSES[system]:
                    for label in (SIGNAL[family], BACKGROUND[family]):
                        a, r = 'abcd_%s_%s' % (label, c), 'recoil_%s_%s' % (label, c)
                        if window is not None:
                            ao, ro = 'abcdown_%s_%s' % (label, c), 'recoilown_%s_%s' % (label, c)
                            if ao in keys:
                                h = f[ao]
                                classes = np.asarray(h.axes[2].edges(), dtype=float)   # the ownership classes of the fill
                                row['abcd'][(label, c)] = owned(np.asarray(h.values(), dtype=float), classes, window, 2) * w
                                row['et_edges'] = np.asarray(h.axes[1].edges(), dtype=float)
                            if ro in keys:
                                # (region A, C) x ownership class on one axis; recoilown is filled only with abcdown
                                h = f[ro]
                                v = np.asarray(h.values(), dtype=float)
                                n_cls = v.shape[0] // 2
                                if ao not in keys or len(classes) - 1 != n_cls:
                                    raise ValueError('%s: %s does not match the ownership classes of %s' % (path, ro, ao))
                                row['recoil'][(label, c)] = np.stack([owned(v[:n_cls], classes, window, 0),
                                                                      owned(v[n_cls:], classes, window, 0)]) * w
                                row['xj_edges'] = np.asarray(h.axes[2].edges(), dtype=float)
                                row['recoil_et_edges'] = np.asarray(h.axes[1].edges(), dtype=float)
                            continue
                        if a in keys:
                            row['abcd'][(label, c)] = np.asarray(f[a].values(), dtype=float) * w
                            row['et_edges'] = np.asarray(f[a].axes[1].edges(), dtype=float)
                        if r in keys:
                            row['recoil'][(label, c)] = np.asarray(f[r].values(), dtype=float) * w
                            row['xj_edges'] = np.asarray(f[r].axes[2].edges(), dtype=float)
                            row['recoil_et_edges'] = np.asarray(f[r].axes[1].edges(), dtype=float)
                out.append(row)
    return out


def ownership_windows(path):
    """{sample: (lo, hi or None)} from a samples configuration (TreeToHists/config/samples.yaml ownership_window_gev)."""
    if path is None:
        return {}
    doc = yaml.safe_load(Path(path).read_text()) or {}
    out = {}
    for name, spec in (doc.get('samples') or {}).items():
        window = (spec or {}).get('ownership_window_gev')
        if window:
            out[name] = (float(window[0]), None if len(window) < 2 or window[1] is None else float(window[1]))
    return out


def total(rows, c, kind, which):
    """Sum over rows of the signal ('signal'), background ('background') or all ('all') labels of class c."""
    acc = None
    for row in rows:
        labels = {'signal': [SIGNAL[row['family']]], 'background': [BACKGROUND[row['family']]],
                  'all': [SIGNAL[row['family']], BACKGROUND[row['family']]]}[which]
        for label in labels:
            v = row[kind].get((label, c))
            if v is not None:
                acc = v.copy() if acc is None else acc + v
    return acc


def regroup(values, edges, target):
    return Q.regrouped(values, edges, target)


def solvable(n, l, r):
    return n[0] > 0 and l[0] > 0 and np.all(np.isfinite(l)) and np.isfinite(r) and r > 0


def solve(counts, leakage, r):
    """Per bin: the leakage-corrected ABCD signal (FinalAnalysis solve_leakage_abcd with R); NaN where it cannot be solved."""
    out = np.full(counts.shape[1], np.nan)
    for k in range(counts.shape[1]):
        n, l = np.maximum(counts[:, k], 0), leakage[:, k]
        if solvable(n, l, r[k]):
            out[k] = CORRECTIONS.solve_leakage_abcd(n, l, float(r[k])).signal_a
    return out


def corrected_pairs(counts, leakage, r, spectra, edges):
    """Per photon-pT bin, the purity-corrected recoil pairs (A - k C)/(1 - k l_C) computed by the measurement itself
    (FinalAnalysis correct_measured, per_pt strategy, background correlation r, no combinatoric subtraction), so the
    closure tests the code the result uses. `spectra` is (region A/C, photon-pT bin, xJ); NaN rows where a bin has no
    physical correction."""
    n_pt, nx = counts.shape[1], spectra.shape[2]
    out = np.full((n_pt, nx), np.nan)
    empty = np.zeros((1, nx))
    for k in range(n_pt):
        n, l = np.maximum(counts[:, k], 0), leakage[:, k]
        if not solvable(n, l, r[k]):
            continue
        data = dict(abcd_counts=n[:, None], abcd_counts_sumw2=np.zeros((4, 1)), recoil_spectra=spectra[:, k:k + 1],
                    recoil_spectra_sumw2=np.zeros((2, 1, nx)))
        response = CORRECTIONS.CorrectionInputs(leakage=l[:, None], leakage_sumw2=np.zeros((4, 1)), combinatoric_reco=empty,
                                              combinatoric_reco_sumw2=empty, combinatoric_normalization_reco=np.ones(1),
                                              combinatoric_normalization_reco_sumw2=np.zeros(1))
        corrections = CORRECTIONS.PurityCorrections(background_correlation=(float(r[k]),))
        try:
            out[k] = CORRECTIONS.correct_measured(data, response, False, pt_edges=edges[k:k + 2],
                                                purity_corrections=corrections)['corrected'][0]
        except ValueError:  # the measurement refuses this bin (no physical leakage denominator)
            pass
    return out


def correlation(background):
    b = background
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where((b[0] > 0) & (b[1] > 0) & (b[2] > 0) & (b[3] > 0), b[0] * b[3] / (b[1] * b[2]), np.nan)


def closure(test_rows, model_rows, c, edges, use_r=True):
    """Photon-purity and per-photon xJ closure of pseudo-data `test_rows` with leakage and R from `model_rows`."""
    if not test_rows or not model_rows:
        return dict(status='insufficient: an empty half')
    parts = [total(test_rows, c, 'abcd', 'all'), total(test_rows, c, 'abcd', 'signal'), total(model_rows, c, 'abcd', 'signal'),
             total(model_rows, c, 'abcd', 'background')]
    if any(p is None for p in parts):
        return dict(status='insufficient: a label has no region leaders')
    with np.errstate(divide='ignore', invalid='ignore'):
        return _closure(test_rows, model_rows, c, edges, use_r)


def _closure(test_rows, model_rows, c, edges, use_r):
    counts = regroup(total(test_rows, c, 'abcd', 'all'), test_rows[0]['et_edges'], edges)
    truth_a = regroup(total(test_rows, c, 'abcd', 'signal'), test_rows[0]['et_edges'], edges)[0]
    signal_model = regroup(total(model_rows, c, 'abcd', 'signal'), model_rows[0]['et_edges'], edges)
    background_model = regroup(total(model_rows, c, 'abcd', 'background'), model_rows[0]['et_edges'], edges)
    r = correlation(background_model) if use_r else np.ones(len(edges) - 1)
    signal = solve(counts, signal_model, r)
    purity_abcd = signal / counts[0]
    purity_truth = truth_a / counts[0]
    out = dict(r=r.tolist(), purity_abcd=purity_abcd.tolist(), purity_truth=purity_truth.tolist(),
               ratio=(purity_abcd / purity_truth).tolist())
    rec_all = total(test_rows, c, 'recoil', 'all')
    rec_sig = total(test_rows, c, 'recoil', 'signal')
    if rec_all is not None and rec_sig is not None:
        re = test_rows[0]['recoil_et_edges']
        xe = test_rows[0]['xj_edges']
        # (region, ET, xJ) -> (region, measurement bin, xJ)
        a = np.stack([regroup(rec_all[0].T, re, edges).T, regroup(rec_all[1].T, re, edges).T])
        sig = regroup(rec_sig[0].T, re, edges).T
        width = np.diff(xe)
        pairs = corrected_pairs(counts, signal_model, r, a, edges)
        per_photon, truth = [], []
        for k in range(len(edges) - 1):
            s = signal[k]
            if not np.isfinite(s) or s <= 0:
                per_photon.append(np.full(len(width), np.nan)); truth.append(np.full(len(width), np.nan))
                continue
            per_photon.append(pairs[k] / s / width)
            truth.append(sig[k] / truth_a[k] / width if truth_a[k] > 0 else np.full(len(width), np.nan))
        out.update(xj_edges=xe.tolist(), per_photon=np.array(per_photon).tolist(), per_photon_truth=np.array(truth).tolist())
    return out


def transfer_terms(rows, c, edges):
    """Why the xJ-level ABCD subtraction (A - k C)/(1 - k l_C) does not close in its own simulation (the identity closure,
    R and leakage from the same rows): per photon-pT bin and xJ bin, per signal photon and unit xJ,
        background   (B_A / S_A) (g_A - g_C) / (1 - k l_C)   background leaders recoil differently in A and in C
        leakage      k l_C (f_A - f_C) / (1 - k l_C)         signal leaking into C recoils unlike the signal in A
    with g the recoil per background leader, f per signal leader, S and B the signal and background leaders, k = B_A / B_C
    and l_C = S_C / S_A. Their sum is the identity closure's per-photon difference to the truth. The background term is
    the input of the ABCD background-transfer systematic."""
    rec_s, rec_b = total(rows, c, 'recoil', 'signal'), total(rows, c, 'recoil', 'background')
    lead_s, lead_b = total(rows, c, 'abcd', 'signal'), total(rows, c, 'abcd', 'background')
    if any(x is None for x in (rec_s, rec_b, lead_s, lead_b)):
        return dict(status='insufficient: a label has no leaders or recoil')
    re, xe = rows[0]['recoil_et_edges'], rows[0]['xj_edges']
    width = np.diff(xe)
    on_bins = lambda v: np.stack([regroup(v[0].T, re, edges).T, regroup(v[1].T, re, edges).T])  # (A/C, pT bin, xJ)
    rs, rb = on_bins(rec_s), on_bins(rec_b)
    s_lead, b_lead = regroup(lead_s, rows[0]['et_edges'], edges), regroup(lead_b, rows[0]['et_edges'], edges)
    with np.errstate(divide='ignore', invalid='ignore'):
        f_a, f_c = rs[0] / s_lead[0][:, None], rs[1] / s_lead[2][:, None]
        g_a, g_c = rb[0] / b_lead[0][:, None], rb[1] / b_lead[2][:, None]
        k, l_c = b_lead[0] / b_lead[2], s_lead[2] / s_lead[0]
        denominator = (1.0 - k * l_c)[:, None]
        background = (b_lead[0] / s_lead[0])[:, None] * (g_a - g_c) / denominator / width
        leakage = (k * l_c)[:, None] * (f_a - f_c) / denominator / width
    return dict(xj_edges=np.asarray(xe).tolist(), background=background.tolist(), leakage=leakage.tolist(),
                recoil_per_background_leader_a=(g_a / width).tolist(), recoil_per_background_leader_c=(g_c / width).tolist(),
                recoil_per_signal_leader_a=(f_a / width).tolist(), recoil_per_signal_leader_c=(f_c / width).tolist(),
                background_over_signal_a=(b_lead[0] / s_lead[0]).tolist(), k=k.tolist(), l_c=l_c.tolist())


def bootstrap_r(rows, c, edges, n, seed):
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n):
        pick = [rows[i] for i in rng.integers(0, len(rows), size=len(rows))]
        background = total(pick, c, 'abcd', 'background')
        if background is not None:
            draws.append(correlation(regroup(background, pick[0]['et_edges'], edges)))
    if len(draws) < 2:
        return np.full(len(edges) - 1, np.nan)
    with np.errstate(invalid='ignore'):
        return np.nanstd(np.array(draws), axis=0, ddof=1)


def data_factorisation(store, system, c):
    """R in background-only sub-bands, DATA against the simulated background: score in [0.1, 0.4) and [0.4, 0.6)
    (non-tight everywhere), isolation in [6, 12) and [12, 40) GeV (non-isolated everywhere), 15-35 GeV."""
    def cells(h):
        v = h.values
        score, iso = h.edges[0], h.edges[1]
        s1 = (score[:-1] >= 0.1 - 1e-9) & (score[1:] <= 0.4 + 1e-9)
        s2 = (score[:-1] >= 0.4 - 1e-9) & (score[1:] <= 0.6 + 1e-9)
        i1 = (iso[:-1] >= 6 - 1e-9) & (iso[1:] <= 12 + 1e-9)
        i2 = (iso[:-1] >= 12 - 1e-9) & (iso[1:] <= 40 + 1e-9)
        n = lambda s, i: v[np.ix_(s, i)].sum()  # noqa: E731
        return n(s1, i1), n(s1, i2), n(s2, i1), n(s2, i2)
    data = store.data(system, 'scoreiso_data_%s' % c)
    sim = None
    for fam in ('phosim', 'jetsim'):
        h = store.sim(system, fam, 'scoreiso_%s_%s' % (BACKGROUND[fam], c))
        sim = h + sim if h is not None else sim
    out = {}
    for name, h in (('data', data), ('simulated_background', sim)):
        if h is None:
            continue
        a, b, cc, d = cells(h)
        out[name] = dict(cells=[a, b, cc, d], r=(a * d / (b * cc)) if b > 0 and cc > 0 else None,
                         r_stat=(np.sqrt(sum(1.0 / x for x in (a, b, cc, d) if x > 0)) * (a * d / (b * cc))) if min(a, b, cc, d) > 0 else None)
    return out


def figure(results, key, text, out):
    res = results[key]
    fine = res['fine']
    if any('status' in fine[v] for v in ('r1', 'half_ab')):
        return
    x = (FINE[:-1] + FINE[1:]) / 2
    f, (a, b) = plt.subplots(2, 1, figsize=(7.1, 6.4), sharex=True, gridspec_kw={'height_ratios': [2.2, 1], 'hspace': .04})
    f.subplots_adjust(left=.14, right=.97, bottom=.1, top=.97)
    R.style(a); R.style(b)
    R.points(a, x, np.array(fine['r1']['purity_truth'], float), None, R.RED, 'o', 'Truth purity (analysis signal)')
    R.points(a, x, np.array(fine['r1']['purity_abcd'], float), None, R.BLACK, 'D', r'ABCD, $R_{\mathrm{bkg}}=1$')
    R.points(a, x + .25, np.array(fine['half_ab']['purity_abcd'], float), None, R.BLUE, 's', r'ABCD, $R_{\mathrm{bkg}}$ and leakage from the other half', open_=True)
    a.set_ylim(0, 1.25); a.set_ylabel('Photon purity (simulation)'); a.legend(loc='lower right', fontsize=9)
    R.brand(a); a.text(.04, .88, text, transform=a.transAxes, va='top', fontsize=10.5)
    R.points(b, x, np.array(fine['r1']['ratio'], float), None, R.BLACK, 'D')
    R.points(b, x + .25, np.array(fine['half_ab']['ratio'], float), None, R.BLUE, 's', open_=True)
    b.axhline(1, ls='--', color=R.BLACK, lw=.8); b.set_ylim(.6, 1.6); b.set_ylabel('ABCD / truth', fontsize=11)
    b.set_xlabel(r'$E_T^{\gamma,\mathrm{reco}}$ [GeV]'); b.set_xlim(15, 35)
    R.OUT = out
    R.save(f, 'purity/%s_purity_closure_rbkg.png' % key)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--hists', type=Path, required=True)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--sample-weights', type=Path)
    ap.add_argument('--unit-weights', action='store_true')
    ap.add_argument('--bootstrap', type=int, default=200)
    ap.add_argument('--samples', type=Path, default=HERE.parent / 'config' / 'samples.yaml',
                    help='ownership windows of the simulated samples (ownership_window_gev)')
    args = ap.parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    store = Q.Store(args.hists, config, json.loads(args.sample_weights.read_text()) if args.sample_weights else None, args.unit_weights)
    measurement = np.asarray(config['measurement']['photon_pt_edges_gev'], dtype=float)
    ownership = ownership_windows(args.samples)
    args.out.mkdir(parents=True, exist_ok=True)
    R.OUT = args.out.resolve()
    doc = dict(schema='PurityClosureV1', hists=str(args.hists), unit_weights=args.unit_weights, measurement_pt_edges=measurement.tolist(),
               method='ATLAS PRC 93 034914 Eq. 4 with R_bkg from the simulated non-signal leaders; closure halves by output group',
               ownership_windows={k: list(v) for k, v in ownership.items()}, entries={}, background_correlation={})
    for key, block in config['systems'].items():
        system, c = block.get('system', key), block['class']
        rows = groups(store, system, ownership)
        rows = [r for r in rows if any(k[1] == c for k in r['abcd'])]
        if len(rows) < 2:
            doc['entries'][key] = dict(status='not enough simulation output groups (%d)' % len(rows))
            continue
        # halves by output group within each sample, so both halves carry every sample
        seen = {}
        half_a, half_b = [], []
        for row in rows:
            k = seen.get(row['sample'], 0)
            (half_a if k % 2 == 0 else half_b).append(row)
            seen[row['sample']] = k + 1
        entry = dict(groups=len(rows))
        for name, edges in (('measurement', measurement), ('fine', FINE)):
            entry[name] = dict(r1=closure(rows, rows, c, edges, use_r=False), same=closure(rows, rows, c, edges),
                               half_ab=closure(half_b, half_a, c, edges), half_ba=closure(half_a, half_b, c, edges))
        entry['transfer'] = transfer_terms(rows, c, measurement)
        entry['r_bkg'] = entry['measurement']['same']['r']
        entry['r_bkg_bootstrap_sd'] = bootstrap_r(rows, c, measurement, args.bootstrap, 3520_0300).tolist()
        entry['data_factorisation'] = data_factorisation(store, system, c)
        doc['entries'][key] = entry
        doc['background_correlation'][key] = entry['r_bkg']
        text = ('pp, ' if system == 'pp' else dict(Q.CLASSES['auau']).get(c, c) + ', ') + Q.energy(system)
        try:
            figure(doc['entries'], key, text, args.out.resolve())
        except Exception as exc:  # a figure must never hide the numbers
            entry['figure_error'] = str(exc)
    (args.out / 'PURITY_CLOSURE.json').write_text(json.dumps(doc, indent=1, default=float) + '\n')
    summary = {k: dict(r_bkg=v.get('r_bkg'), r_sd=v.get('r_bkg_bootstrap_sd'),
                       ratio_r1=(v.get('measurement') or {}).get('r1', {}).get('ratio'),
                       ratio_same=(v.get('measurement') or {}).get('same', {}).get('ratio'),
                       ratio_halves=[(v.get('measurement') or {}).get(h, {}).get('ratio') for h in ('half_ab', 'half_ba')])
               for k, v in doc['entries'].items()}
    print(json.dumps(summary, indent=1, default=float))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
