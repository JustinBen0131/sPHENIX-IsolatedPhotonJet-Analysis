"""The IAN result tables from the result documents make ian writes (results/<entry>.json, FinalAnalysis
measurement.compute): LaTeX (booktabs) and CSV, nothing recomputed.

    tab_xj_<entry>.tex        the reported per-photon xJ distribution of one entry: value, data (statistical) and
                              response half-widths, total
    tab_xj_all.tex            every entry side by side, value +- total, in the reported xJ bins
    tab_measurement.tex       per entry and photon-pT bin: iterations (pairs, photons), purity, purity-corrected
                              photons, jets per photon (DATA / simulation)

    python -m FinalAnalysis.plots.tables --results DIR --out DIR
The two static IAN tables (the Au+Au MC overlay samples, the frozen BDT's XGBoost settings) are records, not results,
and stay as they are.
"""
import argparse
import csv
import json
from pathlib import Path

ORDER = ['pp', 'auau_cent50_80', 'auau_cent20_50', 'auau_cent0_20']
TITLE = {'pp': r'\pp', 'auau_cent0_20': r'0--20\% \auau', 'auau_cent20_50': r'20--50\% \auau', 'auau_cent50_80': r'50--80\% \auau'}


def fmt(x, digits=3):
    return ('%.' + str(digits) + 'f') % x


def reported(result):
    s = result['settings']
    edges = s['xj_edges']
    mask = result['reported_xj_mask']
    return [(i, edges[i], edges[i + 1]) for i in range(len(edges) - 1) if mask[i]]


def xj_table(key, result):
    e = result['entry']
    rows = [(lo, hi, e['central'][i], e['half_width_68_data'][i], e['half_width_68_response'][i], e['sigma_total'][i])
            for i, lo, hi in reported(result)]
    s = result['settings']
    tex = [r'\begin{table}[hbtp!]', r'  \centering', r'  \small',
           r'  \caption{Unfolded per-photon $x_{J\gamma}$ distribution in %s, $%g<p_{T}^{\gamma}<%g$~\GeV, $p_{T}^{\mathrm{jet}}\geq%g$~\GeV: '
           r'value, statistical (data) and response half-widths, and the total uncertainty of the measurement.}' % (
               TITLE.get(key, key), s['pt_edges'][0], s['pt_edges'][-1], s['jet_pt_min_gev']),
           r'  \label{tab:xj_%s}' % key, r'  \begin{tabular}{ccccc}', r'    \toprule',
           r'    $x_{J\gamma}$ interval & $(1/N_{\gamma})\,dN_{\mathrm{pair}}/dx_{J\gamma}$ & data & response & total \\', r'    \midrule']
    for lo, hi, v, d, r, t in rows:
        tex.append(r'    %g--%g & %s & %s & %s & %s \\' % (lo, hi, fmt(v), fmt(d), fmt(r), fmt(t)))
    tex += [r'    \bottomrule', r'  \end{tabular}', r'\end{table}']
    return '\n'.join(tex) + '\n', rows


def all_table(results):
    keys = [k for k in ORDER if k in results] + sorted(k for k in results if k not in ORDER)
    bins = [(lo, hi) for _, lo, hi in reported(results[keys[0]])]
    for k in keys[1:]:
        if [(lo, hi) for _, lo, hi in reported(results[k])] != bins:
            raise ValueError('the entries report different xJ bins')
    tex = [r'\begin{table}[hbtp!]', r'  \centering', r'  \small',
           r'  \caption{Unfolded per-photon $x_{J\gamma}$ distributions, value $\pm$ total uncertainty.}', r'  \label{tab:xj_all}',
           r'  \begin{tabular}{c' + 'c' * len(keys) + '}', r'    \toprule',
           r'    $x_{J\gamma}$ & ' + ' & '.join(TITLE.get(k, k) for k in keys) + r' \\', r'    \midrule']
    for b, (lo, hi) in enumerate(bins):
        cells = []
        for k in keys:
            i = reported(results[k])[b][0]
            cells.append('$%s \\pm %s$' % (fmt(results[k]['entry']['central'][i]), fmt(results[k]['entry']['sigma_total'][i])))
        tex.append(r'    %g--%g & %s \\' % (lo, hi, ' & '.join(cells)))
    tex += [r'    \bottomrule', r'  \end{tabular}', r'\end{table}']
    return '\n'.join(tex) + '\n'


def measurement_table(results):
    keys = [k for k in ORDER if k in results] + sorted(k for k in results if k not in ORDER)
    tex = [r'\begin{table}[hbtp!]', r'  \centering', r'  \small',
           r'  \caption{Measurement summary per collision system and photon-$p_{T}$ interval: Bayesian iterations of the '
           r'joint and photon unfoldings, leakage-corrected ABCD purity, purity-corrected photons, and recoil jets per '
           r'photon in DATA and simulation (reconstructed level).}',
           r'  \label{tab:measurement_summary}', r'  \begin{tabular}{llcccccc}', r'    \toprule',
           r'    System & $p_{T}^{\gamma}$ [\GeV] & $k_{\mathrm{pairs}}$ & $k_{\gamma}$ & purity & $N_{\gamma}$ & jets/$\gamma$ DATA & jets/$\gamma$ sim \\',
           r'    \midrule']
    for k in keys:
        r = results[k]
        it = r.get('iterations') or {}
        gates = (r.get('gates') or {}).get('jets_per_photon') or {}
        none = [None] * len(r['purity'])
        for j, p in enumerate(r['purity']):
            raw, sim = (gates.get('data_raw') or none)[j], (gates.get('simulation_reco') or none)[j]
            tex.append(r'    %s & %g--%g & %s & %s & %s & %s & %s & %s \\' % (
                TITLE.get(k, k) if j == 0 else '', p['pt_low'], p['pt_high'], it.get('pairs', ''), it.get('photons', ''),
                fmt(p['purity']), '%.0f' % r['purity_corrected_photons_ptgamma'][j], fmt(raw) if raw is not None else '--',
                fmt(sim) if sim is not None else '--'))
    tex += [r'    \bottomrule', r'  \end{tabular}', r'\end{table}']
    return '\n'.join(tex) + '\n'


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', type=Path, required=True, help='directory of <entry>.json (or <entry>.result.json)')
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args(argv)
    results = {}
    for p in sorted(args.results.glob('*.json')):
        doc = json.loads(p.read_text())
        if doc.get('schema') == 'PhotonJetPerPhotonResultV1':
            results[doc['system']] = doc
    if not results:
        raise SystemExit('no PhotonJetPerPhotonResultV1 documents in %s' % args.results)
    args.out.mkdir(parents=True, exist_ok=True)
    written = []
    for key, result in results.items():
        tex, rows = xj_table(key, result)
        (args.out / ('tab_xj_%s.tex' % key)).write_text(tex)
        with (args.out / ('tab_xj_%s.csv' % key)).open('w', newline='') as stream:
            w = csv.writer(stream)
            w.writerow(['xj_low', 'xj_high', 'value', 'half_width_68_data', 'half_width_68_response', 'sigma_total'])
            w.writerows(rows)
        written += ['tab_xj_%s.tex' % key, 'tab_xj_%s.csv' % key]
    if len(results) > 1:
        (args.out / 'tab_xj_all.tex').write_text(all_table(results))
        written.append('tab_xj_all.tex')
    (args.out / 'tab_measurement.tex').write_text(measurement_table(results))
    written.append('tab_measurement.tex')
    (args.out / 'TABLE_MANIFEST.json').write_text(json.dumps(dict(schema='IANTableManifestV1', results=sorted(results), tables=written), indent=1) + '\n')
    print('\n'.join(written))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
