"""eventTree quick start: select analysis events, read photons, jets and photon-jet pairs.

    python3 quick_start.py file.eventtree.root

Works for pp and Au+Au, DATA and simulation.
"""
import json
import sys

import awkward as ak
import numpy as np
import uproot

path = sys.argv[1]
with uproot.open(path) as f:
    info = json.loads(str(f['meta/schema']))
    tree = f['eventTree']
    events = tree.arrays(['evt_record', 'evt_nodes_ready', 'vtx_valid', 'evt_terminal_status', 'evt_run', 'vtx_z',
                          'cent_value', 'w_event', 'pho_et', 'pho_eta', 'jet02_pt', 'jet03_pt', 'jet04_pt', 'jet04_eta',
                          'pair_pho', 'pair_jet', 'pair_xjgamma', 'pair_delta_phi'], library='ak')

# The documented reader predicate: retained events whose reconstruction finished with a valid vertex.
analysis = events[(events.evt_record == 0) & (events.evt_nodes_ready == 1) & (events.vtx_valid == 1)
                  & (events.evt_terminal_status == 0)]
print('%s (%s, %s layout): %d entries, %d analysis events'
      % (path, info['system'], info['layout'], len(events), len(analysis)))

# Photons: one entry per candidate in every pho_ branch.
photons = ak.zip({'et': analysis.pho_et, 'eta': analysis.pho_eta})
leading = ak.firsts(photons[ak.argsort(photons.et, ascending=False)])
print('events with a photon candidate: %d' % int(ak.sum(ak.num(photons.et) > 0)))
if ak.any(~ak.is_none(leading)):
    print('leading photon ET median: %.2f GeV' % float(np.median(ak.to_numpy(ak.drop_none(leading.et)))))

# Pairs: pair_pho indexes the event's photons; pair_jet indexes the concatenated reco jet collections
# in meta/collections order (jet02, jet03, jet04, then jetunsub* in Au+Au). Keep the pairs with an R = 0.4 jet:
local04 = analysis.pair_jet - ak.num(analysis.jet02_pt) - ak.num(analysis.jet03_pt)
is04 = (local04 >= 0) & (local04 < ak.num(analysis.jet04_pt))
pt04 = analysis.jet04_pt[local04[is04]]
selected = (pt04 > 10) & (analysis.pair_delta_phi[is04] > 7 * np.pi / 8)
xj = ak.flatten(analysis.pair_xjgamma[is04][selected])
xj = ak.to_numpy(xj[~np.isnan(xj)])
print('back-to-back photon + R=0.4 jet (pT > 10 GeV) pairs: %d' % len(xj))
if len(xj):
    counts, edges = np.histogram(xj, bins=10, range=(0.0, 2.0))
    for low, high, n in zip(edges[:-1], edges[1:], counts):
        print('  xJ %.1f-%.1f  %d' % (low, high, n))
