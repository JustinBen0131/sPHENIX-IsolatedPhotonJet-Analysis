# Reading eventTree

Each ROOT file contains one `eventTree`, with one entry per event in source
order. Event quantities are scalars; photons, jets and other collections are
arrays within that entry. The same structure is used for pp/AuAu and DATA/SIM.

This guide maps the structure. [SCHEMA.md](SCHEMA.md) lists every branch, type
and unit; [config/eventTree_v1.yaml](config/eventTree_v1.yaml) defines the schema.

Jump to [event quantities](#event-quantities), [photons and jets](#photons-and-jets),
[simulation truth](#simulation-truth), or [bookkeeping](#bookkeeping).

## File layout

```text
file.root
├── eventTree          event scalars and collection arrays
├── meta/schema        schema, branch counters and filled optional groups
├── meta/source        input-source identity and event interval
├── meta/collections   jet collection order, radii and reconstruction views
├── meta/constants     file-wide definitions and fixed values
├── meta/codes         dictionaries for integer-coded strings
└── hist/              retained histograms, when present
```

Other metadata record conditions, triggers and completion checks.
See the [file-level reference](SCHEMA.md#file-level-objects) for all keys.

## Event quantities

| Branch group | Contents |
|---|---|
| [`evt_`](SCHEMA.md#event-identity-and-bookkeeping-evt_) | Run/event identity, source entry and event status |
| [`trig_`](SCHEMA.md#trigger-trig_) | Live/scaled trigger bits and decision availability |
| [`gl1_`](SCHEMA.md#gl1-scalers-gl1_) | Beam clock and 64-channel trigger scalers |
| [`vtx_`](SCHEMA.md#vertex-vtx_) | Reconstructed vertex z and validity |
| [`cent_`](SCHEMA.md#centrality-cent_) | Centrality, its validity and PMT inputs |
| [`mbd_`](SCHEMA.md#mbd-mbd_) | MBD charge, timing and detector status |
| [`calo_`](SCHEMA.md#calorimeter-sums-calo_) | Event calorimeter energy sums and availability |
| [`ue_`](SCHEMA.md#underlying-event-ue_) | Underlying-event quantities |
| [`w_`](SCHEMA.md#weights-w_) | Producer event weight and its component factors |

For analysis events, apply:

```text
evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && evt_terminal_status == 0
```

`evt_record == 1` denotes an upstream-rejected DATA event; only its `evt_`
fields are populated. The event predicate does not select photons or jets.

## Photons and jets

| Branch group | Contents and indexing |
|---|---|
| [`pho_`](SCHEMA.md#photons-pho_) | Photon kinematics, shapes, isolation, scores and validity; ordinary fields have one value per candidate |
| [`shv_`](SCHEMA.md#photon-shower-views-stride-7-shv_) | Seven shower views per photon; view-major within each photon |
| [`isow_`](SCHEMA.md#photon-isolation-witnesses-stride-2-isow_) | Two isolation records per photon, for R = 0.3 and 0.4 |
| [`cell_`](SCHEMA.md#shower-cells-parent-index-cell_pho-cell_) | Shower tower cells; `cell_pho` gives the parent photon index |
| [`isoc_`](SCHEMA.md#isolation-constituents-parent-index-isoc_pho-isoc_) | Isolation constituents; `isoc_pho` gives the parent photon index |
| [`jet02_`, `jet03_`, `jet04_`](SCHEMA.md#reco-jets-analysis-view-jetrr_) | Analysis jets at R = 0.2, 0.3 and 0.4; calibrated pp jets or UE-subtracted AuAu jets |
| [`jetunsub02_`, `jetunsub03_`, `jetunsub04_`](SCHEMA.md#reco-jets-auau-unsubtracted-view-jetunsubrr_) | AuAu unsubtracted jets at the same radii |
| [`pair_`](SCHEMA.md#photon-jet-pairs-pair_) | Photon/jet indices, Δφ, xJγ and pair classification |

Ordinary collection arrays use counters such as `n_pho`, `n_jet04` and `n_pair`.
Object links are zero-based and local to the current event; `-1` means no link.
The exact counter for each branch is recorded in `meta/schema`. Strided
members, such as model-input arrays, may have a separate counter.

The shower-view order is `H70, H0, G70, G0, O70, O0, R70`.
For example, `shv_e11[7 * photon + view]` reads that view's energy.
`shv_features` additionally contains the ordered model features of each view.
Use the system's `features.order` in [nominal.json](../PhotonID/models/nominal.json)
and the feature stride in `meta/schema`.

`pair_pho` indexes `pho_`. `pair_jet` indexes the concatenation
`jet02, jet03, jet04`, followed by `jetunsub02, jetunsub03, jetunsub04` in AuAu.
Use `meta/collections` for the recorded order. For an R = 0.4 analysis jet,
its local index is `pair_jet - n_jet02 - n_jet03`; keep indices in `[0, n_jet04)`.
The pair list retains candidate combinations; analysis cuts belong to TreeToHists.

The producer fills nominal `pho_bdt_score`, `pho_bdt_state`, `pho_id_class`
and `pho_iso_class`. Additional named scores appear after PhotonID scoring.
Check validity and filled-group metadata before using optional fields;
jet calorimeter-fraction equivalence remains unqualified.

## Simulation truth

| Branch group | Contents and indexing |
|---|---|
| [`truth_`](SCHEMA.md#truth-event-truth_) | Hard-scatter/minimum-bias truth vertex and denominator completeness |
| [`truthvtx_`](SCHEMA.md#truth-vertices-truthvtx_) | Native vertex IDs, embedding IDs, z and validity |
| [`truthpho_`](SCHEMA.md#truth-photons-truthpho_) | Truth photons, origin, isolation, signal labels and response class |
| [`truthjet02_`, `truthjet03_`, `truthjet04_`](SCHEMA.md#truth-jets-truthjetrr_) | Truth jets at R = 0.2, 0.3 and 0.4 |
| [`plink_`](SCHEMA.md#photon-reco-truth-links-plink_) | Reco/truth photon indices, distance, class and association policy |
| [`jlink_`](SCHEMA.md#jet-reco-truth-links-jlink_) | Reco/truth jet indices, distance and link class |
| [`truthmiss_`](SCHEMA.md#truth-photons-without-a-reco-candidate-truthmiss_) | Established embedded-photon miss records |

`plink_reco` indexes photons; `plink_truth` indexes truth photons.
Jet links use concatenated reco/truth collection indices in their respective
metadata order. Read link classes and capture validity with the indices:
an unavailable association is not evidence of a detector miss. DATA has no
truth-object population.

## Bookkeeping

| Branch group | Contents |
|---|---|
| [`recopho_`](SCHEMA.md#photon-capture-bookkeeping-recopho_) | Inspected, eligible, written and unclassifiable photon counts |
| [`jetview_`](SCHEMA.md#jet-views-jetview_) | Reco jet-view availability |
| [`truthcen_`](SCHEMA.md#truth-centrality-truthcen_) | Truth-photon capture, rejection and signal counts |
| [`truthjetcap_`](SCHEMA.md#truth-jet-capture-truthjetcap_) | Truth-jet container/capture status |
| [`ppg12_`](SCHEMA.md#ppg12-event-quantities-ppg12_) | pp weighting, leading-source and fill-count records |
| [`pdiag_`](SCHEMA.md#ppg12-diagnostic-occurrences-pdiag_) | Individual pp diagnostic occurrences and their weights |
| [`embdiag_`](SCHEMA.md#embedded-photon-diagnostic-occurrences-embdiag_) | Embedded-simulation diagnostic occurrences |
| [`jcon_`](SCHEMA.md#jet-constituents-empty-in-current-production-jcon_) | Jet-constituent fields; empty in current production |
| [`snap_`](SCHEMA.md#event-display-snapshots-empty-in-current-production-snap_) | Event-display fields; empty in current production |

Diagnostic records retain decisions and occurrences used by the analysis;
they are records rather than histogram bins. Integer-coded text is decoded
through `meta/codes`. Unavailable measurements retain NaN or their documented
state; branch presence alone does not establish a valid value.

Energies and momenta use GeV, vertex z uses cm, angles use radians, calibrated
times use ns and centrality uses percent. Raw tower times use ADC samples;
the branch reference gives the conversion and field-specific exceptions.

For working readers, see [quick_start.py](examples/quick_start.py) and
[quick_start.C](examples/quick_start.C).
