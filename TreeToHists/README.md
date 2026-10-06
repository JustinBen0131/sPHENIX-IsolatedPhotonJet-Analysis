# TreeToHists

Fills one ROOT histogram package per DATA collection or simulated sample from eventTree v1 files, with one selection
for spectra, event records, QA and the simulation response. C++/ROOT with bundled MIT-licensed nlohmann/json headers; the inputs are only
read.

```bash
./TreeToHists/tree_to_hists.sh local --settings TreeToHists/config/treetohists.cfg --inputs catalogue.json \
    --work work --out histograms --parallel 8
```

`local` plans, fills and merges on this machine. On SDCC the same plan runs as Condor jobs:

```bash
./TreeToHists/tree_to_hists.sh plan   --settings TreeToHists/config/treetohists.cfg --inputs catalogue.json --work WORK
./TreeToHists/tree_to_hists.sh condor --work WORK --logs LOGS --events-log /tmp/$USER/<campaign>/events.log \
    --setup "/opt/sphenix/core/bin/sphenix_setup.sh -n new"
condor_submit WORK/condor/fill.sub
./TreeToHists/tree_to_hists.sh merge  --work WORK --sample AuAu/data --out histograms     # once per sample
```

## The chain

```
TreeToHists/
  README.md                     this guide
  tree_to_hists.sh              build, plan, fill, merge, local, Condor, centrality, isolation
  config/
    treetohists.cfg             every choice: selection name, packages, triggers, overrides, storage axes
    auau_isolation_wp80.json    the frozen Au+Au isolation line
    auau_embedded_stitching.json  the canonical Au+Au embedded stitching (Jet12/20/30/40, PhotonJet12/20)
  src/
    TreeToHists.cc              the settings and resolved selection, the catalogue and plan, one fill job, the merge
                                into a package, the Au+Au centrality weights, main
    Isolation.cc                the Au+Au isolation package: population, derive, freeze, check
    FillPhotonJet.C             event counts, photons, ABCD regions, recoil jets and their Δφ sidebands, event records
    FillPhotonJetResponse.C     simulation: photon and joint (photon pT, xJ) responses, misses, fakes, leakage, Δφ
    FillPhotonQA.C              photon shapes, isolation, score, timing, event, trigger, calorimeter, generator,
                                matching and photon-ID QA
    PhotonJet.h                 what the fills share: the selection, the eventTree input, the truth definitions
    TreeToHists.h               what the two driver files share
```

`build/` (the compiled driver) is made on the machine that runs it and is not part of the package.

1. **Plan.** The catalogue names every sample and its files, or its published stage-2 groups (`GROUP/CURRENT` with a
   hash-bound `MANIFEST.json`). Each file's system, DATA/simulation family, source identity and score binding is checked
   against the resolved selection; Au+Au simulation must bind its 16 centrality weights, their evidence, the DATA
   trigger and a qualified truth-vertex fit. `WORK/plan.json`, `WORK/selection.json` and `WORK/jobs.txt` result.
2. **Fill.** One job per eventTree file writes its spectra, QA, extension and (simulation) response products into
   `WORK/<system>/<sample>/<shard>/<job>/` and records the input's sha256 (checked against the group manifest when
   there is one). Reuse of `done.json` requires the same input, plan and reduction code. A changed recipe needs a
   new work directory; progress is written every 10,000 events.
3. **Merge.** A sample's jobs are summed into `<out>/<system>/<sample>.root` with a `<sample>.json` sidecar (sha256,
   sources, histograms, event records). A package is never replaced.

## Provenance and incomplete populations

Direct file inputs are SHA256-pinned at planning and rechecked before a fill (including reuse of a completed fill),
after reading and before merging. External vertex-fit files receive the same checks. Plans and fill receipts also
bind the reduction code, so old fills cannot silently be published under a newer implementation.
Source contributions, event records, Sumw2 and source identities remain in the
package; split batches merge through this driver, rather than generic histogram addition.

A catalogue sample may give `population_sources`, the expected number of logical DST sources. Packages record
`COMPLETE`, `PARTIAL` or `UNKNOWN` coverage without changing the denominator. A source-PASS qualification subset is
usable even when its whole run is not published, but must not be presented as complete.

Stage-2 files may have the nine calorimeter-fraction values reserved for later augmentation. With
`Fill.RequireStage3Fractions: 0` they remain usable, and `availability.fractions` records `stage3_not_filled` when
appropriate. Fraction-dependent QA is unavailable; unrelated branches, source identities and selection checks still
apply. Setting the flag to 1 requires the fractions.

`meta/package.render_selection` records the actual Profile and resolved photon-ID/isolation definition used by the
fills, their scientific source-code hashes, QA classes, event cuts and response definitions. FinalAnalysis compiles its
captions from this information. Nominal and optional-model selections remain distinct and hash-bound.

The public entrypoint was exercised with three completed pp DATA sources, six pp/AuAu software-reference sample types
and a pp photon DI sample on ROOT 6.34 (macOS), and with the same pp DATA sources on SDCC ROOT 6.32 (ana.568).
Full-population data and scientific acceptance remain separate gates;
development qualification reports and tests are kept outside the public distribution.

## Settings

`config/treetohists.cfg` (ROOT TEnv, `Key: value`, `#` starts a comment) holds the selection name, optional
replacement packages, the DATA trigger of each system, optional overrides of the working-point lines, isolation radius,
threshold and sideband gap, an augmented score or a sealed PhotonID release (for example the ETcentFine model), and the
storage axes. The nominal selection is the frozen `PhotonID/models/nominal.json` with `config/auau_isolation_wp80.json`:
Au+Au tight above WP70(c), non-tight 0.1 < s < WP80(c), R = 0.3 isolation below 3.22 - 0.0246 c GeV with a 0.8 GeV
sideband; p+p the PPG12 lines with the stored R = 0.4 witness. The selection's sha256 is the profile of every product
and package, so DATA and simulation of one result always share it. A variation is a copy of the file, never new code;
an unknown key is an error.

DATA events need the nominal trigger by exact name in each source run (p+p `Photon 4 GeV + MBD NS >= 1` with |z| < 60 cm
in the run periods 47289-51273 and 51274-53999; Au+Au `Photon 12 GeV + MBD NS >= 2, vtx < 150 cm` with the minimum-bias
classifier and |z| < 10 cm). The package holds only that population: the trigger and the sources are chosen here, once,
and FinalAnalysis reads the package as it is.

## Au+Au centrality weights

Au+Au simulation is weighted per 5% slice of 0-80% to the DATA centrality distribution: `w(c_i) = p_DATA(c_i) /
p_SIM(c_i)`. After the DATA package exists:

```bash
./TreeToHists/tree_to_hists.sh centrality --settings TreeToHists/config/treetohists.cfg --inputs catalogue.json \
    --data histograms/AuAu/data.root --out weights.json
```

DATA is the selected trigger's population in every source where the trigger is enabled (by name in the source's own
trigger table); simulation is `hist/SIM/h_centrality` of the catalogue's Au+Au photon samples, summed with their cross
section per generated event. The 16 weights, their sha256 and every count go into `weights.json`; they become the
catalogue's `centrality_weights` (with `weights.json` as their evidence) before the simulation is planned. An empty
slice on either side is refused.

## Au+Au isolation package

The Au+Au isolation cut (`config/auau_isolation_wp80.json`) is a freeze of this derivation:

```bash
./TreeToHists/tree_to_hists.sh isolation population --out population.root EVENTTREE...    # truth-signal photons
./TreeToHists/tree_to_hists.sh isolation derive --out package.json population.root...
./TreeToHists/tree_to_hists.sh isolation freeze --out freeze.json package.json            # R = 0.3, 80%, 0.8 GeV gap
./TreeToHists/tree_to_hists.sh isolation check  --out check.json freeze.json population.root...
```

The population holds the truth-signal photons (|eta| < 0.7) of minimum-bias embedded events with |z| < 10 cm in 0-80%,
with their reconstructed R = 0.3 and 0.4 isolation. For every photon-pT group and 5% slice the cutoff at 70, 80 and 90%
efficiency is read from the cumulative isolation spectrum; a constant across pT per slice and a line across centrality,
`threshold(c) = intercept + slope c`, follow. The freeze takes the R = 0.3 line at 80% with the sideband gap, bound to
the derivation's sha256; `check` measures a frozen line on new simulation.

## Systematic variations

`Fill.ResponseVariation` fills the simulated response with a varied reco energy scale and resolution (`jet_scale`,
`jet_smear`, `photon_scale`, `photon_smear`, `seed`), for the jet energy scale and resolution systematics of
FinalAnalysis. The variation is part of the plan, so a varied run has its own work directory and packages.

## The catalogue

```json
{"schema": "PhotonJetInputCollectionV1",
 "reference_cross_section_pb": {"pp_photon_signal": 130.4461},
 "samples": [
   {"system": "PP", "name": "data", "family": "data", "groups": ["trees/PP/run47291"]},
   {"system": "PP", "name": "photon20", "family": "photonjet", "files": ["trees/PP/photon20/part0.root"],
    "normalization": {"source_factor": "family_reference", "reference": "pp_photon_signal", "generated_events": 10000000}}]}
```

Paths are relative to the catalogue. Simulation families are `photonjet` and `inclusivejet`, one named source
population per sample. The cross section per generated event is recorded (`pb_per_generated_event`) for FinalAnalysis
to apply once; the fills are not normalised. Au+Au simulation also binds `centrality_weights`, `centrality_evidence`,
`centrality_data_trigger` and `vertex_reweighting` (`path`, `sha256`, `histogram`, `coordinate: truth_vtx_z`,
`selection_sha256`, `centrality_weights_sha256`, `evidence`); its fill weight is `w_event × vertex_fit(truth z) ×
centrality weight`.

For pp simulation, the producer's active PPG12 weight state defines source admission. Inactive events remain in the
tree for accounting but enter neither measurement, response nor QA. `ppg12_weight_lane_code` holds the producer's
lane component code: 1 photon+jet SI, 2 photon+jet DI, 3 inclusive-jet SI, 4 inclusive-jet DI, 0 when inactive. The
expected code follows from the source's `meta/source` dataset (`isSim`, `isSimInclusive`) and its `si_di_role`
(`SI` or `DI`); a pp simulation file without that role is refused at planning. Accepted rows must carry that code,
one weight application and `w_event = w_slice × w_vertex × w_sidi × w_period`. This includes the truth-vertex factor
(both vertices for DI), SI/DI mixture and run-period exposure, once each. The DI component of a slice is its own
catalogued sample with its own generated-event count. The slice factor already contains the PPG12 cross-section
ratio: photon slices relative to Photon20 (130.4461 pb), inclusive-jet slices relative to Jet50 (7.3113 pb). The
catalogue therefore supplies that common family reference over its qualified generated-event count, rather than
multiplying by each slice's cross section again.

Evidence fields identify the analysis's externally qualified weight derivations; a nonempty label does not prove
their scientific validity. The production catalogue must supply actual exposure/count reconciliation and
selection-matched Au+Au centrality and truth-vertex qualification. Historical weights or unit vertex fits are
software fixtures and cannot certify the final physics result.

The Au+Au embedded samples take the canonical stitching, `"normalization": {"source_factor": "auau_embedded_stitching"}`,
from `config/auau_embedded_stitching.json`: Jet12/20/30/40 own the events whose leading R = 0.4 truth jet lies in
12-21, 21-31, 31-41 and above 41 GeV, with ownership-effective cross sections over the sample's events in the same
window; PhotonJet12/20 own the events whose generator-filter photon lies in 12-20 and above 20 GeV, with the accepted
1/ΣL per post-filter event (cross section over all the sample's events). Every fill leaves out the events a sample does
not own (they belong to its neighbour); each fill job counts the events of the denominator and the merge sets
`pb_per_generated_event` from them. The generator spectra keep each sample's whole range, for the overlap at every
boundary. Known events rejected before the photon-ownership decision still count in a photon's all-events
denominator. A retained event without its required ownership witness is an error.

## The package

- `measurement/<class>/`: events, photons, ABCD yields; `jet5`, `jet7`, `jet10`, `jet12` with the recoil jets and the
  near/far Δφ sidebands (A and C regions).
- `response/<class>/jet<t>/` (simulation): photon and pair responses, misses, fakes by category, leakage, truth and
  combinatoric Δφ.
- `qa/`: photon shapes, isolation, score, timing, event, trigger, calorimeter, generator and matching histograms.
- `records/events`, `photons`, `jets`: every region-leading photon and every jet around the A and C leaders with its
  event (64-bit references), for the whole-event bootstrap; `records/sources` and `records/source_contributions`: each
  source's identity, hash, metadata and sparse measurement and response cells, for the source bootstrap.
- `meta/package`, `meta/inventory`: definitions, axes, provenance and contents. Every histogram keeps explicit Sumw2.

## Condor on SDCC

The job environment is not inherited (`getenv = False`): each job sources the `--setup` script. The Condor event log
goes in `/tmp/$USER/<campaign>/`; stdout and stderr stay persistent outside TG, sharded below 1000 files per directory;
there is no file transfer and no runtime limit. Build once on the submit host (`tree_to_hists.sh build`): a job never
compiles. The printed `condor_submit` command is run by the analyst.

## Status

The C++ migration matched the earlier Python-coordinated reference products, including histograms, event records,
source contributions and selection. The subsequent weighting audit corrected pp admission: recorded events rejected
by source stitching no longer enter selected products. Those older pp products are comparison evidence, not a golden
answer for the corrected population. Au+Au centrality and isolation derivations retain their numerical reference
checks; accepted populations, weight products, Sumw2 and source conservation have separate regressions.
Nothing in TreeToHists needs Python.

See `../FinalAnalysis/README.md` for the portable end-to-end command and the distinction between software qualification and complete-population physics results.
