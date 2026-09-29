# TreeToHists

Produce one portable ROOT histogram package per DATA collection or exact SIM sample from final single-`eventTree` inputs. Photon selection, isolation, spectra, QA and response share one resolved definition. Input trees are read only.

## Commands

```bash
./make_hists.sh runPP --inputs collection.json --out ./histograms --plan
./make_hists.sh runAll --inputs collection.json --out ./histograms
./make_hists.sh runAll --inputs collection.json --out ./histograms --resume
```

| Command | Population |
|---|---|
| `runPP` | pp DATA |
| `runAuAu` | AuAu DATA |
| `runAllData` | Both DATA collections |
| `runAllSim` | Each listed SIM sample independently |
| `runAll` | All listed DATA and SIM |
| `--help` | Options and invocation help |

Use `--sample photon20` to restrict a command to an exact sample name. `--plan` validates configuration and input metadata/hashes without producing histograms. `--resume` reuses only unchanged products with matching receipts. A changed input, code version or definition requires a new output directory. The default output is `./histograms`, relative to the caller. Set `PHOTONJET_TREES` to your collection JSON to omit `--inputs`; no personal storage location is embedded.

DATA selection is applied before measurement event counts and selected photon/jet QA. pp uses the exact scaled trigger `Photon 4 GeV + MBD NS >= 1` and |vertex z| < 60 cm; run periods are 47289–51273 (0 mrad) and 51274–53999 (1.5 mrad), restricted to the source's period when bound. AuAu uses `Photon 12 GeV + MBD NS >= 2, vtx < 150 cm`, the minimum-bias classifier pass and the tighter offline |vertex z| < 10 cm cut. Trigger bit positions are resolved per run. Raw trigger QA retains its inclusive population; no parallel Photon10 measurement is produced.

Use the Python environment containing ROOT, NumPy, uproot, awkward and PyYAML. `PHOTONJET_PYTHON` selects that interpreter without modifying its packages. Each worker has a hard deadline (`--timeout`, default 3600 seconds), its own log, progress record and output receipt. The command does not submit farm jobs.

## Input collection

A maintained catalogue binds exact sample identities to published files or group directories. Paths are relative to the catalogue; whitespace in ROOT file paths is currently unsupported by the retained macro interface. Group `CURRENT` pointers are resolved once, and their published hashes are checked.

```json
{
  "schema": "PhotonJetInputCollectionV1",
  "samples": [
    {"system": "PP", "name": "data", "family": "data",
     "groups": ["trees/PP/run47291", "trees/PP/run47292"]},
    {"system": "PP", "name": "photon20", "family": "photonjet",
     "source_sample": "run28_photonjet20", "files": ["trees/PP/photon20/part0.root"],
     "normalization": {"source_factor": "cross_section_over_generated",
                       "cross_section_pb": 133.0, "generated_events": 1000000}}
  ]
}
```

The numbers above illustrate the format; use the production's recorded cross section and generated count. SIM families are `photonjet` and `inclusivejet`; retain distinct names for DI samples. SIM inputs from different source populations cannot be mixed in one package. Each original logical source appears once. Use stage-2 complete input groups for production; loose files are supported for local qualification.

AuAu SIM also requires `centrality_weights` (16 positive values, one per 5% slice), `centrality_evidence` identifying their derivation for that production, and `centrality_data_trigger` containing the exact DATA trigger name. The trigger binding must match the nominal selection. `FinalAnalysis/src/internal/centrality.py` defaults to Photon 12 and retains the established named-trigger DATA/SIM derivation and normalization rules. Never use historical or unit weights as a substitute for current production evidence. For nominal AuAu SIM, the collection also binds `vertex_reweighting`: `path`, `sha256`, `histogram`, `coordinate: truth_vtx_z`, the resolved `selection_sha256`, `centrality_weights_sha256` (the canonical JSON digest of the 16 weights), and `evidence` identifying the vertex correspondence and joint closure qualification. The fit is derived once against the Photon12/minimum-bias DATA vertex distribution, then reused. A pp curve or the historical AuAu reconstructed-vertex ratio is not an AuAu truth-vertex calibration.

The AuAu fill weight is `w_event × vertex_fit(truth_vtx_z) × centrality_weight`. The input's saved vertex factor must be unity; existing vertex-weight metadata also prevents repeat application. The ROOT fit and all input files are hash-checked before and after execution. pp uses its already composed PPG12 `w_event` with no additional vertex factor. Cross-section normalization stays separate for FinalAnalysis to apply exactly once. A missing AuAu fit/closure binding refuses nominal reduction; the software interface alone does not qualify that fit.

## Output contract

```text
histograms/
  manifest.json
  progress.json
  PP/data.root
  PP/photon20.root
  AuAu/data.root
  AuAu/photon12.root
  AuAu/photon20.root
```

Only complete sample candidates are published. The manifest binds each published file's hash and available diagnostics. Hidden `.work/` holds resumable source products, logs and receipts; copying the sample ROOT files and manifest is sufficient for offline reading. Keep `.work/` while a reduction may need resuming. The coordinator holds an exclusive local lock, rejects output collisions and checks hashes after worker completion.

Within each sample:

- `measurement/<centrality>/` contains event/photon counts and ABCD yields; `jet5`, `jet7`, `jet10`, `jet12` contain recoil and near/far angular-sideband histograms.
- AuAu DATA also contains `triggers/photon12/measurement/`, with the same complete inventory and Sumw2. The top-level measurement is an alias of the nominal stream. `records/triggers/<trigger>/` retains its event, photon and jet tables. Empty complete sources remain represented in `records/sources`.
- `response/<centrality>/jet<threshold>/` contains SIM photon/pair responses, misses, fake categories, leakage and truth/combinatoric angular distributions.
- `qa/` contains photon shapes, score/isolation correlations, timing, event/trigger/calorimeter quantities, generator/truth diagnostics and matching. Stage-3 calorimeter fractions use a component axis for EMCal, IHCal and OHCal at each of R=0.2, 0.3 and 0.4.
- `records/events`, `records/photons`, `records/jets` retain event contributions with 64-bit event references. `records/sources` retains source identity, hashes, run-specific trigger metadata and normalization information. `records/source_contributions` retains sparse per-source measurement and SIM response cells and variances, plus source-specific trigger and centrality counts. This allows exact run selection and source bootstrap after packaging.
- `meta/package` and `meta/inventory` describe definitions, axes, provenance, actual contents and diagnostic availability.

Every booked histogram has explicit Sumw2, including DATA and empty selections. Missing stage-3 fraction payloads are marked unavailable and omitted; `require_stage3_fractions: true` in the histogram configuration instead refuses them. Undefined fraction values are omitted from distributions, while valid negative values and histogram flow bins are retained.

The intended inventory is 69/244 objects for pp/AuAu DATA, 170/469 for photon SIM and 168/463 for inclusive SIM. The AuAu DATA count includes 69 measurement objects for Photon12 in addition to the nominal aliases. Current implementation omits the NPB-tagged timing and semantic-trigger diagnostics until their run-calibration/selection contracts are pinned. These omissions and absent fractions are explicit, so the saved object count can be lower than the target. Existing raw trigger-bit and untagged timing QA remain available.

## Definitions and binning

`config/selection.yaml` selects the scored model, matching working points, isolation and photon/jet cuts. Defaults use the original AuAu nominal BDT with centrality-dependent WP70/WP80 and 0.1 non-tight floor, AuAu R=0.3 80% isolation, and the original pp PPG12 definitions. Changes produce a different selection hash shared by all DATA/SIM outputs. Alternative models require their own score provenance and working-point package; they cannot silently replace nominal scores.

For a validated PhotonID release, use a selection YAML like this and pass
`--selection optional.yaml` to the same public command:

```yaml
version: 1
name: optional-model
systems:
  auau:
    photon_id_selection: ../models/AuAu/MLP_v1/selection.json
```

The path is relative to the YAML file. That system's score/state/cut overrides
must be omitted: the sealed release owns its single score, applicability domain
and routed working points. Inapplicable scores are NaN; nominal state never
decides optional applicability. Isolation retains its nominal definition unless
explicitly varied. The same selected release must be present in DATA and SIM.

`config/histograms.yaml` defines storage axes. Photon measurement bins are 1 GeV from 15 to 35, supporting purity/result intervals **17–23 and 23–33 GeV**, with 15–17 and 33–35 guards. QA covers 10–40; response support is 10–40 reconstructed and 5–40 truth. The retained xJ grid supports the candidate result edges **0.35, 0.50, 0.80, 1.00, 1.25, 1.49** and existing reference comparisons. Result xJ bins remain provisional until response statistics, stability and closure qualify them.

The angular sidebands retain both A/C photon regions and near/far delta-eta categories in DATA and SIM. FinalAnalysis owns the combinatoric estimate, purity correction, unfolding, uncertainties and figure production; this package retains their inputs without fitting those results during reduction.

For offline AuAu DATA selection, `package.event_table_arrays(path, xj_edges, all_jets=True, data_trigger="photon12", runs=[...])` returns the matching retained event contributions with remapped photon/jet references. `package.measurement_histograms(path, data_trigger="photon12", runs=[...])` returns all measurement histograms and their variances, plus per-source inclusion records. Omit the trigger to use the package's nominal choice and omit `runs` to include every source where it was enabled. Use these same contributions for ABCD, recoil and combinatoric subtraction; filtering only a final plotted spectrum is insufficient. General detector QA retains its original recorded-event population and is labeled separately.

## Code and build

- `scripts/make_hists.py` coordinates source-bound jobs, receipts, resume and publication.
- `src/PhotonJetHists.h/.cc` is the C++ lifecycle.
- `src/internal/Selection.cc` validates histogram storage settings.
- `src/internal/Histograms.cc` dispatches maintained fillers and additional QA.
- `src/internal/Package.cc` consolidates histograms and remaps event references.
- `src/internal/Centrality.cc` exposes the maintained centrality-sum implementation.
- `src/internal/Types.h` shares the job contract.

Existing `FillPhotonJet.C`, `FillPhotonQA.C`, `FillPhotonJetResponse.C` and shared selection/truth readers remain the single implementations used by the library and compatibility macro entry points. `FinalAnalysis/src/internal/inputs.py` provides array readers for both legacy outputs and new packages. `source_views(path)` exposes each original source within a packaged response for bootstrap without unpacking files. FinalAnalysis owns collection loading and end-to-end analysis qualification.

For a compiled library, from an empty build directory:

```bash
/path/to/TreeToHists/src/autogen.sh --prefix="$PWD/install"
make -j2
make install
export PHOTONJET_HISTS_LIBRARY="$PWD/install/lib/libPhotonJetHists.so"
```

Use `.dylib` if the platform installs that suffix. Without this variable, the command loads the same C++ sources through the installed ROOT interpreter. The compiled library hash is bound into each plan. No library or dependency installation occurs automatically. Autotools generation requires autoreconf/libtool in the build environment; direct native compilation and ROOT interpretation are separately qualified locally. Local development tests and build products are ignored and are not needed to run the package.

## Qualification boundary

Local compilation, reference comparisons and sample-package checks establish software behavior. Production readiness additionally requires recovered input collections, their stage-2 acceptance evidence, current AuAu centrality/normalization evidence, farm execution qualification and downstream analysis regressions. Inclusive/DI runtime coverage, full IAN reproduction and scientific acceptance must be established from their actual inputs; they are not implied by a successful file write.
