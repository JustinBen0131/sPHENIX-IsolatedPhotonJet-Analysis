# TreeProduction

Consumes one physical DST source bundle and writes model-independent ROOT tables. One `PhotonJetTree` SubsysReco serves pp/AuAu and DATA/SIM through typed configuration. This source is uncompiled and unvalidated on sPHENIX.

## Producing trees

From this directory, the only public entrypoint is:

```bash
./produce_trees.sh --help
./produce_trees.sh runPP
./produce_trees.sh runAuAu
./produce_trees.sh runPhotonJetSim
./produce_trees.sh runInclusiveJetSim
./produce_trees.sh runEmbeddedPhotonJetSim
./produce_trees.sh runEmbeddedInclusiveJetSim
./produce_trees.sh runAllData
./produce_trees.sh runAllSim
./produce_trees.sh runAll
```

The explicit Bash shebang works when launched from fish. Invocation also works
from another directory by giving the path to `produce_trees.sh`.

| Command | Lane | Existing profile |
|---|---|---|
| `runPP` | `pp_data` | `pp_data` |
| `runAuAu` | `auau_data` | `auau_data` |
| `runPhotonJetSim` | `pp_photon_sim` | `pp_photon_sim` |
| `runInclusiveJetSim` | `pp_inclusive_sim` | `pp_inclusive_sim` |
| `runEmbeddedPhotonJetSim` | `auau_photon_embedded` | `auau_photon_embedded` |
| `runEmbeddedInclusiveJetSim` | `auau_inclusive_embedded` | `auau_inclusive_embedded` |
| `runAllData` | `runPP`, then `runAuAu` | Same primitive routes |
| `runAllSim` | Four simulation commands above, in order | Same primitive routes |
| `runAll` | `runAllData`, then `runAllSim` | Same primitive routes |

**Current status: all six population bindings are null.** No canonical source
catalogue, input manifests, batch launcher or installed producer binding exists
in this repository. Every production command therefore exits nonzero with the
exact missing binding location. Aggregates check every requested binding before
launching anything. This interface does not supply missing production authority.
The two archived-DI profiles remain unchanged and have no public command here.

```text
produce_trees.sh COMMAND
  -> expand aggregates into primitive lanes
  -> scripts/produce.py: preflight all populations and output locations
  -> one execute_source() for every selected source
  -> macros/Fun4All_PhotonJetTree.C: Job -> Plan -> Fun4All
  -> Production.cc -> PhotonJetTree -> ROOT product
  -> check terminal metadata/accounting before reporting success
```

The one private Python helper avoids brittle shell JSON/quoting and uses only the
standard library for planning. It performs no YAML policy parsing, DST layout
interpretation, reconstruction, weighting or event physics. The macro and C++
implementation are the only production path, including development runs.

### Population versus profile

`config/tree_production.yaml` defines **how** a source is reconstructed.
`config/production_populations.json` binds **which** sources belong to each lane.
The latter is an operational JSON document so help, planning and missing-binding
errors need no YAML/ROOT Python packages. Its fields are:

- `campaign`: reviewed production revision; deliberately null until supplied.
- `output_root`: deterministic storage root, currently `../output/trees`.
- `production_config`: default existing YAML policy, `tree_production.yaml`.
- `producer_library`: exact installed `libPhotonJetTree.so`; currently null.
- `lanes.<lane>.profile`: the centralized lane-to-profile binding shown above.
- `lanes.<lane>.manifest` and `manifest_sha256`: reviewed source manifest and
  its exact byte hash; both currently null. Hashes are populated by the future
  manifest/build preparation, not manually supplied by an operator on the CLI.

Relative paths in this binding file resolve against its own directory. Each
future manifest must contain `schema_version: 1`, its exact `campaign`, its
`lane`, and a nonempty `sources` array. Each source record must explicitly supply:

| Field | Meaning |
|---|---|
| `source_file_ordinal` | Stable, unique nonnegative ordinal in this manifest |
| `input_list`, `input_list_sha256` | Exact existing one-bundle list and its digest |
| `run`, `segment` | Source metadata; never guessed from a filename |
| `first_entry`, `event_count` | Exact starting entry and positive expected count |
| `sample`, `period`, `si_di_role` | Producer provenance; empty period/role is allowed only when explicitly supplied |
| `production_config`, `production_config_sha256` | Optional prepared per-source YAML snapshot and required digest when overriding the default |

Manifest-relative paths resolve against the manifest directory. Input-list
contents remain the existing `Production.cc` contract; use qualified input
paths. If those contents or YAML contain relative paths, the production working
directory is always `TreeProduction`, independent of the caller's directory.
One list represents one physical bundle. This frontend does not invent lists,
split sources, discover runs, calculate metadata or silently deduplicate rows.
Duplicate ordinals/list digests, empty populations and bad hashes are errors.
All manifest rows are checked even with a development source selection.

Source identity uses the existing manifest-hash-plus-ordinal route. It is not
renumbered or rebound to a smaller manifest for development runs. This hashes
the population and input-list bytes, not the entire DST payload. Exact source
catalogue qualification remains a production-preparation responsibility.

Prepared per-source YAML uses the existing parser; this is the binding point
for already-supported run-dependent settings and provenance. The frontend does
not fill missing calibrations, patch science policy, infer a source/library
build relationship, or create a second configuration resolver. Deployment must
provide truthful runtime source/library/macro hashes and required conditions.

### Development restrictions and output

```bash
./produce_trees.sh runAuAu --source 0 --events 10 --dry-run
./produce_trees.sh runAuAu --source 0 --events 10
```

These commands still fail while the population is unbound. `--source N` means
the manifest ordinal, not row position, and is valid only for primitive
commands. `--events N` caps the declared count **per selected source**, never
expands it; without `--source`, all sources remain selected. Zero/negative counts
are rejected. No `--submit`, alternate macro, automatic resume, or standalone
test producer is provided. An entry offset, if required, belongs to the reviewed
manifest; a public `--first-entry` override is intentionally deferred.

`--dry-run` validates operational bindings and prints the plan without creating
directories, importing ROOT/uproot, or launching a child. It does **not** certify
the C++ configuration, conditions, installed environment or input DST contents.

The default output layout is:

```text
output/trees/<campaign>/<lane>/source-<ordinal>/
  trees.root
  production.log
  invocation.json
  completion.json       # written only after terminal checks pass

output/trees/<campaign>/development/<lane>/source-<N-or-all>_events-<N-or-manifest>/
  source-<ordinal>/...  # same route, explicitly a restricted product
```

Each lane execution directory also captures `population.json`. Logs and exact
invocations stay with the output. Existing execution directories are refused;
there is no silent overwrite, skip or retry. Partial failures remain available
for inspection. A fresh campaign/restriction is an explicit operator decision,
not an automatic fallback. Paths alone do not encode source semantics: the
manifest, invocation and producer metadata carry the identities and ranges.

Execution is currently serial/local. It requires ROOT on `PATH`, the configured
installed producer library and headers, and `uproot` for the terminal check.
The configured library directory is prepended to the child loader paths; the
output's library/macro/configuration hashes must agree with the invocation.
Missing calibration or external API bindings fail through the existing C++
path; the frontend preserves the child's exit status and shows a short log tail.

Success requires child exit zero **and** producer `completion_status=complete`,
matching source/provenance, a completed source row, consistent event/rejection
counts and the declared processed count. Merely finding `trees.root` is never
success. These terminal checks are operational completeness, not science
acceptance or proof of complete canonical-population membership. Development
success is explicitly labelled and never reported as full population closure.

### Frontend checks

```bash
bash -n produce_trees.sh
python3 -B -m unittest discover -s tests -v
```

The tests use temporary metadata fixtures and mocked process/terminal-reader
boundaries. They do not execute ROOT, create a fake sPHENIX environment, read
DSTs, or establish runtime/scientific equivalence.

| File | Ownership |
|---|---|
| `produce_trees.sh` | Only public command entrypoint; primitive/aggregate dispatch |
| `scripts/produce.py` | Private population planning, single-source invocation and terminal checks |
| `config/production_populations.json` | Canonical population/profile bindings and operational output/deployment locations |
| `config/tree_production.yaml` | Eight production profiles, input nodes, reconstruction choices, payload bindings and storage policy |
| `macros/Fun4All_PhotonJetTree.C` | Resolve job, configure conditions, create server/producer, register reconstruction/inputs, run, End |
| `src/PhotonJetTree.h` | One producer interface and state; aliases refer to Types.h |
| `src/PhotonJetTree.cc` | Lifecycle and visible event-processing order |
| `src/internal/Types.h` | Plain configurations, records, enums and stable identities |
| `src/internal/Production.h/.cc` | YAML/job/input resolution, CDB, reconstruction registration and input managers |
| `src/internal/Event.cc` | Source/event identity, triggers/scalers, vertex, MBD, centrality, calorimeter/UE state and producer weights |
| `src/internal/Photons.cc` | Candidates, timing, shower views/cells and separate pp/AuAu isolation semantics |
| `src/internal/Jets.cc` | Raw/JES binding, active area, retained encounter order and provisional pairs |
| `src/internal/Truth.cc` | Embedded primary census, HepMC ancestry, truth isolation, dominant-primary photon associations and jet matching |
| `src/internal/Output.cc` | ROOT booking/filling, configuration/provenance and completion only |

## Stored product

Eighteen declared tables: Sources, Events, UpstreamRejectedEvents, TriggerRunInfo, TriggerScalers, Photons, PhotonShowerViews, PhotonCells, Isolation, IsolationConstituents, Jets, PhotonJetPairs, TruthVertices, TruthPhotons, TruthJets, PhotonTruthLinks, JetTruthLinks and WeightComponents.

`metadata` contains file-level key=value provenance, `configuration` the supplied YAML, and `completion` terminal accounting. Only `completion_status=complete` is consumable. Write failures propagate, existing files are refused, and aborted jobs are marked. This is not an atomic publish/validation protocol; close/failure injection still needs runtime checks.

Every producer encounter is retained, including zero-candidate events. Upstream observer accounting is separate. One invocation accepts one physical source bundle; multi-file source transitions are intentionally unsupported. Manifest-based source identity requires an ordinal. Input hashes and exact run/entry ranges must be supplied.

Candidate storage has explicit ET/eta and reconstructed-vertex domains. These remain migration-review items, particularly for native simulation. Event retention alone does not establish complete object acceptance. No model score or working point decides retention.

Reconstructed jets retain finite nonnegative calibrated objects in container encounter order, restarting per view/radius. Raw pT must be finite/nonnegative, area valid, and calibrated ids bound to raw ordinals. JES is applied once upstream, never here or to truth jets. Per-view availability distinguishes valid empty, missing and not applicable. Provisional pair witnesses use absolute wrapped delta-phi, corrected jet pT / photon ET, inclusive 7π/8 recoil, photon_rank=0 and retained jet ordinal.

Truth isolation sums embedded primary **transverse energy**, with the merged core subtracted. A finite stored isolation value and a complete input census are separate facts. Native signal classification retains prompt classes below 3, including -1/0; generator association validity and archived hard-event ownership are separate witnesses. Photon association uses the maximum-energy primary identity, not nearest delta-R. Jet matching keeps its deterministic one-to-one rules and considered edges; incomplete capture cannot create a certified miss.

AuAu UE arrays, flow mode, v2, psi2 and failure witnesses are stored. Isolation axes are retained. Signed SUB1 towers and pp signed topocluster sums follow separate definitions; candidate ET is subtracted once. Cells, isolation constituents and Cartesian pairs stay provisionally enabled until the field ledger proves a lossless alternative. Always-NaN working-point placeholders were removed from this unaccepted base schema.

## Integration still required

- Supply real build/install rules for the library and its exported macro includes; compile all split definitions against the chosen release.
- Reconcile external APIs, especially JetCalib legacy-mode control, centrality, calorimeter status, truth evaluators, tower encoding and TMVA.
- Complete archived G4→waveform/tower registration with its exact random sequence. This lane currently fails explicitly.
- Bind/register the pp topocluster reconstruction when that node is absent from input.
- Bind CDB tags, centrality/ZS/JES payloads, vertex weights and runtime source/library provenance. No ambient or guessed payload is acceptable.
- Validate all eight profiles and compare direct output with the reviewed migrated reference before acceptance.

Never put tight selection, ABCD, purity, leakage, final recoil selection, response construction, unfolding, histogramming or plotting here. Generic optional model support in the upstream builder is disabled by this producer.
