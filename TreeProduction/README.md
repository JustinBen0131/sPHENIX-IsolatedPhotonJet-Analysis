# TreeProduction

Consumes one physical DST source bundle and writes model-independent ROOT tables. One `PhotonJetTree` SubsysReco serves pp/AuAu and DATA/SIM through typed configuration. This source is uncompiled and unvalidated on sPHENIX.

| File | Ownership |
|---|---|
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
