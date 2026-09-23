# Isolated photon–jet analysis

One reconstruction pass writes model-independent facts. Identification and measurement choices remain downstream.

```text
sPHENIX PhotonClusterBuilder dependency
                 ↓
DST → TreeProduction → canonical trees
                              ↓
                           PhotonID → candidate-keyed model scores
                              ↓
                          TreeToHists → additive histogram packages
                              ↓
                         FinalAnalysis → numerical results → plots
```

This is a locally refined implementation for human review. It has **not** been compiled against sPHENIX, exercised on ROOT data, or accepted as scientifically equivalent to the production reference. Configuration bindings deliberately remain unresolved. Do not interpret implemented source or a successful syntax check as production readiness.

| Directory | Owns | Contract/configuration |
|---|---|---|
| [TreeProduction](TreeProduction/README.md) | Reconstruction steering, one producer, identities and stored primitives | `TreeProduction/config/tree_production.yaml` |
| [PhotonID](PhotonID/README.md) | Shared feature definitions, training/import bindings and immutable score sidecars | `PhotonID/model_registry.yaml` |
| [TreeToHists](TreeToHists/README.md) | One batched event loop, measurement selection and additive sufficient statistics | `config/measurement.yaml`, `config/samples.yaml` |
| [FinalAnalysis](FinalAnalysis/README.md) | Purity, corrections, unfolding, normalization and explicitly limited statistical inference | Same measurement configuration; packages only |
| [validation](validation/README.md) | Historical-field review, migration admission and content comparison | `validation/field_mapping.csv` |
| `coresoftware/offline/packages/CaloReco/` | Proposed generic PhotonClusterBuilder dependency change | Two source files; separate upstream review |

The producer remains one `PhotonJetTree` class. Its event flow is in `PhotonJetTree.cc`; methods are grouped in `internal/*.cc`. These are implementation files, not separate stage classes. Production setup never belongs in the macro; physics never belongs in ROOT serialization.

Changing a model uses retained features. Changing working points uses stored scores. Changing unfolding within retained histogram support uses packages. Plotting reads persisted results. Missing reconstruction primitives are the reason to return to DST.

## Review and integration boundary

Read each stage README before its source. First reconcile typed records, configuration and producer event flow; then inspect reconstruction and physics capture; finally inspect downstream contracts and migration. The migration ledger is entirely **REVIEW**, not accepted coverage.

Before integration, resolve build/install registration and external API compatibility, complete archived input reconstruction and topocluster setup, bind conditions/models/normalization, and approve the measurement and uncertainty contracts. Representative lanes must subsequently pass real compile/runtime checks and normalized scientific-content comparison. No ROOT byte-equality requirement is imposed.

The current histogram serialization is a hash-bound **NPZ + JSON** pair, with internal sample/period/SI-DI partitions. A collaborator-facing ROOT histogram container is not implemented. That packaging decision must be resolved before promising that interface.

Development requires Python 3.10+; stages use NumPy, uproot and PyYAML. Training additionally uses XGBoost, scikit-learn and PyROOT/TMVA; plotting uses matplotlib. C++ requires the reviewed sPHENIX release, ROOT, yaml-cpp and FastJet interfaces. No synthetic dependency environment or runtime qualification is supplied here.
