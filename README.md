# Isolated Photon–Jet Analysis

A four-stage analysis pipeline for isolated-photon + jet measurements in sPHENIX.

The design separates reconstruction, photon identification, histogram production, and final statistical analysis so that each stage has one clear responsibility and expensive upstream processing is reused whenever possible.

```text
sPHENIX PhotonClusterBuilder
            ↓
           DST
            ↓
     TreeProduction
            ↓
 model-independent trees
            ↓
         PhotonID
            ↓
 candidate-keyed model scores
            ↓
      TreeToHists
            ↓
 sufficient-statistics packages
            ↓
      FinalAnalysis
            ↓
 numerical results
            ↓
          plots
```

## Repository structure

| Directory | Responsibility | Primary configuration |
|---|---|---|
| [TreeProduction](TreeProduction/README.md) | Reconstruction steering and production of model-independent analysis trees | `TreeProduction/config/tree_production.yaml` |
| [PhotonID](PhotonID/README.md) | Photon-ID feature definitions, model training/import, and score augmentation | `PhotonID/model_registry.yaml` |
| [TreeToHists](TreeToHists/README.md) | Measurement selections, sample weighting, response construction, and additive sufficient statistics | `config/measurement.yaml`, `config/samples.yaml` |
| [FinalAnalysis](FinalAnalysis/README.md) | Purity correction, background subtraction, unfolding, normalization, uncertainties, and final numerical results | `config/measurement.yaml` |
| [validation](validation/README.md) | Schema migration, historical-content accounting, and canonical-content comparison | `validation/field_mapping.csv` |
| `coresoftware/offline/packages/CaloReco/` | PhotonClusterBuilder source used to define the corresponding reusable sPHENIX reconstruction dependency | `PhotonClusterBuilder.h/.cc` |

## Design principles

### TreeProduction stores reconstruction facts

`TreeProduction` performs the expensive DST-level work once and writes a reusable, model-independent representation.

The producer is one `PhotonJetTree` `SubsysReco` serving:

- p+p and Au+Au;
- data and simulation;
- photon+jet and inclusive-jet simulation;
- supported reconstructed and archived input modes.

The main event flow is visible in `PhotonJetTree.cc`, while implementation is grouped by physical responsibility under `TreeProduction/src/internal/`:

```text
Types.h
    canonical typed vocabulary and persisted records

Production.h / Production.cc
    configuration, conditions, inputs, and reconstruction registration

Event.cc
    source/event identity, trigger information, vertex, MBD,
    centrality, calorimeter state, event accounting, and weights

Photons.cc
    photon candidates, timing, shower information, cells, and isolation

Jets.cc
    reconstructed jets, raw/JES-corrected identities, area,
    ordering, and photon–jet relationships

Truth.cc
    truth vertices, photons, jets, truth isolation,
    and reconstruction-to-truth associations

Output.cc
    ROOT serialization, provenance, and completion metadata
```

Photon candidates are retained independently of any analysis BDT or working point. The tree stores the reconstruction primitives needed to evaluate future photon-ID models without returning to DST production.

### PhotonID owns model-dependent information

Photon identification is downstream of tree production:

```text
canonical trees
      ↓
  features.py
      ↓
train.py / imported model
      ↓
model_registry.yaml
      ↓
  augment.py
      ↓
PhotonScores sidecar
```

Training and inference share one feature-definition implementation.

Model outputs are joined to photons through stable source/event/candidate identities rather than row position.

The base tree is never rewritten when a model changes.

The model registry supports distinct reference and canonical classifiers, including:

- PPG12 reference models;
- canonical p+p classifiers;
- canonical Au+Au classifiers;
- future analysis or detector-development variants.

A photon-ID score is a property of a candidate and a model.

Tight/non-tight classification and ABCD-region assignment remain measurement choices and are therefore applied downstream.

### TreeToHists owns the measurement selection

`TreeToHists` combines:

```text
canonical trees
+
PhotonID score products
+
measurement configuration
+
sample configuration
```

and produces one analysis package per logical sample family.

This stage owns:

- photon-ID working points;
- tight/non-tight definitions;
- reconstructed isolation selections;
- ABCD classification;
- photon and recoil selections;
- sample weighting and stitching;
- centrality weighting;
- photon response;
- joint photon-pT × xJγ response;
- fake, miss, combinatoric, and boundary populations;
- additive normalization denominators;
- sum of weights and sum of squared weights;
- retained statistical support required by the final analysis.

The package boundary is designed around sufficient statistics rather than plot-specific histograms so that downstream analysis choices can change without rereading event trees.

### FinalAnalysis owns the final physics result

`FinalAnalysis` consumes histogram packages only.

Its analysis flow is:

```text
compatible packages
      ↓
purity / leakage correction
      ↓
background subtraction
      ↓
response projection
      ↓
iterative-Bayes unfolding
      ↓
photon normalization
      ↓
statistical propagation
      ↓
persisted numerical result
      ↓
plotting
```

The implementation is divided into:

```text
corrections.py
    leakage-aware ABCD and recoil-background corrections

unfolding.py
    response manipulation, feed-in/feed-out treatment,
    and iterative-Bayes unfolding

statistics.py
    statistical summaries, covariance, diagnostics,
    and iteration-quality criteria

run.py
    complete numerical analysis chain

plot.py
    rendering of persisted results only
```

Plotting is intentionally separated from numerical analysis so presentation changes never alter the physics result.

## Configuration ownership

The repository uses separate configuration files for separate kinds of decisions.

### `TreeProduction/config/tree_production.yaml`

Defines reconstruction and storage:

- collision/data mode;
- input nodes;
- calorimeter reconstruction;
- photon capture;
- jet collections;
- Au+Au background subtraction;
- JES configuration;
- centrality reconstruction;
- truth inputs;
- output-table policy;
- source and provenance bindings.

### `PhotonID/model_registry.yaml`

Defines photon-ID models:

- model identity and version;
- collision system;
- ordered feature list;
- shower definition;
- application domain;
- model format;
- model file and SHA256;
- feature-definition identity;
- optional CDB binding.

### `config/measurement.yaml`

Defines the physics measurement:

- nominal photon-ID model;
- working points;
- tight/non-tight definitions;
- isolation selection;
- photon-pT bins and response support;
- jet radius and recoil requirements;
- centrality selection;
- response axes;
- purity and background-correction strategy;
- unfolding configuration;
- uncertainty configuration.

### `config/samples.yaml`

Defines sample relationships and normalization:

- sample identity and role;
- generator slices;
- ownership/stitching regions;
- cross sections;
- generated-event counts;
- SI/DI composition;
- period information;
- centrality reweighting;
- normalization inputs.

## Reuse boundaries

The pipeline is structured so that changes restart from the earliest stage whose information actually changes.

```text
plot style change
    → plot.py only

unfolding/statistical-method change
    → FinalAnalysis

working-point or analysis-selection change
    → TreeToHists

new photon-ID model
    → PhotonID → TreeToHists → FinalAnalysis

new histogram/response support
    → TreeToHists → FinalAnalysis

new reconstruction primitive
    → TreeProduction onward
```

This separation is the central design constraint of the repository.

## PhotonClusterBuilder

The repository includes the corresponding `PhotonClusterBuilder` implementation under:

```text
coresoftware/offline/packages/CaloReco/
```

Its role is reusable detector-level photon reconstruction:

- photon-cluster construction;
- shower-shape quantities;
- configurable isolation;
- p+p and Au+Au support;
- optional generic model evaluation;
- optional local-file or conditions-database model binding.

Analysis-tree retention does not depend on an optional photon-ID model.

The analysis-specific training and model-selection authority remains in `PhotonID`.

## Validation and migration

`validation/` is separate from the normal analysis path.

It exists to establish and preserve the relationship between historical production products and the canonical representation used by this repository.

```text
historical product
        ↓
field mapping / migration
        ↓
canonical representation
        ↕
scientific-content comparison
        ↑
direct canonical production
```

Its responsibilities include:

- field-by-field historical schema accounting;
- classification of information as stored, derived, downstream, or archival;
- migration into the canonical representation;
- identity-based content comparison;
- validation of values, validity states, weights, associations, and provenance.

ROOT-file byte identity is distinct from scientific-content equivalence. Files may have different serialization while representing the same canonical scientific content.

## Development dependencies

Python stages use:

- Python 3.10+
- NumPy
- uproot
- PyYAML

Photon-ID training additionally uses:

- XGBoost
- scikit-learn
- PyROOT / TMVA

Plotting uses:

- matplotlib

Tree production uses the sPHENIX software environment together with ROOT, yaml-cpp, FastJet, and the corresponding detector/reconstruction packages.
