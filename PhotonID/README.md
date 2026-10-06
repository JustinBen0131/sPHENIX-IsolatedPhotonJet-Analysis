# PhotonID

Train photon-identification models from completed eventTree groups and add
their named scores. Trees already contain the nominal `pho_bdt_score`;
optional models preserve that branch and the nominal isolation definition.

## Setup

Use Python with NumPy, SciPy, PyYAML, uproot, awkward and ROOT. BDT fitting
also needs XGBoost. Set `MYINSTALL` to TreeProduction's installed eventTree
library prefix; a `--disable-producer` build is sufficient.
`PHOTONID_PYTHON` can select your Python executable.

## Run

From the repository root:

```bash
./PhotonID/photon_id.sh --help
./PhotonID/photon_id.sh listModels AuAu
./PhotonID/photon_id.sh reTrainNominal AuAu --check
./PhotonID/photon_id.sh reTrainNominal AuAu
./PhotonID/photon_id.sh trainNew centAndETdepModels AuAu
./PhotonID/photon_id.sh trainNew MLP AuAu
./PhotonID/photon_id.sh trainNew logReg AuAu
./PhotonID/photon_id.sh score MLP_v1 AuAu
```

Training prepares the inputs, fits models, validates the saved artifacts,
derives working points and scores complete groups. `score` applies an existing
release. `--check` previews inputs and dependencies; `--no-score` stops after
the model release. Matching completed work is reused.

[config/training.yaml](config/training.yaml) holds the recipes and standard
group bindings. For custom inputs, repeat `--group /path/to/published/group`.
Use `--models DIR` for a writable release root and `--release-version v2`
for an explicitly new release.

## Models

| Model | Score branch |
|---|---|
| `nominal` | `pho_bdt_score` |
| `nominal_retrained_v1` | `pho_bdt_score_nominal_retrained_v1` |
| `centCoarse_v1`, `centFine_v1` | `pho_bdt_score_cent_coarse_v1`, `pho_bdt_score_cent_fine_v1` |
| `ET_v1` | `pho_bdt_score_et_v1` |
| `ETcentCoarse_v1`, `ETcentFine_v1` | `pho_bdt_score_et_cent_coarse_v1`, `pho_bdt_score_et_cent_fine_v1` |
| `MLP_v1`, `logReg_v1` | `pho_mlp_score_v1`, `pho_logreg_score_v1` |

Optional scores are photon-aligned and NaN outside their model's domain.
The recipe fixes the features, source roles, training weights and shared
event-grouped split. Working points belong to their model; they are not
interchangeable between releases.

## Outputs and use in TreeToHists

Releases live in `models/<system>/<name>/` or the `--models` directory.
They contain the manifest, model artifacts, validation, working points and
`selection.json`. Keep that directory intact.

Select an optional model in a copy of TreeToHists' `treetohists.cfg`:

```text
auau.PhotonIDSelection: /path/to/models/AuAu/ETcentFine_v1/selection.json
```

Use `pp.PhotonIDSelection` for a retrained pp model. Omit score, state and
working-point overrides for that system; the release supplies them.
Scoring creates a checked group version and advances `CURRENT` only when
every source is ready. Partial work leaves the published group unchanged.

## Source map

The driver follows preparation → fitting → validation → working points →
scoring. Training labels and weights are defined in `Inputs.py`; histogram
selection belongs to TreeToHists.

| File | Responsibility |
|---|---|
| [photon_id.py](scripts/photon_id.py) | Commands and phase order |
| [Inputs.py](src/Inputs.py) | Features, labels, weights and event split |
| [Training.py](src/Training.py) | BDT, MLP and logistic fitting |
| [Validation.py](src/Validation.py) | Numerical and saved-model checks |
| [WorkingPoints.py](src/WorkingPoints.py) | Efficiency thresholds |
| [Scoring.py](src/Scoring.py) | Photon-aligned scores and group publication |
| [Models.py](src/Models.py) | Recipes, release names and artifact checks |

Unbound source populations are refused; the standard pp training population
is currently unbound. Local model checks do not establish production or
downstream physics qualification.
