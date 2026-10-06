# TreeToHists

Fill histogram packages from eventTree v1 files using one shared selection
for spectra, QA and simulation responses. Inputs are read only.
Requires ROOT 6, its C++ compiler and a POSIX shell.

## Run locally

From the repository root:

```bash
./TreeToHists/tree_to_hists.sh local --settings TreeToHists/config/treetohists.cfg \
    --inputs catalogue.json --work work --out histograms --parallel 8
```

`local` builds the driver, plans one job per input file, fills and merges.
A completed job is reused only when its input and plan match. Use a new work
directory for changed settings. Existing output packages are not replaced.

## Inputs and settings

The input catalogue names produced ROOT files or published tree groups.
TreeProduction's DST `.list` files are used earlier, when making the trees.

```json
{
  "schema": "PhotonJetInputCollectionV1",
  "samples": [
    {"system": "PP", "name": "data", "family": "data",
     "groups": ["trees/PP/run51109"]}
  ]
}
```

Use `files` instead of `groups` for direct ROOT-file paths. Paths are relative
to the catalogue. Simulation families are `photonjet` and `inclusivejet`;
supply their qualified normalization and weighting bindings. AuAu simulation
also requires selection-matched centrality and truth-vertex weights.

[config/treetohists.cfg](config/treetohists.cfg) uses ROOT TEnv `Key: value`
settings. It selects photon ID, isolation, DATA triggers and storage axes.
The defaults use the frozen nominal photon-ID definition and AuAu isolation
line. The DATA triggers are pp Photon4 and AuAu Photon12.
For an optional model, set `auau.PhotonIDSelection` or
`pp.PhotonIDSelection` to its `selection.json`.

The bundled MIT-licensed nlohmann/json headers read catalogues, ROOT metadata
and selection definitions and write package records. No separate JSON-library
installation is needed.

## Run on Condor

```bash
./TreeToHists/tree_to_hists.sh plan --settings TreeToHists/config/treetohists.cfg \
    --inputs catalogue.json --work work
./TreeToHists/tree_to_hists.sh condor --work work --logs logs \
    --events-log /tmp/$USER/photonjet/events.log \
    --setup "/opt/sphenix/core/bin/sphenix_setup.sh -n new"
mkdir -p "/tmp/$USER/photonjet"
condor_submit work/condor/fill.sub
./TreeToHists/tree_to_hists.sh merge --work work --sample PP/data --out histograms
```

`condor` prepares the submit file; the analyst submits it. Merge each sample
after its jobs finish. Jobs use the supplied setup script and prebuilt driver.

## Additional commands

- `centrality` derives AuAu simulation weights from the DATA package.
- `isolation population|derive|freeze|check` builds and checks an isolation definition.
- `Fill.ResponseVariation` in the settings varies simulated response scale/resolution.
- `Fill.RequireStage3Fractions: 1` requires filled jet calorimeter fractions;
  the default 0 permits inputs where they remain unavailable.

Run `tree_to_hists.sh` without arguments for command syntax.
[FinalAnalysis](../FinalAnalysis/README.md) provides the end-to-end `all` command.

## Output

```text
<out>/<system>/data.root       DATA package
<out>/<system>/<sample>.root   simulation package
<out>/<system>/<sample>.json   package sidecar
```

Packages contain measurement and QA histograms, simulation responses,
event records and source contributions. Histograms retain Sumw2; event/source
records support correlated resampling in FinalAnalysis. Simulation's absolute
sample normalization is recorded here and applied by FinalAnalysis.
Incomplete or unknown input coverage remains labeled.

## Source map

Start at `TreeToHists.cc::main`. `plan` resolves inputs and settings;
`fill` calls the common selection and fill macros; `merge` writes one sample package.

| File | Responsibility |
|---|---|
| [TreeToHists.cc](src/TreeToHists.cc) | Settings, catalogue, planning, filling and merge |
| [PhotonJet.h](src/PhotonJet.h) | Shared cuts, regions, weights and truth relations |
| [FillPhotonJet.C](src/FillPhotonJet.C) | Measurement histograms and event records |
| [FillPhotonJetResponse.C](src/FillPhotonJetResponse.C) | Simulation response, misses and fakes |
| [FillPhotonQA.C](src/FillPhotonQA.C) | QA and retained-observable extensions |
| [Isolation.cc](src/Isolation.cc) | Isolation derivation and checks |
| [TreeToHists.h](src/TreeToHists.h) | Shared driver declarations |

Input compatibility is checked automatically. Full-population statistics and
the normalization, weighting and physics-closure qualifications remain
necessary for final scientific use.
