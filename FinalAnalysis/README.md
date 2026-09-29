# FinalAnalysis

Calculate the per-photon xJγ distribution and its whole-event DATA and source-resampled response uncertainties from a TreeToHists collection. This package reads saved histograms and contribution tables; it never reads eventTree. Use the repository environment (numpy, scipy, PyYAML, uproot, matplotlib).

```bash
./FinalAnalysis/final_analysis.sh runPP --manifest hists/manifest.json --out results/pp --check
./FinalAnalysis/final_analysis.sh runPP --manifest hists/manifest.json --out results/pp
./FinalAnalysis/final_analysis.sh runAuAu --manifest hists/manifest.json --out results/auau
./FinalAnalysis/final_analysis.sh runAll --manifest hists/manifest.json --out results/all
./FinalAnalysis/final_analysis.sh render --results results/all --out figures/all
./FinalAnalysis/final_analysis.sh runSystematics --config systematics.yaml --system pp --out results/pp_systematics.json
./FinalAnalysis/final_analysis.sh runSystematics --config systematics.yaml --ratio auau_cent0_20 pp --out results/ratio_systematics.json
```

Set `PHOTONJET_PYTHON` to select the Python interpreter. `--config` defaults to `config/final_analysis.yaml`; copy it to supply production bindings. Its null normalizations and centrality weights are unresolved inputs. Preparation resolves them from the collection only when the exact ingredients and source evidence are recorded. Inconsistent settings, sample identity, hash, selection or normalization fail before calculation. QA-only SIM samples also require recorded normalization and a configuration entry; missing entries stay explicit in the QA notes.

`--check` writes nothing. Execution writes an immutable `PLAN.json`, one JSON/CSV result and numerical figure package per configured system, saved QA arrays and purity/quantile/fit inputs, hash receipts, `RESULTS.json`, and atomic `progress.json`. Each sample's pb/event normalization is applied once; responses retain original source identities independently of physical file count. Fine bins are added before nonlinear corrections. DATA toys resample whole events; response toys resample sources. Neutral correction settings retain the reference numerical behavior.

Use `--resume` only with the identical command, configuration, code and histogram inputs. Completed receipts are verified and reused. Incomplete outputs without a receipt are preserved and refused; choose a new output directory. Changed plans or results are never silently overwritten. Concurrent writers are refused. Rendering requires a fresh or empty directory and reads saved results and QA inputs without reopening ROOT inputs or recomputing the measurement, purity solutions, quantiles or fits.

Rendering writes figure/table manifests and `FIGURE_AVAILABILITY.json` for all 85 configured figure slots. Missing DATA/SIM support and frozen BDT training figures remain `missing_input`; no placeholder image is made. File names carrying an iteration count reflect the actual result, so a fixed slot with a different count remains missing instead of being falsely identified as the same figure.

| Owner | Responsibility |
|---|---|
| `scripts/final_analysis.py` | Prepare, execute/reuse, systematics and render commands |
| `src/analysis.py` | Ordered measurement chain and saved intermediates |
| `src/internal/inputs.py` | Collection validation, event records and source contributions |
| `src/internal/` | Purity, unfolding, prior, bootstrap, combinatoric, closure and systematic kernels; numerical QA preparation |
| `plots/` | Numerical figure packages, QA input preparation, drawing and tables |
| `config/` | Measurement settings, sample definitions and complete figure inventory |

Reference and synthetic qualification test software behavior, not new-production physics acceptance. Nominal model/working-point bindings, contemporary DATA/SIM, cross sections and generated counts, trigger-specific centrality evidence, response/closure support, and the pp SI/DI choice must be qualified on the actual production before scientific use. Measurement edges remain a candidate pending response and closure evidence. Local regression material lives in ignored `tests/`; ordinary commands do not depend on it.
