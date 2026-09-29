# Analysis figures

Run `../final_analysis.sh render --results RESULTS --out FIGURES` from this directory, or use the repository-root command in the package README. Calculation commands save numerical inputs first; rendering consumes them. The measurement figure package must reproduce the recorded per-photon result exactly.

| File | Responsibility |
|---|---|
| `figure_package.py` | Saved measurement intermediates and consistency checks |
| `native_response.py` | Source-weighted finest-grid response and exact restriction checks |
| `render.py` | Measurement, response, closure and ATLAS-comparison drawings |
| `render_qa.py` | Manifest/inventory routing, saved QA arrays and numerical preparation, then saved-input drawing |
| `tables.py` | LaTeX and CSV tables from saved results |
| `purity_closure.py`, `purity_nonclosure.py` | Purity numerical studies using measurement kernels |

QA preparation saves histogram values, explicit variances and edges, source-bound trigger spectra/run tables, the named-trigger centrality target, and purity, signal-efficiency quantile and descriptive-fit results as ordinary NPZ/JSON. Rendering verifies this package and never resolves models or accesses source ROOT files. Projection, visual normalization and evaluation of a recorded fit describe the saved inputs; they do not choose or alter measurement settings.

`config/figures.yaml` inventories every requested file. Missing populations, normalization bindings or training records are reported in the availability manifest and QA notes. Empty histograms are absent support. Trigger turn-ons retain the numerator/denominator semantics of scaled reference events and photon live bits, with each source's enabled-trigger table.
