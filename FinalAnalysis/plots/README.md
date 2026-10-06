# Figures and tables

ROOT macros, compiled by `final_analysis.sh` into `build/`; they draw saved quantities and compute nothing beyond
projections, ratios and normalisations for display.

| File | Draws |
|---|---|
| `IANFigures.C` | The entry point of `draw`, the shared figure helpers, and the measurement figures of `results.root` (purity and combinatoric corrections, responses, iteration stability, closures, prior, results, centrality overlay with systematic boxes, ATLAS overlays) |
| `QAFigures.C` | The photon and event QA figures of `qa.root` (included by `IANFigures.C`) |
| `IANTables.C` | The IAN result tables of `tables` (LaTeX and CSV) |
| `sPhenixStyle.C`, `sPhenixStyle.h` | The sPHENIX plot style |
| `inputs/atlas_reference.json` | The ATLAS p+p and 0-10% Pb+Pb xJγ of PLB 789 (2019) 167 (HEPData ins1694678 Table 1) for the overlays |

A figure whose inputs are absent or empty is not drawn and is listed in `FIGURES.txt` with the reason; nothing is
invented. The classes, samples and triggers are the ones in the files, so a new class or binning needs no change here.
