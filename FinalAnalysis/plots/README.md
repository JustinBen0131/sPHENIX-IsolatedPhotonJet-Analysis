# Figures and tables

ROOT macros read saved analysis products and use the shared sPHENIX style.
They are built and run by `final_analysis.sh`.

| File | Responsibility |
|---|---|
| [IANFigures.C](IANFigures.C) | Measurement figures, the figure layout QA shares, and the `draw` entry point |
| [QAFigures.C](QAFigures.C) | Photon and event QA figures |
| [IANTables.C](IANTables.C) | LaTeX and CSV result tables |
| [sPhenixStyle.C](sPhenixStyle.C), [sPhenixStyle.h](sPhenixStyle.h) | Shared plot style |
| [atlas_reference.json](inputs/atlas_reference.json) | Reference data for ATLAS overlays |

See [../fig/README.md](../fig/README.md) for output names.
Unavailable figures are listed in `FIGURES.txt` with their reason.
Calculations are saved by Analysis and QA before these macros draw them.
