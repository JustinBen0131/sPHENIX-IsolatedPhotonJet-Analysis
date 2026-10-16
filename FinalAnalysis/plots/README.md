# Figures and tables

ROOT macros read saved analysis products and use the shared sPHENIX style.
They are built and run by `final_analysis.sh`.

| File | Responsibility |
|---|---|
| [IANFigures.C](IANFigures.C) | Measurement figures, the figure layout QA shares, and the `draw` entry point |
| [QAFigures.C](QAFigures.C) | Photon, event, simulation-stitching and trigger QA figures, and the figures of the saved photon-ID models |
| [IANTables.C](IANTables.C) | LaTeX and CSV result tables |
| [sPhenixStyle.C](sPhenixStyle.C), [sPhenixStyle.h](sPhenixStyle.h) | Shared plot style |
| [atlas_reference.json](inputs/atlas_reference.json) | Reference data for ATLAS overlays |

See [../fig/README.md](../fig/README.md) for output names. The comment at the head of each macro lists which function
draws which figure of the note, and the sections of the file in order.

A figure whose inputs are absent is listed in `FIGURES.txt` with the reason. A saved-model figure whose optional table
is absent is left out. The jet-radius figure reads the result files of the other radii that stand beside the results
file, `radius_r02.root` and `radius_r03.root`, each with its own receipt, and is drawn when both are there.

Analysis and QA save the quantities these macros draw. What a figure works out itself is for display: projections,
rebinning, normalisations and ratios. The few functions that derive more (a display fit of a stitched spectrum, a
trigger turn-on point, the area under a score curve, the level of a simulated component) say "computed here" in
their comment.

Candidate legends distinguish truth class from sample origin: **Background-tagged**
means truth-tagged background; **Inclusive+jet** means all selected candidates
from that sample without a truth-class restriction. The renderer checks this
population binding before saving and records it in the PNG receipt.

## Layout

Two-dimensional distributions use the shared ROOT `kBird` palette. Renderers inherit it from `sPhenixStyle.C`; avoid local palette overrides. To redraw only the four before/after NCB eta-time maps from an existing QA product, use `final_analysis.sh draw --qa QA.root --fig DIRECTORY --only ncb_eta_time`. This changes no histogram inputs.

The figures write their text inside the frames, in the space the data leave. Only the maps, whose cells fill the
frame, carry the sPHENIX mark and their header above it: the two response maps and the eta-phi maps of the clusters
and the jets.

- **One frame** (`ianFrame`): the sPHENIX mark and the header lines at the top left, the legend at the top right, or
  under the header when the two do not fit side by side, and a y range opened until no drawn point reaches either.
- **Results** (`resultFrame`): the same header, the legend right-aligned, a y axis from zero.
- **Frames on shared axes** (`sharedGrid`, `gridFrame`, `gridLines`, `gridLegend` in QAFigures.C): frames a few
  pixels apart, tick labels on the outer edges, the mark, the lines and the legend inside the frames, each y range
  opened for the text above the points (`gridCeiling`).
- **Event distributions** (`eventDistribution`): for a falling spectrum the header and the legend form one block at
  the top right.
- **Threshold grids** (`workingPoints`): shared axes, the centrality range over each small panel, the mark and the
  sample line in the first panel, the legend in the last panel of the first row, the axis titles once on the canvas.
- **Training diagnostics** (`diagnosticText`, `legendRow`): text placed in canvas pixels; each canvas carries a line
  with the words "Training diagnostic", which the saving code requires.
- **Input note**: a figure drawn from simulation alone takes its note from the simulation samples (`simulationNote`),
  not from the DATA coverage.
- **Words**: as the figures of the note write them; p+p and Au+Au.

## One rule each

- **Text sizes** are in pixels of the canvas. A PNG is written at three times the canvas size; while it is written,
  `save()` widens every line and dash by that factor and sets every text to the exact size of its font step
  (`exactTextSizes`), so the image shows a text as wide as it was measured and placed.
- **A symbol with an index late in a line** is set as one fragment (`apart`), so the index stays on its symbol. An
  absolute value keeps a thin space after its opening bar (`absolute`).
- **Dashed lines** have one pattern, set in the entry point and drawn by `dashedLine`.
- **Axis ticks** have a length in pixels on frames of every shape (`tickLengths`).
