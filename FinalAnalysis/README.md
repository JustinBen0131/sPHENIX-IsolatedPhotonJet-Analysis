# FinalAnalysis

Measures the per-photon xJγ distribution of p+p and Au+Au collisions from TreeToHists histogram packages, and draws
every analysis figure and table of the IAN in sPHENIX style under the IAN's own names. ROOT C++; needs `root-config`
and `root` on PATH.

```bash
./FinalAnalysis/final_analysis.sh run         --config FinalAnalysis/config/analysis.cfg
./FinalAnalysis/final_analysis.sh qa          --config FinalAnalysis/config/analysis.cfg
./FinalAnalysis/final_analysis.sh systematics --config FinalAnalysis/config/systematics.cfg
./FinalAnalysis/final_analysis.sh draw        --results results.root --qa qa.root --systematics systematics.root --fig FinalAnalysis/fig
./FinalAnalysis/final_analysis.sh tables      --results results.root --systematics systematics.root --out FinalAnalysis/tables
```

`run` compiles the C++ on first use into `build/` and writes one `results.root` holding every quantity the measurement
figures draw. `qa` writes the photon and event QA inputs (`qa.root`) from the same packages. `systematics` combines the
nominal results with results of one change each. `draw` reads only those files and writes `fig/`, which mirrors the
IAN's `fig/` folder (review figures before using them in a manuscript; `fig/README.md` names the IAN figure each file is);
`tables` writes the IAN's result tables (LaTeX and CSV). An existing results or QA file is never replaced.

## One command

Copy `config/pipeline.cfg`, `config/analysis.cfg`, `config/pipeline_systematics.cfg` and the TreeToHists reduction
settings into a run directory. Bind the scored-tree catalogue, required figure roster and a fresh output directory:

```bash
./FinalAnalysis/final_analysis.sh all --config /path/to/run/pipeline.cfg
```

Paths are relative to the copied configuration. The command runs reduction, QA, the nominal analysis, declared
variations, systematic combination, figures and tables in order. A reduction override rereads the same pinned trees;
analysis-only variations reuse the nominal histogram collection. `nominal-1` and `nominal+1` bracket the actual selected
iterations. The example declares the maintained prior, iteration, combinatoric, purity and jet-scale variations;
scientific acceptance still requires their input populations and uncertainty choices to be qualified.

To reuse histogram packages, give `Pipeline.HistogramCollection` instead of the tree catalogue and reduction settings.
Use `none` for both variations and systematics for a nominal-only run. `Pipeline.ExpectedFigures` names an exact list
of PNG paths relative to `fig/`; `none` explicitly leaves figure coverage unqualified. Missing requested figures,
changed inputs, unknown/duplicate settings and failed child commands stop the run. `PIPELINE.json` records commands,
input/configuration/code hashes, output hashes and the first failure. Individual commands remain available.

## Figures and tables

The ROOT macros share `plots/sPhenixStyle.C`: Helvetica, the sPHENIX Internal mark, white backgrounds, borderless
legends, fixed colours and readable axes. Curve annotations use a quieter region of the actual plotted data, and matrix
margins reserve only the space their headers require. Response-cell text changes contrast with the colour palette.

Physics captions come from the saved source selection and executed analysis settings. Dataset identity, energy,
centrality, photon windows, vertex/trigger requirements, isolation radius, recoil cuts and iteration choices are tied
to those records. QA projection bounds are the same variables that generate their captions. Full response grids retain
the producer's photon window, while measurement figures use the analysis window. No caption infers a generator or
embedding from a file name.

`results.root` and `qa.root` have hash-bound render metadata and companion `.render.json` receipts. The renderer verifies
both before drawing; systematic bands must bind the same nominal result. Every image has a `.json` receipt containing
its hash, source selections, display projections, visible text and renderer/style identities. Tables use saved axes,
selections and uncertainty payloads and carry equivalent receipts. Legacy inputs without these bindings remain
explicitly marked development references, with `LEGACY_UNVERIFIED_ANNOTATION` status. Partial/unknown input coverage
is visible. A software PASS does not certify final statistics, physics choices or every possible future layout;
representative visual QA and the required scientific gates remain necessary.

These packages need ROOT 6 with its C++ compiler, a POSIX shell and `sha256sum` (or `shasum` on macOS).
The MIT-licensed nlohmann/json headers are bundled under `src/vendor/`; no Python or private workspace is needed.

Clean-source builds and the public entrypoints were exercised on ROOT 6.34 (macOS) and SDCC ROOT 6.32 (ana.568). Real
Stage-2 pp DATA reduction passed; adequate statistics and matching completed AuAu/SIM/DI inputs remain scientific gates.
Detailed development qualification records and tests are kept outside the public distribution.

## The chain

| File | Step |
|---|---|
| `src/FinalAnalysis.h` | The settings, the data types, and one function per step |
| `src/Inputs.cc` | Read and check the settings, the measured histograms or event records and the simulated response (whole, or per source); restrict the response to the measurement window |
| `src/Purity.cc` | Leakage-aware ABCD purity with R_bkg, a purity scale and an efficiency scale; combinatoric jets from the simulated template or the data's Δφ sideband |
| `src/Unfolding.cc` | Iterative-Bayes unfolding with a fake cause; per-photon normalisation; the frozen iteration rule |
| `src/Prior.cc` | The joint unfolding's prior: the simulated truth reweighted in xJ to the corrected data |
| `src/Analysis.cc` | The ordered chain per result, the data and response toys, the scans, the closures, `results.root` |
| `src/Systematics.cc` | Systematic covariances and bands from results of one change each; ratios with coherent sources |
| `src/QA.cc` | The QA inputs: DATA and simulated QA sums, ABCD counts, purity, purity closure, the purity study, stitching, working points, trigger study, centrality target |
| `plots/IANFigures.C`, `plots/QAFigures.C` | Maintained measurement and QA figures, sPHENIX style (`plots/sPhenixStyle.C`), IAN file names |
| `plots/IANTables.C` | The IAN result tables |
| `config/analysis.cfg` | Every physics choice of the measurement |
| `config/systematics.cfg` | The systematic sources: which results files of one change each enter, and how |

For each result (p+p and each Au+Au centrality class):

1. **Inputs.** Photon ABCD counts and region A and C recoil jets from the DATA package, either its histograms or its
   event records (`records/events`, `photons`, `jets`: every photon and jet with its event, the events of the class
   selected by their 5% centrality slices); the photon and joint (photon pT, xJ) responses, misses, fakes,
   combinatoric template and signal leakage from the simulated samples, summed with each sample's cross section per
   generated event (`meta/pb_per_generated_event`, recorded by TreeToHists and applied only here). Coarser analysis
   bins are exact unions of the package bins; a boundary the package lacks is an error.
2. **Purity.** In each photon bin, `(A - s) = R (B - l_B s)(C - l_C s)/(D - l_D s)` gives the prompt-photon signal `s`,
   with the simulated signal leakage `l_r` and R_bkg. The region-A recoil jets are corrected with the transfer
   `k = (A - s)/(C - l_C s)`: `(A_pairs - k C_pairs)/(1 - k l_C)`. A photon bin whose equation has no physical root
   stops the result. A purity scale multiplies `s` and re-derives `k`; an efficiency scale divides the bin's corrected
   photons and jets. Both are 1 in the nominal.
3. **Combinatoric jets.** Either measured in the data (`data_dphi`, the nominal): the Δφ distribution of the jets at
   |Δη| > 0.8 is fitted over 0 < Δφ < 2.5 with `a0 + a2 cos 2Δφ + a3 cos 3Δφ`, normalised per xJ bin to all jets in
   the sideband 0.8 < Δφ < 2.5 and integrated over the recoil window Δφ > 7π/8, and removed from the region A and C
   jets before the purity correction; or the simulated template of recoil jets without a truth jet, scaled to the
   corrected photons and subtracted after it.
4. **Prior.** The simulated truth is reweighted by `w(xJ) = exp(t0 + t1 (xJ - 0.8))` (held constant outside the
   reported range) until the folded simulation, scaled to the measured photons, describes the corrected pairs: a
   generalized least-squares fit with the covariance of the data toys. `P(reco | truth)` is unchanged.
5. **Unfolding.** Photon-yield iterative Bayes from the simulated truth, and joint (photon pT, xJ) iterative Bayes from
   the reweighted truth; migration out of the window counts as misses and into it as fakes. The iterations are either
   configured or chosen by the frozen rule (`rule`): with the simulated prior and its own data and response toys, for
   k = 1..12, Q_k = sqrt(S_data^2 + S_response^2 + D_k^2) from the toys' relative half widths and the change from
   k - 1, gated by toy success, refold, closure, negative fraction and conditioning; the smallest passing Q_k of
   k = 2..4 is taken, and none passing stops the analysis with the failing gates.
6. **Result.** Unfolded pairs summed over photon pT, divided by the unfolded photons and the xJ bin width. Data toys
   are corrected (the data-driven combinatoric estimate refitted), refitted and unfolded like the data; response toys
   resample the simulated sources with replacement within each sample and correct, refit and unfold the data again.
   The statistical uncertainty is the data and response central-68% half widths in quadrature; their robust
   covariances (toys beyond five MAD-sigma trimmed) are saved for the systematics. The nominal toy draws one Poisson(1)
   multiplier per event of the event records, so a photon and its recoil jets move together (`event_bootstrap`);
   `poisson` draws every histogram bin independently, for packages without event records.
7. **Closures.** Photon-yield and joint unfolding of the simulation against its truth, plus prior-fit/unfolding closure
   (`Closure.Toys`) after the correction stage. Pseudo-data are response-column sums plus detector and boundary fakes;
   the prior is fitted with the data covariance in simulation units, followed by both unfoldings. Modes use every
   source (`full`) or source-disjoint halves (`A_to_B`, `B_to_A`: within each sample, sources ranked by hash, even ranks A).
   Source resampling gives the errors and pulls (`estimator_closure_<mode>_*`). This test does not exercise ABCD purity
   or data-driven combinatoric subtraction. Their independent combined-chain closure remains a physics acceptance
   gate; a passing unfolding closure cannot substitute for it. Render metadata records this scope explicitly.

**The reported xJ range** starts where the distribution is fully populated: the first xJ edge at or above the recoil-jet
threshold over the lowest photon pT (5 GeV / 17 GeV gives 0.35). The analysis derives it for every result from the
jet slot (`jet5`) and the photon bins, so a new binning or threshold moves it by itself; only the upper edge is a
setting. Every xJ figure draws only that range, on an axis that starts at 0.

## QA

`qa` writes `qa/<package>/` of the QA file (`Output.QA`): every QA histogram of the DATA package (`data/`) and of the
simulation summed with its cross section per generated event (`phosim/`: the results' photon samples, `jetsim/`: the
`QA.<package>.Inclusive` samples), and the numbers the QA figures show, computed there and never in drawing: the DATA
ABCD counts and purity with the prompt-signal leakage, the purity closure on simulation, the stitched spectra with
their descriptive fits (a sample stitched by ownership enters the sum only inside its window, and keeps its whole range
in its own spectrum), the photon-ID and isolation working points per 5% slice, the trigger study of the runs with
every study trigger enabled, the Au+Au centrality target of the package's trigger, and per result the purity study:
R_bkg of the simulated non-signal candidates, the purity non-closure scale, the ABCD transfer terms and the
background factorisation in DATA and simulation. An empty or missing input leaves its figure and numbers out.

## Settings

`config/analysis.cfg` (ROOT TEnv, `Key: value`, `#` starts a comment) holds every choice: input and outputs, the
recoil-jet slot, xJ bins and the reported upper edge, data toys (model and seed), response toys, the closure's toys,
the iteration rule's toys, the prior, the combinatoric method (shape and scope), the QA's inclusive-jet samples, and per
result the simulated samples, photon bins, response grids, iterations (a number or `rule`), and per photon bin R_bkg,
the purity scale and the efficiency scale, then the combinatoric subtraction. `BackgroundCorrelation` and
`PurityScale` may be `qa`: the values the qa step derived from the simulation, read from the QA file (a purity scale
outside 0.5-1.5 is refused). Every key is required and an unknown key is an error, so a misspelt setting can never fall
back to a default. A study is a copy of this file with different values, never new code.

## Systematics

`config/systematics.cfg` names the nominal results file and every source: Up/Down results (a two-sided nuisance),
one Alternative (mirrored) or several Alternatives with a Representative (an envelope), and whether the source is
Correlated between systems. Each results file comes from a copy of `analysis.cfg` with one change (for example
`Prior.Power: 1.5`, `Prior.Form: exp2`, `Combinatoric.Shape: flat`, `BackgroundCorrelation: qa`, `PurityScale: qa`, or
packages from TreeToHists with `Fill.ResponseVariation` for the jet energy scale). Covariances add, bands add in
quadrature, and ratios such as Au+Au/p+p vary correlated sources together.

## Inputs

One histogram package per collision system and sample, as TreeToHists writes them: `<Package>/data.root` with
`measurement/<class>/abcd`, `measurement/<class>/jet<GeV>/recoil`, the event records and `qa/`, and
`<Package>/<sample>.root` with `response/<class>/jet<GeV>/...` (the photon and joint responses, truth, reco, misses,
detector fakes, combinatoric template and `signal_leakage_truth_matched_abcd`), the per-source records and `qa/`.

## Status

On identical inputs the C++ chain reproduces the earlier Python implementation: purity and photons to machine
precision, the fitted prior to the fit tolerance, the per-photon xJ to 3e-5 of its statistical error; the event-record
inputs and the data-driven combinatoric subtraction to 1e-13; every row of the frozen iteration rule to 1e-12 with the
same selected k; the purity non-closure to 1e-8 (float histograms); R_bkg, the transfer terms and the non-closure to
1e-14 on synthetic inputs. The saved development figures in `fig/` come from legacy reference histograms; they are explicitly diagnostic.
Fresh runs write to the configured output directory and never update a manuscript. The QA figures need packages with
QA fills (the new production). The earlier Python implementation (`scripts/`, `src/*.py`, `src/internal/`,
`plots/*.py`, `config/*.yaml`) is no longer used by anything.
