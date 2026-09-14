# FinalAnalysis: responses, corrections, unfolding, comparison

Entry points:

* `photonjet/cli.py response build` builds the (photon pT, xJ) pair response
  from photon+jet simulation trees, `photonjet/cli.py photon-response build`
  the per-event photon response used for the photon denominator.  Both use
  the nominal truth-signal contract and, with `--sample`, the complete
  analysis weights of `config/samples.yaml`.
* `run_corrections.py` runs the chain for one system and writes the result
  JSON and a CSV of final points.
* `plot_overlay.py` draws final spectra and published reference points
  (`reference/`).

The chain (`photonjet/analysis/chain.py`), in order:

1. Leakage-aware ABCD purity per photon-pT bin from the data counts and the
   simulated prompt-photon leakage; purity-corrected recoil spectrum (region A
   minus the background transfer times region C, region-C prompt leakage
   restored).  Au+Au uses one background transfer from the integrated
   15-35 GeV population, p+p one solve per bin (`corrections.purity_strategy`).
2. Combinatoric subtraction (Au+Au): simulated recoil of matched prompt photons
   whose jet has no truth jet, scaled by the data photon count over the
   simulated Region-A prompt-photon count.
3. Joint iterative-Bayes unfolding in (photon pT, xJ) with the truth spectrum
   as prior; the fake cause is detector plus boundary fakes only, since fake
   photons were removed in 1 and combinatoric recoil in 2.
4. Photon-count unfolding with the per-event photon response; fake cause the
   photons whose truth pT lies outside 15-35 GeV.
5. `(1/N_gamma) dN/dxJ` = unfolded pairs summed over pT divided by unfolded
   photons and the bin width; also per photon-pT bin.
6. Toys of the data sufficient statistics (and of the combinatoric template
   when 2 is on) rerun 1-5; errors are the central 68 percent half-widths,
   the covariance is the winsorized toy covariance rescaled to them.
7. Diagnostics: refolding chi2/ndf, photon refolding chi2/ndf, simulation
   closure and refold, zero-efficiency bins, negative-input fraction, toy
   success fraction.

Iteration count: `unfolding.iterations` in the configuration, or, when null,
a scan of 2..12 with the maintained gate (all diagnostics below their limits,
toy success at least 90 percent, relative errors bounded) and score, keeping
the minimum score.

Not included in the statistical uncertainty: response, leakage and template
statistics (held fixed, as in the maintained analysis).  Systematic
uncertainties are not evaluated in this repository.

`photonjet/analysis/` also holds the reference tree reducers `purity.py` and
`reduce.py`, used by the tests as an independent cross-check of TreeToHists,
the tree validator `io/tree_validation.py`, and the kernels `background.py`,
`response.py`, `response_builder.py`, `photon_response.py`, `unfolding.py`.
