# FinalAnalysis: corrections, unfolding, comparison

* `run_corrections.py` takes the TreeToHists payloads for data and for the
  truth-signal-only photon+jet simulation (prompt-photon leakage), solves the
  leakage-aware ABCD purity per photon-pT bin, forms the purity-corrected
  recoil spectrum (region A minus the background transfer times region C,
  with the region-C prompt leakage restored) and unfolds it with the
  iterative-Bayes kernel on the response bundle.  Output is one JSON file with
  every intermediate quantity and the input hashes.
* `plot_overlay.py` draws corrected or unfolded (1/N_gamma) dN/dxJ for several
  results on one figure.
* `photonjet/` holds the numerical kernels and the reference tree reducers:
  `analysis/background.py` (ABCD solve, purity correction, toys),
  `analysis/response_builder.py` (tree to response bundle),
  `analysis/response.py` (bundle container and boundary categories),
  `analysis/unfolding.py` (iterative Bayes, refold chi2, iteration scan),
  `analysis/purity.py` and `analysis/reduce.py` (independent event-leading
  ABCD counts and recoil histograms used by the tests as a cross-check),
  `io/tree_validation.py` (contract validator) and `cli.py`
  (`response build`, `purity`, `histogram`, `trees validate`).

Response bundle:

```bash
python FinalAnalysis/photonjet/cli.py response build --system auau --dimension 2D \
    --input scored/auau_photonjet_sim.root --output-stem response/auau
```

Not applied in this first pass, but supported by the kernels: the combinatoric
(unmatched-recoil) subtraction, the photon reconstruction-efficiency correction
of the normalisation, and the response-statistics contribution to the
unfolded uncertainty.
