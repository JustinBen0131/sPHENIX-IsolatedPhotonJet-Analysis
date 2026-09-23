# FinalAnalysis

Consumes histogram packages only. It must never read event trees to repair insufficient packages. Numerical results are persisted before `plot.py` renders them.

| File | Ownership |
|---|---|
| `corrections.py` | Leakage-aware ABCD, Region-C restoration and optional combinatoric subtraction |
| `unfolding.py` | Fine-support aggregation, feed-in/out bookkeeping and iterative-Bayes kernels |
| `statistics.py` | Error/covariance summaries, diagnostics and iteration gate |
| `run.py` | Package compatibility, chain execution, scan and numerical JSON/CSV output |
| `plot.py` | Rendering persisted numbers and explicitly labelled references |

The arithmetic follows the collaboration predecessor; this port has not been numerically validated. Fine bins are added before nonlinear operations. Positive normalization denominators are required. Integrated 15–35 GeV transfer is restricted to that window. Combinatoric template/denominator inputs must bind the same normalized simulation population, and the retained response-photon denominator must agree with leakage A.

An iteration scan chooses only finite-scored candidates that **pass** the gate. If none passes, nominal analysis fails. The selected final-toy result and a manually requested iteration must also pass. No least-bad failed result is published as nominal.

## Uncertainty contract

Data toys use coherent stored event-bootstrap replicas, or an explicitly selected independent-Gaussian approximation. Response matrices and leakage remain fixed. When combinatoric subtraction is enabled, template and denominator fluctuate independently as Gaussians, matching the predecessor's executed behavior. Their mutual correlation and correlations with response/leakage are not retained. This is not a complete simulation-statistical or systematic covariance treatment.

Stored bootstrap replicas are used without recycling them as additional independent draws. The central-68 interval and winsorized covariance retain the predecessor definitions. The exact nominal uncertainty/iteration recipe still requires approval in `config/measurement.yaml`.

Packages must be analysis-ready and share grids, exact measurement identity and frozen model/feature bindings. Missing response support or stale/mismatched inputs fail explicitly. JSON diagnostics use null for nonfinite values; failed nominal gates do not produce final points. Existing outputs are protected.

Plotting cannot change physics and does not conceal invalid reported errors. Changing style rerenders existing numerical files. Changing unfolding can reuse packages only within the support and statistical information they actually retain.
