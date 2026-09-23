# TreeToHists

Consumes completed base trees, selected PhotonID score sidecars, one measurement configuration and sample normalization bindings. `build.py` is the single reader/fanout; `histogram_contract.py` owns additive arrays, integrity checks and serialization. Final selections live here, not in TreeProduction or PhotonID.

## Package

One family package is an NPZ array file plus SHA256-bound JSON metadata. It retains ABCD event-leading counts, A/C recoil spectra, photon and event denominators, truth-conditioned leakage, photon response, joint photon-pT × xJ response, fake-photon/combinatoric/detector categories, misses, boundary counts and sumw2. Internal partitions retain sample, period and SI/DI roles before family summation. Values are additive; nonlinear normalization happens later. Absolute luminosity normalization is not implemented.

Fine response bins are summed into downstream measurement bins. Missing intermediate boundaries cause failure; they are never interpolated back into existence. Off-grid reconstructed/truth matches carry boundary witnesses. Zero-candidate events still enter truth-denominator/miss accounting. The current response requires at most one signal truth photon per event and fails explicitly otherwise; that population rule needs approval.

Arrays, shapes, finiteness, nonnegative counts, response partitions, sample sums and receipt hashes are checked. Existing outputs are not overwritten. A successful integrity check is not scientific acceptance.

## Reading and weights

Event batches read only required columns. Objects must follow the producer's contiguous event grouping; orphan/noncontiguous rows fail. Scores may arrive in any row order: a temporary SQLite identity index bounds RAM and verifies exact base/model/feature bindings, duplicates, missing rows and unknown references. Compact event identities span a file; heavy shower/cell tables and Cartesian pairs do not. Recoil kinematics are derived only for selected photons using the producer formula; equivalence closure remains required.

Bootstrap arrays and per-partition fanout have explicit memory limits. A too-large request stops with a sizing error; it does not silently reduce replicas. Temporary score-index disk space must be available.

DATA trigger selection requires an evidence-bound live/scaled bit mask or an explicit inclusive declaration. The historical legacy word is not treated as a raw trigger word.

The event weight is multiplied once by an approved source factor, explicit SI/DI mixture factor and, where required, an approved centrality map. Producer weights, generated/owned counts, cross sections, source-manifest identities and campaign evidence are checked. Missing normalization never becomes unity. Ownership/support exclusions are recorded. Old production numbers have been removed from runtime defaults.

## Statistics and readiness

One source/event-keyed Poisson multiplicity drives all ABCD/recoil fills for that event, independent of input order or chunking. These replicas retain **data ABCD/recoil covariance**. Sumw2 is also stored. They do not retain full joint simulation-response/leakage/template covariance, systematic variations or cross-production physical-event aliases. Do not claim those capabilities.

`config/measurement.yaml` is unresolved, including proposed bins, isolation, centrality, recoil, response support and uncertainty choices. Nominal construction requires approval evidence, model/feature hashes, complete classification support and known truth/capture states. `--diagnostic` explicitly permits missing measurement/working-point bindings, labels the result diagnostic, and cannot bypass missing simulation normalization or corrupt joins. FinalAnalysis rejects diagnostic packages.

The current canonical model domain proposals do not cover the entire proposed response support. Bind the intended application/working-point domain before response production. No ROOT dataset or numerical package has been processed to validate this implementation.
