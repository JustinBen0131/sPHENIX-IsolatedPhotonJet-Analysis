# Au+Au MBD–calorimeter core cut

The canonical `mbd_calo_cubic_20261020` definition keeps the main correlation band independently for DATA, photon+jet SIM and inclusive-jet SIM. It uses the frozen [calibration](../config/mbd_calo_cubic_20261020.json) with [MbdCaloCoreCut.h](../src/MbdCaloCoreCut.h). Loading this calibration reproduces the black boundaries of the accepted before/after figure; a fresh fit creates a different calibration.

## Apply the frozen cut

For MBD charge sum **Q** and total calorimeter energy **E** in GeV, let `x = Q/1000`. Both the mean and core width are independent cubics:

```text
mu(Q)    = a0 + a1*x + a2*x*x + a3*x*x*x
sigma(Q) = b0 + b1*x + b2*x*x + b3*x*x*x
keep: 0 <= Q < 2120, sigma(Q) > 0, abs(E - mu(Q)) <= 3.5*sigma(Q)
```

ROOT and the environment's `nlohmann/json.hpp` are required. From the repository root, load the appropriate family once, then apply it after the baseline event selection:

```cpp
#include "TreeToHists/src/MbdCaloCoreCut.h"

const auto cut = mbdcalo::loadCalibration(
    "TreeToHists/config/mbd_calo_cubic_20261020.json", "data");
// Use "phosim" for photon+jet SIM, "jetsim" for inclusive-jet SIM.
// Inside the event loop, after the baseline selection:
if (!cut.pass(mbdChargeSum, totalCaloEnergyGeV)) continue;
```

Use the stored MBD charge sum and total calorimeter energy used in the diagnostic. This is total energy, not transverse energy. Boundaries are inclusive; non-finite values and charge outside the domain fail. The display stops at Q=2000, but that is not an event cut.

Apply the minimum-bias classifier and **|z_vtx| < 10 cm** in all three samples, with unchanged stored **0–80% centrality**, existing sample ownership and event-quality requirements. DATA additionally requires the scaled trigger **Photon 12 GeV + MBD NS >= 2, vtx < 150 cm**. The calibration was derived with unit event counts, without physics, centrality or vertex weights. Those baseline gates are the caller's responsibility, not implicit in `pass()`.

## Derivation

`mbdcalo::derive(histogram, sampleName)` takes one sample's unweighted pre-cut TH2 with Q on x and total calorimeter energy in GeV on y, each with 150 bins over [0,3000]. Run it separately for each sample.

1. In each native 20-unit charge slice with at least 1000 entries and peak bin count 100, find the main mode and seed the width from FWHM/2.354820045, with a 10 GeV minimum.
2. Fit a Gaussian to the core within approximately 1.3 seed widths, expanding to at least five energy bins. ROOT `TH1::Fit("LISQRN0")` uses bin-integrated Poisson likelihood. Amplitude, mean and width are constrained around the mode seed.
3. Fit the slice means and widths versus Q/1000 with separate `TGraphErrors` cubics, using their fit uncertainties. `SQNR0` obtains the weighted linear least-squares seed; `SQNR0F` refits all four free coefficients with Minuit2/Migrad, tolerance 1e-6 and at most 100000 function calls. Require status 0, a valid fit and covariance status 3.
4. Freeze the resulting coefficients before event processing. The supplied calibration rounds the accepted coefficients to eight decimal places, exactly as used for its plotted selection.

Here sigma is the fitted **Gaussian core width**, not the full-population RMS. It differs from the earlier population-mean/RMS prescription.

## Figure correspondence and scope

| Sample | Rejected, binned estimate | Fitted charge support |
|---|---:|---:|
| DATA | 1.51975749% | 50–2110 |
| Photon+jet SIM | 0.59273226% | 50–1970 |
| Inclusive-jet SIM | 0.54679517% | 50–1970 |

The before/after figure applies the cut at centres of 20-charge by 20-GeV bins. Event-by-event application uses the continuous equation, so exact rejection fractions can differ. Flow bins are not included in these percentages. The full application domain is [0,2120); intervals outside the fitted support are extrapolation. Cubic fits converge but do not have acceptable statistical goodness of fit everywhere, especially near the high-charge bend. These limits remain part of the definition's validation scope.

This release supplies the canonical definition and application interface. Existing production histograms, training products and trees are not retroactively changed. Downstream processing must explicitly load this version; a different vertex selection or a new calibration requires a new derivation and validation.
