# Photon selections for comparisons

Start with [config/treetohists.yaml](config/treetohists.yaml). It selects the
photon-ID models, score working points, DATA triggers and AuAu isolation
package. The resolved selection saved with a histogram package records what
that particular package used; use that record when comparing an existing
measurement.

## Source map

| Definition | Source |
|---|---|
| NPB and shower-shape preselection, evaluated from primitives | [src/Preselection.h](src/Preselection.h), `newPPG12` |
| Resolved event, photon-ID and isolation configuration | [src/PhotonJet.h](src/PhotonJet.h), `configuredProfile` |
| Photon kinematics and ABCD region classification | [src/PhotonJet.h](src/PhotonJet.h), `photonKinematics` and `photonClass` |
| Photon ordering and one leading photon per ABCD region | [src/PhotonJet.h](src/PhotonJet.h), `regionLeadersWithClassification` |
| AuAu isolation threshold and subtraction convention | [config/auau_isolation_wp80_15to35.json](config/auau_isolation_wp80_15to35.json) |
| Where measurement histograms are filled | [src/FillPhotonJet.C](src/FillPhotonJet.C) |

## Nominal photon definition

The common preselection requires finite inputs, `pho_npb_score > 0.5`,
`0.6 < pho_et1 < 1`, `pho_e11_over_e33 < 0.98`, and
`0.8 < pho_e32_over_e35 < 1`. The saved `pho_presel_active` flag is a
cross-check. The implementation recomputes this decision from those fields.

The configured nominal models are `ppg12_original` for pp and
`nominal_retrained_v1` for AuAu. The configuration contains the current
AuAu tight and non-tight score lines; their independent variable is
centrality in percent. The non-tight lower bound is 0.1. The pp score lines
depend on photon transverse energy. `configuredProfile` resolves these
settings before `photonClass` assigns the ABCD regions.

The pp nominal isolation uses radius 0.4 with threshold
`0.490 + 0.037 * photon_ET` GeV and a 0.8 GeV gap to the non-isolated band.
The pp simulation isolation transformation is `1.2 * isolation + 0.1` GeV;
see `analysisIsolation` in `PhotonJet.h`.

AuAu uses the **UE-subtracted SUB1** isolation at radius 0.3. Its nominal
threshold is `3.059196113029283 - 0.02356487430593699 * centrality` GeV for
`0 <= centrality < 80`, with a 0.8 GeV gap to the non-isolated band. The JSON
package records the cone index, subtraction method and threshold together.

The nominal vertex windows are `abs(z) < 60 cm` in pp and `abs(z) < 10 cm`
in AuAu. DATA trigger choices are pp Photon4 and AuAu Photon12; the resolved
event-selection record contains the associated minimum-bias requirement.
Use `photonKinematics` and the resolved package record for the photon energy
and pseudorapidity acceptance and their strict boundary conventions.

## Comparing a photon RAA measurement

Compare the event selection, photon acceptance, preselection, model and
working points, isolation subtraction and cone, and ABCD boundaries. Check
the model release identity as well as its numerical score cut. A different
model can require a different score threshold at the same efficiency.

This photon-jet measurement chooses the highest-energy eligible photon in
each ABCD region. An inclusive-photon RAA analysis must also compare its
photon counting and normalization with that choice. The recoil-jet cuts and
photon-jet angular window are additional selections of the pair measurement.

The source includes systematic variations alongside the nominal settings.
Publishing their implementation does not certify completed numerical
systematic coverage; inspect the result's own coverage record for that.
