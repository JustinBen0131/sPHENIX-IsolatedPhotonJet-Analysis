# The analysis, and why each choice is made

This document is for a physicist who is opening the repository for the first
time.  It explains what is measured, where the method comes from, what every
selection means, and which parts are not yet settled.  Numbers that are still
pending are marked as such rather than filled in.

## 1. Lineage and target

The measurement is the photon-jet transverse momentum balance
`xJgamma = pT(jet) / ET(photon)` for isolated photons in p+p and Au+Au
collisions at 200 GeV.  It follows the ATLAS measurement in 5.02 TeV Pb+Pb
and p+p (Phys. Lett. B 789 (2019) 167, arXiv:1809.07280) as closely as the
sPHENIX detector and kinematics allow: the same observable, the same
inclusive photon-jet pairing with a back-to-back requirement, purity from a
double-sideband ABCD method, and a two-dimensional Bayesian unfolding in
(photon pT, xJgamma).

Within sPHENIX the method rests on three pieces of work:

* the p+p isolated prompt-photon spectrum analysis (internal note PPG12 and
  its paper draft), which supplies the calorimeter clustering, the
  shower-shape variables, the preselection and non-collision-background cut,
  the eleven-input photon-identification BDT, its tight and non-tight windows,
  the isolation working point, the truth photon classes, the ABCD purity
  method, the simulation cross sections and the unfolding precedent;
* the p+p photon-jet studies that fix the observable, the 15 to 35 GeV
  photon window and the recoil definition used here;
* the Au+Au photon-jet program (internal note PPG19, in preparation), of
  which this repository is the analysis code.

The Au+Au extension adds the underlying-event subtraction, a fourteen-input
photon-identification model with centrality as an input, centrality-dependent
working points, and the simulation weighting needed for embedded samples.

## 2. The observable

For every selected photon, every jet passing the jet selection and
`delta phi(photon, jet) > 7 pi / 8` is paired with it (inclusive pairing, not
only the leading jet).  The reported quantity is the per-photon yield

    (1 / N_gamma) dN / dxJgamma

at particle level, in three photon-pT bins (15-20, 20-25, 25-35 GeV) and
summed over them, with xJgamma in 0.1-wide bins from 0 to 2 and one bin from
2 to 3, reported for xJgamma between 0.3 and 1.8.  In Au+Au the nominal
centrality class is 0-20 percent; the trees carry the continuous centrality so
other classes can be formed.

The ATLAS photons are far harder (63 to 200 GeV); the reference points shipped
in `FinalAnalysis/reference/` allow a shape comparison only.

## 3. Data and simulation

Data: Run 24 p+p at 200 GeV and the Run 24 Au+Au dataset, read from
calibrated DSTs with calorimeter and jet nodes.  The p+p photon trigger is the
EMCal 8x8 patch trigger in coincidence with a minimum-bias requirement
(scaled trigger bit 30 in the trigger vector); the Au+Au photon trigger is the
10 GeV patch trigger with the MBD north-south coincidence and vertex
requirement (scaled bit 22), with the offline minimum-bias classifier
required in addition.  The integrated luminosities used for this analysis
are recorded in the production receipts, not in this repository.

Simulation: PYTHIA-8 photon+jet samples (p+p: 5, 10 and 20 GeV photon
slices; Au+Au: 12 and 20 GeV slices embedded in HIJING) and PYTHIA-8
inclusive-jet samples (p+p: 8, 12, 20, 30, 40 GeV jet slices; Au+Au: 12, 20,
30, 40 GeV embedded).  Embedded means the PYTHIA event is overlaid on a
simulated HIJING Au+Au event before reconstruction, so the simulated
underlying event, centrality, and its effect on isolation and jets are
included.

Simulation samples must be combined with the complete per-event weight of
`config/samples.yaml` (section 8).

## 4. Event selection

* Accepted events only (`terminal_status == 0`): the producer marks events it
  could not fully reconstruct, and those are excluded from the physics
  selection while retained in the trees for bookkeeping.
* Vertex: |z| < 30 cm in p+p and |z| < 10 cm in Au+Au for the nominal
  histograms.  The trees are produced with a wider window (60 cm) so the cut
  can be varied offline.
* Au+Au centrality from the MBD charge through the collaboration centrality
  calibration; the nominal class is 0 to 20 percent.  A new calibration is
  expected; the trees carry the inputs needed to recompute the percentile.

## 5. Photon reconstruction and identification

EMCal clusters are formed with the template cluster builder from calibrated
towers above 70 MeV that pass the tower-status check, with the position
correction applied.  Photon candidates are clusters with ET above 5 GeV and
|eta| < 0.7; the nominal analysis uses 15 <= ET < 35 GeV, where the
identification model is trained, and keeps 5 to 40 GeV in the trees for the
unfolding response and for control.

Shower-shape variables are computed on the 7x7 tower grid around the
cluster's maximum tower, with a 70 MeV floor on the cells:

| Input | Meaning |
| --- | --- |
| `cluster_Et`, `cluster_Eta`, `vertexz` | candidate ET, pseudorapidity, event vertex |
| `cluster_weta_cogx`, `cluster_wphi_cogx` | energy-weighted second moments in eta and phi about the centre of gravity, over the towers owned by the cluster, excluding the seed tower |
| `cluster_weta33_cogx`, `cluster_wphi33_cogx` | the same moments restricted to the 3x3 block (Au+Au only) |
| `e11_over_e33` | seed tower over the 3x3 sum |
| `e32_over_e35` | the 3x2 block on the centre-of-gravity side over the 3x5 block |
| `cluster_et1` to `cluster_et4` | fractions built from the four 2x2 blocks around the centre of gravity: total 2x2 fraction, eta asymmetry, phi asymmetry, corner fraction |
| `centrality` | Au+Au only, in percent |

The photon-identification score is a gradient-boosted decision tree
(XGBoost, exported to TMVA for use in ROOT) on these inputs in this order,
trained on the nominal truth-signal photons of the photon+jet simulation
against candidates of the inclusive-jet simulation, with the recipe in
`config/nominal.yaml` and `PhotonID/README.md`.  Candidates are scored once,
by `PhotonID/score_trees.py`, with the same model for data and simulation.

Working points are thresholds on that score.  The tight threshold is set to
a fixed signal efficiency (80 percent in Au+Au, per photon-pT and centrality
bin); the bounded non-tight band is a lower score band used as the ABCD
sideband, chosen from simulation closure rather than tuned on data.  Both are
stored in the configuration and are currently pending derivation for the
new trees.

## 6. Isolation

The isolation energy is the sum of calorimeter tower ET (EMCal, inner and
outer HCal) in a cone of radius R = 0.3 around the candidate, minus the
candidate's own energy.  In Au+Au the cone is evaluated on the
underlying-event-subtracted towers, so the same threshold has the same
meaning across centrality; the shower shapes are evaluated on the
unsubtracted EMCal.  R = 0.3 is the nominal cone for both systems; R = 0.4 is
stored as a control.

The isolated and non-isolated thresholds are pending derivation for the new
trees.  In the p+p reference analysis the isolation working point is a
linear function of ET at 80 percent efficiency, and the non-isolated
sideband starts 0.8 GeV above it; in Au+Au the thresholds are centrality
dependent.  A candidate between the two thresholds is in neither region.

## 7. Truth photon and the training label

In simulation a truth photon is a signal photon when it is a Geant photon
with a valid generator association, its photon class is -1, 0, 1 or 2, its
truth isolation (sum of primary particles in R = 0.3 excluding the photon)
is below 4 GeV, and |eta| < 0.7.  The classes follow the p+p reference
analysis: 1 direct (2 to 2 hard process), 2 fragmentation, 3 decay of a
hadron, 0 unclassifiable, -1 no production vertex found.  Decay photons are
rejected; the others are accepted because the analysis measures isolated
prompt photons irrespective of how the generator labels the hard process.

This definition is used everywhere a truth photon is needed: as the training
signal, for the prompt-photon leakage into the ABCD sidebands, for the
combinatoric template, and for the responses.  Reconstructed candidates in
the signal simulation that are not linked to such a photon are dropped from
training, not relabelled as background.

## 8. Weights

One complete weight per simulated event, in this order:

1. the producer weight stored in the tree (`event_weight`; for p+p simulation
   it already contains the generator-slice cross section relative to the
   family reference, the truth-vertex weight, the interaction-mix weight and
   the period-luminosity weight; otherwise it is 1);
2. the source factor: for p+p simulation the family reference cross section
   over generated events; for embedded photon+jet samples the sample cross
   section over generated events; for embedded inclusive-jet slices the
   ownership stitching, in which an event belongs to the slice whose
   half-open window contains its leading R = 0.4 truth-jet pT (12-21, 21-31,
   31-41, above 41 GeV) with weight sigma_effective over owned events, and
   events outside their own slice's window are dropped;
3. the Au+Au centrality factor, the ratio of the data to the simulated
   centrality distribution in 5 percent bins over 0 to 80 percent, one map
   for photon+jet and one for inclusive-jet simulation.

Sumw2 is accumulated from the square of this complete weight.  Nothing is
applied twice: the histogram stage and the response builds call the same
function.

## 9. Purity and background corrections

Photon candidates are classified into four regions by the two independent
requirements, identification (tight or non-tight) and isolation (isolated or
non-isolated): A tight and isolated, B tight and non-isolated, C non-tight
and isolated, D non-tight and non-isolated.  If identification and isolation
factorise for the background, the signal count in A is
`N_A - N_B N_C / N_D`.  Prompt photons leak into B, C and D; the leakage
fractions are taken from the photon+jet simulation with the nominal truth
definition and enter a leakage-aware solution of the same equation.  The
region-C candidates, with their prompt leakage restored, provide the shape
of the fake-photon recoil that is subtracted from the region-A recoil
spectrum.

In Au+Au a second background is subtracted: recoil jets in events with a
matched prompt photon whose jet has no truth counterpart (combinatoric
recoil from the underlying event).  Its shape and normalisation per prompt
photon come from the embedded simulation and it is scaled to the data photon
count.

## 10. Response and unfolding

The detector response is built from the photon+jet simulation as a matrix
from truth (photon pT, xJgamma) to reconstructed (photon pT, xJgamma), on
grids extended beyond the reported window (truth 5 to 40 GeV, reconstructed
10 to 40 GeV) so that migration across the window edges is explicit.  Truth
pairs without a reconstructed match are misses; reconstructed pairs whose
truth lies outside the window are boundary fakes; reconstructed pairs with a
matched photon and jet that are not the selected truth pair are detector
fakes.  Fake photons and combinatoric recoil are not part of the fake cause
of the unfolding because the purity correction and the combinatoric
subtraction remove them from the data first.

The corrected spectrum is unfolded with the iterative Bayesian method
(D'Agostini), the truth spectrum as prior.  The photon count is unfolded
separately with a per-event photon response (truth photon pT against the
reconstructed ET of the event-leading tight isolated photon).  The number of
iterations is chosen by scanning 2 to 12 and applying the gate and score of
the maintained analysis: finite non-negative result, physical purity,
factorising ABCD, no zero-efficiency bins, at least three reported bins,
refolding and closure chi2 per degree of freedom below 5, at least 90 percent
of toys succeeding, bounded relative uncertainties; among the candidates that
pass, the lowest score wins.

Statistical uncertainties are propagated by repeating the whole correction
and unfolding on Gaussian fluctuations of the data sufficient statistics
(and of the combinatoric template); the reported error is the central 68
percent half-width and the covariance is the rescaled toy covariance.
Response, leakage and template statistics are held fixed, as in the
maintained analysis.  Systematic uncertainties are not evaluated in this
repository.

## 11. Known caveats

* The producer and tree builder in `TreeProduction/` are a frozen snapshot
  from before the current Au+Au production; they are refreshed verbatim from
  the frozen production source when it completes.
* No model, tight threshold, non-tight band or isolation threshold is
  derived for the new trees yet; the stages stop rather than substitute
  values.
* The Au+Au centrality calibration and the EMCal zero-suppression
  cross-calibration are being updated by the collaboration; the centrality
  maps in `config/samples.yaml` were derived for the previous production and
  must be re-derived.
* The number of generated events of the p+p simulation samples is not
  recorded in `config/samples.yaml`; until it is, p+p simulation is combined
  with relative normalisation only.
* The Au+Au integrated luminosity has two published values in the internal
  record that differ by a factor of about three; this repository does not
  compute cross sections, so it does not depend on the choice, but any
  absolute normalisation must state which value it used.

## 12. Names a reader will meet in the trees

* `ReplayFoundationV1` is the directory in every producer output file that
  holds the record trees (`RJEventV1`, `RJPhotonCandidateV1`, `RJJetV1`,
  `RJTruthPhotonV1`, `RJTruthJetV1`, `RJRecoTruthLinkV1`, ...) from which the
  collaborator trees are built.
* `bdt_input_00` ... `bdt_input_13` are the ordered model inputs of section 5
  as stored per candidate; `bdt_input_count` says how many are filled.
* `source_file_index`, `event_id_hi/lo`, `candidate_id_hi/lo`, `jet_id_hi/lo`
  are stable identities; every join in the code uses them, never entry order.
* Identifiers in the producer source that contain a number after `the` are
  historical study labels with no physics meaning.
