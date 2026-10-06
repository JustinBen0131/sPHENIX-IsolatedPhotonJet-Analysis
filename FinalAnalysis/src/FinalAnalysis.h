/* FinalAnalysis: the per-photon xJgamma measurement from TreeToHists histogram packages. */

/* The chain, one function per step (src/):
   Inputs.cc      the settings, the measured histograms or event records, and the simulated response of one class
   Purity.cc      leakage-aware ABCD purity, background correlation R_bkg, combinatoric subtraction
   Unfolding.cc   the iterative-Bayes unfolding and the per-photon normalisation
   Prior.cc       the data-fitted prior of the joint unfolding
   Analysis.cc    the ordered chain, the data and response toys, the iteration rule and results.root
   Systematics.cc the systematic uncertainties from results of one change each
   Every physics choice is a setting in the configuration file (config/analysis.cfg); nothing below chooses one. The
   reported xJ range starts where the distribution is fully populated: the first xJ edge at or above the recoil-jet
   threshold over the lowest photon pT. It follows the binning by itself and is never typed in. */
#pragma once

#include <map>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace fa {

/* Ascending bin boundaries. Energy axes use GeV; xJ is dimensionless. Rebinning uses exact unions of stored bins. */
using Edges = std::vector<double>;

/* One result: a collision system and centrality class with its own grids and unfolding settings. */
struct SystemSettings {
  std::string name;          /* pp, auau_cent0_20, ... */
  std::string package;       /* PP or AuAu: the package directory */
  std::vector<std::string> simulation;  /* the simulated samples of the response (<package>/<sample>.root), summed */
  std::string centrality;    /* all, cent0_20, ... */
  std::string label;         /* plot label, e.g. "Au+Au 0-20%" */
  Edges photonEdges;         /* measurement photon ET bins (GeV) */
  Edges recoEdges;           /* response grid, reconstructed photon ET: package bins, same range */
  Edges truthEdges;          /* response grid, truth photon pT: package bins, same range */
  int pairIterations = 2;    /* iterations of the joint (photon pT, xJ) unfolding; 0: the frozen iteration rule */
  int photonIterations = 2;  /* iterations of the photon-yield unfolding; 0: the frozen iteration rule */
  std::vector<double> backgroundCorrelation;  /* R_bkg per photon bin (1: the background factorises) */
  std::vector<double> purityScale;            /* a factor on the ABCD signal per photon bin (1: none) */
  std::vector<double> efficiencyScale;        /* DATA/simulation photon ID and isolation efficiency per photon bin (1: none) */
  bool backgroundCorrelationFromQA = false;   /* "qa": R_bkg of the simulated non-signal candidates, from the QA file */
  bool purityScaleFromQA = false;             /* "qa": the simulated purity non-closure, from the QA file */
  bool subtractCombinatoric = true;           /* subtract the combinatoric jets (by Combinatoric.Method) */
  double reportedLow = 0;    /* derived, never set: the first xJ edge at or above jet threshold / lowest photon pT */
  bool subtractTemplate = true;               /* derived: subtractCombinatoric with the simulated template */
  std::vector<int> slices;                    /* derived: the 5% centrality slices of the class (-1: p+p) */
};

struct Settings {
  std::string name;
  std::string inputDirectory;      /* holds <package>/data.root and <package>/photonjet.root */
  std::string output;              /* results.root */
  std::string qaOutput;            /* qa.root: the inputs of the photon and event QA figures */
  std::string jet = "jet5";        /* recoil-jet slot of the packages, jet<GeV> */
  double jetThreshold = 5;         /* derived from the slot name */
  Edges xjEdges;
  double reportedHigh = 1.49;      /* upper edge of the reported xJ range; the lower edge is derived per result */
  int dataToys = 500;
  unsigned seed = 1;
  std::string dataToyModel = "poisson";      /* poisson: independent bins; event_bootstrap: Poisson(1) per data event */
  int responseToys = 0;                      /* source bootstrap of the simulated response (0: none) */
  unsigned responseSeed = 1;
  int closureToys = 0;                       /* source resamples of the estimator closure (0: no closure) */
  unsigned closureSeed = 1;
  int ruleDataToys = 2000, ruleResponseToys = 1000;  /* the frozen iteration rule's own toy ensembles */
  unsigned ruleDataSeed = 1, ruleResponseSeed = 1;
  std::string combinatoricMethod = "simulation";  /* simulation template, or data_dphi (Delta-phi sideband of the data) */
  std::string combinatoricShape = "flow";    /* data_dphi: flow (nominal) or flat */
  std::string combinatoricScope = "pt";      /* data_dphi: shape fitted per photon bin (pt) or per photon and xJ bin */
  int scanIterations = 8;          /* iterations shown in the stability scans */
  std::string priorForm = "exp1";  /* joint-unfolding prior: simulation, or the data-fitted exp1 / exp2 reweighting */
  std::string fakeModel = "scaled";  /* detector fakes under the reweighted truth: scaled (fake fraction) or fixed */
  std::string fitRange = "all";    /* prior fit bins: all measured bins with an uncertainty, or reported */
  double priorPower = 1.0;         /* the unfolding starts from w^power (1 nominal; 0.5 and 1.5 bracket it) */
  std::vector<SystemSettings> systems;
  std::map<std::string, std::vector<std::string>> inclusive;  /* QA: the inclusive-jet samples of each package */

  static Settings read(const std::string &path);
};

/* The measured data of one class on the measurement grid. */
struct Measured {
  TMatrixD abcd, abcdVariance;              /* 4 regions (A B C D) x photon bins */
  TMatrixD pairsA, pairsAVariance;          /* photon bins x xJ bins, region A recoil jets */
  TMatrixD pairsC, pairsCVariance;          /* region C (non-tight, isolated) recoil jets */
  TMatrixD backgroundA, backgroundC;        /* data-driven combinatoric jets already removed from pairsA and pairsC */
};

/* The event records of one class on the measurement grid, for the whole-event bootstrap: every photon and jet row keeps
   its event, so one multiplier per event moves that event's photon counts and recoil jets together. */
struct EventTable {
  int nPt = 0, nXj = 0;
  std::vector<double> weight;                               /* per selected event */
  std::vector<int> photonEvent, photonCell;                 /* cell: region * nPt + photon bin */
  std::vector<int> jetEvent, jetCell;                       /* recoil jets: (region * nPt + photon bin) * nXj + xJ bin */
  std::vector<int> dphiEvent, dphiCell;                     /* every jet: (jet cell * nDphi + Delta-phi bin) * 2 + far */
};

/* The simulated response restricted to the measurement window (unfolding.restrict). */
struct Window {
  TMatrixD pairResponse;                    /* truth cells x reco cells, cell = photon bin * nXj + xJ bin */
  TVectorD pairMisses, pairTruth, pairReco, pairFakes;   /* pairFakes: detector plus boundary fakes */
  TMatrixD combinatoric, combinatoricVariance;           /* photon bins x xJ bins */
  TMatrixD photonResponse;                  /* truth x reco photon bins */
  TVectorD photonMisses, photonTruth, photonReco, photonFakes;
  TMatrixD leakage, leakageVariance;        /* 4 x photon bins, truth-matched signal in A B C D */
  TVectorD detector;                        /* reco cells: linked photon and jet that are not the selected truth pair */
  TMatrixD boundaryByXj;                    /* xJ bins x reco cells: matched pairs whose truth photon is outside the window */
  TVectorD boundaryFixed;                   /* boundary migration without a truth xJ bin (expected zero) */
  TVectorD simPhotons;                      /* simulated reco photons in the window (matched plus boundary) */
  TMatrixD responseGridPhoton;              /* photon response on the full response grid (for display) */
  Edges recoEdges, truthEdges;
};

struct Corrected {
  TVectorD photons, purity, transfer;       /* per photon bin: signal, purity, background transfer k */
  TMatrixD pairs, pairsVariance;            /* purity- and combinatoric-corrected recoil jets */
  TMatrixD scaledC, purityCorrected, combinatoric;
  std::vector<bool> physical;               /* per photon bin: the ABCD solution is the physical root */
};

/* Inputs.cc */
std::string inputCoverage(const Settings &s);  /* verify package/selection compatibility and retain population coverage */
std::string renderMetadata(const Settings &s);
void renderProductReceipt(const std::string &path, const std::string &metadata);
Measured readMeasured(const Settings &s, const SystemSettings &sys);
EventTable readEvents(const Settings &s, const SystemSettings &sys);
Measured fromEvents(const EventTable &t, const std::vector<double> *multipliers, const Settings &s);
Window readWindow(const Settings &s, const SystemSettings &sys);
/* The source bootstrap of the response: every simulated source's own window, and the window of weighted occurrences. */
struct SourceWindow { Window window; double weight = 1; int sample = 0; std::string identity; };  /* identity: the source's sha256 */
std::vector<SourceWindow> readSourceWindows(const Settings &s, const SystemSettings &sys);
Window combine(const std::vector<SourceWindow> &sources, const std::vector<double> &multiplicity);
/* Rows are target bins, columns are fine bins; entries select whole bins, with no interpolation. */
std::vector<std::vector<double>> projection(const Edges &fine, const Edges &target);

/* Purity.cc */
struct ABCDSolution { double signal = 0, purity = 0, transfer = 0; bool ok = false; };
ABCDSolution solveLeakageABCD(const double counts[4], const double leakage[4], double backgroundCorrelation = 1.0);
Corrected correct(const Measured &m, const Window &w, const SystemSettings &sys);

/* The data-driven combinatoric-jet estimate from the photon-jet Delta-phi distributions (the sPHENIX Au+Au dijet
   method). dphi[region][((pt * nXj + xj) * nDphi + dphiBin) * 2 + side] holds the jets around the region A (0) and C (1)
   photons, side 0 near (|Delta eta| <= 0.8), 1 far; variance likewise. shape: flow or flat; scope: pt or pt_xj.
   Returns, per region, the combinatoric jets in the recoil window (photon bins x xJ bins) and their variance. */
struct Background { TMatrixD counts[2], variance[2]; };
const Edges &dphiEdges();
Background combinatoricFromDphi(const std::vector<double> dphi[2], const std::vector<double> variance[2], int nPt, int nXj,
                                const std::string &shape, const std::string &scope);

/* Unfolding.cc */
/* Response rows are truth, columns are reco. The returned vector contains truth yields; refolded contains reco yields. */
TVectorD bayes(const TVectorD &measured, const TMatrixD &responseTruthXReco, const TVectorD &misses, int iterations,
               const TVectorD &prior, const TVectorD &fakes, TVectorD *refolded = nullptr);
TVectorD unfoldPhotons(const TVectorD &photons, const Window &w, int iterations);
/* pairCells[photonBin * nXj + xjBin] is summed over photon bins, then divided by total photons and xJ bin width. */
TVectorD perPhoton(const TVectorD &pairCells, const TVectorD &photons, const Edges &xjEdges);

/* The frozen iteration rule (Unfolding.cc). One problem: the photon yield or the joint (photon pT, xJ) pairs. */
struct Problem {
  TMatrixD response;
  TVectorD misses, prior, fakes;
};
Problem photonProblem(const Window &w);
Problem pairProblem(const Window &w);
struct ResponseToy { Problem problem; TVectorD measured; };  /* an accepted response toy: its problem and corrected data */
struct ScanRow {
  int k = 0;
  double dataStat = 0, responseStat = 0, deviation = 0, objective = 0, refold = 0, mcClosure = 0, mcRefold = 0;
  double negativeFraction = 0, dataSuccess = 0, responseSuccess = 0;
  bool passes = false;
  std::string failures;
};
struct IterationScan {
  std::vector<ScanRow> rows;            /* k = 1..12 */
  std::vector<TVectorD> unfolded, refolded;
  int selected = 0;                     /* 0: no k of 2..4 passes */
};
IterationScan scanIterations(const Problem &p, const TVectorD &measured, const std::vector<TVectorD> &dataToys, int dataRequested,
                             const std::vector<ResponseToy> &responseToys, int responseRequested, const std::vector<bool> &mask);
int photonRefoldSelection(const IterationScan &scan, const TVectorD &measured, const std::vector<TVectorD> &dataToys,
                          const std::vector<TVectorD> &responseMeasured);

/* Prior.cc: the data-driven prior of the joint unfolding (the simulated truth reweighted in xJ) */
struct PriorFit { TVectorD theta, weights; double chi2 = 0; int ndf = 0; bool ok = false; };
TVectorD priorWeights(const TVectorD &theta, const Edges &xj, double reportedLow, double reportedHigh);
PriorFit fitPrior(const Window &w, const TVectorD &cells, const TVectorD &photons, const TMatrixD &covariance,
                  const Settings &s, const SystemSettings &sys);
TVectorD unfoldPairsWithPrior(const TVectorD &cells, const Window &w, int iterations, const TVectorD &weights,
                              const std::string &fakeModel, TVectorD *refolded = nullptr);

/* Analysis.cc */
int run(const Settings &s);

/* Systematics.cc: the systematic covariances and bands from results files of one change each */
int systematics(const std::string &configPath);

/* QA.cc: the inputs of the photon and event QA figures (qa/<package>/ of Output.QA) */
void writeQA(const Settings &s, TDirectory *results);
int runQA(const Settings &s);

}  /* namespace fa */
