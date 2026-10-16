/* FinalAnalysis: the per-photon xJgamma measurement from TreeToHists histogram packages. */

/* The chain, one function per step (src/):
   Inputs.cc      the settings, the measured histograms or event records, and the simulated response of one class
   Purity.cc      leakage-aware ABCD purity, background correlation R_bkg, combinatoric subtraction
   Unfolding.cc   the iterative-Bayes unfolding and the per-photon normalisation
   Prior.cc       the data-fitted prior of the joint unfolding
   Analysis.cc    the ordered chain, the data and response toys, the iteration rule and results.root
   Systematics.cc the systematic uncertainties from results of one change each
   Every physics choice is a setting in the configuration file (config/analysis.yaml); nothing below chooses one. The
   reported xJ range starts where the distribution is fully populated: the first xJ edge at or above the recoil-jet
   threshold over the lowest photon pT. It follows the binning by itself and is never typed in. */
#pragma once

#include <set>
#include <map>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "../../TreeToHists/src/ScientificTrace.h"

namespace fa {

/* Scientific ancestry travels with each numerical member, not with a process-wide
   last-read file. It is metadata only: arithmetic and serialized values are unchanged. */
using ScientificLineage = std::map<std::string, photonjet_science::Trace>;
inline photonjet_science::Trace lineageOf(const ScientificLineage &lineage, const std::string &member) {
  const auto found = lineage.find(member);
  return found == lineage.end() ? photonjet_science::Trace::unavailable("missing producer lineage: " + member) : found->second;
}

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
  int pairIterations = 2;    /* iterations of the joint (photon pT, xJ) unfolding; 0: the smallest Q_k of its scan */
  int photonIterations = 2;  /* iterations of the photon-yield unfolding; 0: the smallest Q_k of its scan */
  std::vector<double> backgroundCorrelation;  /* R_bkg per photon bin (1: the background factorises) */
  std::vector<double> purityScale;            /* a factor on the ABCD signal per photon bin (1: none) */
  std::vector<double> efficiencyScale;        /* DATA/simulation photon ID and isolation efficiency per photon bin (1: none) */
  bool backgroundCorrelationFromQA = false;   /* "qa": R_bkg of the simulated non-signal candidates, from the QA file */
  bool purityScaleFromQA = false;             /* "qa": the simulated purity non-closure, from the QA file */
  bool subtractCombinatoric = true;           /* subtract the combinatoric jets (by Combinatoric.Method) */
  double reportedLow = 0;    /* derived, never set: the first xJ edge at or above jet threshold / lowest photon pT */
  bool subtractTemplate = true;               /* derived: subtractCombinatoric with the simulated template */
  bool subtractSideband = false;              /* derived: data_dphi with Combinatoric.Shape simulation */
  std::vector<int> slices;                    /* derived: the 5% centrality slices of the class (-1: p+p) */
  std::vector<int> excludedDataRuns;  /* optional DATA-only run veto; retained event positions keep paired bootstrap draws aligned */
  std::string responseComponent, responseRecipe, responseGridSha256; /* explicit saved response component; empty is nominal */
  ScientificLineage lineage; /* derived QA parameter ancestry, populated when a setting is resolved from QA */
  photonjet_science::Json parameterOrigins=photonjet_science::Json::object();
};

struct Settings {
  photonjet_science::Json parameterOrigins=photonjet_science::Json::object();
  photonjet_science::Json nativeParameterExports=photonjet_science::Json::object();
  bool correctionsOnly = false;   /* explicit diagnostic entry point: serialize pre-unfolding objects, never a final result */
  std::string name;
  std::string inputDirectory;      /* holds <package>/data.root and <package>/photonjet.root */
  std::string inputCollection;     /* required typed collection for named response components */
  std::map<std::string,std::string> namedResponsePins; /* derived package byte pins; never a user setting */
  std::string output;              /* results.root */
  std::string qaOutput;            /* qa.root: the inputs of the photon and event QA figures */
  std::string qaMode = "full";      /* measurement_inputs: varied r_bkg/purity_scale only; diagnostics stay nominal */
  std::string photonIDModels;      /* optional sealed PhotonID models.yaml for saved training QA */
  std::string jet = "jet5";        /* recoil-jet slot of the packages, jet<GeV> */
  double jetThreshold = 5;         /* derived from the slot name */
  Edges xjEdges;
  /* Optional (Association.*): finer photon and xJ bins on which the reconstructed-to-truth association figure
     (IAN Figure 25) is saved, for the listed results. Empty: not saved. Nothing of the measurement uses them. */
  Edges associationPhotonEdges, associationXjEdges;
  std::vector<std::string> associationSystems;
  double reportedHigh = 1.49;      /* upper edge of the reported xJ range; the lower edge is derived per result */
  bool saveToyRecords = false;     /* optional finite unfolded replicas with original attempt IDs for paired comparisons */
  int dataToys = 500;
  unsigned seed = 1;
  std::string dataToyModel = "poisson";      /* poisson: independent bins; event_bootstrap: Poisson(1) per data event */
  int responseToys = 0;                      /* source bootstrap of the simulated response (0: none) */
  unsigned responseSeed = 1;
  int closureToys = 0;                       /* source resamples of the estimator closure (0: no closure) */
  unsigned closureSeed = 1;
  int ruleDataToys = 2000, ruleResponseToys = 1000;  /* the iteration scan's own toy ensembles */
  unsigned ruleDataSeed = 1, ruleResponseSeed = 1;
  std::string rulePrior = "fitted";  /* the prior of the rule's joint scan: fitted (the measurement's) or simulation */
  std::string combinatoricMethod = "simulation";  /* simulation template, or data_dphi (Delta-phi sideband of the data) */
  std::string combinatoricShape = "simulation";  /* data_dphi: simulation (nominal: the simulated transfer), flow or flat */
  std::string combinatoricScope = "pt";      /* flow, flat: shape fitted per photon bin (pt) or per photon and xJ bin */
  double sidebandLow = 0.8, sidebandHigh = 2.5;  /* data_dphi: the Delta-phi sideband [rad], two Delta-phi bin edges */
  int scanIterations = 8;          /* iterations shown in the stability scans */
  std::string priorForm = "exp1";  /* joint-unfolding prior: simulation, or the data-fitted exp1 / exp2 reweighting */
  std::string fakeModel = "scaled";  /* detector fakes under the reweighted truth: scaled (fake fraction) or fixed */
  std::string fitRange = "all";    /* prior fit bins: all measured bins with an uncertainty, or reported */
  double priorFraction = 1.0;      /* the unfolding starts from 1 + f (w - 1): f = 1 nominal, f = 0.5 the prior
                                      systematic (PPG14's half reweighting) */
  std::vector<SystemSettings> systems;
  std::string ppStitching = "nominal";  /* reference: SI-only generator closure and |vz| < 30 cm tag comparison */
  std::map<std::string, std::vector<std::string>> inclusive;  /* QA: the inclusive-jet samples of each package */
  Edges qaAbcdEdges{15,17,19,21,23,25,35};
  Edges qaShowerShapes{15,20}, qaIsolationOverlay{14,16}, qaNcbScore{16,35};
  std::vector<Edges> qaTruthIsolation{{15,20},{25,30}};

  static Settings read(const std::string &path);
};
/* Settings enter the dependency graph only at the operation that consumes them.
   The binding retains the native evaluated source value even if a derived runtime
   value (for example scan-selected k or a QA scale) replaces a struct member. */
inline photonjet_science::Trace withSettings(photonjet_science::Trace trace,
    const photonjet_science::Json &bindings,std::initializer_list<std::string> keys) {
  for(const auto &key:keys) {
    if(!bindings.contains(key)) {
      trace=photonjet_science::Trace::combine("missing_configuration_witness",
        {trace,photonjet_science::Trace::unavailable("configuration not witnessed: "+key)});
      continue;
    }
    const auto &binding=bindings.at(key);
    trace=trace.withParameter(binding.at("reference"),binding.at("category"),binding.at("value"),binding.value("origin",photonjet_science::Json()));
  }
  return trace;
}
inline photonjet_science::Json scientificParameterBindings(const Settings &settings) {
  photonjet_science::Json all={{"global",settings.parameterOrigins},{"systems",photonjet_science::Json::object()}};
  for(const auto &sys:settings.systems)all["systems"][sys.name]=sys.parameterOrigins;
  return all;
}
void validateNamedResponseCollection(Settings &settings); /* existing Pipeline ingestion, no fits or event reduction */

/* The measured data of one class on the measurement grid. */
struct Measured {
  TMatrixD abcd, abcdVariance;              /* 4 regions (A B C D) x photon bins */
  TMatrixD pairsA, pairsAVariance;          /* photon bins x xJ bins, region A recoil jets */
  TMatrixD pairsC, pairsCVariance;          /* region C (non-tight, isolated) recoil jets */
  TMatrixD backgroundA, backgroundC;        /* data-driven combinatoric jets already removed from pairsA and pairsC */
  TMatrixD sidebandA, sidebandC;            /* Shape simulation: every jet around the region A and C photons in the
                                               Delta-phi sideband, photon bins x xJ bins; correct() subtracts */
  ScientificLineage lineage;
};

/* The event records of one class on the measurement grid, for the whole-event bootstrap: every photon and jet row keeps
   its event, so one multiplier per event moves that event's photon counts and recoil jets together. */
struct EventTable {
  int nPt = 0, nXj = 0;
  size_t excludedEvents = 0;
  std::vector<double> weight;                               /* per selected event */
  std::vector<int> photonEvent, photonCell;                 /* cell: region * nPt + photon bin */
  std::vector<int> jetEvent, jetCell;                       /* recoil jets: (region * nPt + photon bin) * nXj + xJ bin */
  std::vector<int> dphiEvent, dphiCell;                     /* every jet: (jet cell * nDphi + Delta-phi bin) * 2 + far */
  ScientificLineage lineage;
};

/* The simulated response restricted to the measurement window (unfolding.restrict). */
struct Window {
  TMatrixD pairResponse;                    /* truth cells x reco cells, cell = photon bin * nXj + xJ bin */
  TVectorD pairMisses, pairTruth, pairReco, pairFakes;   /* pairFakes: detector plus boundary fakes */
  TMatrixD combinatoric, combinatoricVariance;           /* photon bins x xJ bins */
  TMatrixD sidebandJets;                    /* photon bins x xJ bins: every jet of a matched region-A photon in the
                                               Delta-phi sideband (with or without a truth jet nearby) */
  TMatrixD photonResponse;                  /* truth x reco photon bins */
  TVectorD photonMisses, photonTruth, photonReco, photonFakes;
  TMatrixD leakage, leakageVariance;        /* 4 x photon bins, truth-matched signal in A B C D */
  TVectorD detector;                        /* reco cells: linked photon and jet that are not the selected truth pair */
  TMatrixD boundaryByXj;                    /* xJ bins x reco cells: matched pairs whose truth photon is outside the window */
  TVectorD boundaryFixed;                   /* boundary migration without a truth xJ bin (expected zero) */
  TVectorD simPhotons;                      /* simulated reco photons in the window (matched plus boundary) */
  TMatrixD responseGridPhoton;              /* photon response on the full response grid (for display) */
  Edges recoEdges, truthEdges;
  ScientificLineage lineage;
};

struct Corrected {
  TVectorD photons, purity, transfer;       /* per photon bin: signal, purity, background transfer k */
  TMatrixD pairs, pairsVariance;            /* purity- and combinatoric-corrected recoil jets */
  TMatrixD scaledC, scaledCVariance;        /* the stages the figures show: region C times the transfer, */
  TMatrixD purityCorrected, purityCorrectedVariance;  /* the recoil jets after the purity correction, */
  TMatrixD combinatoric;                    /* and the simulated combinatoric template taken off after it */
  std::vector<bool> physical;               /* per photon bin: the ABCD solution is the physical root */
};

/* Inputs.cc */
std::string inputCoverage(const Settings &s, bool allowMissingData = false);  /* verify package/selection compatibility and retain population coverage */
std::string renderMetadata(const Settings &s, bool allowMissingData = false);
void renderProductReceipt(const std::string &path, const std::string &metadata);
Measured readMeasured(const Settings &s, const SystemSettings &sys);
EventTable readEvents(const Settings &s, const SystemSettings &sys);
Measured fromEvents(const EventTable &t, const std::vector<double> *multipliers, const Settings &s);
Window readWindow(const Settings &s, const SystemSettings &sys);
/* The simulated photon-jet Delta-phi distributions of a result, for the Delta-phi figure: every jet above the
   threshold around a truth-matched region-A photon inside the result's photon and xJ bins, per bin of dphiEdges().
   truthJet: jets with a truth jet nearby; combinatoric: jets with none. Simulation units, with their variance; the
   photons they belong to are row A of the saved leakage. */
struct SimulatedDphi { TVectorD truthJet, truthJetVariance, combinatoric, combinatoricVariance; ScientificLineage lineage; };
SimulatedDphi readSimulatedDphi(const Settings &s, const SystemSettings &sys);
/* The source bootstrap of the response. A simulated source is one produced file of a sample, and its window is what
   that file alone contributes. Windows add, so every source is kept as one sparse row of numbers (the additive
   members of Window, flattened in the order of `blocks` in Inputs.cc) and the window of any set of occurrences is a
   single pass over the filled entries. Reading the rows out of the packages is the slowest step of an analysis, about
   a minute, so readSources saves them and reads the packages again only when they, the grids or the input code change. */
struct Sources {
  Window layout;                      /* an empty window of the result's dimensions */
  std::vector<int> sample;            /* per source: its sample, an index of SystemSettings::simulation */
  std::vector<double> weight;         /* per source: its sample's cross section per generated event */
  std::vector<std::string> identity;  /* per source: the sha256 of the produced file */
  std::vector<int> begin;             /* per source: its first entry; begin[size()] ends the last source */
  std::vector<int> column;            /* per entry: its position in the flattened window */
  std::vector<double> value;          /* per entry: the source's own count or variance, before any weight */
  int counts = 0;                     /* columns below this hold counts, the others variances */
  size_t size() const { return sample.size(); }
};
/* `window`: the result's window from the summed samples (readWindow). The rows must add up to it in every number. */
Sources readSources(const Settings &s, const SystemSettings &sys, const Window &window);
/* The window of weighted source occurrences: counts times multiplicity x weight, variances times multiplicity x weight^2. */
Window combine(const Sources &sources, const std::vector<double> &multiplicity);
/* Rows are target bins, columns are fine bins; entries select whole bins, with no interpolation. */
std::vector<std::vector<double>> projection(const Edges &fine, const Edges &target);

/* Purity.cc */
struct ABCDSolution { double signal = 0, purity = 0, transfer = 0; bool ok = false; };
ABCDSolution solveLeakageABCD(const double counts[4], const double leakage[4], double backgroundCorrelation = 1.0);
Corrected correct(const Measured &m, const Window &w, const SystemSettings &sys);

/* The data-driven combinatoric-jet estimate from the photon-jet Delta-phi distributions (the sPHENIX Au+Au dijet
   method). dphi[region][((pt * nXj + xj) * nDphi + dphiBin) * 2 + side] holds the jets around the region A (0) and C (1)
   photons, side 0 near (|Delta eta| <= 0.8), 1 far; variance likewise. shape: flow or flat; scope: pt or pt_xj;
   the sideband is sidebandLow < Delta phi < sidebandHigh and the flow fit runs over 0 < Delta phi < sidebandHigh.
   Returns, per region, the combinatoric jets in the recoil window (photon bins x xJ bins) and their variance.
   With Combinatoric.Shape simulation the sideband jets are kept in Measured and correct() applies the simulated transfer. */
struct Background { TMatrixD counts[2], variance[2]; };
const Edges &dphiEdges();
Background combinatoricFromDphi(const std::vector<double> dphi[2], const std::vector<double> variance[2], int nPt, int nXj,
                                const std::string &shape, const std::string &scope, double sidebandLow, double sidebandHigh);

/* Unfolding.cc */
/* Response rows are truth, columns are reco. The returned vector contains truth yields; refolded contains reco yields. */
TVectorD bayes(const TVectorD &measured, const TMatrixD &responseTruthXReco, const TVectorD &misses, int iterations,
               const TVectorD &prior, const TVectorD &fakes, TVectorD *refolded = nullptr);
TVectorD unfoldPhotons(const TVectorD &photons, const Window &w, int iterations);
/* pairCells[photonBin * nXj + xjBin] is summed over photon bins, then divided by total photons and xJ bin width. */
TVectorD perPhoton(const TVectorD &pairCells, const TVectorD &photons, const Edges &xjEdges);

/* The iteration scan (Unfolding.cc). One problem: the photon yield or the joint (photon pT, xJ) pairs. */
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
  int selected = 0;                     /* the k >= 2 with the smallest Q_k; 0: none is finite */
};
IterationScan scanIterations(const Problem &p, const TVectorD &measured, const std::vector<TVectorD> &dataToys, int dataRequested,
                             const std::vector<ResponseToy> &responseToys, int responseRequested, const std::vector<bool> &mask);

/* Prior.cc: the data-driven prior of the joint unfolding (the simulated truth reweighted in xJ) */
struct PriorFit { TVectorD theta, weights; double chi2 = 0; int ndf = 0; bool ok = false; };
/* Copy every object of a directory (sub-directories and trees included) into another; `skip` holds paths relative
   to `from` that are left out. */
void copyDirectory(TDirectory *from, TDirectory *to, const std::set<std::string> &skip = {}, const std::string &prefix = "");
/* The studies that need a DATA package's per-source records: the photon-trigger study (out/triggers) and the Au+Au
   centrality target (out/centrality_target). A compact package carries their result in derived/. */
void sourceStudies(TFile *data, const std::string &package, TDirectory *out);
extern const char *const kSourceStudies;  /* the version of sourceStudies a compact package was made with */
/* A twin of a validated collection without the DATA packages' per-source records (final_analysis.sh compact). */
int compactCollection(const std::string &collection, const std::string &output);
TVectorD priorWeights(const TVectorD &theta, const Edges &xj, double reportedLow, double reportedHigh);
/* The joint unfolding problem that starts from the simulated truth times weights per xJ bin (as unfoldPairsWithPrior). */
Problem pairProblemWithPrior(const Window &w, const TVectorD &weights, const std::string &fakeModel);
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
/* Import a sealed training-diagnostics collection without histogram inputs. */
int runModelQA(const std::string &collection, const std::string &output);

}  /* namespace fa */
