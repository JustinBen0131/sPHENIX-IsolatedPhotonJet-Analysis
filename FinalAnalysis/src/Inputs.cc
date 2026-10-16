/* Inputs: the settings file, the measured histograms and the simulated response of one centrality class. */

/* A histogram package keeps every quantity on its own fine bins. Coarser analysis bins are exact unions of those bins
   (projection); a boundary that the package does not have is an error, never an interpolation. */
#include "FinalAnalysis.h"
#include "NativeProducerEvidence.h"
#include "RenderMetadata.h"
#include "MeasurementReplay.h"
#include "NamedResponse.h"
#include "RootScientificTrace.h"
#include "ConfigurationScientificWitness.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include <unistd.h>

#include "TEnv.h"
#include "../../TreeToHists/src/Configuration.h"
#include "TFile.h"
#include "TMD5.h"
#include "THashList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TObjString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

/* A hash of this file and FinalAnalysis.h, passed by final_analysis.sh: it names the saved source rows (readSources),
   so they are made again when the code that makes them changes. Without it nothing is saved. */
#ifndef FINALANALYSIS_INPUT_CODE
#define FINALANALYSIS_INPUT_CODE ""
#endif

namespace fa {
namespace science = photonjet_science;
constexpr char kPpStitchingDefault[] = "nominal";
std::string inputCoverage(const Settings &s, bool allowMissingData) {
  using Json = nlohmann::json;
  Json coverage = Json::object(); std::set<std::string> seen;
  Json replay; bool sawNominal=false, sawReplay=false;
  std::map<std::string, Json> selections; bool partial = false, unknown = false, legacy = false;
  for (const auto &sys : s.systems) {
    std::vector<std::string> names{"data"}; names.insert(names.end(), sys.simulation.begin(), sys.simulation.end());
    const auto inclusive = s.inclusive.find(sys.package);
    if (inclusive != s.inclusive.end()) names.insert(names.end(), inclusive->second.begin(), inclusive->second.end());
    for (const auto &name : names) {
      const std::string key = sys.package + "/" + name;
      if (!seen.insert(key).second) continue;
      const std::string path = s.inputDirectory + "/" + key + ".root";
      if (allowMissingData && name == "data" && !std::filesystem::exists(path)) {
        coverage[key] = {{"status", "MISSING"}};
        partial = true;
        continue;
      }
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
      auto *saved = dynamic_cast<TObjString *>(file->Get("meta/package"));
      if (!saved) throw std::runtime_error(path + ": missing sample-package metadata");
      const Json meta = Json::parse(saved->GetString().Data());
      if (meta.value("schema", "") != "PhotonJetSampleV1" || meta.value("system", "") != sys.package || meta.value("name", "") != name)
        throw std::runtime_error(path + ": package identity differs from its analysis binding");
      if (meta.value("legacy", false)) { legacy = true; coverage[key] = {{"status", "LEGACY_REFERENCE"}}; continue; }
      const Json recipe=measurementReplayRecipe(meta);
      validateMeasurementReplayConsumerScope(recipe,s.combinatoricMethod,s.combinatoricShape);
      if(recipe.is_null())sawNominal=true;
      else {
        if(file->GetDirectory("records/candidates"))throw std::runtime_error("replay candidate archive must remain in nominal parent");
        if(sawReplay && replay!=recipe)throw std::runtime_error("analysis inputs mix measurement replay recipes");
        replay=recipe;sawReplay=true;
      }
      if(sawReplay && sawNominal)throw std::runtime_error("analysis inputs mix replayed and nominal packages");
      /* DATA is read first, so a package's simulation is checked against the selection its DATA evaluated. */
      const Json &selection = meta.at("selection");
      if (selection.value("selection_sha256", "").empty() ||
          (selections.count(sys.package) && !photonjet_render::selectionsAgree(selection, selections.at(sys.package), false)))
        throw std::runtime_error(path + ": DATA/SIM selections differ");
      selections.emplace(sys.package, selection);
      if (name != "data") {
        auto *weight = dynamic_cast<TParameter<double> *>(file->Get("meta/pb_per_generated_event"));
        const double recorded = meta.at("normalization").at("pb_per_generated_event");
        if (!weight || !(recorded > 0) || !std::isfinite(recorded) || weight->GetVal() != recorded)
          throw std::runtime_error(path + ": SIM normalization metadata differs from its applied weight");
      }
      coverage[key] = meta.value("coverage", Json{{"status", "UNKNOWN"}});
      /* A combined DATA package counts its declared excluded runs in the roster it is compared with. When every run
         it combines is complete and together they are all of its sources, the runs analysed are completely covered. */
      if (coverage[key].value("status", "") == "PARTIAL" && meta.contains("run_coverages") && !meta.at("run_coverages").empty() &&
          meta.contains("excluded_runs") && !meta.at("excluded_runs").empty()) {
        long long sources = 0;
        bool complete = true;
        for (const auto &run : meta.at("run_coverages")) {
          complete = complete && run.value("status", "") == "COMPLETE" && run.at("population_sources") == run.at("selected_sources");
          sources += run.value("selected_sources", 0LL);
        }
        if (complete && sources == coverage[key].value("selected_sources", -1LL) && sources == meta.value("source_count", -2LL)) {
          coverage[key]["roster_status"] = "PARTIAL";
          coverage[key]["excluded_runs"] = meta.at("excluded_runs");
          coverage[key]["status"] = "COMPLETE";
        }
      }
      const std::string status = coverage[key].value("status", "UNKNOWN");
      partial = partial || status == "PARTIAL"; unknown = unknown || status != "PARTIAL" && status != "COMPLETE";
    }
  }
  Json result={{"status", legacy ? "REFERENCE_INPUTS" : partial ? "PARTIAL" : unknown ? "UNKNOWN" : "COMPLETE"}, {"samples", coverage}};
  if(sawReplay)result["measurement_replay"]=replay;
  for(const auto &sys:s.systems)if(named_response::active(sys))result["named_response_components"][sys.name]=named_response::selection(sys);
  return result.dump();
}

namespace {

/* A setting's value: the text after "Key:", up to a trailing "# comment". Every setting is required. */
std::string text(TEnv &env, const std::string &key) {
  std::string value = env.GetValue(key.c_str(), "");
  value = value.substr(0, value.find('#'));
  const auto first = value.find_first_not_of(" \t"), last = value.find_last_not_of(" \t");
  if (first == std::string::npos) throw std::runtime_error("setting " + key + " is missing");
  return value.substr(first, last - first + 1);
}

Edges numbers(TEnv &env, const std::string &key) {
  std::istringstream in(text(env, key));
  Edges out;
  for (std::string item; in >> item;) {
    size_t used = 0;
    double x = 0;
    try { x = std::stod(item, &used); } catch (const std::exception &) { used = 0; }
    if (used != item.size() || !std::isfinite(x)) throw std::runtime_error("setting " + key + ": '" + item + "' is not a number");
    out.push_back(x);
  }
  return out;
}

double number(TEnv &env, const std::string &key) {
  const Edges v = numbers(env, key);
  if (v.size() != 1) throw std::runtime_error("setting " + key + " takes one number");
  return v[0];
}

int integer(TEnv &env, const std::string &key) {
  const double v = number(env, key);
  if (v != std::floor(v)) throw std::runtime_error("setting " + key + " takes an integer");
  return static_cast<int>(v);
}

std::string choice(TEnv &env, const std::string &key, const std::set<std::string> &allowed) {
  const std::string v = text(env, key);
  if (!allowed.count(v)) throw std::runtime_error("setting " + key + ": '" + v + "' is not one of the allowed values");
  return v;
}

Edges increasing(TEnv &env, const std::string &key) {
  const Edges e = numbers(env, key);
  if (e.size() < 2) throw std::runtime_error("setting " + key + " needs at least two edges");
  for (size_t i = 1; i < e.size(); ++i)
    if (!(e[i] > e[i - 1])) throw std::runtime_error("setting " + key + " must increase");
  return e;
}

bool isEdge(const Edges &edges, double x) {
  for (double e : edges)
    if (std::fabs(e - x) < 1e-9) return true;
  return false;
}

/* The recoil-jet threshold is the slot's name: the package contract stores recoil jets above t GeV under jet<t>/. */
double slotThreshold(const std::string &slot) {
  size_t used = 0;
  double t = 0;
  try { t = slot.rfind("jet", 0) == 0 ? std::stod(slot.substr(3), &used) : 0; } catch (const std::exception &) { used = 0; }
  if (used == 0 || used != slot.size() - 3 || !(t > 0)) throw std::runtime_error("Recoil.Jet must name a package slot jet<GeV>, not " + slot);
  return t;
}

/* The 5% centrality slices of a class: all (p+p) is slice -1; centA_B holds slices A/5 .. B/5 - 1. */
std::vector<int> slicesOf(const std::string &centrality) {
  if (centrality == "all") return {-1};
  const auto u = centrality.find('_');
  int low = -1, high = -1;
  try { if (centrality.rfind("cent", 0) == 0 && u != std::string::npos) { low = std::stoi(centrality.substr(4, u - 4)); high = std::stoi(centrality.substr(u + 1)); } } catch (const std::exception &) {}
  if (low < 0 || high <= low || low % 5 || high % 5) throw std::runtime_error("centrality class " + centrality + " is not all or centA_B in 5% steps");
  std::vector<int> out;
  for (int k = low / 5; k < high / 5; ++k) out.push_back(k);
  return out;
}

/* The fully populated xJ region starts at the first xJ edge at or above (jet threshold) / (lowest photon pT): below it
   the jet threshold removes recoil jets of the lowest photons, and the distribution is a kinematic turn-on. */
double populatedLow(const Edges &xj, double threshold, double photonLow) {
  const double turnOn = threshold / photonLow;
  for (double e : xj)
    if (e >= turnOn - 1e-9) return e;
  throw std::runtime_error("no xJ edge lies above the jet-threshold turn-on");
}

Edges axisEdges(const TAxis *a) {
  Edges e;
  for (int i = 1; i <= a->GetNbins() + 1; ++i) e.push_back(a->GetBinLowEdge(i));
  return e;
}

template <class T> T *get(TFile &f, const std::string &path) {
  auto *h = dynamic_cast<T *>(f.Get(path.c_str()));
  if (!h) throw std::runtime_error(std::string(f.GetName()) + ": missing " + path);
  return h;
}

std::unique_ptr<TFile> open(const std::string &path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  observeNativeInput(*f);
  return f;
}

/* Pair cells c = photon bin * nX + xJ bin: the matrix that takes fine cells to window cells. */
TMatrixD cellProjection(const std::vector<std::vector<double>> &pt, const std::vector<std::vector<double>> &xj) {
  const int nPt = pt.size(), nFinePt = pt[0].size(), nX = xj.size(), nFineX = xj[0].size();
  TMatrixD m(nPt * nX, nFinePt * nFineX);
  for (int i = 0; i < nPt; ++i)
    for (int f = 0; f < nFinePt; ++f)
      if (pt[i][f] > 0)
        for (int x = 0; x < nX; ++x)
          for (int g = 0; g < nFineX; ++g) m(i * nX + x, f * nFineX + g) = xj[x][g];
  return m;
}

TMatrixD toMatrix(const std::vector<std::vector<double>> &p) {
  TMatrixD m(p.size(), p[0].size());
  for (size_t i = 0; i < p.size(); ++i)
    for (size_t j = 0; j < p[0].size(); ++j) m(i, j) = p[i][j];
  return m;
}

TVectorD column(const TMatrixD &m) {  /* 1 - (sum of each column): the fine cells outside the window */
  TVectorD out(m.GetNcols());
  for (int j = 0; j < m.GetNcols(); ++j) {
    double s = 0;
    for (int i = 0; i < m.GetNrows(); ++i) s += m(i, j);
    out[j] = 1.0 - s;
  }
  return out;
}

TVectorD values1(const TH1 *h, bool variance = false) {
  TVectorD v(h->GetNbinsX());
  for (int i = 0; i < v.GetNrows(); ++i) {
    const double e = h->GetBinError(i + 1);
    v[i] = variance ? e * e : h->GetBinContent(i + 1);
  }
  return v;
}

TMatrixD values2(const TH2 *h, bool variance = false) {  /* x bins x y bins */
  TMatrixD m(h->GetNbinsX(), h->GetNbinsY());
  for (int i = 0; i < m.GetNrows(); ++i)
    for (int j = 0; j < m.GetNcols(); ++j) {
      const double e = h->GetBinError(i + 1, j + 1);
      m(i, j) = variance ? e * e : h->GetBinContent(i + 1, j + 1);
    }
  return m;
}

TVectorD cellsOf(const TH2 *h, bool variance = false) {  /* (photon bin, xJ bin) flattened, photon bin outer */
  TMatrixD m = values2(h, variance);
  TVectorD v(m.GetNrows() * m.GetNcols());
  for (int i = 0; i < m.GetNrows(); ++i)
    for (int j = 0; j < m.GetNcols(); ++j) v[i * m.GetNcols() + j] = m(i, j);
  return v;
}

TMatrixD reshape(const TVectorD &v, int rows, int cols) {
  TMatrixD m(rows, cols);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) m(i, j) = v[i * cols + j];
  return m;
}

/* The cells of cellsOf for a (photon pT, xJ, Delta phi) histogram summed over the Delta-phi sideband. */
TVectorD sidebandCells(const TH3 *h, double low, double high) {
  const int nx = h->GetNbinsX(), ny = h->GetNbinsY();
  TVectorD v(nx * ny);
  for (int z = 1; z <= h->GetNbinsZ(); ++z) {
    if (h->GetZaxis()->GetBinLowEdge(z) < low - 1e-9 || h->GetZaxis()->GetBinUpEdge(z) > high + 1e-9) continue;
    for (int i = 1; i <= nx; ++i)
      for (int j = 1; j <= ny; ++j) v[(i - 1) * ny + j - 1] += h->GetBinContent(i, j, z);
  }
  return v;
}

}  /* namespace */

Settings Settings::read(const std::string &path) {
  TEnv env;
  photonjet_config::readSettings(path, env);
  Settings s;
  s.name = text(env, "Analysis.Name");
  s.inputDirectory = text(env, "Input.Directory");
  if(env.Defined("Input.Collection"))s.inputCollection=std::filesystem::absolute(std::filesystem::path(path).parent_path()/text(env,"Input.Collection")).lexically_normal().string();
  s.output = text(env, "Output.Results");
  s.qaOutput = text(env, "Output.QA");
  if(env.Defined("Output.ToyRecords")) {
    const int enabled=integer(env,"Output.ToyRecords");
    if(enabled!=0 && enabled!=1)throw std::runtime_error("Output.ToyRecords is 0 or 1");
    s.saveToyRecords=enabled==1;
  }
  if(env.Defined("QA.Mode"))s.qaMode=choice(env,"QA.Mode",{"full","measurement_inputs"});
  if(env.Defined("QA.AbcdEdges"))s.qaAbcdEdges=increasing(env,"QA.AbcdEdges");
  auto window=[&](const std::string &key,Edges &range){
    if(env.Defined(key.c_str()))range=increasing(env,key);
    if(range.size()!=2 || range.front()<0)throw std::runtime_error(key+" requires two nonnegative increasing bounds in GeV");
  };
  if(s.qaAbcdEdges.front()<0)throw std::runtime_error("QA.AbcdEdges must be nonnegative");
  window("QA.Windows.ShowerShapes",s.qaShowerShapes);
  window("QA.Windows.IsolationOverlay",s.qaIsolationOverlay);
  window("QA.Windows.NcbScore",s.qaNcbScore);
  if(env.Defined("QA.Windows.TruthIsolation")){
    const auto pairs=numbers(env,"QA.Windows.TruthIsolation");
    if(pairs.empty() || pairs.size()%2)throw std::runtime_error("QA.Windows.TruthIsolation requires low high pairs in GeV");
    s.qaTruthIsolation.clear();
    for(size_t i=0;i<pairs.size();i+=2){
      if(pairs[i]<0 || !(pairs[i]<pairs[i+1]))throw std::runtime_error("QA.Windows.TruthIsolation bounds must increase within each pair");
      s.qaTruthIsolation.push_back({pairs[i],pairs[i+1]});
    }
  }
  s.photonIDModels=env.GetValue("QA.PhotonIDModels","");
  if (s.photonIDModels=="none") s.photonIDModels.clear();
  if (!s.photonIDModels.empty()) s.photonIDModels=std::filesystem::absolute(std::filesystem::path(path).parent_path()/s.photonIDModels).string();
  s.jet = text(env, "Recoil.Jet");
  s.jetThreshold = slotThreshold(s.jet);
  s.xjEdges = increasing(env, "Measurement.XjEdges");
  s.reportedHigh = number(env, "Measurement.ReportedXjMax");
  if (!isEdge(s.xjEdges, s.reportedHigh)) throw std::runtime_error("Measurement.ReportedXjMax must be an xJ bin edge");
  if (env.Defined("Association.Systems")) {  /* all three keys, or none */
    s.associationPhotonEdges = increasing(env, "Association.PhotonEdges");
    s.associationXjEdges = increasing(env, "Association.XjEdges");
    std::istringstream listed(text(env, "Association.Systems"));
    for (std::string name; listed >> name;) s.associationSystems.push_back(name);
  }
  s.dataToys = integer(env, "Uncertainty.DataToys");
  s.seed = integer(env, "Uncertainty.Seed");
  s.dataToyModel = choice(env, "Uncertainty.DataToyModel", {"poisson", "event_bootstrap"});
  s.responseToys = integer(env, "Uncertainty.ResponseToys");
  s.responseSeed = integer(env, "Uncertainty.ResponseSeed");
  s.closureToys = integer(env, "Closure.Toys");
  s.closureSeed = integer(env, "Closure.Seed");
  if (s.closureToys < 0 || s.closureToys == 1) throw std::runtime_error("Closure.Toys is 0 or at least 2");
  s.ruleDataToys = integer(env, "IterationRule.DataToys");
  s.ruleDataSeed = integer(env, "IterationRule.DataSeed");
  s.ruleResponseToys = integer(env, "IterationRule.ResponseToys");
  s.ruleResponseSeed = integer(env, "IterationRule.ResponseSeed");
  s.rulePrior = choice(env, "IterationRule.Prior", {"fitted", "simulation"});
  if (s.ruleDataToys < 2 || s.ruleResponseToys < 2) throw std::runtime_error("the iteration rule needs at least two data and two response toys");
  if (s.responseToys < 0 || s.responseToys == 1) throw std::runtime_error("Uncertainty.ResponseToys is 0 or at least 2");
  s.combinatoricMethod = choice(env, "Combinatoric.Method", {"simulation", "data_dphi"});
  s.combinatoricShape = choice(env, "Combinatoric.Shape", {"simulation", "flow", "flat"});
  {  /* the sideband: two Delta-phi bin edges below the recoil window, as written or to three decimals */
    const Edges sideband = numbers(env, "Combinatoric.Sideband");
    auto edge = [](double x) {
      for (double e : dphiEdges()) if (std::fabs(e - x) < 2e-3) return e;
      throw std::runtime_error("Combinatoric.Sideband names a value that is not a Delta-phi bin edge (0.8, 2.5 or a multiple of pi/32)");
    };
    if (sideband.size() != 2) throw std::runtime_error("Combinatoric.Sideband is [low, high] in rad");
    s.sidebandLow = edge(sideband[0]), s.sidebandHigh = edge(sideband[1]);
    if (!(s.sidebandLow < s.sidebandHigh) || s.sidebandHigh > 7 * M_PI / 8 + 1e-9)
      throw std::runtime_error("Combinatoric.Sideband lies below the recoil window, low < high");
  }
  s.combinatoricScope = choice(env, "Combinatoric.Scope", {"pt", "pt_xj"});
  if (s.combinatoricMethod == "data_dphi" && s.dataToyModel != "event_bootstrap")
    throw std::runtime_error("Combinatoric.Method data_dphi is refitted per event: it needs Uncertainty.DataToyModel event_bootstrap");
  s.scanIterations = integer(env, "Unfolding.ScanIterations");
  s.priorForm = choice(env, "Prior.Form", {"simulation", "exp1", "exp2"});
  s.fakeModel = choice(env, "Prior.FakeModel", {"scaled", "fixed"});
  s.fitRange = choice(env, "Prior.FitRange", {"all", "reported"});
  s.priorFraction = number(env, "Prior.Fraction");
  if (s.dataToys < 2 || s.scanIterations < 1 || !(s.priorFraction > 0 && s.priorFraction <= 1))
    throw std::runtime_error("DataToys >= 2, ScanIterations >= 1, 0 < Prior.Fraction <= 1");
  std::set<std::string> known = {"Analysis.Name", "Input.Directory", "Input.Collection", "Output.Results", "Output.QA", "Output.ToyRecords", "QA.Mode", "QA.PhotonIDModels", "Recoil.Jet", "Measurement.XjEdges",
                                 "Measurement.ReportedXjMax", "Association.Systems", "Association.PhotonEdges", "Association.XjEdges",
                                 "Uncertainty.DataToys", "Uncertainty.Seed",
                                 "Uncertainty.DataToyModel", "Uncertainty.ResponseToys", "Uncertainty.ResponseSeed", "Closure.Toys", "Closure.Seed",
                                 "IterationRule.DataToys", "IterationRule.DataSeed", "IterationRule.ResponseToys", "IterationRule.ResponseSeed", "IterationRule.Prior",
                                 "Combinatoric.Method", "Combinatoric.Shape", "Combinatoric.Scope", "Combinatoric.Sideband",
                                 "Unfolding.ScanIterations", "Prior.Form", "Prior.FakeModel", "Prior.FitRange", "Prior.Fraction",
                                 "Systems", "QA.AbcdEdges", "QA.Windows.ShowerShapes", "QA.Windows.IsolationOverlay",
                                 "QA.Windows.TruthIsolation", "QA.Windows.NcbScore"};
  std::istringstream names(text(env, "Systems"));
  for (std::string name; names >> name;) {
    for (const char *key : {"Package", "Centrality", "Label", "Simulation", "PhotonEdges", "RecoEdges", "TruthEdges", "PairIterations",
                            "PhotonIterations", "BackgroundCorrelation", "PurityScale", "EfficiencyScale", "SubtractCombinatoric", "ResponseComponent", "ResponseRecipe", "ResponseGridSHA256", "ExcludeDataRuns"})
      known.insert(name + "." + key);
    SystemSettings sys;
    sys.name = name;
    if(env.Defined((name+".ResponseComponent").c_str()))sys.responseComponent=text(env,name+".ResponseComponent");
    if(env.Defined((name+".ResponseRecipe").c_str()))sys.responseRecipe=text(env,name+".ResponseRecipe");
    if(env.Defined((name+".ResponseGridSHA256").c_str()))sys.responseGridSha256=text(env,name+".ResponseGridSHA256");
    sys.package = text(env, name + ".Package");
    if(env.Defined((name+".ExcludeDataRuns").c_str())) {
      for(double run:numbers(env,name+".ExcludeDataRuns")) {
        if(run<=0 || run>std::numeric_limits<int>::max() || run!=std::floor(run))
          throw std::runtime_error(name+".ExcludeDataRuns needs positive integer run numbers");
        sys.excludedDataRuns.push_back(static_cast<int>(run));
      }
      std::sort(sys.excludedDataRuns.begin(),sys.excludedDataRuns.end());
      if(std::adjacent_find(sys.excludedDataRuns.begin(),sys.excludedDataRuns.end())!=sys.excludedDataRuns.end())
        throw std::runtime_error(name+".ExcludeDataRuns contains duplicate runs");
      if(s.dataToyModel!="event_bootstrap")throw std::runtime_error("ExcludeDataRuns needs retained event_bootstrap records");
    }
    sys.centrality = text(env, name + ".Centrality");
    sys.label = text(env, name + ".Label");
    sys.photonEdges = increasing(env, name + ".PhotonEdges");
    sys.recoEdges = increasing(env, name + ".RecoEdges");
    sys.truthEdges = increasing(env, name + ".TruthEdges");
    std::istringstream samples(text(env, name + ".Simulation"));
    for (std::string sample; samples >> sample;) sys.simulation.push_back(sample);
    /* a number, or "scan": the k with the smallest Q_k of the result's iteration scan */
    sys.pairIterations = text(env, name + ".PairIterations") == "scan" ? 0 : integer(env, name + ".PairIterations");
    sys.photonIterations = text(env, name + ".PhotonIterations") == "scan" ? 0 : integer(env, name + ".PhotonIterations");
    /* one positive value per photon bin; BackgroundCorrelation and PurityScale may instead be "qa": derived by the qa
       step from the simulation and read from the QA file */
    auto perBin = [&](const std::string &key, std::vector<double> &values, bool *fromQA) {
      if (fromQA && (*fromQA = text(env, name + "." + key) == "qa")) return;
      values = numbers(env, name + "." + key);
      if (values.size() + 1 != sys.photonEdges.size()) throw std::runtime_error(name + "." + key + " needs one value per photon bin");
      for (double v : values) if (!(v > 0) || !std::isfinite(v)) throw std::runtime_error(name + "." + key + " values are positive");
    };
    perBin("BackgroundCorrelation", sys.backgroundCorrelation, &sys.backgroundCorrelationFromQA);
    perBin("PurityScale", sys.purityScale, &sys.purityScaleFromQA);
    perBin("EfficiencyScale", sys.efficiencyScale, nullptr);
    sys.subtractCombinatoric = integer(env, name + ".SubtractCombinatoric") != 0;
    if (sys.pairIterations < 0 || sys.photonIterations < 0) throw std::runtime_error(name + ": iterations are a positive number or rule");
    sys.reportedLow = populatedLow(s.xjEdges, s.jetThreshold, sys.photonEdges.front());
    sys.subtractTemplate = sys.subtractCombinatoric && s.combinatoricMethod == "simulation";
    sys.subtractSideband = sys.subtractCombinatoric && s.combinatoricMethod == "data_dphi" && s.combinatoricShape == "simulation";
    if (s.combinatoricMethod == "data_dphi" && !sys.subtractCombinatoric)
      throw std::runtime_error(name + ": the data-driven combinatoric estimate is always subtracted");
    sys.slices = slicesOf(sys.centrality);
    named_response::validateChoice(sys);
    if (!(sys.reportedLow < s.reportedHigh)) throw std::runtime_error(name + ": the fully populated xJ region is empty");
    s.systems.push_back(sys);
  }
  if (s.systems.empty()) throw std::runtime_error("no Systems");
  /* the QA's inclusive-jet samples of each package ("none": no inclusive simulation) */
  std::set<std::string> packages;
  for (const auto &sys : s.systems) packages.insert(sys.package);
  for (const auto &package : packages) {
    known.insert("QA." + package + ".Inclusive");
    std::istringstream names(text(env, "QA." + package + ".Inclusive"));
    for (std::string name; names >> name;) if (name != "none") s.inclusive[package].push_back(name);
  }
  known.insert("QA.PP.Stitching");
  s.ppStitching = env.GetValue("QA.PP.Stitching", kPpStitchingDefault);
  if (s.ppStitching != "nominal" && s.ppStitching != "reference")
    throw std::runtime_error("QA.PP.Stitching is nominal or reference");
  /* a misspelt or leftover setting is an error, never silently ignored */
  TIter next(env.GetTable());
  while (auto *record = static_cast<TEnvRec *>(next()))
    if (!known.count(record->GetName())) throw std::runtime_error(path + ": unknown setting " + record->GetName());
  if(!s.inputCollection.empty() || std::any_of(s.systems.begin(),s.systems.end(),named_response::active))validateNamedResponseCollection(s);
  ConfigurationScientificWitness witness;witness.load(path);
  auto bind=[&](science::Json &bindings,const std::string &key,const std::string &sourceKey,const science::Json &value,const std::string &category="models") {
    bindings[key]={{"reference","FinalAnalysis.setting/"+sourceKey},{"category",category},{"value",value},{"origin",witness.evaluated(sourceKey,value)}};
  };
  auto global=[&](const std::string &key,const science::Json &value,const std::string &category="models"){bind(s.parameterOrigins,key,key,value,category);};
  global("Recoil.Jet",s.jet,"selection");global("Measurement.XjEdges",s.xjEdges,"projection");
  global("Measurement.ReportedXjMax",s.reportedHigh,"projection");
  global("Association.Systems",s.associationSystems,"selection");global("Association.PhotonEdges",s.associationPhotonEdges,"projection");
  global("Association.XjEdges",s.associationXjEdges,"projection");
  global("Uncertainty.DataToys",s.dataToys);global("Uncertainty.Seed",s.seed);global("Uncertainty.DataToyModel",s.dataToyModel);
  global("Uncertainty.ResponseToys",s.responseToys);global("Uncertainty.ResponseSeed",s.responseSeed);
  global("Closure.Toys",s.closureToys);global("Closure.Seed",s.closureSeed);
  global("IterationRule.DataToys",s.ruleDataToys);global("IterationRule.DataSeed",s.ruleDataSeed);
  global("IterationRule.ResponseToys",s.ruleResponseToys);global("IterationRule.ResponseSeed",s.ruleResponseSeed);global("IterationRule.Prior",s.rulePrior);
  global("Combinatoric.Method",s.combinatoricMethod);global("Combinatoric.Shape",s.combinatoricShape);global("Combinatoric.Scope",s.combinatoricScope);
  global("Combinatoric.Sideband",Edges{s.sidebandLow,s.sidebandHigh},"selection");
  global("Unfolding.ScanIterations",s.scanIterations);global("Prior.Form",s.priorForm);global("Prior.FakeModel",s.fakeModel);
  global("Prior.FitRange",s.fitRange);global("Prior.Fraction",s.priorFraction);
  global("QA.AbcdEdges",s.qaAbcdEdges,"projection");global("QA.Windows.ShowerShapes",s.qaShowerShapes,"projection");
  global("QA.Windows.IsolationOverlay",s.qaIsolationOverlay,"projection");global("QA.Windows.NcbScore",s.qaNcbScore,"projection");
  global("QA.Windows.TruthIsolation",s.qaTruthIsolation,"projection");global("QA.PP.Stitching",s.ppStitching,"normalization");
  for(const auto &package:packages)global("QA."+package+".Inclusive",s.inclusive[package],"selection");
  for(auto &sys:s.systems) {
    auto local=[&](const std::string &key,const science::Json &value,const std::string &category="models"){bind(sys.parameterOrigins,key,sys.name+"."+key,value,category);};
    local("Centrality",sys.centrality,"selection");local("Simulation",sys.simulation,"selection");
    if(!sys.excludedDataRuns.empty())local("ExcludeDataRuns",sys.excludedDataRuns,"selection");
    local("PhotonEdges",sys.photonEdges,"projection");local("RecoEdges",sys.recoEdges,"projection");local("TruthEdges",sys.truthEdges,"projection");
    local("PairIterations",sys.pairIterations);local("PhotonIterations",sys.photonIterations);
    local("BackgroundCorrelation",sys.backgroundCorrelationFromQA?science::Json("qa"):science::Json(sys.backgroundCorrelation));
    local("PurityScale",sys.purityScaleFromQA?science::Json("qa"):science::Json(sys.purityScale));
    local("EfficiencyScale",sys.efficiencyScale,"normalization");local("SubtractCombinatoric",sys.subtractCombinatoric);
    for(const auto *key:{"Combinatoric.Method","Combinatoric.Shape"})sys.parameterOrigins[key]=s.parameterOrigins.at(key);
  }
  // A genuinely omitted setting is still a concrete native default. Bind the
  // actual constant used above, rather than pretending the YAML contains it.
  if(s.parameterOrigins.at("QA.PP.Stitching").at("origin").is_null() && !env.Defined("QA.PP.Stitching")) {
    ConfigurationScientificWitness defaults;defaults.load(__FILE__);
    const std::string marker="constexpr char kPpStitchingDefault[] = ";
    const std::string literal=science::Json(kPpStitchingDefault).dump();
    const auto begin=defaults.source.find(marker);
    if(begin!=std::string::npos && defaults.source.compare(begin+marker.size(),literal.size(),literal)==0) {
      defaults.span("QA.PP.Stitching",begin+marker.size(),begin+marker.size()+literal.size());
      auto origin=defaults.sites.value("QA.PP.Stitching",science::Json());
      if(!origin.is_null())origin["encoding"]="string";
      s.parameterOrigins["QA.PP.Stitching"]["origin"]=origin;
    }
  }
  s.nativeParameterExports=witness.exports;
  return s;
}

/* Fine bins -> target bins: 1 where a fine bin lies inside a target bin. Every target edge must be a fine edge. */
std::vector<std::vector<double>> projection(const Edges &fine, const Edges &target) {
  for (double edge : target) {
    bool found = false;
    for (double f : fine) found = found || std::fabs(f - edge) < 1e-9;
    if (!found) throw std::runtime_error("boundary " + std::to_string(edge) + " is not a package bin edge");
  }
  std::vector<std::vector<double>> p(target.size() - 1, std::vector<double>(fine.size() - 1, 0.0));
  for (size_t i = 0; i + 1 < target.size(); ++i)
    for (size_t j = 0; j + 1 < fine.size(); ++j)
      p[i][j] = (fine[j] >= target[i] - 1e-9 && fine[j + 1] <= target[i + 1] + 1e-9) ? 1.0 : 0.0;
  return p;
}

namespace {

/* Reads the chosen columns of a record tree row by row, as doubles whatever their stored type. */
void readRows(TFile &f, const std::string &path, const std::vector<std::string> &names,
              const std::function<void(const std::vector<double> &)> &row) {
  auto *tree = dynamic_cast<TTree *>(f.Get(path.c_str()));
  if (!tree) throw std::runtime_error(std::string(f.GetName()) + ": missing " + path);
  tree->SetBranchStatus("*", false);
  std::vector<TLeaf *> leaves;
  for (const auto &n : names) {
    tree->SetBranchStatus(n.c_str(), true);
    TLeaf *leaf = tree->GetLeaf(n.c_str());
    if (!leaf) throw std::runtime_error(path + ": missing column " + n);
    leaves.push_back(leaf);
  }
  std::vector<double> values(names.size());
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    for (size_t k = 0; k < leaves.size(); ++k) values[k] = leaves[k]->GetValue();
    row(values);
  }
}

/* Measurement bin of a native ROOT bin number: the bin holding its lower edge; -1 outside the measurement window. */
int measurementBin(int rootBin, const Edges &native, const Edges &target) {
  if (rootBin < 1 || rootBin > static_cast<int>(native.size()) - 1) return -1;
  const double low = native[rootBin - 1];
  const int index = std::upper_bound(target.begin(), target.end(), low + 1e-9) - target.begin() - 1;
  return index >= 0 && index < static_cast<int>(target.size()) - 1 ? index : -1;
}

}  /* namespace */

/* Rebin retained DATA records while preserving event membership. One bootstrap multiplier must move an event's
   photons and jets together; independently resampling histogram cells would lose their correlations. */
EventTable readEvents(const Settings &s, const SystemSettings &sys) {
  /* TreeToHists selects the trigger and the sources once: the top-level records are the analysis population */
  auto f = open(s.inputDirectory + "/" + sys.package + "/data.root");
  const Edges native = axisEdges(get<TH2>(*f, "measurement/" + sys.centrality + "/abcd")->GetYaxis());
  /* the package's recoil-jet thresholds, in the order of the thresholds bits: its jet<GeV> slots, ascending */
  std::vector<double> thresholds;
  TIter next(get<TDirectory>(*f, "measurement/" + sys.centrality)->GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    const std::string n = key->GetName();
    if (n.rfind("jet", 0) == 0) thresholds.push_back(slotThreshold(n));
  }
  std::sort(thresholds.begin(), thresholds.end());
  const int threshold = std::find(thresholds.begin(), thresholds.end(), s.jetThreshold) - thresholds.begin();
  if (threshold == static_cast<int>(thresholds.size())) throw std::runtime_error("the package has no " + s.jet + " slot");
  EventTable t;
  /* Records have their own producer definitions. An old package without these
     witnesses remains unknown; a nearby histogram is not evidence for a tree. */
  for (const auto &name : {"photons", "jets", "events"}) {
    auto *definition = dynamic_cast<TObjString *>(f->Get((std::string("meta/record_definitions/") + name).c_str()));
    auto trace = definition ? science::Trace::source(science::Json::parse(definition->GetString().Data()))
                            : science::Trace::unavailable(std::string("record definition absent: ") + name);
    science::Json parameters={{"centrality_slices",sys.slices}};
    if (std::string(name) != "events") parameters["photon_edges"] = sys.photonEdges;
    if (std::string(name) == "jets") { parameters["xj_edges"] = s.xjEdges; parameters["jet_threshold_gev"] = s.jetThreshold; }
    trace=withSettings(trace,sys.parameterOrigins,{"Centrality"});
    if(std::string(name)!="events")trace=withSettings(trace,sys.parameterOrigins,{"PhotonEdges"});
    if(std::string(name)=="jets")trace=withSettings(trace,s.parameterOrigins,{"Measurement.XjEdges","Recoil.Jet"});
    if(!sys.excludedDataRuns.empty()) {
      parameters["excluded_data_runs"]=sys.excludedDataRuns;
      parameters["run_veto_algorithm"]="zero_event_weight_preserve_bootstrap_positions";
      trace=withSettings(trace,sys.parameterOrigins,{"ExcludeDataRuns"});
    }
    t.lineage[name] = trace.then("select_event_records", parameters);
  }
  t.nPt = sys.photonEdges.size() - 1;
  t.nXj = s.xjEdges.size() - 1;
  const int nDphi = dphiEdges().size() - 1;
  std::vector<int> renumber;
  std::vector<std::string> eventColumns={"slice","weight"};
  if(!sys.excludedDataRuns.empty())eventColumns.push_back("run");
  readRows(*f, "records/events", eventColumns, [&](const std::vector<double> &v) {
    const bool selected = std::find(sys.slices.begin(), sys.slices.end(), static_cast<int>(v[0])) != sys.slices.end();
    renumber.push_back(selected ? static_cast<int>(t.weight.size()) : -1);
    if (selected) {
      bool excluded=false;
      if(!sys.excludedDataRuns.empty()) {
        if(!std::isfinite(v[2]) || v[2]<=0 || v[2]>std::numeric_limits<int>::max() || v[2]!=std::floor(v[2]))
          throw std::runtime_error("run veto encountered an invalid retained event run");
        excluded=std::binary_search(sys.excludedDataRuns.begin(),sys.excludedDataRuns.end(),static_cast<int>(v[2]));
      }
      // Keep every selected event in place: nominal and subset consume the same
      // Poisson draw for each original event, including events with zero weight.
      t.weight.push_back(excluded ? 0.0 : std::isfinite(v[1]) ? v[1] : 1.0);
      t.excludedEvents+=excluded;
    }
  });
  auto event = [&](double index) { const long i = static_cast<long>(index); return i >= 0 && i < static_cast<long>(renumber.size()) ? renumber[i] : -1; };
  readRows(*f, "records/photons", {"event", "region", "et_bin"}, [&](const std::vector<double> &v) {
    const int e = event(v[0]), pt = measurementBin(static_cast<int>(v[2]), native, sys.photonEdges);
    if (e < 0 || pt < 0) return;
    t.photonEvent.push_back(e);
    t.photonCell.push_back(static_cast<int>(v[1]) * t.nPt + pt);
  });
  readRows(*f, "records/jets", {"event", "region", "et_bin", "xj", "dphi", "far", "window", "thresholds"}, [&](const std::vector<double> &v) {
    const int e = event(v[0]), pt = measurementBin(static_cast<int>(v[2]), native, sys.photonEdges);
    const int x = std::upper_bound(s.xjEdges.begin(), s.xjEdges.end(), v[3]) - s.xjEdges.begin() - 1;
    int passed = 0;
    for (long mask = static_cast<long>(v[7]); mask > 0; mask >>= 1) ++passed;  /* bit length: the thresholds are nested */
    if (e < 0 || pt < 0 || x < 0 || x >= t.nXj || passed <= threshold) return;
    const int cell = (static_cast<int>(v[1]) * t.nPt + pt) * t.nXj + x;
    if (static_cast<int>(v[6]) == 1) { t.jetEvent.push_back(e); t.jetCell.push_back(cell); }
    const int d = std::upper_bound(dphiEdges().begin(), dphiEdges().end(), v[4]) - dphiEdges().begin() - 1;
    if (d >= 0 && d < nDphi) { t.dphiEvent.push_back(e); t.dphiCell.push_back((cell * nDphi + d) * 2 + static_cast<int>(v[5])); }
  });
  return t;
}

/* The measured statistics of the (resampled) events; the variance of every count is the count, as for unit-weight
   Poisson fills. With Combinatoric.Method data_dphi the data-driven combinatoric jets are removed from the recoil jets. */
Measured fromEvents(const EventTable &t, const std::vector<double> *multipliers, const Settings &s) {
  auto weight = [&](int e) { return t.weight[e] * (multipliers ? (*multipliers)[e] : 1.0); };
  Measured m{TMatrixD(4, t.nPt), TMatrixD(4, t.nPt), TMatrixD(t.nPt, t.nXj), TMatrixD(t.nPt, t.nXj),
             TMatrixD(t.nPt, t.nXj), TMatrixD(t.nPt, t.nXj), TMatrixD(t.nPt, t.nXj), TMatrixD(t.nPt, t.nXj)};
  const auto photons = science::Trace::combine("count_record_photons", {lineageOf(t.lineage,"photons"),lineageOf(t.lineage,"events")},
    {{"bootstrap",multipliers != nullptr},{"variance","count"}});
  const auto jets = science::Trace::combine("count_record_recoil_jets", {lineageOf(t.lineage,"jets"),lineageOf(t.lineage,"events")},
    {{"bootstrap",multipliers != nullptr},{"variance","count"}});
  m.lineage["abcd"] = photons;
  for (const auto &region : {"A","C"}) {
    m.lineage[std::string("rawPairs")+region] = jets.then("select_region",{{"region",region}});
    m.lineage[std::string("pairs")+region] = m.lineage[std::string("rawPairs")+region];
    if (s.combinatoricMethod == "data_dphi") {
      auto background = jets.then("dphi_background",{{"region",region},{"shape",s.combinatoricShape},
        {"scope",s.combinatoricScope},{"sideband_low",s.sidebandLow},{"sideband_high",s.sidebandHigh}});
      background=withSettings(background,s.parameterOrigins,{"Combinatoric.Method","Combinatoric.Shape","Combinatoric.Sideband"});
      if(s.combinatoricShape!="simulation")background=withSettings(background,s.parameterOrigins,{"Combinatoric.Scope"});
      m.lineage[std::string("sideband")+region] = background;
      m.lineage[std::string("background")+region] = background;
      if (s.combinatoricShape != "simulation")
        m.lineage[std::string("pairs")+region] = science::Trace::combine("subtract_data_dphi_background",
          {m.lineage[std::string("rawPairs")+region],background});
    }
  }
  for (size_t i = 0; i < t.photonEvent.size(); ++i) m.abcd(t.photonCell[i] / t.nPt, t.photonCell[i] % t.nPt) += weight(t.photonEvent[i]);
  for (size_t i = 0; i < t.jetEvent.size(); ++i) {
    const int region = t.jetCell[i] / (t.nPt * t.nXj), rest = t.jetCell[i] % (t.nPt * t.nXj);
    (region == 0 ? m.pairsA : m.pairsC)(rest / t.nXj, rest % t.nXj) += weight(t.jetEvent[i]);
  }
  m.abcdVariance = m.abcd;
  m.pairsAVariance = m.pairsA;
  m.pairsCVariance = m.pairsC;
  if (s.combinatoricMethod == "data_dphi") {
    const int nDphi = dphiEdges().size() - 1, perRegion = t.nPt * t.nXj * nDphi * 2;
    std::vector<double> dphi[2] = {std::vector<double>(perRegion), std::vector<double>(perRegion)};
    for (size_t i = 0; i < t.dphiEvent.size(); ++i) dphi[t.dphiCell[i] / perRegion][t.dphiCell[i] % perRegion] += weight(t.dphiEvent[i]);
    if (s.combinatoricShape == "simulation") {
      /* the jets of the Delta-phi sideband at every Delta eta; correct() applies the simulated transfer to them */
      for (TMatrixD *x : {&m.sidebandA, &m.sidebandC, &m.backgroundA, &m.backgroundC}) x->ResizeTo(t.nPt, t.nXj);
      const Edges &e = dphiEdges();
      for (int r = 0; r < 2; ++r)
        for (int cell = 0; cell < t.nPt * t.nXj; ++cell)
          for (int d = 0; d < nDphi; ++d)
            if (e[d] >= s.sidebandLow - 1e-12 && e[d + 1] <= s.sidebandHigh + 1e-12)
              (r == 0 ? m.sidebandA : m.sidebandC)(cell / t.nXj, cell % t.nXj) += dphi[r][(cell * nDphi + d) * 2] + dphi[r][(cell * nDphi + d) * 2 + 1];
      return m;
    }
    const Background b = combinatoricFromDphi(dphi, dphi, t.nPt, t.nXj, s.combinatoricShape, s.combinatoricScope, s.sidebandLow, s.sidebandHigh);
    m.backgroundA = b.counts[0];
    m.backgroundC = b.counts[1];
    m.pairsA -= b.counts[0];
    m.pairsC -= b.counts[1];
    m.pairsAVariance += b.variance[0];
    m.pairsCVariance += b.variance[1];
  }
  return m;
}

Measured readMeasured(const Settings &s, const SystemSettings &sys) {
  if(!sys.excludedDataRuns.empty())return fromEvents(readEvents(s,sys),nullptr,s);
  auto f = open(s.inputDirectory + "/" + sys.package + "/data.root");
  const std::string base = "measurement/" + sys.centrality + "/";
  auto *abcd = get<TH2>(*f, base + "abcd");
  auto *recoil = get<TH3>(*f, base + s.jet + "/recoil");
  const auto pAbcd = projection(axisEdges(abcd->GetYaxis()), sys.photonEdges);
  const auto pReco = projection(axisEdges(recoil->GetYaxis()), sys.photonEdges);
  const Edges fineXj = axisEdges(recoil->GetZaxis());
  if (std::fabs(fineXj.front() - s.xjEdges.front()) > 1e-9 || std::fabs(fineXj.back() - s.xjEdges.back()) > 1e-9)
    throw std::runtime_error("xJ bins must span the package xJ range");
  const auto pXj = projection(fineXj, s.xjEdges);
  const int nPt = sys.photonEdges.size() - 1, nX = s.xjEdges.size() - 1;
  Measured m{TMatrixD(4, nPt), TMatrixD(4, nPt), TMatrixD(nPt, nX), TMatrixD(nPt, nX), TMatrixD(nPt, nX), TMatrixD(nPt, nX),
             TMatrixD(nPt, nX), TMatrixD(nPt, nX)};
  m.lineage["abcd"] = science::projectionOf(abcd, {{{"axis",1},{"low",sys.photonEdges.front()},{"high",sys.photonEdges.back()}}})
    .then("exact_union_rebin",{{"photon_edges",sys.photonEdges}});
  m.lineage["abcd"]=withSettings(m.lineage["abcd"],sys.parameterOrigins,{"Centrality","PhotonEdges"});
  for (int region=0;region<2;++region) {
    const std::string name = region == 0 ? "pairsA" : "pairsC";
    m.lineage[name] = science::projectionOf(recoil, {{{"axis",0},{"low",recoil->GetXaxis()->GetBinLowEdge(region+1)},{"high",recoil->GetXaxis()->GetBinUpEdge(region+1)}},
      {{"axis",1},{"low",sys.photonEdges.front()},{"high",sys.photonEdges.back()}},{{"axis",2},{"low",s.xjEdges.front()},{"high",s.xjEdges.back()}}})
      .then("exact_union_rebin",{{"photon_edges",sys.photonEdges},{"xj_edges",s.xjEdges}});
    m.lineage[name]=withSettings(withSettings(m.lineage[name],sys.parameterOrigins,{"Centrality","PhotonEdges"}),s.parameterOrigins,{"Recoil.Jet","Measurement.XjEdges"});
    m.lineage[region == 0 ? "rawPairsA" : "rawPairsC"] = m.lineage[name];
  }
  for (int r = 0; r < 4; ++r)
    for (int i = 0; i < nPt; ++i)
      for (size_t j = 0; j < pAbcd[i].size(); ++j)
        if (pAbcd[i][j] > 0) {
          const double e = abcd->GetBinError(r + 1, j + 1);
          m.abcd(r, i) += abcd->GetBinContent(r + 1, j + 1);
          m.abcdVariance(r, i) += e * e;
        }
  for (int region = 0; region < 2; ++region) {
    TMatrixD &v = region == 0 ? m.pairsA : m.pairsC;
    TMatrixD &w = region == 0 ? m.pairsAVariance : m.pairsCVariance;
    for (int i = 0; i < nPt; ++i)
      for (size_t j = 0; j < pReco[i].size(); ++j)
        if (pReco[i][j] > 0)
          for (int x = 0; x < nX; ++x)
            for (size_t g = 0; g < pXj[x].size(); ++g)
              if (pXj[x][g] > 0) {
                const double e = recoil->GetBinError(region + 1, j + 1, g + 1);
                v(i, x) += recoil->GetBinContent(region + 1, j + 1, g + 1);
                w(i, x) += e * e;
              }
  }
  return m;
}

/* The response restricted to the measurement window (FinalAnalysis unfolding.restrict): reco migration out of the
   window joins the misses, truth migration into it joins the boundary fakes; counts are added before anything
   nonlinear. Restricting from the fine bins equals restricting from the response grid, which is a union of them. */
namespace {

/* The simulated samples of a result (<Package>/<sample>.root), each weighted by its cross section per generated event
   (meta/pb_per_generated_event, which TreeToHists records and never applies): every response histogram is their sum. */
/* A filled cell of the joint response: truth cell, reco cell (both from zero) and its content. */
struct ResponseCell { int truth, reco; double value; };

class Simulation {
 public:
  std::vector<ResponseCell> responseCells(const std::string &path) {
    auto *h = get<TH2>(path);
    std::vector<ResponseCell> cells;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
      for (int j = 1; j <= h->GetNbinsY(); ++j)
        if (const double v = h->GetBinContent(i, j); v != 0) cells.push_back({i - 1, j - 1, v});
    return cells;
  }
  Simulation(const Settings &s, const SystemSettings &sys) {
    for (const auto &name : sys.simulation) {
      files_.push_back(open(s.inputDirectory + "/" + sys.package + "/" + name + ".root"));
      named_response::requirePinned(s,sys,name,*files_.back());
      auto *w = dynamic_cast<TParameter<double> *>(files_.back()->Get("meta/pb_per_generated_event"));
      if (!w || !(w->GetVal() > 0) || !std::isfinite(w->GetVal()))
        throw std::runtime_error(std::string(files_.back()->GetName()) + ": no positive meta/pb_per_generated_event");
      weights_.push_back(w->GetVal());
    }
  }
  template <class T> T *get(const std::string &path) {
    auto &sum = cache_[path];
    for (size_t k = 0; !cached_.count(path) && k < files_.size(); ++k) {
      T *h = fa::get<T>(*files_[k], path);
      if (!sum) {
        sum.reset(static_cast<TH1 *>(h->Clone()));
        sum->SetDirectory(nullptr);
        sum->Scale(weights_[k]);
        science::sumDefinition(sum.get(),h,weights_[k],true);
      } else if (!sum->Add(h, weights_[k])) {
        throw std::runtime_error("the simulated samples disagree in the binning of " + path);
      } else science::sumDefinition(sum.get(),h,weights_[k],false);
    }
    cached_.insert(path);
    return static_cast<T *>(sum.get());
  }
  science::Trace trace(const std::string &path) { return science::traceOf(get<TH1>(path)); }

 private:
  std::vector<std::unique_ptr<TFile>> files_;
  std::vector<double> weights_;
  std::map<std::string, std::unique_ptr<TH1>> cache_;
  std::set<std::string> cached_;
};

/* One source's response histograms rebuilt from its sparse cells (records/source_contributions of its package). */
class SourceHistograms {
 public:
  /* `sparse`: the joint response, kept as its filled cells (a dense copy per source costs more than everything else) */
  SourceHistograms(TFile &package, const std::vector<std::string> &paths, const std::string &sparse) : sparse_(sparse) {
    for (const auto &path : paths) {
      traces_[path] = science::traceOf(fa::get<TH1>(package,path));
      if (path == sparse) { stride_ = fa::get<TH2>(package, path)->GetNbinsX() + 2; continue; }
      auto *h = static_cast<TH1 *>(fa::get<TH1>(package, path)->Clone());
      h->SetDirectory(nullptr);
      h->Reset();
      hists_[path].reset(h);
    }
  }
  std::vector<ResponseCell> responseCells(const std::string &path) const {
    if (path != sparse_) throw std::runtime_error("no per-source response cells for " + path);
    return cells_;
  }
  template <class T> T *get(const std::string &path) {
    auto it = hists_.find(path);
    if (it == hists_.end()) throw std::runtime_error("no per-source histogram " + path);
    return static_cast<T *>(it->second.get());
  }
  science::Trace trace(const std::string &path) const { return lineageOf(traces_,path); }
  void clear() { for (auto &item : hists_) item.second->Reset(); cells_.clear(); }
  void set(const std::string &path, int bin, double value, double variance) {
    if (path == sparse_) {  /* a global bin of a TH2: x + (nx + 2) y, flow bins included */
      const int x = bin % stride_, y = bin / stride_;
      if (x >= 1 && y >= 1 && value != 0) cells_.push_back({x - 1, y - 1, value});
      return;
    }
    auto it = hists_.find(path);
    if (it == hists_.end()) return;  /* a cell of another histogram */
    it->second->SetBinContent(bin, value);
    it->second->GetSumw2()->SetAt(variance, bin);
  }

 private:
  ScientificLineage traces_;
  std::map<std::string, std::unique_ptr<TH1>> hists_;
  std::string sparse_;
  int stride_ = 1;
  std::vector<ResponseCell> cells_;
};

/* Everything of the settings that a window depends on. windowOf takes this and nothing else of the settings, so the
   description that names a saved source table (readSources) cannot miss one. */
struct Grid {
  std::string centrality, jet;       /* the class and the recoil-jet slot of the packages */
  std::string responsePrefix, responseIdentity;
  science::Json settings,system;
  Edges photon, xj;                  /* the measurement bins */
  Edges reco, truth;                 /* the display grid of the photon response */
  double sidebandLow, sidebandHigh;  /* the Delta-phi sideband of the simulated transfer [rad] */
  Grid(const Settings &s, const SystemSettings &sys)
      : centrality(sys.centrality), jet(s.jet), photon(sys.photonEdges), xj(s.xjEdges), reco(sys.recoEdges),
        truth(sys.truthEdges), sidebandLow(s.sidebandLow), sidebandHigh(s.sidebandHigh) {
    settings=s.parameterOrigins;system=sys.parameterOrigins;
    responsePrefix=named_response::prefix(sys);if(named_response::active(sys))responseIdentity=named_response::selection(sys).dump();
  }
  std::string text() const {
    std::ostringstream out;
    out.precision(17);
    out << centrality << " " << jet << " sideband " << sidebandLow << " " << sidebandHigh;
    for (const auto &[name, edges] : {std::pair<const char *, const Edges *>{"photon", &photon}, {"xj", &xj},
                                      {"reco", &reco}, {"truth", &truth}}) {
      out << " " << name;
      for (double edge : *edges) out << " " << edge;
    }
    if(!responseIdentity.empty())out<<" named_response "<<responseIdentity;
    return out.str();
  }
};

template <class Provider> Window windowOf(Provider &f, const Grid &g);

/* Every additive member of a window in one fixed order, the counts first and the two variances last. Flattening a
   window, rebuilding one and the saved source table walk this list and nothing else does, so a new member of Window
   is added here alone. visit(its numbers, how many, whether they are variances). */
template <class W, class Visit> void blocks(W &w, Visit visit) {
  /* Every member of Window by position, so the build stops here when a member is added and has no place below.
     recoEdges and truthEdges are the display grid, not sums. */
  [[maybe_unused]] auto &[pairResponse, pairMisses, pairTruth, pairReco, pairFakes, combinatoric, combinatoricVariance,
                          sidebandJets, photonResponse, photonMisses, photonTruth, photonReco, photonFakes, leakage,
                          leakageVariance, detector, boundaryByXj, boundaryFixed, simPhotons, responseGridPhoton,
                          recoEdges, truthEdges, lineage] = w;
  for (auto *m : {&pairResponse, &combinatoric, &sidebandJets, &photonResponse, &leakage, &boundaryByXj, &responseGridPhoton})
    visit(m->GetMatrixArray(), m->GetNoElements(), false);
  for (auto *v : {&pairMisses, &pairTruth, &pairReco, &pairFakes, &photonMisses, &photonTruth, &photonReco, &photonFakes,
                  &detector, &boundaryFixed, &simPhotons})
    visit(v->GetMatrixArray(), v->GetNoElements(), false);
  for (auto *m : {&combinatoricVariance, &leakageVariance}) visit(m->GetMatrixArray(), m->GetNoElements(), true);
}

/* A window as one row of numbers; `counts` receives how many of them are counts. */
std::vector<double> flatten(const Window &w, int *counts = nullptr) {
  std::vector<double> row;
  int n = 0;
  blocks(w, [&](const double *x, int size, bool variance) {
    row.insert(row.end(), x, x + size);
    if (!variance) n += size;
  });
  if (counts) *counts = n;
  return row;
}

/* The rows of every source of the result's samples, read from the packages (records/sources and
   records/source_contributions). The sample's cross section per generated event stays a separate factor; combine
   applies it once. */
Sources sourcesFromPackages(const Settings &s, const SystemSettings &sys, const Grid &grid, const Window &layout) {
  const std::string base = grid.responsePrefix + "response/" + grid.centrality + "/" + grid.jet + "/";
  std::vector<std::string> paths;
  for (const char *n : {"xj_reco", "xj_truth", "xj_response_truth_global_x_reco_global", "xj_misses", "xj_detector_fakes_reco",
                        "xj_combinatoric_reco", "photon_response_truth_x_reco", "photon_misses", "photon_truth", "photon_reco",
                        "signal_leakage_truth_matched_abcd", "dphi_truthjet_near", "dphi_truthjet_far", "dphi_combinatoric_near",
                        "dphi_combinatoric_far"})
    paths.push_back(base + n);
  Sources out{layout};
  for (size_t k = 0; k < sys.simulation.size(); ++k) {
    auto file = open(s.inputDirectory + "/" + sys.package + "/" + sys.simulation[k] + ".root");
    named_response::requirePinned(s,sys,sys.simulation[k],*file);
    auto *weight = dynamic_cast<TParameter<double> *>(file->Get("meta/pb_per_generated_event"));
    if (!weight) throw std::runtime_error(std::string(file->GetName()) + ": no meta/pb_per_generated_event");
    auto *sources = dynamic_cast<TTree *>(file->Get("records/sources"));
    auto *cells = dynamic_cast<TTree *>(file->Get("records/source_contributions"));
    if (!sources || !cells) throw std::runtime_error(std::string(file->GetName()) + ": no per-source records for the source bootstrap");
    SourceHistograms h(*file, paths, base + "xj_response_truth_global_x_reco_global");
    TTreeReader rs(sources), rc(cells);
    TTreeReaderValue<Long64_t> begin(rs, "contribution_begin"), end(rs, "contribution_end");
    TTreeReaderValue<std::string> json(rs, "json"), histogram(rc, "histogram");
    TTreeReaderValue<Long64_t> bin(rc, "bin");
    TTreeReaderValue<double> value(rc, "value"), variance(rc, "variance");
    while (rs.Next()) {
      h.clear();
      for (Long64_t i = *begin; i < *end; ++i) {
        if (rc.SetEntry(i) != TTreeReader::kEntryValid) throw std::runtime_error("cannot read source contribution " + std::to_string(i));
        h.set(*histogram, static_cast<int>(*bin), *value, *variance);
      }
      const std::vector<double> row = flatten(windowOf(h, grid), &out.counts);
      out.begin.push_back(out.column.size());
      for (size_t c = 0; c < row.size(); ++c)
        if (row[c] != 0) out.column.push_back(c), out.value.push_back(row[c]);
      out.sample.push_back(k);
      out.weight.push_back(weight->GetVal());
      out.identity.push_back(nlohmann::json::parse(*json).at("sha256"));
    }
  }
  if (out.size() == 0) throw std::runtime_error(sys.name + ": no simulated sources");
  if (out.column.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error(sys.name + ": more source entries than an int counts");
  out.begin.push_back(out.column.size());
  return out;
}

/* The saved rows, if the file was made from exactly this description. */
bool readSaved(const std::string &path, const std::string &description, Sources &out) {
  if (!std::filesystem::is_regular_file(path)) return false;
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) return false;
  auto *made = dynamic_cast<TObjString *>(file->Get("description"));
  auto *counts = dynamic_cast<TParameter<int> *>(file->Get("counts"));
  if (!made || !counts || description != made->GetString().Data()) return false;
  auto read = [&](const char *name, auto &values) {
    std::remove_reference_t<decltype(values)> *saved = nullptr;
    file->GetObject(name, saved);  /* a copy this function owns */
    if (!saved) return false;
    values = std::move(*saved);
    delete saved;
    return true;
  };
  out.counts = counts->GetVal();
  if (!(read("sample", out.sample) && read("weight", out.weight) && read("identity", out.identity) &&
        read("begin", out.begin) && read("column", out.column) && read("value", out.value)))
    return false;
  /* the vectors describe one table: a begin per source and one more, a value per column, every column inside the layout */
  const size_t n = out.sample.size(), width = flatten(out.layout).size();
  if (out.weight.size() != n || out.identity.size() != n || out.begin.size() != n + 1 || out.column.size() != out.value.size() ||
      out.begin.front() != 0 || out.begin.back() != static_cast<int>(out.column.size()) ||
      !std::is_sorted(out.begin.begin(), out.begin.end()))
    return false;
  const bool valid=std::all_of(out.column.begin(), out.column.end(), [&](int c) { return c >= 0 && c < static_cast<int>(width); });
  if(valid) observeNativeInput(*file);
  return valid;
}

void save(const Sources &t, const std::string &path, const std::string &description) {
  std::filesystem::create_directories(std::filesystem::path(path).parent_path());
  const std::string partial = path + ".partial" + std::to_string(getpid());  /* a result runs in its own process */
  {
    std::unique_ptr<TFile> file(TFile::Open(partial.c_str(), "RECREATE"));
    if (!file || file->IsZombie()) throw std::runtime_error("cannot write " + partial);
    TObjString(description.c_str()).Write("description");
    TParameter<int>("counts", t.counts).Write();
    file->WriteObject(&t.sample, "sample");
    file->WriteObject(&t.weight, "weight");
    file->WriteObject(&t.identity, "identity");
    file->WriteObject(&t.begin, "begin");
    file->WriteObject(&t.column, "column");
    file->WriteObject(&t.value, "value");
    writeNativeProducerWitness(*file);
  }
  std::filesystem::rename(partial, path);
}

/* The largest difference between two windows, each number relative to the largest of its block in either window. A
   block of counts that is zero up to rounding in both (below 1e-12 of the largest count; boundaryFixed is one by
   construction) has nothing to compare. A number that is not finite, or windows of other sizes, is a full mismatch. */
double mismatch(const Window &a, const Window &b) {
  const std::vector<double> x = flatten(a), y = flatten(b);
  const double everything = std::numeric_limits<double>::infinity();
  if (x.size() != y.size()) return everything;
  double largest = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return everything;
    largest = std::max({largest, std::fabs(x[i]), std::fabs(y[i])});
  }
  double worst = 0;
  size_t at = 0;
  blocks(a, [&](const double *, int size, bool variance) {
    double scale = 0, difference = 0;
    for (int i = 0; i < size; ++i, ++at) {
      scale = std::max({scale, std::fabs(x[at]), std::fabs(y[at])});
      difference = std::max(difference, std::fabs(x[at] - y[at]));
    }
    if (scale > (variance ? 0.0 : 1e-12 * largest)) worst = std::max(worst, difference / scale);
  });
  return worst;
}

}  /* namespace */

Window readWindow(const Settings &s, const SystemSettings &sys) {
  Simulation f(s, sys);
  return windowOf(f, Grid(s, sys));
}

SimulatedDphi readSimulatedDphi(const Settings &s, const SystemSettings &sys) {
  Simulation f(s, sys);
  const std::string base = named_response::prefix(sys) + "response/" + sys.centrality + "/" + s.jet + "/dphi_";
  const Edges &dphi = dphiEdges();
  const int n = dphi.size() - 1;
  SimulatedDphi out;
  for (auto *v : {&out.truthJet, &out.truthJetVariance, &out.combinatoric, &out.combinatoricVariance}) v->ResizeTo(n);
  /* a (photon E_T, xJ, Delta phi) histogram summed over the photon and xJ bins of the result */
  auto add = [&](const std::string &name, TVectorD &values, TVectorD &variance) {
    auto *h = f.get<TH3>(base + name);
    out.lineage[name] = science::projectionOf(h, {{{"axis",0},{"low",sys.photonEdges.front()},{"high",sys.photonEdges.back()}},
      {{"axis",1},{"low",s.xjEdges.front()},{"high",s.xjEdges.back()}}})
      .then("sum_photon_xj_bins",{{"dphi_edges",dphi}});
    if (h->GetNbinsZ() != n) throw std::runtime_error(base + name + ": other Delta-phi bins than the analysis");
    for (int z = 1; z <= n; ++z) {
      if (std::fabs(h->GetZaxis()->GetBinLowEdge(z) - dphi[z - 1]) > 1e-9) throw std::runtime_error(base + name + ": other Delta-phi edges");
      for (int x = 1; x <= h->GetNbinsX(); ++x) {
        if (h->GetXaxis()->GetBinLowEdge(x) < sys.photonEdges.front() - 1e-9 || h->GetXaxis()->GetBinUpEdge(x) > sys.photonEdges.back() + 1e-9) continue;
        for (int y = 1; y <= h->GetNbinsY(); ++y) {
          if (h->GetYaxis()->GetBinLowEdge(y) < s.xjEdges.front() - 1e-9 || h->GetYaxis()->GetBinUpEdge(y) > s.xjEdges.back() + 1e-9) continue;
          values[z - 1] += h->GetBinContent(x, y, z);
          variance[z - 1] += h->GetBinError(x, y, z) * h->GetBinError(x, y, z);
        }
      }
    }
  };
  for (const char *side : {"near", "far"}) {
    add(std::string("truthjet_") + side, out.truthJet, out.truthJetVariance);
    add(std::string("combinatoric_") + side, out.combinatoric, out.combinatoricVariance);
  }
  out.lineage["truthJet"] = science::Trace::combine("sum_near_far",{out.lineage.at("truthjet_near"),out.lineage.at("truthjet_far")});
  out.lineage["combinatoric"] = science::Trace::combine("sum_near_far",{out.lineage.at("combinatoric_near"),out.lineage.at("combinatoric_far")});
  return out;
}

/* The rows are a function of the packages, the grid and the code that makes them, and of nothing else. They are saved
   under that description in $FINALANALYSIS_CACHE (final_analysis.sh sets it to build/cache) and read back while the
   description is the same; a package is named by the SHA-256 of its bytes.
   Read or new, the rows must add up to the result's window in every number: saved rows that do not are made again,
   new rows that do not are an error in the packages. Deleting the folder costs only the time to read the packages. */
Sources readSources(const Settings &s, const SystemSettings &sys, const Window &window) {
  const Grid grid(s, sys);
  std::string description = "source rows " FINALANALYSIS_INPUT_CODE " | " + grid.text() + " |";
  for (const auto &name : sys.simulation) {
    description += " " + name + " " + photonjet_render::sha256(s.inputDirectory + "/" + sys.package + "/" + name + ".root", true);
  }
  const char *folder = std::getenv("FINALANALYSIS_CACHE");
  const bool saving = folder && *folder && *FINALANALYSIS_INPUT_CODE;
  TMD5 name;
  name.Update(reinterpret_cast<const UChar_t *>(description.data()), description.size());
  name.Final();
  const std::string path = saving ? std::string(folder) + "/sources_" + name.AsString() + ".root" : "";
  Window layout(window);  /* its dimensions, emptied (a ROOT matrix is sized by construction, never by assignment) */
  blocks(layout, [](double *x, int size, bool) { std::fill(x, x + size, 0.0); });
  auto addsUp = [&](const Sources &rows) {
    return mismatch(combine(rows, std::vector<double>(rows.size(), 1.0)), window) <= 1e-9;
  };
  if (saving) {
    Sources saved{layout};
    if (readSaved(path, description, saved) && addsUp(saved)) return saved;
  }
  Sources rows = sourcesFromPackages(s, sys, grid, layout);
  if (!addsUp(rows)) throw std::runtime_error(sys.name + ": the per-source contributions do not add up to the package response");
  if (saving) save(rows, path, description);
  return rows;
}

Window combine(const Sources &sources, const std::vector<double> &multiplicity) {
  std::vector<double> row(flatten(sources.layout).size());
  for (size_t i = 0; i < sources.size(); ++i) {
    const double c = multiplicity[i] * sources.weight[i], c2 = c * sources.weight[i];
    if (c == 0) continue;
    for (int j = sources.begin[i]; j < sources.begin[i + 1]; ++j)
      row[sources.column[j]] += (sources.column[j] < sources.counts ? c : c2) * sources.value[j];
  }
  Window w = sources.layout;
  const double *next = row.data();
  blocks(w, [&](double *x, int size, bool) { std::copy(next, next + size, x); next += size; });
  return w;
}

namespace {

template <class Provider> Window windowOf(Provider &f, const Grid &g) {
  const std::string base = g.responsePrefix + "response/" + g.centrality + "/" + g.jet + "/";
  auto *pairReco = f.template get<TH2>(base + "xj_reco");
  auto *pairTruth = f.template get<TH2>(base + "xj_truth");
  const Edges fineReco = axisEdges(pairReco->GetXaxis()), fineTruth = axisEdges(pairTruth->GetXaxis());
  const Edges fineXj = axisEdges(pairReco->GetYaxis());
  const TMatrixD rp = toMatrix(projection(fineReco, g.photon));
  const TMatrixD tp = toMatrix(projection(fineTruth, g.photon));
  const auto pXj = projection(fineXj, g.xj);
  const TMatrixD rpair = cellProjection(projection(fineReco, g.photon), pXj);
  const TMatrixD tpair = cellProjection(projection(fineTruth, g.photon), pXj);
  const int nPt = g.photon.size() - 1, nX = g.xj.size() - 1;

  /* The joint response as its filled cells. The projections hold only 0 and 1, so a fine cell lies in at most one
     window cell and every product below is one sum over the filled cells. */
  const std::vector<ResponseCell> M = f.responseCells(base + "xj_response_truth_global_x_reco_global");
  const TVectorD misses = cellsOf(f.template get<TH2>(base + "xj_misses"));
  const TVectorD rOut = column(rpair), tOut = column(tpair);
  auto owner = [](const TMatrixD &p) {
    std::vector<int> cell(p.GetNcols(), -1);
    for (int c = 0; c < p.GetNrows(); ++c)
      for (int j = 0; j < p.GetNcols(); ++j)
        if (p(c, j) != 0) {
          if (p(c, j) != 1 || cell[j] >= 0) throw std::runtime_error("a fine cell lies in more than one window cell");
          cell[j] = c;
        }
    return cell;
  };
  const std::vector<int> tCell = owner(tpair), rCell = owner(rpair);
  const int nFineX = fineXj.size() - 1;
  Window w;
  w.recoEdges = g.reco;
  w.truthEdges = g.truth;
  w.pairResponse.ResizeTo(nPt * nX, nPt * nX);
  TVectorD recoOutside(tpair.GetNcols());   /* per fine truth cell: matched pairs whose reco cell lies outside the window */
  TVectorD truthOutside(rpair.GetNcols());  /* per fine reco cell: matched pairs whose truth lies outside the window */
  TMatrixD truthOutsideByXj(nX, rpair.GetNcols());
  for (const ResponseCell &c : M) {
    if (c.truth >= tpair.GetNcols() || c.reco >= rpair.GetNcols()) throw std::runtime_error("a response cell lies outside the package grid");
    if (tCell[c.truth] >= 0 && rCell[c.reco] >= 0) w.pairResponse(tCell[c.truth], rCell[c.reco]) += c.value;
    recoOutside[c.truth] += c.value * rOut[c.reco];
    if (tOut[c.truth] > 0) {
      truthOutside[c.reco] += tOut[c.truth] * c.value;
      for (int x = 0; x < nX; ++x)
        if (pXj[x][c.truth % nFineX] > 0) truthOutsideByXj(x, c.reco) += tOut[c.truth] * c.value;
    }
  }
  w.pairMisses.ResizeTo(nPt * nX);
  w.pairMisses = tpair * (misses + recoOutside);
  TVectorD boundary = rpair * truthOutside;
  w.detector.ResizeTo(nPt * nX);
  w.detector = rpair * cellsOf(f.template get<TH2>(base + "xj_detector_fakes_reco"));
  w.pairFakes.ResizeTo(nPt * nX);
  w.pairFakes = boundary + w.detector;
  /* the boundary migration per truth xJ bin: it scales with the prior weight of that bin */
  w.boundaryByXj.ResizeTo(nX, nPt * nX);
  for (int x = 0; x < nX; ++x) {
    TVectorD part(rpair.GetNcols());
    for (int j = 0; j < rpair.GetNcols(); ++j) part[j] = truthOutsideByXj(x, j);
    const TVectorD projected = rpair * part;
    for (int j = 0; j < nPt * nX; ++j) w.boundaryByXj(x, j) = projected[j];
  }
  w.boundaryFixed.ResizeTo(nPt * nX);
  w.boundaryFixed = boundary;
  for (int x = 0; x < nX; ++x)
    for (int j = 0; j < nPt * nX; ++j) w.boundaryFixed[j] -= w.boundaryByXj(x, j);
  w.pairTruth.ResizeTo(nPt * nX);
  w.pairTruth = tpair * cellsOf(pairTruth);
  w.pairReco.ResizeTo(nPt * nX);
  w.pairReco = rpair * cellsOf(pairReco);
  auto *comb = f.template get<TH2>(base + "xj_combinatoric_reco");
  w.combinatoric.ResizeTo(nPt, nX);
  w.combinatoric = reshape(rpair * cellsOf(comb), nPt, nX);
  w.combinatoricVariance.ResizeTo(nPt, nX);
  w.combinatoricVariance = reshape(rpair * cellsOf(comb, true), nPt, nX);
  /* every jet of a matched region-A photon in the Delta-phi sideband: the denominator of the simulated transfer */
  TVectorD sideband(rpair.GetNcols());
  for (const char *n : {"dphi_truthjet_near", "dphi_truthjet_far", "dphi_combinatoric_near", "dphi_combinatoric_far"})
    sideband += sidebandCells(f.template get<TH3>(base + n), g.sidebandLow, g.sidebandHigh);
  w.sidebandJets.ResizeTo(nPt, nX);
  w.sidebandJets = reshape(rpair * sideband, nPt, nX);

  const TMatrixD PM = values2(f.template get<TH2>(base + "photon_response_truth_x_reco"));  /* truth x reco, fine */
  const TVectorD photonMisses = values1(f.template get<TH1>(base + "photon_misses"));
  TVectorD rpOut(rp.GetNcols()), tpOut(tp.GetNcols());
  for (int j = 0; j < rp.GetNcols(); ++j) { double c = 0; for (int i = 0; i < nPt; ++i) c += rp(i, j); rpOut[j] = 1 - c; }
  for (int j = 0; j < tp.GetNcols(); ++j) { double c = 0; for (int i = 0; i < nPt; ++i) c += tp(i, j); tpOut[j] = 1 - c; }
  w.photonResponse.ResizeTo(nPt, nPt);
  w.photonResponse = tp * PM * TMatrixD(TMatrixD::kTransposed, rp);
  w.photonMisses.ResizeTo(nPt);
  w.photonMisses = tp * (photonMisses + PM * rpOut);
  TVectorD photonTruthOutside(PM.GetNcols());
  for (int j = 0; j < PM.GetNcols(); ++j)
    for (int i = 0; i < PM.GetNrows(); ++i) photonTruthOutside[j] += tpOut[i] * PM(i, j);
  w.photonFakes.ResizeTo(nPt);
  w.photonFakes = rp * photonTruthOutside;
  w.photonTruth.ResizeTo(nPt);
  w.photonTruth = tp * values1(f.template get<TH1>(base + "photon_truth"));
  w.photonReco.ResizeTo(nPt);
  w.photonReco = rp * values1(f.template get<TH1>(base + "photon_reco"));
  w.simPhotons.ResizeTo(nPt);
  for (int j = 0; j < nPt; ++j) {
    for (int i = 0; i < nPt; ++i) w.simPhotons[j] += w.photonResponse(i, j);
    w.simPhotons[j] += w.photonFakes[j];
  }

  auto *leak = f.template get<TH2>(base + "signal_leakage_truth_matched_abcd");  /* 4 x fine reco */
  const TMatrixD rpT(TMatrixD::kTransposed, rp);
  w.leakage.ResizeTo(4, nPt);
  w.leakage = values2(leak) * rpT;
  w.leakageVariance.ResizeTo(4, nPt);
  w.leakageVariance = values2(leak, true) * rpT;

  /* the photon response on the configured response grid, for display */
  const TMatrixD gridT = toMatrix(projection(fineTruth, g.truth)), gridR = toMatrix(projection(fineReco, g.reco));
  w.responseGridPhoton.ResizeTo(gridT.GetNrows(), gridR.GetNrows());
  w.responseGridPhoton = gridT * PM * TMatrixD(TMatrixD::kTransposed, gridR);
  /* Follow the exact additive members. Photon-only quantities do not inherit
     joint-pair binning, sideband settings or prior/unfolding choices. */
  const science::Json pairGrid={{"photon_edges",g.photon},{"xj_edges",g.xj}};
  auto pair = [&](const std::string &name) { return withSettings(withSettings(f.trace(base+name),g.system,{"Centrality","Simulation","PhotonEdges"}),g.settings,{"Measurement.XjEdges","Recoil.Jet"}).then("restrict_pair_window",pairGrid); };
  auto photon = [&](const std::string &name) { return withSettings(f.trace(base+name),g.system,{"Centrality","Simulation","PhotonEdges"}).then("restrict_photon_window",{{"photon_edges",g.photon}}); };
  auto join = [&](const std::string &operation,std::initializer_list<science::Trace> inputs) {
    return science::Trace::combine(operation,std::vector<science::Trace>(inputs));
  };
  const auto response=pair("xj_response_truth_global_x_reco_global");
  w.lineage["pairResponse"]=response;
  w.lineage["pairMisses"]=join("add_reco_outside_misses",{pair("xj_misses"),response});
  w.lineage["detector"]=pair("xj_detector_fakes_reco");
  w.lineage["pairFakes"]=join("add_truth_outside_fakes",{response,w.lineage.at("detector")});
  w.lineage["boundaryByXj"]=response.then("truth_outside_by_xj");
  w.lineage["boundaryFixed"]=response.then("truth_outside_without_xj");
  w.lineage["pairTruth"]=pair("xj_truth");
  w.lineage["pairReco"]=pair("xj_reco");
  w.lineage["combinatoric"]=pair("xj_combinatoric_reco");
  std::vector<science::Trace> sidebandInputs;
  for (const char *name : {"dphi_truthjet_near","dphi_truthjet_far","dphi_combinatoric_near","dphi_combinatoric_far"})
    sidebandInputs.push_back(pair(name));
  w.lineage["sidebandJets"]=science::Trace::combine("sum_dphi_sideband",sidebandInputs,
    {{"low",g.sidebandLow},{"high",g.sidebandHigh}});
  w.lineage["sidebandJets"]=withSettings(w.lineage["sidebandJets"],g.settings,{"Combinatoric.Sideband"});
  const auto photonResponse=photon("photon_response_truth_x_reco");
  w.lineage["photonResponse"]=photonResponse;
  w.lineage["photonMisses"]=join("add_reco_outside_misses",{photon("photon_misses"),photonResponse});
  w.lineage["photonFakes"]=photonResponse.then("truth_outside_fakes");
  w.lineage["photonTruth"]=photon("photon_truth");
  w.lineage["photonReco"]=photon("photon_reco");
  w.lineage["simPhotons"]=join("column_sum_plus_fakes",{photonResponse,w.lineage.at("photonFakes")});
  w.lineage["leakage"]=photon("signal_leakage_truth_matched_abcd");
  w.lineage["responseGridPhoton"]=f.trace(base+"photon_response_truth_x_reco").then("response_display_grid",{{"truth_edges",g.truth},{"reco_edges",g.reco}});
  w.lineage["responseGridPhoton"]=withSettings(w.lineage["responseGridPhoton"],g.system,{"Centrality","Simulation","RecoEdges","TruthEdges"});
  return w;
}

}  /* namespace */

}  /* namespace fa */
