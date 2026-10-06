/* Inputs: the settings file, the measured histograms and the simulated response of one centrality class. */

/* A histogram package keeps every quantity on its own fine bins. Coarser analysis bins are exact unions of those bins
   (projection); a boundary that the package does not have is an error, never an interpolation. */
#include "FinalAnalysis.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "TEnv.h"
#include "TFile.h"
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

namespace fa {
std::string inputCoverage(const Settings &s) {
  using Json = nlohmann::json;
  Json coverage = Json::object(); std::set<std::string> seen;
  std::map<std::string, std::string> selections; bool partial = false, unknown = false, legacy = false;
  for (const auto &sys : s.systems) {
    std::vector<std::string> names{"data"}; names.insert(names.end(), sys.simulation.begin(), sys.simulation.end());
    const auto inclusive = s.inclusive.find(sys.package);
    if (inclusive != s.inclusive.end()) names.insert(names.end(), inclusive->second.begin(), inclusive->second.end());
    for (const auto &name : names) {
      const std::string key = sys.package + "/" + name;
      if (!seen.insert(key).second) continue;
      const std::string path = s.inputDirectory + "/" + key + ".root";
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
      auto *saved = dynamic_cast<TObjString *>(file->Get("meta/package"));
      if (!saved) throw std::runtime_error(path + ": missing sample-package metadata");
      const Json meta = Json::parse(saved->GetString().Data());
      if (meta.value("schema", "") != "PhotonJetSampleV1" || meta.value("system", "") != sys.package || meta.value("name", "") != name)
        throw std::runtime_error(path + ": package identity differs from its analysis binding");
      if (meta.value("legacy", false)) { legacy = true; coverage[key] = {{"status", "LEGACY_REFERENCE"}}; continue; }
      const std::string selection = meta.at("selection").at("selection_sha256");
      if (selection.empty() || (selections.count(sys.package) && selections.at(sys.package) != selection))
        throw std::runtime_error(path + ": DATA/SIM selection hashes differ");
      selections[sys.package] = selection;
      if (name != "data") {
        auto *weight = dynamic_cast<TParameter<double> *>(file->Get("meta/pb_per_generated_event"));
        const double recorded = meta.at("normalization").at("pb_per_generated_event");
        if (!weight || !(recorded > 0) || !std::isfinite(recorded) || weight->GetVal() != recorded)
          throw std::runtime_error(path + ": SIM normalization metadata differs from its applied weight");
      }
      coverage[key] = meta.value("coverage", Json{{"status", "UNKNOWN"}});
      const std::string status = coverage[key].value("status", "UNKNOWN");
      partial = partial || status == "PARTIAL"; unknown = unknown || status != "PARTIAL" && status != "COMPLETE";
    }
  }
  return Json{{"status", legacy ? "REFERENCE_INPUTS" : partial ? "PARTIAL" : unknown ? "UNKNOWN" : "COMPLETE"}, {"samples", coverage}}.dump();
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

}  /* namespace */

Settings Settings::read(const std::string &path) {
  TEnv env;
  if (env.ReadFile(path.c_str(), kEnvUser) != 0) throw std::runtime_error("cannot read settings " + path);
  Settings s;
  s.name = text(env, "Analysis.Name");
  s.inputDirectory = text(env, "Input.Directory");
  s.output = text(env, "Output.Results");
  s.qaOutput = text(env, "Output.QA");
  s.jet = text(env, "Recoil.Jet");
  s.jetThreshold = slotThreshold(s.jet);
  s.xjEdges = increasing(env, "Measurement.XjEdges");
  s.reportedHigh = number(env, "Measurement.ReportedXjMax");
  if (!isEdge(s.xjEdges, s.reportedHigh)) throw std::runtime_error("Measurement.ReportedXjMax must be an xJ bin edge");
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
  if (s.ruleDataToys < 2 || s.ruleResponseToys < 2) throw std::runtime_error("the iteration rule needs at least two data and two response toys");
  if (s.responseToys < 0 || s.responseToys == 1) throw std::runtime_error("Uncertainty.ResponseToys is 0 or at least 2");
  s.combinatoricMethod = choice(env, "Combinatoric.Method", {"simulation", "data_dphi"});
  s.combinatoricShape = choice(env, "Combinatoric.Shape", {"flow", "flat"});
  s.combinatoricScope = choice(env, "Combinatoric.Scope", {"pt", "pt_xj"});
  if (s.combinatoricMethod == "data_dphi" && s.dataToyModel != "event_bootstrap")
    throw std::runtime_error("Combinatoric.Method data_dphi is refitted per event: it needs Uncertainty.DataToyModel event_bootstrap");
  s.scanIterations = integer(env, "Unfolding.ScanIterations");
  s.priorForm = choice(env, "Prior.Form", {"simulation", "exp1", "exp2"});
  s.fakeModel = choice(env, "Prior.FakeModel", {"scaled", "fixed"});
  s.fitRange = choice(env, "Prior.FitRange", {"all", "reported"});
  s.priorPower = number(env, "Prior.Power");
  if (s.dataToys < 2 || s.scanIterations < 1 || !(s.priorPower > 0)) throw std::runtime_error("DataToys >= 2, ScanIterations >= 1, Prior.Power > 0");
  std::set<std::string> known = {"Analysis.Name", "Input.Directory", "Output.Results", "Output.QA", "Recoil.Jet", "Measurement.XjEdges",
                                 "Measurement.ReportedXjMax", "Uncertainty.DataToys", "Uncertainty.Seed",
                                 "Uncertainty.DataToyModel", "Uncertainty.ResponseToys", "Uncertainty.ResponseSeed", "Closure.Toys", "Closure.Seed",
                                 "IterationRule.DataToys", "IterationRule.DataSeed", "IterationRule.ResponseToys", "IterationRule.ResponseSeed",
                                 "Combinatoric.Method", "Combinatoric.Shape", "Combinatoric.Scope",
                                 "Unfolding.ScanIterations", "Prior.Form", "Prior.FakeModel", "Prior.FitRange", "Prior.Power",
                                 "Systems"};
  std::istringstream names(text(env, "Systems"));
  for (std::string name; names >> name;) {
    for (const char *key : {"Package", "Centrality", "Label", "Simulation", "PhotonEdges", "RecoEdges", "TruthEdges", "PairIterations",
                            "PhotonIterations", "BackgroundCorrelation", "PurityScale", "EfficiencyScale", "SubtractCombinatoric"})
      known.insert(name + "." + key);
    SystemSettings sys;
    sys.name = name;
    sys.package = text(env, name + ".Package");
    sys.centrality = text(env, name + ".Centrality");
    sys.label = text(env, name + ".Label");
    sys.photonEdges = increasing(env, name + ".PhotonEdges");
    sys.recoEdges = increasing(env, name + ".RecoEdges");
    sys.truthEdges = increasing(env, name + ".TruthEdges");
    std::istringstream samples(text(env, name + ".Simulation"));
    for (std::string sample; samples >> sample;) sys.simulation.push_back(sample);
    /* a number, or "rule": the frozen iteration rule selects k */
    sys.pairIterations = text(env, name + ".PairIterations") == "rule" ? 0 : integer(env, name + ".PairIterations");
    sys.photonIterations = text(env, name + ".PhotonIterations") == "rule" ? 0 : integer(env, name + ".PhotonIterations");
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
    if (s.combinatoricMethod == "data_dphi" && !sys.subtractCombinatoric)
      throw std::runtime_error(name + ": the data-driven combinatoric estimate is always subtracted");
    sys.slices = slicesOf(sys.centrality);
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
  /* a misspelt or leftover setting is an error, never silently ignored */
  TIter next(env.GetTable());
  while (auto *record = static_cast<TEnvRec *>(next()))
    if (!known.count(record->GetName())) throw std::runtime_error(path + ": unknown setting " + record->GetName());
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
  t.nPt = sys.photonEdges.size() - 1;
  t.nXj = s.xjEdges.size() - 1;
  const int nDphi = dphiEdges().size() - 1;
  std::vector<int> renumber;
  readRows(*f, "records/events", {"slice", "weight"}, [&](const std::vector<double> &v) {
    const bool selected = std::find(sys.slices.begin(), sys.slices.end(), static_cast<int>(v[0])) != sys.slices.end();
    renumber.push_back(selected ? static_cast<int>(t.weight.size()) : -1);
    if (selected) t.weight.push_back(std::isfinite(v[1]) ? v[1] : 1.0);
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
  for (size_t i = 0; i < t.photonEvent.size(); ++i) m.abcd(t.photonCell[i] / t.nPt, t.photonCell[i] % t.nPt) += weight(t.photonEvent[i]);
  for (size_t i = 0; i < t.jetEvent.size(); ++i) {
    const int region = t.jetCell[i] / (t.nPt * t.nXj), rest = t.jetCell[i] % (t.nPt * t.nXj);
    (region == 0 ? m.pairsA : m.pairsC)(rest / t.nXj, rest % t.nXj) += weight(t.jetEvent[i]);
  }
  m.abcdVariance = m.abcd;
  m.pairsAVariance = m.pairsA;
  m.pairsCVariance = m.pairsC;
  if (s.combinatoricMethod == "data_dphi") {
    const int perRegion = t.nPt * t.nXj * (dphiEdges().size() - 1) * 2;
    std::vector<double> dphi[2] = {std::vector<double>(perRegion), std::vector<double>(perRegion)};
    for (size_t i = 0; i < t.dphiEvent.size(); ++i) dphi[t.dphiCell[i] / perRegion][t.dphiCell[i] % perRegion] += weight(t.dphiEvent[i]);
    const Background b = combinatoricFromDphi(dphi, dphi, t.nPt, t.nXj, s.combinatoricShape, s.combinatoricScope);
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
class Simulation {
 public:
  Simulation(const Settings &s, const SystemSettings &sys) {
    for (const auto &name : sys.simulation) {
      files_.push_back(open(s.inputDirectory + "/" + sys.package + "/" + name + ".root"));
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
      } else if (!sum->Add(h, weights_[k])) {
        throw std::runtime_error("the simulated samples disagree in the binning of " + path);
      }
    }
    cached_.insert(path);
    return static_cast<T *>(sum.get());
  }

 private:
  std::vector<std::unique_ptr<TFile>> files_;
  std::vector<double> weights_;
  std::map<std::string, std::unique_ptr<TH1>> cache_;
  std::set<std::string> cached_;
};

/* One source's response histograms rebuilt from its sparse cells (records/source_contributions of its package). */
class SourceHistograms {
 public:
  SourceHistograms(TFile &package, const std::vector<std::string> &paths) {
    for (const auto &path : paths) {
      auto *h = static_cast<TH1 *>(fa::get<TH1>(package, path)->Clone());
      h->SetDirectory(nullptr);
      h->Reset();
      hists_[path].reset(h);
    }
  }
  template <class T> T *get(const std::string &path) {
    auto it = hists_.find(path);
    if (it == hists_.end()) throw std::runtime_error("no per-source histogram " + path);
    return static_cast<T *>(it->second.get());
  }
  void clear() { for (auto &item : hists_) item.second->Reset(); }
  void set(const std::string &path, int bin, double value, double variance) {
    auto it = hists_.find(path);
    if (it == hists_.end()) return;  /* a cell of another histogram */
    it->second->SetBinContent(bin, value);
    it->second->GetSumw2()->SetAt(variance, bin);
  }

 private:
  std::map<std::string, std::unique_ptr<TH1>> hists_;
};

template <class Provider> Window windowOf(Provider &f, const Settings &s, const SystemSettings &sys);

}  /* namespace */

Window readWindow(const Settings &s, const SystemSettings &sys) {
  Simulation f(s, sys);
  return windowOf(f, s, sys);
}

/* Every source of the result's simulated samples with its own window (records/sources, records/source_contributions)
   and its sample's cross section per generated event. */
std::vector<SourceWindow> readSourceWindows(const Settings &s, const SystemSettings &sys) {
  const std::string base = "response/" + sys.centrality + "/" + s.jet + "/";
  std::vector<std::string> paths;
  for (const char *n : {"xj_reco", "xj_truth", "xj_response_truth_global_x_reco_global", "xj_misses", "xj_detector_fakes_reco",
                        "xj_combinatoric_reco", "photon_response_truth_x_reco", "photon_misses", "photon_truth", "photon_reco",
                        "signal_leakage_truth_matched_abcd"})
    paths.push_back(base + n);
  std::vector<SourceWindow> out;
  for (size_t k = 0; k < sys.simulation.size(); ++k) {
    auto file = open(s.inputDirectory + "/" + sys.package + "/" + sys.simulation[k] + ".root");
    auto *weight = dynamic_cast<TParameter<double> *>(file->Get("meta/pb_per_generated_event"));
    if (!weight) throw std::runtime_error(std::string(file->GetName()) + ": no meta/pb_per_generated_event");
    auto *sources = dynamic_cast<TTree *>(file->Get("records/sources"));
    auto *cells = dynamic_cast<TTree *>(file->Get("records/source_contributions"));
    if (!sources || !cells) throw std::runtime_error(std::string(file->GetName()) + ": no per-source records for the source bootstrap");
    SourceHistograms h(*file, paths);
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
      out.push_back(SourceWindow{windowOf(h, s, sys), weight->GetVal(), static_cast<int>(k), nlohmann::json::parse(*json).at("sha256")});
    }
  }
  if (out.empty()) throw std::runtime_error(sys.name + ": no simulated sources");
  return out;
}

/* The window of weighted source occurrences: counts times multiplicity x weight, variances times multiplicity x weight^2. */
Window combine(const std::vector<SourceWindow> &sources, const std::vector<double> &multiplicity) {
  Window w = sources.front().window;
  auto zero = [](auto &x) { x.Zero(); };
  for (auto *m : {&w.pairResponse, &w.combinatoric, &w.combinatoricVariance, &w.photonResponse, &w.leakage, &w.leakageVariance,
                  &w.boundaryByXj, &w.responseGridPhoton}) zero(*m);
  for (auto *v : {&w.pairMisses, &w.pairTruth, &w.pairReco, &w.pairFakes, &w.photonMisses, &w.photonTruth, &w.photonReco,
                  &w.photonFakes, &w.detector, &w.boundaryFixed, &w.simPhotons}) zero(*v);
  for (size_t i = 0; i < sources.size(); ++i) {
    const Window &x = sources[i].window;
    const double c = multiplicity[i] * sources[i].weight, c2 = c * sources[i].weight;
    if (c == 0) continue;
    for (auto [target, source] : {std::pair<TMatrixD *, const TMatrixD *>{&w.pairResponse, &x.pairResponse}, {&w.combinatoric, &x.combinatoric},
                                  {&w.photonResponse, &x.photonResponse}, {&w.leakage, &x.leakage}, {&w.boundaryByXj, &x.boundaryByXj},
                                  {&w.responseGridPhoton, &x.responseGridPhoton}})
      Add(*target, c, *source);
    Add(w.combinatoricVariance, c2, x.combinatoricVariance);
    Add(w.leakageVariance, c2, x.leakageVariance);
    for (auto [target, source] : {std::pair<TVectorD *, const TVectorD *>{&w.pairMisses, &x.pairMisses}, {&w.pairTruth, &x.pairTruth},
                                  {&w.pairReco, &x.pairReco}, {&w.pairFakes, &x.pairFakes}, {&w.photonMisses, &x.photonMisses},
                                  {&w.photonTruth, &x.photonTruth}, {&w.photonReco, &x.photonReco}, {&w.photonFakes, &x.photonFakes},
                                  {&w.detector, &x.detector}, {&w.boundaryFixed, &x.boundaryFixed}, {&w.simPhotons, &x.simPhotons}})
      Add(*target, c, *source);
  }
  return w;
}

namespace {

template <class Provider> Window windowOf(Provider &f, const Settings &s, const SystemSettings &sys) {
  const std::string base = "response/" + sys.centrality + "/" + s.jet + "/";
  auto *pairReco = f.template get<TH2>(base + "xj_reco");
  auto *pairTruth = f.template get<TH2>(base + "xj_truth");
  const Edges fineReco = axisEdges(pairReco->GetXaxis()), fineTruth = axisEdges(pairTruth->GetXaxis());
  const Edges fineXj = axisEdges(pairReco->GetYaxis());
  const TMatrixD rp = toMatrix(projection(fineReco, sys.photonEdges));
  const TMatrixD tp = toMatrix(projection(fineTruth, sys.photonEdges));
  const auto pXj = projection(fineXj, s.xjEdges);
  const TMatrixD rpair = cellProjection(projection(fineReco, sys.photonEdges), pXj);
  const TMatrixD tpair = cellProjection(projection(fineTruth, sys.photonEdges), pXj);
  const int nPt = sys.photonEdges.size() - 1, nX = s.xjEdges.size() - 1;

  const TMatrixD M = values2(f.template get<TH2>(base + "xj_response_truth_global_x_reco_global"));
  const TVectorD misses = cellsOf(f.template get<TH2>(base + "xj_misses"));
  const TVectorD rOut = column(rpair), tOut = column(tpair);
  Window w;
  w.recoEdges = sys.recoEdges;
  w.truthEdges = sys.truthEdges;
  w.pairResponse.ResizeTo(nPt * nX, nPt * nX);
  w.pairResponse = tpair * M * TMatrixD(TMatrixD::kTransposed, rpair);
  w.pairMisses.ResizeTo(nPt * nX);
  w.pairMisses = tpair * (misses + M * rOut);
  TVectorD truthOutside(M.GetNcols());  /* matched reco cells whose truth lies outside the window */
  for (int j = 0; j < M.GetNcols(); ++j)
    for (int i = 0; i < M.GetNrows(); ++i) truthOutside[j] += tOut[i] * M(i, j);
  TVectorD boundary = rpair * truthOutside;
  w.detector.ResizeTo(nPt * nX);
  w.detector = rpair * cellsOf(f.template get<TH2>(base + "xj_detector_fakes_reco"));
  w.pairFakes.ResizeTo(nPt * nX);
  w.pairFakes = boundary + w.detector;
  /* the boundary migration per truth xJ bin: it scales with the prior weight of that bin */
  const int nFineX = fineXj.size() - 1;
  w.boundaryByXj.ResizeTo(nX, nPt * nX);
  for (int x = 0; x < nX; ++x) {
    TVectorD part(M.GetNcols());
    for (int i = 0; i < M.GetNrows(); ++i)
      if (tOut[i] > 0 && pXj[x][i % nFineX] > 0)
        for (int j = 0; j < M.GetNcols(); ++j) part[j] += tOut[i] * M(i, j);
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
  const TMatrixD gridT = toMatrix(projection(fineTruth, sys.truthEdges)), gridR = toMatrix(projection(fineReco, sys.recoEdges));
  w.responseGridPhoton.ResizeTo(gridT.GetNrows(), gridR.GetNrows());
  w.responseGridPhoton = gridT * PM * TMatrixD(TMatrixD::kTransposed, gridR);
  return w;
}

}  /* namespace */

}  /* namespace fa */
