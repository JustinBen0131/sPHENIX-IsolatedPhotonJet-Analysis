/* Isolation: the Au+Au isolation package, the reconstructed-isolation cutoffs of truth-signal photons in embedded
   simulation (config/auau_isolation_wp80.json is one of its freezes). */

/* tree_to_hists isolation population --out population.root EVENTTREE...
   tree_to_hists isolation derive     --out package.json POPULATION...
   tree_to_hists isolation freeze     --out freeze.json PACKAGE [--efficiency 0.8] [--gap 0.8]
   tree_to_hists isolation check      --out check.json FREEZE POPULATION... */

/* population: the truth-signal photons of Au+Au embedded simulation with their reconstructed isolation (below).
   derive: for every photon-pT group and 5% centrality slice (native bin b, slice
   (b - 1) / 5) the reconstructed isolation is histogrammed (700 bins from -20 to 50 GeV, unweighted, finite values) and
   the cutoff where the cumulative spectrum reaches 70, 80 and 90% is read inside its bin by linear interpolation; its
   error is the binomial error of the efficiency over the local density, and at least half a bin. Per slice a constant
   is fitted across the pT groups, per efficiency a straight line across the slice centres,
   threshold(c) = intercept + slope * c        (GeV; c the centrality in percent),
   for the R = 0.3 and 0.4 cones, with the sha256 of every input. freeze: the R = 0.3 line at the nominal efficiency with
   the sideband gap, bound to the package's sha256: the package the eventTree conversion applies to every Au+Au photon.
   check: a frozen line on new simulation (15 <= ET < 35 GeV): its realized efficiency per slice against the target, and
   the line refitted there against the frozen one. An existing output is never replaced. */
#include "TreeToHists.h"

#include "PhotonJet.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <limits>
#include <map>
#include <string>
#include <vector>

/* === The population ==================================================================================================
   The truth-signal photons of Au+Au embedded simulation with their reconstructed isolation: the input of the package. */

/* ./TreeToHists/tree_to_hists.sh isolation population --out population.root a.eventtree.root b.eventtree.root */

/* Events: the corrected minimum-bias decision (cent_mb_decision == 1), a finite vertex with |z| < 10 cm, centrality in
   [0, 80). Photons: matched to an analysis truth-signal photon (pho_truth_match_state == 1), finite E_T and eta,
   |eta| < 0.7. No preselection, identification or E_T cut. One row per photon in the TTree "isolation_population":
   et, eta, centrality, native_bin (cent_native_bin), iso_r03, iso_r04 (the reconstructed isolation witnesses, R = 0.3
   and 0.4; NaN when the witness is not reconstructed), method (their subtraction method), weight (w_event), file
   plus "sources" (as FillPhotonJetResponse) and "counts" (events and photons at each gate, JSON). */
namespace photonjet {
namespace isolation {

/* The event and photon gates of the isolation population. */
constexpr int kMinimumBias = 1;            /* cent_mb_decision of a minimum-bias event */
constexpr double kVertexAbsZMax = 10.0;    /* cm */
constexpr double kCentralityMax = 80.0;    /* percent */
constexpr int kTruthSignalMatch = 1;       /* pho_truth_match_state of a photon matched to an analysis truth-signal photon */
constexpr double kPhotonAbsEtaMax = 0.7;
constexpr int kReconstructed = 0;          /* isow_reconstructed_or_truth of a reconstructed witness */

struct Row {
  double et, eta, centrality, isoR03, isoR04, weight;
  int nativeBin, method, file;
};

inline void fillFile(const std::string& path, int fileIndex, TTree& table, Row& row, std::map<std::string, long>& counts) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  if (!isSimulation(*file) || fileSystem(*file) != System::auau)
    throw std::runtime_error("the isolation population needs Au+Au simulation: " + path);
  TTreeReader reader("eventTree", file.get());
  TTreeReaderValue<int> record(reader, "evt_record"), minimumBias(reader, "cent_mb_decision"),
      nativeBin(reader, "cent_native_bin");
  TTreeReaderValue<double> vz(reader, "vtx_z"), centrality(reader, "cent_value"), weight(reader, "w_event");
  TTreeReaderArray<double> et(reader, "pho_et"), eta(reader, "pho_eta"), cone(reader, "isow_cone_sum");
  TTreeReaderArray<int> match(reader, "pho_truth_match_state"), kind(reader, "isow_reconstructed_or_truth"),
      method(reader, "isow_subtraction_method");
  while (reader.Next()) {
    if (*record != 0) continue;
    ++counts["events"];
    if (*minimumBias != kMinimumBias) continue;
    ++counts["events_minimum_bias"];
    if (!(std::isfinite(*vz) && std::fabs(*vz) < kVertexAbsZMax)) continue;
    ++counts["events_vertex"];
    if (!(std::isfinite(*centrality) && *centrality >= 0.0 && *centrality < kCentralityMax)) continue;
    ++counts["events_selected"];
    const size_t n = et.GetSize();
    if (eta.GetSize() != n || match.GetSize() != n || cone.GetSize() != 2 * n || kind.GetSize() != 2 * n)
      throw std::runtime_error("photon groups or isolation witnesses (2 per photon) disagree in " + path);
    for (size_t i = 0; i < n; ++i) {
      ++counts["photons"];
      if (match[i] != kTruthSignalMatch) continue;
      ++counts["photons_truth_signal"];
      if (!(std::isfinite(et[i]) && std::isfinite(eta[i]) && std::fabs(eta[i]) < kPhotonAbsEtaMax)) continue;
      ++counts["photons_selected"];
      const double nan = std::numeric_limits<double>::quiet_NaN();
      if (method[2 * i] != method[2 * i + 1])
        throw std::runtime_error("the two isolation witnesses of one photon have different methods in " + path);
      row = Row{et[i], eta[i], *centrality, kind[2 * i] == kReconstructed ? cone[2 * i] : nan,
                kind[2 * i + 1] == kReconstructed ? cone[2 * i + 1] : nan, *weight, *nativeBin, method[2 * i], fileIndex};
      table.Fill();
    }
  }
  /* a missing branch or a read error ends the loop like the end of the tree; only the end is accepted */
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
}

inline int run(const char* inputs, const char* output) {
  const std::vector<std::string> paths = inputPaths(inputs);
  const std::string sources = sourceRecords(paths);  /* before the output is opened (sourceRecord) */
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  auto* table = new TTree("isolation_population", "truth-signal photons of Au+Au embedded simulation");  /* owned by out */
  Row row{};
  table->Branch("et", &row.et, "et/D");
  table->Branch("eta", &row.eta, "eta/D");
  table->Branch("centrality", &row.centrality, "centrality/D");
  table->Branch("native_bin", &row.nativeBin, "native_bin/I");
  table->Branch("iso_r03", &row.isoR03, "iso_r03/D");
  table->Branch("iso_r04", &row.isoR04, "iso_r04/D");
  table->Branch("method", &row.method, "method/I");
  table->Branch("weight", &row.weight, "weight/D");
  table->Branch("file", &row.file, "file/I");
  std::map<std::string, long> counts;
  for (size_t i = 0; i < paths.size(); ++i) fillFile(paths[i], static_cast<int>(i), *table, row, counts);
  out->cd();
  table->Write("", TObject::kOverwrite);  /* one key: an autosave cycle is replaced */
  std::string json = "{";
  for (const auto& [key, value] : counts) json += (json.size() > 1 ? "," : "") + ("\"" + key + "\":" + std::to_string(value));
  TObjString((json + "}").c_str()).Write("counts");
  TObjString(sources.c_str()).Write("sources");
  out->Close();
  return 0;
}

}  /* namespace isolation */
}  /* namespace photonjet */


/* === The package ===================================================================================================== */
namespace treetohists {
namespace fs = std::filesystem;
namespace {

using Ordered = nlohmann::ordered_json;
const int kBins = 700;
const double kLow = -20.0, kHigh = 50.0;
/* photon-pT groups: the group and the ET ranges it collects */
const std::vector<std::pair<std::pair<int, int>, std::vector<std::pair<int, int>>>> kGroups = {
    {{15, 17}, {{15, 17}}}, {{17, 19}, {{17, 19}}}, {{19, 21}, {{19, 21}}}, {{21, 23}, {{21, 23}}}, {{23, 26}, {{23, 26}}},
    {{26, 35}, {{26, 30}, {30, 35}}}};
const std::vector<double> kEfficiencies = {0.70, 0.80, 0.90};

std::vector<double> isolationEdges() {  /* numpy.linspace(-20, 50, 701) */
  std::vector<double> e(kBins + 1);
  const double step = (kHigh - kLow) / kBins;
  for (int i = 0; i <= kBins; ++i) e[i] = i * step + kLow;
  e[kBins] = kHigh;
  return e;
}

struct Population {
  std::vector<double> et, centrality, isoR03, isoR04;
  std::vector<long> native;
  std::vector<int> method;
};

/* The rows of every population file, in order; one subtraction method. */
Population readPopulation(const std::vector<std::string> &paths) {
  Population p;
  for (const auto &path : paths) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    auto *tree = f ? dynamic_cast<TTree *>(f->Get("isolation_population")) : nullptr;
    if (!tree) throw std::runtime_error(path + ": no isolation_population");
    TTreeReader r(tree);
    TTreeReaderValue<double> et(r, "et"), centrality(r, "centrality"), r03(r, "iso_r03"), r04(r, "iso_r04");
    TTreeReaderValue<int> native(r, "native_bin"), method(r, "method");
    while (r.Next()) {
      p.et.push_back(*et), p.centrality.push_back(*centrality), p.isoR03.push_back(*r03), p.isoR04.push_back(*r04);
      p.native.push_back(*native), p.method.push_back(*method);
    }
  }
  for (int m : p.method)
    if (m != p.method.front()) throw std::runtime_error("more than one reconstructed isolation subtraction method");
  return p;
}

/* Counts of the finite values on the isolation edges, with under- and overflow (numpy.histogram: the last bin holds its
   upper edge; a value at the upper edge is also counted as overflow). */
struct Histogram { std::vector<double> counts; double under = 0, over = 0; };
Histogram histogram(const std::vector<double> &values) {
  static const std::vector<double> edges = isolationEdges();
  Histogram h;
  h.counts.assign(kBins, 0.0);
  for (double x : values) {
    if (!std::isfinite(x)) continue;
    if (x < kLow) h.under += 1;
    if (x >= kHigh) h.over += 1;
    if (x < edges.front() || x > edges.back()) continue;
    const int b = x == edges.back() ? kBins - 1 : static_cast<int>(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
    h.counts[b] += 1;
  }
  return h;
}

/* The cutoff at the probability of the cumulative spectrum, with its error and the counts behind it. */
Ordered quantile(const Histogram &h, double p) {
  static const std::vector<double> edges = isolationEdges();
  double inRange = 0;
  for (double v : h.counts) inRange += v;
  const double total = inRange + h.under + h.over;
  if (!std::isfinite(total) || total <= 0) throw std::runtime_error("isolation histogram has no positive in-range weight");
  const double absolute = p * total;
  if (absolute <= h.under) throw std::runtime_error("the quantile lies in the underflow below -20 GeV");
  if (absolute > h.under + inRange) throw std::runtime_error("the quantile lies in the overflow above 50 GeV");
  const double target = absolute - h.under;
  std::vector<double> cumulative(kBins);
  double running = 0;
  for (int i = 0; i < kBins; ++i) cumulative[i] = running += h.counts[i];
  int index = static_cast<int>(std::lower_bound(cumulative.begin(), cumulative.end(), target) - cumulative.begin());
  index = std::min(std::max(index, 0), kBins - 1);
  const double previous = index > 0 ? cumulative[index - 1] : 0.0, weight = h.counts[index];
  const double fraction = std::min(std::max(weight <= 0 ? 0.5 : (target - previous) / weight, 0.0), 1.0);
  const double width = edges[index + 1] - edges[index], threshold = edges[index] + fraction * width;
  const double sumw2 = inRange + h.under + h.over;  /* unweighted: the variances are the counts */
  const double effective = sumw2 > 0 ? total * total / sumw2 : total, density = width > 0 ? weight / (total * width) : 0.0;
  const double statistical = effective > 0 && density > 0 ? std::sqrt(p * (1 - p) / effective) / density : width;
  return Ordered{{"threshold_gev", threshold}, {"threshold_error_gev", std::max(0.5 * width, statistical)}, {"sumw", total},
                 {"sumw2", sumw2}, {"effective_entries", effective}, {"underflow", h.under}, {"overflow", h.over},
                 {"flow_fraction", (h.under + h.over) / std::max(total, 1e-12)}};
}

Ordered constantFit(const std::vector<double> &y, const std::vector<double> &error) {
  double w = 0, wy = 0;
  std::vector<std::pair<double, double>> valid;
  for (size_t i = 0; i < y.size(); ++i)
    if (std::isfinite(y[i]) && std::isfinite(error[i]) && error[i] > 0) valid.push_back({y[i], error[i]});
  if (valid.size() < 2) throw std::runtime_error("a constant fit needs at least two valid points");
  for (const auto &[v, e] : valid) w += 1 / (e * e), wy += v / (e * e);
  const double value = wy / w;
  double chi2 = 0;
  for (const auto &[v, e] : valid) chi2 += (v - value) * (v - value) / (e * e);
  return Ordered{{"value_gev", value}, {"error_gev", std::sqrt(1 / w)}, {"chi2", chi2}, {"ndf", static_cast<int>(valid.size()) - 1},
                 {"n_points", static_cast<int>(valid.size())}};
}

Ordered linearFit(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &error) {
  std::vector<std::array<double, 3>> valid;
  for (size_t i = 0; i < x.size(); ++i)
    if (std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(error[i]) && error[i] > 0) valid.push_back({x[i], y[i], error[i]});
  if (valid.size() < 3) throw std::runtime_error("a linear fit needs at least three valid points");
  double s = 0, sx = 0, sxx = 0, sy = 0, sxy = 0;
  for (const auto &[a, b, e] : valid) {
    const double w = 1 / (e * e);
    s += w, sx += w * a, sxx += w * a * a, sy += w * b, sxy += w * a * b;
  }
  const double determinant = s * sxx - sx * sx;
  const double c00 = sxx / determinant, c01 = -sx / determinant, c11 = s / determinant;  /* (D^T W D)^-1 */
  const double intercept = c00 * sy + c01 * sxy, slope = c01 * sy + c11 * sxy;
  double chi2 = 0;
  for (const auto &[a, b, e] : valid) chi2 += (b - intercept - slope * a) * (b - intercept - slope * a) / (e * e);
  return Ordered{{"intercept_gev", intercept}, {"slope_gev_per_percent", slope}, {"intercept_error_gev", std::sqrt(c00)},
                 {"slope_error_gev_per_percent", std::sqrt(c11)}, {"covariance", {{c00, c01}, {c01, c11}}}, {"chi2", chi2},
                 {"ndf", static_cast<int>(valid.size()) - 2}, {"n_points", static_cast<int>(valid.size())}};
}

/* (cutoff points, constant fits per slice, lines per efficiency) of one cone. */
std::array<Ordered, 3> derive(const std::vector<double> &et, const std::vector<long> &native, const std::vector<double> &iso) {
  Ordered points = Ordered::array(), flat = Ordered::array(), lines = Ordered::array();
  for (int k = 0; k < 16; ++k) {
    const int centMin = 5 * k, centMax = 5 * k + 5;
    for (const auto &[group, sources] : kGroups) {
      std::vector<double> selected;
      for (size_t i = 0; i < et.size(); ++i) {
        if (native[i] < 1 || (native[i] - 1) / 5 != k) continue;
        bool inPt = false;
        for (const auto &[lo, hi] : sources) inPt = inPt || (et[i] >= lo && et[i] < hi);
        if (inPt) selected.push_back(iso[i]);
      }
      const Histogram h = histogram(selected);
      for (double efficiency : kEfficiencies) {
        Ordered point = {{"efficiency", efficiency}, {"cent_min", centMin}, {"cent_max", centMax}, {"cent_center", 0.5 * (centMin + centMax)},
                         {"pt_min", group.first}, {"pt_max", group.second}, {"pt_center", 0.5 * (group.first + group.second)},
                         {"photons", static_cast<long>(selected.size())}};
        point.update(quantile(h, efficiency));
        points.push_back(point);
      }
    }
  }
  for (int k = 0; k < 16; ++k)
    for (double efficiency : kEfficiencies) {
      std::vector<double> y, e;
      for (const auto &p : points)
        if (p.at("cent_min") == 5 * k && p.at("efficiency") == efficiency) y.push_back(p.at("threshold_gev")), e.push_back(p.at("threshold_error_gev"));
      Ordered fit = {{"efficiency", efficiency}, {"cent_min", 5 * k}, {"cent_max", 5 * k + 5}, {"cent_center", 5 * k + 2.5}};
      fit.update(constantFit(y, e));
      flat.push_back(fit);
    }
  for (double efficiency : kEfficiencies) {
    std::vector<double> x, y, e;
    for (const auto &f : flat)
      if (f.at("efficiency") == efficiency) x.push_back(f.at("cent_center")), y.push_back(f.at("value_gev")), e.push_back(f.at("error_gev"));
    Ordered line = {{"efficiency", efficiency}};
    line.update(linearFit(x, y, e));
    lines.push_back(line);
  }
  return {points, flat, lines};
}

Ordered inputsOf(const std::vector<std::string> &paths) {
  Ordered out = Ordered::array();
  for (const auto &p : paths) out.push_back(Ordered{{"path", p}, {"sha256", sha256File(p)}});
  return out;
}

Ordered package(const std::vector<std::string> &paths) {
  const Population p = readPopulation(paths);
  Ordered groups = Ordered::array();
  for (const auto &[group, sources] : kGroups) {
    Ordered s = Ordered::array();
    for (const auto &[lo, hi] : sources) s.push_back({lo, hi});
    groups.push_back({{group.first, group.second}, s});
  }
  Ordered out = {{"schema", "AuAuIsolationPackageV1"},
                 {"method", "per pT group and 5% centrality slice: the 70/80/90% cutoff of the reconstructed-isolation spectrum of "
                            "truth-signal photons; a constant across pT per slice; a line across centrality"},
                 {"histogram", {{"nbins", kBins}, {"xmin", kLow}, {"xmax", kHigh}, {"weighted", false}}},
                 {"pt_groups", groups}, {"efficiencies", kEfficiencies},
                 {"centrality_convention", "native bin b (1-based); slice k = bins 5k+1..5k+5"},
                 {"photons", static_cast<long>(p.et.size())}, {"subtraction_method", p.method.empty() ? Ordered() : Ordered(p.method.front())},
                 {"inputs", inputsOf(paths)}, {"cones", Ordered::object()}};
  for (const auto &[cone, iso] : std::vector<std::pair<std::string, const std::vector<double> *>>{{"R30", &p.isoR03}, {"R40", &p.isoR04}}) {
    const auto [points, flat, lines] = derive(p.et, p.native, *iso);
    out["cones"][cone] = Ordered{{"lines", lines}, {"constant_fits", flat}, {"cutoff_points", points}};
  }
  return out;
}

Ordered freeze(const std::string &derivation, double efficiency, double gap) {
  std::ifstream in(derivation);
  const Ordered doc = Ordered::parse(in);
  if (doc.value("schema", "") != "AuAuIsolationPackageV1") throw std::runtime_error("not an isolation derivation package");
  const Ordered *line = nullptr;
  for (const auto &l : doc.at("cones").at("R30").at("lines"))
    if (std::fabs(l.at("efficiency").get<double>() - efficiency) < 1e-12) line = &l;
  if (!line) throw std::runtime_error("no R = 0.3 line at the efficiency");
  return Ordered{{"schema", "CanonicalIsolationFreezeV1"}, {"system", "auau"},
                 {"witness", {{"radius", 0.3}, {"stride_index", 0}, {"subtraction_method", doc.at("subtraction_method")}, {"reconstructed", true}}},
                 {"threshold", {{"intercept_gev", line->at("intercept_gev")}, {"slope_gev_per_percent", line->at("slope_gev_per_percent")},
                                {"variable", "cent_value"}, {"efficiency", efficiency}}},
                 {"sideband_gap_gev", gap}, {"domain", {{"centrality_min", 0.0}, {"centrality_max_exclusive", 80.0}}},
                 {"derivation", {{"file", fs::path(derivation).filename().string()}, {"package_sha256", sha256File(derivation)}}}};
}

Ordered check(const Ordered &frozen, const Population &p) {
  if (frozen.value("schema", "") != "CanonicalIsolationFreezeV1") throw std::runtime_error("not an isolation freeze");
  const Ordered &line = frozen.at("threshold");
  const double target = line.at("efficiency"), a = line.at("intercept_gev"), b = line.at("slope_gev_per_percent");
  std::vector<bool> inside(p.et.size()), isolated(p.et.size());
  for (size_t i = 0; i < p.et.size(); ++i) {
    inside[i] = p.et[i] >= 15.0 && p.et[i] < 35.0 && std::isfinite(p.isoR03[i]) && std::isfinite(p.centrality[i]) && p.native[i] >= 1;
    isolated[i] = p.isoR03[i] < a + b * p.centrality[i];
  }
  auto fraction = [&](const std::function<bool(size_t)> &selected) {
    long n = 0, pass = 0;
    for (size_t i = 0; i < p.et.size(); ++i)
      if (selected(i)) ++n, pass += isolated[i];
    if (n == 0) return Ordered{{"photons", 0}, {"efficiency", nullptr}, {"error", nullptr}, {"pull", nullptr}};
    const double f = static_cast<double>(pass) / n, error = std::sqrt(std::max(f * (1 - f), 0.0) / n);
    return Ordered{{"photons", n}, {"efficiency", f}, {"error", error}, {"pull", error > 0 ? Ordered((f - target) / error) : Ordered()}};
  };
  Ordered slices = Ordered::array();
  for (int k = 0; k < 16; ++k) {
    Ordered s = {{"cent_min", 5 * k}, {"cent_max", 5 * k + 5}};
    s.update(fraction([&](size_t i) { return inside[i] && (p.native[i] - 1) / 5 == k; }));
    slices.push_back(s);
  }
  Ordered out = {{"schema", "AuAuIsolationCheckV1"}, {"target_efficiency", target},
                 {"frozen_line", {{"intercept_gev", a}, {"slope_gev_per_percent", b}}}, {"slices", slices},
                 {"overall", fraction([&](size_t i) { return static_cast<bool>(inside[i]); })}};
  std::vector<double> et, iso;
  std::vector<long> native;
  for (size_t i = 0; i < p.et.size(); ++i)
    if (inside[i]) et.push_back(p.et[i]), native.push_back(p.native[i]), iso.push_back(p.isoR03[i]);
  try {
    const auto lines = derive(et, native, iso)[2];
    const Ordered *refit = nullptr;
    for (const auto &l : lines)
      if (std::fabs(l.at("efficiency").get<double>() - target) < 1e-12) refit = &l;
    const double da = refit->at("intercept_gev").get<double>() - a, db = refit->at("slope_gev_per_percent").get<double>() - b;
    const auto &c = refit->at("covariance");
    Ordered difference = Ordered::array(), error = Ordered::array(), centres = Ordered::array();
    for (int k = 0; k < 16; ++k) {
      const double x = 5 * k + 2.5;
      difference.push_back(da + db * x);
      error.push_back(std::sqrt(c[0][0].get<double>() + 2 * x * c[0][1].get<double>() + x * x * c[1][1].get<double>()));
      centres.push_back(x);
    }
    out["refit"] = {{"line", *refit}, {"threshold_minus_frozen_gev", difference}, {"error_gev", error}, {"slice_centres", centres}};
  } catch (const std::runtime_error &e) {  /* too few photons in a slice and pT group to refit */
    out["refit"] = {{"error", e.what()}};
  }
  return out;
}

void write(const std::string &path, const Ordered &doc) {
  if (fs::exists(path)) throw std::runtime_error(path + " exists");
  std::ofstream(path) << doc.dump(1) << "\n";
}

}  /* namespace */

/* tree_to_hists isolation population|derive|freeze|check ...: argv from the subcommand on. */
/* population stores eligible truth-signal witnesses; derive fits cutoffs,
   freeze binds the chosen line, and check evaluates that fixed definition on
   a supplied population. Ordinary fills consume the resulting freeze. */
int isolation(const std::vector<std::string> &args) {
  if (args.size() < 3 || args[1] != "--out") throw std::runtime_error("usage: tree_to_hists isolation population|derive|freeze|check --out FILE INPUTS...");
  const std::string command = args[0], out = args[2];
  std::vector<std::string> inputs;
  double efficiency = 0.80, gap = 0.8;
  for (size_t i = 3; i < args.size(); ++i) {
    if (args[i] == "--efficiency" && i + 1 < args.size()) efficiency = std::stod(args[++i]);
    else if (args[i] == "--gap" && i + 1 < args.size()) gap = std::stod(args[++i]);
    else inputs.push_back(args[i]);
  }
  if (inputs.empty()) throw std::runtime_error("isolation " + command + ": no inputs");
  if (command == "population") {
    if (fs::exists(out)) throw std::runtime_error(out + " exists");
    std::string list;
    for (const auto &i : inputs) list += i + " ";
    return photonjet::isolation::run(list.c_str(), out.c_str());
  }
  if (command == "derive") {
    const Ordered doc = package(inputs);
    write(out, doc);
    for (const auto &[cone, block] : doc.at("cones").items())
      for (const auto &l : block.at("lines"))
        std::cout << cone << " " << l.at("efficiency").get<double>() << ": " << l.at("intercept_gev").get<double>() << " + "
                  << l.at("slope_gev_per_percent").get<double>() << " c (chi2/ndf " << l.at("chi2").get<double>() << "/" << l.at("ndf") << ")\n";
    return 0;
  }
  if (command == "freeze") {
    write(out, freeze(inputs.at(0), efficiency, gap));
    return 0;
  }
  if (command == "check") {
    std::ifstream in(inputs.at(0));
    const Ordered frozen = Ordered::parse(in);
    const std::vector<std::string> populations(inputs.begin() + 1, inputs.end());
    const Population p = readPopulation(populations);
    if (!p.method.empty() && p.method.front() != frozen.at("witness").at("subtraction_method").get<int>())
      throw std::runtime_error("the population's isolation method differs from the freeze's witness");
    Ordered result = check(frozen, p);
    result["freeze"] = {{"path", inputs.at(0)}, {"sha256", sha256File(inputs.at(0))}};
    result["inputs"] = inputsOf(populations);
    write(out, result);
    for (const auto &s : result.at("slices"))
      if (s.at("photons").get<long>() > 0)
        std::cout << s.at("cent_min") << "-" << s.at("cent_max") << "% " << s.at("photons") << " photons: " << s.at("efficiency").get<double>()
                  << " +- " << s.at("error").get<double>() << "\n";
    return 0;
  }
  throw std::runtime_error("unknown isolation command " + command);
}

}  /* namespace treetohists */
