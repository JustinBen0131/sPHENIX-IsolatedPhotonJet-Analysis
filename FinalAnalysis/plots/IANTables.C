/* IANTables: the IAN result tables, LaTeX (booktabs) and CSV, from results.root and, when given, systematics.root. */

/* ./FinalAnalysis/final_analysis.sh tables --results results.root [--systematics systematics.root] --out DIR */

/* tab_xj_<result>.tex, .csv  the unfolded per-photon xJ distribution of one result over its reported bins: value,
   data and response half-widths, the systematic uncertainty (with systematics) and the total
   tab_xj_all.tex             every result side by side, value +- total, over the xJ bins every result reports
   tab_measurement.tex        per result and photon-pT bin: the iterations of the joint and photon unfoldings, the
   purity, the purity-corrected photons and the recoil jets per photon in the reported xJ
   range (raw region-A DATA, simulated reconstruction)
   Nothing is recomputed beyond sums over the reported xJ bins and sums in quadrature. */
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TObjString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "../src/RenderMetadata.h"

namespace {
using Json = photonjet_render::Json;
Json gTableMetadata, gTableCoverage, gTableInputs;
std::string coverageNote() {
  const auto status = gTableCoverage.value("status",std::string("UNKNOWN"));
  if (status=="REFERENCE_INPUTS") return " Reference inputs for software development.";
  if (status=="PARTIAL") return " Partial input population.";
  if (status=="UNKNOWN") return " Input population coverage is unverified.";
  return "";
}

struct Entry {
  std::string key, title;          /* the result's directory and its LaTeX name */
  std::map<std::string, std::string> settings;
  double order = 0;                /* p+p first, then Au+Au from peripheral to central */
};

std::string format(const char *f, double x) {
  char b[64];
  snprintf(b, sizeof b, f, x);
  return b;
}
std::string fixed(double x) { return std::isfinite(x) ? format("%.3f", x) : "--"; }

std::vector<Entry> entriesOf(TFile &f) {
  std::vector<Entry> out;
  TIter next(f.GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next())) {
    auto *s = dynamic_cast<TObjString *>(f.Get((std::string(k->GetName()) + "/settings").c_str()));
    if (!s) continue;
    Entry e;
    e.key = k->GetName();
    std::istringstream in(s->GetString().Data());
    for (std::string item; std::getline(in, item, ';');)
      if (item.find('=') != std::string::npos) e.settings[item.substr(0, item.find('='))] = item.substr(item.find('=') + 1);
    const std::string c = e.settings["centrality"];
    const auto u = c.find('_');
    if (e.settings["package"] == "PP") e.title = "\\pp", e.order = -1;
    else if (c.rfind("cent", 0) == 0 && u != std::string::npos)
      e.title = c.substr(4, u - 4) + "--" + c.substr(u + 1) + "\\% \\auau", e.order = -std::atof(c.c_str() + 4) + 1000;
    else e.title = e.key, e.order = 2000;
    if (!gTableMetadata.is_null()) {
      const auto &saved = gTableMetadata.at("results").at(e.key);
      if (saved.at("package").get<std::string>()!=e.settings.at("package") || saved.at("centrality").get<std::string>()!=c)
        throw std::runtime_error("table identity differs from the saved result metadata");
      const auto &range = saved.at("centrality_percent");
      if (!range.is_null()) {
        e.title = format("%g",range.at(0).get<double>())+"--"+format("%g",range.at(1).get<double>())+"\\% \\auau";
        e.order = 1000-range.at(0).get<double>();
      }
    }
    out.push_back(e);
  }
  std::stable_sort(out.begin(), out.end(), [](const Entry &a, const Entry &b) { return a.order < b.order; });
  return out;
}

TH1 *get(TFile *f, const std::string &key, const std::string &name) { return f ? dynamic_cast<TH1 *>(f->Get((key + "/" + name).c_str())) : nullptr; }

/* The reported xJ bins of a result: [low, high] of settings reported_xj. */
std::vector<int> reportedBins(const TH1 *xj, const Entry &e) {
  double low = 0, high = 0;
  char comma;
  std::istringstream(e.settings.at("reported_xj")) >> low >> comma >> high;
  std::vector<int> bins;
  for (int b = 1; b <= xj->GetNbinsX(); ++b)
    if (xj->GetXaxis()->GetBinLowEdge(b) >= low - 1e-9 && xj->GetXaxis()->GetBinUpEdge(b) <= high + 1e-9) bins.push_back(b);
  return bins;
}

/* The systematic sigma of a bin of a result (systematics.root holds the reported bins only); NaN without it. */
double systematicOf(TFile *s, const std::string &key, double low, double high) {
  TH1 *h = get(s, key, "xj_systematic_sigma");
  if (!h) return std::nan("");
  for (int b = 1; b <= h->GetNbinsX(); ++b)
    if (std::fabs(h->GetXaxis()->GetBinLowEdge(b) - low) < 1e-9 && std::fabs(h->GetXaxis()->GetBinUpEdge(b) - high) < 1e-9) return h->GetBinContent(b);
  return std::nan("");
}

}  /* namespace */

/* Export the saved reported-bin values and uncertainties to LaTeX and CSV, using the same input checks as figures.
   No correction, fit or unfolding is performed during table export. */
void IANTables(const char *results, const char *out, const char *systematics = "") {
  std::unique_ptr<TFile> f(TFile::Open(results, "READ")), s;
  if (!f || f->IsZombie()) throw std::runtime_error(std::string("IANTables: cannot open ") + results);
  gTableInputs = Json::object(); gTableMetadata = Json(); gTableCoverage = Json::object();
  gTableInputs[results] = photonjet_render::sha256(results,true);
  if (auto *saved = dynamic_cast<TObjString *>(f->Get("render_metadata"))) {
    const auto sealed = Json::parse(saved->GetString().Data());
    const auto receipt = photonjet_render::read(std::string(results)+".render.json");
    if (receipt.at("product_sha256")!=gTableInputs.at(results) || receipt.at("metadata_sha256")!=sealed.at("sha256"))
      throw std::runtime_error("table input or annotation metadata changed");
    gTableMetadata = photonjet_render::verify(sealed);
  }
  if (auto *coverage = dynamic_cast<TObjString *>(f->Get("input_coverage"))) gTableCoverage = Json::parse(coverage->GetString().Data());
  if (systematics && *systematics) {
    s.reset(TFile::Open(systematics, "READ"));
    if (!s || s->IsZombie()) throw std::runtime_error(std::string("IANTables: cannot open ") + systematics);
    gTableInputs[systematics] = photonjet_render::sha256(systematics,true);
    if (auto *nominal = dynamic_cast<TObjString *>(s->Get("nominal_input_sha256")))
      if (std::string(nominal->GetString().Data())!=gTableInputs.at(results).get<std::string>())
        throw std::runtime_error("table uncertainties belong to a different nominal result");
  }
  const bool withSystematics = s != nullptr;
  gSystem->mkdir(out, true);
  const std::string dir = out;
  const std::vector<Entry> entries = entriesOf(*f);
  if (entries.empty()) throw std::runtime_error("IANTables: no results in the supplied file");
  std::vector<std::string> written;
  /* value and total per entry and reported bin, for the side-by-side table */
  std::map<std::string, std::map<std::pair<double, double>, std::pair<double, double>>> totals;
  for (const auto &e : entries) {
    TH1 *xj = get(f.get(), e.key, "xj_per_photon"), *data = get(f.get(), e.key, "xj_half_width_data"),
        *response = get(f.get(), e.key, "xj_half_width_response"), *purity = get(f.get(), e.key, "purity");
    if (!xj || !data || !response || !purity) throw std::runtime_error("IANTables: missing required result for " + e.key);
    if (!gTableMetadata.is_null()) {
      const auto &edges = gTableMetadata.at("results").at(e.key).at("photon_edges_gev");
      if (edges.front().get<double>()!=purity->GetXaxis()->GetXmin() || edges.back().get<double>()!=purity->GetXaxis()->GetXmax())
        throw std::runtime_error("table photon bins differ from the executed selection");
    }
    std::ostringstream tex, csv;
    tex << "\\begin{table}[hbtp!]\n  \\centering\n  \\small\n"
        << "  \\caption{Unfolded per-photon $x_{J\\gamma}$ distribution in " << e.title << ", $"
        << format("%g", purity->GetXaxis()->GetXmin()) << "\\leq E_{T}^{\\gamma}<" << format("%g", purity->GetXaxis()->GetXmax())
        << "$~\\GeV, $p_{T}^{\\mathrm{jet}}\\geq" << e.settings.at("jet_threshold_gev") << "$~\\GeV: value, statistical (data) and response half-widths, "
        << (withSystematics ? "systematic uncertainty, " : "") << "and the total uncertainty of the measurement." << coverageNote() << "}\n"
        << "  \\label{tab:xj_" << e.key << "}\n  \\begin{tabular}{" << (withSystematics ? "cccccc" : "ccccc") << "}\n    \\toprule\n"
        << "    $x_{J\\gamma}$ interval & $(1/N_{\\gamma})\\,dN_{\\mathrm{pair}}/dx_{J\\gamma}$ & data & response & "
        << (withSystematics ? "systematic & " : "") << "total \\\\\n    \\midrule\n";
    csv << "xj_low,xj_high,value,half_width_68_data,half_width_68_response," << (withSystematics ? "systematic," : "") << "sigma_total\n";
    for (int b : reportedBins(xj, e)) {
      const double low = xj->GetXaxis()->GetBinLowEdge(b), high = xj->GetXaxis()->GetBinUpEdge(b);
      const double v = xj->GetBinContent(b), d = data->GetBinContent(b), r = response->GetBinContent(b);
      const double syst = withSystematics ? systematicOf(s.get(), e.key, low, high) : 0.0;
      if (!std::isfinite(v) || !std::isfinite(d) || !std::isfinite(r) || !std::isfinite(syst) || d < 0 || r < 0 || syst < 0)
        throw std::runtime_error("IANTables: invalid or missing uncertainty for " + e.key);
      const double total = std::sqrt(d * d + r * r + (std::isfinite(syst) ? syst * syst : 0.0));
      tex << "    " << format("%g", low) << "--" << format("%g", high) << " & " << fixed(v) << " & " << fixed(d) << " & " << fixed(r) << " & "
          << (withSystematics ? fixed(syst) + " & " : "") << fixed(total) << " \\\\\n";
      csv << format("%.10g", low) << "," << format("%.10g", high) << "," << format("%.10g", v) << "," << format("%.10g", d) << ","
          << format("%.10g", r) << "," << (withSystematics ? format("%.10g", syst) + "," : "") << format("%.10g", total) << "\n";
      totals[e.key][{low, high}] = {v, total};
    }
    tex << "    \\bottomrule\n  \\end{tabular}\n\\end{table}\n";
    std::ofstream(dir + "/tab_xj_" + e.key + ".tex") << tex.str();
    std::ofstream(dir + "/tab_xj_" + e.key + ".csv") << csv.str();
    written.push_back("tab_xj_" + e.key + ".tex"), written.push_back("tab_xj_" + e.key + ".csv");
  }
  /* every result side by side over the bins all of them report */
  if (totals.size() > 1) {
    std::vector<std::pair<double, double>> bins;
    for (const auto &[bin, value] : totals.begin()->second) {
      bool everywhere = true;
      for (const auto &[key, rows] : totals) everywhere = everywhere && rows.count(bin);
      if (everywhere) bins.push_back(bin);
    }
    std::ostringstream tex;
    std::string columns = "c", header;
    for (const auto &e : entries)
      if (totals.count(e.key)) columns += "c", header += " & " + e.title;
    tex << "\\begin{table}[hbtp!]\n  \\centering\n  \\small\n  \\caption{Unfolded per-photon $x_{J\\gamma}$ distributions, value $\\pm$ total uncertainty." << coverageNote() << "}\n"
        << "  \\label{tab:xj_all}\n  \\begin{tabular}{" << columns << "}\n    \\toprule\n    $x_{J\\gamma}$" << header << " \\\\\n    \\midrule\n";
    for (const auto &bin : bins) {
      tex << "    " << format("%g", bin.first) << "--" << format("%g", bin.second);
      for (const auto &e : entries)
        if (totals.count(e.key)) tex << " & $" << fixed(totals[e.key][bin].first) << " \\pm " << fixed(totals[e.key][bin].second) << "$";
      tex << " \\\\\n";
    }
    tex << "    \\bottomrule\n  \\end{tabular}\n\\end{table}\n";
    std::ofstream(dir + "/tab_xj_all.tex") << tex.str();
    written.push_back("tab_xj_all.tex");
  }
  /* the measurement summary */
  std::ostringstream tex;
  tex << "\\begin{table}[hbtp!]\n  \\centering\n  \\small\n"
      << "  \\caption{Measurement summary per collision system and photon-$p_{T}$ interval: Bayesian iterations of the joint and photon "
         "unfoldings, leakage-corrected ABCD purity, purity-corrected photons, and recoil jets per photon in DATA and simulation "
         "(reconstructed level)." << coverageNote() << "}\n"
      << "  \\label{tab:measurement_summary}\n  \\begin{tabular}{llcccccc}\n    \\toprule\n"
      << "    System & $p_{T}^{\\gamma}$ [\\GeV] & $k_{\\mathrm{pairs}}$ & $k_{\\gamma}$ & purity & $N_{\\gamma}$ & jets/$\\gamma$ DATA & jets/$\\gamma$ sim \\\\\n"
      << "    \\midrule\n";
  for (const auto &e : entries) {
    TH1 *purity = get(f.get(), e.key, "purity"), *photons = get(f.get(), e.key, "photons_corrected"), *xj = get(f.get(), e.key, "xj_per_photon");
    auto *pairsA = dynamic_cast<TH2 *>(get(f.get(), e.key, "pairs_A")), *abcd = dynamic_cast<TH2 *>(get(f.get(), e.key, "abcd"));
    auto *pairsSim = dynamic_cast<TH2 *>(get(f.get(), e.key, "pairs_reco_simulation"));
    TH1 *photonsSim = get(f.get(), e.key, "photon_reco");
    if (!purity || !photons || !xj) continue;
    const std::vector<int> shown = reportedBins(xj, e);
    auto perPhoton = [&](TH2 *pairs, double n, int pt) {
      if (!pairs || !(n > 0)) return std::nan("");
      double sum = 0;
      for (int b : shown) sum += pairs->GetBinContent(pt, b);
      return sum / n;
    };
    for (int i = 1; i <= purity->GetNbinsX(); ++i)
      tex << "    " << (i == 1 ? e.title : "") << " & " << format("%g", purity->GetXaxis()->GetBinLowEdge(i)) << "--"
          << format("%g", purity->GetXaxis()->GetBinUpEdge(i)) << " & " << e.settings.at("pair_iterations") << " & "
          << e.settings.at("photon_iterations") << " & " << fixed(purity->GetBinContent(i)) << " & " << format("%.0f", photons->GetBinContent(i))
          << " & " << fixed(perPhoton(pairsA, abcd ? abcd->GetBinContent(1, i) : 0.0, i)) << " & "
          << fixed(perPhoton(pairsSim, photonsSim ? photonsSim->GetBinContent(i) : 0.0, i)) << " \\\\\n";
  }
  tex << "    \\bottomrule\n  \\end{tabular}\n\\end{table}\n";
  std::ofstream(dir + "/tab_measurement.tex") << tex.str();
  written.push_back("tab_measurement.tex");
  std::ofstream list((dir + "/TABLES.txt").c_str());
  list << "# Tables written by FinalAnalysis/plots/IANTables.C from " << results << (withSystematics ? std::string(" and ") + systematics : "") << "\n";
  for (const auto &w : written) list << w << "\n";
  if (!list) throw std::runtime_error("IANTables: cannot save TABLES.txt");
  for (const auto &w : written) {
    std::ifstream saved(dir + "/" + w, std::ios::binary | std::ios::ate);
    if (!saved || saved.tellg() <= 0) throw std::runtime_error("IANTables: table was not saved: " + w);
    photonjet_render::write(dir+"/"+w+".json", {{"schema","PhotonJetTableRenderReceiptV1"},
      {"payload_sha256",photonjet_render::sha256(dir+"/"+w,true)}, {"inputs",gTableInputs}, {"input_coverage",gTableCoverage},
      {"result_metadata",gTableMetadata}, {"renderer_sha256",photonjet_render::sha256(__FILE__,true)}});
  }
  printf("IANTables: %zu tables in %s\n", written.size(), out);
}
