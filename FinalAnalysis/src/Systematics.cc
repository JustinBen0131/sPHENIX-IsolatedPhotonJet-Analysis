/* Systematics: the systematic uncertainties of the per-photon results, from results files that differ by one change each. */

/* final_analysis --systematics config/systematics.cfg */

/* The settings name the nominal results file and every source; a variation is a copy of the analysis settings with one
   change, run like the nominal. Per source and result, the signed shifts (varied minus nominal, per reported xJ bin):
   Up, Down                 a two-sided nuisance: covariance d d^T with d = (up - down)/2; the band takes, per bin, the
   larger positive and the larger negative of the two shifts
   Alternative              one-sided, mirrored: d = alternative - nominal, band +-|d|
   Alternatives, Representative   an envelope: the band is the largest |shift| of the alternatives, the covariance
   d d^T of the named representative
   Sources are independent: covariances add, bands add in quadrature. The statistical covariance is the nominal's own
   (data plus response toys). A ratio A/B varies a source present in both results together when it is Correlated (the
   same physical input), and each system alone otherwise; the two systems' statistics are independent. */
#include "FinalAnalysis.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THashList.h"
#include "TKey.h"
#include "TObjString.h"
#include "RenderMetadata.h"
#include "TSystem.h"

namespace fa {
namespace {

std::string setting(TEnv &env, const std::string &key, bool required = true) {
  std::string v = env.GetValue(key.c_str(), "");
  v = v.substr(0, v.find('#'));
  const auto first = v.find_first_not_of(" \t"), last = v.find_last_not_of(" \t");
  v = first == std::string::npos ? "" : v.substr(first, last - first + 1);
  if (v.empty() && required) throw std::runtime_error("setting " + key + " is missing");
  return v;
}

std::vector<std::string> words(const std::string &text) {
  std::istringstream in(text);
  std::vector<std::string> out;
  for (std::string w; in >> w;) out.push_back(w);
  return out;
}

/* The reported part of one result of a results file. */
struct Result {
  Edges edges;         /* the reported xJ bins */
  TVectorD central;    /* per-photon value per reported bin */
  TMatrixD statistical;  /* data plus response toy covariance */
};

std::map<std::string, Result> readResults(const std::string &path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  std::map<std::string, Result> out;
  TIter next(f->GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    if (std::string(key->GetClassName()) != "TDirectoryFile") continue;
    const std::string name = key->GetName();
    auto *settings = dynamic_cast<TObjString *>(f->Get((name + "/settings").c_str()));
    auto *h = dynamic_cast<TH1 *>(f->Get((name + "/xj_per_photon").c_str()));
    auto *cd = dynamic_cast<TH2 *>(f->Get((name + "/xj_covariance_data").c_str()));
    auto *cr = dynamic_cast<TH2 *>(f->Get((name + "/xj_covariance_response").c_str()));
    if (!settings || !h || !cd || !cr) throw std::runtime_error(path + ": " + name + " lacks its result or covariance");
    const std::string text = settings->GetString().Data(), tag = "reported_xj=";
    const auto at = text.find(tag);
    std::istringstream range(text.substr(at + tag.size(), text.find(';', at) - at - tag.size()));
    double low = 0, high = 0;
    char comma;
    range >> low >> comma >> high;
    std::vector<int> bins;
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (h->GetXaxis()->GetBinLowEdge(b) >= low - 1e-9 && h->GetXaxis()->GetBinUpEdge(b) <= high + 1e-9) bins.push_back(b);
    Result r;
    r.central.ResizeTo(bins.size());
    r.statistical.ResizeTo(bins.size(), bins.size());
    for (size_t i = 0; i < bins.size(); ++i) {
      r.edges.push_back(h->GetXaxis()->GetBinLowEdge(bins[i]));
      r.central[i] = h->GetBinContent(bins[i]);
      for (size_t j = 0; j < bins.size(); ++j) r.statistical(i, j) = cd->GetBinContent(bins[i], bins[j]) + cr->GetBinContent(bins[i], bins[j]);
    }
    r.edges.push_back(h->GetXaxis()->GetBinUpEdge(bins.back()));
    out.emplace(name, r);
  }
  return out;
}

struct Shift {
  std::string kind;                        /* updown, one_sided, envelope */
  std::map<std::string, TVectorD> templates;
  TVectorD covarianceTemplate, bandUp, bandDown;
};

TVectorD difference(const Result &varied, const Result &nominal, const std::string &label) {
  if (varied.edges != nominal.edges) throw std::runtime_error(label + ": reported bins differ from the nominal's");
  return varied.central - nominal.central;
}

TH1D *hist(const std::string &name, const Edges &edges, const TVectorD &v, const TVectorD *errors = nullptr) {
  auto *h = new TH1D(name.c_str(), "", edges.size() - 1, edges.data());
  for (int i = 0; i < v.GetNrows(); ++i) h->SetBinContent(i + 1, v[i]), h->SetBinError(i + 1, errors ? (*errors)[i] : 0.0);
  return h;
}

TH2D *hist(const std::string &name, const Edges &edges, const TMatrixD &m) {
  auto *h = new TH2D(name.c_str(), "", edges.size() - 1, edges.data(), edges.size() - 1, edges.data());
  for (int i = 0; i < m.GetNrows(); ++i) for (int j = 0; j < m.GetNcols(); ++j) h->SetBinContent(i + 1, j + 1, m(i, j));
  return h;
}

TVectorD diagonalSqrt(const TMatrixD &m) {
  TVectorD v(m.GetNrows());
  for (int i = 0; i < v.GetNrows(); ++i) v[i] = std::sqrt(std::max(m(i, i), 0.0));
  return v;
}

}  /* namespace */

/* Compare saved nominal and variation results on compatible axes, then build nuisance covariances and ratios.
   Variation production is a separate step; this command does not reread or reprocess eventTrees. */
int systematics(const std::string &configPath) {
  TEnv env;
  if (env.ReadFile(configPath.c_str(), kEnvUser) != 0) throw std::runtime_error("cannot read " + configPath);
  const std::string base = gSystem->GetDirName(configPath.c_str()).Data();
  auto path = [&](const std::string &p) { return p.empty() || p[0] == '/' ? p : base + "/" + p; };
  std::set<std::string> known = {"Systematics.Nominal", "Systematics.Output", "Systematics.Sources", "Systematics.Ratios"};
  photonjet_render::Json inputHashes = photonjet_render::Json::object();
  auto pinnedResults = [&](const std::string &file) {
    inputHashes[file] = photonjet_render::sha256(file,true);
    return readResults(file);
  };
  const std::string nominalPath = path(setting(env, "Systematics.Nominal"));
  const auto nominal = pinnedResults(nominalPath);
  const std::string output = path(setting(env, "Systematics.Output"));
  if (!gSystem->AccessPathName(output.c_str())) throw std::runtime_error(output + " exists; choose a new output");
  std::map<std::string, std::map<std::string, Shift>> shifts;  /* system -> source -> shift */
  std::map<std::string, bool> correlated;
  for (const auto &source : words(setting(env, "Systematics.Sources"))) {
    for (const char *k : {"Family", "Correlated", "Up", "Down", "Alternative", "Alternatives", "Representative"}) known.insert(source + "." + k);
    setting(env, source + ".Family");
    correlated[source] = setting(env, source + ".Correlated") == "1";
    const std::string up = setting(env, source + ".Up", false), alternative = setting(env, source + ".Alternative", false),
                      alternatives = setting(env, source + ".Alternatives", false);
    if ((!up.empty()) + (!alternative.empty()) + (!alternatives.empty()) != 1)
      throw std::runtime_error(source + ": give Up and Down, an Alternative, or Alternatives");
    std::map<std::string, std::map<std::string, Result>> files;  /* label -> results */
    if (!up.empty()) files["up"] = pinnedResults(path(up)), files["down"] = pinnedResults(path(setting(env, source + ".Down")));
    else if (!alternative.empty()) files["alternative"] = pinnedResults(path(alternative));
    else for (const auto &item : words(alternatives)) {
      const auto eq = item.find('=');
      if (eq == std::string::npos) throw std::runtime_error(source + ".Alternatives: name=path items");
      files[item.substr(0, eq)] = pinnedResults(path(item.substr(eq + 1)));
    }
    for (const auto &[system, reference] : nominal) {
      bool present = true;
      for (const auto &[label, results] : files) present = present && results.count(system);
      if (!present) continue;  /* the source does not vary this result */
      Shift s;
      const int n = reference.central.GetNrows();
      s.bandUp.ResizeTo(n), s.bandDown.ResizeTo(n), s.covarianceTemplate.ResizeTo(n);
      for (const auto &[label, results] : files) s.templates[label].ResizeTo(n), s.templates[label] = difference(results.at(system), reference, source + "/" + label);
      if (!up.empty()) {
        s.kind = "updown";
        const TVectorD &u = s.templates["up"], &d = s.templates["down"];
        for (int i = 0; i < n; ++i) {
          s.covarianceTemplate[i] = 0.5 * (u[i] - d[i]);
          s.bandUp[i] = std::max({u[i], d[i], 0.0});
          s.bandDown[i] = std::min({u[i], d[i], 0.0});
        }
      } else if (!alternative.empty()) {
        s.kind = "one_sided";
        s.covarianceTemplate = s.templates["alternative"];
        for (int i = 0; i < n; ++i) s.bandUp[i] = std::fabs(s.covarianceTemplate[i]), s.bandDown[i] = -s.bandUp[i];
      } else {
        s.kind = "envelope";
        const std::string representative = setting(env, source + ".Representative");
        if (!s.templates.count(representative)) throw std::runtime_error(source + ": the representative is one of the alternatives");
        s.covarianceTemplate = s.templates[representative];
        for (int i = 0; i < n; ++i) {
          for (const auto &[label, t] : s.templates) s.bandUp[i] = std::max(s.bandUp[i], std::fabs(t[i]));
          s.bandDown[i] = -s.bandUp[i];
        }
      }
      shifts[system].emplace(source, s);
    }
  }
  const auto ratios = words(setting(env, "Systematics.Ratios", false));
  TIter next(env.GetTable());
  while (auto *record = static_cast<TEnvRec *>(next()))
    if (!known.count(record->GetName())) throw std::runtime_error(configPath + ": unknown setting " + record->GetName());

  const std::string partial = output + ".partial";
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot write " + partial);
  for (const auto &[system, reference] : nominal) {
    const int n = reference.central.GetNrows();
    TMatrixD covariance(n, n);
    TVectorD up(n), down(n);
    for (const auto &[source, s] : shifts[system]) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) covariance(i, j) += s.covarianceTemplate[i] * s.covarianceTemplate[j];
        up[i] += s.bandUp[i] * s.bandUp[i];
        down[i] += s.bandDown[i] * s.bandDown[i];
      }
    }
    for (int i = 0; i < n; ++i) up[i] = std::sqrt(up[i]), down[i] = -std::sqrt(down[i]);
    TMatrixD total = covariance + reference.statistical, correlation(n, n);
    const TVectorD sigma = diagonalSqrt(total), statistical = diagonalSqrt(reference.statistical), systematic = diagonalSqrt(covariance);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) correlation(i, j) = sigma[i] > 0 && sigma[j] > 0 ? total(i, j) / (sigma[i] * sigma[j]) : 0.0;
    out->mkdir(system.c_str())->cd();
    for (auto *h : std::vector<TObject *>{hist("xj_per_photon", reference.edges, reference.central, &statistical),
                                          hist("xj_systematic_sigma", reference.edges, systematic), hist("xj_systematic_up", reference.edges, up),
                                          hist("xj_systematic_down", reference.edges, down), hist("systematic_covariance", reference.edges, covariance),
                                          hist("statistical_covariance", reference.edges, reference.statistical),
                                          hist("total_correlation", reference.edges, correlation)})
      h->Write();
    for (const auto &[source, s] : shifts[system])
      for (const auto &[label, t] : s.templates) hist("shift_" + source + "_" + label, reference.edges, t)->Write();
    std::cout << "[FinalAnalysis] systematics " << system << ": " << shifts[system].size() << " sources" << std::endl;
  }
  /* ratios A/B: coherent variation of correlated sources present in both; independent, linearised statistics */
  for (const auto &item : ratios) {
    const auto slash = item.find('/');
    const std::string a = item.substr(0, slash), b = item.substr(slash + 1);
    if (!nominal.count(a) || !nominal.count(b)) throw std::runtime_error("ratio " + item + ": both results must be in the nominal file");
    const Result &ra = nominal.at(a), &rb = nominal.at(b);
    if (ra.edges != rb.edges) throw std::runtime_error("ratio " + item + ": the reported bins differ");
    const int n = ra.central.GetNrows();
    TVectorD central(n);
    for (int i = 0; i < n; ++i) central[i] = ra.central[i] / rb.central[i];
    std::map<std::string, std::pair<TVectorD, TVectorD>> grouped;
    std::set<std::string> sources;
    for (const auto &[s, x] : shifts[a]) sources.insert(s);
    for (const auto &[s, x] : shifts[b]) sources.insert(s);
    for (const auto &source : sources) {
      const bool coherent = correlated[source] && shifts[a].count(source) && shifts[b].count(source);
      for (const auto &[side, system] : {std::pair<const char *, std::string>{"a", a}, {"b", b}}) {
        if (!shifts[system].count(source)) continue;
        const std::string key = coherent ? source : source + "[" + system + "]";
        auto &g = grouped[key];
        g.first.ResizeTo(n), g.second.ResizeTo(n);
        (std::string(side) == "a" ? g.first : g.second) = shifts[system].at(source).covarianceTemplate;
      }
    }
    TMatrixD covariance(n, n), statistical(n, n);
    std::map<std::string, TVectorD> templates;
    for (const auto &[key, g] : grouped) {
      TVectorD d(n);
      for (int i = 0; i < n; ++i) d[i] = (ra.central[i] + g.first[i]) / (rb.central[i] + g.second[i]) - central[i];
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) covariance(i, j) += d[i] * d[j];
      templates[key].ResizeTo(n), templates[key] = d;
    }
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        statistical(i, j) = ra.statistical(i, j) / (rb.central[i] * rb.central[j]) +
                            rb.statistical(i, j) * ra.central[i] * ra.central[j] / (rb.central[i] * rb.central[i] * rb.central[j] * rb.central[j]);
    const TVectorD statSigma = diagonalSqrt(statistical), systSigma = diagonalSqrt(covariance);
    out->mkdir(("ratio_" + a + "_over_" + b).c_str())->cd();
    hist("ratio", ra.edges, central, &statSigma)->Write();
    hist("ratio_systematic_sigma", ra.edges, systSigma)->Write();
    hist("systematic_covariance", ra.edges, covariance)->Write();
    hist("statistical_covariance", ra.edges, statistical)->Write();
    for (const auto &[key, d] : templates) hist("shift_" + key, ra.edges, d)->Write();
  }
  out->cd();
  TObjString(configPath.c_str()).Write("config");
  TObjString(inputHashes.at(nominalPath).get<std::string>().c_str()).Write("nominal_input_sha256");
  TObjString(inputHashes.dump().c_str()).Write("input_sha256");
  out->Close();
  for (const auto &[file,hash] : inputHashes.items()) if (photonjet_render::sha256(file,true)!=hash.get<std::string>())
    throw std::runtime_error("systematics input changed while computing uncertainties");
  if (gSystem->Rename(partial.c_str(), output.c_str()) != 0) throw std::runtime_error("cannot move systematics into place");
  std::cout << "[FinalAnalysis] wrote " << output << std::endl;
  return 0;
}

}  /* namespace fa */
