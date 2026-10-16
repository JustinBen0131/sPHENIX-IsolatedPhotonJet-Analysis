#include "NativeProducerEvidence.h"
/* Systematics: the systematic uncertainties of the per-photon results, from results files that differ by one change each. */

/* final_analysis --systematics config/systematics.yaml */

/* The settings name the nominal results file and every source; a variation is a copy of the analysis settings with one
   change, run like the nominal. Per source and result, the signed shifts (varied minus nominal, per reported xJ bin):
   Up, Down                 a two-sided nuisance: covariance d d^T with d = (up - down)/2; the band takes, per bin, the
   larger positive and the larger negative of the two shifts
   Alternative              one-sided, mirrored: d = alternative - nominal, band +-|d|
   Alternatives, Representative   an envelope: the band is the largest |shift| of the alternatives, the covariance
   d d^T of the named representative
   Components                  independent boundary prescriptions under one displayed source; each component uses
   the rules above, and their covariances and squared bands add. Systems optionally scopes a source/component.
   Sources are independent: covariances add, bands add in quadrature. The statistical covariance is the nominal's own
   (data plus response toys). A ratio A/B varies a source present in both results together when it is Correlated (the
   same physical input), and each system alone otherwise; the two systems' statistics are independent. */
#include "FinalAnalysis.h"
#include "HistogramErrors.h"
#include "DerivedScientificTrace.h"
#include "ConfigurationScientificWitness.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "TEnv.h"
#include "../../TreeToHists/src/Configuration.h"
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
namespace sc=photonjet_science;
namespace ds=derived_science;

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
  sc::Trace centralTrace,statisticalTrace;
};

std::map<std::string, Result> readResults(const std::string &path,const std::string &role,sc::Json &exports) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  observeNativeInput(*f);
  if(auto *saved=dynamic_cast<TObjString*>(f->Get("native_parameter_exports"))) {
    const auto values=sc::Json::parse(saved->GetString().Data());
    for(const auto &[key,value]:values.items()) {
      if(exports.contains(key) && exports.at(key)!=value)throw std::runtime_error("conflicting native parameter export");
      exports[key]=value;
    }
  }
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
    r.centralTrace=sc::traceOf(h).inRole(role).then("reported_xj_bins",{{"edges",r.edges}});
    r.statisticalTrace=ds::join("sum_independent_data_and_response_covariance",
      {sc::traceOf(cd).inRole(role),sc::traceOf(cr).inRole(role)},{{"edges",r.edges}});
    out.emplace(name, r);
  }
  return out;
}

struct Shift {
  std::string kind;                        /* updown, one_sided, envelope */
  std::map<std::string, TVectorD> templates;
  TVectorD covarianceTemplate, bandUp, bandDown;
  std::map<std::string,sc::Trace> templateTraces;
  sc::Trace covarianceTrace,bandTrace;
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
  h->Sumw2(); // Matrix entries have no independent histogram counting error.
  for (int i = 0; i < m.GetNrows(); ++i) for (int j = 0; j < m.GetNcols(); ++j) h->SetBinContent(i + 1, j + 1, m(i, j));
  return h;
}

TVectorD diagonalSqrt(const TMatrixD &m) {
  TVectorD v(m.GetNrows());
  for (int i = 0; i < v.GetNrows(); ++i) v[i] = std::sqrt(std::max(m(i, i), 0.0));
  return v;
}

/* A coverage table distinguishes an unavailable uncertainty from an inapplicable
   one. Every nominal result and named source must be declared; an included
   source is accepted only when its actual variation files contain that result. */
photonjet_render::Json readCoverage(const std::string &file, const std::string &nominalHash,
                                   const std::map<std::string, Result> &nominal,
                                   const std::vector<std::string> &configured) {
  using Json = photonjet_render::Json;
  const auto keys = [](const Json &object) {
    if (!object.is_object()) throw std::runtime_error("systematic coverage requires an object");
    std::set<std::string> names;
    for (const auto &item : object.items()) names.insert(item.key());
    return names;
  };
  std::ifstream input(file);
  if (!input) throw std::runtime_error("cannot read systematic coverage " + file);
  Json table; input >> table;
  if (keys(table) != std::set<std::string>{"schema", "nominal_sha256", "sources", "results"} ||
      table.at("schema") != "PhotonJetSystematicCoverageV1" || table.at("nominal_sha256") != nominalHash)
    throw std::runtime_error("systematic coverage schema or nominal hash differs");
  const auto names = table.at("sources").get<std::vector<std::string>>();
  const std::set<std::string> sources(names.begin(), names.end()), active(configured.begin(), configured.end());
  if (sources.empty() || sources.size() != names.size() || sources.count("") || active.size() != configured.size())
    throw std::runtime_error("systematic coverage sources must be distinct and nonempty");
  for (const auto &source : active)
    if (!sources.count(source)) throw std::runtime_error("systematic coverage omits configured source " + source);
  std::set<std::string> expected;
  for (const auto &[name, result] : nominal) expected.insert(name);
  if (keys(table.at("results")) != expected)
    throw std::runtime_error("systematic coverage result keys differ from nominal");
  std::set<std::string> included;
  for (const auto &[name, result] : nominal) {
    const auto &row = table.at("results").at(name);
    if (keys(row) != sources) throw std::runtime_error("systematic coverage source keys differ for " + name);
    for (const auto &source : sources) {
      const auto &cell = row.at(source);
      if (keys(cell) != std::set<std::string>{"status", "reason"} || !cell.at("reason").is_string() ||
          cell.at("reason").get<std::string>().empty())
        throw std::runtime_error("systematic coverage needs a reason for " + name + "/" + source);
      const auto status = cell.at("status").get<std::string>();
      if (status != "included" && status != "not_yet_included" && status != "not_applicable")
        throw std::runtime_error("unknown systematic coverage status for " + name + "/" + source);
      if (status == "included") {
        if (!active.count(source)) throw std::runtime_error("systematic coverage includes unconfigured source " + source);
        included.insert(source);
      }
    }
  }
  if (included != active) throw std::runtime_error("configured sources must be included in at least one result");
  return table;
}

}  /* namespace */

/* Compare saved nominal and variation results on compatible axes, then build nuisance covariances and ratios.
   Variation production is a separate step; this command does not reread or reprocess eventTrees. */
int systematics(const std::string &configPath) {
  resetNativeInputs();
  TEnv env;
  photonjet_config::readSettings(configPath, env);
  ConfigurationScientificWitness witness;witness.load(configPath);
  sc::Json parameterExports=sc::Json::object();
  auto parameter=[&](sc::Trace trace,const std::string &key,const sc::Json &value) {
    return trace.withParameter("FinalAnalysis.systematics/"+key,"models",value,witness.evaluated(key,value));
  };
  const std::string base = gSystem->GetDirName(configPath.c_str()).Data();
  auto path = [&](const std::string &p) { return p.empty() || p[0] == '/' ? p : base + "/" + p; };
  std::set<std::string> known = {"Systematics.Nominal", "Systematics.Output", "Systematics.Sources", "Systematics.Ratios", "Systematics.Coverage"};
  photonjet_render::Json inputHashes = photonjet_render::Json::object();
  auto pinnedResults = [&](const std::string &file,const std::string &role) {
    inputHashes[file] = photonjet_render::sha256(file,true);
    return readResults(file,role,parameterExports);
  };
  const std::string nominalPath = path(setting(env, "Systematics.Nominal"));
  const auto nominal = pinnedResults(nominalPath,"nominal");
  const auto sources = words(setting(env, "Systematics.Sources"));
  photonjet_render::Json coverage;
  const auto coverageSetting = setting(env, "Systematics.Coverage", false);
  if (!coverageSetting.empty()) {
    const auto coveragePath = path(coverageSetting);
    inputHashes[coveragePath] = photonjet_render::sha256(coveragePath, true);
    coverage = readCoverage(coveragePath, inputHashes.at(nominalPath), nominal, sources);
  }
  const std::string output = path(setting(env, "Systematics.Output"));
  if (!gSystem->AccessPathName(output.c_str())) throw std::runtime_error(output + " exists; choose a new output");
  std::map<std::string, std::map<std::string, Shift>> shifts;  /* system -> source -> shift */
  std::map<std::string, bool> correlated;
  photonjet_render::Json sourceReferences = photonjet_render::Json::object();
  // Components keep separately defined ID boundaries independent without publishing them as unrelated sources.
  std::vector<std::string> nuisanceSources;
  std::map<std::string,std::string> parentSource,correlationKeys;
  photonjet_render::Json sourceGroups=photonjet_render::Json::object();
  for(const auto &source:sources) {
    known.insert(source+".Components");known.insert(source+".Systems");
    const auto components=words(setting(env,source+".Components",false));
    if(components.empty()) {nuisanceSources.push_back(source);parentSource[source]=source;continue;}
    known.insert(source+".Family");known.insert(source+".Correlated");
    const auto family=setting(env,source+".Family"),correlation=setting(env,source+".Correlated");
    std::set<std::string> unique;
    for(const auto &component:components) {
      if(component.empty() || component.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos ||
          !unique.insert(component).second)throw std::runtime_error(source+": invalid/duplicate component name");
      const auto child=source+"."+component;
      if(std::find(sources.begin(),sources.end(),child)!=sources.end())throw std::runtime_error("component is also an independent source: "+child);
      nuisanceSources.push_back(child);parentSource[child]=source;sourceGroups[source].push_back(child);
      correlationKeys[child]=env.Defined((child+".Correlated").c_str())?child+".Correlated":source+".Correlated";
      for(const auto &field:{std::string("Family"),std::string("Correlated"),std::string("Systems")})
        if(!env.Defined((child+"."+field).c_str()) && env.Defined((source+"."+field).c_str()))
          env.SetValue((child+"."+field).c_str(),setting(env,source+"."+field).c_str());
    }
  }
  std::map<std::string,std::set<std::string>> sourceScopes;
  for(const auto &source:nuisanceSources) {
    const auto scoped=words(setting(env,source+".Systems",false));
    for(const auto &system:scoped)
      if(!nominal.count(system) || !sourceScopes[source].insert(system).second)
        throw std::runtime_error(source+": unknown/duplicate result in Systems");
    if(env.Defined((source+".Systems").c_str()) && scoped.empty())throw std::runtime_error(source+": Systems is empty");
  }
  if(!coverage.is_null())for(const auto &[system,result]:nominal)for(const auto &source:sources)
    if(coverage.at("results").at(system).at(source).at("status")=="included") {
      bool covered=false;
      for(const auto &child:nuisanceSources)if(parentSource.at(child)==source)
        covered=covered || sourceScopes[child].empty() || sourceScopes[child].count(system);
      if(!covered)throw std::runtime_error(source+": included result has no component scope: "+system);
    }
  for (const auto &source : nuisanceSources) {
    for (const char *k : {"Family", "Correlated", "Systems", "Up", "Down", "Alternative", "Alternatives", "Representative", "Reference", "ReferenceReason"}) known.insert(source + "." + k);
    setting(env, source + ".Family");
    correlated[source] = setting(env, source + ".Correlated") == "1";
    const auto referenceSetting = setting(env, source + ".Reference", false);
    const auto referenceReason = setting(env, source + ".ReferenceReason", false);
    if (referenceSetting.empty() != referenceReason.empty())
      throw std::runtime_error(source + ": Reference and nonempty ReferenceReason must be supplied together");
    const auto referencePath = referenceSetting.empty() ? nominalPath : path(referenceSetting);
    const auto sourceReference = referenceSetting.empty() ? nominal : pinnedResults(referencePath,"source/"+source+"/reference");
    sourceReferences[source] = {{"path", referencePath}, {"sha256", inputHashes.at(referencePath)},
      {"treatment", referenceSetting.empty() ? "nominal" : "additive_delta_transfer"}, {"reason", referenceReason}};
    const std::string up = setting(env, source + ".Up", false), alternative = setting(env, source + ".Alternative", false),
                      alternatives = setting(env, source + ".Alternatives", false);
    if ((!up.empty()) + (!alternative.empty()) + (!alternatives.empty()) != 1)
      throw std::runtime_error(source + ": give Up and Down, an Alternative, or Alternatives");
    std::map<std::string, std::map<std::string, Result>> files;  /* label -> results */
    if (!up.empty()) files["up"] = pinnedResults(path(up),"source/"+source+"/up"), files["down"] = pinnedResults(path(setting(env, source + ".Down")),"source/"+source+"/down");
    else if (!alternative.empty()) files["alternative"] = pinnedResults(path(alternative),"source/"+source+"/alternative");
    else for (const auto &item : words(alternatives)) {
      const auto eq = item.find('=');
      if (eq == std::string::npos) throw std::runtime_error(source + ".Alternatives: name=path items");
      files[item.substr(0, eq)] = pinnedResults(path(item.substr(eq + 1)),"source/"+source+"/"+item.substr(0,eq));
    }
    for (const auto &[system, reference] : nominal) {
      if(!sourceScopes[source].empty() && !sourceScopes[source].count(system))continue;
      bool present = true, anyPresent = false;
      for (const auto &[label, results] : files) {
        present = present && results.count(system);
        anyPresent = anyPresent || results.count(system);
      }
      if (!coverage.is_null()) {
        const auto status = coverage.at("results").at(system).at(parentSource.at(source)).at("status").get<std::string>();
        if (status != "included") {
          if (anyPresent) throw std::runtime_error("systematic coverage excludes supplied result " + system + "/" + source);
          continue;
        }
        if (!present) throw std::runtime_error("systematic coverage included result is missing: " + system + "/" + source);
      }
      if (!present) continue;  /* the source does not vary this result */
      if (!sourceReference.count(system)) throw std::runtime_error(source + ": reference result is missing: " + system);
      const auto &deltaReference = sourceReference.at(system);
      if (deltaReference.edges != reference.edges)
        throw std::runtime_error(source + ": reference reported bins differ from nominal for " + system);
      Shift s;
      const int n = reference.central.GetNrows();
      s.bandUp.ResizeTo(n), s.bandDown.ResizeTo(n), s.covarianceTemplate.ResizeTo(n);
      for (const auto &[label, results] : files) s.templates[label].ResizeTo(n), s.templates[label] = difference(results.at(system), deltaReference, source + "/" + label);
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
      std::vector<sc::Trace> alternativesTrace;
      for(const auto &[label,results]:files) {
        s.templateTraces[label]=ds::join("variation_minus_source_reference",
          {results.at(system).centralTrace,deltaReference.centralTrace},{{"source",source},{"direction",label}});
        alternativesTrace.push_back(s.templateTraces.at(label));
      }
      s.bandTrace=sc::Trace::combine("systematic_source_band",alternativesTrace,{{"mode",s.kind}});
      if(s.kind=="envelope") {
        const auto representative=setting(env,source+".Representative");
        s.covarianceTrace=parameter(s.templateTraces.at(representative),source+".Representative",representative);
      } else s.covarianceTrace=s.bandTrace;
      s.covarianceTrace=s.covarianceTrace.then("outer_product_nuisance_covariance",{{"mode",s.kind}});
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
    // Bind only the source definitions actually used by this output quantity.
    std::vector<sc::Trace> bandInputs{reference.centralTrace},covInputs{reference.centralTrace};
    for(const auto &[source,shift]:shifts[system]) {bandInputs.push_back(shift.bandTrace);covInputs.push_back(shift.covarianceTrace);}
    const auto bandTrace=sc::Trace::combine("quadrature_of_independent_source_bands",bandInputs);
    const auto covTrace=sc::Trace::combine("sum_independent_source_covariances",covInputs);
    const auto centralTrace=ds::join("central_with_statistical_variance",{reference.centralTrace,reference.statisticalTrace});
    ScientificLineage lineages={{"xj_per_photon",centralTrace},{"xj_systematic_sigma",covTrace},
      {"xj_systematic_up",bandTrace},{"xj_systematic_down",bandTrace},{"systematic_covariance",covTrace},
      {"statistical_covariance",reference.statisticalTrace},
      {"total_correlation",ds::join("total_covariance_to_correlation",{covTrace,reference.statisticalTrace})}};
    const std::map<std::string,sc::Json> covariances={{"xj_per_photon",ds::covarianceRecord(reference.statistical,"nominal data and response covariance")}};
    out->mkdir(system.c_str())->cd();
    for (auto *h : std::vector<TObject *>{hist("xj_per_photon", reference.edges, reference.central, &statistical),
                                          hist("xj_systematic_sigma", reference.edges, systematic), hist("xj_systematic_up", reference.edges, up),
                                          hist("xj_systematic_down", reference.edges, down), hist("systematic_covariance", reference.edges, covariance),
                                          hist("statistical_covariance", reference.edges, reference.statistical),
                                          hist("total_correlation", reference.edges, correlation)})
      ds::write(h,lineages,covariances);
    for (const auto &[source, s] : shifts[system])
      for (const auto &[label, t] : s.templates) {
        auto *h=hist("shift_"+source+"_"+label,reference.edges,t);
        ds::write(h,{{h->GetName(),s.templateTraces.at(label)}});
      }
    /* The integrals of the reported distribution: the mean xJ and the pairs per photon (value x bin width, summed). A
       single unfolded bin follows the unfolding's prior when the xJ resolution exceeds the bin width; these do far
       less. Statistical errors by linear propagation of the statistical covariance. Per source, the integrals of every
       varied result (nominal + shift) minus the nominal's, with the band rules of the bins. */
    auto integrals = [&](const TVectorD &y) {
      double yield = 0, first = 0;
      for (int i = 0; i < n; ++i) {
        const double width = reference.edges[i + 1] - reference.edges[i];
        yield += y[i] * width, first += y[i] * width * 0.5 * (reference.edges[i] + reference.edges[i + 1]);
      }
      TVectorD out(2);
      out[0] = yield != 0 ? first / yield : 0.0, out[1] = yield;
      return out;
    };
    const TVectorD centralIntegrals = integrals(reference.central);
    TVectorD statisticalIntegrals(2), upIntegrals(2), downIntegrals(2);
    for (int q = 0; q < 2; ++q) {
      TVectorD jacobian(n);
      for (int i = 0; i < n; ++i) {
        const double width = reference.edges[i + 1] - reference.edges[i], centre = 0.5 * (reference.edges[i] + reference.edges[i + 1]);
        jacobian[i] = q == 0 ? (centralIntegrals[1] != 0 ? width * (centre - centralIntegrals[0]) / centralIntegrals[1] : 0.0) : width;
      }
      double variance = 0;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) variance += jacobian[i] * reference.statistical(i, j) * jacobian[j];
      statisticalIntegrals[q] = std::sqrt(std::max(variance, 0.0));
    }
    const Edges two = {0, 1, 2};
    auto labelled = [&](TH1D *h) {
      h->GetXaxis()->SetBinLabel(1, "mean_xj"), h->GetXaxis()->SetBinLabel(2, "pairs_per_photon");
      return h;
    };
    for (const auto &[source, s] : shifts[system]) {
      std::map<std::string, TVectorD> moved;
      for (const auto &[label, t] : s.templates) {
        moved[label].ResizeTo(2);
        moved[label] = integrals(reference.central + t) - centralIntegrals;
        auto *h=labelled(hist("integrals_shift_"+source+"_"+label,two,moved[label]));
        ds::write(h,{{h->GetName(),ds::join("shift_of_bin_width_integrals",{reference.centralTrace,s.templateTraces.at(label)},{{"edges",reference.edges}})}});
      }
      for (int q = 0; q < 2; ++q) {
        double high = 0, low = 0;
        if (s.kind == "updown") high = std::max({moved["up"][q], moved["down"][q], 0.0}), low = std::min({moved["up"][q], moved["down"][q], 0.0});
        else for (const auto &[label, d] : moved) high = std::max(high, std::fabs(d[q])), low = -high;
        upIntegrals[q] += high * high, downIntegrals[q] += low * low;
      }
    }
    for (int q = 0; q < 2; ++q) upIntegrals[q] = std::sqrt(upIntegrals[q]), downIntegrals[q] = -std::sqrt(downIntegrals[q]);
    for(auto *h:{labelled(hist("integrals",two,centralIntegrals,&statisticalIntegrals)),
      labelled(hist("integrals_systematic_up",two,upIntegrals)),labelled(hist("integrals_systematic_down",two,downIntegrals))})
      ds::write(h,{{h->GetName(),(std::string(h->GetName())=="integrals"?centralTrace:bandTrace).then("bin_width_integrals",{{"edges",reference.edges}})}});
    std::cout << "[FinalAnalysis] systematics " << system << ": " << shifts[system].size() << " sources; mean xJ " << centralIntegrals[0] << " +- "
              << statisticalIntegrals[0] << " (stat) +" << upIntegrals[0] << " " << downIntegrals[0] << " (syst); pairs per photon " << centralIntegrals[1]
              << " +- " << statisticalIntegrals[1] << " (stat) +" << upIntegrals[1] << " " << downIntegrals[1] << " (syst)" << std::endl;
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
    std::map<std::string,std::vector<sc::Trace>> groupedTraces;
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
        auto trace=shifts[system].at(source).covarianceTrace;
        if(shifts[a].count(source) && shifts[b].count(source))trace=parameter(trace,correlationKeys.count(source)?correlationKeys.at(source):source+".Correlated",correlated[source]);
        groupedTraces[key].push_back(trace);
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
    const auto centralTrace=ds::join("ratio_with_independent_statistics",{ra.centralTrace,rb.centralTrace,ra.statisticalTrace,rb.statisticalTrace});
    std::vector<sc::Trace> ratioInputs{ra.centralTrace,rb.centralTrace};
    std::map<std::string,sc::Trace> ratioShifts;
    for(const auto &[key,inputs]:groupedTraces) {
      auto all=inputs;all.push_back(ra.centralTrace);all.push_back(rb.centralTrace);
      ratioShifts[key]=sc::Trace::combine("correlated_ratio_source_shift",all,{{"source_group",key}});
      ratioInputs.push_back(ratioShifts.at(key));
    }
    const auto systematicTrace=sc::Trace::combine("ratio_independent_source_covariance",ratioInputs);
    out->mkdir(("ratio_" + a + "_over_" + b).c_str())->cd();
    ScientificLineage ratioLineage={{"ratio",centralTrace},{"ratio_systematic_sigma",systematicTrace},
      {"systematic_covariance",systematicTrace},{"statistical_covariance",centralTrace}};
    for(auto *h:std::vector<TObject*>{hist("ratio",ra.edges,central,&statSigma),hist("ratio_systematic_sigma",ra.edges,systSigma),
      hist("systematic_covariance",ra.edges,covariance),hist("statistical_covariance",ra.edges,statistical)})
      ds::write(h,ratioLineage,{{"ratio",ds::covarianceRecord(statistical,"independent numerator and denominator statistics")}});
    for(const auto &[key,d]:templates) {auto *h=hist("shift_"+key,ra.edges,d);ds::write(h,{{h->GetName(),ratioShifts.at(key)}});}
  }
  out->cd();
  TObjString(configPath.c_str()).Write("config");
  TObjString(inputHashes.at(nominalPath).get<std::string>().c_str()).Write("nominal_input_sha256");
  TObjString(inputHashes.dump().c_str()).Write("input_sha256");
  TObjString(sourceReferences.dump().c_str()).Write("systematic_source_references");
  TObjString(sourceGroups.dump().c_str()).Write("systematic_source_groups");
  if (!coverage.is_null()) TObjString(coverage.dump().c_str()).Write("systematic_coverage");
  for(const auto &[key,value]:witness.exports.items())parameterExports[key]=value;
  TObjString(parameterExports.dump().c_str()).Write("native_parameter_exports");
  writeNativeProducerWitness(*out);
  out->cd();
  requireExplicitHistogramErrors(out.get());
  out->Close();
  for (const auto &[file,hash] : inputHashes.items()) if (photonjet_render::sha256(file,true)!=hash.get<std::string>())
    throw std::runtime_error("systematics input changed while computing uncertainties");
  if (gSystem->Rename(partial.c_str(), output.c_str()) != 0) throw std::runtime_error("cannot move systematics into place");
  std::cout << "[FinalAnalysis] wrote " << output << std::endl;
  return 0;
}

}  /* namespace fa */
