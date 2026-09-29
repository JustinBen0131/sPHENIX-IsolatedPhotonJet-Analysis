// Resolved YAML input and event-level validation shared by all three reducers.
// Python selection.py validates YAML and pins packages; no YAML parser is added to ROOT.
#pragma once
#include "EventTreeInput.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <TH1.h>
#include <TTreeReaderValue.h>

namespace photonjet {
using SelectionJson = nlohmann::json;

inline Profile configuredProfile(const char* legacyName, System system, const char* resolvedPath) {
  if (!resolvedPath || !*resolvedPath) return namedProfile(legacyName, system);
  if (std::string(legacyName) != "nominal")
    throw std::runtime_error("use one resolved selection; do not combine it with a legacy profile variant");
  std::ifstream stream(resolvedPath);
  if (!stream) throw std::runtime_error("cannot open resolved selection");
  SelectionJson j; stream >> j;
  if (j.at("schema") != "PhotonJetSelectionV1" || j.at("version") != 1)
    throw std::runtime_error("unsupported resolved selection");
  Profile p = nominal(system);
  const auto& s = j.at("systems").at(system == System::pp ? "pp" : "auau");
  auto line = [](const SelectionJson& v) {
    const double a = v.at("intercept").get<double>(), b = v.at("slope").get<double>();
    if (!std::isfinite(a) || !std::isfinite(b)) throw std::runtime_error("non-finite selection line");
    return Line{a,b};
  };
  p.configured = true;
  p.augmentedClasses = false;
  p.name = j.at("profile").get<std::string>();
  p.selectionId = j.at("selection_sha256").get<std::string>();
  if (p.selectionId.size()!=64 || p.name != "selection:"+p.selectionId)
    throw std::runtime_error("resolved selection identity differs");
  p.selectionJson = j.dump();
  p.vertexAbsZMax = s.at("event_selection").at("vertex_abs_z_max_cm").get<double>();
  if (!std::isfinite(p.vertexAbsZMax) || p.vertexAbsZMax <= 0)
    throw std::runtime_error("invalid event vertex window");
  p.tight = line(s.at("tight"));
  p.nonTightLow = line(s.at("non_tight_low"));
  p.nonTightHigh = line(s.at("non_tight_high"));
  p.scoreBranch = s.at("score_branch").get<std::string>();
  p.optionalID = s.contains("optional_model");
  // The nominal state remains part of every eventTree, but is never used to
  // infer optional applicability. Optional scores encode that with NaN.
  p.stateBranch = p.optionalID ? "pho_bdt_state" : s.at("state_branch").get<std::string>();
  if (p.optionalID) for (const auto& cell : s.at("optional_model").at("cells")) {
    const auto& et = cell.at("et"); const auto& cent = cell.at("centrality");
    p.idCells.push_back(IDCell{et.at(0), et.at(1), cent.at(0), cent.at(1),
                             line(cell.at("tight")), line(cell.at("non_tight_high"))});
  }
  if (p.optionalID && p.idCells.empty()) throw std::runtime_error("optional photon-ID has no routing cells");
  p.verifyNominalID = s.at("verify_nominal_id").get<bool>();
  p.verifyNominalIsolation = s.at("verify_nominal_isolation").get<bool>();
  const auto& iso = s.at("isolation");
  const double radius = iso.at("radius").get<double>();
  if (radius != 0.3 && radius != 0.4) throw std::runtime_error("isolation radius is not a retained witness");
  p.isolationRadiusIndex = radius == 0.3 ? 0 : 1;
  p.isolationSubtraction = iso.at("subtraction_method").get<int>();
  p.isolationFromWitness = iso.at("threshold").is_string() && iso.at("threshold") == "stored";
  if (!p.isolationFromWitness) p.isolationThreshold = line(iso.at("threshold"));
  p.sidebandGap = iso.at("sideband_gap_gev").is_null() ? std::numeric_limits<double>::quiet_NaN()
                 : iso.at("sideband_gap_gev").get<double>();
  if (!std::isnan(p.sidebandGap) && (!std::isfinite(p.sidebandGap) || p.sidebandGap < 0))
    throw std::runtime_error("invalid isolation sideband gap");
  p.isolationCentMin = iso.at("centrality_min").get<double>();
  p.isolationCentMax = iso.at("centrality_max").get<double>();
  for (double x : (system==System::pp ? std::vector<double>{15,35} : std::vector<double>{0,80}))
    if (!(p.nonTightLow.at(x) < p.nonTightHigh.at(x) && p.nonTightHigh.at(x) <= p.tight.at(x)))
      throw std::runtime_error("selection bands overlap or are empty");
  return p;
}

inline Profile validateSelectionInput(TFile& file, Profile p) {
  if (!p.configured) return p;
  const auto resolved = SelectionJson::parse(p.selectionJson);
  const auto& s = resolved.at("systems").at(p.system==System::pp ? "pp" : "auau");
  const std::string key = s.at("score_metadata").get<std::string>();
  const auto id = SelectionJson::parse(metaText(file,key));
  if (p.optionalID) {
    const auto& requested = s.at("optional_model");
    const auto& receipt = id.at(requested.at("model_name").get<std::string>());
    if (receipt.value("status", std::string()) != "PASS" ||
        receipt.at("release_sha256") != requested.at("release_sha256") ||
        receipt.at("score_branch") != requested.at("score_branch") || receipt.at("domain") != requested.at("domain"))
      throw std::runtime_error("optional score receipt differs from selected release/domain");
  } else {
    if (id.at("model_sha256") != s.at("model_sha256") ||
        id.value("freeze_sha256", std::string()) != s.at("photon_id_package_sha256").get<std::string>())
      throw std::runtime_error("saved score differs from selected feature/model package");
    if (p.scoreBranch != "pho_bdt_score" || p.stateBranch != "pho_bdt_state")
      if (id.value("score_branch", std::string()) != p.scoreBranch || id.value("state_branch", std::string()) != p.stateBranch)
        throw std::runtime_error("augmentation receipt does not bind requested score/state branches");
  }
  if (p.verifyNominalIsolation && p.system==System::auau) {
    const auto iso = SelectionJson::parse(metaText(file,"meta/isolation"));
    if (iso.at("freeze_sha256") != s.at("isolation").at("package_sha256"))
      throw std::runtime_error("saved nominal isolation definition differs from the configured package");
  }
  // Older pp deliveries have no filled isolation class; their witness remains
  // sufficient to evaluate PPG12. New pp deliveries additionally prove label parity.
  if (p.system==System::pp && !file.Get("meta/isolation")) p.verifyNominalIsolation=false;
  return p;
}

template <class Values, class Ints>
void validateSelectionEvent(const Profile& p, double centrality, const Values& et, const Values& score, const Ints& state,
                            const Values& cone, const Values& threshold, const Values& sideband, const Ints& method,
                            const Ints& role, const Ints& idClass, const Ints& isoClass) {
  if (!p.configured) return;
  const std::size_t n=et.GetSize();
  if (score.GetSize()!=n || (!p.optionalID && state.GetSize()!=n) || cone.GetSize()!=2*n || threshold.GetSize()!=2*n ||
      sideband.GetSize()!=2*n || method.GetSize()!=2*n || role.GetSize()!=2*n ||
      (p.verifyNominalID && idClass.GetSize()!=n) || (p.verifyNominalIsolation && isoClass.GetSize()!=n))
    throw std::runtime_error("selection branches are not aligned with photons");
  for (std::size_t i=0;i<n;++i) {
    if (p.optionalID && (!std::isnan(score[i]) && (!std::isfinite(score[i]) || score[i] < 0. || score[i] > 1.)))
      throw std::runtime_error("optional score must be finite in [0,1] or NaN outside its domain");
    if (!p.optionalID && (state[i]<0 || state[i]>2 || (state[i]==1 && (!std::isfinite(score[i]) ||
        (p.system==System::auau && !std::isfinite(centrality))))))
      throw std::runtime_error("invalid score applicability/finite-value witness");
    const auto c=photonClass(p,i,et[i],centrality,score,state,cone,threshold,sideband,method,role,idClass,isoClass);
    if (p.verifyNominalID && idClass[i]!=c.idCode)
      throw std::runtime_error("saved nominal photon-ID class disagrees with its pinned score/definition");
    if (p.verifyNominalIsolation && isoClass[i]!=c.isoCode)
      throw std::runtime_error("saved nominal isolation class disagrees with its pinned witness/definition");
  }
}

inline void writeSelection(const Profile& p) {
  if (p.configured) TObjString(p.selectionJson.c_str()).Write("selection");
}

// The public job binds one qualified AuAu fit to its exact selection and
// centrality weights. pp already contains its vertex factor in w_event.
inline std::unique_ptr<TH1> auauVertexFit;
inline void configureVertexFit(const SelectionJson& spec) {
  auauVertexFit.reset();
  if(spec.is_null() || spec.empty()) return;  // historical macro/reference mode
  if(spec.at("coordinate")!="truth_vtx_z") throw std::runtime_error("unsupported vertex-fit coordinate");
  TDirectory::TContext context;
  std::unique_ptr<TFile> file(TFile::Open(spec.at("path").get<std::string>().c_str(),"READ"));
  if(!file || file->IsZombie()) throw std::runtime_error("cannot open AuAu vertex fit");
  auto* h=dynamic_cast<TH1*>(file->Get(spec.at("histogram").get<std::string>().c_str()));
  if(!h || h->GetDimension()!=1 || h->GetNbinsX()<2) throw std::runtime_error("invalid AuAu vertex histogram");
  for(int i=1;i<=h->GetNbinsX();++i)
    if(!std::isfinite(h->GetBinContent(i)) || h->GetBinContent(i)<0.)
      throw std::runtime_error("invalid AuAu vertex-fit weight");
  auauVertexFit.reset(static_cast<TH1*>(h->Clone("auau_vertex_fit")));auauVertexFit->SetDirectory(nullptr);
}

class VertexEventWeight {
  std::unique_ptr<TTreeReaderValue<double>> z_, previous_;
 public:
  VertexEventWeight(TFile& file,TTreeReader& reader) {
    if(!auauVertexFit) return;
    if(fileSystem(file)!=System::auau || !isSimulation(file))
      throw std::runtime_error("additional vertex fit applies only to AuAu SIM");
    if(file.Get("meta/weighting") && SelectionJson::parse(metaText(file,"meta/weighting")).value("vertex_enabled",false))
      throw std::runtime_error("input already includes a vertex weight");
    z_=std::make_unique<TTreeReaderValue<double>>(reader,"truth_vtx_z");
    previous_=std::make_unique<TTreeReaderValue<double>>(reader,"w_vertex");
  }
  double operator()() const {
    if(!z_) return 1.;
    if(**previous_!=1.) throw std::runtime_error("input vertex factor is not unity; refusing double weighting");
    const double z=**z_;
    if(!std::isfinite(z) || z<=-9000.) throw std::runtime_error("missing truth vertex for nominal AuAu weight");
    if(z<auauVertexFit->GetXaxis()->GetXmin() || z>auauVertexFit->GetXaxis()->GetXmax()) return 0.;
    const double w=auauVertexFit->Interpolate(z);
    if(!std::isfinite(w) || w<0.) throw std::runtime_error("invalid interpolated AuAu vertex weight");
    return w;
  }
};

// Exact trigger names and enabled prescales are resolved per source, never by
// a fixed bit index. Absence/disablement is a valid empty stream; missing or
// ambiguous metadata is an error. The mask is applied to scaled event bits.
inline ULong64_t dataTriggerMask(TFile& file, const std::string& key) {
  const std::string name = key=="photon4" ? "Photon 4 GeV + MBD NS >= 1" :
                          key=="photon12" ? "Photon 12 GeV + MBD NS >= 2, vtx < 150 cm" : "";
  if(name.empty()) throw std::runtime_error("unknown nominal DATA trigger: "+key);
  const auto table = SelectionJson::parse(metaText(file,file.Get("meta/trigger_selection") ?
                                                   "meta/trigger_selection" : "meta/trigger_run_info"));
  const auto& names=table.at("name");const auto& prescales=table.at("initial_prescale");
  if(!names.is_array()||names.size()!=64||prescales.size()!=64)
    throw std::runtime_error("trigger table must contain 64 names/prescales");
  int bit=-1;
  for(int i=0;i<64;++i) if(names[i]==name) {
    if(bit!=-1) throw std::runtime_error("ambiguous exact trigger name");bit=i;
  }
  return bit<0||prescales[bit].get<int>()<0 ? 0 : ULong64_t(1)<<bit;
}

// DATA admission is shared by measurement and selected-object QA. Raw trigger
// QA remains inclusive. Run windows are half-open; a declared source period
// restricts the file to that period, otherwise the two supported periods apply.
struct DataEventSelection {
  bool active=false, minimumBias=false;
  ULong64_t mask=0;
  std::vector<std::pair<int,int>> runWindows;
  DataEventSelection(TFile& file, const Profile& p) {
    active=p.configured && !isSimulation(file);
    if(!active) return;
    const auto definition=SelectionJson::parse(p.selectionJson);
    const auto& s=definition.at("systems").at(p.system==System::pp?"pp":"auau");
    const auto& e=s.at("event_selection");
    minimumBias=e.at("require_minimum_bias").get<bool>();
    mask=dataTriggerMask(file,s.at("data_trigger").get<std::string>());
    const auto source=SelectionJson::parse(metaText(file,"meta/source"));
    const auto period=source.value("period",std::string());
    const auto& periods=e.at("data_run_periods");
    if(p.system==System::pp && !period.empty() && !periods.contains(period))
      throw std::runtime_error("unrecognized pp source run period: "+period);
    for(auto it=periods.begin();it!=periods.end();++it)
      if(period.empty() || it.key()==period) runWindows.emplace_back(it.value().at(0),it.value().at(1));
  }
  bool accepts(int run, ULong64_t scaled, int mb) const {
    if(!active) return true;
    if(!(scaled&mask) || (minimumBias && mb!=1)) return false;
    if(runWindows.empty()) return true;
    return std::any_of(runWindows.begin(),runWindows.end(),[&](auto w){return run>=w.first && run<w.second;});
  }
};
}  // namespace photonjet
