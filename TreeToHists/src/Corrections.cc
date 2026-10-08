/* AuAu truth-vertex correction: a bounded tree pass retains the reco/truth
   vertex population; the PPG12 Gaussian scan and residual fit run on those
   saved sums. Ownership and cross-section denominators match TreeToHists.
   Derivation writes a candidate, never a production-qualified freeze. */
#include "TreeToHists.h"
#include "Configuration.h"
#include "PhotonJet.h"

#include <TH1D.h>
#include <TH3D.h>
#include <TObjString.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <limits>
#include <regex>
#include <set>
#include <stdexcept>

namespace treetohists {
Json stitchingNormalization(const std::string&,const std::string&,const std::string&,const std::string&);
photonjet::Ownership ownershipOf(const Json&);
namespace {
namespace fs=std::filesystem;
constexpr const char* trigger="Photon 12 GeV + MBD NS >= 2, vtx < 150 cm";
std::string option(const std::vector<std::string>& args,const std::string& key) {
  for(size_t i=0;i+1<args.size();++i) if(args[i]==key) return args[i+1];
  throw std::runtime_error("vertex: missing "+key);
}
std::string relative(const fs::path& base,const Json& path) {return fs::canonical(base/path.get<std::string>()).string();}
void write(const fs::path& path,const Json& value) {
  std::ofstream out(path);out<<value.dump(1)<<'\n';out.close();
  if(!out) throw std::runtime_error("cannot write "+path.string());
}
struct Configuration {
  Json settings,fit,binding;
  std::string inputs,centrality,home;
  int truthBins=480,recoBins=80;
  double truthLow=-60,truthHigh=60,recoLow=-10,recoHigh=10;
};
Configuration configuration(const std::vector<std::string>& args) {
  const fs::path path=fs::absolute(option(args,"--config")),base=path.parent_path();
  const Json doc=photonjet_config::readYaml(path.string()).at("VertexCorrection");
  const std::set<std::string> known={"Settings","Inputs","CentralityWeights","TruthRangeCM","RecoRangeCM","BinWidthCM",
    "SigmaTargetRangeCM","SigmaTargetStepCM","ResidualSmoothingPasses","ResidualClip","ResidualEdgeClamp"};
  for(const auto& [key,value]:doc.items()) if(!known.count(key)) throw std::runtime_error("unknown VertexCorrection setting: "+key);
  Configuration c;c.home=option(args,"--home");c.inputs=relative(base,doc.at("Inputs"));
  c.settings=settingsSnapshot(relative(base,doc.at("Settings")),c.home).at("selection");
  const auto& au=c.settings.at("systems").at("auau");
  if(au.at("event_selection").at("trigger_name").get<std::string>()!=trigger || !au.at("event_selection").at("require_minimum_bias").get<bool>())
    throw std::runtime_error("AuAu vertex target requires nominal Photon12 and minimum-bias classification");
  const auto tr=doc.value("TruthRangeCM",Json::array({-60.,60.}));
  const auto rr=doc.value("RecoRangeCM",Json::array({-10.,10.}));
  if(tr.size()!=2 || rr.size()!=2) throw std::runtime_error("vertex ranges need two edges");
  c.truthLow=tr.at(0);c.truthHigh=tr.at(1);c.recoLow=rr.at(0);c.recoHigh=rr.at(1);
  const double width=doc.value("BinWidthCM",0.25);
  const double cut=au.at("event_selection").at("vertex_abs_z_max_cm");
  if(!(width>0 && std::isfinite(width) && c.truthLow<0 && c.truthHigh>0 && c.recoLow==-cut && c.recoHigh==cut))
    throw std::runtime_error("vertex binning must cover truth z and exactly match the reco-vertex selection");
  auto bins=[&](double low,double high) {
    const double n=(high-low)/width;
    if(!std::isfinite(n) || n<4 || n>4000 || std::abs(n-std::round(n))>1.e-8) throw std::runtime_error("vertex bin width does not divide its range");
    return static_cast<int>(std::round(n));
  };
  c.truthBins=bins(c.truthLow,c.truthHigh);c.recoBins=bins(c.recoLow,c.recoHigh);
  if(1LL*c.truthBins*c.recoBins*16>16000000) throw std::runtime_error("vertex population exceeds the bounded histogram size");
  c.fit={{"sigma_target_range_cm",doc.value("SigmaTargetRangeCM",Json::array({0.5,30.}))},
    {"sigma_target_step_cm",doc.value("SigmaTargetStepCM",0.125)},
    {"residual_smoothing_passes",doc.value("ResidualSmoothingPasses",1)},
    {"residual_clip",doc.value("ResidualClip",Json::array({0.1,10.}))},
    {"residual_edge_clamp",doc.value("ResidualEdgeClamp",true)}};
  c.centrality=doc.contains("CentralityWeights")?fs::weakly_canonical(base/doc.at("CentralityWeights").get<std::string>()).string():"";
  c.binding={{"selection_sha256",c.settings.at("selection_sha256")},{"data_trigger",trigger},
    {"truth_range_cm",{c.truthLow,c.truthHigh}},{"reco_range_cm",{c.recoLow,c.recoHigh}},
    {"bin_width_cm",width},{"stitching_sha256",sha256File(c.home+"/config/auau_embedded_stitching.json")}};
  return c;
}
std::unique_ptr<TH1D> h1(const char* name,int bins,double low,double high) {
  auto h=std::make_unique<TH1D>(name,name,bins,low,high);h->SetDirectory(nullptr);h->Sumw2();return h;
}
std::unique_ptr<TH3D> joint(const Configuration& c) {
  auto h=std::make_unique<TH3D>("truth_reco_centrality",";truth z [cm];reco z [cm];5% centrality slice",
    c.truthBins,c.truthLow,c.truthHigh,c.recoBins,c.recoLow,c.recoHigh,16,0,16);
  h->SetDirectory(nullptr);h->Sumw2();return h;
}

/* Population sums remain unweighted: derive combines denominator counts
   before applying the canonical sample cross section, once per sample. */
void population(const std::vector<std::string>& args,const Configuration& c) {
  const fs::path out=fs::absolute(option(args,"--out"));if(fs::exists(out)) throw std::runtime_error("vertex population exists");
  const Json snapshot=freezeInputs(c.inputs);
  long long total=0;
  std::set<std::string> identities;
  for(const auto& sample:snapshot.at("samples")) if(sample.at("system")=="AuAu") for(const auto& pin:sample.at("files")) {
    const std::string path=pin.at("path"),hash=pin.at("sha256");
    if(sha256File(path)!=hash) throw std::runtime_error("changed vertex population source");
    TFile file(path.c_str());auto* tree=file.Get<TTree>("eventTree");
    const Json source=Json::parse(photonjet::metaText(file,"meta/source"));
    const std::string identity=source.value("logical_source_id",source.value("input_file_sha256",std::string()));
    if(identity.empty() || !identities.insert(identity).second) throw std::runtime_error("missing or duplicate vertex population original source");
    if(!tree) throw std::runtime_error("vertex population source has no eventTree");total+=tree->GetEntries();
  }
  if(total<=0) throw std::runtime_error("empty AuAu vertex population");
  fs::create_directories(out.parent_path());
  photonjet::progressPath=out.string()+".progress.json";photonjet::progressTotal=total;photonjet::progressOffset=0;
  photonjet::progressStart=std::chrono::steady_clock::now();photonjet::progressPhase="vertex_population";
  photonjet::entryBegin=0;photonjet::progressIdentity={{"campaign",digest(snapshot)},{"row","VERTEX_POPULATION"},
    {"sample","AuAu"},{"source",digest(snapshot)}};
  Json metadata={{"schema","AuAuVertexPopulationV1"},{"binding",c.binding},{"samples",Json::object()},
    {"input_snapshot_sha256",digest(snapshot)},{"publication_bindings",snapshot.at("publication_bindings")},
    {"sources",Json::array()},{"processed_entries",0},{"producer_code",reductionCode()}};
  const fs::path partial=out.string()+".partial";
  TFile result(partial.c_str(),"CREATE");if(result.IsZombie()) throw std::runtime_error("cannot create vertex population");
  for(const auto& sample:snapshot.at("samples")) {
    if(sample.at("system")!="AuAu") continue;
    const std::string family=sample.at("family"),name=sample.at("name");const bool data=family=="data";
    if(!data && family!="photonjet" && family!="inclusivejet") throw std::runtime_error("invalid vertex population family");
    if(!std::regex_match(name,std::regex("[A-Za-z0-9_-]+"))) throw std::runtime_error("invalid vertex sample name");
    Json norm=data?Json():stitchingNormalization(c.home,"AuAu",name,family);
    const auto ownership=data?photonjet::Ownership{}:ownershipOf(norm);
    auto reco=h1("data_reco",c.recoBins,c.recoLow,c.recoHigh);
    auto truth=h1("truth_z",c.truthBins,c.truthLow,c.truthHigh);
    auto pairs=joint(c);
    long long denominator=0,selected=0,missing=0,outside=0,entries=0;
    double truthCount=0,truthSum=0,truthSquareSum=0;
    for(const auto& pin:sample.at("files")) {
      const std::string path=pin.at("path"),hash=pin.at("sha256");
      TFile file(path.c_str());if(file.IsZombie() || photonjet::fileSystem(file)!=photonjet::System::auau || photonjet::isSimulation(file)==data)
        throw std::runtime_error("vertex source system/family differs from catalogue");
      const Json source=Json::parse(photonjet::metaText(file,"meta/source"));
      const std::string identity=source.value("logical_source_id",source.value("input_file_sha256",std::string()));
      if(!data && source.value("dataset",std::string())!=(family=="photonjet"?"isSimEmbedded":"isSimEmbeddedInclusive"))
        throw std::runtime_error("vertex source simulation family differs from catalogue");
      auto* tree=file.Get<TTree>("eventTree");if(!tree) throw std::runtime_error("missing vertex population tree");
      photonjet::progressRowTotal=tree->GetEntries();photonjet::entryEnd=tree->GetEntries();
      photonjet::progressIdentity["sample"]="AuAu/"+name;photonjet::progressIdentity["source"]=identity;photonjet::reportProgress(0);
      TTreeReader r(tree);
      TTreeReaderValue<int> record(r,"evt_record"),ready(r,"evt_nodes_ready"),terminal(r,"evt_terminal_status"),valid(r,"vtx_valid"),mb(r,"cent_mb_decision");
      TTreeReaderValue<double> z(r,"vtx_z"),cent(r,"cent_value");
      std::unique_ptr<TTreeReaderValue<ULong64_t>> scaled;
      std::unique_ptr<TTreeReaderValue<double>> tz;
      ULong64_t mask=0;
      if(data) {scaled=std::make_unique<TTreeReaderValue<ULong64_t>>(r,"trig_scaled");mask=photonjet::dataTriggerMask(file,c.settings.at("systems").at("auau").at("data_trigger").get<std::string>());}
      else tz=std::make_unique<TTreeReaderValue<double>>(r,"truth_vtx_z");
      photonjet::OwnershipGate owner(r,ownership,true);
      long long seen=0;
      while(r.Next()) {
        ++seen;++entries;photonjet::reportProgress(seen);
        if(!data && (norm.at("ownership").at("denominator").get<std::string>()=="all_events" || ownership.contains(owner.value()))) ++denominator;
        if(*record!=0 || *ready!=1 || *terminal!=0 || *valid!=1 || !std::isfinite(*z) || !std::isfinite(*cent) || *cent<0 || *cent>=80) continue;
        if(data) {
          if(*mb==1 && (**scaled&mask) && c.recoLow<*z && *z<c.recoHigh) {reco->Fill(*z);++selected;}
          continue;
        }
        if(!owner.owned()) continue;
        const double t=**tz;
        if(!std::isfinite(t) || t<=-9000.) {++missing;continue;}
        ++truthCount;truthSum+=t;truthSquareSum+=t*t;
        truth->Fill(t);
        if(t<c.truthLow || t>=c.truthHigh) {++outside;continue;}
        if(c.recoLow<*z && *z<c.recoHigh) {pairs->Fill(t,*z,std::floor(*cent/5)+0.5);++selected;}
      }
      if(seen!=tree->GetEntries() || sha256File(path)!=hash) throw std::runtime_error("vertex population read incomplete or input changed");
      metadata["sources"].push_back({{"sample",name},{"path",path},{"sha256",hash},{"logical_source_id",identity},{"entries",seen}});
      photonjet::progressOffset+=seen;
    }
    auto* directory=result.mkdir(name.c_str());directory->cd();
    if(data) reco->Write("data_reco");else {truth->Write("truth_z");pairs->Write("truth_reco_centrality");}
    metadata["samples"][name]={{"family",family},{"normalization",norm},{"denominator",denominator},{"selected",selected},
      {"missing_truth_vertex",missing},{"truth_outside_range",outside},{"entries",entries},
      {"truth_moments",{truthCount,truthSum,truthSquareSum}}};
  }
  metadata["processed_entries"]=total;metadata["sha256"]=digest(metadata);
  result.cd();TObjString(metadata.dump().c_str()).Write("population");result.Close();
  for(const auto& pin:snapshot.at("publication_bindings")) if(sha256File(pin.at("path"))!=pin.at("sha256").get<std::string>())
    throw std::runtime_error("publication manifest changed during vertex population");
  fs::rename(partial,out);photonjet::progressPhase="vertex_population_complete";
  photonjet::progressOffset=total-photonjet::progressRowTotal;
  photonjet::reportProgress(photonjet::progressRowTotal);photonjet::progressPath.clear();
  write(out.string()+".json",{{"schema","AuAuVertexPopulationReceiptV1"},{"status","POPULATION_ONLY"},
    {"path",out.string()},{"sha256",sha256File(out.string())},{"population",metadata}});
}

struct Sample {
  Json metadata;
  std::unique_ptr<TH1D> truth;
  std::unique_ptr<TH3D> pairs;
};
double sum(const std::vector<double>& a) {return std::accumulate(a.begin(),a.end(),0.);}
double chi2(const std::vector<double>& d,const std::vector<double>& m) {
  double value=0;for(size_t i=0;i<d.size();++i)
    value+=(d[i]-m[i])*(d[i]-m[i])/((d[i]>0?d[i]:1.)+(m[i]>0?m[i]:1.));return value;
}
void normalize(std::vector<double>& a,double total) {
  const double integral=sum(a);if(!(integral>0 && std::isfinite(integral))) throw std::runtime_error("empty/nonfinite vertex model");
  for(auto& x:a) x*=total/integral;
}

/* Dense source sums bound memory independently of the tree count. A fit
   candidate records binning; exact event-level reweighted closure is still
   required before the production weight receipt can be qualified. */
void derive(const std::vector<std::string>& args,const Configuration& c) {
  const fs::path out=fs::absolute(option(args,"--out"));if(fs::exists(out)) throw std::runtime_error("vertex fit output exists");
  if(c.centrality.empty()) throw std::runtime_error("vertex derive needs CentralityWeights");
  const Json centrality=readJson(c.centrality);const auto weights=centrality.at("per_slice").get<std::vector<double>>();
  if(weights.size()!=16 || centrality.at("centrality_weights_sha256").get<std::string>()!=digest(weights) ||
     centrality.at("selection_sha256")!=c.settings.at("selection_sha256")) throw std::runtime_error("centrality weights do not match vertex selection");
  for(double w:weights) if(!(w>0 && std::isfinite(w))) throw std::runtime_error("vertex needs positive finite centrality weights");
  std::map<std::string,Sample> samples;std::set<std::string> sources;Json pins=Json::array();
  auto data=h1("data_reco",c.recoBins,c.recoLow,c.recoHigh);
  for(size_t i=0;i+1<args.size();++i) if(args[i]=="--population") {
    const std::string path=fs::canonical(args[++i]).string(),hash=sha256File(path);TFile file(path.c_str());
    const Json pin=readJson(path+".json");
    if(pin.at("schema").get<std::string>()!="AuAuVertexPopulationReceiptV1" || pin.at("sha256").get<std::string>()!=hash)
      throw std::runtime_error("vertex population changed from its saved receipt");
    auto* object=file.Get<TObjString>("population");if(!object) throw std::runtime_error("missing vertex population metadata");
    Json m=Json::parse(object->GetString().Data());const std::string claimed=m.at("sha256");m.erase("sha256");
    if(digest(m)!=claimed || pin.at("population").at("sha256").get<std::string>()!=claimed || m.at("schema").get<std::string>()!="AuAuVertexPopulationV1" || m.at("binding")!=c.binding ||
       m.at("producer_code")!=reductionCode()) throw std::runtime_error("vertex population binding changed");
    for(const auto& s:m.at("sources")) {
      const std::string identity=s.value("logical_source_id",std::string());
      if(identity.empty() || !sources.insert(identity).second) throw std::runtime_error("missing or duplicate vertex population original source");
    }
    for(const auto& [name,s]:m.at("samples").items()) {
      if(s.at("family").get<std::string>()=="data") {
        auto* h=file.Get<TH1D>((name+"/data_reco").c_str());
        if(!h || h->GetNbinsX()!=c.recoBins || h->GetXaxis()->GetXmin()!=c.recoLow || h->GetXaxis()->GetXmax()!=c.recoHigh)
          throw std::runtime_error("invalid DATA vertex population");data->Add(h);continue;
      }
      auto& target=samples[name];
      if(!target.pairs) {target.pairs=joint(c);target.truth=h1("truth_z",c.truthBins,c.truthLow,c.truthHigh);target.metadata=s;}
      else {
        if(target.metadata.at("normalization")!=s.at("normalization") || target.metadata.at("family")!=s.at("family")) throw std::runtime_error("vertex sample normalization changed");
        for(const char* key:{"denominator","selected","missing_truth_vertex","truth_outside_range","entries"})
          target.metadata[key]=target.metadata.at(key).get<long long>()+s.at(key).get<long long>();
        for(int k=0;k<3;++k) target.metadata["truth_moments"][k]=target.metadata.at("truth_moments").at(k).get<double>()+s.at("truth_moments").at(k).get<double>();
      }
      auto* h=file.Get<TH3D>((name+"/truth_reco_centrality").c_str());auto* t=file.Get<TH1D>((name+"/truth_z").c_str());
      if(!h || !t || h->GetNbinsX()!=c.truthBins || h->GetNbinsY()!=c.recoBins || h->GetNbinsZ()!=16 ||
         h->GetXaxis()->GetXmin()!=c.truthLow || h->GetXaxis()->GetXmax()!=c.truthHigh ||
         h->GetYaxis()->GetXmin()!=c.recoLow || h->GetYaxis()->GetXmax()!=c.recoHigh ||
         t->GetNbinsX()!=c.truthBins || t->GetXaxis()->GetXmin()!=c.truthLow || t->GetXaxis()->GetXmax()!=c.truthHigh)
        throw std::runtime_error("invalid SIM vertex population");
      target.pairs->Add(h);target.truth->Add(t);
    }
    if(sha256File(path)!=hash) throw std::runtime_error("vertex population changed during reading");pins.push_back({{"path",path},{"sha256",hash}});
  }
  if(pins.empty() || samples.empty() || !(data->Integral()>0)) throw std::runtime_error("vertex fit needs DATA and simulation populations");
  auto truth=h1("truth_z",c.truthBins,c.truthLow,c.truthHigh);
  auto pairs=joint(c);auto inclusive=joint(c);Json sampleInfo=Json::object();
  double total=0,mean=0,second=0;
  for(auto& [name,s]:samples) {
    const long long denominator=s.metadata.at("denominator");
    if(denominator<=0 || s.metadata.at("missing_truth_vertex").get<long long>() || s.metadata.at("truth_outside_range").get<long long>())
      throw std::runtime_error("vertex sample has no denominator or missing/out-of-range truth vertices: "+name);
    const double factor=s.metadata.at("normalization").at("cross_section_pb").get<double>()/denominator;
    if(s.metadata.at("family").get<std::string>()=="photonjet") {
      truth->Add(s.truth.get(),factor);pairs->Add(s.pairs.get(),factor);
      const auto& moments=s.metadata.at("truth_moments");
      total+=factor*moments.at(0).get<double>();mean+=factor*moments.at(1).get<double>();second+=factor*moments.at(2).get<double>();
    }
    else inclusive->Add(s.pairs.get(),factor);
    sampleInfo[name]=s.metadata;sampleInfo[name]["pb_per_generated_event"]=factor;
  }
  if(!(total>0)) throw std::runtime_error("no photon-family truth vertex population");
  const double sigma=std::sqrt(std::max(0.,second/total-(mean/total)*(mean/total)));
  if(!(sigma>0 && std::isfinite(sigma))) throw std::runtime_error("degenerate truth-vertex distribution");
  std::vector<double> d(c.recoBins);for(int b=1;b<=c.recoBins;++b)d[b-1]=data->GetBinContent(b);
  auto model=[&](const TH3D& joint,const TH1D& fit) {
    std::vector<double> m(c.recoBins);
    for(int x=1;x<=c.truthBins;++x) for(int y=1;y<=c.recoBins;++y) for(int z=1;z<=16;++z)
      m[y-1]+=joint.GetBinContent(x,y,z)*weights[z-1]*fit.Interpolate(joint.GetXaxis()->GetBinCenter(x));
    return m;
  };
  auto fit=h1("h_w_iterative",c.truthBins,c.truthLow,c.truthHigh);
  auto parametric=[&](double target) {for(int b=1;b<=c.truthBins;++b) {const double z=fit->GetBinCenter(b);
    const double value=std::exp(0.5*z*z*(1/(sigma*sigma)-1/(target*target)));
    if(!std::isfinite(value)) throw std::runtime_error("vertex Gaussian scan overflows; narrow its configured range");fit->SetBinContent(b,value);}};
  const auto range=c.fit.at("sigma_target_range_cm");const double lo=range.at(0),hi=range.at(1),step=c.fit.at("sigma_target_step_cm");
  if(range.size()!=2 || !(lo>0 && lo<hi && step>0 && std::isfinite(step)) || (hi-lo)/step>2000) throw std::runtime_error("invalid bounded vertex sigma scan");
  fs::create_directories(out);photonjet::progressPath=(out/"progress.json").string();
  const long long trials=static_cast<long long>(std::floor((hi-lo+step/100)/step))+1;
  photonjet::entryBegin=0;photonjet::entryEnd=trials+1;photonjet::progressTotal=photonjet::progressRowTotal=trials+1;photonjet::progressOffset=0;
  photonjet::progressStart=std::chrono::steady_clock::now();photonjet::progressPhase="vertex_gaussian_scan";
  photonjet::progressIdentity={{"campaign",digest(pins)},{"row","VERTEX_FIT"},
    {"sample","AuAu/photonjet"},{"source",digest(pins)}};photonjet::reportProgress(0);
  double best=std::numeric_limits<double>::infinity(),target=0;Json scan=Json::array();
  for(double s=lo;s<=hi+step/100;s+=step) {parametric(s);auto m=model(*pairs,*fit);normalize(m,sum(d));const double score=chi2(d,m);
    scan.push_back({{"sigma_target_cm",s},{"chi2",score}});if(score<best) {best=score;target=s;}photonjet::reportProgress(scan.size());}
  photonjet::progressPhase="vertex_residual";photonjet::reportProgress(trials);
  // PPG12 refines an interior grid minimum with its three-point parabola.
  const auto bestIt=std::min_element(scan.begin(),scan.end(),[](const Json& a,const Json& b) {return a.at("chi2").get<double>()<b.at("chi2").get<double>();});
  const size_t index=static_cast<size_t>(std::distance(scan.begin(),bestIt));
  if(index>0 && index+1<scan.size()) {
    const double left=scan.at(index-1).at("chi2"),right=scan.at(index+1).at("chi2"),curvature=left-2*best+right;
    if(curvature>0) {const double refined=target+step*(left-right)/(2*curvature);
      if(refined>target-step && refined<target+step) target=refined;}
  }
  parametric(target);auto m1=model(*pairs,*fit);normalize(m1,sum(d));best=chi2(d,m1);
  auto residual=h1("residual",c.recoBins,c.recoLow,c.recoHigh);std::vector<double> values(c.recoBins,1.);
  for(int b=0;b<c.recoBins;++b) {
    if(m1[b]<=0 && d[b]>0) throw std::runtime_error("DATA reco vertex bin has no simulated support");
    if(m1[b]>0) values[b]=d[b]/m1[b];
  }
  const int smoothing=c.fit.at("residual_smoothing_passes");if(smoothing<0 || smoothing>20) throw std::runtime_error("invalid residual smoothing passes");
  for(int pass=0;pass<smoothing;++pass) {auto copy=values;for(int b=0;b<c.recoBins;++b)
    values[b]=0.25*copy[std::max(0,b-1)]+0.5*copy[b]+0.25*copy[std::min(c.recoBins-1,b+1)];}
  const auto clip=c.fit.at("residual_clip");if(clip.size()!=2 || !(clip.at(0).get<double>()>0 && clip.at(0).get<double>()<=clip.at(1).get<double>())) throw std::runtime_error("invalid residual clip");
  for(auto& v:values) v=std::clamp(v,clip.at(0).get<double>(),clip.at(1).get<double>());
  if(c.fit.at("residual_edge_clamp").get<bool>()) {values.front()=values[1];values.back()=values[c.recoBins-2];}
  for(int b=1;b<=c.recoBins;++b) residual->SetBinContent(b,values[b-1]);
  // PPG12 applies the reco residual as a function of truth z. Outside the
  // analysis interval it is unity; within it interpolation clamps at centers.
  for(int b=1;b<=c.truthBins;++b) {const double z=fit->GetBinCenter(b);if(z<=c.recoLow || z>=c.recoHigh) continue;
    const double clamped=std::clamp(z,residual->GetBinCenter(1),residual->GetBinCenter(c.recoBins));fit->SetBinContent(b,fit->GetBinContent(b)*residual->Interpolate(clamped));}
  auto m2=model(*pairs,*fit);normalize(m2,sum(d));
  auto before=h1("model_stage1",c.recoBins,c.recoLow,c.recoHigh),after=h1("model_stage2",c.recoBins,c.recoLow,c.recoHigh);
  for(int b=1;b<=c.recoBins;++b) {before->SetBinContent(b,m1[b-1]);after->SetBinContent(b,m2[b-1]);}
  fs::create_directories(out);const auto product=out/"vertex.root";TFile result(product.c_str(),"CREATE");
  fit->Write();residual->Write();data->Write();truth->Write();before->Write();after->Write();result.Close();
  Json receipt={{"schema","AuAuTruthVertexDerivationV1"},{"status","DERIVED_UNQUALIFIED"},{"binding",c.binding},
    {"method","PPG12 Gaussian-ratio scan followed by smoothed/clipped reco residual applied at truth z"},
    {"fit_settings",c.fit},{"sigma_generated_cm",sigma},{"sigma_target_cm",target},{"stage1_chi2",best},{"stage2_chi2",chi2(d,m2)},
    {"gaussian_scan",scan},{"centrality_path",c.centrality},{"centrality_sha256",sha256File(c.centrality)},
    {"centrality_weights_sha256",digest(weights)},{"population_files",pins},{"sample_normalization",sampleInfo},
    {"vertex",{{"path",product.string()},{"sha256",sha256File(product.string())},{"histogram","h_w_iterative"},{"coordinate","truth_vtx_z"}}},
    {"acceptance",{{"binned_fit_only",true},{"exact_event_level_vertex_closure","REQUIRED"},{"centrality_photon_and_inclusive_closure","REQUIRED"},
      {"production_qualified",false}}},{"producer_code",reductionCode()}};
  write(out/"DERIVATION.json",receipt);
  photonjet::progressPhase="vertex_derivation_complete";photonjet::reportProgress(trials+1);photonjet::progressPath.clear();
}
}  // namespace
int vertex(const std::vector<std::string>& args) {
  if(args.empty()) throw std::runtime_error("vertex population|derive ...");
  const auto c=configuration(args);
  if(args.front()=="population") {population(args,c);return 0;}
  if(args.front()=="derive") {derive(args,c);return 0;}
  throw std::runtime_error("vertex supports population or derive");
}
}  // namespace treetohists
