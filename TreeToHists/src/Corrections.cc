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

// BEGIN HISTOGRAM INPUT ADAPTER
/* Native DATA input to the existing vertex derive path. This adapter reads
   selected-event histogram packages and original-source cells only. A sealed
   contract pins package bytes, the producing code and the complete ordered
   source roster. It never reads eventTree or changes a population's provenance.
   This block is inserted in Corrections.cc, after Sample and before derive. */

Json vertexObject(TFile& file,const char* name) {
  const auto* object=file.Get<TObjString>(name);
  if(!object) throw std::runtime_error(std::string("vertex input missing ")+name);
  return Json::parse(object->GetString().Data());
}
double vertexUnitCount(double count,double variance) {
  if(!std::isfinite(count) || count<0 || count>9007199254740991. || count!=std::floor(count) || variance!=count)
    throw std::runtime_error("vertex DATA needs finite integer unit counts and equal Sumw2");
  return count;
}
Json vertexSourceIdentity(const Json& s) {
  const auto& metadata=s.at("metadata");const auto& original=metadata.at("source");
  return {{"identity",s.at("identity")},{"sha256",s.at("sha256")},{"metadata_sha256",digest(metadata)},
    {"logical_source_id",original.value("logical_source_id",original.value("input_file_sha256",std::string()))},
    {"run",original.at("run")},{"entries",metadata.at("entries")}};
}

/* The manifest is an exact parent-selected input contract, not a scientific
   qualification. source_records_sha256 hashes the ordered vertexSourceIdentity
   objects. Record the manifest hash on the invocation, not inside that manifest. */
Json importDataHistograms(const std::string& manifest,const std::string& expectedHash,
                          const Configuration& c,TH1D& target,std::set<std::string>& originals) {
  if(sha256File(manifest)!=expectedHash) throw std::runtime_error("DATA histogram contract hash changed");
  const Json contract=readJson(manifest);
  if(contract.at("schema")!="AuAuVertexDataHistogramInputV1" || contract.at("selection")!=c.settings ||
     c.recoBins!=80 || c.recoLow!=-10 || c.recoHigh!=10)
    throw std::runtime_error("DATA histogram selection or vertex binning differs");
  const auto& au=c.settings.at("systems").at("auau");const auto& event=au.at("event_selection");
  if(au.at("data_trigger")!="photon12" || event.at("trigger_name")!=trigger ||
     event.at("trigger_word")!="scaled" || !event.at("data_run_periods").empty() ||
     !event.at("require_minimum_bias").get<bool>() || event.at("vertex_abs_z_max_cm")!=10)
    throw std::runtime_error("DATA target requires nominal Photon12, MB and abs(z)<10");
  if(!contract.at("producer_code").is_object() || contract.at("producer_code").empty() ||
     !contract.at("packages").is_array() || contract.at("packages").empty())
    throw std::runtime_error("DATA histogram contract has no producer or packages");
  long long expectedSources=0;
  for(const auto& pin:contract.at("packages")) {
    const auto count=pin.at("source_count").get<long long>();
    if(count<=0 || count>1000000000LL-expectedSources)
      throw std::runtime_error("DATA histogram source total is invalid");
    expectedSources+=count;
  }
  ReductionProgress importProgress("vertex-data-import",expectedSources+1);
  Json pins=Json::array();std::set<std::string> packagePaths,identities;long long sourceCount=0;
  constexpr const char* key="qa/events/vertex_selected";
  for(const auto& pin:contract.at("packages")) {
    const std::string path=fs::canonical(pin.at("path").get<std::string>()).string();
    const std::string hash=pin.at("sha256");
    if(!packagePaths.insert(path).second || sha256File(path)!=hash)
      throw std::runtime_error("duplicate DATA package or changed package hash");
    fs::path sidecar=path;sidecar.replace_extension(".json");
    if(sha256File(sidecar.string())!=pin.at("receipt_sha256").get<std::string>())
      throw std::runtime_error("DATA package receipt hash changed");
    const Json receipt=readJson(sidecar.string());TFile file(path.c_str(),"READ");
    if(file.IsZombie()) throw std::runtime_error("cannot read DATA histogram package");
    const Json m=vertexObject(file,"meta/package");
    if(m.at("schema")!="PhotonJetSampleV1" || m.at("system")!="AuAu" || m.at("family")!="data" ||
       m.at("package_kind")!="FINAL" || m.at("selection")!=c.settings || m.at("name")!=pin.at("name") ||
       m.at("plan_sha256")!=pin.at("plan_sha256") || m.at("package_binding_sha256")!=pin.at("package_binding_sha256") ||
       m.at("render_selection").at("producer_code")!=contract.at("producer_code") ||
       m.at("data_trigger")!="photon12" || m.at("trigger_selection")!="applied_in_package")
      throw std::runtime_error("DATA package source, selection or producer binding differs");
    if(receipt.at("sha256")!=hash || receipt.at("plan_sha256")!=m.at("plan_sha256") ||
       receipt.at("package_binding_sha256")!=m.at("package_binding_sha256") ||
       receipt.at("selection_sha256")!=c.settings.at("selection_sha256") ||
       receipt.at("sources")!=m.at("source_count") || receipt.at("validation").at("schema")!="PhotonJetPackageValidationV1" ||
       receipt.at("validation").at("status")!="PASS" ||
       receipt.at("validation").at("producer_code_sha256")!=digest(contract.at("producer_code")))
      throw std::runtime_error("DATA package validation receipt differs");
    auto* histogram=file.Get<TH1>(key);auto* sources=file.Get<TTree>("records/sources");
    auto* cells=file.Get<TTree>("records/source_contributions");
    if(!histogram || !sources || !cells || histogram->GetDimension()!=1 || histogram->GetNbinsX()!=80 ||
       histogram->GetSumw2N()!=82 || sources->GetEntries()!=m.at("source_count").get<long long>() ||
       sources->GetEntries()!=pin.at("source_count").get<long long>())
      throw std::runtime_error("missing selected DATA histogram/source records or wrong dimensions");
    for(int b=1;b<=81;++b) if(histogram->GetXaxis()->GetBinLowEdge(b)!=-10.+0.25*(b-1))
      throw std::runtime_error("selected DATA vertex bin edges differ");
    for(int b=0;b<82;++b) vertexUnitCount(histogram->GetBinContent(b),histogram->GetSumw2()->At(b));
    if(histogram->GetBinContent(0)!=0 || histogram->GetBinContent(81)!=0)
      throw std::runtime_error("selected DATA vertex has flow counts");
    auto summed=h1("source_data_reco",80,-10,10);
    const bool aggregate=finalDataAggregatesValidated(receipt,m,hash);
    TTreeReader rs(sources),rc(cells);
    TTreeReaderValue<Long64_t> index(rs,"source"),begin(rs,"contribution_begin"),end(rs,"contribution_end");
    TTreeReaderValue<std::string> json(rs,"json"),name(rc,"histogram");
    TTreeReaderValue<Long64_t> owner(rc,"source"),bin(rc,"bin");
    TTreeReaderValue<double> value(rc,"value"),variance(rc,"variance");
    Long64_t previous=0,ordinal=0;
    const std::string rosterHash=digestStream([&](std::ostream& out) {
      out<<'[';
      while(rs.Next()) {
        const Json s=Json::parse(*json),identity=vertexSourceIdentity(s);
        const std::string id=identity.at("logical_source_id"),recordIdentity=identity.at("identity");
        const long long entries=identity.at("entries");
        if(*index!=ordinal || id.empty() || recordIdentity.empty() || entries<0 ||
           !originals.insert(id).second || !identities.insert(recordIdentity).second)
          throw std::runtime_error("duplicate or invalid DATA original source");
        if(identity.at("run").get<int>()<=0 || identity.at("sha256").get<std::string>().size()!=64)
          throw std::runtime_error("DATA original run or hash missing");
        long long next=0;
        for(const auto& range:s.at("entry_ranges")) {
          if(range.size()!=2 || range.at(0).get<long long>()!=next || range.at(1).get<long long>()<=next)
            throw std::runtime_error("DATA source range gap or overlap");
          next=range.at(1).get<long long>();
        }
        if(next!=entries || s.at("reduced_entries")!=entries || *begin!=previous || *end<*begin || *end>cells->GetEntries())
          throw std::runtime_error("incomplete DATA source or contribution ranges");
        const auto& sm=s.at("metadata");
        const Json table=sm.contains("trigger_selection")?sm.at("trigger_selection"):sm.at("trigger_run_info");
        if(table.at("name").size()!=64 || table.at("initial_prescale").size()!=64)
          throw std::runtime_error("invalid source trigger table");
        int bit=-1;
        for(int b=0;b<64;++b) if(table.at("name").at(b)==trigger) {
          if(bit>=0) throw std::runtime_error("ambiguous source Photon12 trigger");bit=b;
        }
        const bool enabled=bit>=0 && table.at("initial_prescale").at(bit).get<int>()>=0;
        std::set<long long> bins;double selected=0;
        // Full native validation has already checked the aggregate against
        // these source cells. Still inspect disabled-trigger sources for zero
        // selected entries; older or unvalidated inputs keep the full check.
        for(Long64_t row=*begin;(!aggregate || !enabled) && row<*end;++row) {
          if(rc.SetEntry(row)!=TTreeReader::kEntryValid || *owner!=ordinal)
            throw std::runtime_error("DATA contribution belongs to wrong source");
          if(*name!=key) continue;
          if(*bin<1 || *bin>80 || !bins.insert(*bin).second)
            throw std::runtime_error("duplicate or flow DATA source bin");
          const double count=vertexUnitCount(*value,*variance);selected+=count;
          summed->AddBinContent(*bin,count);summed->GetSumw2()->AddAt(summed->GetSumw2()->At(*bin)+count,*bin);
        }
        if(selected>entries || (!enabled && selected!=0))
          throw std::runtime_error("DATA source selected counts contradict entries or disabled trigger");
        if(ordinal) out<<',';out<<identity.dump(-1,' ',true);
        previous=*end;++ordinal;
        importProgress.advance(sourceCount+ordinal);
      }
      out<<']';
    });
    if(ordinal!=sources->GetEntries() || previous!=cells->GetEntries() || rosterHash!=pin.at("source_records_sha256").get<std::string>())
      throw std::runtime_error("DATA exact source roster or contribution coverage differs");
    if(aggregate) {
      for(int b=0;b<82;++b) {
        summed->SetBinContent(b,histogram->GetBinContent(b));
        summed->GetSumw2()->AddAt(histogram->GetSumw2()->At(b),b);
      }
    }
    for(int b=0;b<82;++b) if(summed->GetBinContent(b)!=histogram->GetBinContent(b) ||
                              summed->GetSumw2()->At(b)!=histogram->GetSumw2()->At(b))
      throw std::runtime_error("selected DATA histogram differs from original-source counts or Sumw2");
    if(sha256File(path)!=hash || sha256File(sidecar.string())!=pin.at("receipt_sha256").get<std::string>())
      throw std::runtime_error("DATA input changed during import");
    target.Add(summed.get());sourceCount+=ordinal;
    Json saved=pin;saved["selected_events"]=summed->Integral();
    saved["aggregation"]=aggregate?"validated_final_histogram":"per_source_cells";pins.push_back(saved);
  }
  if(sha256File(manifest)!=expectedHash || sourceCount!=contract.at("source_count").get<long long>())
    throw std::runtime_error("DATA contract changed or total source count differs");
  importProgress.finish();
  return {{"schema","AuAuVertexDataHistogramImportV1"},{"status","STRUCTURE_VERIFIED_NOT_SCIENTIFIC_ACCEPTANCE"},
    {"manifest",manifest},{"manifest_sha256",expectedHash},{"packages",pins},{"source_count",sourceCount},
    {"selected_events",target.Integral()},{"producer_code",contract.at("producer_code")},
    {"selection_sha256",c.settings.at("selection_sha256")},{"histogram",key},{"event_tree_reads",0}};
}

/* Seal selected DATA package metadata for derive without reopening eventTree.
   The importer still verifies the counts and source cells when it consumes
   this contract. Constructing it does not qualify a physical correction. */
void dataHistogramInput(const std::vector<std::string>& args,const Configuration& c) {
  const fs::path out=fs::absolute(option(args,"--out"));
  if(fs::exists(out)) throw std::runtime_error("DATA histogram input already exists");
  Json contract={{"schema","AuAuVertexDataHistogramInputV1"},{"selection",c.settings},
    {"producer_code",Json()},{"packages",Json::array()},{"source_count",0}};
  // Read only the small headers to bind progress before any full-file hash.
  long long expectedTotal=0;
  for(size_t i=0;i<args.size();++i) if(args[i]=="--data") {
    if(i+1==args.size()) throw std::runtime_error("vertex: missing --data path");
    const std::string path=fs::canonical(args[++i]).string();TFile file(path.c_str(),"READ");
    if(file.IsZombie()) throw std::runtime_error("cannot read DATA histogram package");
    const long long count=vertexObject(file,"meta/package").at("source_count");
    if(count<=0 || count>1000000000LL-expectedTotal) throw std::runtime_error("invalid DATA source count");
    expectedTotal+=count;
  }
  if(!expectedTotal) throw std::runtime_error("vertex data-input needs at least one DATA package");
  ReductionProgress progress("vertex-data-input",expectedTotal+1);
  std::set<std::string> paths,originals;
  long long total=0;
  for(size_t i=0;i<args.size();++i) if(args[i]=="--data") {
    if(i+1==args.size()) throw std::runtime_error("vertex: missing --data path");
    const std::string path=fs::canonical(args[++i]).string();
    if(!paths.insert(path).second || paths.size()>256) throw std::runtime_error("duplicate or excessive DATA packages");
    fs::path sidecar=path;sidecar.replace_extension(".json");
    const std::string hash=sha256File(path),receiptHash=sha256File(sidecar.string());
    const Json receipt=readJson(sidecar.string());TFile file(path.c_str(),"READ");
    if(file.IsZombie()) throw std::runtime_error("cannot read DATA histogram package");
    const Json m=vertexObject(file,"meta/package");const auto& code=m.at("render_selection").at("producer_code");
    const long long count=m.at("source_count");
    if(m.at("schema")!="PhotonJetSampleV1" || m.at("system")!="AuAu" || m.at("family")!="data" ||
       m.at("package_kind")!="FINAL" || m.at("selection")!=c.settings || m.at("data_trigger")!="photon12" ||
       m.at("trigger_selection")!="applied_in_package" || count<=0 || count>1000000000LL-total ||
       !code.is_object() || code.empty() || receipt.at("sha256")!=hash || receipt.at("sources")!=count ||
       receipt.at("plan_sha256")!=m.at("plan_sha256") || receipt.at("package_binding_sha256")!=m.at("package_binding_sha256") ||
       receipt.at("selection_sha256")!=c.settings.at("selection_sha256") ||
       receipt.at("validation").at("schema")!="PhotonJetPackageValidationV1" ||
       receipt.at("validation").at("status")!="PASS" || receipt.at("validation").at("producer_code_sha256")!=digest(code))
      throw std::runtime_error("DATA package is incomplete or its selection/receipt differs");
    if(contract.at("producer_code").is_null()) contract["producer_code"]=code;
    else if(contract.at("producer_code")!=code) throw std::runtime_error("DATA package producers differ");
    contract["packages"].push_back({{"path",path},{"sha256",hash},{"receipt_sha256",receiptHash},
      {"name",m.at("name")},{"plan_sha256",m.at("plan_sha256")},
      {"package_binding_sha256",m.at("package_binding_sha256")},{"source_count",count}});
    total+=count;
  }
  if(total!=expectedTotal) throw std::runtime_error("DATA source count changed during input verification");
  long long processed=0;
  for(auto& pin:contract.at("packages")) {
    const std::string path=pin.at("path");TFile file(path.c_str(),"READ");
    auto* sources=file.Get<TTree>("records/sources");
    if(file.IsZombie() || !sources || sources->GetEntries()!=pin.at("source_count").get<long long>())
      throw std::runtime_error("DATA package source records are missing or incomplete");
    TTreeReader reader(sources);TTreeReaderValue<Long64_t> index(reader,"source");
    TTreeReaderValue<std::string> text(reader,"json");long long ordinal=0;
    pin["source_records_sha256"]=digestStream([&](std::ostream& stream) {
      stream<<'[';
      while(reader.Next()) {
        const Json identity=vertexSourceIdentity(Json::parse(*text));
        const std::string original=identity.at("logical_source_id");
        if(*index!=ordinal || original.empty() || !originals.insert(original).second ||
           identity.at("run").get<int>()<=0 || identity.at("entries").get<long long>()<0 ||
           identity.at("sha256").get<std::string>().size()!=64)
          throw std::runtime_error("DATA source identity is duplicated or invalid");
        if(ordinal) stream<<',';stream<<identity.dump(-1,' ',true);++ordinal;
        progress.advance(++processed);
      }
      stream<<']';
    });
    fs::path sidecar=path;sidecar.replace_extension(".json");
    if(ordinal!=sources->GetEntries() || sha256File(path)!=pin.at("sha256").get<std::string>() ||
       sha256File(sidecar.string())!=pin.at("receipt_sha256").get<std::string>())
      throw std::runtime_error("DATA package changed while its source roster was sealed");
  }
  contract["source_count"]=total;fs::create_directories(out.parent_path());
  const fs::path partial=out.string()+".partial";
  if(fs::exists(partial)) throw std::runtime_error("DATA histogram input partial already exists");
  std::ofstream stream(partial);stream<<contract.dump(1)<<'\n';stream.close();
  if(!stream) throw std::runtime_error("cannot write DATA histogram input");
  fs::rename(partial,out);
  progress.finish();
}

/* A changed producer bundle can retain exactly the same unweighted population
   kernel. Compare its saved source bytes with the current scientific slices;
   do not replace producer_code or trust a PASS label. The SIM branch reads no
   photon-ID/isolation selection; a changed full selection hash is permitted
   only after the exact kernel and every other binding field compare equal.
   Unrelated fills/plots are absent from population(). Config parser changes and
   any changed population/ownership/normalization/metadata slice fail closed. */
std::string populationCodeRegion(const std::string& path,const std::string& begin,const std::string& end) {
  std::ifstream in(path);if(!in) throw std::runtime_error("missing archived population code");
  const std::string text((std::istreambuf_iterator<char>(in)),{});
  const std::string start="\n"+begin,stop="\n"+end;
  const auto first=text.find(start),last=first==std::string::npos?std::string::npos:text.find(stop,first+start.size());
  if(first==std::string::npos || last==std::string::npos || text.find(start,first+start.size())!=std::string::npos)
    throw std::runtime_error("ambiguous population implementation slice");
  return text.substr(first,last-first);
}
Json checkPopulationProducer(const Json& m,const Configuration& c,const std::vector<std::string>& args) {
  const Json current=reductionCode();
  const bool sameSelection=m.at("binding")==c.binding;
  if(!sameSelection) {
    Json before=m.at("binding"),after=c.binding;before.erase("selection_sha256");after.erase("selection_sha256");
    if(before!=after || !m.contains("samples") || m.at("samples").empty())
      throw std::runtime_error("population event-selection or binning compatibility is unproved");
    for(const auto& [name,sample]:m.at("samples").items())
      if(sample.at("family")!="photonjet" && sample.at("family")!="inclusivejet")
        throw std::runtime_error("only SIM populations may use the proven event-population selection projection");
  }
  if(sameSelection && m.at("producer_code")==current) return {{"status","EXACT_PRODUCER_AND_BINDING"}};
  const std::string path=option(args,"--population-code"),hash=option(args,"--population-code-sha256");
  if(sha256File(path)!=hash) throw std::runtime_error("population code evidence hash changed");
  const Json proof=readJson(path);
  if(proof.at("schema")!="AuAuVertexPopulationCodeInputV1" || proof.at("producer_code")!=m.at("producer_code"))
    throw std::runtime_error("population code evidence belongs to another producer");
  for(const char* key:{"Configuration.h","TreeProduction/Configuration.h"})
    if(m.at("producer_code").at(key)!=current.at(key))
      throw std::runtime_error("population configuration compatibility is unproved");
  const std::vector<std::array<std::string,3>> regions={
    {{"Corrections.cc","constexpr const char* trigger=","struct Sample {"}},
    {{"PhotonJet.h","struct Ownership {","struct CentralityClass {"}},
    {{"PhotonJet.h","inline std::string metaText(TFile&","/* The identity of a file's production input:"}},
    {{"PhotonJet.h","inline ULong64_t dataTriggerMask(","/* DATA admission is shared"}},
    {{"TreeToHists.cc","Json stitchingNormalization(const std::string &home, const std::string &system, const std::string &name, const std::string &family) {",
      "/* The events of one eventTree file in a sample's stitching denominator:"}}
  };
  Json checked=Json::array();
  for(const auto& r:regions) {
    const Json& pin=proof.at("sources").at(r[0]);const std::string before=pin.at("path"),after=c.home+"/src/"+r[0];
    if(pin.at("sha256")!=m.at("producer_code").at(r[0]) || sha256File(before)!=pin.at("sha256").get<std::string>() ||
       sha256File(after)!=current.at(r[0]).get<std::string>())
      throw std::runtime_error("population code source hash changed");
    const std::string oldBody=populationCodeRegion(before,r[1],r[2]),newBody=populationCodeRegion(after,r[1],r[2]);
    if(oldBody!=newBody) throw std::runtime_error("population implementation compatibility is unproved: "+r[0]);
    checked.push_back({{"source",r[0]},{"begin",r[1]},{"end",r[2]},{"canonical_string_sha256",digest(oldBody)}});
  }
  // These declarations affect branch system detection; they lie outside the function slices.
  for(const std::string declaration:{"enum class System { pp, auau };","using SelectionJson = nlohmann::json;"}) {
    const auto& pin=proof.at("sources").at("PhotonJet.h");
    std::ifstream old(pin.at("path").get<std::string>()),now(c.home+"/src/PhotonJet.h");
    const std::string a((std::istreambuf_iterator<char>(old)),{}),b((std::istreambuf_iterator<char>(now)),{});
    if(a.find(declaration)==std::string::npos || b.find(declaration)==std::string::npos)
      throw std::runtime_error("population type compatibility is unproved");
  }
  if(sha256File(path)!=hash) throw std::runtime_error("population code evidence changed during verification");
  return {{"status",sameSelection?"IDENTICAL_POPULATION_KERNEL_AND_BINDING":"IDENTICAL_SIM_EVENT_POPULATION_KERNEL_AND_BINDING"},
    {"full_selection_identity_equal",sameSelection},{"original_binding",m.at("binding")},{"consumer_binding",c.binding},
    {"selection_projection",sameSelection?"none":"SIM branch ignores full photon-ID/isolation selection; every non-selection binding field is identical"},
    {"evidence",path},{"sha256",hash},
    {"original_producer_code",m.at("producer_code")},{"consumer_producer_code",current},{"checked_regions",checked},
    {"scope","unweighted population selection, ownership, normalization and binning only"}};
}

void checkSimulationSums(const Json& sample,const Configuration& c,const std::string& name,const TH1D& truth,const TH3D& pairs) {
  const std::string family=sample.at("family");
  if((family!="photonjet" && family!="inclusivejet") || sample.at("normalization")!=stitchingNormalization(c.home,"AuAu",name,family))
    throw std::runtime_error("retained SIM family or normalization differs");
  const auto wantedTruth=h1("expected_truth",c.truthBins,c.truthLow,c.truthHigh);const auto wantedPairs=joint(c);
  const std::vector<std::pair<const TAxis*,const TAxis*>> axes={{truth.GetXaxis(),wantedTruth->GetXaxis()},
    {pairs.GetXaxis(),wantedPairs->GetXaxis()},{pairs.GetYaxis(),wantedPairs->GetYaxis()},{pairs.GetZaxis(),wantedPairs->GetZaxis()}};
  for(const auto& axesPair:axes) {
    if(axesPair.first->GetNbins()!=axesPair.second->GetNbins()) throw std::runtime_error("retained SIM axis differs");
    for(int b=1;b<=axesPair.first->GetNbins()+1;++b) if(axesPair.first->GetBinLowEdge(b)!=axesPair.second->GetBinLowEdge(b))
      throw std::runtime_error("retained SIM bin edges differ");
  }
  for(const TH1* h:{static_cast<const TH1*>(&truth),static_cast<const TH1*>(&pairs)}) {
    if(h->GetSumw2N()!=h->GetNcells()) throw std::runtime_error("retained SIM has no unweighted Sumw2");
    for(int b=0;b<h->GetNcells();++b) vertexUnitCount(h->GetBinContent(b),h->GetSumw2()->At(b));
  }
  const long long entries=sample.at("entries"),denominator=sample.at("denominator"),selected=sample.at("selected");
  if(entries<=0 || denominator<=0 || denominator>entries || selected<0 || selected>entries ||
     pairs.Integral()!=selected || truth.Integral(0,truth.GetNbinsX()+1)!=sample.at("truth_moments").at(0).get<double>())
    throw std::runtime_error("retained SIM counts contradict population metadata");
}

// END HISTOGRAM INPUT ADAPTER

double sum(const std::vector<double>& a) {return std::accumulate(a.begin(),a.end(),0.);}
double chi2(const std::vector<double>& d,const std::vector<double>& m) {
  double value=0;for(size_t i=0;i<d.size();++i)
    value+=(d[i]-m[i])*(d[i]-m[i])/((d[i]>0?d[i]:1.)+(m[i]>0?m[i]:1.));return value;
}
void normalize(std::vector<double>& a,double total) {
  const double integral=sum(a);if(!(integral>0 && std::isfinite(integral))) throw std::runtime_error("empty/nonfinite vertex model");
  for(auto& x:a) x*=total/integral;
}

/* Match the DATA vertex selection before forming p(DATA)/p(SIM). Normalize
   each centrality distribution to unity; cross sections and generated-event
   denominators still set the relative simulation mixture, not training weights. */
Json selectedCentrality(const Json& input,const std::map<std::string,Sample>& samples,
                        const Json& sourceFiles) {
  const auto data=input.at("data").at("counts").get<std::vector<double>>();
  if(data.size()!=16) throw std::runtime_error("centrality DATA needs 16 slices");
  std::vector<double> simulation(16,0.);Json sampleInfo=Json::object();
  for(const auto& [name,s]:samples) {
    if(s.metadata.at("family")!="photonjet") continue;
    const double denominator=s.metadata.at("denominator").get<double>();
    const double cross=s.metadata.at("normalization").at("cross_section_pb").get<double>();
    if(!(denominator>0 && cross>0 && std::isfinite(denominator) && std::isfinite(cross)) ||
       !s.pairs || s.pairs->GetNbinsZ()!=16)
      throw std::runtime_error("invalid selected centrality normalization");
    std::vector<double> counts(16,0.);const double factor=cross/denominator;
    for(int z=1;z<=16;++z) for(int x=1;x<=s.pairs->GetNbinsX();++x)
      for(int y=1;y<=s.pairs->GetNbinsY();++y) {
        const double n=s.pairs->GetBinContent(x,y,z);
        if(!std::isfinite(n) || n<0) throw std::runtime_error("invalid selected centrality population");
        counts[z-1]+=n;
      }
    for(int z=0;z<16;++z) simulation[z]+=factor*counts[z];
    sampleInfo[name]={{"counts",counts},{"files",sourceFiles.at(name)},
      {"weight_pb_per_event",factor}};
  }
  const double nd=sum(data),ns=sum(simulation);
  if(!(nd>0 && ns>0 && std::isfinite(nd) && std::isfinite(ns)))
    throw std::runtime_error("empty selected centrality target");
  std::vector<double> weights(16);
  for(int z=0;z<16;++z) {
    if(!(data[z]>0 && simulation[z]>0 && std::isfinite(data[z]) && std::isfinite(simulation[z])))
      throw std::runtime_error("selected centrality has an unsupported slice");
    weights[z]=(data[z]/nd)/(simulation[z]/ns);
  }
  Json result=input;
  result["status"]="DERIVED_UNQUALIFIED";
  result["formula"]="p(DATA selected centrality) / p(SIM selected centrality); unit-normalized shapes";
  result["simulation"]={{"counts",simulation},{"samples",sampleInfo},
    {"histogram","truth_reco_centrality: centrality projection after event/vertex selection"}};
  result["per_slice"]=weights;result["centrality_weights_sha256"]=digest(weights);
  return result;
}

/* Dense source sums bound memory independently of the tree count. A fit
   candidate records binning; exact event-level reweighted closure is still
   required before the production weight receipt can be qualified. */
void derive(const std::vector<std::string>& args,const Configuration& c) {
  const fs::path out=fs::absolute(option(args,"--out"));if(fs::exists(out)) throw std::runtime_error("vertex fit output exists");
  checkReductionCode({{"producer_code",reductionCode()}});
  if(c.centrality.empty()) throw std::runtime_error("vertex derive needs CentralityWeights");
  Json centrality=readJson(c.centrality);auto weights=centrality.at("per_slice").get<std::vector<double>>();
  if(weights.size()!=16 || centrality.at("centrality_weights_sha256").get<std::string>()!=digest(weights) ||
     centrality.at("selection_sha256")!=c.settings.at("selection_sha256")) throw std::runtime_error("centrality weights do not match vertex selection");
  for(double w:weights) if(!(w>0 && std::isfinite(w))) throw std::runtime_error("vertex needs positive finite centrality weights");
  std::map<std::string,Sample> samples;std::set<std::string> sources;Json pins=Json::array();Json sourceFiles=Json::object();
  auto data=h1("data_reco",c.recoBins,c.recoLow,c.recoHigh);
  Json dataImport=nullptr;
  if(std::find(args.begin(),args.end(),"--data-histograms")!=args.end())
    dataImport=importDataHistograms(option(args,"--data-histograms"),option(args,"--data-histograms-sha256"),c,*data,sources);
  for(size_t i=0;i+1<args.size();++i) if(args[i]=="--population") {
    const std::string path=fs::canonical(args[++i]).string(),hash=sha256File(path);TFile file(path.c_str());
    const Json pin=readJson(path+".json");
    if(pin.at("schema").get<std::string>()!="AuAuVertexPopulationReceiptV1" || pin.at("sha256").get<std::string>()!=hash)
      throw std::runtime_error("vertex population changed from its saved receipt");
    auto* object=file.Get<TObjString>("population");if(!object) throw std::runtime_error("missing vertex population metadata");
    Json m=Json::parse(object->GetString().Data());
    if(pin.at("population")!=m) throw std::runtime_error("population receipt metadata differs");
    const std::string claimed=m.at("sha256");m.erase("sha256");
    if(digest(m)!=claimed || pin.at("population").at("sha256").get<std::string>()!=claimed || m.at("schema").get<std::string>()!="AuAuVertexPopulationV1") throw std::runtime_error("vertex population binding changed");
    const Json compatibility=checkPopulationProducer(m,c,args);
    std::map<std::string,long long> sourceEntries;
    for(const auto& s:m.at("sources")) {
      if(s.at("sha256").get<std::string>().size()!=64 || s.at("entries").get<long long>()<0 ||
         !m.at("samples").contains(s.at("sample").get<std::string>())) throw std::runtime_error("invalid population source metadata");
      sourceEntries[s.at("sample").get<std::string>()]+=s.at("entries").get<long long>();
      const std::string sample=s.at("sample");
      if(!sourceFiles.contains(sample)) sourceFiles[sample]=Json::array();
      sourceFiles[sample].push_back({{"path",s.at("path")},{"sha256",s.at("sha256")}});
      const std::string identity=s.value("logical_source_id",std::string());
      if(identity.empty() || !sources.insert(identity).second) throw std::runtime_error("missing or duplicate vertex population original source");
    }
    for(const auto& [name,s]:m.at("samples").items()) {
      if(sourceEntries[name]!=s.at("entries").get<long long>()) throw std::runtime_error("population source counts differ");
      if(s.at("family").get<std::string>()=="data") {
        if(!dataImport.is_null()) throw std::runtime_error("cannot mix histogram DATA target with population DATA");
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
      checkSimulationSums(s,c,name,*t,*h);
      target.pairs->Add(h);target.truth->Add(t);
    }
    if(sha256File(path)!=hash) throw std::runtime_error("vertex population changed during reading");pins.push_back({{"path",path},{"sha256",hash},{"receipt_sha256",sha256File(path+".json")},
      {"original_producer_code",m.at("producer_code")},{"original_binding",m.at("binding")},{"compatibility",compatibility}});
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
  centrality=selectedCentrality(centrality,samples,sourceFiles);
  centrality["data_target_input_snapshot_sha256"]=centrality.at("input_snapshot_sha256");
  centrality["simulation"]["population_files"]=pins;
  centrality["input_snapshot_sha256"]=digest(Json{{"data_package_sha256",centrality.at("data").at("package_sha256")},{"population_files",pins}});
  weights=centrality.at("per_slice").get<std::vector<double>>();
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
  const auto centralityProduct=out/"centrality.json";write(centralityProduct,centrality);
  Json receipt={{"schema","AuAuTruthVertexDerivationV1"},{"status","DERIVED_UNQUALIFIED"},{"binding",c.binding},
    {"method","PPG12 Gaussian-ratio scan followed by smoothed/clipped reco residual applied at truth z"},
    {"fit_settings",c.fit},{"sigma_generated_cm",sigma},{"sigma_target_cm",target},{"stage1_chi2",best},{"stage2_chi2",chi2(d,m2)},
    {"gaussian_scan",scan},{"centrality_path",centralityProduct.string()},{"centrality_sha256",sha256File(centralityProduct.string())},
    {"data_centrality_input",{{"path",c.centrality},{"sha256",sha256File(c.centrality)}}},
    {"centrality_weights_sha256",digest(weights)},{"population_files",pins},{"data_histogram_import",dataImport},{"sample_normalization",sampleInfo},
    {"vertex",{{"path",product.string()},{"sha256",sha256File(product.string())},{"histogram","h_w_iterative"},{"coordinate","truth_vtx_z"}}},
    {"acceptance",{{"binned_fit_only",true},{"exact_event_level_vertex_closure","REQUIRED"},{"centrality_photon_and_inclusive_closure","REQUIRED"},
      {"production_qualified",false}}},{"producer_code",reductionCode()}};
  write(out/"DERIVATION.json",receipt);
  photonjet::progressPhase="vertex_derivation_complete";photonjet::reportProgress(trials+1);photonjet::progressPath.clear();
}
}  // namespace
int vertex(const std::vector<std::string>& args) {
  if(args.empty()) throw std::runtime_error("vertex population|data-input|derive ...");
  const auto c=configuration(args);
  if(args.front()=="population") {population(args,c);return 0;}
  if(args.front()=="data-input") {dataHistogramInput(args,c);return 0;}
  if(args.front()=="derive") {derive(args,c);return 0;}
  throw std::runtime_error("vertex supports population, data-input or derive");
}
}  // namespace treetohists
