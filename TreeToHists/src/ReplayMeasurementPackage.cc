// Source-aware post-read measurement, response and physics-purity replay.
// Numerical fixture admission is explicit; production physics acceptance is separate.
#define main treetohists_original_main
#include "TreeToHists.cc"
#undef main
#include "CandidateResponseReplay.h"
#include "CandidateCachePackage.h"
#include "SelectionReplay.h"
#include "ResponseOnlyResolution.h"
#include "CentralityReplay.h"
#include <TMemFile.h>
#include "ReplayBuildInfo.h"

namespace replay_package {
using namespace photonjet;
using namespace treetohists;
using Range=std::pair<Long64_t,Long64_t>;
using Cell=std::pair<double,double>;
using Cells=std::map<std::pair<std::string,Long64_t>,Cell>;

Json object(TFile& f,const char* name) {
  auto* p=f.Get<TObjString>(name);
  if(!p)throw std::runtime_error(std::string("missing ")+name);
  return Json::parse(p->GetString().Data());
}
TTree& table(TFile& f,const std::string& name) {
  auto* p=f.Get<TTree>(name.c_str());
  if(!p)throw std::runtime_error("missing "+name);
  return *p;
}
// One linear scan per table, never one full scan per source.
std::vector<Range> sourceRanges(TTree& t,size_t count) {
  std::vector<Range> ranges(count,{0,0}); Long64_t previous=-1;
  t.SetBranchStatus("*",false);t.SetBranchStatus("source",true);
  auto* leaf=t.GetLeaf("source");if(!leaf)throw std::runtime_error("missing source leaf");
  for(Long64_t i=0;i<t.GetEntries();++i) {
    if(t.GetEntry(i)<=0)throw std::runtime_error("unreadable source row");
    const auto s=leaf->GetValueLong64();
    if(s<previous || s<0 || s>=static_cast<Long64_t>(count))throw std::runtime_error("unordered or foreign source");
    if(s!=previous)ranges[s].first=i;
    ranges[s].second=i+1;previous=s;
    if(i%10000==0)ReductionProgress::pulse("indexing_candidate_sources");
  }
  t.SetBranchStatus("*",true);return ranges;
}

Json consumerCode() {
  const auto expected=Json::parse(TREETOHISTS_REPLAY_COMPILED_CODE);Json observed=Json::object();
  const fs::path here=fs::path(__FILE__).parent_path().parent_path().parent_path();
  for(const auto& [name,value]:expected.items())observed[name]=sha256File((here/name).string());
  if(observed!=expected)throw std::runtime_error("consumer source differs from compiled code");
  return observed;
}
Json counts(const energy_selection::Counters& c){return {{"candidates",c.candidates},{"incoming_id",c.incomingId},{"outgoing_id",c.outgoingId},
  {"incoming_npb",c.incomingNpb},{"outgoing_npb",c.outgoingNpb},{"preselection_changed",c.preselectionChanged}};}
TDirectory* directory(TFile& f,const std::string& name) {
  if(!f.GetDirectory(name.c_str()))f.mkdir(name.c_str(),"",true);
  auto* d=f.GetDirectory(name.c_str());if(!d)throw std::runtime_error("missing output directory "+name);return d;
}
int run(const Json& request) {
  const std::string inputPath=request.at("package"),target=request.at("out");
  const bool admissionOnly=request.value("admission_only",false);
  const fs::path partial=target+".partial.root";
  fs::path side=target;side.replace_extension(".json");
  if(fs::exists(target)||fs::exists(partial)||fs::exists(side))throw std::runtime_error("output already exists; no overwrite");
  const auto plan=readJson(request.at("plan"));
  std::unique_ptr<TFile> input(TFile::Open(inputPath.c_str(),"READ"));
  if(!input||input->IsZombie())throw std::runtime_error("unreadable package");
  auto meta=object(*input,"meta/package");const auto rows=object(*input,"meta/inventory");
  if(meta.contains("replay_provenance"))throw std::runtime_error("replay must start from an original nominal package");
  const auto code=consumerCode();
  if(sha256File(inputPath)!=request.at("parent_package_sha256"))throw std::runtime_error("parent bytes differ from requested input");
  const auto& recipe=plan.at("histograms").at("measurement_replay");
  if(admissionOnly && !recipe.contains("centrality_weight"))
    throw std::runtime_error("read-only replay admission is qualified only for the centrality source");
  if(recipe.at("schema")!="MeasurementReplayRecipeV1" || recipe.at("consumer_code")!=code || recipe.at("admission")!="MECHANICAL_FIXTURE_ONLY")throw std::runtime_error("unqualified or foreign replay recipe");
  if(recipe.at("candidate_archive")!="parent_only")throw std::runtime_error("candidate archive must remain in nominal parent");
  if(!recipe.contains("selection_change") && !recipe.contains("response_only_resolution") && !recipe.contains("centrality_weight") && recipe.at("cone_prescription")!="fixed-measured-cone-mechanical-fixture")throw std::runtime_error("unsupported cone prescription");
  photonjet_selection_replay::validate(recipe,meta.at("selection"),plan.at("selection"),[](const Json& j){return digest(j);});
  if(recipe.contains("selection_change"))for(const auto& source:recipe.at("selection_change").at("definition_sources"))
    if(sha256File(source.at("path"))!=source.at("sha256"))throw std::runtime_error("selection definition source changed");
  Json sample;for(const auto& item:plan.at("samples"))
    if(item.at("system")==meta.at("system")&&item.at("name")==meta.at("name"))sample=item;
  if(sample.is_null())throw std::runtime_error("sample absent from requested plan");
  const bool responseOnly=photonjet_response_only::validate(recipe,sample.at("key"),sample.at("family"));
  const bool centralityOff=photonjet_centrality_replay::validate(recipe,sample.at("key"),sample.at("family"));
  const auto centralityFreezes=photonjet_centrality_replay::admit(recipe,meta.at("selection"),sha256File);
  photonjet_centrality_replay::validatePackage(recipe,meta,centralityFreezes);
  if(recipe.contains("centrality_weight")) {
    const auto& c=recipe.at("centrality_weight");
    Json roster=Json::object();for(const auto& row:plan.at("samples"))roster[row.at("key").get<std::string>()]=row.at("family");
    if(roster.size()!=plan.at("samples").size() || roster!=c.at("sample_families") ||
       request.at("parent_package_sha256")!=c.at("parent_packages").at(sample.at("key").get<std::string>()))
      throw std::runtime_error("centrality replay plan or package differs from complete AuAu parent roster");
  }
  const bool nominalMeasurement=recipe.contains("response_only_resolution")||recipe.contains("centrality_weight");
  if(recipe.contains("response_only_resolution"))for(const auto& source:recipe.at("response_only_resolution").at("definition_sources"))
    if(sha256File(source.at("path"))!=source.at("sha256"))throw std::runtime_error("resolution definition source changed");
  auto originalPlan=plan;originalPlan["histograms"]=meta.at("histograms");originalPlan["selection"]=meta.at("selection");originalPlan["plan_sha256"]=meta.at("plan_sha256");
  auto allowed=meta.at("histograms");allowed["measurement_replay"]=recipe;
  if(allowed!=plan.at("histograms"))throw std::runtime_error("only the explicit record recipe may change");
  if(admissionOnly) {
    fs::path receipt=inputPath;receipt.replace_extension(".json");
    if(!fs::exists(receipt) || !readJson(receipt.string()).contains("validation"))
      throw std::runtime_error("read-only replay admission requires an existing validated parent receipt; no fallback scan or write");
  }
  verifyReplayParent(inputPath,originalPlan,sample);
  const auto count=sample.at("sources").size();
  if(count==0||count>100)throw std::runtime_error("bounded prototype permits 1..100 original sources");
  if(admissionOnly) {
    std::cout<<Json{{"status","PASS_MEASUREMENT_REPLAY_ADMISSION_ONLY"},
      {"admission","MECHANICAL_FIXTURE_ONLY"},{"sample",sample.at("key")},
      {"sources",count},{"consumer_code",code},{"recipe_sha256",digest(recipe)},
      {"event_loops_executed",0},{"package_written",false}}.dump()<<'\n';
    return 0;
  }
  ReductionProgress progress("replay_measurement_package",count+1);
  const auto& h=plan.at("histograms");configureJetRadius(h.at("jet_radius"));
  configuredXjEdges=h.at("native_xj_edges").get<std::vector<double>>();
  configuredMeasurementEtEdges=h.at("native_measurement_et_edges_gev").get<std::vector<double>>();
  fs::create_directories(fs::path(target).parent_path());
  const std::string selectionPath=target+".selection.json";
  if(fs::exists(selectionPath))throw std::runtime_error("selection output already exists");
  std::ofstream(selectionPath)<<plan.at("selection").dump(2)<<'\n';
  const auto system=meta.at("system")=="PP"?System::pp:System::auau;
  auto profile=configuredProfile("nominal",system,selectionPath.c_str());profile.simulation=sample.at("family")!="data";
  configuredDataJetScale=h.at("data_jet_scale").at(system==System::pp?"pp":"auau");
  const auto classes=centralityClasses(system);
  response::ResponseGrids grids;grids.xj=configuredXjEdges;
  grids.reco=h.at("native_reco_response_edges_gev").get<std::vector<double>>();
  grids.truth=h.at("native_truth_response_edges_gev").get<std::vector<double>>();
  const std::string effect=recipe.at("sample_effects").at(sample.at("key").get<std::string>());
  const auto variation=response::parseVariation(effect.c_str());
  if(responseOnly!=response::responseOnly(variation))throw std::runtime_error("response scope and applied effect differ");
  if(centralityOff!=variation.centralityOff)throw std::runtime_error("centrality scope and applied effect differ");
  const bool reevaluatePhotons=energy_selection::changesPhotons(variation)&&!responseOnly;
  const auto& pp=recipe.at("models").at("pp");const auto& npb=recipe.at("models").at("npb");
  energy_selection::Evaluator evaluator(profile,{pp.at("path"),"myBDT",pp.at("sha256")},{npb.at("path"),"myBDT",npb.at("sha256")},
      {{.490,.037},.8},sha256File);
  const auto unchangedCone=[](const energy_selection::Input& in,double){return in.prescribedCone;};
  const Json activeRows=measurementReplayInventory(sample.at("system"),sample.at("family"));
  const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(180);
  entryBegin=0;entryEnd=0;progressPath=target+".progress.json";
  progressTotal=table(*input,"records/candidates/events").GetEntries()*(profile.simulation?2:1);progressOffset=0;
  progressIdentity={{"schema","LongRunningCanaryProgressV1"},{"task","THE-482"},{"workstream","the474_systematics_augmentation_20261014"},
    {"campaign","measurement_package_replay"},{"row_id",sample.at("key")}};
  if(recipe.contains("centrality_weight")) {
    progressIdentity["task"]="THE-485";
    progressIdentity["workstream"]="the482_parallel_reweighting_closure_20261015";
    progressIdentity["campaign"]="integrated_centrality_v1";
  }
  progressStart=progressLast=std::chrono::steady_clock::now();progressPhase="source_index";reportProgress(0);
  std::map<std::string,std::vector<Range>> ranges;
  for(const auto& name:candidateCacheNames())ranges[name]=sourceRanges(table(*input,"records/candidates/"+name),count);
  std::set<std::string> responsePaths;std::map<std::string,std::unique_ptr<TH1>> sums;
  std::map<std::string,std::unique_ptr<THnSparse>> sparseSums;
  for(const auto& row:activeRows)responsePaths.insert(row.at("path"));

  std::unique_ptr<TFile> output(TFile::Open(partial.c_str(),"CREATE"));
  if(!output||output->IsZombie())throw std::runtime_error("cannot open partial package");
  directory(*output,"records")->cd();
  Long64_t source=0,bin=0,begin=0,end=0;double value=0,variance=0;std::string histogram,sourceJson;
  TTree contributions("source_contributions","Source-aware replay cells, including ROOT flow bins");
  contributions.Branch("source",&source);contributions.Branch("histogram",&histogram);
  contributions.Branch("bin",&bin);contributions.Branch("value",&value);contributions.Branch("variance",&variance);
  TTree sources("sources","Original source identities and replay contribution ranges");
  sources.Branch("source",&source);sources.Branch("json",&sourceJson);
  sources.Branch("contribution_begin",&begin);sources.Branch("contribution_end",&end);
  RecordTable savedEvents,savedPhotons,savedJets;Long64_t eventOffset=0;
  auto& oldSources=table(*input,"records/sources");
  std::string* oldJson=nullptr;oldSources.SetBranchAddress("json",&oldJson);
  std::vector<std::string> identities;
  for(Long64_t i=0;i<static_cast<Long64_t>(count);++i) {
    if(oldSources.GetEntry(i)<=0)throw std::runtime_error("source identity read failed");
    const auto origin=Json::parse(*oldJson).at("metadata").at("source");
    identities.push_back(origin.contains("logical_source_id")?origin.at("logical_source_id").get<std::string>():origin.at("input_file_sha256").get<std::string>());
  }
  Json evidence=Json::array();
  for(source=0;source<static_cast<Long64_t>(count);++source) {
    if(std::chrono::steady_clock::now()>deadline)throw std::runtime_error("local replay package runtime ceiling");
    if(oldSources.GetEntry(source)<=0)throw std::runtime_error("source read failed");sourceJson=*oldJson;
    const auto src=Json::parse(sourceJson);const auto& origin=src.at("metadata").at("source");
    const std::string identity=origin.contains("logical_source_id")?origin.at("logical_source_id"):origin.at("input_file_sha256");
    Cells cells;
    candidate_replay::Slice slice; slice.source=source;
    for(const auto& name:candidateCacheNames()) {
      const auto range=ranges.at(name).at(source);
      if(range.second-range.first>10000000)throw std::runtime_error("source table exceeds qualified row ceiling");
      slice.rows[name]=range;
    }
    const auto cached=slice.rows.at("events").second-slice.rows.at("events").first;
    progressRowTotal=cached;
    progressPhase="record_measurement_and_purity";reportProgress(0);
    Histograms measurement(classes);EventTable eventTable;purity_inputs::Book purity(sample.at("family"),classes);
    energy_selection::Context measurementContext{evaluator,unchangedCone,"fixed-measured-cone-mechanical-fixture"};
    const auto summary=candidate_replay::measurement(*input,"records/candidates/",profile,classes,measurement,eventTable,
        1.,count,false,10000000,nominalMeasurement?nullptr:&variation,reevaluatePhotons?&measurementContext:nullptr,identities,&purity,&slice);
    TMemFile responseFile("response.root","RECREATE");responseFile.cd();measurement.write();eventTable.write();purity.write();
    auto* recordDest=directory(*output,"records");
    savedEvents.append(table(responseFile,"table_events"),"events",source,eventOffset,recordDest);
    savedPhotons.append(table(responseFile,"table_photons"),"photons",source,eventOffset,recordDest);
    savedJets.append(table(responseFile,"table_jets"),"jets",source,eventOffset,recordDest);
    eventOffset+=eventTable.run.size();progressOffset+=cached;
    Json responseCounters=nullptr,resolutionCounts=nullptr;
    if(profile.simulation) {
      progressPhase="record_response";reportProgress(0);
      std::vector<std::vector<response::Accumulators>> acc(classes.size());
      for(size_t k=0;k<classes.size();++k)for(double threshold:kJetThresholds)
        acc[k].emplace_back("_"+classes[k].name+"_"+thresholdName(threshold),grids,system==System::pp);
      energy_selection::Context responseContext{evaluator,unchangedCone,"fixed-measured-cone-mechanical-fixture"};
      response::ResolutionMigration migration;
      const Json* admittedWeights=recipe.contains("centrality_weight")?&centralityFreezes.at(sample.at("family").get<std::string>()).at("centrality").at("weights"):nullptr;
      candidate_replay::responseMeasurement(*input,"records/candidates/",profile,classes,grids,acc,identities,1.,variation,10000000,reevaluatePhotons?&responseContext:nullptr,&slice,responseOnly?&migration:nullptr,admittedWeights);
      if(responseOnly)resolutionCounts=response::resolutionMigrationJson(migration);
      responseFile.cd();for(auto& group:acc)for(auto& a:group)a.write();responseCounters=counts(responseContext.counters);progressOffset+=cached;
    }
    for(const auto& row:activeRows) {
      const std::string path=row.at("path"),key=row.at("key");auto* object=responseFile.Get(key.c_str());
      if(auto* sparse=dynamic_cast<THnSparse*>(object)) {
        if(!sparseSums.count(path))sparseSums[path].reset(static_cast<THnSparse*>(sparse->Clone()));
        else sparseSums[path]->Add(sparse);
        // Native source_contributions stores dense cells only. Preserve that
        // format; sparse diagnostics still retain their full aggregate Sumw2.
        continue;
      }
      auto* response=dynamic_cast<TH1*>(object);
      if(!response)throw std::runtime_error("missing or unsupported replay histogram "+key);
      if(!sums.count(path)){sums[path].reset(static_cast<TH1*>(response->Clone()));sums[path]->SetDirectory(nullptr);}
      else { sums[path]->Add(response);photonjet_science::mergeDefinitions(*sums[path],*response); }
      for(bin=0;bin<response->GetNcells();++bin) {
        value=response->GetBinContent(bin);variance=response->GetSumw2()->At(bin);
        if(value!=0||variance!=0)cells[{path,bin}]={value,variance};
      }
    }
    begin=contributions.GetEntries();
    for(const auto& [key,cell]:cells){histogram=key.first;bin=key.second;value=cell.first;variance=cell.second;contributions.Fill();}
    end=contributions.GetEntries();sources.Fill();
    evidence.push_back({{"source",source},{"logical_source_id",identity},{"candidate_events",summary.cacheEvents},{"contribution_cells",end-begin},{"measurement_counts",counts(measurementContext.counters)},{"response_counts",responseCounters},{"response_only_resolution_counts",resolutionCounts}});
    progress.advance(source+1);
  }
  oldSources.ResetBranchAddresses();
  for(const auto& [path,hist]:sums){const auto slash=path.rfind('/');directory(*output,path.substr(0,slash))->cd();hist->Write(path.substr(slash+1).c_str(),TObject::kOverwrite);}
  for(const auto& [path,hist]:sparseSums){const auto slash=path.rfind('/');directory(*output,path.substr(0,slash))->cd();hist->Write(path.substr(slash+1).c_str(),TObject::kOverwrite);}
  output->GetDirectory("records")->cd();contributions.Write("",TObject::kOverwrite);sources.Write("",TObject::kOverwrite);
  contributions.SetDirectory(nullptr);sources.SetDirectory(nullptr);
  savedEvents.write();savedPhotons.write();savedJets.write();
  // Candidate archives remain immutable in the original nominal parent.
  meta["parent_candidate_records"]=meta.at("candidate_records");
  meta.erase("candidate_records");
  if(recipe.contains("selection_change")) {
    meta["parent_selection"]=meta.at("selection");meta["selection"]=plan.at("selection");
    meta["render_selection"]["selection_sha256"]=plan.at("selection").at("selection_sha256");
    meta["render_selection"]["photon_id"]=plan.at("selection").at("systems").at(system==System::pp?"pp":"auau");
  }
  meta["histograms"]=h;meta["render_selection"]["histogram_settings"]=h;
  meta["plan_sha256"]=plan.at("plan_sha256");meta["package_binding_sha256"]=packageBinding(plan,sample);
  meta["histograms_written"]=activeRows.size();meta["event_records"]=eventOffset;
  meta["replay_provenance"]={{"schema","PhotonJetMeasurementReplayV1"},{"parent_package_sha256",sha256File(inputPath)},
    {"consumer_code",code},{"original_producer_code",meta.at("render_selection").at("producer_code")},
    {"recipe",recipe},{"recipe_sha256",digest(recipe)},{"applied_effect",effect},{"consumer_root_version",gROOT->GetVersion()},
    {"normalization_applied",false},{"candidate_archive","parent_only"},{"scope","measurement_response_purity_inputs"},{"nominal_diagnostics","parent_only"},{"sources",evidence}};
  directory(*output,"meta")->cd();TObjString(meta.dump().c_str()).Write("package",TObject::kOverwrite);
  TObjString(activeRows.dump().c_str()).Write("inventory",TObject::kOverwrite);
  if(auto* normalization=input->Get("meta/pb_per_generated_event"))normalization->Write("pb_per_generated_event");
  output->Close();output.reset();progressPhase="validation_and_publication";reportProgress(0);
  publishMeasurementReplayPackage(partial.string(),target,plan,sample,inputPath,code);
  progress.finish();progressPhase="complete";reportProgress(0);
  std::cout<<Json{{"status","PASS_MEASUREMENT_REPLAY_PACKAGE"},{"package",target},{"sources",evidence}}.dump()<<'\n';return 0;
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::runtime_error("one request JSON required");return replay_package::run(treetohists::readJson(argv[1]));}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
