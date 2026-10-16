/* Bounded native response-component fanout. The owner resolves its normal
   selection/weight context first. This API neither plans nor dispatches jobs,
   fills nominal QA/spectra, nor claims complete varied sample packages. */
#pragma once
#include "TreeToHists.h"
#include "CollectionContract.h"
#include "FillPhotonJetResponse.C"
#include <TSystem.h>
#include <filesystem>
#include <set>

namespace photonjet { namespace response { namespace native_grid {
using Json=nlohmann::json;
namespace fs=std::filesystem;
inline void require(bool ok,const char *why) {if(!ok)throw std::runtime_error(std::string("native response grid: ")+why);}
inline bool pin(const Json &j) {
  if(!j.is_string())return false;const auto s=j.get<std::string>();
  return s.size()==64 && s.find_first_not_of("0123456789abcdef")==std::string::npos;
}
inline std::uint64_t number(const Json &j) {
  require(j.is_number_integer() && j.get<long long>()>=0,"expected a nonnegative integer");return j.get<std::uint64_t>();
}
inline Json fit(const TH1 *h) {
  if(!h)return nullptr;
  Json edges=Json::array(),bins=Json::array();
  for(int i=1;i<=h->GetNbinsX()+1;++i)edges.push_back(h->GetXaxis()->GetBinLowEdge(i));
  for(int i=0;i<h->GetNcells();++i)bins.push_back(h->GetBinContent(i));
  return {{"edges",edges},{"contents",bins}};
}
inline Json profileState(const Profile &p) {
  const auto line=[](const Line &l){return Json::array({l.intercept,l.slope});};Json cells=Json::array();
  for(const auto &c:p.idCells)cells.push_back({{"range",{c.etMin,c.etMax,c.centMin,c.centMax}},
    {"tight",line(c.tight)},{"nontight_high",line(c.nonTightHigh)}});
  return {{"photon",{p.photonEtMin,p.photonEtMax,p.photonAbsEtaMax}},{"vertex_abs_z",p.vertexAbsZMax},
    {"tight",line(p.tight)},{"nontight_low",line(p.nonTightLow)},{"nontight_high",line(p.nonTightHigh)},
    {"nontight_low_constant",{p.nonTightLowIsConstant,p.nonTightLowConstant}},
    {"isolation_indices",{p.isolationRadiusIndex,p.isolationSubtraction}},
    {"jet",{p.jetRadius,p.jetPtMin,p.jetAbsEtaMax,p.overlapDeltaR,p.deltaPhiMin}},
    {"augmented_classes",p.augmentedClasses},{"sideband",{p.nonTightMin,p.nonTightLowShift,p.sidebandGap}},
    {"flags",{p.configured,p.isolationFromWitness,p.verifyNominalID,p.verifyNominalIsolation,p.simulation,p.ppIsolationParity,p.optionalID,p.nativeID}},
    {"isolation_threshold",line(p.isolationThreshold)},{"isolation_centrality",{p.isolationCentMin,p.isolationCentMax}},
    {"score_branch",p.scoreBranch},{"state_branch",p.stateBranch},{"selection_id",p.selectionId},{"id_cells",cells}};
}
inline void encodeNonfinite(Json &value) {
  if(value.is_number_float() && !std::isfinite(value.get<double>())) {
    const double number=value.get<double>();value=std::isnan(number) ? "NaN" : (number>0 ? "+Inf" : "-Inf");
  } else if(value.is_object() || value.is_array())for(auto &child:value)encodeNonfinite(child);
}
/* Capture effective callback state, including actual fits, rather than trusting
   a caller-supplied context label. No numerical parameter is supplied here. */
inline Json context(const Profile &p,const ResponseGrids &g,const CentralityWeights &w) {
  Json periods=Json::array();
  for(const auto &period:ppPeriods)periods.push_back({{"name",period.name},{"share",period.share},
    {"double_fraction",period.doubleFraction},{"fit",fit(period.fit.get())}});
  const char *owner=configuredOwnership.owner==Ownership::Owner::none ? "none" :
    configuredOwnership.owner==Ownership::Owner::leadingTruthJet ? "leading_truth_jet" : "filter_photon";
  Json value={{"profile",p.name},{"profile_state",profileState(p)},{"system",p.system==System::pp ? "PP" : "AuAu"},{"selection",p.selectionJson},
    {"jet_radius",p.jetRadius},{"configured_jet_radius",configuredJetRadius},{"strict_reduction",strictReduction},
    {"xj_edges",configuredXjEdges},{"measurement_et_edges",configuredMeasurementEtEdges},
    {"ownership",{{"owner",owner},{"low",configuredOwnership.low},
      {"high",std::isfinite(configuredOwnership.high) ? Json(configuredOwnership.high) : Json()}}},
    {"auau_vertex_fit",fit(auauVertexFit.get())},{"auau_receipt_sha256",auauReceiptSha},
    {"auau_stitching_sha256",auauStitchingSha},{"pp_periods",periods},
    {"grids",{{"reco",g.reco},{"truth",g.truth},{"xj",g.xj}}},{"centrality_weights",w.perSlice}};
  encodeNonfinite(value);return value; // preserve distinct sanctioned sentinels in the context hash
}
inline void atomic(const fs::path &path,const Json &value) {
  const auto next=path.string()+".next";std::ofstream stream(next);stream<<value.dump(2)<<'\n';stream.close();
  require(bool(stream),"cannot write atomic progress/receipt");fs::rename(next,path);
}
inline std::uint64_t residentBytes() {
  ProcInfo_t info;require(gSystem->GetProcInfo(&info)==0 && info.fMemResident>0,"resident-memory observation unavailable");
  return std::uint64_t(info.fMemResident)*1024;
}
inline std::uint64_t denseCells(const ResponseGrids &g,bool pp) {
  const std::uint64_t r=g.reco.size()-1,t=g.truth.size()-1,x=g.xj.size()-1;
  const auto a=associationEtEdges().size()+1,b=associationXjEdges().size()+1,d=deltaPhiEdges().size()+1;
  return 2*(r+2)+2*(t+2)+(r+2)*(t+2)+8*(r+2)*(x+2)+2*(t+2)*(x+2)+
    (r*x+2)*(t*x+2)+(2+(pp?3:0))*a*b+6*(r+2)+4*(r+2)*(x+2)*d;
}
struct Restore {
  Long64_t begin=entryBegin,end=entryEnd;std::string progress=progressPath;
  ~Restore(){entryBegin=begin;entryEnd=end;progressPath=progress;}
};
inline Json resolvedPlan(const Json &spec,const Json &variant) {
  Json plan={{"schema","PhotonJetNativeResponsePlanV1"},{"producer_code",spec.at("producer_code")},
    {"runtime_context",spec.at("runtime_context")},{"sources",spec.at("sources")},
    {"name",variant.at("name")},{"variation",variant.at("variation")}};
  if(spec.contains("ems_selection_contract"))plan["ems_selection_contract"]=spec.at("ems_selection_contract");
  return plan;
}

// EMS is a calibration shift: evaluate the existing whole-selection contract
// at scaled original ET. Named EMR/JER still use response-only classification.
// Keep the optional field absent on old grids so their identity is unchanged.
inline bool settingsKeys(Json config) {
  config.erase("ems_selection_contract");
  return photonjet_collection::keys(config,{"schema","variants","limits"});
}
inline void copyEmsContract(const Json &config,Json &spec) {
  if(config.contains("ems_selection_contract"))spec["ems_selection_contract"]=config.at("ems_selection_contract");
}
inline bool supportedRecipe(const Variation &v) {
  if(responseOnly(v))return true;
  if(emsSelectionResponse(v))return !v.centralityOff && v.photonResolution=="ppg18_nominal" &&
    v.jetResolution=="sphenix_pp2026_nominal" && v.emScale!=1.;
  return v.centralityOff && v.photonResolution=="off" && v.jetResolution=="off" &&
    v.emScale==1. && v.photonScale==1. && v.jetScale==1. && v.photonSmear==0. && v.jetSmear==0.;
}
inline void validateEmsContract(const Json &config) {
  bool ems=false;
  for(const auto &v:config.at("variants"))ems=ems || emsSelectionResponse(parseVariation(v.at("variation").get<std::string>().c_str()));
  require(ems==config.contains("ems_selection_contract"),"EMS components require exactly one explicit selection contract; other grids must omit it");
  if(ems)ems_native::validateContract(config.at("ems_selection_contract"));
}

/* An opt-in setting on the existing campaign planner. It declares response
   components only; the nominal measurement/QA and DATA paths keep their owners. */
inline void validateSettings(const Json &config) {
  using photonjet_collection::keys;
  require(settingsKeys(config) && config.at("schema")=="PhotonJetNativeResponseGridSettingsV1","invalid settings schema/keys");
  const auto &variants=config.at("variants"),&limits=config.at("limits");
  require(variants.is_array() && !variants.empty() && variants.size()<=16,"settings require 1..16 variants");
  require(keys(limits,{"max_histogram_bytes","max_resident_bytes","max_output_bytes","max_sparse_bins_per_histogram","max_seconds"}),"settings require explicit closed resource bounds");
  for(const auto &[name,value]:limits.items())require(number(value)>0,"settings resource bounds must be positive integers");
  const auto hist=number(limits.at("max_histogram_bytes")),resident=number(limits.at("max_resident_bytes")),output=number(limits.at("max_output_bytes")),sparse=number(limits.at("max_sparse_bins_per_histogram")),seconds=number(limits.at("max_seconds"));
  require(hist<=64ULL*1024*1024*1024 && resident>=hist && resident<=64ULL*1024*1024*1024 && output<=64ULL*1024*1024*1024 && sparse<=1000000 && seconds<=43200,"settings resource bounds exceed native envelope");
  std::set<std::string> names,recipes;bool nominal=false;
  for(const auto &variant:variants) {
    require(keys(variant,{"name","variation"}) && variant.at("name").is_string() && variant.at("variation").is_string(),"invalid settings variant");
    const auto name=variant.at("name").get<std::string>(),recipe=variant.at("variation").get<std::string>();
    require(!name.empty() && name.size()<=64 && name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")==std::string::npos && names.insert(name).second,"duplicate/unsafe settings name");
    const auto parsed=parseVariation(recipe.c_str());
    require(recipe.empty() || supportedRecipe(parsed),"unsupported response component recipe");
    require(recipes.insert(variationJson(parsed)).second,"duplicate resolved settings recipe");
    if(recipe.empty()){require(name=="nominal","empty recipe must be named nominal");nominal=true;}
    else require(name!="nominal","nominal label must have the empty recipe");
  }
  require(nominal,"settings require the nominal response component");
  validateEmsContract(config);
}

/* `spec` is a create-once, closed, fully resolved finite grid. Each source is
   checksum-verified once, then opened once by the native reader. Original
   source/entry/object indices feed the existing shared random draws. */
inline Json run(const Json &spec,const Profile &profile,const ResponseGrids &grids,const CentralityWeights &weights,
    const std::vector<ResponseView>& borrowedSelectionViews={}) {
  require(borrowedSelectionViews.size()<=16,"borrowed selection sink roster exceeds bound");
  for(const auto& view:borrowedSelectionViews)require(view.acc && view.measurement && view.measurementRecords && selection_views::cutComposition(view.variation),"borrowed sink must be a complete nominal-recipe measurement/response view");
  using photonjet_collection::keys;
  auto closed=spec;closed.erase("ems_selection_contract");
  require(keys(closed,{"schema","identity","sources","variants","producer_code","runtime_context","limits","output_directory"}),"unknown/missing grid key");
  require(spec.at("schema")=="PhotonJetNativeResponseGridV1","unknown grid schema");
  require(keys(spec.at("identity"),{"task","workstream","campaign","row"}),"complete identity is required");
  for(const auto &[key,value]:spec.at("identity").items())require(value.is_string() && !value.get<std::string>().empty() && value.get<std::string>().size()<=256,"empty/oversized grid identity");
  const auto code=treetohists::reductionCode();treetohists::checkReductionCode({{"producer_code",code}});
  require(spec.at("producer_code")==code,"grid code differs from actual source/binary");
  const auto effective=context(profile,grids,weights);
  require(spec.at("runtime_context")==effective,"resolved selector/weight context differs");
  require((profile.system==System::pp && weights.perSlice.empty()) ||
    (profile.system==System::auau && weights.perSlice.size()==16),"resolved centrality map has the wrong system/dimensions");
  for(double weight:weights.perSlice)require(std::isfinite(weight) && weight>0.,"resolved centrality map is not finite/positive");
  const auto &sources=spec.at("sources"),&variants=spec.at("variants"),&limits=spec.at("limits");
  require(sources.is_array() && !sources.empty() && sources.size()<=64,"source roster must be finite (1..64)");
  require(variants.is_array() && !variants.empty() && variants.size()<=16,"variant roster must be finite (1..16)");
  validateEmsContract(spec);
  if(spec.contains("ems_selection_contract"))require(spec.at("ems_selection_contract").at("selection")==Json::parse(profile.selectionJson),"EMS contract selection differs from resolved grid profile");
  require(keys(limits,{"max_histogram_bytes","max_resident_bytes","max_output_bytes","max_sparse_bins_per_histogram","max_seconds"}),"unknown/missing resource bound");
  const auto histogramLimit=number(limits.at("max_histogram_bytes")),residentLimit=number(limits.at("max_resident_bytes")),
    outputLimit=number(limits.at("max_output_bytes")),sparseLimit=number(limits.at("max_sparse_bins_per_histogram")),secondsLimit=number(limits.at("max_seconds"));
  require(histogramLimit>0 && histogramLimit<=64ULL*1024*1024*1024 && residentLimit>=histogramLimit &&
          residentLimit<=64ULL*1024*1024*1024 && outputLimit>0 && outputLimit<=64ULL*1024*1024*1024 &&
          sparseLimit>0 && sparseLimit<=1000000 && secondsLimit>0 && secondsLimit<=43200,"invalid finite resource bound");
  require(grids.reco.size()>=2 && grids.truth.size()>=2 && grids.xj.size()>=2 &&
    grids.reco.size()<=128 && grids.truth.size()<=128 && grids.xj.size()<=128,"grid axes exceed bounded native dimensions");
  for(const auto *axis:{&grids.reco,&grids.truth,&grids.xj})for(size_t i=0;i<axis->size();++i)
    require(std::isfinite((*axis)[i]) && (!i || (*axis)[i]>(*axis)[i-1]),"invalid native axis");
  const auto classes=centralityClasses(profile.system);
  const std::uint64_t accumulatorCount=(variants.size()+borrowedSelectionViews.size())*classes.size()*kJetThresholds.size();
  const std::uint64_t sparseReserved=((sparseLimit+16383)/16384)*16384*128;
  const std::uint64_t estimated=accumulatorCount*(denseCells(grids,profile.system==System::pp)*16+1048576+sparseReserved);
  require(estimated<=histogramLimit && estimated<=outputLimit,"declared histogram/output budget cannot hold this grid");
  require(residentBytes()+estimated+64*1024*1024<=residentLimit,"resident budget lacks histogram/runtime reserve");
  std::set<std::string> names,recipes,paths,originals;std::vector<Variation> parsed;
  for(const auto &variant:variants) {
    require(keys(variant,{"name","variation","plan","plan_sha256"}) && variant.at("name").is_string() &&
      variant.at("variation").is_string() && pin(variant.at("plan_sha256")),"invalid variant declaration");
    require(variant.at("plan")==resolvedPlan(spec,variant) && treetohists::digest(variant.at("plan"))==variant.at("plan_sha256").get<std::string>(),"variant plan/source/code/recipe binding differs");
    const auto name=variant.at("name").get<std::string>(),recipe=variant.at("variation").get<std::string>();
    require(!name.empty() && name.size()<=64 && name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")==std::string::npos &&
      names.insert(name).second,"duplicate/unsafe grid name");
    auto variation=parseVariation(recipe.c_str());
    require(recipe.empty() || supportedRecipe(variation),"unsupported grid response component recipe");
    require(recipes.insert(variationJson(variation)).second,"duplicate resolved variation recipe");
    validateVariationSystem(variation,profile.system);parsed.push_back(variation);
  }
  std::uint64_t total=0;
  for(const auto &source:sources) {
    require(keys(source,{"path","sha256","input_file_sha256","entries","entry_begin","entry_end"}) &&
      source.at("path").is_string() && pin(source.at("sha256")) && source.at("input_file_sha256").is_string(),"invalid source declaration");
    const auto path=source.at("path").get<std::string>(),original=source.at("input_file_sha256").get<std::string>();
    require(fs::path(path).is_absolute() && fs::path(path).lexically_normal()==fs::path(path) && paths.insert(path).second &&
      !original.empty() && originals.insert(original).second,"source identity duplicated or not absolute");
    const auto begin=number(source.at("entry_begin")),end=number(source.at("entry_end")),entries=number(source.at("entries"));
    require(begin<=end && end<=entries && (begin<end || entries==0),"invalid finite source range");
    total+=end-begin;require(total<=100000,"grid exceeds 100000 original entries");
  }
  const fs::path output=spec.at("output_directory").get<std::string>();
  require(output.is_absolute() && output.lexically_normal()==output && !fs::exists(output),"output must be a new absolute directory");
  require(fs::create_directory(output),"cannot create output directory");
  const auto gridSha=treetohists::digest(spec);atomic(output/"grid.json",spec);
  const auto start=std::chrono::steady_clock::now();auto last=start;
  std::uint64_t processed=0,rowProcessed=0,rowTotal=0,opens=0,checksums=0,peak=0;Json sourceRecords=Json::array();
  std::string phase="booking",error;bool terminal=false;
  auto progress=[&](bool force=false) {
    const auto now=std::chrono::steady_clock::now();const double seconds=std::chrono::duration<double>(now-start).count();
    // This is a planning estimate, not a production kill timer. Bounded
    // foreground/canary callers own their separate process deadline.
    const auto memory=residentBytes();peak=std::max(peak,memory);require(memory<=residentLimit || phase=="FIRST_BAD" || terminal,"resident-memory bound exceeded");
    if(!force && processed%1000 && std::chrono::duration<double>(now-last).count()<2)return;
    last=now;const double rate=seconds>0 ? processed/seconds : 0.;Json record=spec.at("identity");
    record.update({{"schema","PhotonJetNativeResponseGridProgressV1"},{"grid_sha256",gridSha},{"phase",phase},
      {"processed",processed},{"total",total},{"remaining",total-processed},{"row_processed",rowProcessed},{"row_total",rowTotal},
      {"row_remaining",rowTotal-rowProcessed},{"event_variation_processed",processed*(variants.size()+borrowedSelectionViews.size())},
      {"event_variation_total",total*(variants.size()+borrowedSelectionViews.size())},{"event_variation_remaining",(total-processed)*(variants.size()+borrowedSelectionViews.size())},
      {"deterministic_skip_processed",0},{"deterministic_skip_total",0},{"deterministic_skip_remaining",0},
      {"rate",rate},{"eta_seconds",rate>0 ? (total-processed)/rate : -1.},{"elapsed_seconds",seconds},
      {"planned_wall_seconds",secondsLimit},{"over_plan",seconds>secondsLimit},
      {"updated_at_unix",std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},
      {"latest_error",error.empty() ? Json() : Json(error)},{"source_file_opens",opens},{"input_checksum_passes",checksums},
      {"estimated_histogram_bound_bytes",estimated},{"resident_bytes",memory},{"peak_observed_resident_bytes",peak},
      {"work_definition","K declared variation decisions per original entry; shared skips complete all K slots"}});
    atomic(output/"progress.json",record);if(progressObserver)progressObserver(record);
  };
  Restore restore;progressPath.clear();
  try {
    progress(true);std::vector<ResponseAccumulators> accumulators(variants.size());std::vector<ResponseView> views;
    for(size_t v=0;v<variants.size();++v) {
      auto &acc=accumulators[v];acc.resize(classes.size());
      for(size_t c=0;c<classes.size();++c)for(double threshold:kJetThresholds)
        acc[c].emplace_back("_"+classes[c].name+"_"+thresholdName(threshold),grids,profile.system==System::pp);
      views.push_back({parsed[v],&acc});
    }
    views.insert(views.end(),borrowedSelectionViews.begin(),borrowedSelectionViews.end());
    auto checkSparse=[&](){for(const auto &view:views)for(const auto &acc:std::vector<const ResponseAccumulators*>{view.acc})for(const auto &perClass:*acc)for(const auto &h:perClass)
      require(std::uint64_t(h.pairReason->GetNbins())<=sparseLimit,"sparse-bin bound exceeded");};
    for(const auto &source:sources) {
      phase="reading";rowProcessed=0;rowTotal=number(source.at("entry_end"))-number(source.at("entry_begin"));
      entryBegin=number(source.at("entry_begin"));entryEnd=number(source.at("entry_end"));progress(true);
      const auto path=source.at("path").get<std::string>();
      require(treetohists::sha256File(path)==source.at("sha256").get<std::string>(),"input file hash differs");++checksums;
      const auto size=fs::file_size(path);const auto mtime=fs::last_write_time(path);
      ResponseReadObserver observer;
      observer.opened=[&](TFile &file,TTreeReader &reader){
        ++opens;require(reader.GetEntries()==source.at("entries").get<Long64_t>(),"input entry total differs");
        const auto original=sourceSha(file);require(original==source.at("input_file_sha256").get<std::string>(),"original source identity differs");
        sourceRecords.push_back({{"path",path},{"input_file_sha256",original}});
      };
      observer.completed=[&](long long){++processed;++rowProcessed;require(processed<=total && rowProcessed<=rowTotal,"progress exceeds declared work");checkSparse();progress();};
      const auto emsText=spec.contains("ems_selection_contract")?spec.at("ems_selection_contract").dump():std::string();
      fillFileViews(path,profile,classes,grids,weights,views,emsText.c_str(),treetohists::sha256File,observer);
      require(rowProcessed==rowTotal,"reader ended before declared entry range");
      require(fs::file_size(path)==size && fs::last_write_time(path)==mtime,"source changed during native read");progress(true);
    }
    require(processed==total && opens==sources.size(),"incomplete original-source pass");phase="writing";progress(true);
    const std::uint64_t controlReserve=65536;
    std::uint64_t bytes=fs::file_size(output/"grid.json")+controlReserve;Json products=Json::array();
    for(size_t v=0;v<variants.size();++v) {
      const auto name=variants[v].at("name").get<std::string>();const auto partial=output/(name+".partial.root");
      Json provenance={{"schema","PhotonJetNativeResponseGridProductV1"},{"scope","RESPONSE_COMPONENT_ONLY"},
        {"grid_sha256",gridSha},{"grid",spec},{"variant",variants[v]},{"variant_index",v}};
      const auto emsText=spec.contains("ems_selection_contract")?spec.at("ems_selection_contract").dump():std::string();
      writeProduct(partial.c_str(),profile,grids,weights,parsed[v],accumulators[v],sourceRecords.dump(),emsText.c_str(),provenance.dump());
      bytes+=fs::file_size(partial);require(bytes<=outputLimit,"output-byte bound exceeded");
      products.push_back({{"name",name},{"path",name+".root"},{"bytes",fs::file_size(partial)},
        {"sha256",treetohists::sha256File(partial.string())},{"plan_sha256",variants[v].at("plan_sha256")},
        {"variation",Json::parse(variationJson(parsed[v]))}});progress(true);
    }
    phase="publishing";progress(true);
    Json receipt={{"schema","PhotonJetNativeResponseGridReceiptV1"},{"status","PASS_RESPONSE_COMPONENT_ONLY"},
      {"grid_sha256",gridSha},{"original_entries",processed},{"event_variation_work",processed*(variants.size()+borrowedSelectionViews.size())},
      {"source_file_opens",opens},{"input_checksum_passes",checksums},{"borrowed_selection_sink_count",borrowedSelectionViews.size()},{"products",products},{"output_bytes_including_control_reserve",bytes},
      {"control_reserve_bytes",controlReserve},
      {"peak_observed_resident_bytes",peak},{"ems_measurement_qa_fanout","UNQUALIFIED"},{"science_acceptance",false}};
    require(receipt.dump(2).size()+1+fs::file_size(output/"progress.json")<=controlReserve,"control-output reserve exceeded");
    for(const auto &product:products)fs::rename(output/(product.at("name").get<std::string>()+".partial.root"),output/product.at("path").get<std::string>());
    phase="complete";terminal=true;progress(true);atomic(output/"receipt.json",receipt);return receipt;
  } catch(const std::exception &e) {
    error=e.what();phase="FIRST_BAD";try{progress(true);}catch(...){}
    throw;
  }
}
} } } // photonjet::response::native_grid
