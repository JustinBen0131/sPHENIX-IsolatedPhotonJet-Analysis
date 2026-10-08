/* Source-aware histogram merging and package validation. Histograms remain
   unscaled; FinalAnalysis applies each sample's normalization once. */
#include "TreeToHists.h"
#include "Configuration.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <TFile.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TROOT.h>

namespace treetohists {
namespace fs = std::filesystem;
thread_local ReductionProgress *ReductionProgress::active_=nullptr;
ReductionProgress::ReductionProgress(const std::string &operation,long long total) {
  const char *path=std::getenv("TREETOHISTS_PROGRESS_PATH");
  if (!path || !*path || active_) return;
  if (total<=0 || !fs::path(path).is_absolute()) throw std::runtime_error("invalid reduction progress binding");
  path_=path; operation_=phase_=operation; total_=total; exceptions_=std::uncaught_exceptions();
  start_=last_=std::chrono::steady_clock::now();
  emit(true); active_=this;
}
ReductionProgress::~ReductionProgress() noexcept {
  if (path_.empty()) return;
  if (!finished_) {
    phase_=std::uncaught_exceptions()>exceptions_ ? "failed" : "incomplete";
    try {emit(true);} catch (const std::exception &e) {std::cerr<<"[TreeToHists] progress: "<<e.what()<<'\n';}
  }
  active_=nullptr;
}
void ReductionProgress::emit(bool force) {
  if (path_.empty()) return;
  const auto now=std::chrono::steady_clock::now();
  if (!force && processed_-emitted_<10000 && std::chrono::duration<double>(now-last_).count()<2) return;
  const double elapsed=std::chrono::duration<double>(now-start_).count();
  const double rate=elapsed>0 ? processed_/elapsed : 0;
  Json row={{"schema","LongRunningCanaryProgressV1"},{"producer","TreeToHists"},
    {"row_id",operation_},{"unit","work units"},{"phase",phase_},
    {"processed",processed_},{"total",total_},{"remaining",total_-processed_},
    {"row_processed",processed_},{"row_total",total_},{"row_remaining",total_-processed_},
    {"deterministic_skip_processed",0},{"deterministic_skip_total",0},{"deterministic_skip_remaining",0},
    {"elapsed_seconds",elapsed},{"rate",rate},{"eta_seconds",rate>0 ? (total_-processed_)/rate : -1},
    {"updated_at_unix",std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},
    {"latest_error",phase_=="failed" || phase_=="incomplete" ? Json("native reduction did not complete") : Json()}};
  const fs::path next=path_+".next";
  std::ofstream out(next);out<<row.dump()<<'\n';out.close();
  if (!out) throw std::runtime_error("cannot write native reduction progress");
  fs::rename(next,path_);last_=now;emitted_=processed_;
}
void ReductionProgress::advance(long long processed) {
  if (path_.empty()) return;
  if (processed<processed_ || processed>=total_) throw std::runtime_error("native reduction progress regressed or completed early");
  processed_=processed;phase_=operation_;emit(processed==1);
}
void ReductionProgress::finish() {
  if (path_.empty()) return;
  processed_=total_;phase_="complete";emit(true);finished_=true;
}
void ReductionProgress::pulse(const std::string &phase) {
  if (!active_) return;
  active_->phase_=phase;active_->emit(false);
}

namespace {
Json object(TFile &file, const char *name) {
  auto *o = dynamic_cast<TObjString *>(file.Get(name));
  if (!o) throw std::runtime_error(std::string(file.GetName()) + ": missing " + name);
  return Json::parse(o->GetString().Data());
}
TTree &tree(TFile &file, const char *name) {
  auto *t = dynamic_cast<TTree *>(file.Get(name));
  if (!t) throw std::runtime_error(std::string(file.GetName()) + ": missing " + name);
  return *t;
}
TDirectory *directory(TFile &file, const std::string &path) {
  if (!file.GetDirectory(path.c_str())) file.mkdir(path.c_str(), "", true);
  auto *d = file.GetDirectory(path.c_str());
  if (!d) throw std::runtime_error("cannot create " + path);
  return d;
}
void writeJson(const fs::path &path, const Json &value) {
  const fs::path next = path.string() + ".partial";
  std::ofstream out(next); out << value.dump(1) << '\n'; out.close();
  if (!out) throw std::runtime_error("cannot write " + path.string());
  fs::rename(next, path);
}
Json sidecar(TFile &file, const std::string &path) {
  const Json meta = object(file, "meta/package");
  Json out={{"package",fs::path(path).filename().string()}, {"sha256",sha256File(path)},
    {"sources",meta.at("source_count")}, {"histograms",meta.at("histograms_written")},
    {"event_records",meta.at("event_records")}, {"availability",meta.at("availability")},
    {"plan_sha256",meta.at("plan_sha256")}, {"selection_sha256",meta.at("selection").at("selection_sha256")}};
  if (meta.contains("package_binding_sha256")) out["package_binding_sha256"]=meta.at("package_binding_sha256");
  out["validation"]={{"schema","PhotonJetPackageValidationV1"},{"status","PASS"},
    {"producer_code_sha256",digest(meta.at("render_selection").at("producer_code"))},
    {"root_version",meta.at("runtime").at("root_version")},
    {"checks",Json::array({"identity_and_settings","histograms_and_sumw2","source_ranges_and_contributions","event_references"})}};
  return out;
}
void sameAxis(const TAxis &a, const TAxis &b) {
  if (a.GetNbins() != b.GetNbins()) throw std::runtime_error("histogram axis differs");
  for (int i=1; i<=a.GetNbins()+1; ++i)
    if (a.GetBinLowEdge(i) != b.GetBinLowEdge(i)) throw std::runtime_error("histogram bin edges differ");
}
void finiteHistogram(TObject &o, const std::string &path) {
  if (auto *h = dynamic_cast<TH1 *>(&o)) {
    if (h->GetSumw2N() != h->GetNcells()) throw std::runtime_error("missing Sumw2: " + path);
    for (int i=0; i<h->GetNcells(); ++i) {
      if (i%10000==0) ReductionProgress::pulse("validating_histograms");
      if (!std::isfinite(h->GetBinContent(i)) || !std::isfinite(h->GetSumw2()->At(i)) || h->GetSumw2()->At(i)<0)
        throw std::runtime_error("invalid histogram cell: " + path);
    }
  } else if (auto *h = dynamic_cast<THnSparse *>(&o)) {
    if (!h->GetCalculateErrors()) throw std::runtime_error("missing sparse Sumw2: " + path);
    for (Long64_t i=0; i<h->GetNbins(); ++i) {
      if (i%10000==0) ReductionProgress::pulse("validating_histograms");
      if (!std::isfinite(h->GetBinContent(i)) || !std::isfinite(h->GetBinError2(i)) || h->GetBinError2(i)<0)
        throw std::runtime_error("invalid sparse cell: " + path);
    }
  } else throw std::runtime_error("unsupported histogram object: " + path);
}
}  // namespace

/* A collection grows; an unchanged sample does not. Bind reusable packages
   to their own original inputs and science settings, independently of the
   collection's location, job partitioning and other samples. */
namespace {
size_t sourceCount(const Json &sample) {
  return sample.contains("_source_count") ? sample.at("_source_count").get<size_t>() : sample.at("sources").size();
}
std::string bindingWithSources(const Json &plan,const Json &sample,
                              const std::function<void(std::ostream &)> &writeSources) {
  Json body={{"schema","PhotonJetPackageBindingV1"},{"system",sample.at("system")},{"name",sample.at("name")},
    {"family",sample.at("family")},{"sources",nullptr},{"selection",plan.at("selection")},
    {"histograms",plan.at("histograms")},{"producer_code",plan.at("producer_code")},{"root_version",plan.at("root_version")},
    {"normalization",sample.at("normalization")},{"coverage",sample.at("coverage")},
    {"weights",sample.at("weights")},{"vertex_reweighting",sample.at("vertex_reweighting")},
    {"intermediate",sample.value("intermediate",false)}};
  for (const char *key : {"centrality_evidence","centrality_data_trigger","run","data_run","run_packages","run_coverages","run_records"})
    body[key]=sample.value(key,Json());
  /* Emit the same sorted canonical JSON as the original binding, one source
     at a time. The digest and reuse contract are independent of merge size. */
  return digestStream([&](std::ostream &out) {
    out<<'{'; bool first=true;
    for (const auto &[key,value] : body.items()) {
      if (!first) out<<','; first=false;
      out<<Json(key).dump(-1,' ',true)<<':';
      if (key!="sources") {out<<value.dump(-1,' ',true); continue;}
      writeSources(out);
    }
    out<<'}';
  });
}
}  // namespace

std::string packageBindingFromSources(const Json &plan,const Json &sample,
                                     const std::function<void(std::ostream &)> &writeSources) {
  return bindingWithSources(plan,sample,writeSources);
}

std::string packageBinding(const Json &plan,const Json &sample) {
  /* Compact merge descriptors carry the already planned ordered-source
     digest. Full validation recomputes it from the output's source table. */
  if (sample.contains("_source_count")) return sample.at("_package_binding_sha256");
  return bindingWithSources(plan,sample,[&](std::ostream &out) {
    out<<'['; size_t count=0;
    for (const auto &src : sample.at("sources")) {
      if (count%1000==0) ReductionProgress::pulse("hashing_source_bindings");
      if (count++) out<<',';
      out<<Json{{"identity",src.at("identity")},{"sha256",src.at("expected_sha256")},
        {"metadata",src.at("metadata")}}.dump(-1,' ',true);
    }
    out<<']';
  });
}

Json sampleHeader(const Json &sample) {
  Json out=Json::object();
  for (const auto &[key,value] : sample.items()) if (key!="sources") out[key]=value;
  out["sources"]=Json::array(); return out;
}

void RecordTable::append(TTree &input, const std::string &name, Long64_t source, Long64_t offset,
                         TDirectory *dest, const std::vector<Long64_t> *sourceMap) {
  if (!output) {
    dest->cd(); output.reset(new TTree(name.c_str(), input.GetTitle()));
    output->SetAutoFlush(-8*1024*1024);
    for (auto *obj : *input.GetListOfLeaves()) {
      auto *leaf = static_cast<TLeaf *>(obj);
      const std::string n=leaf->GetName(), t=leaf->GetTypeName();
      if (t=="Double_t" || t=="Float_t") output->Branch(n.c_str(), &(doubles[n]=0), (n+"/D").c_str());
      else output->Branch(n.c_str(), &(ints[n]=0), (n+"/L").c_str());
    }
    if (!ints.count("source")) output->Branch("source", &(ints["source"]=0), "source/L");
  }
  for (Long64_t row=0; row<input.GetEntries(); ++row) {
    if (row%10000==0) ReductionProgress::pulse("copying_"+name);
    if (input.GetEntry(row)<=0) throw std::runtime_error("cannot read event record");
    for (auto &[n,v] : doubles) {
      auto *leaf=input.GetLeaf(n.c_str());
      if (!leaf) throw std::runtime_error("event-record schema changed");
      v=leaf->GetValue();
    }
    for (auto &[n,v] : ints) {
      if (n=="source" && !sourceMap) {v=source; continue;}
      auto *leaf=input.GetLeaf(n.c_str());
      if (!leaf) throw std::runtime_error("event-record schema changed");
      v=leaf->GetValueLong64();
      if (n=="event") v+=offset;
      if (n=="file" || n=="source") {
        if (sourceMap) {
          if (v<0 || static_cast<size_t>(v)>=sourceMap->size()) throw std::runtime_error("invalid event source reference");
          v=sourceMap->at(v);
        } else v=source;
      }
    }
    output->Fill();
  }
}
void RecordTable::write() {
  if (!output) return;
  output->GetDirectory()->cd(); output->Write("",TObject::kOverwrite); output->SetDirectory(nullptr);
}

/* Validate against the immutable input snapshot, including every original
   source and its full, gap-free range union. Recreate only a missing sidecar
   after a valid ROOT publication; existing sidecars must match exactly. */
namespace {
void checkPackageMetadata(const Json &meta,const Json &plan,const Json &sample,const std::string &path) {
  const std::string binding=packageBinding(plan,sample);
  if (meta.at("schema")!="PhotonJetSampleV1" || meta.at("system")!=sample.at("system") ||
      meta.at("name")!=sample.at("name") || meta.at("family")!=sample.at("family") ||
      meta.at("selection")!=plan.at("selection") || meta.at("histograms")!=plan.at("histograms") ||
      (meta.at("plan_sha256")!=plan.at("plan_sha256") &&
       (sample.value("intermediate",false) || sample.contains("run_packages") ||
        meta.value("package_binding_sha256","")!=binding)) || meta.at("coverage")!=sample.at("coverage") ||
      meta.at("source_count")!=sourceCount(sample) ||
      meta.at("package_kind")!=(sample.value("intermediate",false) ? "INTERMEDIATE" : "FINAL"))
    throw std::runtime_error("package identity, configuration, coverage or snapshot differs: " + path);
  if (meta.contains("package_binding_sha256") && meta.at("package_binding_sha256")!=binding)
    throw std::runtime_error("package input/settings binding differs: " + path);
  if (meta.at("centrality_weights")!=sample.at("weights") || meta.at("vertex_reweighting")!=sample.at("vertex_reweighting") ||
      meta.at("render_selection").at("producer_code")!=plan.at("producer_code"))
    throw std::runtime_error("package correction or producer binding differs");
  if (meta.value("run",Json())!=sample.value("run",Json()) ||
      meta.value("run_packages",Json())!=sample.value("run_packages",Json()) ||
      meta.value("data_run",Json())!=sample.value("data_run",Json()) ||
      meta.value("run_coverages",Json())!=sample.value("run_coverages",Json()) ||
      meta.value("run_records",Json())!=sample.value("run_records",Json()))
    throw std::runtime_error("DATA run or combined-run binding differs");
  Json normal=meta.at("normalization"), expected=sample.at("normalization");
  if (expected.contains("ownership")) {
    const long long denominator=normal.at("denominator_events");
    const Json weight=denominator>0 ? Json(expected.at("cross_section_pb").get<double>()/denominator) : Json();
    if (normal.at("pb_per_generated_event")!=weight || (denominator<=0 && !sample.value("intermediate",false)))
      throw std::runtime_error("stitching normalization differs");
    normal.erase("denominator_events"); normal["pb_per_generated_event"]=nullptr;
  }
  if (normal!=expected || normal.value("applied_in_histograms",false)) throw std::runtime_error("sample normalization differs");
}

/* A completed package has already passed the full record scan. Verify its
   exact bytes and binding before reuse; old receipts require a full check. */
void checkValidatedInput(const std::string &path,const Json &plan,const Json &sample) {
  fs::path jsonPath=path; jsonPath.replace_extension(".json");
  if (!fs::exists(jsonPath) || !readJson(jsonPath.string()).contains("validation")) {
    validatePackage(path,plan,sample); return;
  }
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open package " + path);
  checkPackageMetadata(object(*file,"meta/package"),plan,sample,path);
  if (readJson(jsonPath.string())!=sidecar(*file,path))
    throw std::runtime_error("validated package bytes or receipt changed: " + path);
}
}  // namespace

void validatePackage(const std::string &path, const Json &plan, const Json &sample) {
  ReductionProgress progress("validate",sourceCount(sample)+1);
  checkReductionCode(plan); checkVertexFit(sample); checkPpPeriodFits(ppPeriodsOf(plan,sample));
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot open package " + path);
  const Json meta=object(*f,"meta/package");
  checkPackageMetadata(meta,plan,sample,path);
  size_t histograms=0;
  std::map<std::string,TH1*> dense;
  for (const auto &row : inventory(sample.at("system"),sample.at("family"))) {
    const std::string name=row.at("path"); TObject *o=f->Get(name.c_str());
    if (!o && row.at("kind")=="extensions") continue;
    if (!o) throw std::runtime_error("missing histogram " + name);
    finiteHistogram(*o,name); ++histograms;
    if (auto *h=dynamic_cast<TH1*>(o)) dense[name]=h;
  }
  if (meta.at("histograms_written")!=histograms) throw std::runtime_error("histogram inventory differs");
  TTree &sources=tree(*f,"records/sources"), &cells=tree(*f,"records/source_contributions"), &events=tree(*f,"records/events");
  if (sources.GetEntries()!=static_cast<Long64_t>(sourceCount(sample)) || meta.at("event_records")!=events.GetEntries())
    throw std::runtime_error("source/event counts differ");
  std::string *json=nullptr, *histogram=nullptr;
  sources.SetBranchAddress("json",&json); cells.SetBranchAddress("histogram",&histogram);
  Long64_t previous=0;
  const std::string observedBinding=bindingWithSources(plan,sample,[&](std::ostream &out) {
  out<<'[';
  for (Long64_t i=0; i<sources.GetEntries(); ++i) {
    if (sources.GetEntry(i)<=0 || sources.GetLeaf("source")->GetValueLong64()!=i) throw std::runtime_error("invalid source record");
    const Json src=Json::parse(*json);
    if (!sample.contains("_source_count")) {
      const Json &wanted=sample.at("sources").at(i);
      if (src.at("identity")!=wanted.at("identity") || src.at("sha256")!=wanted.at("expected_sha256") || src.at("metadata")!=wanted.at("metadata"))
        throw std::runtime_error("original source identity, input hash or metadata differs");
    }
    if (i) out<<',';
    out<<Json{{"identity",src.at("identity")},{"sha256",src.at("sha256")},
      {"metadata",src.at("metadata")}}.dump(-1,' ',true);
    long long next=0;
    for (const auto &range : src.at("entry_ranges")) {
      if (range.size()!=2 || range[0].get<long long>()!=next || range[1].get<long long>()<=next) throw std::runtime_error("invalid source range union");
      next=range[1].get<long long>();
    }
    if (next!=src.at("metadata").at("entries").get<long long>() || src.at("reduced_entries")!=next)
      throw std::runtime_error("incomplete original-source coverage");
    const Long64_t begin=sources.GetLeaf("contribution_begin")->GetValueLong64(), end=sources.GetLeaf("contribution_end")->GetValueLong64();
    if (begin!=previous || end<begin || end>cells.GetEntries()) throw std::runtime_error("invalid contribution range");
    std::pair<std::string,Long64_t> prior; bool first=true;
    for (Long64_t c=begin;c<end;++c) {
      if (c%10000==0) ReductionProgress::pulse("validating_contributions");
      if (cells.GetEntry(c)<=0 || cells.GetLeaf("source")->GetValueLong64()!=i)
        throw std::runtime_error("contribution source differs");
      const Long64_t bin=cells.GetLeaf("bin")->GetValueLong64();
      const auto found=dense.find(*histogram);
      auto *h=found==dense.end() ? nullptr : found->second;
      const double value=cells.GetLeaf("value")->GetValue(), variance=cells.GetLeaf("variance")->GetValue();
      const std::pair<std::string,Long64_t> key{*histogram,bin};
      if (!h || bin<0 || bin>=h->GetNcells() || !std::isfinite(value) || !std::isfinite(variance) || variance<0 || (!first && key<=prior))
        throw std::runtime_error("invalid, duplicate or unordered source cell");
      prior=key; first=false;
    }
    previous=end; progress.advance(i+1);
  }
  out<<']';
  });
  sources.ResetBranchAddresses(); cells.ResetBranchAddresses();
  if (observedBinding!=packageBinding(plan,sample)) throw std::runtime_error("ordered source identities or metadata differ from the planned binding");
  if (previous!=cells.GetEntries()) throw std::runtime_error("unattributed source contributions");
  std::vector<std::pair<Long64_t,Long64_t>> eventRanges(sources.GetEntries(),{0,0});
  Long64_t lastSource=-1;
  for (const char *name : {"records/events","records/photons","records/jets"}) {
    TTree &t=tree(*f,name);
    if (!t.GetLeaf("source")) throw std::runtime_error("missing event source reference");
    t.SetBranchStatus("*",false); t.SetBranchStatus("source",true);
    if (t.GetLeaf("event")) t.SetBranchStatus("event",true);
    for (Long64_t i=0; i<t.GetEntries(); ++i) {
      if (i%10000==0) ReductionProgress::pulse("validating_event_references");
      if (t.GetEntry(i)<=0) throw std::runtime_error("unreadable event record");
      const long long id=t.GetLeaf("source")->GetValueLong64();
      if (id<0 || id>=sources.GetEntries()) throw std::runtime_error("event source outside source table");
      if (std::string(name)=="records/events") {
        if (id<lastSource) throw std::runtime_error("unordered event source records");
        if (id!=lastSource) eventRanges.at(id).first=i;
        eventRanges.at(id).second=i+1; lastSource=id;
      } else if (auto *e=t.GetLeaf("event")) {
        const auto range=eventRanges.at(id);
        if (e->GetValueLong64()<range.first || e->GetValueLong64()>=range.second)
          throw std::runtime_error("event reference belongs to another source or is outside event table");
      } else throw std::runtime_error("missing object event reference");
    }
  }
  const Json receipt=sidecar(*f,path);
  fs::path jsonPath=path; jsonPath.replace_extension(".json");
  if (fs::exists(jsonPath) && readJson(jsonPath.string())!=receipt) throw std::runtime_error("package sidecar differs: " + jsonPath.string());
  if (!fs::exists(jsonPath)) writeJson(jsonPath,receipt);
  progress.finish();
}

void publishValidatedPackage(const std::string &partial,const std::string &path,
                             const Json &plan,const Json &sample) {
  fs::path temporary=partial; temporary.replace_extension(".json");
  fs::remove(temporary);  // derived receipt of an unpublished attempt
  validatePackage(partial,plan,sample);
  Json receipt=readJson(temporary.string()); receipt["package"]=fs::path(path).filename().string();
  fs::path side=path; side.replace_extension(".json");
  if (fs::exists(side) && readJson(side.string())!=receipt) throw std::runtime_error("existing publication receipt differs");
  /* Atomic rename preserves the validated bytes. A crash before the sidecar
     write is recovered by the ordinary full validation of the final file. */
  fs::rename(partial,path); writeJson(side,receipt); fs::remove(temporary);
}

/* Only final packages with an exact binding may survive an expanded input
   snapshot. The combined DATA product is rebuilt from the selected runs.
   Old packages without this binding remain readable but are never reused. */
void reusePackages(Json &plan, const std::string &collection) {
  const fs::path manifest=fs::canonical(collection), base=manifest.parent_path();
  const Json previous=photonjet_config::readYaml(manifest.string());
  if (previous.value("schema","")!="PhotonJetHistogramCollectionV1" || previous.value("status","")!="VALIDATED")
    throw std::runtime_error("reuse_collection must name a validated collection.yaml");
  std::map<std::string,Json> samples;
  std::map<std::string,std::string> hashes;
  for (const auto &row : previous.at("samples")) {
    const std::string key=row.at("key");
    if (!samples.emplace(key,row).second) throw std::runtime_error("duplicate sample in reuse collection");
  }
  for (const auto &file : previous.at("files"))
    if (!hashes.emplace(file.at("path"),file.at("sha256")).second) throw std::runtime_error("duplicate reuse payload");
  std::set<std::string> reused;
  for (auto &sample : plan.at("samples")) {
    const std::string key=sample.at("key"); const auto found=samples.find(key);
    if (sample.contains("run_packages") || found==samples.end() ||
        found->second.value("package_binding_sha256","")!=packageBinding(plan,sample)) continue;
    const Json &row=found->second; const std::string root=key+".root", side=key+".json";
    const fs::path path=fs::canonical(base/root), json=fs::canonical(base/side);
    if (path.lexically_relative(base).empty() || *path.lexically_relative(base).begin()==".." ||
        json.lexically_relative(base).empty() || *json.lexically_relative(base).begin()==".." ||
        row.at("path")!=root || !hashes.count(root) || !hashes.count(side) ||
        hashes.at(root)!=row.at("sha256").get<std::string>() || sha256File(path.string())!=hashes.at(root) || sha256File(json.string())!=hashes.at(side))
      throw std::runtime_error("reused collection payload changed or escaped its directory: " + key);
    validatePackage(path.string(),plan,sample);
    sample["reused_package"]={{"path",path.string()},{"sha256",row.at("sha256")},
      {"package_plan_sha256",row.value("package_plan_sha256",previous.at("input_snapshot_sha256"))}};
    reused.insert(key);
  }
  Json jobs=Json::array();
  for (const auto &job : plan.at("jobs")) if (!reused.count(job.at("sample"))) jobs.push_back(job);
  plan["jobs"]=jobs;
  plan["reuse_collection"]={{"path",manifest.string()},{"sha256",sha256File(manifest.string())},{"reused_samples",reused.size()}};
}

/* Merge a bounded number of disjoint packages. Record/source offsets move
   together; no event or source payload is materialized in memory. */
void combinePackages(const Json &plan, const Json &sample, const std::vector<std::string> &inputs,
                     const std::string &output,const Json &inputSamples) {
  ReductionProgress progress("merge",sourceCount(sample)+1);
  checkReductionCode(plan); checkVertexFit(sample); checkPpPeriodFits(ppPeriodsOf(plan,sample));
  if (inputs.empty()) throw std::runtime_error("no packages to combine");
  if (!inputSamples.empty() && inputSamples.size()!=inputs.size()) throw std::runtime_error("merge parent bindings differ");
  if (fs::exists(output)) {validatePackage(output,plan,sample); progress.finish(); return;}
  fs::create_directories(fs::path(output).parent_path());
  const std::string partial=output+".partial";
  fs::remove(partial);
  std::unique_ptr<TFile> out(TFile::Open(partial.c_str(),"CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("cannot create merged package");
  TDirectory *records=directory(*out,"records"); records->cd();
  Long64_t source=0, bin=0, begin=0, end=0, eventOffset=0;
  double value=0,variance=0; std::string histogram,json;
  TTree sources("sources","Original source identities and contribution ranges"), cells("source_contributions","Sparse cells of each original source");
  sources.Branch("source",&source); sources.Branch("json",&json);
  sources.Branch("contribution_begin",&begin); sources.Branch("contribution_end",&end);
  cells.Branch("source",&source); cells.Branch("histogram",&histogram); cells.Branch("bin",&bin);
  cells.Branch("value",&value); cells.Branch("variance",&variance);
  sources.SetAutoFlush(-8*1024*1024); cells.SetAutoFlush(-8*1024*1024);
  RecordTable events,photons,jets;
  std::map<std::string,std::unique_ptr<TH1>> sums;
  std::map<std::string,std::unique_ptr<THnSparse>> sparse;
  Json meta,available=Json::object(); std::string grids;
  long long denominator=0; size_t cursor=0;
  for (size_t inputIndex=0; inputIndex<inputs.size(); ++inputIndex) {
    const std::string &path=inputs[inputIndex];
    std::unique_ptr<TFile> in(TFile::Open(path.c_str(),"READ"));
    if (!in || in->IsZombie()) throw std::runtime_error("cannot open merge input " + path);
    const Json m=object(*in,"meta/package");
    const size_t count=m.at("source_count");
    if (!count || cursor+count>sourceCount(sample)) throw std::runtime_error("merge source count differs");
    Json subset;
    if (!inputSamples.empty()) {
      subset=inputSamples.at(inputIndex);
      if (sourceCount(subset)!=count) throw std::runtime_error("merge parent source count differs");
    } else {
      subset=sampleHeader(sample); subset["intermediate"]=true;
      for (size_t i=0; i<count; ++i) subset["sources"].push_back(sample.at("sources").at(cursor+i));
    }
    if (m.contains("run") && inputSamples.empty()) {
      const std::string key=m.at("system").get<std::string>()+"/"+m.at("name").get<std::string>();
      if (sample.contains("run_packages")) {
        if (std::find(sample.at("run_packages").begin(),sample.at("run_packages").end(),Json(key))==sample.at("run_packages").end())
          throw std::runtime_error("unplanned DATA run package " + key);
        subset["key"]=key; subset["name"]=m.at("name"); subset["run"]=m.at("run"); subset.erase("run_packages");
        subset["coverage"]=sample.at("run_coverages").at(key); subset.erase("run_coverages");
        subset["data_run"]=sample.at("run_records").at(key); subset.erase("run_records");
      }
      subset["intermediate"]=m.at("package_kind")=="INTERMEDIATE";
      const int run=subset.at("run");
      for (const auto &src : subset.at("sources")) if (src.at("metadata").at("source").value("run",-1)!=run)
        throw std::runtime_error("DATA source belongs to another run");
    }
    checkValidatedInput(path,plan,subset);
    if (meta.is_null()) meta=m;
    if (m.at("render_selection")!=meta.at("render_selection")) throw std::runtime_error("merge render/selection bindings differ");
    if (m.at("normalization").contains("ownership")) denominator+=m.at("normalization").at("denominator_events").get<long long>();
    for (const auto &[key,state] : m.at("availability").items()) {
      if (!available.contains(key)) available[key]=state;
      else if (available.at(key)!=state) available[key]="PARTIAL";
    }
    for (const auto &row : inventory(sample.at("system"),sample.at("family"))) {
      const std::string name=row.at("path"); TObject *o=in->Get(name.c_str());
      if (!o && row.at("kind")=="extensions") continue;
      if (!o) throw std::runtime_error("missing merge histogram " + name);
      if (auto *h=dynamic_cast<TH1 *>(o)) {
        if (!sums.count(name)) {sums[name].reset(static_cast<TH1 *>(h->Clone())); sums[name]->SetDirectory(nullptr);}
        else {
          sameAxis(*sums[name]->GetXaxis(),*h->GetXaxis());
          if (h->GetDimension()>1) sameAxis(*sums[name]->GetYaxis(),*h->GetYaxis());
          if (h->GetDimension()>2) sameAxis(*sums[name]->GetZaxis(),*h->GetZaxis());
          if (sums[name]->GetDimension()!=h->GetDimension() || !sums[name]->Add(h)) throw std::runtime_error("histogram merge failed");
        }
      } else if (auto *h=dynamic_cast<THnSparse *>(o)) {
        if (!sparse.count(name)) sparse[name].reset(static_cast<THnSparse *>(h->Clone()));
        else {
          if (sparse[name]->GetNdimensions()!=h->GetNdimensions()) throw std::runtime_error("sparse dimensions differ");
          for (int d=0;d<h->GetNdimensions();++d) sameAxis(*sparse[name]->GetAxis(d),*h->GetAxis(d));
          sparse[name]->Add(h);
        }
      }
    }
    if (auto *g=dynamic_cast<TObjString *>(in->Get("grids"))) {
      const std::string current=g->GetString().Data();
      if (!grids.empty() && grids!=current) throw std::runtime_error("response grids differ"); grids=current;
    }
    TTree &is=tree(*in,"records/sources"), &ic=tree(*in,"records/source_contributions");
    std::string *sj=nullptr,*ch=nullptr; is.SetBranchAddress("json",&sj); ic.SetBranchAddress("histogram",&ch);
    std::vector<Long64_t> map;
    for (Long64_t i=0;i<is.GetEntries();++i) {
      if (is.GetEntry(i)<=0) throw std::runtime_error("unreadable source record");
      source=cursor+i; map.push_back(source); json=*sj; begin=cells.GetEntries();
      const Long64_t cb=is.GetLeaf("contribution_begin")->GetValueLong64(), ce=is.GetLeaf("contribution_end")->GetValueLong64();
      for (Long64_t c=cb;c<ce;++c) {
        if (c%10000==0) ReductionProgress::pulse("copying_contributions");
        if (ic.GetEntry(c)<=0 || ic.GetLeaf("source")->GetValueLong64()!=i) throw std::runtime_error("contribution source differs");
        histogram=*ch; bin=ic.GetLeaf("bin")->GetValueLong64(); value=ic.GetLeaf("value")->GetValue(); variance=ic.GetLeaf("variance")->GetValue(); cells.Fill();
      }
      end=cells.GetEntries(); sources.Fill();
    }
    is.ResetBranchAddresses(); ic.ResetBranchAddresses();
    TTree &ev=tree(*in,"records/events"), &ph=tree(*in,"records/photons"), &jt=tree(*in,"records/jets");
    events.append(ev,"events",0,eventOffset,records,&map); photons.append(ph,"photons",0,eventOffset,records,&map); jets.append(jt,"jets",0,eventOffset,records,&map);
    eventOffset+=ev.GetEntries(); cursor+=count; progress.advance(cursor);
  }
  if (cursor!=sourceCount(sample)) throw std::runtime_error("merge did not cover the planned sources");
  for (auto &[name,h] : sums) {const auto slash=name.rfind('/'); directory(*out,name.substr(0,slash))->cd(); h->Write(name.substr(slash+1).c_str());}
  for (auto &[name,h] : sparse) {const auto slash=name.rfind('/'); directory(*out,name.substr(0,slash))->cd(); h->Write(name.substr(slash+1).c_str());}
  records->cd(); sources.Write(); cells.Write(); events.write(); photons.write(); jets.write();
  meta["source_count"]=cursor; meta["event_records"]=eventOffset; meta["histograms_written"]=sums.size()+sparse.size(); meta["availability"]=available;
  meta["system"]=sample.at("system"); meta["name"]=sample.at("name"); meta["family"]=sample.at("family"); meta["coverage"]=sample.at("coverage");
  meta["plan_sha256"]=plan.at("plan_sha256"); meta["package_binding_sha256"]=packageBinding(plan,sample);
  meta.erase("run"); meta.erase("run_packages"); meta.erase("data_run"); meta.erase("run_coverages"); meta.erase("run_records");
  if (sample.contains("run")) meta["run"]=sample.at("run");
  if (sample.contains("data_run")) meta["data_run"]=sample.at("data_run");
  if (sample.contains("run_packages")) {meta["run_packages"]=sample.at("run_packages");meta["run_coverages"]=sample.at("run_coverages");meta["run_records"]=sample.at("run_records");meta["excluded_runs"]=sample.value("excluded_runs",Json::array());}
  meta["package_kind"]=sample.value("intermediate",false) ? "INTERMEDIATE" : "FINAL";
  if (meta.at("normalization").contains("ownership")) {
    if (denominator<=0 && !sample.value("intermediate",false)) throw std::runtime_error("no stitching denominator events");
    meta["normalization"]["denominator_events"]=denominator;
    meta["normalization"]["pb_per_generated_event"]=denominator>0 ? Json(meta.at("normalization").at("cross_section_pb").get<double>()/denominator) : Json();
  }
  out->cd(); TObjString(plan.at("selection").at("profile").get<std::string>().c_str()).Write("profile");
  TObjString(plan.at("selection").dump().c_str()).Write("selection"); if (!grids.empty()) TObjString(grids.c_str()).Write("grids");
  directory(*out,"meta")->cd(); TObjString(meta.dump().c_str()).Write("package");
  TObjString(inventory(sample.at("system"),sample.at("family")).dump().c_str()).Write("inventory");
  if (sample.at("family")!="data" && !meta.at("normalization").at("pb_per_generated_event").is_null())
    TParameter<double>("pb_per_generated_event",meta.at("normalization").at("pb_per_generated_event").get<double>()).Write();
  sources.SetDirectory(nullptr); cells.SetDirectory(nullptr); out->Close();
  sums.clear(); sparse.clear(); events.output.reset(); photons.output.reset(); jets.output.reset();
  publishValidatedPackage(partial,output,plan,sample); progress.finish();
}
}  // namespace treetohists
