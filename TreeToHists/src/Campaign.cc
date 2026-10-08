/* One farm launch: freeze the published input inventory, plan on a compute
   node, fill bounded rows, merge disjoint sources, then publish a collection.
   This coordinator never runs FinalAnalysis or transfers files. */
#include "TreeToHists.h"
#include "Configuration.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <TFile.h>
#include <TObjString.h>
#include <TROOT.h>

namespace treetohists {
namespace fs=std::filesystem;
namespace {
std::string arg(const std::vector<std::string> &args,const std::string &key) {
  for (size_t i=0;i+1<args.size();++i) if (args[i]==key) return args[i+1];
  throw std::runtime_error("missing " + key);
}
void atomic(const fs::path &path,const std::string &body) {
  fs::create_directories(path.parent_path()); const fs::path next=path.string()+".partial";
  std::ofstream stream(next); stream<<body; stream.close();
  if (!stream) throw std::runtime_error("cannot write " + path.string()); fs::rename(next,path);
}
void yaml(std::ostream &out,const Json &j,int indent=0) {
  const std::string pad(indent,' ');
  if (j.is_object()) for (const auto &[key,value] : j.items()) {
    out<<pad<<Json(key).dump()<<":";
    if ((value.is_object() || value.is_array()) && !value.empty()) {out<<'\n'; yaml(out,value,indent+2);}
    else out<<" "<<value.dump()<<'\n';
  } else if (j.is_array()) for (const auto &value : j) {
    out<<pad<<"-";
    if ((value.is_object() || value.is_array()) && !value.empty()) {out<<'\n'; yaml(out,value,indent+2);}
    else out<<" "<<value.dump()<<'\n';
  }
}
void writeYaml(const fs::path &path,const Json &j) {std::ostringstream body; yaml(body,j); atomic(path,body.str());}
std::string utc() {
  const std::time_t now=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm time{}; gmtime_r(&now,&time); std::ostringstream out; out<<std::put_time(&time,"%Y-%m-%dT%H:%M:%SZ"); return out.str();
}
std::string safePath(const fs::path &base,const Json &value) {
  const fs::path path=fs::absolute(base/value.get<std::string>()).lexically_normal();
  const std::string text=path.string();
  if (!std::regex_match(text,std::regex("/[A-Za-z0-9_./+@-]+"))) throw std::runtime_error("campaign paths need ordinary absolute Unix path names: " + text);
  return text;
}
Json sampleOf(const Json &p,const std::string &key) {
  for (const auto &s : p.at("samples")) if (s.at("key")==key) return s;
  throw std::runtime_error("unknown sample " + key);
}
Json filePin(const fs::path &path) {
  return {{"path",path.string()},{"sha256",sha256File(path.string())},{"size_bytes",fs::file_size(path)}};
}
Json inspectionParts(const Json &c) {
  const fs::path path=fs::path(c.at("work").get<std::string>())/"inspection/index.json";
  if (sha256File(path.string())!=c.at("inspection_index_sha256").get<std::string>())
    throw std::runtime_error("input inspection index changed");
  return readJson(path.string());
}
Json inspectionSnapshot(const Json &c) {
  Json snapshot=c.at("resolved_settings");snapshot["producer_code"]=c.at("producer_code");
  snapshot["root_version"]=gROOT->GetVersion();return snapshot;
}
size_t sourceCount(const Json &s) {
  return s.contains("_source_count") ? s.at("_source_count").get<size_t>() : s.at("sources").size();
}
Json workerHeader(const Json &s) {
  Json out=sampleHeader(s);
  for (const char *key : {"_source_count","_package_binding_sha256","_source_file","_source_files"}) out.erase(key);
  return out;
}
struct SourceSpan {std::string path; long long offset=0; size_t count=0;};
void spanSources(const std::vector<SourceSpan> &spans,const std::function<void(const Json &)> &visit) {
  for (const auto &span : spans) {
    std::ifstream in(span.path); in.seekg(span.offset); std::string line;
    for (size_t i=0;i<span.count;++i) {
      if (!std::getline(in,line)) throw std::runtime_error("planned source inventory is truncated");
      visit(Json::parse(line));
    }
  }
}
Json spanSample(const Json &plan,const Json &sample,const std::vector<SourceSpan> &spans,bool intermediate) {
  Json out=workerHeader(sample); out["intermediate"]=intermediate;
  size_t count=0; for (const auto &span : spans) count+=span.count;
  out["_source_count"]=count;
  out["_package_binding_sha256"]=packageBindingFromSources(plan,out,[&](std::ostream &stream) {
    stream<<'['; size_t seen=0;
    spanSources(spans,[&](const Json &src) {
      if (seen%1000==0) ReductionProgress::pulse("hashing_source_bindings");
      if (seen++) stream<<',';
      stream<<Json{{"identity",src.at("identity")},{"sha256",src.at("expected_sha256")},
        {"metadata",src.at("metadata")}}.dump(-1,' ',true);
    });
    if (seen!=count) throw std::runtime_error("planned source span count differs");
    stream<<']';
  });
  return out;
}
long fillPasses(const Json &sample) {
  return sample.at("family")=="data" ? 3 : 4+(sample.at("normalization").contains("ownership") ? 1 : 0);
}
void bindings(const Json &c,bool checkInputs=true) {
  if (c.at("producer_code")!=reductionCode() || c.at("binary_sha256").get<std::string>()!=sha256File(c.at("binary")) ||
      c.at("entrypoint_sha256").get<std::string>()!=sha256File(c.at("entrypoint")) ||
      c.at("setup_sha256").get<std::string>()!=sha256File(c.at("setup_script")) ||
      c.at("settings_sha256").get<std::string>()!=sha256File(c.at("settings")))
    throw std::runtime_error("campaign code, build, setup, settings or frozen inputs changed");
  if (c.at("resolved_settings")!=settingsSnapshot(c.at("settings"),c.at("home")))
    throw std::runtime_error("frozen physics definitions changed");
  /* Planning verifies the complete input snapshot once. A merge reads only
     its sealed descriptor and validated parent packages, whose bindings and
     payload hashes are checked by the reducer. Rehashing all manifests for
     every merge would multiply shared-filesystem work without checking any
     additional bytes used by that worker. */
  if (checkInputs) {
    if (c.at("inputs_sha256").get<std::string>()!=sha256File(c.at("inputs"))) throw std::runtime_error("frozen inputs changed");
    for (const auto &pin : c.at("publication_bindings"))
      if (sha256File(pin.at("path"))!=pin.at("sha256").get<std::string>()) throw std::runtime_error("published tree manifest changed");
    if (!c.at("timing").is_null() && sha256File(c.at("timing_path"))!=c.at("timing_sha256").get<std::string>())
      throw std::runtime_error("timing evidence changed");
  }
}
Json timing(const std::string &path) {
  if (path=="none") return nullptr;
  Json t=photonjet_config::readYaml(path);
  if (t.at("schema")!="TreeToHistsTimingV1" || t.at("root_version").get<std::string>()!=gROOT->GetVersion() ||
      t.at("producer_code_sha256")!=digest(reductionCode()) || t.at("classes").empty())
    throw std::runtime_error("timing receipt does not bind this ROOT build and reduction code");
  for (auto &[name,row] : t.at("classes").items()) {
    if (row.at("status")!="MEASURED_TERMINAL" || row.at("processed_entries").get<long long>()<=0 ||
        !(row.at("seconds").get<double>()>0) || !std::isfinite(row.at("seconds").get<double>()))
      throw std::runtime_error("terminal measured timing required for " + name);
    const std::string receipt=row.at("receipt");
    if (sha256File(receipt)!=row.at("receipt_sha256").get<std::string>()) throw std::runtime_error("timing source receipt changed");
    const Json measured=readJson(receipt);
    if (measured.value("sample_class", "")!=name)
      throw std::runtime_error("timing receipt belongs to another sample class: " + name);
    for (const char *field : {"producer_code_sha256","root_version","processed_entries","seconds","reduction_settings_sha256"})
      if (measured.at(field)!=(row.contains(field) ? row.at(field) : t.at(field))) throw std::runtime_error("timing does not match the terminal fill receipt");
  }
  return t;
}
std::string submit(const Json &c,const std::string &arguments,const std::string &label) {
  const std::string logs=c.at("logs"), work=c.at("work"), setup=c.at("setup");
  std::ostringstream out;
  out<<"universe = vanilla\nexecutable = "<<c.at("entrypoint").get<std::string>()
    <<"\narguments = "<<arguments<<"\ngetenv = False\nenvironment = \"TREETOHISTS_SETUP='"<<setup<<"'\""
    <<"\ninitialdir = "<<work<<"\nshould_transfer_files = NO\nrequest_cpus = 1\nrequest_memory = "
    <<c.at(label=="fill" || label=="inspect" ? "fill_memory_mb" : label=="merge" ? "merge_memory_mb" : label=="collection" ? "collection_memory_mb" : "plan_memory_mb")
    <<"MB\nrequest_disk = "<<c.at("disk_mb")<<"MB\non_exit_hold = False\non_exit_remove = True"
    <<"\nlog = "<<c.at("events").get<std::string>()<<"\noutput = "<<logs<<"/$(shard)/"<<label<<"_$(row).out"
    <<"\nerror = "<<logs<<"/$(shard)/"<<label<<"_$(row).err\nqueue 1\n";
  return out.str();
}

/* Each final merge validates its ROOT package before publication. The last
   job needs only these small expected bindings and the validated receipts,
   not the campaign's complete source/range inventory in memory. */
Json collectionIndex(const Json &plan) {
  Json samples=Json::array();
  for (const auto &s : plan.at("samples")) {
    Json row={{"key",s.at("key")},{"system",s.at("system")},{"name",s.at("name")},
      {"family",s.at("family")},{"coverage",s.at("coverage")},
      {"source_count",sourceCount(s)},{"package_binding_sha256",packageBinding(plan,s)}};
    for (const char *key : {"run","run_packages"}) if (s.contains(key)) row[key]=s.at(key);
    if (s.at("system")=="PP" && s.at("family")!="data") {
      std::set<std::string> roles;
      if (s.contains("si_di_role")) roles.insert(s.at("si_di_role"));
      else for (const auto &src : s.at("sources")) roles.insert(src.at("metadata").at("source").at("si_di_role"));
      if (roles.size()!=1) throw std::runtime_error("mixed p+p SI/DI source roles in one sample");
      row["si_di_role"]=*roles.begin(); row["simulation_periods"]=ppPeriodsOf(plan,s).at("mode");
    }
    samples.push_back(std::move(row));
  }
  Json index={{"schema","TreeToHistsCollectionIndexV1"},{"plan_sha256",plan.at("plan_sha256")},
    {"producer_code",plan.at("producer_code")},{"root_version",plan.at("root_version")},
    {"selection",plan.at("selection")},{"histograms",plan.at("histograms")},{"samples",samples}};
  index["descriptor_sha256"]=digest(index); return index;
}
void prepare(const std::vector<std::string> &args) {
  const fs::path path=fs::absolute(arg(args,"--config")); const fs::path base=path.parent_path();
  const Json cfg=photonjet_config::readYaml(path.string()).at("Campaign");
  const std::set<std::string> known={"Settings","Inputs","WorkDirectory","OutputDirectory","LogsDirectory","EventLog",
    "SetupScript","SetupArguments","TimingFile","MaxEventsPerJob","MaxSourcesPerJob","MergeFanIn","FillConcurrency","MergeConcurrency",
    "MemoryMB","PlanMemoryMB","FillMemoryMB","MergeMemoryMB","CollectionMemoryMB","DiskMB","InspectSourcesPerJob","InspectConcurrency"};
  for (const auto &[key,value] : cfg.items()) if (!known.count(key)) throw std::runtime_error("unknown Campaign setting: " + key);
  Json c={{"schema","TreeToHistsCampaignV1"},{"home",safePath(base,arg(args,"--home"))},
    {"work",safePath(base,cfg.at("WorkDirectory"))},{"output",safePath(base,cfg.at("OutputDirectory"))},
    {"logs",safePath(base,cfg.at("LogsDirectory"))},{"events",safePath(base,cfg.at("EventLog"))},
    {"setup_script",safePath(base,cfg.at("SetupScript"))}, {"producer_code",reductionCode()},
    {"max_events",cfg.value("MaxEventsPerJob",50000L)},{"max_sources",cfg.value("MaxSourcesPerJob",32)},{"fan_in",cfg.value("MergeFanIn",200)},
    {"inspect_sources",cfg.value("InspectSourcesPerJob",500)},{"inspect_concurrency",cfg.value("InspectConcurrency",128)},
    {"fill_concurrency",cfg.value("FillConcurrency",20000)},{"merge_concurrency",cfg.value("MergeConcurrency",128)},
    {"plan_memory_mb",cfg.value("PlanMemoryMB",cfg.value("MemoryMB",4096))},
    {"fill_memory_mb",cfg.value("FillMemoryMB",cfg.value("MemoryMB",2048))},
    {"merge_memory_mb",cfg.value("MergeMemoryMB",cfg.value("MemoryMB",2048))},
    {"collection_memory_mb",cfg.value("CollectionMemoryMB",cfg.value("MemoryMB",1024))},{"disk_mb",cfg.value("DiskMB",8192)}};
  if (c.at("max_events").get<long>()<1 || c.at("max_events").get<long>()>100000 ||
      c.at("max_sources").get<int>()<1 || c.at("max_sources").get<int>()>256 ||
      c.at("fan_in").get<int>()<2 || c.at("fan_in").get<int>()>1024 || c.at("fill_concurrency").get<int>()<1 ||
      c.at("merge_concurrency").get<int>()<1 || c.at("disk_mb").get<int>()<8192 ||
      c.at("inspect_sources").get<int>()<1 || c.at("inspect_sources").get<int>()>10000 || c.at("inspect_concurrency").get<int>()<1)
    throw std::runtime_error("invalid campaign resource bounds");
  for (const char *key : {"plan_memory_mb","fill_memory_mb","merge_memory_mb","collection_memory_mb"})
    if (c.at(key).get<int>()<512 || c.at(key).get<int>()>65536) throw std::runtime_error("invalid campaign memory request");
  const fs::path work=c.at("work").get<std::string>();
  const std::string settings=safePath(base,cfg.at("Settings")), inputs=safePath(base,cfg.at("Inputs"));
  Json frozen=freezeInputs(inputs);
  c["resolved_settings"]=settingsSnapshot(settings,c.at("home"));
  c["publication_bindings"]=frozen.at("publication_bindings");
  c["setup"]=c.at("setup_script");
  for (const auto &a : cfg.value("SetupArguments",Json::array())) {
    const std::string token=a;
    if (!std::regex_match(token,std::regex("[A-Za-z0-9_./+@=-]+"))) throw std::runtime_error("invalid setup argument");
    c["setup"]=c.at("setup").get<std::string>()+" "+token;
  }
  c["setup_sha256"]=sha256File(c.at("setup_script"));
  c["binary"]=c.at("home").get<std::string>()+"/build/tree_to_hists"; c["binary_sha256"]=sha256File(c.at("binary"));
  c["entrypoint"]=c.at("home").get<std::string>()+"/tree_to_hists.sh"; c["entrypoint_sha256"]=sha256File(c.at("entrypoint"));
  const std::string timingFile=cfg.value("TimingFile",std::string("none"));
  c["timing_path"]=timingFile=="none" ? "none" : safePath(base,timingFile);
  c["timing"]=timing(c.at("timing_path"));
  if (!c.at("timing").is_null() && c.at("timing").at("reduction_settings_sha256")!=digest(c.at("resolved_settings")))
    throw std::runtime_error("timing receipt belongs to different reduction settings");
  c["timing_sha256"]=c.at("timing").is_null() ? Json() : Json(sha256File(c.at("timing_path")));
  c["settings_source"]=settings; c["settings_sha256"]=sha256File(settings); c["settings"]=settings;
  c["inputs"]=(work/"inputs.yaml").string();
  std::ostringstream frozenText; yaml(frozenText,frozen);
  if (fs::exists(work/"campaign.json")) {
    Json old=readJson((work/"campaign.json").string());
    c["inputs_sha256"]=old.at("inputs_sha256");
    c["inspection_index_sha256"]=old.at("inspection_index_sha256");inspectionParts(c);
    if (old!=c || sha256File(c.at("inputs"))!=c.at("inputs_sha256").get<std::string>())
      throw std::runtime_error("campaign differs; choose a new work directory");
    std::cout<<"[TreeToHists] existing frozen campaign: "<<work<<'\n'; return;
  }
  if (fs::exists(work) && !fs::is_empty(work)) throw std::runtime_error("campaign work directory is not empty");
  fs::create_directories(work/"condor"); atomic(work/"inputs.yaml",frozenText.str()); c["inputs_sha256"]=sha256File(c.at("inputs"));
  /* Small header jobs precede graph assembly. One result per group avoids
     both a serial scan of the campaign and one metadata file per source. */
  Json parts=Json::array(),files=Json::array();
  auto flushInspection=[&]() {
    if (files.empty()) return;
    const std::string id=std::to_string(parts.size());const fs::path path=work/"inspection/rows"/(id+".json");
    atomic(path,Json{{"inputs_sha256",c.at("inputs_sha256")},{"files",files}}.dump()+"\n");
    Json pin=filePin(path);pin["result"]=(work/"inspection/results"/(id+".json")).string();pin["count"]=files.size();
    parts.push_back(pin);files=Json::array();
  };
  for (const auto &sample : frozen.at("samples")) for (const auto &file : sample.at("files")) {
    files.push_back(file);if (files.size()==c.at("inspect_sources").get<size_t>()) flushInspection();
  }
  flushInspection();atomic(work/"inspection/index.json",parts.dump()+"\n");
  c["inspection_index_sha256"]=sha256File((work/"inspection/index.json").string());
  atomic(work/"campaign.json",c.dump(1)+"\n");
  const bool ready=!c.at("timing").is_null() && c.at("timing").value("environment","")=="SDCC";
  atomic(work/"admission.json",Json{{"status",ready ? "READY_FOR_SITE_PREFLIGHT" : "NEEDS_SDCC_TIMING"},
    {"campaign_sha256",sha256File((work/"campaign.json").string())},{"submission_authorized",false}}.dump(1)+"\n");
  fs::create_directories(fs::path(c.at("logs").get<std::string>())/"0");
  atomic(work/"condor/plan.sub",submit(c,"campaign-plan --work "+work.string(),"plan"));
  atomic(work/"condor/inspect.sub",submit(c,"campaign-inspect --work "+work.string()+" --part $(row)","inspect"));
  atomic(work/"condor/fill.sub",submit(c,"fill --work "+work.string()+" --job $(row)","fill"));
  atomic(work/"condor/fill_group.sub",submit(c,"fill-group --work "+work.string()+" --group $(row)","fill"));
  atomic(work/"condor/merge.sub",submit(c,"merge-group --work "+work.string()+" --node $(node)","merge"));
  atomic(work/"condor/collection.sub",submit(c,"collection --work "+work.string()+" --out "+c.at("output").get<std::string>(),"collection"));
  atomic(work/"condor/refuse.sub","universe = vanilla\nexecutable = /usr/bin/false\ngetenv = False\non_exit_hold = False\non_exit_remove = True\nqueue 1\n");
  atomic(work/"condor/reduce.dag","# Replaced atomically by the compute-node planner before this child starts.\nJOB UNPLANNED refuse.sub\n");
  std::ostringstream launch;
  for (size_t i=0;i<parts.size();++i) launch<<"JOB I"<<i<<" inspect.sub\nVARS I"<<i<<" row=\""<<i<<"\" shard=\"0\"\nCATEGORY I"<<i<<" INSPECT\n";
  launch<<"JOB PLAN plan.sub\nVARS PLAN row=\"0\" shard=\"0\"\nPARENT";
  for (size_t i=0;i<parts.size();++i) launch<<" I"<<i;
  launch<<" CHILD PLAN\nMAXJOBS INSPECT "<<c.at("inspect_concurrency")<<"\nSUBDAG EXTERNAL REDUCE reduce.dag\nPARENT PLAN CHILD REDUCE\n";
  atomic(work/"condor/campaign.dag",launch.str());
  std::cout<<"[TreeToHists] frozen campaign: "<<work<<"\n[TreeToHists] "<<(ready ? "timing bound; site preflight still required" : "preparation only; measured SDCC timing is required before submission")<<'\n';
}

void makeDag(const fs::path &work) {
  const Json c=readJson((work/"campaign.json").string()); bindings(c);
  if (c.at("timing").is_null()) throw std::runtime_error("measured timing is required before farm planning");
  long cap=c.at("max_events");
  for (const auto &[name,row] : c.at("timing").at("classes").items())
    cap=std::min(cap,static_cast<long>(std::floor(8*3600/(1.5*row.at("seconds").get<double>()/row.at("processed_entries").get<long long>()))));
  if (cap<1) throw std::runtime_error("one event exceeds the wall-time target");
  Json p=preparePlan(c.at("settings"),c.at("inputs"),c.at("home"),{},cap,work.string(),true,inspectionParts(c));
  if (p.at("selection")!=c.at("resolved_settings").at("selection")) throw std::runtime_error("planned selection differs from frozen settings");
  for (const auto &s : p.at("samples"))
    if (!c.at("timing").at("classes").contains(s.at("system").get<std::string>()+"/"+s.at("family").get<std::string>()))
      throw std::runtime_error("missing measured timing for a selected sample class");
  savePlan(p,work.string(),false); // Only range rows need individual GPFS descriptors.
  ReductionProgress graphProgress("campaign-graph",p.at("samples").size()+1);
  long long plannedSamples=0;
  std::ostringstream dag;
  /* The site dispatcher consumes this same graph, not a second implementation
     of grouping or merging. Descriptor hashes bind every exported node to the
     native plan; the ordinary Condor DAG remains the public launch path. */
  Json nodes=Json::array();
  long fillCount=0,groupId=0;
  std::vector<std::string> finals; long nodeId=0;
  std::map<std::string,Json> packageSamples;
  auto compactSample=[&](const Json &s) {
    Json small=workerHeader(s); small.erase("sources");
    small["_source_count"]=sourceCount(s); small["_package_binding_sha256"]=packageBinding(p,s);
    return small;
  };
  for (const auto &s : p.at("samples")) if (s.contains("_source_file")) {
    const auto &pin=s.at("_source_file");
    if (sha256File(pin.at("path"))!=pin.at("sha256").get<std::string>()) throw std::runtime_error("planned source inventory changed");
  }
  std::map<std::string,std::pair<std::string,std::string>> runFinals;
  std::map<std::string,std::map<size_t,std::vector<std::string>>> sourceJobs;
  for (size_t i=0;i<p.at("jobs").size();++i) {
    const auto &job=p.at("jobs")[i]; sourceJobs[job.at("sample")][job.at("source").get<size_t>()].push_back("F"+std::to_string(i));
  }
  auto node=[&](const Json &s,const std::vector<std::string> &parents,const std::vector<std::string> &inputs) {
    const std::string id="M"+std::to_string(nodeId++);
    const std::string out=(s.value("intermediate",false) ? work/"merge"/id : fs::path(c.at("output").get<std::string>())) .string();
    Json small={{"producer_code",p.at("producer_code")},{"selection",p.at("selection")},{"histograms",p.at("histograms")},
      {"root_version",p.at("root_version")},{"plan_sha256",p.at("plan_sha256")},{"jobs",Json::array()}};
    if (inputs.empty() && !s.contains("reused_package")) for (const auto &src : s.at("sources"))
      for (const auto &name : sourceJobs[s.at("key")][src.at("index").get<size_t>()])
        small["jobs"].push_back(p.at("jobs").at(std::stoul(name.substr(1))));
    Json compact=compactSample(s), inputSamples=Json::array();
    if (!inputs.empty()) for (const auto &parent : parents) inputSamples.push_back(packageSamples.at(parent));
    Json descriptor={{"sample",inputs.empty() && !s.contains("reused_package") ? s : compact},{"plan",small},{"inputs",inputs},
      {"input_samples",inputSamples},{"output_directory",out},{"plan_sha256",p.at("plan_sha256")}};
    packageSamples[id]=compact;
    descriptor["descriptor_sha256"]=digest(descriptor);
    const fs::path descriptorPath=work/"merge_rows"/(id+".json");
    atomic(descriptorPath,descriptor.dump()+"\n");
    dag<<"JOB "<<id<<" merge.sub\nVARS "<<id<<" row=\""<<id<<"\" shard=\"0\" node=\""<<id<<"\"\nCATEGORY "<<id<<" MERGE\n";
    if (!parents.empty()) {
      dag<<"PARENT"; for (const auto &parent : parents) dag<<' '<<parent;
      dag<<" CHILD "<<id<<'\n';
    }
    nodes.push_back({{"id",id},{"operation","merge-group"},{"sample",s.at("key")},
      {"sample_class",s.at("system").get<std::string>()+"/"+s.at("family").get<std::string>()},
      {"argv",Json::array({"merge-group","--work",work.string(),"--node",id})},
      {"parents",parents},{"descriptor",filePin(descriptorPath)},{"source_count",sourceCount(s)},
      {"progress",{{"row_id","merge"},{"total",sourceCount(s)+1},{"unit","work units"},{"path","environment"}}},
      {"output",(fs::path(out)/(s.at("key").get<std::string>()+".root")).string()}});
    return std::make_pair(id,(fs::path(out)/(s.at("key").get<std::string>()+".root")).string());
  };
  for (const auto &s : p.at("samples")) {
    if (s.contains("reused_package")) {
      const auto result=node(s,{},{}); finals.push_back(result.first);
      if (s.contains("run")) runFinals[s.at("key")]=result;
      graphProgress.advance(++plannedSamples);
      continue;
    }
    struct Part {std::vector<SourceSpan> sources; std::string node,path;}; std::vector<Part> level;
    const size_t fan=c.at("fan_in");
    if (s.contains("run_packages")) {
      for (const auto &key : s.at("run_packages")) {
        const auto result=runFinals.at(key.get<std::string>());
        const auto &perRun=sampleOf(p,key);
        const auto &pin=perRun.at("_source_file");
        level.push_back({{{pin.at("path"),0,sourceCount(perRun)}},result.first,result.second});
      }
      /* Even a one-run collection is a distinct combined analysis product. */
      if (level.size()==1) {
        auto result=node(s,{level[0].node},{level[0].path});
        level={{level[0].sources,result.first,result.second}};
      }
    }
    if (!s.contains("run_packages")) {
      Json batch=Json::array(); long long batchEntries=0,batchOffset=0;
      const std::string sourcePath=s.at("_source_file").at("path");
      auto flush=[&]() {
        if (batch.empty()) return;
        const std::string shard=std::to_string(groupId/1000), id="G"+std::to_string(groupId++);
        Json sub=workerHeader(s); sub["sources"]=batch; sub["intermediate"]=true;
        Json local=sub, jobs=Json::array();
        for (size_t i=0;i<batch.size();++i) {
          const size_t original=batch.at(i).at("index"); local["sources"][i]["index"]=i;
          for (const auto &name : sourceJobs[s.at("key")][original]) {
            Json job=p.at("jobs").at(std::stoul(name.substr(1))); job["source"]=i; jobs.push_back(job);
          }
        }
        Json small={{"schema",p.at("schema")},{"producer_code",p.at("producer_code")},{"selection",p.at("selection")},
          {"histograms",p.at("histograms")},{"root_version",p.at("root_version")},{"plan_sha256",p.at("plan_sha256")},
          {"samples",Json::array({local})},{"jobs",jobs}};
        Json descriptor={{"plan",small}}; descriptor["descriptor_sha256"]=digest(descriptor);
        const fs::path descriptorPath=work/"fill_groups"/shard/(id+".json");
        atomic(descriptorPath,descriptor.dump()+"\n");
        fs::create_directories(fs::path(c.at("logs").get<std::string>())/std::to_string(fillCount/1000));
        dag<<"JOB "<<id<<" fill_group.sub\nVARS "<<id<<" row=\""<<id<<"\" shard=\""<<fillCount/1000<<"\"\nCATEGORY "<<id<<" FILL\n";
        nodes.push_back({{"id",id},{"operation","fill-group"},{"sample",s.at("key")},
          {"sample_class",s.at("system").get<std::string>()+"/"+s.at("family").get<std::string>()},
          {"argv",Json::array({"fill-group","--work",work.string(),"--group",id})},{"parents",Json::array()},
          {"descriptor",filePin(descriptorPath)},{"source_count",batch.size()},{"physical_entries",batchEntries},
          {"progress",{{"row_id",id},{"total",batchEntries*fillPasses(s)},{"unit","event_passes"},
            {"path",(work/"chunks"/shard/id/"progress.json").string()}}},
          {"output",(work/"chunks"/shard/id/(s.at("key").get<std::string>()+".root")).string()}});
        packageSamples[id]=compactSample(sub);
        level.push_back({{{sourcePath,batchOffset,batch.size()}},id,(work/"chunks"/shard/id/(s.at("key").get<std::string>()+".root")).string()});
        ++fillCount; batch=Json::array(); batchEntries=0;
      };
      std::ifstream sourceStream(sourcePath); std::string sourceLine;
      if (!sourceStream) throw std::runtime_error("cannot read planned source inventory");
      while (true) {
        const long long offset=static_cast<long long>(sourceStream.tellg());
        if (!std::getline(sourceStream,sourceLine)) break;
        const Json src=Json::parse(sourceLine);
        const long long entries=src.at("metadata").at("entries");
        if (entries>cap) {
          flush(); // Keep every range of a large original in the same source-aware merge.
          std::vector<std::string> parents;
          for (const auto &name : sourceJobs[s.at("key")][src.at("index").get<size_t>()]) {
            const size_t i=std::stoul(name.substr(1)); savePlanRow(p,i,work.string(),src);
            fs::create_directories(fs::path(c.at("logs").get<std::string>())/std::to_string(i/1000));
            dag<<"JOB "<<name<<" fill.sub\nVARS "<<name<<" row=\""<<i<<"\" shard=\""<<i/1000<<"\"\nCATEGORY "<<name<<" FILL\n";
            const Json &job=p.at("jobs").at(i);
            const long long events=job.at("range_end").get<long long>()-job.at("range_begin").get<long long>();
            nodes.push_back({{"id",name},{"operation","fill"},{"sample",s.at("key")},
              {"sample_class",s.at("system").get<std::string>()+"/"+s.at("family").get<std::string>()},
              {"argv",Json::array({"fill","--work",work.string(),"--job",std::to_string(i)})},{"parents",Json::array()},
              {"descriptor",filePin(work/"rows"/std::to_string(i/1000)/(std::to_string(i)+".json"))},
              {"source_count",1},{"physical_entries",events},
              {"progress",{{"row_id",name},{"total",events*fillPasses(s)},{"unit","event_passes"},
                {"path",(fs::path(jobDirectory(work.string(),job))/"progress.json").string()}}},
              {"output",jobDirectory(work.string(),job)}});
            parents.push_back(name); ++fillCount;
          }
          Json sub=workerHeader(s); sub["sources"]=Json::array({src}); sub["intermediate"]=true;
          const auto result=node(sub,parents,{}); level.push_back({{{sourcePath,offset,1}},result.first,result.second});
        } else {
          if (!batch.empty() && (batchEntries+entries>cap || batch.size()>=c.at("max_sources").get<size_t>())) flush();
          if (batch.empty()) batchOffset=offset;
          batch.push_back(src); batchEntries+=entries;
        }
      }
      if (!sourceStream.eof()) throw std::runtime_error("planned source inventory read failed");
      flush();
      if (level.size()==1) {
        const auto result=node(s,{level[0].node},{level[0].path});
        level={{level[0].sources,result.first,result.second}};
      }
    }
    while (level.size()>1) {
      std::vector<Part> next;
      for (size_t begin=0;begin<level.size();begin+=fan) {
        std::vector<SourceSpan> spans;
        std::vector<std::string> parents,inputs;
        for (size_t i=begin;i<std::min(level.size(),begin+fan);++i) {
          spans.insert(spans.end(),level[i].sources.begin(),level[i].sources.end());
          parents.push_back(level[i].node); inputs.push_back(level[i].path);
        }
        Json sub=spanSample(p,s,spans,level.size()>fan);
        auto result=node(sub,parents,inputs); next.push_back({std::move(spans),result.first,result.second});
      }
      level=std::move(next);
    }
    finals.push_back(level.at(0).node);
    if (s.contains("run")) runFinals[s.at("key")]=std::make_pair(level.at(0).node,level.at(0).path);
    graphProgress.advance(++plannedSamples);
  }
  dag<<"JOB COLLECTION collection.sub\nVARS COLLECTION row=\"0\" shard=\"0\"\nPARENT";
  for (const auto &f : finals) dag<<' '<<f;
  dag<<" CHILD COLLECTION\nMAXJOBS FILL "<<c.at("fill_concurrency")<<"\nMAXJOBS MERGE "<<c.at("merge_concurrency")<<'\n';
  atomic(work/"condor/reduce.dag",dag.str());
  atomic(work/"collection_index.json",collectionIndex(p).dump()+"\n");
  nodes.push_back({{"id","COLLECTION"},{"operation","collection"},{"argv",Json::array({"collection","--work",work.string(),"--out",c.at("output")})},
    {"parents",finals},{"descriptor",filePin(work/"collection_index.json")},
    {"progress",{{"row_id","collection"},{"total",p.at("samples").size()+1},{"unit","work units"},{"path","environment"}}},
    {"output",(fs::path(c.at("output").get<std::string>())/"collection.yaml").string()}});
  Json graph={{"schema","TreeToHistsCampaignGraphV1"},{"campaign",filePin(work/"campaign.json")},
    {"plan",filePin(work/"plan.json")},{"plan_sha256",p.at("plan_sha256")},{"dag",filePin(work/"condor/reduce.dag")},
    {"producer_code_sha256",digest(p.at("producer_code"))},{"root_version",p.at("root_version")},
    {"fill_event_cap",cap},{"fill_source_cap",c.at("max_sources")},{"merge_fan_in",c.at("fan_in")},
    {"fill_concurrency",c.at("fill_concurrency")},{"merge_concurrency",c.at("merge_concurrency")},{"nodes",nodes}};
  graph["graph_sha256"]=digest(graph);
  atomic(work/"campaign_graph.json",graph.dump()+"\n");
  graphProgress.finish();
  std::cout<<"[TreeToHists] "<<fillCount<<" bounded fill workers ("<<groupId<<" source groups), "<<nodeId
    <<" merge nodes; "<<p.at("jobs").size()<<" exact source ranges; collection is the terminal stage\n";
}
void collection(const fs::path &work,const fs::path &output) {
  Json p;
  if (fs::exists(work/"collection_index.json")) p=readJson((work/"collection_index.json").string());
  else {
    /* Interactive/local reductions have no campaign graph. Construct the
       same small index from their plan; farm campaigns always write it once. */
    p=collectionIndex(readJson((work/"plan.json").string()));
  }
  const std::string expected=p.at("descriptor_sha256"); p.erase("descriptor_sha256");
  if (p.at("schema")!="TreeToHistsCollectionIndexV1" || digest(p)!=expected ||
      p.at("plan_sha256")!=readJson((work/"binding.json").string()).at("plan_sha256"))
    throw std::runtime_error("collection index differs from its immutable plan");
  checkReductionCode(p);
  Json files=Json::array(),samples=Json::array();
  ReductionProgress progress("collection",p.at("samples").size()+1); long long completed=0;
  for (const auto &s : p.at("samples")) {
    const std::string key=s.at("key"), path=(output/(key+".root")).string();
    fs::path side=path; side.replace_extension(".json"); const Json receipt=readJson(side.string());
    const auto &validation=receipt.at("validation");
    if (receipt.at("package")!=fs::path(path).filename().string() || receipt.at("sha256")!=sha256File(path) ||
        receipt.at("selection_sha256")!=p.at("selection").at("selection_sha256") ||
        receipt.at("sources")!=s.at("source_count") || receipt.at("package_binding_sha256")!=s.at("package_binding_sha256") ||
        validation.at("schema")!="PhotonJetPackageValidationV1" || validation.at("status")!="PASS" ||
        validation.at("producer_code_sha256")!=digest(p.at("producer_code")) || validation.at("root_version")!=p.at("root_version"))
      throw std::runtime_error("package is incomplete, changed or belongs to another collection: " + path);
    Json row={{"key",key},{"path",key+".root"},{"system",s.at("system")},{"name",s.at("name")},{"family",s.at("family")},
      {"sha256",receipt.at("sha256")},{"coverage",s.at("coverage")},{"source_count",receipt.at("sources")},
      {"event_records",receipt.at("event_records")},{"histograms",receipt.at("histograms")},
      {"availability",receipt.at("availability")},{"validated_at",utc()}};
    row["package_plan_sha256"]=receipt.at("plan_sha256");
    if (receipt.contains("package_binding_sha256")) row["package_binding_sha256"]=receipt.at("package_binding_sha256");
    if (s.contains("run")) row["run"]=s.at("run");
    if (s.contains("run_packages")) row["run_packages"]=s.at("run_packages");
    for (const char *field : {"si_di_role","simulation_periods"}) if (s.contains(field)) row[field]=s.at(field);
    samples.push_back(row);
    files.push_back({{"path",key+".root"},{"sha256",receipt.at("sha256")}});
    files.push_back({{"path",key+".json"},{"sha256",sha256File(side.string())}});
    progress.advance(++completed);
  }
  const Json configuration={{"selection.yaml",p.at("selection")},{"histograms.yaml",p.at("histograms")}};
  for (const auto &[name,value] : configuration.items()) {
    writeYaml(output/name,value); files.push_back({{"path",name},{"sha256",sha256File((output/name).string())}});
  }
  Json manifest={{"schema","PhotonJetHistogramCollectionV1"},{"status","VALIDATED"},{"created_at",utc()},
    {"input_snapshot_sha256",p.at("plan_sha256")},{"selection_sha256",p.at("selection").at("selection_sha256")},
    {"samples",samples},{"files",files}};
  const fs::path path=output/"collection.yaml";
  if (fs::exists(path)) {
    const Json old=photonjet_config::readYaml(path.string());
    manifest["created_at"]=old.at("created_at");
    for (size_t i=0;i<samples.size();++i) manifest["samples"][i]["validated_at"]=old.at("samples").at(i).at("validated_at");
    if (old!=manifest) throw std::runtime_error("existing collection differs; publication is never replaced");
  } else writeYaml(path,manifest);
  progress.finish();
  std::cout<<"[TreeToHists] validated collection: "<<path<<'\n';
}
}  // namespace
int campaign(const std::string &command,const std::vector<std::string> &args) {
  if (command=="campaign") {prepare(args); return 0;}
  if (command=="campaign-submit-check") {
    const fs::path config=fs::absolute(arg(args,"--config"));
    const Json cfg=photonjet_config::readYaml(config.string()).at("Campaign");
    const fs::path work=safePath(config.parent_path(),cfg.at("WorkDirectory"));
    const Json c=readJson((work/"campaign.json").string()); bindings(c);
    const Json admission=readJson((work/"admission.json").string());
    if (admission.at("status")!="READY_FOR_SITE_PREFLIGHT" || admission.at("campaign_sha256")!=sha256File((work/"campaign.json").string()))
      throw std::runtime_error("campaign is preparation-only; measured SDCC timing and site preflight are required");
    const char *account=std::getenv("USER");
    if (!account || !*account) throw std::runtime_error("USER is required for the SDCC namespace checks");
    const std::string user=account;
    for (const char *key : {"work","output","logs"}) {
      const std::string path=c.at(key);
      if (path.rfind("/sphenix/u/"+user+"/",0)!=0 && path.rfind("/sphenix/user/"+user+"/",0)!=0)
        throw std::runtime_error("campaign directories must be in the normal GPFS user namespace");
    }
    const fs::path events=c.at("events").get<std::string>();
    if (events.string().rfind("/tmp/"+user+"/",0)!=0) throw std::runtime_error("Condor event log must be under /tmp/USER/");
    fs::create_directories(events.parent_path());
    std::cout<<work.string()<<'\n'; return 0;
  }
  const fs::path work=fs::absolute(arg(args,"--work"));
  if (command=="campaign-inspect") {
    const Json c=readJson((work/"campaign.json").string());bindings(c,false);
    const std::string id=arg(args,"--part");
    if (!std::regex_match(id,std::regex("0|[1-9][0-9]*"))) throw std::runtime_error("invalid inspection part");
    const Json parts=inspectionParts(c);const auto &pin=parts.at(std::stoul(id));
    if (!std::getenv("TREETOHISTS_PROGRESS_PATH")) {
      const fs::path progress=work/"inspection/progress"/(id+".json");fs::create_directories(progress.parent_path());
      setenv("TREETOHISTS_PROGRESS_PATH",progress.c_str(),1);
    }
    if (sha256File(pin.at("path"))!=pin.at("sha256").get<std::string>()) throw std::runtime_error("inspection descriptor changed");
    const Json descriptor=readJson(pin.at("path")),snapshot=inspectionSnapshot(c);
    if (descriptor.at("inputs_sha256")!=c.at("inputs_sha256") || descriptor.at("files").size()!=pin.at("count"))
      throw std::runtime_error("inspection descriptor belongs to another input snapshot");
    const fs::path out=pin.at("result").get<std::string>();
    if (fs::exists(out)) {
      Json saved=readJson(out.string());const std::string sha=saved.at("sha256");saved.erase("sha256");
      if (digest(saved)!=sha || saved.at("snapshot")!=snapshot || saved.at("descriptor_sha256")!=pin.at("sha256") ||
          saved.at("rows").size()!=descriptor.at("files").size()) throw std::runtime_error("existing inspection result differs");
      ReductionProgress progress("inspect",saved.at("rows").size());progress.finish();return 0;
    }
    Json result={{"schema","TreeToHistsSourceInspectionV1"},{"descriptor_sha256",pin.at("sha256")},{"snapshot",snapshot},
      {"rows",inspectSources(descriptor,snapshot)}};
    result["sha256"]=digest(result);atomic(out,result.dump()+"\n");return 0;
  }
  if (command=="campaign-plan") {makeDag(work); return 0;}
  if (command=="collection") {collection(work,arg(args,"--out")); return 0;}
  if (command=="merge-group") {
    const Json c=readJson((work/"campaign.json").string()); bindings(c,false);
    const std::string id=arg(args,"--node");
    if (!std::regex_match(id,std::regex("M[0-9]+"))) throw std::runtime_error("invalid merge node");
    Json node=readJson((work/"merge_rows"/(id+".json")).string()); const std::string sha=node.at("descriptor_sha256");
    node.erase("descriptor_sha256");
    if (digest(node)!=sha || node.at("plan_sha256")!=readJson((work/"binding.json").string()).at("plan_sha256")) throw std::runtime_error("merge node differs from its immutable plan");
    const Json &p=node.at("plan");
    const Json &s=node.at("sample");
    if (node.at("inputs").empty()) merge(p,s,work.string(),node.at("output_directory"));
    else combinePackages(p,s,node.at("inputs").get<std::vector<std::string>>(),
      (fs::path(node.at("output_directory").get<std::string>())/(s.at("key").get<std::string>()+".root")).string(),node.value("input_samples",Json::array()));
    return 0;
  }
  throw std::runtime_error("unknown campaign command");
}
}  // namespace treetohists
