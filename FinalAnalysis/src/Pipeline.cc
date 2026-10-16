/* Pipeline: `final_analysis.sh all` and `final_analysis.sh compact`. */

/* `all` runs the package's own commands in order (QA, the nominal measurement, every variation, the systematic
   uncertainties, the figures and the tables) into one new run folder, and PIPELINE.json binds what went in and what
   came out. No physics is computed here: every number comes from a command this file starts.

   In reading order:
     1. Settings files and small helpers
     2. Starting the package's commands
     3. The histogram collection: checking a delivered one, choosing its samples, making its compact twin
     4. The pipeline: its plan, each variation's settings, the run and its record */
#include "FinalAnalysis.h"
#include "RenderMetadata.h"
#include "MeasurementReplay.h"
#include "NamedResponse.h"
#include "PhysicsGeneration.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cerrno>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include "TEnv.h"
#include "../../TreeToHists/src/Configuration.h"
#include "../../TreeToHists/src/CollectionContract.h"
#include "TFile.h"
#include "THashList.h"
#include "TObjString.h"
#include "TROOT.h"

namespace fa {
namespace {

namespace fs = std::filesystem;
using Json = nlohmann::json;
using Config = std::map<std::string, std::string>;  /* a settings file: key -> value */

const char *const kNameCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";

/* ---------------------------------------------------------------------------------------------------------------
   1. Settings files and small helpers
   --------------------------------------------------------------------------------------------------------------- */

std::vector<std::string> words(const std::string &s) {
  std::istringstream in(s);
  std::vector<std::string> out;
  for (std::string v; in >> v;) out.push_back(v);
  return out;
}

/* A settings file as key -> value. A key that appears twice is an error (TEnv would keep the first and only warn). */
Config readSettings(const fs::path &path) {
  if (!photonjet_config::isYaml(path.string())) {
    std::ifstream source(path);
    std::set<std::string> keys;
    for (std::string line; std::getline(source, line);) {
      line = line.substr(0, line.find('#'));
      const auto colon = line.find(':');
      if (colon == std::string::npos) continue;
      const auto names = words(line.substr(0, colon));
      if (names.size() == 1 && !keys.insert(names[0]).second)
        throw std::runtime_error("duplicate configuration key " + names[0] + " in " + path.string());
    }
  }
  TEnv env;
  photonjet_config::readSettings(path.string(), env);
  Config out;
  TIter next(env.GetTable());
  while (auto *r = static_cast<TEnvRec *>(next())) out[r->GetName()] = r->GetValue();
  return out;
}

void writeSettings(const fs::path &path, const Config &settings) {
  std::ofstream out(path);
  for (const auto &[key, value] : settings) out << key << ": " << Json(value).dump() << '\n';
  out.close();
  if (!out) throw std::runtime_error("cannot save " + path.string());
}

/* Written beside its final name and moved into place, so a reader never sees half a file. */
void writeJson(const fs::path &path, const Json &value) {
  std::ofstream out(path.string() + ".partial");
  out << value.dump(2) << '\n';
  out.close();
  if (!out) throw std::runtime_error("cannot save " + path.string());
  fs::rename(path.string() + ".partial", path);
}

std::string required(const Config &c, const std::string &key) {
  const auto it = c.find(key);
  if (it == c.end() || it->second.empty()) throw std::runtime_error("missing " + key);
  return it->second;
}

/* A path of a settings file: relative paths are relative to that file. */
fs::path resolve(const fs::path &base, const std::string &v) {
  return fs::absolute(fs::path(v).is_absolute() ? fs::path(v) : base / v).lexically_normal();
}

std::string hash(const fs::path &p) {
  if (!fs::is_regular_file(p)) throw std::runtime_error("missing pinned file " + p.string());
  return photonjet_render::sha256(p.string(), true);
}

/* Applies the settings of `run` that start with `prefix` (<variation>.Analysis. or <variation>.Reduction.) to `c`. */
void overlay(Config &c, const Config &run, const std::string &prefix) {
  for (const auto &[k, v] : run)
    if (k.rfind(prefix, 0) == 0) {
      const auto key = k.substr(prefix.size());
      if (!c.count(key)) throw std::runtime_error("override names an absent setting: " + k);
      c[key] = v;
    }
}

Config reductionSettings(const fs::path &path) {
  Config c = readSettings(path);
  for (const char *k : {"Selection.PhotonIDPackage", "Selection.AuAuIsolationPackage"})
    if (c.count(k) && !c[k].empty()) c[k] = resolve(path.parent_path(), c[k]).string();
  return c;
}

/* ---------------------------------------------------------------------------------------------------------------
   2. Starting the package's commands
   --------------------------------------------------------------------------------------------------------------- */

/* One command: the name it is reported under, its exact arguments (no shell reads them) and the log of its output. */
struct Command {
  std::string name;
  std::vector<std::string> arguments;
  fs::path log;
};

/* Every command runs in its own process group, so an interrupt of the pipeline reaches everything a command started. */
volatile sig_atomic_t interrupted = 0;
volatile sig_atomic_t groups[64] = {0};  /* the process groups running now */
void terminate(int signal) {
  interrupted = signal;
  for (int i = 0; i < 64; ++i)
    if (groups[i] > 0) kill(-groups[i], signal);
}

/* Runs the commands, at most `limit` at once, and returns their exit codes in order (-1: never started, the pipeline
   was interrupted first). */
std::vector<int> run(const std::vector<Command> &commands, size_t limit) {
  limit = std::max<size_t>(1, std::min<size_t>(limit, 64));
  std::vector<int> codes(commands.size(), -1);
  std::map<pid_t, std::pair<size_t, int>> running;  /* process -> (its command, its place in `groups`) */
  size_t next = 0;
  while (next < commands.size() || !running.empty()) {
    while (next < commands.size() && running.size() < limit && !interrupted) {
      const Command &command = commands[next];
      const int fd = open(command.log.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
      if (fd < 0) throw std::runtime_error("cannot create log " + command.log.string());
      const pid_t pid = fork();
      if (pid == 0) {
        setpgid(0, 0);
        dup2(fd, STDOUT_FILENO);
        dup2(fd, STDERR_FILENO);
        close(fd);
        std::vector<char *> argv;
        for (const auto &a : command.arguments) argv.push_back(const_cast<char *>(a.c_str()));
        argv.push_back(nullptr);
        execvp(argv[0], argv.data());
        _exit(127);
      }
      close(fd);
      if (pid < 0) throw std::runtime_error("cannot start " + command.name);
      setpgid(pid, pid);
      int place = 0;
      while (place < 64 && groups[place] > 0) ++place;
      groups[place] = pid;  /* at most 64 run at once, so a place is free */
      if (interrupted) kill(-pid, interrupted);
      running[pid] = {next++, place};
    }
    if (running.empty()) break;  /* interrupted before every command started */
    int status = 0;
    const pid_t pid = waitpid(-1, &status, 0);
    if (pid < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error("cannot wait for a command");
    }
    const auto found = running.find(pid);
    if (found == running.end()) continue;
    codes[found->second.first] = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    groups[found->second.second] = 0;
    running.erase(found);
  }
  return codes;
}

/* ---------------------------------------------------------------------------------------------------------------
   3. The histogram collection
   --------------------------------------------------------------------------------------------------------------- */

/* A transferred collection uses paths relative to its manifest. Only its declared payload is verified, including
   package identity and recorded coverage; no eventTrees or farm paths are needed. */
struct HistogramCollection {
  fs::path directory;
  Json manifest;                         /* collection.yaml; null for a plain directory of packages */
  std::vector<fs::path> files;
  std::map<std::string, std::string> pins;  /* file -> sha256, for the run record */
  Json otherSelections = Json::array();  /* samples analysed under the collection's selection but filled under another */
};

HistogramCollection histogramCollection(const fs::path &input) {
  if (fs::is_directory(input)) return {input, Json(), {}, {}, Json::array()};  /* a plain directory of packages */
  if (!photonjet_config::isYaml(input.string()))
    throw std::runtime_error("histogram collection must be a directory or collection.yaml");
  const Json c = photonjet_config::readYaml(input.string());
  if (c.value("schema", "") != "PhotonJetHistogramCollectionV1" || c.value("status", "") != "VALIDATED" ||
      c.value("samples", Json::array()).empty() || c.value("files", Json::array()).empty())
    throw std::runtime_error("histogram collection is incomplete or has not passed package validation");
  HistogramCollection out{fs::canonical(input.parent_path()), c, {input}, {{input.string(), hash(input)}}, Json::array()};

  /* a declared file, which must lie inside the collection's directory */
  auto payload = [&](const std::string &name) {
    const fs::path relative(name);
    if (name.empty() || relative.is_absolute() || relative.lexically_normal() != relative)
      throw std::runtime_error("unsafe collection path " + name);
    for (const auto &part : relative)
      if (part == ".." || part == ".") throw std::runtime_error("unsafe collection path " + name);
    const fs::path path = fs::canonical(out.directory / relative);
    const fs::path inside = path.lexically_relative(out.directory);
    if (inside.empty() || *inside.begin() == ".." || !fs::is_regular_file(path))
      throw std::runtime_error("collection file escapes its directory: " + name);
    return path;
  };

  /* every declared file has its declared hash */
  std::map<std::string, std::string> hashes;
  for (const auto &file : c.at("files")) {
    const std::string name = file.at("path"), expected = file.at("sha256");
    if (!hashes.emplace(name, expected).second) throw std::runtime_error("duplicate collection file " + name);
    const auto path = payload(name);
    const std::string actual = hash(path);
    if (actual != expected) {
      /* A compact twin (final_analysis.sh compact): the delivered package without its per-source records, bound to
         the delivered hash by <package>.compact.json beside it. Nothing else may differ from the manifest. */
      fs::path binding = path;
      binding.replace_extension(".compact.json");
      std::ifstream twin(binding);
      const Json record = twin ? Json::parse(twin, nullptr, false) : Json();
      if (!record.is_object() || record.value("schema", "") != "PhotonJetCompactPackageV1" ||
          record.value("source_sha256", "") != expected || record.value("compact_sha256", "") != actual)
        throw std::runtime_error("transferred collection file has changed: " + name);
      out.pins[binding.string()] = hash(binding);
    }
    out.files.push_back(path);
    out.pins[path.string()] = actual;
  }
  if (!hashes.count("selection.yaml") || !hashes.count("histograms.yaml"))
    throw std::runtime_error("collection configuration files are missing");
  const Json selection = photonjet_config::readYaml(payload("selection.yaml").string());
  const Json histograms = photonjet_config::readYaml(payload("histograms.yaml").string());
  const Json responseContract = c.value("simulation_response", Json());
  photonjet_collection::requireContract(responseContract, histograms);
  photonjet_collection::requireSelection(responseContract, selection);
  const bool replay=histograms.contains("measurement_replay");
  if(replay!=c.contains("measurement_replay"))throw std::runtime_error("collection replay declaration differs");
  if(replay) {
    const auto& binding=c.at("measurement_replay");
    const bool centrality=histograms.at("measurement_replay").contains("centrality_weight");
    if(centrality!=binding.contains("system_subset") || (centrality && binding.at("system_subset")!="AuAu"))
      throw std::runtime_error("centrality collection system subset differs");
    if(binding.at("schema")!="PhotonJetMeasurementReplayCollectionV1" ||
       binding.at("recipe_sha256")!=photonjet_render::sha256(histograms.at("measurement_replay").dump()) ||
       binding.at("admission")!=histograms.at("measurement_replay").at("admission") ||
       binding.at("nominal_diagnostics")!="parent_collection" ||
       binding.at("parent_packages").size()!=c.at("samples").size() ||
       c.contains("compact_of") || c.contains("component_collections"))
      throw std::runtime_error("invalid measurement replay collection binding");
  }
  if (selection.at("selection_sha256") != c.at("selection_sha256"))
    throw std::runtime_error("collection selection binding differs");
  std::set<std::string> allowed = {"selection.yaml", "histograms.yaml"};

  /* Component manifests preserve each production plan and the retained DATA run lineage without requiring a second
     copy of every run ROOT file locally. */
  std::map<std::string, Json> componentSamples, runLineage;
  std::map<std::string, std::string> componentSelection;  /* each final sample's own production selection */
  std::set<std::string> componentCompactPackages;  /* hash-bound native DATA components retain compact provenance */
  if (c.contains("component_collections")) {
    if (!c.at("component_collections").is_array() || c.at("component_collections").empty())
      throw std::runtime_error("empty component collection declaration");
    for (const auto &pin : c.at("component_collections")) {
      const std::string name = pin.at("path");
      if (!hashes.count(name) || hashes.at(name) != pin.at("sha256").get<std::string>() || !allowed.insert(name).second)
        throw std::runtime_error("component manifest is missing or duplicated");
      const Json component = photonjet_config::readYaml(payload(name).string());
      if (component.value("schema", "") != "PhotonJetHistogramCollectionV1" || component.value("status", "") != "VALIDATED" ||
          component.contains("component_collections") || !component.at("selection_sha256").is_string())
        throw std::runtime_error("invalid component collection");
      std::set<std::string> compactPaths;
      if (component.contains("compact_of")) {
        const auto &listed = component.at("compact_of").at("compact_packages");
        if (!listed.is_array()) throw std::runtime_error("invalid component compact package declaration");
        for (const auto &path : listed)
          if (!path.is_string() || !compactPaths.insert(path.get<std::string>()).second)
            throw std::runtime_error("invalid or duplicate component compact package");
      }
      for (const auto &row : component.at("samples")) {
        const std::string key = row.at("system").get<std::string>() + "/" + row.at("name").get<std::string>();
        if (compactPaths.erase(row.at("path").get<std::string>())) {
          if (row.at("family") != "data" || row.contains("run"))
            throw std::runtime_error("component compaction is restricted to combined DATA packages");
          componentCompactPackages.insert(row.at("path").get<std::string>());
        }
        auto &index = row.contains("run") ? runLineage : componentSamples;
        if (row.at("key") != key || !index.emplace(key, row).second)
          throw std::runtime_error("duplicate component sample or run lineage");
        if (!row.contains("run")) componentSelection[key] = component.at("selection_sha256");
      }
      if (!compactPaths.empty()) throw std::runtime_error("component compact declaration names an undeclared package");
    }
  }

  /* every sample: its package and sidecar are declared, and the package's own metadata says what the manifest says */
  std::set<std::string> samples;
  for (const auto &sample : c.at("samples")) {
    const std::string system = sample.at("system"), name = sample.at("name"), family = sample.at("family");
    const std::string key = system + "/" + name, root = key + ".root", sidecar = key + ".json";
    const bool dataRun = family == "data" && sample.contains("run") && sample.at("run").is_number_integer() &&
                         sample.at("run").get<int>() > 0 && name == "run_" + std::to_string(sample.at("run").get<int>());
    if ((system != "PP" && system != "AuAu") || name.empty() || name.find_first_not_of(kNameCharacters) != std::string::npos ||
        (family != "data" && family != "photonjet" && family != "inclusivejet") ||
        ((family == "data") != (name == "data" || dataRun)) || !samples.insert(key).second)
      throw std::runtime_error("invalid or duplicate sample in histogram collection: " + key);
    if (sample.at("path") != root || !hashes.count(root) || !hashes.count(sidecar) ||
        hashes.at(root) != sample.at("sha256").get<std::string>())
      throw std::runtime_error("collection is missing the ROOT package or sidecar for " + key);
    allowed.insert(root);
    allowed.insert(sidecar);
    if (c.contains("component_collections") && (!componentSamples.count(key) || componentSamples.at(key) != sample))
      throw std::runtime_error("final package differs from its native component collection: " + key);
    TFile file(payload(root).c_str(), "READ");
    auto *saved = dynamic_cast<TObjString *>(file.Get("meta/package"));
    if (file.IsZombie() || !saved) throw std::runtime_error("cannot read package metadata for " + key);
    const Json meta = Json::parse(saved->GetString().Data());
    const Json recipe=measurementReplayRecipe(meta);
    if(replay==recipe.is_null())throw std::runtime_error("collection mixes nominal and replay packages");
    if(replay && (recipe!=histograms.at("measurement_replay") ||
        meta.at("replay_provenance").at("parent_package_sha256")!=c.at("measurement_replay").at("parent_packages").at(key)))
      throw std::runtime_error("collection replay parent or recipe differs for "+key);
    /* A component keeps the selection it was produced with. Only a p+p DATA trigger name that the sample never
       evaluates may differ from the collection's selection (photonjet_render::selectionsAgree). */
    const Json &filled = meta.at("selection");
    if (!photonjet_render::selectionsAgree(filled, selection, system == "PP" && family == "data") ||
        (componentSelection.count(key) && componentSelection.at(key) != filled.at("selection_sha256").get<std::string>()))
      throw std::runtime_error("package selection differs from its collection for " + key);
    if (filled != selection) out.otherSelections.push_back({{"key", key}, {"selection_sha256", filled.at("selection_sha256")}});
    if (meta.value("schema", "") != "PhotonJetSampleV1" || meta.at("system") != system || meta.at("name") != name ||
        meta.at("family") != family || meta.at("coverage") != sample.at("coverage") ||
        meta.at("plan_sha256") != sample.value("package_plan_sha256", c.at("input_snapshot_sha256")) ||
        !photonjet_collection::histogramsAgree(meta.at("histograms"), histograms, responseContract, family, system) ||
        meta.at("package_kind") != "FINAL" ||
        meta.at("source_count") != sample.at("source_count") || meta.at("event_records") != sample.at("event_records") ||
        meta.at("histograms_written") != sample.at("histograms"))
      throw std::runtime_error("package identity, coverage or input snapshot differs for " + key);
    if (!c.contains("component_collections") && meta.at("plan_sha256") != c.at("input_snapshot_sha256") &&
        (!sample.contains("package_binding_sha256") || sample.contains("run_packages") ||
         meta.value("package_binding_sha256", Json()) != sample.at("package_binding_sha256")))
      throw std::runtime_error("reused package has no matching input/settings binding for " + key);
    /* A compact collection written as a collection of its own (compact_of at its top, made on the farm by
       redeclare_compact.py so that the transfer can carry it): its DATA packages are listed with the hash of the
       compact file, carry their own meta/compact record and come without their run packages. */
    bool compactPackage = componentCompactPackages.count(root) != 0;
    if (c.contains("compact_of"))
      for (const auto &listed : c.at("compact_of").at("compact_packages")) compactPackage = compactPackage || listed == root;
    if (compactPackage && !file.Get("meta/compact")) throw std::runtime_error("declared compact package has no compact record: " + key);
    if (meta.value("run", Json()) != sample.value("run", Json()) ||
        (!compactPackage && meta.value("run_packages", Json()) != sample.value("run_packages", Json())))
      throw std::runtime_error("collection DATA run binding differs for " + key);
    if (system == "PP" && family != "data" &&
        (meta.at("si_di_role") != sample.at("si_di_role") ||
         meta.at("selection").at("systems").at("pp").at("simulation_periods").at("mode") != sample.at("simulation_periods")))
      throw std::runtime_error("p+p period or SI/DI binding differs for " + key);
    std::ifstream sidecarFile(payload(sidecar));
    const Json record = Json::parse(sidecarFile);
    if (record.at("sha256") != hashes.at(root) || record.at("plan_sha256") != meta.at("plan_sha256") ||
        record.value("package_binding_sha256", Json()) != meta.value("package_binding_sha256", Json()) ||
        record.at("selection_sha256") != filled.at("selection_sha256"))
      throw std::runtime_error("package sidecar does not bind " + key);
  }
  const auto variedSystem=photonjet_collection::dataScaleSystem(responseContract);
  if (photonjet_collection::recipeV2(responseContract)) {
    bool haveData=false,haveSimulation=false;
    for (const auto &sample:c.at("samples")) {
      haveData=haveData || sample.at("family")=="data";
      haveSimulation=haveSimulation || sample.at("family")=="photonjet" || sample.at("family")=="inclusivejet";
    }
    if (!haveData || !haveSimulation) throw std::runtime_error("composed collection requires DATA and simulation");
  }
  if (!variedSystem.empty() && std::none_of(c.at("samples").begin(),c.at("samples").end(),[&](const Json &sample) {
        return sample.at("family")=="data" && sample.at("system")== (variedSystem=="pp" ? "PP" : "AuAu");
      })) throw std::runtime_error("DATA scale recipe has no DATA sample of its declared system");
  if (allowed.size() != hashes.size()) throw std::runtime_error("collection declares an unrelated payload");
  if (c.contains("component_collections") && samples.size() != componentSamples.size())
    throw std::runtime_error("portable collection omitted a final component sample");

  /* a combined DATA package accounts for exactly its run packages */
  for (const auto &sample : c.at("samples")) {
    if (!sample.contains("run_packages")) continue;
    std::set<std::string> runs;
    long long sourceCount = 0, eventCount = 0;
    if (sample.at("family") != "data" || sample.at("name") != "data" || sample.at("run_packages").empty())
      throw std::runtime_error("invalid combined DATA lineage");
    for (const auto &key : sample.at("run_packages")) {
      const auto found = std::find_if(c.at("samples").begin(), c.at("samples").end(), [&](const Json &s) {
        return s.at("system").get<std::string>() + "/" + s.at("name").get<std::string>() == key.get<std::string>();
      });
      const Json *run = found != c.at("samples").end() ? &*found
                        : runLineage.count(key.get<std::string>()) ? &runLineage.at(key.get<std::string>()) : nullptr;
      if (!run || !run->contains("run") || run->at("system") != sample.at("system") || !runs.insert(key).second)
        throw std::runtime_error("combined DATA is missing a declared run package");
      sourceCount += run->at("source_count").get<long long>();
      eventCount += run->at("event_records").get<long long>();
    }
    if (sourceCount != sample.at("source_count").get<long long>() || eventCount != sample.at("event_records").get<long long>())
      throw std::runtime_error("combined DATA source/event accounting differs from its run packages");
  }
  return out;
}

/* A varied collection must be derived from this exact nominal collection,
   with every sample and original package accounted for. Diagnostics still
   refer to that nominal input; only measurement-consumed inputs are varied. */
void bindMeasurementReplay(const HistogramCollection& varied,const HistogramCollection& nominal) {
  if(!varied.manifest.contains("measurement_replay"))return;
  if(nominal.manifest.is_null() || nominal.manifest.contains("measurement_replay"))
    throw std::runtime_error("measurement replay requires an original nominal collection");
  const auto& b=varied.manifest.at("measurement_replay");
  if(b.at("parent_collection_sha256")!=hash(nominal.directory/"collection.yaml"))
    throw std::runtime_error("replay collection belongs to another nominal collection");
  Json parents=Json::object();
  const auto recipe=photonjet_config::readYaml((varied.directory/"histograms.yaml").string()).at("measurement_replay");
  const bool centrality=recipe.contains("centrality_weight");
  for(const auto& sample:nominal.manifest.at("samples"))
    if(!centrality || sample.at("system")=="AuAu")parents[sample.at("key").get<std::string>()]=sample.at("sha256");
  if(b.at("parent_packages")!=parents)
    throw std::runtime_error("replay collection source roster differs from nominal");
  if(centrality) {
    Json variedKeys=Json::object();
    for(const auto& sample:varied.manifest.at("samples"))variedKeys[sample.at("key").get<std::string>()]=sample.at("family");
    const auto& c=recipe.at("centrality_weight");
    if(c.at("collection_subset")!="all_auau_parent_samples" || c.at("parent_packages")!=parents ||
       c.at("parent_collection").at("sha256")!=b.at("parent_collection_sha256") ||
       c.at("sample_families")!=variedKeys || variedKeys.size()!=parents.size())
      throw std::runtime_error("centrality replay collection has an incomplete or mixed AuAu roster");
  }
  if(recipe.contains("selection_change")) {
    const auto original=photonjet_config::readYaml((nominal.directory/"selection.yaml").string());
    const auto effective=photonjet_config::readYaml((varied.directory/"selection.yaml").string());
    photonjet_selection_replay::validate(recipe,original,effective,
      [](const Json& value){return photonjet_render::sha256(value.dump(-1,' ',true));});
  } else if(varied.manifest.at("selection_sha256")!=nominal.manifest.at("selection_sha256"))
    throw std::runtime_error("replay collection has an undeclared selection change");
}

/* Pipeline.UseCollectionSamples: the response and the inclusive-jet QA use exactly the samples of the verified
   collection. Analysis windows and toy settings stay in the analyst's settings. Missing DATA or photon simulation is
   an error, and a p+p sample of the two-period sum needs its single- or double-interaction partner. */
void useCollectionSamples(Config &analysis, const Json &collection) {
  if (collection.is_null()) throw std::runtime_error("Pipeline.UseCollectionSamples requires collection.yaml");
  if(collection.contains("measurement_replay") && collection.at("measurement_replay").contains("system_subset")) {
    if(collection.at("measurement_replay").at("system_subset")!="AuAu")
      throw std::runtime_error("unqualified collection system subset");
    std::string selected;std::set<std::string> kept,removed;
    for(const auto& system:words(required(analysis,"Systems"))) {
      if(required(analysis,system+".Package")=="AuAu") {
        selected+=(selected.empty()?"":" ")+system;kept.insert(system);
      } else removed.insert(system);
    }
    if(selected.empty())throw std::runtime_error("centrality alternative needs at least one configured AuAu result");
    analysis["Systems"]=selected;
    for(auto it=analysis.begin();it!=analysis.end();) {
      const auto dot=it->first.find('.');
      if(removed.count(it->first.substr(0,dot)) || it->first=="QA.PP.Inclusive")it=analysis.erase(it);
      else ++it;
    }
    if(analysis.count("Association.Systems")) {
      std::string association;
      for(const auto& system:words(analysis.at("Association.Systems"))) if(kept.count(system))
        association+=(association.empty()?"":" ")+system;
      if(association.empty())for(const char* key:{"Association.Systems","Association.PhotonEdges","Association.XjEdges"})analysis.erase(key);
      else analysis["Association.Systems"]=association;
    }
  }
  std::map<std::string, std::vector<std::string>> photon, inclusive;
  std::set<std::string> data;
  std::map<std::string, Json> pp;
  for (const auto &sample : collection.at("samples")) {
    const std::string system = sample.at("system"), name = sample.at("name"), family = sample.at("family");
    if (sample.contains("run")) continue;  /* retained for offline QA; the combined DATA is the analysis input */
    if (family == "data") data.insert(system);
    else (family == "photonjet" ? photon : inclusive)[system].push_back(name);
    if (system == "PP" && family != "data") pp[name] = sample;
  }
  const std::string suffix = "_double";
  for (const auto &[name, sample] : pp) {
    if (sample.at("simulation_periods") != "all") continue;
    const bool di = sample.at("si_di_role") == "DI";
    if (!di && sample.at("si_di_role") != "SI") throw std::runtime_error("unknown p+p interaction role");
    const bool named = name.size() > suffix.size() && name.substr(name.size() - suffix.size()) == suffix;
    const std::string partner = !di ? name + suffix : named ? name.substr(0, name.size() - suffix.size()) : "";
    if (partner.empty() || !pp.count(partner) || pp.at(partner).at("family") != sample.at("family") ||
        pp.at(partner).at("si_di_role") != (di ? "SI" : "DI") || pp.at(partner).at("simulation_periods") != "all")
      throw std::runtime_error("p+p two-period analysis needs both SI and DI packages for " + name);
  }
  auto joined = [](const std::vector<std::string> &names) {
    std::string out;
    for (const auto &name : names) out += (out.empty() ? "" : " ") + name;
    return out.empty() ? "none" : out;
  };
  for (const auto &system : words(required(analysis, "Systems"))) {
    const auto package = required(analysis, system + ".Package");
    if (!data.count(package) || photon[package].empty())
      throw std::runtime_error("collection needs DATA and photon simulation for " + package);
    analysis[system + ".Simulation"] = joined(photon[package]);
    analysis["QA." + package + ".Inclusive"] = joined(inclusive[package]);
  }
}

}  /* namespace */

/* The twin of a validated collection for everyday running: every file is first checked against the collection
   (the one full read of the delivered DATA packages), then copied, except that a DATA package is written without
   records/source_contributions and records/sources and with derived/ (QA.cc sourceStudies) in their place.
   <package>.compact.json binds each twin to the delivered sha256; histogramCollection accepts exactly that. */
int compactCollection(const std::string &collectionPath, const std::string &outputPath) {
  const fs::path input = fs::canonical(collectionPath), out = fs::absolute(outputPath).lexically_normal();
  if (fs::exists(out)) throw std::runtime_error("the output directory exists; choose a new one");
  fs::create_directories(out);
  const HistogramCollection c = histogramCollection(input);
  if(c.manifest.contains("measurement_replay"))throw std::runtime_error("measurement replay collection is already a restricted payload; do not compact as nominal QA");
  if (c.manifest.is_null()) throw std::runtime_error("compact needs a collection.yaml");
  std::map<std::string, std::string> data;  /* path -> system */
  for (const auto &sample : c.manifest.at("samples"))
    if (sample.at("family") == "data") data[sample.at("path").get<std::string>()] = sample.at("system").get<std::string>();
  const std::set<std::string> dropped = {"records/source_contributions", "records/sources"};
  for (const auto &file : c.manifest.at("files")) {
    const std::string name = file.at("path"), expected = file.at("sha256");
    const fs::path from = c.directory / name, to = out / name;
    fs::create_directories(to.parent_path());
    if (!data.count(name)) {
      fs::copy_file(from, to);
      continue;
    }
    std::unique_ptr<TFile> source(TFile::Open(from.c_str(), "READ"));
    if (!source || source->IsZombie()) throw std::runtime_error("cannot open " + from.string());
    if (source->GetDirectory("derived")) {  /* already a twin: its binding travels with it */
      source.reset();
      fs::copy_file(from, to);
      fs::path binding = from, copy = to;
      binding.replace_extension(".compact.json");
      copy.replace_extension(".compact.json");
      fs::copy_file(binding, copy);
      continue;
    }
    std::unique_ptr<TFile> target(TFile::Open(to.c_str(), "NEW"));
    if (!target || target->IsZombie()) throw std::runtime_error("cannot write " + to.string());
    copyDirectory(source.get(), target.get(), dropped);
    TDirectory *derived = target->mkdir("derived");
    sourceStudies(source.get(), data.at(name), derived);
    derived->cd();
    TObjString(kSourceStudies).Write("version");
    Json record = {{"schema", "PhotonJetCompactPackageV1"}, {"source_sha256", expected}, {"source_bytes", fs::file_size(from)},
                   {"dropped", dropped}, {"derived", kSourceStudies}};
    TDirectory *meta = target->GetDirectory("meta");
    if (!meta) throw std::runtime_error(from.string() + ": no meta directory");
    meta->cd();
    TObjString(record.dump().c_str()).Write("compact");
    target->Close();
    source.reset();
    record["compact_sha256"] = hash(to);
    record["compact_bytes"] = fs::file_size(to);
    fs::path binding = to;
    binding.replace_extension(".compact.json");
    writeJson(binding, record);
    std::cout << "[FinalAnalysis compact] " << name << ": " << record["source_bytes"] << " -> " << record["compact_bytes"]
              << " bytes" << std::endl;
  }
  fs::copy_file(input, out / input.filename());
  std::cout << "[FinalAnalysis compact] wrote " << (out / input.filename()).string() << std::endl;
  return 0;
}

/* ---------------------------------------------------------------------------------------------------------------
   4. The pipeline
   --------------------------------------------------------------------------------------------------------------- */
namespace {

/* What `all` was asked to do: the pipeline settings, read and checked before anything is written. */
struct Plan {
  fs::path file, base;             /* the pipeline settings file and its folder */
  fs::path package, treeToHists;   /* FinalAnalysis/ and TreeToHists/ */
  Config settings;                 /* the pipeline settings */
  bool trees = false;              /* from eventTrees (a local reduction first) instead of a histogram collection */
  fs::path input;                  /* the eventTree collection, or the histogram collection */
  fs::path nominalPath;            /* the nominal analysis settings */
  Config nominal;
  fs::path reductionPath;          /* trees: the nominal reduction settings */
  Config nominalReduction;
  std::string systematics, expected;  /* the systematics settings and the list of figures that must exist, or "none" */
  Config systematicsSettings;      /* read with the plan, so a broken file stops the run before it starts */
  std::vector<std::string> variations, expectedFigures;
  std::string parallel = "1";      /* how many variations run at once */
  bool collectionSamples = false;  /* Pipeline.UseCollectionSamples */
  std::vector<std::pair<std::string, fs::path>> radii;  /* Pipeline.RadiusCollections: label and collection.yaml */
  fs::path output;

  /* the settings of this variation alone start with <name>.Analysis. or <name>.Reduction. */
  bool sets(const std::string &name, const char *kind) const {
    const std::string prefix = name + "." + kind + ".";
    return std::any_of(settings.begin(), settings.end(), [&](const auto &item) { return item.first.rfind(prefix, 0) == 0; });
  }
};

Plan readPlan(const std::string &file, const std::string &home) {
  Plan p;
  p.file = fs::absolute(file);
  p.base = p.file.parent_path();
  p.package = fs::absolute(home);
  p.treeToHists = p.package.parent_path() / "TreeToHists";
  p.settings = readSettings(p.file);
  const Config &s = p.settings;
  p.trees = s.count("Pipeline.InputCollection");
  if (p.trees == bool(s.count("Pipeline.HistogramCollection")))
    throw std::runtime_error("give one input collection or histogram collection");
  p.input = resolve(p.base, required(s, p.trees ? "Pipeline.InputCollection" : "Pipeline.HistogramCollection"));
  p.nominalPath = resolve(p.base, required(s, "Pipeline.AnalysisSettings"));
  p.nominal = readSettings(p.nominalPath);
  if (p.trees) {
    p.reductionPath = resolve(p.base, required(s, "Pipeline.ReductionSettings"));
    p.nominalReduction = reductionSettings(p.reductionPath);
  }
  p.systematics = required(s, "Pipeline.SystematicsSettings");
  p.expected = required(s, "Pipeline.ExpectedFigures");
  p.variations = words(required(s, "Pipeline.Variations") == "none" ? "" : required(s, "Pipeline.Variations"));
  std::set<std::string> unique;
  for (const auto &v : p.variations)
    if (v == "nominal" || v.find_first_not_of(kNameCharacters) != std::string::npos || !unique.insert(v).second)
      throw std::runtime_error("invalid or duplicate variation " + v);
  if (p.systematics == "none" && !p.variations.empty()) throw std::runtime_error("variations require systematics settings");
  if (p.systematics != "none") p.systematicsSettings = readSettings(resolve(p.base, p.systematics));

  const std::set<std::string> known = {
      "Pipeline.InputCollection", "Pipeline.HistogramCollection", "Pipeline.ReductionSettings", "Pipeline.AnalysisSettings",
      "Pipeline.SystematicsSettings", "Pipeline.ExpectedFigures", "Pipeline.Variations", "Pipeline.OutputDirectory",
      "Pipeline.Parallel", "Pipeline.UseCollectionSamples", "Pipeline.PhotonIDModels", "Pipeline.RadiusCollections"};
  for (const auto &[k, v] : s) {
    bool ok = known.count(k);
    for (const auto &name : p.variations) ok = ok || k.rfind(name + ".Analysis.", 0) == 0 || k.rfind(name + ".Reduction.", 0) == 0 || k == name + ".HistogramCollection";
    if (!ok) throw std::runtime_error("unknown pipeline setting " + k);
  }
  /* Refuse before creating output or launching the nominal reduction. */
  auto hasNamed=[](const Config &config,const std::string &prefix) {
    for(const auto &[key,value]:config)if(key.rfind(prefix,0)==0)
      for(const std::string suffix:{".ResponseComponent",".ResponseRecipe",".ResponseGridSHA256"})
        if(key.size()>suffix.size() && key.compare(key.size()-suffix.size(),suffix.size(),suffix)==0)return true;
    return false;
  };
  const bool nominalNamed=hasNamed(p.nominal,"");
  if(nominalNamed && p.trees)throw std::runtime_error("named response selection requires saved histogram input without eventTree reduction");
  for(const auto &name:p.variations)
    if((nominalNamed || hasNamed(s,name+".Analysis.")) && (p.trees || p.sets(name,"Reduction")))
      throw std::runtime_error("named response selection requires saved histogram input without eventTree reduction");
  if (s.count("Pipeline.Parallel")) p.parallel = s.at("Pipeline.Parallel");
  if (p.parallel.empty() || p.parallel.find_first_not_of("0123456789") != std::string::npos || std::stol(p.parallel) < 1 ||
      std::stol(p.parallel) > 256)
    throw std::runtime_error("Pipeline.Parallel must be an integer in 1..256");
  const std::string samples = s.count("Pipeline.UseCollectionSamples") ? s.at("Pipeline.UseCollectionSamples") : "0";
  if (samples != "0" && samples != "1") throw std::runtime_error("Pipeline.UseCollectionSamples must be 0 or 1");
  if (p.trees && samples == "1") throw std::runtime_error("Pipeline.UseCollectionSamples applies only to saved histogram collections");
  p.collectionSamples = samples == "1";
  /* Pipeline.RadiusCollections: "r02=COLLECTION.yaml r03=COLLECTION.yaml", the histogram collections of the other jet
     radii. The nominal analysis is repeated on each (results/radius_<label>.root) for the radius overlay. */
  if (s.count("Pipeline.RadiusCollections") && s.at("Pipeline.RadiusCollections") != "none")
    for (const auto &item : words(s.at("Pipeline.RadiusCollections"))) {
      const auto equals = item.find('=');
      const std::string label = item.substr(0, equals);
      if (equals == std::string::npos || label.empty() || label.find_first_not_of(kNameCharacters) != std::string::npos || p.trees)
        throw std::runtime_error("Pipeline.RadiusCollections needs label=COLLECTION.yaml items and a histogram collection");
      for (const auto &[other, path] : p.radii) if (other == label) throw std::runtime_error("duplicate radius label " + label);
      p.radii.push_back({label, resolve(p.base, item.substr(equals + 1))});
    }

  if (p.expected != "none") {
    std::ifstream list(resolve(p.base, p.expected));
    if (!list) throw std::runtime_error("cannot read expected figure list");
    for (std::string line; std::getline(list, line);) {
      const auto tokens = words(line.substr(0, line.find('#')));
      if (tokens.empty()) continue;
      if (tokens.size() != 1 || fs::path(tokens[0]).is_absolute() || tokens[0].find("..") != std::string::npos ||
          fs::path(tokens[0]).extension() != ".png")
        throw std::runtime_error("expected figures need safe relative PNG paths");
      p.expectedFigures.push_back(tokens[0]);
    }
    if (p.expectedFigures.empty()) throw std::runtime_error("expected figure list is empty");
  }
  if (!fs::exists(p.input)) throw std::runtime_error("input collection missing");
  p.output = resolve(p.base, required(s, "Pipeline.OutputDirectory"));
  return p;
}

/* A pipeline coverage template pins the just-produced nominal, while the
   numerical combiner still receives its ordinary exact-hash coverage contract. */
Json bindCoverageNominal(Json coverage,const std::string &nominalHash) {
  if(coverage.value("schema","")!="PhotonJetSystematicCoverageV1" || !coverage.contains("nominal_sha256") ||
      !coverage.at("nominal_sha256").is_string())throw std::runtime_error("invalid systematic coverage template");
  if(coverage.at("nominal_sha256")=="@nominal")coverage["nominal_sha256"]=nominalHash;
  return coverage;
}

/* The analysis settings of the nominal run or of one variation, before the paths of this run are filled in.

   A variation is run for its central values, which is all the systematics step takes from it. The statistical
   errors, the response toys and the closure are the nominal's, so they are switched off here unless the variation sets
   them itself. Its data toys stay: the prior fit uses their covariance, so the central values are those of a full
   run. Its iterations are the nominal result's choice (or one less, one more): the iteration scan is never run again
   under a variation, where it could select another k. */
Config analysisSettings(const Plan &plan, const std::string &name, const Json &manifest, const fs::path &nominalResults) {
  Config analysis = plan.nominal;
  if (plan.collectionSamples) useCollectionSamples(analysis, manifest);
  /* These opt-in fields have no nominal placeholder. Admit them only on an
     existing result; the ordinary overlay still rejects misspelled keys. */
  for(const auto &system:words(required(analysis,"Systems")))
    for(const char *field:{"ResponseComponent","ResponseRecipe","ResponseGridSHA256"}) {
      const std::string key=system+"."+field;
      if(plan.settings.count(name+".Analysis."+key))analysis.emplace(key,"");
    }
  // Saved measurement views can omit nominal model diagnostics through the existing QA interface.
  for(const auto *field:{"QA.Mode","QA.PhotonIDModels"})
    if(plan.settings.count(name+".Analysis."+field))analysis.emplace(field,"");
  overlay(analysis, plan.settings, name + ".Analysis.");
  if (name != "nominal") {
    if (!plan.settings.count(name + ".Analysis.Uncertainty.ResponseToys")) analysis["Uncertainty.ResponseToys"] = "0";
    if (!plan.settings.count(name + ".Analysis.Closure.Toys")) analysis["Closure.Toys"] = "0";
  }
  for (auto &[key, value] : analysis) {
    const auto dot = key.find('.');
    const std::string field = dot == std::string::npos ? std::string() : key.substr(dot + 1);
    const bool iterations = field == "PairIterations" || field == "PhotonIterations";
    if (name != "nominal" && iterations && value == "scan") value = "nominal";
    if (value != "nominal-1" && value != "nominal+1" && !(iterations && value == "nominal")) continue;
    if (name == "nominal" || !iterations)
      throw std::runtime_error("relative nominal iterations apply only to variation iteration settings");
    /* the k the nominal result used, from its saved settings */
    TFile reference(nominalResults.c_str(), "READ");
    auto *saved = dynamic_cast<TObjString *>(reference.Get((key.substr(0, dot) + "/settings").c_str()));
    if (!saved) throw std::runtime_error("nominal iteration result missing for " + key);
    const std::string marker = field == "PairIterations" ? "pair_iterations=" : "photon_iterations=";
    const std::string text = saved->GetString().Data();
    const auto at = text.find(marker);
    if (at == std::string::npos) throw std::runtime_error("nominal iteration metadata missing");
    const int k = std::stoi(text.substr(at + marker.size())) + (value == "nominal-1" ? -1 : value == "nominal+1" ? 1 : 0);
    if (k < 1) throw std::runtime_error("iteration variation is outside positive iterations");
    value = std::to_string(k);
  }
  return analysis;
}

/* Whether a variation can use the nominal QA file. The QA step reads the packages, each result's samples and photon
   bins, and the recoil-jet slot. A variation that changes none of them (prior, unfolding, toys, combinatoric
   estimate, iterations, purity inputs) has the nominal QA, and its QA step is not repeated. */
bool usesNominalQA(const Plan &plan, const std::string &name) {
  if (name == "nominal" || plan.sets(name, "Reduction") || plan.settings.count(name+".HistogramCollection")) return false;
  const std::set<std::string> sections = {"Prior", "Unfolding", "IterationRule", "Uncertainty", "Closure", "Combinatoric"};
  const std::set<std::string> fields = {"PairIterations", "PhotonIterations", "BackgroundCorrelation", "PurityScale", "EfficiencyScale", "ResponseComponent", "ResponseRecipe", "ResponseGridSHA256"};
  const std::string prefix = name + ".Analysis.";
  for (const auto &[key, value] : plan.settings) {
    if (key.rfind(prefix, 0) != 0) continue;
    const std::string rest = key.substr(prefix.size());
    if (!sections.count(rest.substr(0, rest.find('.'))) && !fields.count(rest.substr(rest.rfind('.') + 1))) return false;
  }
  return true;
}

/* Bind the existing collection to the generated analysis recipe. This helper
   runs before any analysis/QA command and refuses an eventTree reduction route. */
Json bindNamedResponseInput(Config &analysis,const HistogramCollection &selected,bool trees,bool reduction) {
  const bool active=std::any_of(analysis.begin(),analysis.end(),[](const auto &entry){for(const std::string suffix:{".ResponseComponent",".ResponseRecipe",".ResponseGridSHA256"})if(entry.first.size()>suffix.size() && entry.first.compare(entry.first.size()-suffix.size(),suffix.size(),suffix)==0)return true;return false;});
  if(!active)return Json();
  if(trees || reduction || selected.manifest.is_null())throw std::runtime_error("named response selection requires an existing typed histogram collection, without eventTree reduction");
  const auto manifest=fs::canonical(selected.files.front());
  if(analysis.count("Input.Collection") && fs::canonical(analysis.at("Input.Collection"))!=manifest)throw std::runtime_error("analysis Input.Collection differs from its pipeline collection");
  analysis["Input.Collection"]=manifest.string();
  return {{"collection",manifest.string()},{"collection_sha256",hash(manifest)},{"scope","RESPONSE_COMPONENT_ONLY"}};
}

}  /* namespace */

/* Runs the package's commands in dependency order. Each has a log and stops the chain on failure; PIPELINE.json binds
   the settings, the inputs and the products of a completed run. */
void validateNamedResponseCollection(Settings &settings) {
  const bool active=std::any_of(settings.systems.begin(),settings.systems.end(),named_response::active);
  if(!active){if(!settings.inputCollection.empty())throw std::runtime_error("Input.Collection is only used with explicit named response selection");return;}
  named_response::require(!settings.inputCollection.empty(),"Input.Collection is required");
  named_response::require(settings.combinatoricMethod=="data_dphi" && settings.combinatoricShape=="simulation","consumer requires Combinatoric.Method=data_dphi and Combinatoric.Shape=simulation; other correction modes are unqualified");
  const auto collection=histogramCollection(settings.inputCollection);
  named_response::require(!collection.manifest.is_null() && fs::canonical(settings.inputDirectory)==collection.directory,"component input directory differs from its typed collection");
  const auto contract=collection.manifest.value("simulation_response",Json());
  named_response::require(!collection.manifest.contains("measurement_replay") &&
    (contract.is_null() || contract.at("simulation_response_variation")==""),"component selection cannot overlay a collection response/replay recipe");
  settings.namedResponsePins.clear();
  for(const auto &sys:settings.systems)if(named_response::active(sys)) {
    named_response::validateChoice(sys);named_response::require(!sys.simulation.empty(),"selected result has no simulation");
    std::set<std::string> names;
    for(const auto &name:sys.simulation) {
      named_response::require(names.insert(name).second,"duplicate selected simulation");
      const std::string key=sys.package+"/"+name;
      const auto sample=std::find_if(collection.manifest.at("samples").begin(),collection.manifest.at("samples").end(),[&](const Json &row){return row.at("system")==sys.package && row.at("name")==name && row.at("family")=="photonjet";});
      named_response::require(sample!=collection.manifest.at("samples").end(),"selected simulation is not a signal-SIM sample of the typed collection: "+key);
      const auto path=fs::canonical(collection.directory/(key+".root"));TFile file(path.c_str(),"READ");
      named_response::require(!file.IsZombie(),"cannot open selected package");named_response::verifyPackage(file,sys);
      settings.namedResponsePins[key]=collection.pins.at(path.string());
    }
  }
}

int pipeline(const std::string &file, const std::string &home) {
  const Plan plan = readPlan(file, home);
  const fs::path &output = plan.output, &package = plan.package;
  const std::string command = (package / "final_analysis.sh").string();
  fs::create_directories(output.parent_path());
  if (!fs::create_directory(output)) throw std::runtime_error("output directory already exists; choose a new version");
  for (const char *d : {"config", "logs", "results", "fig", "tables"}) fs::create_directory(output / d);
  std::signal(SIGINT, terminate);
  std::signal(SIGTERM, terminate);

  Json record = {{"schema", "PhotonJetPipelineRunV1"}, {"status", "RUNNING"}, {"root_version", gROOT->GetVersion()},
                 {"physics_acceptance", "NOT_ESTABLISHED_BY_SOFTWARE_RUN"}, {"steps", Json::array()},
                 {"input_pins", Json::object()}};
  auto save = [&]() { writeJson(output / "PIPELINE.json", record); };
  auto pin = [&](const fs::path &p) { record["input_pins"][p.string()] = hash(p); };
  /* The code that computes (src) is pinned at launch; the code that draws (plots) when drawing starts, so figures can
     be worked on while results are computed. Both are checked again at the end. */
  const std::vector<fs::path> codeRoots = plan.trees ? std::vector<fs::path>{package, plan.treeToHists} : std::vector<fs::path>{package};
  auto pinCode = [&](const char *folder) {
    const std::set<std::string> code = {".cc", ".h", ".C", ".json", ".py"};
    for (const auto &root : codeRoots) {
      if (!fs::is_directory(root / folder)) continue;
      for (const auto &e : fs::recursive_directory_iterator(root / folder))
        if (e.is_regular_file() && code.count(e.path().extension().string())) pin(e.path());
    }
  };
  /* Runs steps at once (their names are their log names), records each and stops at the first that failed. */
  auto steps = [&](const std::vector<std::pair<std::string, std::vector<std::string>>> &named, size_t limit) {
    if (interrupted) throw std::runtime_error("interrupted");
    save();
    std::vector<Command> commands;
    for (const auto &[name, arguments] : named) commands.push_back({name, arguments, output / "logs" / (name + ".log")});
    const std::vector<int> codes = run(commands, limit);
    std::string failure;
    for (size_t i = 0; i < commands.size(); ++i) {
      std::cout << "[FinalAnalysis all] " << commands[i].name << " exit=" << codes[i] << " log=" << commands[i].log << std::endl;
      record["steps"].push_back({{"name", commands[i].name}, {"exit_code", codes[i]}, {"log", commands[i].log.filename().string()}});
      if (codes[i] && failure.empty()) failure = commands[i].name + " failed; see " + commands[i].log.string();
    }
    save();
    if (!failure.empty()) throw std::runtime_error(failure);
  };
  auto step = [&](const std::string &name, const std::vector<std::string> &arguments) { steps({{name, arguments}}, 1); };

  try {
    /* what goes in */
    pin(plan.file);
    pin(plan.nominalPath);
    if (plan.trees) pin(plan.input), pin(plan.reductionPath);
    const HistogramCollection collection =
        plan.trees ? HistogramCollection{plan.input, Json(), {}, {}, Json::array()} : histogramCollection(plan.input);
    if(collection.manifest.contains("measurement_replay"))throw std::runtime_error("the nominal analysis input cannot be a replay variation");
    for (const auto &[path, digest] : collection.pins) record["input_pins"][path] = digest;
    if (!collection.manifest.is_null()) record["histogram_collection"] = collection.manifest;
    if (!collection.otherSelections.empty())
      record["selection_compatibility"] = {{"rule", "identical except the p+p DATA trigger name, which these samples do not evaluate"},
                                           {"collection_selection_sha256", collection.manifest.at("selection_sha256")},
                                           {"samples", collection.otherSelections}};
    if (plan.systematics != "none") pin(resolve(plan.base, plan.systematics));
    if (plan.expected != "none") pin(resolve(plan.base, plan.expected));
    for (const auto &root : codeRoots) pin(root / (root == package ? "final_analysis.sh" : "tree_to_hists.sh"));
    pinCode("src");
    for (const auto &file : measurementReplayCodeFiles()) pin(file);
    pin(named_response::codeFile());
    if (!plan.trees) {
      pin(plan.treeToHists / "src/Configuration.h");
      pin(package.parent_path() / "TreeProduction/src/internal/Configuration.h");
    }

    /* the nominal run, then every variation: its settings, its QA unless the nominal's serves, its measurement. The
       nominal measurement runs first because a variation's iterations follow it; the variations' then run at once. */
    std::map<std::string, fs::path> results;
    std::vector<std::pair<std::string, std::vector<std::string>>> variationRuns;
    fs::path nominalPackages = collection.directory;
    std::vector<std::string> names = {"nominal"};
    names.insert(names.end(), plan.variations.begin(), plan.variations.end());
    for (const auto &name : names) {
      const bool ownReduction = plan.sets(name, "Reduction");
      if (ownReduction && !plan.trees) throw std::runtime_error("response variations require an input eventTree collection");
      const bool ownCollection=plan.settings.count(name+".HistogramCollection");
      if (ownCollection && (ownReduction || plan.trees))
        throw std::runtime_error("a saved variation collection cannot also request eventTree reduction");
      HistogramCollection varied;
      if (ownCollection) {
        const fs::path manifest=resolve(plan.base,plan.settings.at(name+".HistogramCollection"));
        if (fs::is_directory(manifest)) throw std::runtime_error("a variation requires a validated collection.yaml");
        varied=histogramCollection(manifest);
        bindMeasurementReplay(varied,collection);
        for (const auto &[file,digest]:varied.pins) record["input_pins"][file]=digest;
        record["variation_collections"][name]={{"path",manifest.string()},{"selection_sha256",varied.manifest.at("selection_sha256")},{"other_selections",varied.otherSelections}};
        if(varied.manifest.contains("measurement_replay")) {
          record["variation_collections"][name]["measurement_replay"]=varied.manifest.at("measurement_replay");
          record["variation_collections"][name]["nominal_diagnostics"]=(output/"results/nominal_qa.root").string();
        }
      }
      Config analysis = analysisSettings(plan, name, ownCollection ? varied.manifest : collection.manifest, name == "nominal" ? fs::path() : results.at("nominal"));
      fs::path packages = nominalPackages;
      if (ownCollection) packages=varied.directory;
      if (plan.trees && (name == "nominal" || ownReduction)) {
        Config reduction = plan.nominalReduction;
        overlay(reduction, plan.settings, name + ".Reduction.");
        const auto settings = output / "config" / (name + "_reduction.yaml");
        writeSettings(settings, reduction);
        packages = output / ("packages_" + name);
        step(name + "_reduction", {(plan.treeToHists / "tree_to_hists.sh").string(), "local", "--settings", settings.string(),
                                   "--inputs", plan.input.string(), "--work", (output / ("work_" + name)).string(),
                                   "--out", packages.string(), "--parallel", plan.parallel});
        if (name == "nominal") nominalPackages = packages;
      }
      if ((name == "nominal" || ownReduction) && collection.manifest.is_null())
        for (const auto &e : fs::recursive_directory_iterator(packages))
          if (e.is_regular_file() && e.path().extension() == ".root") pin(e.path());

      const bool nominalQA = usesNominalQA(plan, name);
      const fs::path settings = output / "config" / (name + "_analysis.yaml");
      results[name] = output / "results" / (name + ".root");
      analysis["Input.Directory"] = packages.string();
      const auto componentInput=bindNamedResponseInput(analysis,ownCollection?varied:collection,plan.trees,ownReduction);
      if(!componentInput.is_null())record["named_response_components"][name]=componentInput;
      analysis["Output.Results"] = results[name].string();
      analysis["Output.QA"] = (output / "results" / ((nominalQA ? "nominal" : name) + "_qa.root")).string();
      const bool replayInputs=ownCollection && varied.manifest.contains("measurement_replay");
      if(replayInputs) {
        analysis["QA.Mode"]="measurement_inputs";
        analysis["QA.PhotonIDModels"]="none";
      }
      if (!replayInputs && plan.settings.count("Pipeline.PhotonIDModels") && plan.settings.at("Pipeline.PhotonIDModels") != "none") {
        const auto models = resolve(plan.base, plan.settings.at("Pipeline.PhotonIDModels"));
        pin(models);
        analysis["QA.PhotonIDModels"] = models.string();
      }
      writeSettings(settings, analysis);
      if (!nominalQA) step(name + "_qa", {command, "qa", "--config", settings.string()});
      const std::vector<std::string> measurement = {command, "run", "--config", settings.string()};
      if (name == "nominal" || ownReduction) step(name + "_analysis", measurement);
      else variationRuns.push_back({name + "_analysis", measurement});
    }
    if (!variationRuns.empty()) steps(variationRuns, std::stol(plan.parallel));

    /* the nominal analysis on the collection of each other jet radius: its own QA file and results beside the
       nominal's, with every toy and its own iteration choice, so each radius is a complete measurement */
    std::vector<std::pair<std::string, std::vector<std::string>>> radiusQA, radiusRuns;
    for (const auto &[label, path] : plan.radii) {
      pin(path);
      const HistogramCollection other = histogramCollection(path);
      for (const auto &[file, digest] : other.pins) record["input_pins"][file] = digest;
      Config analysis = analysisSettings(plan, "nominal", other.manifest, fs::path());
      const std::string name = "radius_" + label;
      const fs::path settings = output / "config" / (name + "_analysis.yaml");
      analysis["Input.Directory"] = other.directory.string();
      analysis["Output.Results"] = (output / "results" / (name + ".root")).string();
      analysis["Output.QA"] = (output / "results" / (name + "_qa.root")).string();
      if (plan.settings.count("Pipeline.PhotonIDModels") && plan.settings.at("Pipeline.PhotonIDModels") != "none")
        analysis["QA.PhotonIDModels"] = resolve(plan.base, plan.settings.at("Pipeline.PhotonIDModels")).string();
      writeSettings(settings, analysis);
      radiusQA.push_back({name + "_qa", {command, "qa", "--config", settings.string()}});
      radiusRuns.push_back({name + "_analysis", {command, "run", "--config", settings.string()}});
    }
    if (!radiusQA.empty()) steps(radiusQA, std::stol(plan.parallel)), steps(radiusRuns, std::stol(plan.parallel));

    /* the systematic uncertainties: "@name" in the systematics settings is the result of that variation */
    fs::path systematics;
    if (!plan.systematicsSettings.empty()) {
      const fs::path source = resolve(plan.base, plan.systematics);
      Config sys = plan.systematicsSettings;
      systematics = output / "results/systematics.root";
      sys["Systematics.Nominal"] = results.at("nominal").string();
      sys["Systematics.Output"] = systematics.string();
      auto result = [&](const std::string &value) {
        if (value.empty()) throw std::runtime_error("empty systematic result binding");
        if (value[0] == '@') {
          const auto name = value.substr(1);
          if (name == "nominal" || !results.count(name)) throw std::runtime_error("systematic names an undeclared variation: " + name);
          return results.at(name).string();
        }
        const auto external = resolve(source.parent_path(), value);
        pin(external);
        return external.string();
      };
      for (auto &[key, value] : sys) {
        const auto suffix = key.substr(key.rfind('.') + 1);
        if (suffix == "Up" || suffix == "Down" || suffix == "Alternative") value = result(value);
        if (suffix == "Alternatives") {  /* label=result label=result ... */
          std::string rebuilt;
          for (const auto &item : words(value)) {
            const auto equals = item.find('=');
            if (equals == std::string::npos) throw std::runtime_error("systematic alternative needs label=result");
            rebuilt += (rebuilt.empty() ? "" : " ") + item.substr(0, equals + 1) + result(item.substr(equals + 1));
          }
          value = rebuilt;
        }
      }
      if(sys.count("Systematics.Coverage")) {
        const auto templatePath=resolve(source.parent_path(),sys.at("Systematics.Coverage"));
        pin(templatePath);
        std::ifstream input(templatePath);Json coverage;input>>coverage;
        coverage=bindCoverageNominal(coverage,hash(results.at("nominal")));
        const auto coveragePath=output/"config/systematics_coverage.json";
        writeJson(coveragePath,coverage);sys["Systematics.Coverage"]=coveragePath.string();
        record["systematics_coverage"]={{"template",templatePath.string()},{"path",coveragePath.string()},
          {"nominal_sha256",coverage.at("nominal_sha256")}};
      }
      const auto settings = output / "config/systematics.yaml";
      writeSettings(settings, sys);
      step("systematics", {command, "systematics", "--config", settings.string()});
    }

    /* the figures and the tables, from the saved results alone */
    std::vector<std::string> draw = {command, "draw", "--results", results.at("nominal").string(),
                                     "--qa", (output / "results/nominal_qa.root").string(), "--fig", (output / "fig").string()};
    std::vector<std::string> tables = {command, "tables", "--results", results.at("nominal").string(), "--out", (output / "tables").string()};
    if (!systematics.empty()) {
      draw.insert(draw.end(), {"--systematics", systematics.string()});
      tables.insert(tables.end(), {"--systematics", systematics.string()});
    }
    pinCode("plots");
    pinCode("scripts");
    pin(package / "plot_physics_coverage.json");
    if (fs::is_regular_file(package / "plot_physics_qualification.json"))
      pin(package / "plot_physics_qualification.json");
    step("draw", draw);
    record["plot_physics_generation"] = physicsGeneration(output / "fig", hash);
    save();
    step("tables", tables);

    /* what came out: every listed figure exists, every expected one too, and nothing that went in has changed */
    record["figures"] = Json::array();
    record["unavailable_figures"] = Json::array();
    std::ifstream list(output / "fig/FIGURES.txt");
    if (!list) throw std::runtime_error("missing figure coverage record");
    auto drawn = [&](const fs::path &figure) { return fs::is_regular_file(figure) && fs::file_size(figure) > 0; };
    for (std::string line; std::getline(list, line);) {
      if (line.rfind("fig/", 0) == 0) {
        if (!drawn(output / line)) throw std::runtime_error("missing saved figure " + line);
        record["figures"].push_back(line);
      }
      if (line.rfind("# not drawn:", 0) == 0) record["unavailable_figures"].push_back(line.substr(13));
    }
    if (record["figures"].empty()) throw std::runtime_error("no figures generated");
    record["missing_expected_figures"] = Json::array();
    for (const auto &figure : plan.expectedFigures)
      if (!drawn(output / "fig" / figure)) record["missing_expected_figures"].push_back(figure);
    if (!record["missing_expected_figures"].empty()) throw std::runtime_error("required figure coverage incomplete");
    for (const auto &[path, original] : record["input_pins"].items())
      if (hash(path) != original.get<std::string>()) throw std::runtime_error("input changed during run: " + path);
    record["output_pins"] = Json::object();
    for (const char *folder : {"results", "fig", "tables", "config"})
      for (const auto &e : fs::recursive_directory_iterator(output / folder))
        if (e.is_regular_file()) record["output_pins"][fs::relative(e.path(), output).string()] = hash(e.path());
    record["status"] = plan.expected == "none" ? "COMPONENTS_PASSED_COVERAGE_UNQUALIFIED" : "PASS_EXPECTED_OUTPUTS";
    save();
    return 0;
  } catch (const std::exception &e) {
    record["status"] = "FAILED";
    record["first_error"] = e.what();
    save();
    throw;
  }
}

}  /* namespace fa */
