/* Pipeline: portable orchestration of the existing reduction, analysis and rendering commands.
   Physics is evaluated only by TreeToHists and the existing FinalAnalysis components. */
#include "FinalAnalysis.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>
#include <cerrno>
#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#include "TEnv.h"
#include "THashList.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObjString.h"

namespace fa {
namespace {
namespace fs = std::filesystem;
using Json = nlohmann::json;
using Config = std::map<std::string, std::string>;
volatile sig_atomic_t interrupted = 0, child = 0;
void terminate(int signal) { interrupted = signal; if (child > 0) kill(-child, signal); }
std::vector<std::string> words(const std::string &s) {
  std::istringstream in(s); std::vector<std::string> out;
  for (std::string v; in >> v;) out.push_back(v);
  return out;
}
Config read(const fs::path &path) {
  /* TEnv otherwise warns and silently keeps the first of two settings. */
  std::ifstream source(path); std::set<std::string> keys;
  for (std::string line; std::getline(source, line);) {
    line = line.substr(0, line.find('#')); const auto colon = line.find(':');
    if (colon == std::string::npos) continue;
    const auto names = words(line.substr(0, colon));
    if (names.size() == 1 && !keys.insert(names[0]).second)
      throw std::runtime_error("duplicate configuration key " + names[0] + " in " + path.string());
  }
  TEnv env;
  if (env.ReadFile(path.c_str(), kEnvUser) != 0) throw std::runtime_error("cannot read " + path.string());
  Config out; TIter next(env.GetTable());
  while (auto *r = static_cast<TEnvRec *>(next())) out[r->GetName()] = r->GetValue();
  return out;
}
std::string required(const Config &c, const std::string &key) {
  auto it = c.find(key);
  if (it == c.end() || it->second.empty()) throw std::runtime_error("missing " + key);
  return it->second;
}
fs::path resolve(const fs::path &base, const std::string &v) {
  return fs::absolute(fs::path(v).is_absolute() ? fs::path(v) : base / v).lexically_normal();
}
void atomic(const fs::path &p, const Json &j) {
  std::ofstream out(p.string() + ".partial"); out << j.dump(2) << '\n'; out.close();
  if (!out) throw std::runtime_error("cannot save " + p.string());
  fs::rename(p.string() + ".partial", p);
}
void config(const fs::path &p, const Config &c) {
  std::ofstream out(p);
  for (const auto &[k, v] : c) out << k << ": " << v << '\n';
  out.close(); if (!out) throw std::runtime_error("cannot save " + p.string());
}
/* Exact argv, no shell evaluation. A signal reaches the complete child process group. */
int execute(const std::vector<std::string> &args, const fs::path &log, bool quiet = false) {
  if (interrupted) throw std::runtime_error("interrupted before dispatch");
  const int fd = open(log.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
  if (fd < 0) throw std::runtime_error("cannot create log " + log.string());
  const pid_t pid = fork();
  if (pid == 0) {
    setpgid(0, 0); dup2(fd, STDOUT_FILENO); dup2(fd, STDERR_FILENO); close(fd);
    std::vector<char *> argv; for (const auto &a : args) argv.push_back(const_cast<char *>(a.c_str())); argv.push_back(nullptr);
    execvp(argv[0], argv.data()); _exit(127);
  }
  close(fd); if (pid < 0) throw std::runtime_error("cannot start component");
  setpgid(pid, pid); child = pid; if (interrupted) kill(-pid, interrupted);
  int status = 0;
  while (waitpid(pid, &status, 0) < 0) if (errno != EINTR) { child = 0; throw std::runtime_error("cannot reap component"); }
  child = 0;
  const int code = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
  if (!quiet) std::cout << "[FinalAnalysis all] " << args[1] << " exit=" << code << " log=" << log << std::endl;
  return code;
}
std::string hash(const fs::path &p, const fs::path &scratch) {
  if (!fs::is_regular_file(p)) throw std::runtime_error("missing pinned file " + p.string());
  fs::path log = scratch / ("hash_" + std::to_string(getpid()) + ".txt");
  int code = execute({"sha256sum", p.string()}, log, true);
  if (code == 127) { fs::remove(log); code = execute({"shasum", "-a", "256", p.string()}, log, true); }
  std::ifstream in(log); std::string value; in >> value; fs::remove(log);
  if (code || value.size() != 64 || value.find_first_not_of("0123456789abcdef") != std::string::npos)
    throw std::runtime_error("SHA256 failed for " + p.string());
  return value;
}
void overlays(Config &c, const Config &run, const std::string &prefix) {
  for (const auto &[k, v] : run) if (k.rfind(prefix, 0) == 0) {
    const auto key = k.substr(prefix.size());
    if (!c.count(key)) throw std::runtime_error("override names an absent setting: " + k);
    c[key] = v;
  }
}
Config reduction(const fs::path &path) {
  Config c = read(path);
  for (const char *k : {"Selection.PhotonIDPackage", "Selection.AuAuIsolationPackage"})
    if (c.count(k) && !c[k].empty()) c[k] = resolve(path.parent_path(), c[k]).string();
  return c;
}
}  /* namespace */

int pipeline(const std::string &file, const std::string &home) {
  const fs::path path = fs::absolute(file), base = path.parent_path(), package = fs::absolute(home), tth = package.parent_path() / "TreeToHists";
  const Config settings = read(path);
  const bool trees = settings.count("Pipeline.InputCollection");
  if (trees == bool(settings.count("Pipeline.HistogramCollection"))) throw std::runtime_error("give one input collection or histogram collection");
  const fs::path input = resolve(base, required(settings, trees ? "Pipeline.InputCollection" : "Pipeline.HistogramCollection"));
  const fs::path nominalPath = resolve(base, required(settings, "Pipeline.AnalysisSettings"));
  const Config nominal = read(nominalPath);
  const fs::path reductionPath = trees ? resolve(base, required(settings, "Pipeline.ReductionSettings")) : fs::path();
  const Config nominalReduction = trees ? reduction(reductionPath) : Config();
  const std::string systematic = required(settings, "Pipeline.SystematicsSettings"), expected = required(settings, "Pipeline.ExpectedFigures");
  const auto variations = words(required(settings, "Pipeline.Variations") == "none" ? "" : required(settings, "Pipeline.Variations"));
  std::set<std::string> unique;
  for (const auto &v : variations) if (v == "nominal" || v.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos || !unique.insert(v).second)
    throw std::runtime_error("invalid or duplicate variation " + v);
  if (systematic == "none" && !variations.empty()) throw std::runtime_error("variations require systematics settings");
  Config sys = systematic == "none" ? Config() : read(resolve(base, systematic));
  std::set<std::string> known = {"Pipeline.InputCollection", "Pipeline.HistogramCollection", "Pipeline.ReductionSettings", "Pipeline.AnalysisSettings",
    "Pipeline.SystematicsSettings", "Pipeline.ExpectedFigures", "Pipeline.Variations", "Pipeline.OutputDirectory", "Pipeline.Parallel"};
  for (const auto &[k, v] : settings) {
    bool ok = known.count(k);
    for (const auto &name : variations) ok = ok || k.rfind(name + ".Analysis.", 0) == 0 || k.rfind(name + ".Reduction.", 0) == 0;
    if (!ok) throw std::runtime_error("unknown pipeline setting " + k);
  }
  const std::string parallel = settings.count("Pipeline.Parallel") ? settings.at("Pipeline.Parallel") : "1";
  if (parallel.empty() || parallel.find_first_not_of("0123456789") != std::string::npos || std::stol(parallel) < 1 || std::stol(parallel) > 256)
    throw std::runtime_error("Pipeline.Parallel must be an integer in 1..256");
  std::vector<std::string> expectedFiles;
  if (expected != "none") {
    std::ifstream list(resolve(base, expected)); if (!list) throw std::runtime_error("cannot read expected figure list");
    for (std::string l; std::getline(list, l);) {
      l = l.substr(0, l.find('#')); const auto tokens = words(l); if (tokens.empty()) continue;
      if (tokens.size() != 1 || fs::path(tokens[0]).is_absolute() || tokens[0].find("..") != std::string::npos || fs::path(tokens[0]).extension() != ".png")
        throw std::runtime_error("expected figures need safe relative PNG paths");
      expectedFiles.push_back(tokens[0]);
    }
    if (expectedFiles.empty()) throw std::runtime_error("expected figure list is empty");
  }
  if (!fs::exists(input)) throw std::runtime_error("input collection missing");
  const fs::path output = resolve(base, required(settings, "Pipeline.OutputDirectory"));
  fs::create_directories(output.parent_path());
  if (!fs::create_directory(output)) throw std::runtime_error("output directory already exists; choose a new version");
  for (const char *d : {"config", "logs", "results", "fig", "tables"}) fs::create_directory(output / d);
  std::signal(SIGINT, terminate); std::signal(SIGTERM, terminate);
  Json record = {{"schema", "PhotonJetPipelineRunV1"}, {"status", "RUNNING"}, {"root_version", gROOT->GetVersion()},
    {"physics_acceptance", "NOT_ESTABLISHED_BY_SOFTWARE_RUN"}, {"steps", Json::array()}, {"input_pins", Json::object()}};
  auto pin = [&](const fs::path &p) { record["input_pins"][p.string()] = hash(p, output / "logs"); };
  auto step = [&](const std::string &name, const std::vector<std::string> &args) {
    atomic(output / "PIPELINE.json", record);
    const auto log = output / "logs" / (name + ".log"); const int code = execute(args, log);
    record["steps"].push_back({{"name", name}, {"exit_code", code}, {"log", log.filename().string()}});
    atomic(output / "PIPELINE.json", record);
    if (code) throw std::runtime_error(name + " failed; see " + log.string());
  };
  try {
    pin(path); pin(nominalPath);
    if (trees) { pin(input); pin(reductionPath); }
    if (systematic != "none") pin(resolve(base, systematic));
    if (expected != "none") pin(resolve(base, expected));
    for (const auto &root : {package, tth}) {
      pin(root / (root == package ? "final_analysis.sh" : "tree_to_hists.sh"));
      for (const auto &dir : {"src", "plots"}) if (fs::is_directory(root / dir)) for (const auto &e : fs::recursive_directory_iterator(root / dir))
        if (e.is_regular_file() && (e.path().extension() == ".cc" || e.path().extension() == ".h" || e.path().extension() == ".C" || e.path().extension() == ".json")) pin(e.path());
    }
    std::map<std::string, fs::path> results;
    fs::path nominalPackages = input;
    for (const auto &name : [&]() { std::vector<std::string> n{"nominal"}; n.insert(n.end(), variations.begin(), variations.end()); return n; }()) {
      Config analysis = nominal, fill = nominalReduction;
      const bool differentReduction = std::any_of(settings.begin(), settings.end(), [&](const auto &item) { return item.first.rfind(name + ".Reduction.", 0) == 0; });
      if (differentReduction && !trees) throw std::runtime_error("response variations require an input eventTree collection");
      overlays(analysis, settings, name + ".Analysis."); overlays(fill, settings, name + ".Reduction.");
      /* Iteration variations use the nominal result's choice, never rerun its rule with a different seed. */
      for (auto &[key, value] : analysis) if (value == "nominal-1" || value == "nominal+1") {
        const auto dot = key.find('.'); const auto field = key.substr(dot + 1);
        if (name == "nominal" || dot == std::string::npos || (field != "PairIterations" && field != "PhotonIterations"))
          throw std::runtime_error("relative nominal iterations apply only to variation iteration settings");
        TFile reference(results.at("nominal").c_str(), "READ");
        auto *saved = dynamic_cast<TObjString *>(reference.Get((key.substr(0, dot) + "/settings").c_str()));
        if (!saved) throw std::runtime_error("nominal iteration result missing for " + key);
        const std::string marker = field == "PairIterations" ? "pair_iterations=" : "photon_iterations=", text = saved->GetString().Data();
        const auto at = text.find(marker); if (at == std::string::npos) throw std::runtime_error("nominal iteration metadata missing");
        const int k = std::stoi(text.substr(at + marker.size())) + (value == "nominal-1" ? -1 : 1);
        if (k < 1) throw std::runtime_error("iteration variation is outside positive iterations");
        value = std::to_string(k);
      }
      fs::path packages = nominalPackages;
      if (trees && (name == "nominal" || differentReduction)) {
        const auto cfg = output / "config" / (name + "_reduction.cfg"); config(cfg, fill);
        packages = output / ("packages_" + name);
        step(name + "_reduction", {(tth / "tree_to_hists.sh").string(), "local", "--settings", cfg.string(), "--inputs", input.string(), "--work", (output / ("work_" + name)).string(), "--out", packages.string(), "--parallel", parallel});
        if (name == "nominal") nominalPackages = packages;
      }
      if (name == "nominal" || differentReduction) {
        for (const auto &e : fs::recursive_directory_iterator(packages)) if (e.is_regular_file() && e.path().extension() == ".root") pin(e.path());
      }
      const auto result = output / "results" / (name + ".root"), qa = output / "results" / (name + "_qa.root"), cfg = output / "config" / (name + "_analysis.cfg");
      analysis["Input.Directory"] = packages.string(); analysis["Output.Results"] = result.string(); analysis["Output.QA"] = qa.string();
      config(cfg, analysis); results[name] = result;
      step(name + "_qa", {(package / "final_analysis.sh").string(), "qa", "--config", cfg.string()});
      step(name + "_analysis", {(package / "final_analysis.sh").string(), "run", "--config", cfg.string()});
    }
    fs::path systematicResult;
    if (!sys.empty()) {
      sys["Systematics.Nominal"] = results.at("nominal").string(); systematicResult = output / "results/systematics.root"; sys["Systematics.Output"] = systematicResult.string();
      auto bind = [&](std::string value) {
        if (value.empty()) throw std::runtime_error("empty systematic result binding");
        if (value[0] == '@') {
          const auto name = value.substr(1); if (!unique.count(name)) throw std::runtime_error("systematic names an undeclared variation: " + name);
          return results.at(name).string();
        }
        const auto external = resolve(resolve(base, systematic).parent_path(), value); pin(external); return external.string();
      };
      for (auto &[key, value] : sys) {
        const auto suffix = key.substr(key.rfind('.') + 1);
        if (suffix == "Up" || suffix == "Down" || suffix == "Alternative") value = bind(value);
        if (suffix == "Alternatives") {
          std::string rebuilt;
          for (const auto &item : words(value)) { const auto eq = item.find('='); if (eq == std::string::npos) throw std::runtime_error("systematic alternative needs label=result"); rebuilt += (rebuilt.empty() ? "" : " ") + item.substr(0, eq + 1) + bind(item.substr(eq + 1)); }
          value = rebuilt;
        }
      }
      const auto cfg = output / "config/systematics.cfg"; config(cfg, sys);
      step("systematics", {(package / "final_analysis.sh").string(), "systematics", "--config", cfg.string()});
    }
    std::vector<std::string> draw = {(package / "final_analysis.sh").string(), "draw", "--results", results.at("nominal").string(), "--qa", (output / "results/nominal_qa.root").string(), "--fig", (output / "fig").string()};
    std::vector<std::string> tables = {(package / "final_analysis.sh").string(), "tables", "--results", results.at("nominal").string(), "--out", (output / "tables").string()};
    if (!systematicResult.empty()) { draw.insert(draw.end(), {"--systematics", systematicResult.string()}); tables.insert(tables.end(), {"--systematics", systematicResult.string()}); }
    step("draw", draw); step("tables", tables);
    record["figures"] = Json::array(); record["unavailable_figures"] = Json::array();
    std::ifstream list(output / "fig/FIGURES.txt"); if (!list) throw std::runtime_error("missing figure coverage record");
    for (std::string l; std::getline(list, l);) {
      if (l.rfind("fig/", 0) == 0) { if (!fs::is_regular_file(output / l) || fs::file_size(output / l) == 0) throw std::runtime_error("missing saved figure " + l); record["figures"].push_back(l); }
      if (l.rfind("# not drawn:", 0) == 0) record["unavailable_figures"].push_back(l.substr(13));
    }
    if (record["figures"].empty()) throw std::runtime_error("no figures generated");
    record["missing_expected_figures"] = Json::array();
    for (const auto &f : expectedFiles) if (!fs::is_regular_file(output / "fig" / f) || fs::file_size(output / "fig" / f) == 0) record["missing_expected_figures"].push_back(f);
    if (!record["missing_expected_figures"].empty()) throw std::runtime_error("required figure coverage incomplete");
    for (const auto &[p, original] : record["input_pins"].items()) if (hash(p, output / "logs") != original.get<std::string>()) throw std::runtime_error("input changed during run: " + p);
    record["output_pins"] = Json::object();
    for (const auto &dir : {"results", "fig", "tables", "config"}) for (const auto &e : fs::recursive_directory_iterator(output / dir)) if (e.is_regular_file()) record["output_pins"][fs::relative(e.path(), output).string()] = hash(e.path(), output / "logs");
    record["status"] = expected == "none" ? "COMPONENTS_PASSED_COVERAGE_UNQUALIFIED" : "PASS_EXPECTED_OUTPUTS";
    atomic(output / "PIPELINE.json", record); return 0;
  } catch (const std::exception &e) {
    record["status"] = "FAILED"; record["first_error"] = e.what(); atomic(output / "PIPELINE.json", record); throw;
  }
}
}  /* namespace fa */
