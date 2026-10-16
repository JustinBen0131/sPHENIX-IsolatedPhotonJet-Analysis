/* Phases: how long each named phase of an analysis took. A phase that lasts a second or more is reported on standard
   output when it ends, so a log shows where the time went. */
#pragma once

#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <condition_variable>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <thread>
#include <nlohmann/json.hpp>

namespace fa {

class Phases {
 public:
  Phases() = default;
  /* Opt-in operational telemetry; one atomic record per collision-system calculation. It never changes a draw or
     random stream. Overall units are complete system calculations, phase units are explicitly named attempts. */
  Phases(const std::string &system, const std::string &analysis, size_t photonBins, const std::string &centrality)
      : system_(system), analysis_(analysis), centrality_(centrality), photonBins_(photonBins) {
    const char *directory = std::getenv("FINALANALYSIS_PROGRESS");
    if (!directory || !*directory) return;
    std::filesystem::create_directories(directory);
    path_ = std::string(directory) + "/" + system + ".json";
    publish();
    reporter_ = std::thread([this] {
      std::unique_lock<std::mutex> lock(mutex_);
      while (!changed_.wait_for(lock, std::chrono::seconds(1), [this] { return stop_; })) {
        try { publish(); } catch (const std::exception &e) { telemetryError_ = e.what(); }
      }
    });
  }
  ~Phases() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
      if (!path_.empty()) {
        if (!complete_ && error_.empty()) error_ = "system calculation ended without completion";
        try { publish(); } catch (...) {}
      }
    }
    changed_.notify_all();
    if (reporter_.joinable()) reporter_.join();
    finishDuration();
  }
  /* Ends the phase in progress, if any, and starts `name`. */
  void start(const std::string &name, size_t total = 1, const std::string &unit = "phase_completion") {
    std::lock_guard<std::mutex> lock(mutex_);
    check();
    finishDuration();
    name_ = name;
    since_ = std::chrono::steady_clock::now();
    phaseTotal_ = total; phaseProcessed_ = 0; unit_ = unit;
    if (!path_.empty()) publish();
  }
  void tick() {
    if (path_.empty()) return;
    std::lock_guard<std::mutex> lock(mutex_);
    check();
    if (phaseProcessed_ >= phaseTotal_) throw std::runtime_error("progress exceeds phase total");
    ++phaseProcessed_;
  }
  /* A scoped attempt counts even when a toy is rejected or an exception is handled by its existing caller. */
  struct Attempt { Phases &owner; ~Attempt() {
    try { owner.tick(); } catch (const std::exception &e) {
      std::lock_guard<std::mutex> lock(owner.mutex_); owner.telemetryError_ = e.what();
    }
  } };
  Attempt attempt() { return Attempt{*this}; }
  void complete() {
    std::lock_guard<std::mutex> lock(mutex_); check();
    phaseProcessed_ = phaseTotal_; complete_ = true;
    if (!path_.empty()) publish();
  }
  void fail(const std::string &error) {
    std::lock_guard<std::mutex> lock(mutex_); error_ = error;
    if (!path_.empty()) publish();
  }

 private:
  void finishDuration() {
    const auto now = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(now - since_).count();
    if (!name_.empty() && seconds >= 1)
      std::cout << "[FinalAnalysis] " << name_ << ": " << std::round(10 * seconds) / 10 << " s" << std::endl;
  }
  void check() const { if (!telemetryError_.empty()) throw std::runtime_error(telemetryError_); }
  void publish() const {
    using Json = nlohmann::json;
    const auto now = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(now - began_).count();
    const double phaseElapsed = std::chrono::duration<double>(now - since_).count();
    const bool measuredPhase = unit_ == "toy_attempt" || unit_ == "iteration_attempt";
    const Json rate = measuredPhase && phaseProcessed_ > 0 && phaseElapsed > 0 ? Json(phaseProcessed_ / phaseElapsed) : Json();
    const Json eta = complete_ ? Json(0) : !rate.is_null() ? Json((phaseTotal_ - phaseProcessed_) / rate.get<double>()) : Json();
    const std::time_t stamp = std::time(nullptr); std::tm utc{}; gmtime_r(&stamp, &utc);
    char updated[32]; std::strftime(updated, sizeof(updated), "%Y-%m-%dT%H:%M:%SZ", &utc);
    const auto context = [](const char *name) { const char *v = std::getenv(name); return v ? v : "UNBOUND"; };
    Json record = {{"schema", "FinalAnalysisProgressV1"}, {"task", context("FINALANALYSIS_TASK")},
      {"workstream", context("FINALANALYSIS_WORKSTREAM")}, {"campaign", analysis_}, {"row", system_},
      {"centrality", centrality_}, {"photon_bins", photonBins_}, {"phase", name_},
      {"work_unit", "collision_system_calculation"}, {"processed", complete_ ? 1 : 0}, {"total", 1}, {"remaining", complete_ ? 0 : 1},
      {"row_processed", complete_ ? 1 : 0}, {"row_total", 1}, {"row_remaining", complete_ ? 0 : 1},
      {"phase_work_unit", unit_}, {"phase_processed", phaseProcessed_}, {"phase_total", phaseTotal_},
      {"phase_remaining", phaseTotal_ - phaseProcessed_},
      {"deterministic_skip_processed", 0}, {"deterministic_skip_total", 0}, {"deterministic_skip_remaining", 0},
      {"rate", complete_ && elapsed > 0 ? Json(1 / elapsed) : Json()}, {"eta", complete_ ? Json(0) : Json()},
      {"phase_rate", rate}, {"phase_eta", eta}, {"updated_at", updated}, {"latest_error", error_},
      {"state", !error_.empty() ? "FAILED" : complete_ ? "COMPLETE" : "RUNNING"}};
    const std::string partial = path_ + ".partial";
    { std::ofstream out(partial); out << record.dump(2) << '\n'; if (!out) throw std::runtime_error("cannot write progress: " + partial); }
    std::filesystem::rename(partial, path_);
  }
  std::string name_;
  std::chrono::steady_clock::time_point since_;
  std::chrono::steady_clock::time_point began_ = std::chrono::steady_clock::now();
  std::string path_, system_, analysis_, centrality_, unit_ = "phase_completion", error_, telemetryError_;
  size_t photonBins_ = 0, phaseTotal_ = 1, phaseProcessed_ = 0;
  bool stop_ = false, complete_ = false;
  std::mutex mutex_;
  std::condition_variable changed_;
  std::thread reporter_;
};

}  /* namespace fa */
