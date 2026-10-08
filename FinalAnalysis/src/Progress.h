/* Optional local analysis progress. Counts completed toy attempts and analysis
   phases without touching the random stream or scientific products. Set
   FINALANALYSIS_PROGRESS_PATH to an absolute JSON path to enable it. */
#pragma once

#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <nlohmann/json.hpp>

namespace fa {

class AnalysisProgress {
 public:
  AnalysisProgress(const std::string &campaign, long long total) : total_(total) {
    const char *path = std::getenv("FINALANALYSIS_PROGRESS_PATH");
    if (!path || !*path) return;
    path_ = path;
    if (total <= 0 || !std::filesystem::path(path_).is_absolute())
      throw std::runtime_error("invalid local analysis progress binding");
    campaign_ = environment("FINALANALYSIS_CAMPAIGN", campaign);
    task_ = environment("FINALANALYSIS_TASK_ID", "local");
    workstream_ = environment("FINALANALYSIS_WORKSTREAM", "local_analysis");
    row_ = environment("FINALANALYSIS_ROW_ID", campaign);
    start_ = last_ = std::chrono::steady_clock::now();
    emit(true);
  }
  ~AnalysisProgress() noexcept {
    if (path_.empty() || finished_) return;
    phase_ = "failed";
    error_ = "native analysis did not complete; see its error log";
    try { emit(true); } catch (const std::exception &e) { std::cerr << "[FinalAnalysis] progress: " << e.what() << '\n'; }
  }
  void phase(const std::string &phase, long long total) {
    if (path_.empty()) return;
    if (total < 0 || rowProcessed_ != rowTotal_) throw std::runtime_error("incomplete analysis progress phase");
    phase_ = phase; rowTotal_ = total; rowProcessed_ = 0;
    emit(true);
  }
  void advance(long long count = 1, bool skipped = false) {
    if (path_.empty()) return;
    if (count < 0 || count > rowTotal_ - rowProcessed_ || count > total_ - processed_)
      throw std::runtime_error("analysis progress exceeds its work bound");
    rowProcessed_ += count; processed_ += count;
    if (skipped) skipped_ += count;
    emit(rowProcessed_ == rowTotal_);
  }
  void finish() {
    if (path_.empty()) return;
    if (processed_ != total_ || rowProcessed_ != rowTotal_) throw std::runtime_error("incomplete analysis progress total");
    phase_ = "complete"; emit(true); finished_ = true;
  }

 private:
  static std::string environment(const char *key, const std::string &fallback) {
    const char *value = std::getenv(key); return value && *value ? value : fallback;
  }
  void emit(bool force) {
    const auto now = std::chrono::steady_clock::now();
    if (!force && processed_ - emitted_ < 25 && std::chrono::duration<double>(now - last_).count() < 2) return;
    const double elapsed = std::chrono::duration<double>(now - start_).count();
    const double rate = elapsed > 0 ? processed_ / elapsed : 0;
    using Json = nlohmann::json;
    Json record = {{"schema", "LongRunningCanaryProgressV1"}, {"producer", "FinalAnalysis"},
      {"task_id", task_}, {"workstream_id", workstream_}, {"campaign_tag", campaign_}, {"row_id", row_},
      {"unit", "toy attempts and analysis phases"}, {"phase", phase_},
      {"processed", processed_}, {"total", total_}, {"remaining", total_ - processed_},
      {"row_processed", rowProcessed_}, {"row_total", rowTotal_}, {"row_remaining", rowTotal_ - rowProcessed_},
      {"deterministic_skip_processed", skipped_}, {"deterministic_skip_total", skipped_}, {"deterministic_skip_remaining", 0},
      {"elapsed_seconds", elapsed}, {"rate", rate}, {"eta_seconds", rate > 0 ? (total_ - processed_) / rate : -1},
      {"updated_at_unix", std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},
      {"latest_error", error_.empty() ? Json() : Json(error_)}};
    const std::filesystem::path temporary = path_ + ".next";
    std::ofstream out(temporary); out << record.dump() << '\n'; out.close();
    if (!out) throw std::runtime_error("cannot write local analysis progress");
    std::filesystem::rename(temporary, path_);
    last_ = now; emitted_ = processed_;
  }
  std::string path_, campaign_, task_, workstream_, row_, phase_ = "starting", error_;
  long long total_ = 0, processed_ = 0, emitted_ = 0, rowTotal_ = 0, rowProcessed_ = 0, skipped_ = 0;
  bool finished_ = false;
  std::chrono::steady_clock::time_point start_, last_;
};

}  // namespace fa
