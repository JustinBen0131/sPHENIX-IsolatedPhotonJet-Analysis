/* Shared interfaces for planning, fills, package merging and correction preparation. */
#pragma once

#include <nlohmann/json.hpp>

#include <chrono>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <ostream>
#include <TTree.h>

namespace treetohists {
using Json = nlohmann::json;
std::string sha256File(const std::string &path);     /* TreeToHists.cc */
std::string digest(const Json &value);
std::string digestStream(const std::function<void(std::ostream &)> &write);
Json readJson(const std::string &path);
Json reductionCode();
void checkReductionCode(const Json &plan);
void checkVertexFit(const Json &sample);
void checkPpPeriodFits(const Json &periods);
Json ppPeriodsOf(const Json &plan, const Json &sample);
Json inventory(const std::string &system, const std::string &family);
void merge(const Json &plan, const Json &sample, const std::string &work, const std::string &output);
Json preparePlan(const std::string &settings, const std::string &inputs, const std::string &home,
                 const std::vector<std::string> &samples, long long maxEvents, const std::string &work = "",
                 bool bounded = false, const Json &inspections = Json::array());
Json inspectSources(const Json &descriptor, const Json &snapshot);
void savePlan(const Json &plan, const std::string &work, bool writeRows = true);
void savePlanRow(const Json &plan, long index, const std::string &work,
                 const Json &source = Json());
Json freezeInputs(const std::string &inputs);
Json settingsSnapshot(const std::string &settings, const std::string &home);
void fillRow(const std::string &work, long index);
void fillGroup(const std::string &work, const std::string &group, const std::string &scratch);
std::string jobDirectory(const std::string &work, const Json &job);
void validatePackage(const std::string &path, const Json &plan, const Json &sample);
void publishValidatedPackage(const std::string &partial, const std::string &path,
                             const Json &plan, const Json &sample);
std::string packageBinding(const Json &plan, const Json &sample);
std::string packageBindingFromSources(const Json &plan, const Json &sample,
                                     const std::function<void(std::ostream &)> &writeSources);
Json sampleHeader(const Json &sample);  /* configuration without a copied source inventory */
void reusePackages(Json &plan, const std::string &collection);
void combinePackages(const Json &plan, const Json &sample, const std::vector<std::string> &inputs,
                     const std::string &output, const Json &inputSamples = Json::array());
int campaign(const std::string &command, const std::vector<std::string> &args);

/* Merge and validation own this counter. Nested record/hash loops refresh
   the active operation without counting another source or starting a poller.
   Unset TREETOHISTS_PROGRESS_PATH keeps ordinary interactive use quiet. */
class ReductionProgress {
 public:
  ReductionProgress(const std::string &operation, long long total);
  ~ReductionProgress() noexcept;
  void advance(long long processed);
  void finish();
  static void pulse(const std::string &phase);
 private:
  void emit(bool force);
  std::string path_, operation_, phase_;
  long long total_=0, processed_=0, emitted_=0;
  bool finished_=false;
  int exceptions_=0;
  std::chrono::steady_clock::time_point start_, last_;
  static thread_local ReductionProgress *active_;
};

/* Stream record tables through both raw and hierarchical merges. A package
   merge supplies the old-to-new source map; event references shift together. */
struct RecordTable {
  std::unique_ptr<TTree> output;
  std::map<std::string, Long64_t> ints;
  std::map<std::string, Double_t> doubles;
  void append(TTree &input, const std::string &name, Long64_t source, Long64_t offset, TDirectory *dest,
              const std::vector<Long64_t> *sourceMap = nullptr);
  void write();
};
int isolation(const std::vector<std::string> &args);  /* Isolation.cc: tree_to_hists isolation ... */
int vertex(const std::vector<std::string> &args);     /* Corrections.cc: AuAu truth-vertex population and fit */
}  /* namespace treetohists */

namespace treetohists { Json treeCopySummary(const Json &sample); }
