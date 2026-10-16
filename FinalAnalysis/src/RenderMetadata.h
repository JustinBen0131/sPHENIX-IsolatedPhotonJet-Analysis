/* Shared metadata checks for C++ analysis and ROOT rendering.
   seal/verify bind the JSON description; product receipts separately bind the actual ROOT-file bytes.
   Renderers use both checks before labeling a figure with the saved selection and binning. */
#ifndef PHOTONJET_RENDER_METADATA_H
#define PHOTONJET_RENDER_METADATA_H
#include <nlohmann/json.hpp>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <system_error>
#include <unistd.h>

namespace photonjet_render {
using Json = nlohmann::json;
/* Sample origin and selected truth class are different facts. These labels
   describe the candidates plotted; the source family stays in the receipt. */
enum class CandidatePopulation { TruthTaggedBackground, AllInclusiveJetCandidates, OtherInclusiveJetCandidates };
inline Json candidatePopulation(CandidatePopulation population) {
  switch (population) {
    case CandidatePopulation::TruthTaggedBackground:
      return {{"population", "truth_tagged_background"}, {"label", "Background-tagged"}};
    case CandidatePopulation::AllInclusiveJetCandidates:
      return {{"population", "all_inclusive_jet_candidates"}, {"label", "Inclusive+jet"}};
    case CandidatePopulation::OtherInclusiveJetCandidates:
      return {{"population", "other_inclusive_jet_candidates"}, {"label", "Other inclusive+jet candidates"}};
  }
  throw std::runtime_error("unknown candidate population");
}
inline std::string candidateLegendLabel(const Json &population) {
  const std::string base = population.at("label");
  if (!population.contains("source_slices_gev")) return base;
  const auto &slices = population.at("source_slices_gev");
  if (!slices.is_array() || slices.empty()) throw std::runtime_error("candidate source slices are required");
  std::string ranges;
  int previous = 0;
  for (const auto &slice : slices) {
    if (!slice.is_number_integer() || slice.get<int>() <= previous)
      throw std::runtime_error("candidate source slices must be positive and ordered");
    previous = slice.get<int>();
    ranges += (ranges.empty() ? "" : "/") + std::to_string(previous);
  }
  const std::string prefix = population.at("population") == "all_inclusive_jet_candidates" ? "Inclusive" : base;
  return prefix + ": " + population.at("sample_family").get<std::string>() + " " + ranges;
}
inline void requireCandidateLabel(const Json &population, const std::string &label) {
  const auto kind = population.at("population").get<std::string>();
  const auto expected = candidatePopulation(kind == "truth_tagged_background"
      ? CandidatePopulation::TruthTaggedBackground
      : kind == "all_inclusive_jet_candidates" ? CandidatePopulation::AllInclusiveJetCandidates
      : kind == "other_inclusive_jet_candidates" ? CandidatePopulation::OtherInclusiveJetCandidates
      : throw std::runtime_error("unknown candidate population: " + kind));
  if (population.value("sample_family", "").empty())
    throw std::runtime_error("candidate source family is required");
  if (kind != "truth_tagged_background" && population.at("sample_family") != "inclusive+jet")
    throw std::runtime_error("inclusive+jet candidates require the inclusive+jet sample family");
  if (population.at("label") != expected.at("label") || label != candidateLegendLabel(population))
    throw std::runtime_error("candidate legend differs from plotted population: " + label);
}
/* Physics distributions use the producer event weight and nominal source
   normalization. The producer records the complete weight-array digest and
   each source's actual normalization; unit display normalization changes none
   of this ancestry. Training weights and absent/unknown origins never render. */
inline void requirePhysicsPlotWeights(const Json &weights, const std::string &plot) {
  const auto fail = [&]() { throw std::runtime_error("physics plot weight provenance is missing, training, or unknown: " + plot); };
  const auto hash = [](const Json &value) {
    if (!value.is_string()) return false;
    const auto text = value.get<std::string>();
    return text.size() == 64 && text.find_first_not_of("0123456789abcdef") == std::string::npos;
  };
  // Native histogram products already contain their producer event weights.
  // Bind the exact products and the normalization applied by FinalAnalysis;
  // a training diagnostic has neither this package schema nor this contract.
  if (weights.is_object() && weights.value("schema", "") == "PhotonJetHistogramPlotWeightsV1") {
    if (weights.value("role", "") != "physics_event_source" ||
        !weights.contains("training_weights_used") || weights.at("training_weights_used") != false ||
        !weights.contains("packages") || !weights.at("packages").is_array() || weights.at("packages").empty()) fail();
    std::set<std::string> packages;
    for (const auto &p : weights.at("packages")) {
      if (p.value("schema", "") != "PhotonJetSampleV1" || p.value("source_count", 0LL) <= 0 ||
          !hash(p.at("sha256")) || !hash(p.at("selection_sha256")) || !hash(p.at("producer_sha256")) ||
          !packages.insert(p.at("name").get<std::string>()).second) fail();
      const auto family = p.value("family", "");
      const auto &n = p.at("normalization");
      if (family == "data") { if (!n.is_null()) fail(); continue; }
      if (family != "photonjet" && family != "inclusivejet") fail();
      const auto rule = n.value("source_factor", "");
      /* Embedded Au+Au stitching records its numerator as the sample's cross section. */
      const double factor = n.value("pb_per_generated_event", -1.0),
                   numerator = n.value("numerator_pb", rule == "auau_embedded_stitching" ? n.value("cross_section_pb", -1.0) : -1.0),
                   denominator = n.value("denominator_events", -1.0);
      if (!std::isfinite(factor) || !(factor > 0) || !std::isfinite(numerator) || !(numerator > 0) ||
          !std::isfinite(denominator) || !(denominator > 0) ||
          n.value("applied_in_histograms", true) ||
          std::abs(factor - numerator / denominator) > 1e-12 * factor ||
          (rule != "family_reference" && rule != "cross_section_over_generated" && rule != "auau_embedded_stitching")) fail();
    }
    return;
  }
  /* The native score collector records the actual event-weighted sources and
     applies its declared cross section / ownership denominator before writing. */
  if (weights.is_object() && weights.value("schema", "") == "PhotonIDWorkingPointReviewV1") {
    if (weights.value("status", "") != "PASS" || weights.value("weighting", "") != "physics_event_times_sample_normalization" ||
        !weights.contains("training_weights_used") || !weights.at("training_weights_used").is_boolean() ||
        weights.at("training_weights_used").get<bool>() || !weights.contains("sample_normalization") ||
        !weights.at("sample_normalization").is_object() || weights.at("sample_normalization").empty() ||
        !weights.contains("sources") || !weights.at("sources").is_array() || weights.at("sources").empty()) fail();
    std::set<std::string> seen;
    for (const auto &source : weights.at("sources")) {
      if (!source.is_object() || !source.contains("sha256") || !hash(source.at("sha256")) ||
          !source.contains("logical_id") || !hash(source.at("logical_id")) ||
          !seen.insert(source.at("logical_id").get<std::string>()).second) fail();
    }
    for (const auto &sample : weights.at("sample_normalization")) {
      const double cross = sample.value("cross_section_pb", -1.0), denominator = sample.value("denominator_events", -1.0),
                   factor = sample.value("pb_per_event", -1.0);
      if (!std::isfinite(cross) || cross <= 0 || !std::isfinite(denominator) || denominator <= 0 ||
          std::floor(denominator) != denominator || !std::isfinite(factor) || factor <= 0 ||
          std::abs(factor - cross / denominator) > 1e-12 * std::max(factor, cross / denominator) ||
          !sample.contains("definition") || !sample.at("definition").is_object()) fail();
      const auto &definition = sample.at("definition");
      const auto family = definition.value("family", ""), owner = definition.value("owner", ""), basis = definition.value("denominator", "");
      if (definition.value("cross_section_pb", -1.0) != cross ||
          !((family == "photonjet" && owner == "generator_filter_photon_pt" && basis == "all_events") ||
            (family == "inclusivejet" && owner == "leading_truth_jet_r04_pt" && basis == "events_in_window"))) fail();
    }
    return;
  }
  if (!weights.is_object() || weights.value("schema", "") != "PhotonJetPlotWeightProvenanceV1" ||
      weights.value("role", "") != "physics_event_source") fail();
  const auto field = weights.value("event_weight_field", "");
  if ((field != "event_weight" && field != "w_event" && field != "m_mcEventWeight") ||
      !weights.contains("complete_weights_sha256") || !hash(weights.at("complete_weights_sha256")) ||
      !weights.contains("sources") || !weights.at("sources").is_array() || weights.at("sources").empty()) fail();
  std::set<std::string> seen;
  for (const auto &source : weights.at("sources")) {
    if (!source.is_object() || !source.contains("source_sha256") || !hash(source.at("source_sha256")) ||
        !source.contains("normalization_sha256") || !hash(source.at("normalization_sha256")) ||
        !source.contains("normalization_factor") || !source.at("normalization_factor").is_number()) fail();
    const auto rule = source.value("source_factor", ""), id = source.value("source_id", "");
    const double factor = source.at("normalization_factor");
    if (id.empty() || !seen.insert(id).second || !std::isfinite(factor) || factor <= 0 ||
        (rule != "data_event_weight" && rule != "cross_section_over_generated" &&
         rule != "family_reference" && rule != "auau_embedded_stitching")) fail();
  }
}
/* Training preparation and learned-feature summaries have a separate, narrow
   contract. Matched learner objectives and unweighted paired-model mappings
   are diagnostics; this never authorizes a physical score spectrum or ROC. */
inline void requireModelTrainingDiagnostic(const Json &proof, const std::string &plot) {
  const auto fail = [&]() { throw std::runtime_error("invalid training diagnostic provenance: " + plot); };
  if (!proof.is_object() || proof.value("schema", "") != "PhotonIDTrainingDiagnosticV1" ||
      proof.value("role", "") != "model_training_diagnostic" ||
      !proof.contains("physics_distribution") || proof.at("physics_distribution") != false) fail();
  const auto name = plot.substr(plot.find_last_of('/') + 1), kind = proof.value("diagnostic", "");
  const bool matched=kind=="model_routing" || kind=="training_capacity";
  const bool pairs=kind=="paired_score_mapping" || kind=="paired_score_difference";
  const std::string weighting = pairs ? "unit_counts" : (kind == "feature_separation" || matched ? "training_weighted" :
    (kind == "split_gain" ? "learned_model_split_gain" : "raw_and_training_weighted"));
  if (proof.value("weighting", "") != weighting) fail();
  const bool allowed = (kind == "training_reweighting" &&
    (name == "auau_training_et_reweighting.png" || name == "auau_training_eta_reweighting.png" ||
     name == "pp_training_et_reweighting.png" || name == "pp_training_eta_reweighting.png")) ||
    (kind == "feature_separation" && (name == "auau_bdt_single_feature_separation.png" || name == "pp_bdt_single_feature_separation.png")) ||
    (kind == "split_gain" && (name == "auau_bdt_feature_split_gain.png" || name == "pp_bdt_feature_split_gain.png")) ||
    (kind == "model_routing" && name == "auau_bdt_routing_auc_by_et_centrality.png") ||
    (kind == "training_capacity" && name == "auau_bdt_tree_depth_train_holdout_diagnostic.png") ||
    (kind == "paired_score_mapping" && name == "pp_bdt_original_vs_retrained_scores.png") ||
    (kind == "paired_score_difference" && name == "pp_bdt_retraining_score_difference.png");
  if (!allowed || !proof.contains("models") || !proof.at("models").is_array() || proof.at("models").empty()) fail();
  std::set<std::string> modelNames;
  for (const auto &model : proof.at("models")) {
    if (!model.is_object() || model.value("name", "").empty() || !modelNames.insert(model.at("name").get<std::string>()).second) fail();
    for (const char *key : {"population_sha256", "model_artifacts_sha256", "training_weights_sha256"}) {
      const auto value = model.value(key, "");
      if (value.size() != 64 || value.find_first_not_of("0123456789abcdef") != std::string::npos) fail();
    }
    if(matched)for(const char *key:{"population_sha256","training_weights_sha256","event_split_sha256","split_sha256"}) {
      const auto value=model.value(key,"");
      if(value.size()!=64 || value.find_first_not_of("0123456789abcdef")!=std::string::npos ||
         value!=proof.at("models").at(0).value(key,""))fail();
    }
  }
  if(matched && proof.at("models").size()<2)fail();
  if(pairs) {
    const auto reference=proof.value("reference_model_sha256","");
    if(reference.size()!=64 || reference.find_first_not_of("0123456789abcdef")!=std::string::npos ||
       !proof.contains("paired_rows") || !proof.at("paired_rows").is_number_integer() || proof.at("paired_rows").get<long long>()<=0)fail();
  }
}
/* The SHA-256 of a text, or of a file's bytes. Exact argv: neither file names nor scientific metadata are evaluated
   by a shell.
   A file's hash is kept with the other saved work ($FINALANALYSIS_CACHE/sha256/, set by final_analysis.sh), so a
   package of gigabytes is read once and not at the start and end of every step that binds it. The record is believed
   while the file is the same file: its path, inode, size, the time it was last written and the time its status last
   changed. No program can set the last one back, so a file that was rewritten or replaced and given its old time
   again (cp -p, a time-keeping copy, a rename back) is still read again. */
inline std::string sha256(const std::string &value, bool file = false) {
  const auto identity = [](const std::string &path) {
    struct stat s;
    if (::stat(path.c_str(), &s) != 0 || !S_ISREG(s.st_mode)) return std::string();
#ifdef __APPLE__
    const auto &written = s.st_mtimespec, &changed = s.st_ctimespec;
#else
    const auto &written = s.st_mtim, &changed = s.st_ctim;
#endif
    return path + "\n" + std::to_string(s.st_ino) + " " + std::to_string(s.st_size) + " " +
           std::to_string(written.tv_sec) + " " + std::to_string(written.tv_nsec) + " " +
           std::to_string(changed.tv_sec) + " " + std::to_string(changed.tv_nsec) + "\n";
  };
  std::string path, record, stamp;
  const char *cache = std::getenv("FINALANALYSIS_CACHE");
  if (file && cache && *cache) {
    path = std::filesystem::absolute(value).string();
    stamp = identity(path);
    record = std::string(cache) + "/sha256/" + std::to_string(std::hash<std::string>{}(path));
    std::ifstream saved(record);
    const std::string text((std::istreambuf_iterator<char>(saved)), std::istreambuf_iterator<char>());
    if (!stamp.empty() && text.size() == stamp.size() + 64 && text.compare(0, stamp.size(), stamp) == 0 &&
        text.find_first_not_of("0123456789abcdef", stamp.size()) == std::string::npos)
      return text.substr(stamp.size());
  }
  int input[2], output[2];
  if (pipe(input) || pipe(output)) throw std::runtime_error("SHA256 pipe failed");
  const pid_t pid = fork();
  if (pid == 0) {
    dup2(input[0], STDIN_FILENO); dup2(output[1], STDOUT_FILENO);
    close(input[0]); close(input[1]); close(output[0]); close(output[1]);
    if (file) { execlp("sha256sum", "sha256sum", "--", value.c_str(), nullptr); execlp("shasum", "shasum", "-a", "256", "--", value.c_str(), nullptr); }
    else { execlp("sha256sum", "sha256sum", nullptr); execlp("shasum", "shasum", "-a", "256", nullptr); }
    _exit(127);
  }
  close(input[0]); close(output[1]);
  if (pid < 0) { close(input[1]); close(output[0]); throw std::runtime_error("SHA256 child failed"); }
  size_t offset = 0;
  while (!file && offset < value.size()) {
    const auto count = ::write(input[1], value.data()+offset, value.size()-offset);
    if (count < 0 && errno == EINTR) continue;
    if (count <= 0) break;
    offset += count;
  }
  close(input[1]); std::string text; char buffer[1024]; ssize_t count;
  while ((count = ::read(output[0], buffer, sizeof buffer)) != 0) {
    if (count < 0 && errno == EINTR) continue;
    if (count < 0) break;
    text.append(buffer, count);
  }
  close(output[0]); int status = 0;
  while (waitpid(pid, &status, 0) < 0) if (errno != EINTR) throw std::runtime_error("SHA256 child not reaped");
  const auto hash = text.substr(0, text.find_first_of(" \t\n"));
  if (!WIFEXITED(status) || WEXITSTATUS(status) || (!file && offset != value.size()) || hash.size()!=64 ||
      hash.find_first_not_of("0123456789abcdef")!=std::string::npos) throw std::runtime_error("SHA256 failed");
  /* Saved only when the file is, after the read, exactly what it was before it; a file that changed while it was read
     is read again next time. Saving is a convenience: a folder that cannot be written never fails the hash. */
  if (!stamp.empty() && identity(path) == stamp) {
    std::error_code ignored;
    std::filesystem::create_directories(std::filesystem::path(record).parent_path(), ignored);
    const std::string partial = record + ".partial" + std::to_string(getpid());
    std::ofstream out(partial);
    out << stamp << hash;
    out.close();
    if (out) std::filesystem::rename(partial, record, ignored);
    if (!out || ignored) std::filesystem::remove(partial, ignored);
  }
  return hash;
}
/* The hash TreeToHists gives a resolved selection: its content without the identifiers derived from it. */
inline std::string selectionContentHash(Json selection) {
  for (const char *derived : {"selection_sha256", "profile", "name"}) selection.erase(derived);
  return sha256(selection.dump(-1, ' ', true));
}
/* Whether a sample's saved selection agrees with the selection it is analysed under. A resolved selection names the
   p+p DATA trigger, which only p+p DATA evaluates: any other sample may carry another p+p DATA trigger name, with the
   identifiers derived from it (selection_sha256, profile). Every other setting must be identical, and both selections
   must carry the hash of their own content, so a relabelled selection never passes. */
inline bool selectionsAgree(const Json &sample, const Json &reference, bool ppData) {
  if (sample == reference) return true;
  if (ppData) return false;
  const auto evaluated = [](Json s) {
    for (const char *derived : {"selection_sha256", "profile"}) s.erase(derived);
    if (s.contains("systems") && s.at("systems").contains("pp") && s.at("systems").at("pp").contains("event_selection"))
      s["systems"]["pp"]["event_selection"].erase("trigger_name");
    return s;
  };
  if (evaluated(sample) != evaluated(reference)) return false;
  for (const Json *s : {&sample, &reference}) {
    const std::string id = selectionContentHash(*s);
    if (s->value("selection_sha256", "") != id || s->value("profile", "") != "selection:" + id) return false;
  }
  return true;
}
inline Json seal(const Json &body) { return {{"body", body}, {"sha256", sha256(body.dump())}}; }
inline Json verify(const Json &sealed) {
  if (sha256(sealed.at("body").dump()) != sealed.at("sha256").get<std::string>()) throw std::runtime_error("render metadata hash differs");
  return sealed.at("body");
}
inline void write(const std::string &path, const Json &body) {
  std::ofstream out(path + ".partial"); out << body.dump(2) << '\n'; out.close();
  if (!out || std::rename((path + ".partial").c_str(), path.c_str())) throw std::runtime_error("cannot write render receipt " + path);
}
inline Json read(const std::string &path) {
  std::ifstream in(path); if (!in) throw std::runtime_error("cannot read render receipt " + path);
  Json value; in >> value; return value;
}
}  /* namespace photonjet_render */
#endif
