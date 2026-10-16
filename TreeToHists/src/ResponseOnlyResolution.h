#ifndef PHOTONJET_RESPONSE_ONLY_RESOLUTION_H
#define PHOTONJET_RESPONSE_ONLY_RESOLUTION_H

#include "SelectionReplay.h"
#include <set>

/* Shared publisher/consumer declaration for a response-only nuisance.
   This admits local mechanics only, never numerical calibration adoption. */
namespace photonjet_response_only {
using Json = nlohmann::json;
inline void require(bool ok, const char* why) {
  if (!ok) throw std::runtime_error(std::string("response-only resolution: ") + why);
}
inline bool validate(const Json& recipe, const std::string& sample, const std::string& family) {
  const auto& effects = recipe.at("sample_effects");
  if (!recipe.contains("response_only_resolution")) {
    for (const auto& effect : effects)
      require(effect.get<std::string>().find("photon_selection=") == std::string::npos,
              "selection policy lacks its scope declaration");
    return false;
  }
  const auto& r = recipe.at("response_only_resolution");
  require(!recipe.contains("selection_change"), "combined cut and resolution variation is not qualified");
  require(recipe.at("admission") == "MECHANICAL_FIXTURE_ONLY" &&
          recipe.at("cone_prescription") == "unchanged-nominal-measurement-and-isolation",
          "unqualified admission or cone policy");
  require(r.at("schema") == "PhotonJetResponseOnlyEMResolutionV1" &&
          r.at("measurement") == "nominal" && r.at("classification") == "original_unsmeared" &&
          r.at("purity_leakage") == "original_unsmeared" &&
          r.at("response_kinematics") == "varied" &&
          r.at("leader") == "highest_varied_et_among_nominal_id_isolation",
          "measurement, classification or leader policy differs");
  const std::string mode = r.at("mode");
  require(mode == "ppg18_low" || mode == "ppg18_nominal" || mode == "ppg18_high", "unknown resolution mode");
  require(r.at("seed").is_number_integer(), "seed must be an exact integer");
  const auto seed = r.at("seed").get<long long>();
  require(seed >= 0 && seed <= 9007199254740991LL, "seed outside exact integer domain");
  require(r.at("affected_samples").is_array() && !r.at("affected_samples").empty(), "missing affected samples");
  std::set<std::string> affected;
  for (const auto& item : r.at("affected_samples")) {
    const std::string key = item;
    require(effects.contains(key) && affected.insert(key).second, "unknown or duplicated affected sample");
  }
  const std::string expected = "photon_resolution=" + mode + " photon_selection=response_only seed=" + std::to_string(seed);
  for (auto it = effects.begin(); it != effects.end(); ++it)
    require(it.value() == (affected.count(it.key()) ? expected : ""), "undeclared or combined sample effect");
  require(family != "data" || !affected.count(sample), "resolution cannot act on DATA");
  require(r.at("definition_sources").is_array() && !r.at("definition_sources").empty(), "missing definition source");
  for (const auto& source : r.at("definition_sources"))
    require(photonjet_selection_replay::digestString(source.at("sha256")) &&
            source.at("path").is_string() && !source.at("path").get<std::string>().empty() &&
            source.at("locator").is_string() && !source.at("locator").get<std::string>().empty(),
            "invalid definition identity");
  return affected.count(sample) != 0;
}
}
#endif
