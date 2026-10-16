/* Per-output scientific identity. ROOT rendering remains in IANFigures/QAFigures.
   No policy is inferred from canvas text or source-file hashes. Incomplete is the
   default until the drawing function has registered every covered category. */
#ifndef PHOTONJET_PLOT_PHYSICS_H
#define PHOTONJET_PLOT_PHYSICS_H
#include "RenderMetadata.h"
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace photonjet_render {
inline std::string physicsDouble(double value) {
  static_assert(sizeof(double) == sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559,
                "PlotPhysicsV1 requires IEEE754 binary64");
  if (!std::isfinite(value)) throw std::runtime_error("nonfinite physics value needs a validity mask");
  if (value == 0) value = 0;  // signed zero has the same scientific value
  std::uint64_t bits;
  std::memcpy(&bits, &value, sizeof bits);
  std::ostringstream out;
  out << "f64:" << std::hex << std::setw(16) << std::setfill('0') << bits;
  return out.str();
}

inline void requirePortablePhysics(const Json &value) {
  if (value.is_number_float()) throw std::runtime_error("PlotPhysicsV1: encode real values with physicsDouble");
  if (value.is_structured()) for (const auto &child : value) requirePortablePhysics(child);
}

/* Construct afresh for every independent output, including each supplemental
   output. The key of every quantity is its stable scientific series identity,
   never its colour, legend wording, ROOT counter, filename or draw order. */
class PlotPhysics {
  Json definitions_ = Json::object(), quantities_ = Json::object(), dependencies_ = Json::object();
  Json coverage_ = Json::object();
  std::string family_;
public:
  explicit PlotPhysics(const std::string &family = "unregistered") : family_(family) {
    for (const char *name : {"quantities", "errors", "covariance", "scalar_annotations", "selection",
                            "weights", "normalization", "models", "projection", "dependencies"})
      coverage_[name] = {{"state", "unknown"}, {"evidence", "renderer not yet instrumented"}};
  }
  void definition(const std::string &id, const Json &value) { add(definitions_, id, value); }
  void quantity(const std::string &id, const Json &value) { add(quantities_, id, value); }
  void dependency(const std::string &id, const Json &value) { add(dependencies_, id, value); }
  void coverage(const std::string &name, const std::string &state, const std::string &evidence) {
    if (!coverage_.contains(name) || evidence.empty() ||
        (state != "complete" && state != "not_applicable" && state != "unknown"))
      throw std::runtime_error("invalid plot physics coverage");
    coverage_[name] = {{"state", state}, {"evidence", evidence}};
  }
  Json finish(const std::string &output) const {
    if (output.empty()) throw std::runtime_error("plot output identity is required");
    Json body = {{"schema", "PlotPhysicsV1"}, {"output", output}, {"renderer_family", family_},
                 {"contract_version", 1}, {"definitions", definitions_},
                 {"quantities", quantities_}, {"dependencies", dependencies_}, {"coverage", coverage_}};
    requirePortablePhysics(body);
    const std::string text = body.dump(-1, ' ', true);
    return {{"canonical", text}, {"sha256", sha256(text)}};
  }
private:
  static void add(Json &to, const std::string &id, const Json &value) {
    if (id.empty() || to.contains(id)) throw std::runtime_error("duplicate/empty scientific identity: " + id);
    requirePortablePhysics(value);
    to[id] = value;
  }
};

/* Adapter for the existing Bins values; scientific arrays before cosmetic
   line smoothing/marker offsets. This neither computes nor redraws any bin.
   The renderer must explicitly supply its existing omitted-empty-bin mask. */
inline Json physicsBins(const std::vector<double> &edges, const std::vector<double> &values,
                        const std::vector<double> &variances, const std::vector<bool> &valid) {
  if (edges.size() != values.size()+1 || variances.size() != values.size() || valid.size() != values.size())
    throw std::runtime_error("physics bin shape mismatch");
  Json x = Json::array(), y = Json::array(), e2 = Json::array(), mask = Json::array();
  for (double v : edges) x.push_back(physicsDouble(v));
  for (size_t i = 0; i < values.size(); ++i) {
    if (variances[i] < 0 || !std::isfinite(variances[i])) throw std::runtime_error("invalid variance");
    mask.push_back(valid[i]);
    y.push_back(valid[i] ? Json(physicsDouble(values[i])) : Json(nullptr));
    e2.push_back(physicsDouble(variances[i]));
  }
  return {{"edges", x}, {"values", y}, {"variances", e2}, {"valid", mask}};
}
}  // namespace photonjet_render
#endif
