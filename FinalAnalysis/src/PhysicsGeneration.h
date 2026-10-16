/* Read back the exact native inventory emitted by the draw command. */
#ifndef FINALANALYSIS_PHYSICS_GENERATION_H
#define FINALANALYSIS_PHYSICS_GENERATION_H
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

namespace fa {
template<class Hash>
nlohmann::json physicsGeneration(const std::filesystem::path& figures, Hash hash) {
  namespace fs = std::filesystem;
  using J = nlohmann::json;
  const auto base = fs::canonical(figures);
  const auto pointerPath = base / "PHYSICS_CURRENT.json";
  std::ifstream pointerFile(pointerPath);
  if (!pointerFile) throw std::runtime_error("missing native physics generation pointer");
  J pointer; pointerFile >> pointer;
  if (pointer.value("schema", "") != "PlotPhysicsGenerationPointerV1")
    throw std::runtime_error("invalid native physics generation pointer");
  auto path = fs::path(pointer.at("inventory").at("path").get<std::string>());
  if (!path.is_absolute()) path = base / path;
  path = fs::canonical(path);
  const auto relative = path.lexically_relative(base);
  if (relative.empty() || *relative.begin() == ".." || !fs::is_regular_file(path))
    throw std::runtime_error("native physics inventory escapes figure generation");
  const auto digest = hash(path);
  if (digest != pointer.at("inventory").at("sha256").get<std::string>())
    throw std::runtime_error("native physics inventory hash changed");
  std::ifstream inventoryFile(path);
  J inventory; inventoryFile >> inventory;
  if (inventory.value("schema", "") != "PhysicsNativeInventoryV1" ||
      !inventory.contains("owner_admissible") || !inventory["owner_admissible"].is_boolean() ||
      !pointer.contains("owner_admissible") || !pointer["owner_admissible"].is_boolean() ||
      pointer["owner_admissible"] != inventory["owner_admissible"])
    throw std::runtime_error("native physics completeness differs from generation pointer");
  // Incomplete upstream metadata does not invalidate a scientifically useful
  // run. It explicitly withholds automatic canonical propagation.
  pointer["pointer"] = {{"path", pointerPath.string()}, {"sha256", hash(pointerPath)}};
  pointer["inventory"]["path"] = path.string();
  return pointer;
}
} // namespace fa
#endif
