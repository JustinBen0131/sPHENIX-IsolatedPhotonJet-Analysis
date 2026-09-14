#ifndef RJ_CEMC_STATUS_RULES_H
#define RJ_CEMC_STATUS_RULES_H

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace rj_cemc_status
{
constexpr std::size_t channel_count = 24576;
constexpr std::uint8_t hot_bit = 1;

inline void require(bool valid, const std::string& reason)
{
  if (!valid) throw std::runtime_error("CEMC_STATUS_CONTRACT: " + reason);
}

inline std::vector<std::uint8_t> validated_mask(const std::vector<int>& status)
{
  require(status.size() == channel_count, "incomplete map channel coverage");
  std::vector<std::uint8_t> mask;
  mask.reserve(status.size());
  for (const auto value : status)
  {
    // Missing CDB integer fields return INT_MIN. Never treat them as GOOD.
    require(value >= 0 && value <= 255, "missing or invalid map status field");
    mask.push_back(value > 0);
  }
  return mask;
}

inline void validate_flags(const std::vector<std::uint8_t>& flags,
                           const std::vector<std::uint8_t>& mask)
{
  require(flags.size() == channel_count && mask.size() == channel_count,
          "incomplete tower container or map");
  for (std::size_t i = 0; i < flags.size(); ++i)
    require(!mask[i] || (flags[i] & hot_bit),
            "map-rejected tower lacks HOT flag at channel " + std::to_string(i));
}

inline void add_missing_hot_flags(std::vector<std::uint8_t>& flags,
                                  const std::vector<std::uint8_t>& mask)
{
  require(flags.size() == channel_count && mask.size() == channel_count,
          "incomplete recovery inputs");
  for (std::size_t i = 0; i < flags.size(); ++i)
    if (mask[i]) flags[i] |= hot_bit;
  validate_flags(flags, mask);
}
}
#endif
