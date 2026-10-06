#ifndef PHOTONJET_CONFIGURATION_H
#define PHOTONJET_CONFIGURATION_H

/**
 * @file Configuration.h
 * Read-only access to the production configuration, which the runner parses
 * from YAML into JSON once; the native plan (Production.cc) reads it through
 * this node interface and owns every physics decision.
 */

#include <nlohmann/json.hpp>
#include <cstddef>
#include <string>
#include <vector>

namespace photonjet::production::detail
{
// Configuration is parsed once by the launch layer, without selecting or
// changing any physics values. Missing nodes remain distinguishable from
// false, zero and empty containers in the native plan's validation.
class ConfigurationNode
{
 public:
  ConfigurationNode() = default;
  ConfigurationNode(const nlohmann::json& value) : m_value(value) {}

  explicit operator bool() const { return !m_value.is_null(); }
  bool IsNull() const { return m_value.is_null(); }
  bool IsSequence() const { return m_value.is_array(); }
  std::size_t size() const { return m_value.size(); }

  // Keys of a mapping (empty for any other node).
  std::vector<std::string> keys() const
  {
    std::vector<std::string> result;
    if (m_value.is_object())
      for (auto it = m_value.cbegin(); it != m_value.cend(); ++it) result.push_back(it.key());
    return result;
  }

  ConfigurationNode operator[](const std::string& key) const
  {
    if (!m_value.is_object() || !m_value.contains(key)) return {};
    return m_value.at(key);
  }
  ConfigurationNode operator[](std::size_t index) const
  {
    if (!m_value.is_array() || index >= m_value.size()) return {};
    return m_value.at(index);
  }
  template <class T> T as() const { return m_value.get<T>(); }
  ConfigurationNode operator[](int index) const
  {
    return index < 0 ? ConfigurationNode{} : (*this)[static_cast<std::size_t>(index)];
  }

  class Iterator
  {
   public:
    explicit Iterator(nlohmann::json::const_iterator position) : m_position(position) {}
    ConfigurationNode operator*() const { return *m_position; }
    Iterator& operator++() { ++m_position; return *this; }
    bool operator!=(const Iterator& other) const { return m_position != other.m_position; }
   private:
    nlohmann::json::const_iterator m_position;
  };
  Iterator begin() const { return Iterator(m_value.cbegin()); }
  Iterator end() const { return Iterator(m_value.cend()); }

 private:
  nlohmann::json m_value;
};
}  // namespace photonjet::production::detail
#endif
