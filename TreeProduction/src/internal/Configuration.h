#ifndef PHOTONJET_CONFIGURATION_H
#define PHOTONJET_CONFIGURATION_H

/**
 * @file Configuration.h
 * @brief Read-only access to the parsed TreeProduction configuration.
 *
 * The production frontend parses the user-facing YAML configuration once and
 * passes the resulting JSON document to the native production layer.
 *
 * ConfigurationNode is a deliberately thin view over that JSON. It provides
 * the small YAML-like access interface used by Production.cc without assigning
 * defaults, interpreting physics values, or changing configuration semantics.
 *
 * Production.cc remains authoritative for:
 *   - required versus optional fields,
 *   - profile interpretation,
 *   - defaults,
 *   - consistency checks, and
 *   - every production or physics decision derived from the configuration.
 *
 * Missing nodes remain distinct from explicit false, zero, empty-string and
 * empty-container values so native validation can distinguish absence from a
 * deliberately configured value.
 */

#include <nlohmann/json.hpp>
#include <cstddef>
#include <string>
#include <vector>

namespace photonjet::production::detail
{

/**
 * Lightweight read-only node wrapper for the parsed configuration tree.
 *
 * A default-constructed node represents a missing lookup. Existing JSON values
 * retain their native type and are converted only when the caller explicitly
 * requests a type through as<T>().
 *
 * This class intentionally contains no TreeProduction policy. It only exposes
 * object lookup, array lookup, iteration and typed value access to the native
 * plan resolver.
 */
class ConfigurationNode
{
 public:
  ConfigurationNode() = default;
  ConfigurationNode(const nlohmann::json& value) : m_value(value) {}

  /**
   * A node is present when its stored JSON value is not null.
   *
   * This permits Production.cc to distinguish a missing/null field from valid
   * false, zero or empty-container values.
   */
  explicit operator bool() const { return !m_value.is_null(); }

  /** Return whether this node contains JSON null. */
  bool IsNull() const { return m_value.is_null(); }

  /** Return whether this node is a JSON array. */
  bool IsSequence() const { return m_value.is_array(); }

  /** Return the underlying JSON container/value size. */
  std::size_t size() const { return m_value.size(); }

  /**
   * Return the keys of a JSON object.
   *
   * Non-object nodes return an empty vector. No ordering or interpretation is
   * added beyond the underlying nlohmann::json object traversal.
   */
  std::vector<std::string> keys() const
  {
    std::vector<std::string> result;
    if (m_value.is_object())
      for (auto it = m_value.cbegin(); it != m_value.cend(); ++it) result.push_back(it.key());
    return result;
  }

  /**
   * Read one object member.
   *
   * A missing key, or lookup on a non-object node, returns an empty
   * ConfigurationNode rather than manufacturing a value.
   */
  ConfigurationNode operator[](const std::string& key) const
  {
    if (!m_value.is_object() || !m_value.contains(key)) return {};
    return m_value.at(key);
  }

  /**
   * Read one array element.
   *
   * An out-of-range index, or lookup on a non-array node, returns an empty
   * ConfigurationNode.
   */
  ConfigurationNode operator[](std::size_t index) const
  {
    if (!m_value.is_array() || index >= m_value.size()) return {};
    return m_value.at(index);
  }

  /**
   * Convert the current JSON value to the type explicitly requested by the
   * native plan resolver.
   *
   * Type validation and any resulting conversion exception come directly from
   * nlohmann::json; this wrapper does not coerce incompatible configuration.
   */
  template <class T> T as() const { return m_value.get<T>(); }

  /**
   * Signed-index convenience overload.
   *
   * Negative indices represent an invalid/missing lookup; nonnegative indices
   * use the bounded size_t overload above.
   */
  ConfigurationNode operator[](int index) const
  {
    return index < 0 ? ConfigurationNode{} : (*this)[static_cast<std::size_t>(index)];
  }

  /**
   * Minimal iterator adapter used for range-based traversal of configuration
   * sequences. Dereferencing produces another read-only ConfigurationNode.
   */
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

  /** Begin iteration over the underlying JSON container. */
  Iterator begin() const { return Iterator(m_value.cbegin()); }

  /** End iteration over the underlying JSON container. */
  Iterator end() const { return Iterator(m_value.cend()); }

 private:
  /**
   * Parsed configuration value owned by this node.
   *
   * The wrapper stores JSON structure only; semantic interpretation belongs to
   * the callers in the native production-plan implementation.
   */
  nlohmann::json m_value;
};

}  // namespace photonjet::production::detail

#endif
