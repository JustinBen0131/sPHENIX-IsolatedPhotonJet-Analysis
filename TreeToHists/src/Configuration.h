/* Shared YAML configuration reader for TreeToHists and FinalAnalysis.
   libyaml parses the document; settings are passed to the existing ROOT TEnv
   readers, which retain ownership of every physics choice and validation. */
#pragma once

#include <yaml.h>
#include <nlohmann/json.hpp>
#include "TEnv.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace photonjet_config {
using Json = nlohmann::json;

inline bool isYaml(const std::string &path) {
  const auto extension = std::filesystem::path(path).extension();
  return extension == ".yaml" || extension == ".yml";
}

class YamlReader {
 public:
  explicit YamlReader(const std::string &path) : path_(path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot read configuration " + path);
    text_.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
    if (!yaml_parser_initialize(&parser_)) throw std::runtime_error("cannot initialize YAML reader");
    yaml_parser_set_input_string(&parser_, reinterpret_cast<const unsigned char *>(text_.data()), text_.size());
  }
  ~YamlReader() { if (eventReady_) yaml_event_delete(&event_); yaml_parser_delete(&parser_); }
  YamlReader(const YamlReader &) = delete;
  YamlReader &operator=(const YamlReader &) = delete;

  Json read() {
    next(); expect(YAML_STREAM_START_EVENT); next(); expect(YAML_DOCUMENT_START_EVENT); next();
    Json value = node(0);
    expect(YAML_DOCUMENT_END_EVENT); next();
    if (event_.type != YAML_STREAM_END_EVENT) fail("expected one YAML document");
    return value;
  }

 private:
  std::string path_, text_;
  yaml_parser_t parser_{};
  yaml_event_t event_{};
  bool eventReady_ = false;
  std::map<std::string, Json> anchors_;

  [[noreturn]] void fail(const std::string &message) const {
    throw std::runtime_error(path_ + ":" + std::to_string(event_.start_mark.line + 1) + ": " + message);
  }
  static std::string text(const yaml_char_t *value) { return value ? reinterpret_cast<const char *>(value) : ""; }
  void next() {
    if (eventReady_) yaml_event_delete(&event_);
    eventReady_ = false;
    if (!yaml_parser_parse(&parser_, &event_))
      throw std::runtime_error(path_ + ":" + std::to_string(parser_.problem_mark.line + 1) + ": " +
                               (parser_.problem ? parser_.problem : "invalid YAML"));
    eventReady_ = true;
  }
  void expect(yaml_event_type_t type) const { if (event_.type != type) fail("invalid YAML structure"); }

  Json scalar() const {
    const std::string value(reinterpret_cast<const char *>(event_.data.scalar.value), event_.data.scalar.length);
    const std::string tag = text(event_.data.scalar.tag);
    if (tag == YAML_STR_TAG || (tag.empty() && event_.data.scalar.style != YAML_PLAIN_SCALAR_STYLE)) return value;
    if (!tag.empty() && tag != YAML_INT_TAG && tag != YAML_FLOAT_TAG && tag != YAML_BOOL_TAG && tag != YAML_NULL_TAG)
      fail("unsupported YAML scalar tag " + tag);
    if (value == "null" || value == "Null" || value == "NULL" || value == "~" || value.empty()) {
      if (!tag.empty() && tag != YAML_NULL_TAG) fail("invalid tagged scalar " + value);
      return nullptr;
    }
    if (value == "true" || value == "True" || value == "TRUE" || value == "false" || value == "False" || value == "FALSE") {
      if (!tag.empty() && tag != YAML_BOOL_TAG) fail("invalid tagged scalar " + value);
      return value == "true" || value == "True" || value == "TRUE";
    }
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lower == ".nan" || lower == ".inf" || lower == "+.inf" || lower == "-.inf") fail("number must be finite");
    std::string number = value;
    number.erase(std::remove(number.begin(), number.end(), '_'), number.end());
    static const std::regex numeric("[+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?");
    if (std::regex_match(number, numeric)) {
      if (number.front() == '+') number.erase(0, 1);
      const size_t first = number.front() == '-' ? 1 : 0;
      if (number[first] == '.') number.insert(first, "0");
      while (number.size() > first + 1 && number[first] == '0' && number[first + 1] >= '0' && number[first + 1] <= '9')
        number.erase(first, 1);
      const auto exponent = number.find_first_of("eE");
      const size_t end = exponent == std::string::npos ? number.size() : exponent;
      if (number[end - 1] == '.') number.insert(end, "0");
      const Json parsed = Json::parse(number, nullptr, false);
      if (!parsed.is_number() || !std::isfinite(parsed.get<double>())) fail("number must be finite and representable");
      if (tag == YAML_BOOL_TAG || tag == YAML_NULL_TAG || (tag == YAML_INT_TAG && !parsed.is_number_integer()))
        fail("invalid tagged scalar " + value);
      if (tag == YAML_FLOAT_TAG) return parsed.get<double>();
      return parsed;
    }
    if (!tag.empty()) fail("invalid tagged scalar " + value);
    return value;
  }

  Json node(size_t depth) {
    if (depth > 64) fail("configuration nesting is too deep");
    Json value;
    std::string anchor;
    if (event_.type == YAML_SCALAR_EVENT) {
      anchor = text(event_.data.scalar.anchor); value = scalar(); next();
    } else if (event_.type == YAML_ALIAS_EVENT) {
      const auto name = text(event_.data.alias.anchor);
      const auto it = anchors_.find(name);
      if (it == anchors_.end()) fail("unknown or recursive YAML alias " + name);
      value = it->second; next();
    } else if (event_.type == YAML_SEQUENCE_START_EVENT) {
      anchor = text(event_.data.sequence_start.anchor);
      const std::string tag = text(event_.data.sequence_start.tag);
      if (!tag.empty() && tag != YAML_SEQ_TAG) fail("unsupported YAML sequence tag");
      value = Json::array(); next();
      while (event_.type != YAML_SEQUENCE_END_EVENT) value.push_back(node(depth + 1));
      next();
    } else if (event_.type == YAML_MAPPING_START_EVENT) {
      anchor = text(event_.data.mapping_start.anchor);
      const std::string tag = text(event_.data.mapping_start.tag);
      if (!tag.empty() && tag != YAML_MAP_TAG) fail("unsupported YAML mapping tag");
      value = Json::object(); std::set<std::string> keys; next();
      while (event_.type != YAML_MAPPING_END_EVENT) {
        expect(YAML_SCALAR_EVENT);
        const std::string key(reinterpret_cast<const char *>(event_.data.scalar.value), event_.data.scalar.length);
        if (key.empty() || !keys.insert(key).second) fail("empty or duplicate YAML key " + key);
        next(); Json item = node(depth + 1);
        if (key == "<<") {
          const Json maps = item.is_array() ? item : Json::array({item});
          for (const auto &map : maps) {
            if (!map.is_object()) fail("YAML merge needs a mapping");
            for (const auto &[k, v] : map.items()) if (!value.contains(k)) value[k] = v;
          }
        } else value[key] = item;
      }
      next();
    } else fail("expected a YAML value");
    if (!anchor.empty()) anchors_[anchor] = value;
    return value;
  }
};

inline Json readYaml(const std::string &path) { return YamlReader(path).read(); }

/* Nested keys become the existing dotted setting names. Scalar sequences are
   passed in their original order; no setting is inferred or selected here. */
inline void flatten(const Json &value, const std::string &prefix, std::map<std::string, std::string> &settings) {
  if (value.is_object()) {
    for (const auto &[key, child] : value.items()) flatten(child, prefix.empty() ? key : prefix + "." + key, settings);
    return;
  }
  auto scalar = [](const Json &v) {
    if (v.is_structured()) throw std::runtime_error("settings sequences must contain scalars");
    return v.is_string() ? v.get<std::string>() : v.is_null() ? std::string() : v.is_boolean() ? (v.get<bool>() ? "1" : "0") : v.dump();
  };
  std::string text;
  if (value.is_array()) {
    for (const auto &v : value) { if (!text.empty()) text += " "; text += scalar(v); }
  } else text = scalar(value);
  if (prefix.empty() || !settings.emplace(prefix, text).second)
    throw std::runtime_error("empty or duplicate flattened setting " + prefix);
}

inline void readSettings(const std::string &path, TEnv &env) {
  if (!isYaml(path)) {
    if (env.ReadFile(path.c_str(), kEnvUser) != 0) throw std::runtime_error("cannot read settings " + path);
    return;
  }
  const Json value = readYaml(path);
  if (!value.is_object()) throw std::runtime_error(path + ": settings must be a YAML mapping");
  std::map<std::string, std::string> settings;
  flatten(value, "", settings);
  for (const auto &[key, text] : settings) env.SetValue(key.c_str(), text.c_str());
}
}  // namespace photonjet_config
