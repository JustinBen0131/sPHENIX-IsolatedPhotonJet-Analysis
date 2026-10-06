/* TreeToHists.h: what the driver's two files share (TreeToHists.cc, Isolation.cc). */
#pragma once

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace treetohists {
using Json = nlohmann::json;
std::string sha256File(const std::string &path);     /* TreeToHists.cc */
int isolation(const std::vector<std::string> &args);  /* Isolation.cc: tree_to_hists isolation ... */
}  /* namespace treetohists */
