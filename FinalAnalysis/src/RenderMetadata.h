/* Shared metadata checks for C++ analysis and ROOT rendering.
   seal/verify bind the JSON description; product receipts separately bind the actual ROOT-file bytes.
   Renderers use both checks before labeling a figure with the saved selection and binning. */
#ifndef PHOTONJET_RENDER_METADATA_H
#define PHOTONJET_RENDER_METADATA_H
#include <nlohmann/json.hpp>
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

namespace photonjet_render {
using Json = nlohmann::json;
/* Exact argv: neither file names nor scientific metadata are evaluated by a shell. */
inline std::string sha256(const std::string &value, bool file = false) {
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
  return hash;
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
