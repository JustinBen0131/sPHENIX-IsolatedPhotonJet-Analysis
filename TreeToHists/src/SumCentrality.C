// One centrality histogram summed over many files, written under the same key, so TreeToHists/centrality.py can read
// the per-group simulation sums instead of every eventTree part; sums of sums work the same way.
//
//   root -l -q 'TreeToHists/src/SumCentrality.C("a.eventtree.root b.eventtree.root", "sum.root",
//                                               "hist/SIM/h_centrality")'
//
// Every input must carry the histogram with identical bin edges, compared explicitly (TH1::Add alone accepts other
// binnings in some ROOT versions). The output holds the sum at the input key plus "sources" (the input paths, JSON) and
// "counts" (files summed and total entries, JSON). Failures return non-zero before the output is created; no
// exceptions, so interpreted and compiled runs behave the same.
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TObjString.h>

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {
bool sameEdges(const TH1& a, const TH1& b) {
  const int n = a.GetNbinsX();
  if (n != b.GetNbinsX()) return false;
  for (int i = 1; i <= n + 1; ++i)
    if (a.GetXaxis()->GetBinLowEdge(i) != b.GetXaxis()->GetBinLowEdge(i)) return false;
  return true;
}
}  // namespace

// inputs: whitespace-separated paths (eventTree files or earlier sums); key: the histogram to sum.
int SumCentrality(const char* inputs, const char* output, const char* key) {
  std::vector<std::string> paths;
  std::istringstream stream(inputs);
  for (std::string path; stream >> path;) paths.push_back(path);
  if (paths.empty()) {
    std::cerr << "[SumCentrality] no input files" << std::endl;
    return 1;
  }
  std::unique_ptr<TH1> sum;
  std::string sources = "[";
  for (const auto& path : paths) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::cerr << "[SumCentrality] cannot open " << path << std::endl;
      return 2;
    }
    auto* h = dynamic_cast<TH1*>(f->Get(key));
    if (!h) {
      std::cerr << "[SumCentrality] " << path << " has no " << key << std::endl;
      return 3;
    }
    if (!sum) {
      sum.reset(static_cast<TH1*>(h->Clone()));
      sum->SetDirectory(nullptr);
    } else if (!sameEdges(*sum, *h) || !sum->Add(h)) {
      std::cerr << "[SumCentrality] " << path << ": " << key << " has other bin edges" << std::endl;
      return 4;
    }
    sources += (sources.size() > 1 ? ",\"" : "\"") + path + "\"";
  }
  const std::string name(key);
  const auto slash = name.rfind('/');
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "[SumCentrality] cannot create " << output << std::endl;
    return 5;
  }
  TDirectory* dir = out.get();
  if (slash != std::string::npos) {
    // mkdir("a/b") returns the top directory "a"; the leaf is fetched by its path
    out->mkdir(name.substr(0, slash).c_str(), "", true);
    dir = out->GetDirectory(name.substr(0, slash).c_str());
  }
  if (!dir) {
    std::cerr << "[SumCentrality] cannot create the directory of " << key << std::endl;
    return 6;
  }
  dir->cd();
  sum->Write(name.substr(slash == std::string::npos ? 0 : slash + 1).c_str());
  out->cd();
  TObjString((sources + "]").c_str()).Write("sources");
  TObjString(("{\"files\":" + std::to_string(paths.size()) + ",\"entries\":" + std::to_string(static_cast<long long>(sum->GetEntries())) + "}").c_str())
      .Write("counts");
  out->Close();
  return 0;
}
