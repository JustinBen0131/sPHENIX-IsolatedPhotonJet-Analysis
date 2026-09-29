// The truth-signal photons of Au+Au embedded simulation with their reconstructed isolation: the input of the Au+Au
// isolation package (TreeToHists/isolation.py).
//
//   root -l -q 'TreeToHists/src/FillIsolationPopulation.C("a.eventtree.root b.eventtree.root", "population.root")'
//
// Events: the corrected minimum-bias decision (cent_mb_decision == 1), a finite vertex with |z| < 10 cm, centrality in
// [0, 80). Photons: matched to an analysis truth-signal photon (pho_truth_match_state == 1), finite E_T and eta,
// |eta| < 0.7. No preselection, identification or E_T cut. One row per photon in the TTree "isolation_population":
//   et, eta, centrality, native_bin (cent_native_bin), iso_r03, iso_r04 (the reconstructed isolation witnesses, R = 0.3
//   and 0.4; NaN when the witness is not reconstructed), method (their subtraction method), weight (w_event), file
// plus "sources" (as FillPhotonJetResponse) and "counts" (events and photons at each gate, JSON).
#include "EventTreeInput.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libTreePlayer)
#endif

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace photonjet {
namespace isolation {

// The event and photon gates of the isolation population.
constexpr int kMinimumBias = 1;            // cent_mb_decision of a minimum-bias event
constexpr double kVertexAbsZMax = 10.0;    // cm
constexpr double kCentralityMax = 80.0;    // percent
constexpr int kTruthSignalMatch = 1;       // pho_truth_match_state of a photon matched to an analysis truth-signal photon
constexpr double kPhotonAbsEtaMax = 0.7;
constexpr int kReconstructed = 0;          // isow_reconstructed_or_truth of a reconstructed witness

struct Row {
  double et, eta, centrality, isoR03, isoR04, weight;
  int nativeBin, method, file;
};

inline void fillFile(const std::string& path, int fileIndex, TTree& table, Row& row, std::map<std::string, long>& counts) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  if (!isSimulation(*file) || fileSystem(*file) != System::auau)
    throw std::runtime_error("the isolation population needs Au+Au simulation: " + path);
  TTreeReader reader("eventTree", file.get());
  TTreeReaderValue<int> record(reader, "evt_record"), minimumBias(reader, "cent_mb_decision"),
      nativeBin(reader, "cent_native_bin");
  TTreeReaderValue<double> vz(reader, "vtx_z"), centrality(reader, "cent_value"), weight(reader, "w_event");
  TTreeReaderArray<double> et(reader, "pho_et"), eta(reader, "pho_eta"), cone(reader, "isow_cone_sum");
  TTreeReaderArray<int> match(reader, "pho_truth_match_state"), kind(reader, "isow_reconstructed_or_truth"),
      method(reader, "isow_subtraction_method");
  while (reader.Next()) {
    if (*record != 0) continue;
    ++counts["events"];
    if (*minimumBias != kMinimumBias) continue;
    ++counts["events_minimum_bias"];
    if (!(std::isfinite(*vz) && std::fabs(*vz) < kVertexAbsZMax)) continue;
    ++counts["events_vertex"];
    if (!(std::isfinite(*centrality) && *centrality >= 0.0 && *centrality < kCentralityMax)) continue;
    ++counts["events_selected"];
    const size_t n = et.GetSize();
    if (eta.GetSize() != n || match.GetSize() != n || cone.GetSize() != 2 * n || kind.GetSize() != 2 * n)
      throw std::runtime_error("photon groups or isolation witnesses (2 per photon) disagree in " + path);
    for (size_t i = 0; i < n; ++i) {
      ++counts["photons"];
      if (match[i] != kTruthSignalMatch) continue;
      ++counts["photons_truth_signal"];
      if (!(std::isfinite(et[i]) && std::isfinite(eta[i]) && std::fabs(eta[i]) < kPhotonAbsEtaMax)) continue;
      ++counts["photons_selected"];
      const double nan = std::numeric_limits<double>::quiet_NaN();
      if (method[2 * i] != method[2 * i + 1])
        throw std::runtime_error("the two isolation witnesses of one photon have different methods in " + path);
      row = Row{et[i], eta[i], *centrality, kind[2 * i] == kReconstructed ? cone[2 * i] : nan,
                kind[2 * i + 1] == kReconstructed ? cone[2 * i + 1] : nan, *weight, *nativeBin, method[2 * i], fileIndex};
      table.Fill();
    }
  }
  // a missing branch or a read error ends the loop like the end of the tree; only the end is accepted
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
}

inline int run(const char* inputs, const char* output) {
  const std::vector<std::string> paths = inputPaths(inputs);
  const std::string sources = sourceRecords(paths);  // before the output is opened (sourceRecord)
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  auto* table = new TTree("isolation_population", "truth-signal photons of Au+Au embedded simulation");  // owned by out
  Row row{};
  table->Branch("et", &row.et, "et/D");
  table->Branch("eta", &row.eta, "eta/D");
  table->Branch("centrality", &row.centrality, "centrality/D");
  table->Branch("native_bin", &row.nativeBin, "native_bin/I");
  table->Branch("iso_r03", &row.isoR03, "iso_r03/D");
  table->Branch("iso_r04", &row.isoR04, "iso_r04/D");
  table->Branch("method", &row.method, "method/I");
  table->Branch("weight", &row.weight, "weight/D");
  table->Branch("file", &row.file, "file/I");
  std::map<std::string, long> counts;
  for (size_t i = 0; i < paths.size(); ++i) fillFile(paths[i], static_cast<int>(i), *table, row, counts);
  out->cd();
  table->Write("", TObject::kOverwrite);  // one key: an autosave cycle is replaced
  std::string json = "{";
  for (const auto& [key, value] : counts) json += (json.size() > 1 ? "," : "") + ("\"" + key + "\":" + std::to_string(value));
  TObjString((json + "}").c_str()).Write("counts");
  TObjString(sources.c_str()).Write("sources");
  out->Close();
  return 0;
}

}  // namespace isolation
}  // namespace photonjet

// inputs: whitespace-separated eventTree paths of Au+Au embedded simulation.
int FillIsolationPopulation(const char* inputs, const char* output) { return photonjet::isolation::run(inputs, output); }
