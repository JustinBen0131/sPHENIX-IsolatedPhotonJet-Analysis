// eventTree quick start in C++: select analysis events and histogram back-to-back photon-jet xJ.
//
//   root -l -b -q 'quick_start.C("file.eventtree.root")'
//
// TTreeReaderArray reads the collection branches (C arrays sized by their group counter, e.g. n_pho).
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iostream>

void quick_start(const char* path) {
  TFile file(path);
  TTreeReader reader("eventTree", &file);
  TTreeReaderValue<int> record(reader, "evt_record");
  TTreeReaderValue<int> nodes_ready(reader, "evt_nodes_ready");
  TTreeReaderValue<int> vertex_valid(reader, "vtx_valid");
  TTreeReaderValue<int> terminal(reader, "evt_terminal_status");
  TTreeReaderArray<double> photon_et(reader, "pho_et");
  TTreeReaderArray<double> pair_xj(reader, "pair_xjgamma");
  TTreeReaderArray<double> pair_dphi(reader, "pair_delta_phi");
  TTreeReaderArray<int> pair_photon(reader, "pair_pho");
  TTreeReaderArray<int> pair_jet(reader, "pair_jet");
  TTreeReaderArray<double> jet02_pt(reader, "jet02_pt");
  TTreeReaderArray<double> jet03_pt(reader, "jet03_pt");
  TTreeReaderArray<double> jet04_pt(reader, "jet04_pt");

  TH1D xj("xj", "back-to-back photon + R=0.4 jet (pT > 10 GeV);x_{J#gamma};pairs", 20, 0.0, 2.0);
  long long entries = 0, analysis = 0, with_photon = 0;
  while (reader.Next()) {
    ++entries;
    // The documented reader predicate: reconstruction finished with a valid vertex, and the established
    // event classification accepted the event (evt_terminal_status 0).
    if (*record != 0 || *nodes_ready != 1 || *vertex_valid != 1 || *terminal != 0) continue;
    ++analysis;
    if (photon_et.GetSize() > 0) ++with_photon;
    // pair_photon[i] indexes this event's photons; pair_jet[i] indexes the concatenated reco jets in
    // meta/collections order (jet02, jet03, jet04, then jetunsub* in Au+Au). Keep R = 0.4 jets above 10 GeV.
    const int start04 = jet02_pt.GetSize() + jet03_pt.GetSize();
    for (size_t i = 0; i < pair_xj.GetSize(); ++i) {
      const int local04 = pair_jet[i] - start04;
      if (local04 < 0 || local04 >= static_cast<int>(jet04_pt.GetSize()) || jet04_pt[local04] <= 10.0) continue;
      if (pair_dphi[i] > 7.0 * TMath::Pi() / 8.0 && !std::isnan(pair_xj[i]) && photon_et[pair_photon[i]] > 0) xj.Fill(pair_xj[i]);
    }
  }
  std::cout << path << ": " << entries << " entries, " << analysis << " analysis events, " << with_photon
            << " with a photon candidate, " << xj.GetEntries() << " back-to-back pairs, mean xJ " << xj.GetMean() << "\n";
}
