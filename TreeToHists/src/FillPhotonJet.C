/* Measurement reduction from eventTree v1.

   Fills event counts, photons, ABCD leaders, recoil spectra and Delta-phi
   sidebands, plus event/photon/jet records for whole-event resampling.
   PhotonJet.h owns all cuts, stitching admission and event weights.
   Histograms retain Sumw2; public plans enforce finite weights and input
   provenance. The inventory below defines the saved cells and records.
   Run through tree_to_hists.sh; README.md describes the package contract. */
#include "PhotonJet.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libTreePlayer)
#endif

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TAxis.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/* Fill the photon + recoil-jet histograms from eventTree v1 files. */

/* root -l -q 'TreeToHists/src/FillPhotonJet.C("a.eventtree.root b.eventtree.root", "out.root", "nominal")' */

/* Analysis events: evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && evt_terminal_status == 0.
   Per centrality class <c> (Au+Au: cent0_20, cent20_50, cent50_80; pp: all) and recoil-jet threshold <t>
   (jet5, jet7, jet10, jet12 GeV):
   events_<c>                TH1D  one bin                            all analysis events of the class
   photons_<c>               TH1D  photon ET                          photons passing the kinematic selection (no ID)
   abcd_<c>                  TH2D  region (A,B,C,D) x photon ET       the event-leading photon of each region
   recoil_<c>_<t>            TH3D  region (A,C) x photon ET x xJ      every recoil jet of the A and C leaders
   dphi<R>_<c>_<t>_<eta>     TH3D  photon ET x xJ x |dphi|           R = A, C; every jet of the leader passing
   the recoil cuts except dphi; eta = near
   (|deta| <= 0.8) or far (|deta| > 0.8)
   and one event table per output, for whole-event resampling (bootstrap):
   table_events   run, seq, file, entry, centrality, slice (5% slices of 0-80; -1 in pp), weight
   table_photons  event, region (0-3 = A-D), et, et_bin              the event-leading photon of each region
   table_jets     event, region (0 = A, 1 = C), et_bin, xj, xj_bin, pt, dphi, far, window, thresholds (bit k: pt >= threshold k)
   every jet around the A and C leaders passing the recoil cuts at the lowest threshold except dphi;
   far = 1 when |deta| > 0.8; window = 1 when it also passes the dphi cut (a recoil jet); bins are ROOT
   bin numbers
   Weights are w_event (1 when it is not finite). All histograms keep Sumw2. */

namespace photonjet {

class Histograms {
 public:
  explicit Histograms(const std::vector<CentralityClass>& classes) {
    const auto et = photonEtEdges(), xj = xjEdges(), dphi = deltaPhiEdges();
    const double regions[5] = {0, 1, 2, 3, 4}, recoilRegions[3] = {0, 1, 2};
    for (const auto& c : classes) {
      add1(new TH1D(("events_" + c.name).c_str(), "analysis events;;events", 1, 0.0, 1.0));
      add1(new TH1D(("photons_" + c.name).c_str(), "photons passing the kinematic selection;photon E_{T} [GeV]", et.size() - 1,
                    et.data()));
      add2(new TH2D(("abcd_" + c.name).c_str(), "event-leading photons;region (A,B,C,D);photon E_{T} [GeV]", 4, regions,
                    et.size() - 1, et.data()));
      for (double t : kJetThresholds) {
        const std::string tag = c.name + "_" + thresholdName(t);
        add3(new TH3D(("recoil_" + tag).c_str(), "recoil jets;region (A,C);photon E_{T} [GeV];x_{J#gamma}", 2, recoilRegions,
                      et.size() - 1, et.data(), xj.size() - 1, xj.data()));
        for (const char* r : {"A", "C"})
          for (const char* eta : {"near", "far"})
            add3(new TH3D((std::string("dphi") + r + "_" + tag + "_" + eta).c_str(),
                          "jets around the leader;photon E_{T} [GeV];x_{J#gamma};|#Delta#phi|", et.size() - 1, et.data(),
                          xj.size() - 1, xj.data(), dphi.size() - 1, dphi.data()));
      }
    }
  }
  TH1D& h1(const std::string& name) { return *h1_.at(name); }
  TH2D& h2(const std::string& name) { return *h2_.at(name); }
  TH3D& h3(const std::string& name) { return *h3_.at(name); }
  void write() {
    for (auto& [name, h] : h1_) h->Write();
    for (auto& [name, h] : h2_) h->Write();
    for (auto& [name, h] : h3_) h->Write();
  }

 private:
  void add1(TH1D* h) {
    h->Sumw2();
    h->SetDirectory(nullptr);
    h1_[h->GetName()].reset(h);
  }
  void add2(TH2D* h) {
    h->Sumw2();
    h->SetDirectory(nullptr);
    h2_[h->GetName()].reset(h);
  }
  void add3(TH3D* h) {
    h->Sumw2();
    h->SetDirectory(nullptr);
    h3_[h->GetName()].reset(h);
  }
  std::map<std::string, std::unique_ptr<TH1D>> h1_;
  std::map<std::string, std::unique_ptr<TH2D>> h2_;
  std::map<std::string, std::unique_ptr<TH3D>> h3_;
};

struct EventTable {
  std::vector<Int_t> run, fileIndex;
  std::vector<Long64_t> seq, entry;
  std::vector<Double_t> centrality, weight;
  std::vector<Int_t> slice;
  std::vector<Int_t> photonEvent, photonRegion, photonEtBin;
  std::vector<Double_t> photonEt;
  std::vector<Int_t> jetEvent, jetRegion, jetEtBin, jetXjBin, jetThresholds, jetWindow, jetFar;
  std::vector<Double_t> jetXj, jetPt, jetDeltaPhi;

  void write() {
    TTree events("table_events", "analysis events with at least one region-leading photon");
    Int_t r, f, s;
    Long64_t q, e;
    Double_t c, w;
    events.Branch("run", &r, "run/I");
    events.Branch("seq", &q, "seq/L");
    events.Branch("file", &f, "file/I");
    events.Branch("entry", &e, "entry/L");
    events.Branch("centrality", &c, "centrality/D");
    events.Branch("slice", &s, "slice/I");
    events.Branch("weight", &w, "weight/D");
    for (size_t k = 0; k < run.size(); ++k) {
      r = run[k], q = seq[k], f = fileIndex[k], e = entry[k], c = centrality[k], s = slice[k], w = weight[k];
      events.Fill();
    }
    events.Write("", TObject::kOverwrite);  /* one key: an autosave cycle is replaced */
    TTree photons("table_photons", "event-leading photon of each region");
    Int_t pe, pr, pb;
    Double_t pet;
    photons.Branch("event", &pe, "event/I");
    photons.Branch("region", &pr, "region/I");
    photons.Branch("et", &pet, "et/D");
    photons.Branch("et_bin", &pb, "et_bin/I");
    for (size_t k = 0; k < photonEvent.size(); ++k) {
      pe = photonEvent[k], pr = photonRegion[k], pet = photonEt[k], pb = photonEtBin[k];
      photons.Fill();
    }
    photons.Write("", TObject::kOverwrite);  /* one key: an autosave cycle is replaced */
    TTree jets("table_jets", "jets around the A and C leaders at the lowest threshold; window = 1 for recoil jets");
    Int_t je, jr, jb, jx, jt, jw, jf;
    Double_t xj, pt, dphi;
    jets.Branch("event", &je, "event/I");
    jets.Branch("region", &jr, "region/I");
    jets.Branch("et_bin", &jb, "et_bin/I");
    jets.Branch("xj", &xj, "xj/D");
    jets.Branch("xj_bin", &jx, "xj_bin/I");
    jets.Branch("pt", &pt, "pt/D");
    jets.Branch("dphi", &dphi, "dphi/D");
    jets.Branch("far", &jf, "far/I");
    jets.Branch("window", &jw, "window/I");
    jets.Branch("thresholds", &jt, "thresholds/I");
    for (size_t k = 0; k < jetEvent.size(); ++k) {
      je = jetEvent[k], jr = jetRegion[k], jb = jetEtBin[k], xj = jetXj[k], jx = jetXjBin[k], pt = jetPt[k];
      dphi = jetDeltaPhi[k], jt = jetThresholds[k], jw = jetWindow[k], jf = jetFar[k];
      jets.Fill();
    }
    jets.Write("", TObject::kOverwrite);  /* one key: an autosave cycle is replaced */
  }
};

inline void fillFile(const std::string& path, int fileIndex, const Profile& requested,
                     const std::vector<CentralityClass>& classes, Histograms& h, EventTable& table,
                     const CentralityWeights& weights = {}, const std::string& dataTrigger = "",
                     Long64_t progressOffset = 0) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  const Profile profile = validateSelectionInput(*file, requested);
  if (fileSystem(*file) != profile.system) throw std::runtime_error("profile system differs from " + path);
  TTreeReader reader("eventTree", file.get());
  const VertexEventWeight vertexWeight(*file,reader);
  const PpSimulationGate ppSelection(*file,reader);
  const DataEventSelection dataSelection(*file,profile);
  const ULong64_t triggerMask = dataTrigger.empty() ? 0 : dataTriggerMask(*file,dataTrigger);
  std::unique_ptr<TTreeReaderValue<ULong64_t>> scaled;
  if(!dataTrigger.empty() || dataSelection.active) scaled.reset(new TTreeReaderValue<ULong64_t>(reader,"trig_scaled"));
  std::unique_ptr<TTreeReaderValue<int>> mb;
  if(dataSelection.minimumBias) mb.reset(new TTreeReaderValue<int>(reader,"cent_mb_decision"));
  TTreeReaderValue<int> record(reader, "evt_record"), ready(reader, "evt_nodes_ready"), vtxValid(reader, "vtx_valid"),
      terminal(reader, "evt_terminal_status"), run(reader, "evt_run");
  TTreeReaderValue<Long64_t> seq(reader, "evt_seq");
  TTreeReaderValue<double> vz(reader, "vtx_z"), centrality(reader, "cent_value"), weightValue(reader, "w_event");
  TTreeReaderArray<double> et(reader, "pho_et"), eta(reader, "pho_eta"), phi(reader, "pho_phi"),
      score(reader, profile.scoreBranch.c_str());
  TTreeReaderArray<int> presel(reader, "pho_presel_active"), state(reader, profile.stateBranch.c_str());
  TTreeReaderArray<double> cone(reader, "isow_cone_sum"), threshold(reader, "isow_threshold"),
      sideband(reader, "isow_sideband_threshold");
  TTreeReaderArray<int> method(reader, "isow_subtraction_method"), role(reader, "isow_reconstructed_or_truth");
  TTreeReaderArray<int> idClass(reader, "pho_id_class"), isoClass(reader, "pho_iso_class");
  TTreeReaderArray<double> jetPt(reader, "jet04_pt"), jetEta(reader, "jet04_eta"), jetPhi(reader, "jet04_phi");
  const OwnershipGate owner(reader);  /* a stitched sample's events outside its window belong to another sample */

  while (reader.Next()) {
    reportProgress(progressOffset+reader.GetCurrentEntry()+1);
    validateSelectionEvent(profile, *centrality, et, score, state, cone, threshold, sideband, method, role, idClass, isoClass);
    if (*record != 0 || *ready != 1 || *vtxValid != 1 || *terminal != 0) continue;
    if (!ppSelection.accepts()) continue;
    if (!owner.owned()) continue;
    if(profile.configured && (!std::isfinite(*vz) || std::fabs(*vz)>=profile.vertexAbsZMax)) continue;
    if(!dataSelection.accepts(*run,scaled?**scaled:0,mb?**mb:0)) continue;
    if(!dataTrigger.empty() && !(**scaled & triggerMask)) continue;
    const double product = *weightValue * vertexWeight() * weights(*centrality);
    if (strictReduction && !std::isfinite(product)) throw std::runtime_error("nonfinite event weight in " + path);
    const double w = std::isfinite(product) ? product : 1.0;
    int k = 0;  /* the event's centrality class (classes are disjoint) */
    if (profile.system == System::auau) {
      while (k < static_cast<int>(classes.size()) && !(classes[k].low <= *centrality && *centrality < classes[k].high)) ++k;
      if (k == static_cast<int>(classes.size())) continue;
    }
    const CentralityClass& c = classes[k];
    h.h1("events_" + c.name).Fill(0.5, w);
    const size_t nPhoton = et.GetSize();
    if (nPhoton == 0) continue;
    if (cone.GetSize() != 2 * nPhoton) throw std::runtime_error("isolation witnesses are not 2 per photon in " + path);
    if (score.GetSize() != nPhoton || state.GetSize() != nPhoton)
      throw std::runtime_error("photon-ID branches pho_bdt_score, pho_bdt_state are not filled in " + path);
    if (profile.augmentedClasses && (idClass.GetSize() != nPhoton || isoClass.GetSize() != nPhoton))
      throw std::runtime_error("the augmented classes pho_id_class, pho_iso_class are not filled in " + path);
    TH1D& kinematic = h.h1("photons_" + c.name);
    for (size_t i = 0; i < nPhoton; ++i)
      if (photonKinematics(profile, *vz, et[i], eta[i], phi[i])) kinematic.Fill(et[i], w);
    const std::array<long, 4> leader = regionLeaders(profile, *vz, *centrality, nPhoton, et, eta, phi, score, presel, state,
                                                     cone, threshold, sideband, method, role, idClass, isoClass);
    if (leader[kA] < 0 && leader[kB] < 0 && leader[kC] < 0 && leader[kD] < 0) continue;
    TH2D& abcd = h.h2("abcd_" + c.name);
    const int event = static_cast<int>(table.run.size());
    table.run.push_back(*run);
    table.seq.push_back(*seq);
    table.fileIndex.push_back(fileIndex);
    table.entry.push_back(reader.GetCurrentEntry());
    table.centrality.push_back(*centrality);
    table.slice.push_back(profile.system == System::auau ? centralitySlice(*centrality) : -1);
    table.weight.push_back(w);
    for (int r = 0; r < 4; ++r) {
      if (leader[r] < 0) continue;
      abcd.Fill(r + 0.5, et[leader[r]], w);
      table.photonEvent.push_back(event);
      table.photonRegion.push_back(r);
      table.photonEt.push_back(et[leader[r]]);
      table.photonEtBin.push_back(abcd.GetYaxis()->FindFixBin(et[leader[r]]));
    }
    for (int r : {kA, kC}) {
      if (leader[r] < 0) continue;
      const size_t i = leader[r];
      const std::string name = r == kA ? "A" : "C";
      for (size_t j = 0; j < jetPt.GetSize(); ++j) {
        const double xj = jetPt[j] / et[i];
        if (!(xj >= 0.0 && xj < 3.0)) continue;
        const double dphi = wrappedAbsDeltaPhi(jetPhi[j], phi[i]);
        int passed = 0, open = 0;
        for (size_t t = 0; t < kJetThresholds.size(); ++t) {
          Profile p = profile;
          p.jetPtMin = kJetThresholds[t];
          const std::string tag = c.name + "_" + thresholdName(kJetThresholds[t]);
          if (recoilJet(p, eta[i], phi[i], jetPt[j], jetEta[j], jetPhi[j])) {
            h.h3("recoil_" + tag).Fill(r == kA ? 0.5 : 1.5, et[i], xj, w);
            passed |= 1 << t;
          }
          /* the same cuts without the dphi requirement, for a data-driven fake-jet estimate */
          Profile anyAngle = p;
          anyAngle.deltaPhiMin = 0.0;
          if (recoilJet(anyAngle, eta[i], phi[i], jetPt[j], jetEta[j], jetPhi[j])) {
            open |= 1 << t;
            const char* side = std::fabs(jetEta[j] - eta[i]) > kFarDeltaEta ? "far" : "near";
            h.h3("dphi" + name + "_" + tag + "_" + side).Fill(et[i], xj, dphi, w);
          }
        }
        if (open & 1) {
          TH3D& recoil = h.h3("recoil_" + c.name + "_" + thresholdName(kJetThresholds[0]));
          table.jetEvent.push_back(event);
          table.jetRegion.push_back(r == kA ? 0 : 1);
          table.jetEtBin.push_back(recoil.GetYaxis()->FindFixBin(et[i]));
          table.jetXj.push_back(xj);
          table.jetXjBin.push_back(recoil.GetZaxis()->FindFixBin(xj));
          table.jetPt.push_back(jetPt[j]);
          table.jetDeltaPhi.push_back(dphi);
          table.jetFar.push_back(std::fabs(jetEta[j] - eta[i]) > kFarDeltaEta ? 1 : 0);
          table.jetWindow.push_back(passed & 1);
          table.jetThresholds.push_back(open);
        }
      }
    }
  }
  /* a missing branch or a read error ends the loop like the end of the tree; only the end is accepted */
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
}

}  /* namespace photonjet */

/* inputs: whitespace-separated eventTree paths. profile: "nominal" (the analysis), a variant of it (PhotonJet.h
   namedProfile), or "historical" (the reference reproduction only). */
int FillPhotonJet(const char* inputs, const char* output, const char* profileName = "nominal", const char* selection = "",
                  const char* centralityWeightText = "") {
  using namespace photonjet;
  const std::vector<std::string> paths = inputPaths(inputs);
  const Profile profile = configuredProfile(profileName, fileSystem(paths.front()), selection);
  const std::vector<CentralityClass> classes = centralityClasses(profile.system);
  Histograms h(classes);
  EventTable table;
  const CentralityWeights weights = std::string(centralityWeightText).empty() ? CentralityWeights{} :
      centralityWeights(centralityWeightText, profile.system);
  std::unique_ptr<TFile> first(TFile::Open(paths.front().c_str(),"READ"));
  const bool retainTriggers=profile.configured && profile.system==System::auau && !isSimulation(*first);
  std::map<std::string,std::unique_ptr<Histograms>> triggerHists;
  std::map<std::string,EventTable> triggerTables;
  std::string nominalTrigger;
  if(retainTriggers) {
    nominalTrigger=SelectionJson::parse(profile.selectionJson).at("systems").at("auau").at("data_trigger").get<std::string>();
    for(const auto& key:{nominalTrigger}) triggerHists[key].reset(new Histograms(classes));
    Long64_t offset=0;
    for(const auto& key:{nominalTrigger}) for(size_t f=0;f<paths.size();++f) {
      fillFile(paths[f],static_cast<int>(f),profile,classes,*triggerHists.at(key),triggerTables[key],weights,key,offset);
      std::unique_ptr<TFile> source(TFile::Open(paths[f].c_str(),"READ"));
      offset+=static_cast<TTree*>(source->Get("eventTree"))->GetEntries();
    }
  } else {
    for (size_t f = 0; f < paths.size(); ++f) fillFile(paths[f], static_cast<int>(f), profile, classes, h, table, weights);
  }
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  if(retainTriggers) {
    triggerHists.at(nominalTrigger)->write();triggerTables.at(nominalTrigger).write();
    for(const auto& key:{nominalTrigger}) {
      const std::string folder=std::string("triggers/")+key;
      out->mkdir(folder.c_str(),"",true);out->GetDirectory(folder.c_str())->cd();
      triggerHists.at(key)->write();triggerTables.at(key).write();
    }
    out->cd();TObjString(nominalTrigger.c_str()).Write("data_trigger");
  } else { h.write();table.write(); }
  out->cd();
  TObjString(profile.name.c_str()).Write("profile");
  writeSelection(profile);
  out->Close();
  return 0;
}
