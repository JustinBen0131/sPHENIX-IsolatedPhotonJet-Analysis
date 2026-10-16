#pragma once
#include "PhotonJet.h"
#include "LeadingAssociation.h"
#include "MeasurementScientificDefinitions.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>
#include <map>
#include <memory>
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
        const auto ae=associationEtEdges(), ax=associationXjEdges();
        add2(new TH2D(("leading_global_"+tag).c_str(),"global leading jet;reco photon E_{T};reco xJ",ae.size()-1,ae.data(),ax.size()-1,ax.data()));
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
  void scientificContext(const photonjet_science::Json& context) {
    for(auto& [name,h]:h1_)photonjet_science::observeMeasurementDefinition(*h,context);
    for(auto& [name,h]:h2_)photonjet_science::observeMeasurementDefinition(*h,context);
    for(auto& [name,h]:h3_)photonjet_science::observeMeasurementDefinition(*h,context);
  }
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
  std::map<std::string,std::set<std::string>> scientificDefinitions;
  void scientificContext(const photonjet_science::Json& context) {
    using namespace photonjet_science;
    for(const std::string name:{"events","photons","jets"}) {
      // Unlike events_all counts, these events require a region leader.
      auto d=measurementDefinition(name=="jets"?"dphiA_all_jet5_near":"abcd_all",context);
      d["observable"]="records/"+name;
      d["terms"].push_back(term("projection","record_retention",name=="events"?
        Json("accepted event with at least one region leader"):name=="photons"?
        Json("one highest-ET photon per ABCD region"):Json("jets around region A/C leaders; all angles; threshold bitmask")));
      if(name=="jets") {
        auto& terms=d["terms"];
        terms.erase(std::remove_if(terms.begin(),terms.end(),[](const Json& x){return x.at("id")=="eta_side"||x.at("id")=="jet_threshold_gev";}),terms.end());
        terms.push_back(term("selection","jet_thresholds_gev",photonjet::kJetThresholds));
      }
      scientificDefinitions[name].insert(d.dump());
    }
  }
  void scientificTag(TTree& tree,const std::string& name) const {
    using namespace photonjet_science;
    Json definitions=Json::array();const auto found=scientificDefinitions.find(name);
    if(found!=scientificDefinitions.end())for(const auto& value:found->second)definitions.push_back(Json::parse(value));
    if(definitions.empty())definitions.push_back(nullptr);
    tree.GetUserInfo()->Add(new TNamed("physics_effective_definitions",definitions.dump().c_str()));
  }
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
    scientificTag(events,"events");
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
    scientificTag(photons,"photons");
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
    scientificTag(jets,"jets");
    jets.Write("", TObject::kOverwrite);  /* one key: an autosave cycle is replaced */
  }
};

/* Shared measurement arithmetic for native trees and candidate records. */
template<class Values, class Ints, class Scores, class Energies, class States, class WitnessValues, class JetMomenta>
inline void accumulateMeasurement(const Profile& profile, const CentralityClass& c,
    int run, Long64_t seq, Long64_t entry, int fileIndex, double centrality, double vz, double w,
    Histograms& h, EventTable* saved, const Energies& et, const Values& eta, const Values& phi,
    const Scores& score, const std::vector<int>& accepted, const std::vector<int>& visible,
    const States& state, const WitnessValues& cone, const WitnessValues& threshold, const WitnessValues& sideband,
    const Ints& method, const Ints& role, const Ints& idClass, const Ints& isoClass,
    const JetMomenta& storedJetPt, const Values& jetEta, const Values& jetPhi, double jetScale, const std::array<long,4>* sharedLeaders=nullptr) {
    const size_t nPhoton=et.GetSize();
    EventTable discarded;
    EventTable& table=saved?*saved:discarded;
    TH1D& kinematic = h.h1("photons_" + c.name);
    for (size_t i = 0; i < nPhoton; ++i)
      if (visible[i] && photonKinematics(profile, vz, et[i], eta[i], phi[i])) kinematic.Fill(et[i], w);
    const std::array<long, 4> leader = sharedLeaders ? *sharedLeaders : regionLeaders(profile, vz, centrality, nPhoton, et, eta, phi, score, accepted, state,
                                                     cone, threshold, sideband, method, role, idClass, isoClass);
    if (leader[kA] < 0 && leader[kB] < 0 && leader[kC] < 0 && leader[kD] < 0) return;
    TH2D& abcd = h.h2("abcd_" + c.name);
    const int event = static_cast<int>(table.run.size());
    table.run.push_back(run);
    table.seq.push_back(seq);
    table.fileIndex.push_back(fileIndex);
    table.entry.push_back(entry);
    table.centrality.push_back(centrality);
    table.slice.push_back(profile.system == System::auau ? centralitySlice(centrality) : -1);
    table.weight.push_back(w);
    for (int r = 0; r < 4; ++r) {
      if (leader[r] < 0) continue;
      abcd.Fill(r + 0.5, et[leader[r]], w);
      table.photonEvent.push_back(event);
      table.photonRegion.push_back(r);
      table.photonEt.push_back(et[leader[r]]);
      table.photonEtBin.push_back(abcd.GetYaxis()->FindFixBin(et[leader[r]]));
    }
    if(leader[kA]>=0){const size_t a=leader[kA];const ScaledValues pts(storedJetPt,jetScale);
      for(double threshold:kJetThresholds){Profile p=profile;p.jetPtMin=threshold;
        const long j=associationLeader(p,eta[a],phi[a],pts,jetEta,jetPhi);
        if(j>=0)h.h2("leading_global_"+c.name+"_"+thresholdName(threshold)).Fill(et[a],pts[j]/et[a],w);
      }
    }
    for (int r : {kA, kC}) {
      if (leader[r] < 0) continue;
      const size_t i = leader[r];
      const std::string name = r == kA ? "A" : "C";
      const ScaledValues jetPt(storedJetPt, jetScale);
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
} // namespace photonjet
