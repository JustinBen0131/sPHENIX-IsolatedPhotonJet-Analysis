#pragma once
/* Physics inputs used by FinalAnalysis's background-correlation and purity
   non-closure calculations. No diagnostic shapes, timing or raw isolation.
   Axes, float storage and fill order match the established native QA inputs. */
#include "PhotonJet.h"
#include <TH2F.h>
#include <TH3F.h>
#include <memory>
#include <map>

namespace photonjet { namespace purity_inputs {
class Book {
 public:
  const std::string family;
  Book(const std::string& kind,const std::vector<CentralityClass>& classes):family(kind) {
    if(kind!="data" && kind!="photonjet" && kind!="inclusivejet")
      throw std::runtime_error("unknown purity-input family");
    const std::vector<double> region{0.,1.,2.};
    for(const auto& c:classes)for(const auto& label:labels()) {
      const auto tag=label+"_"+c.name;
      auto a=std::make_unique<TH2F>(("abcd_"+tag).c_str(),"event-leading photons;region (A,B,C,D);E_{T} [GeV]",4,0.,4.,20,15.,35.);
      auto r=std::make_unique<TH3F>(("recoil_"+tag).c_str(),"recoil jets (5 GeV);region (A,C);photon E_{T} [GeV];x_{J#gamma}",
          2,region.data(),photonEtEdges().size()-1,photonEtEdges().data(),xjEdges().size()-1,xjEdges().data());
      a->SetDirectory(nullptr);r->SetDirectory(nullptr);a->Sumw2();r->Sumw2();
      abcd_.emplace(tag,std::move(a));recoil_.emplace(tag,std::move(r));
    }
  }
  std::vector<std::string> labels()const {
    if(family=="data")return {"data"};
    return family=="photonjet"?std::vector<std::string>{"sig","phoother"}:std::vector<std::string>{"bkg","jetprompt"};
  }
  std::string label(bool signal)const {
    if(family=="data")return "data";
    return family=="photonjet"?(signal?"sig":"phoother"):(signal?"jetprompt":"bkg");
  }
  template<class Values,class Ints,class Scores,class Preselection>
  void accumulate(const Profile& profile,const CentralityClass& centralityClass,double vz,double centrality,double weight,
      const Values& et,const Values& eta,const Values& phi,const Scores& score,const Preselection& accepted,
      const Ints& state,const Values& cone,const Values& threshold,const Values& sideband,const Ints& method,const Ints& role,
      const Ints& id,const Ints& iso,const Values& jetPt,const Values& jetEta,const Values& jetPhi,
      const std::vector<int>& signal,double jetScale) {
    if(signal.size()!=et.GetSize() || (family=="data")==profile.simulation)
      throw std::runtime_error("purity label/profile cardinality differs");
    const auto leaders=regionLeaders(profile,vz,centrality,et.GetSize(),et,eta,phi,score,accepted,state,cone,threshold,sideband,method,role,id,iso);
    Profile selected=profile;selected.jetPtMin=kJetThresholds[0];
    for(int region=0;region<4;++region) {
      const auto i=leaders[region];if(i<0)continue;
      const auto tag=label(signal.at(i))+"_"+centralityClass.name;
      abcd_.at(tag)->Fill(region+.5,et[i],weight);
      if(region!=kA && region!=kC)continue;
      for(size_t j=0;j<jetPt.GetSize();++j) {
        const double pt=jetPt[j]*jetScale,xj=pt/et[i];
        if(!(xj>=0. && xj<3.) || !recoilJet(selected,eta[i],phi[i],pt,jetEta[j],jetPhi[j]))continue;
        recoil_.at(tag)->Fill(region==kA?.5:1.5,et[i],xj,weight);
      }
    }
  }
  void write()const {for(const auto& [_,h]:abcd_)h->Write();for(const auto& [_,h]:recoil_)h->Write();}
 private:
  std::map<std::string,std::unique_ptr<TH2F>> abcd_;
  std::map<std::string,std::unique_ptr<TH3F>> recoil_;
};
} }
