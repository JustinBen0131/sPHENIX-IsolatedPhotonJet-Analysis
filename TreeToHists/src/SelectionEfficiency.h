/* The stage efficiency arithmetic from FillPhotonQA.C, on the same native
   score, preselection and subtracted isolation witnesses as the recoil fill.
   Photon-family SIM only. Diagnostics and physics correction are distinct. */
#pragma once
#include "PhotonJet.h"
#include <TH3F.h>
namespace photonjet { namespace selection_efficiency {
class Book {
 std::map<std::string,std::unique_ptr<TH3F>> histograms_;
 public:
 explicit Book(const std::vector<CentralityClass>& classes){
  for(const auto& c:classes){auto h=std::make_unique<TH3F>(("eff_"+c.name).c_str(),
   "analysis-signal truth photons;p_{T}^{truth} [GeV];stage;matched reco E_{T} [GeV]",35,5.,40.,6,0.,6.,21,-2.5,50.);
   h->SetDirectory(nullptr);h->Sumw2();histograms_[c.name]=std::move(h);}
 }
 template<class Values,class Ints,class Scores,class Preselection>
 void accumulate(const Profile& profile,const CentralityClass& c,double vz,double centrality,double w,
  const Values& et,const Values& eta,const Values& phi,const Scores& score,const Preselection& presel,
  const Ints& state,const Values& cone,const Values& threshold,const Values& sideband,
  const Ints& method,const Ints& role,const Ints& id,const Ints& iso,const Values& truthPt,const Values& truthEta,
  const Ints& truthRole,const std::vector<long>& truthOf){
  if(!profile.simulation || truthOf.size()!=et.GetSize())throw std::runtime_error("efficiency requires aligned SIM");
  auto& h=*histograms_.at(c.name);
  for(size_t t=0;t<truthPt.GetSize();++t){
   if(!std::isfinite(truthPt[t]) || !std::isfinite(truthEta[t]) || truthRole[t]!=kSignalSourceRole)continue;
   h.Fill(truthPt[t],.5,-1.,w);
   long r=-1;
   for(size_t i=0;i<et.GetSize();++i)if(truthOf[i]==long(t) && std::isfinite(vz) && std::isfinite(et[i]) &&
    std::isfinite(eta[i]) && std::isfinite(phi[i]) && std::abs(vz)<profile.vertexAbsZMax && std::abs(eta[i])<profile.photonAbsEtaMax)r=i;
   if(r<0)continue;
   const auto d=photonClass(profile,r,et[r],centrality,score,state,cone,threshold,sideband,method,role,id,iso);
   const bool pre=presel[r]==1;
   const std::array<bool,5> reached{true,pre,pre&&d.isolated,pre&&d.tight,pre&&d.tight&&d.isolated};
   for(int stage=0;stage<5;++stage)if(reached[stage])h.Fill(truthPt[t],stage+1.5,et[r],w);
  }
 }
 void write()const{for(const auto& [_,h]:histograms_)h->Write();}
};
} }
