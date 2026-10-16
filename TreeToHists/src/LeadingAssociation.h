#ifndef PHOTONJET_LEADING_ASSOCIATION_H
#define PHOTONJET_LEADING_ASSOCIATION_H
#include "PhotonJet.h"
namespace photonjet {
inline std::vector<double> associationEtEdges(){std::vector<double> v;for(int i=10;i<=40;++i)v.push_back(i);return v;}
inline std::vector<double> associationXjEdges(){std::vector<double> v;for(int i=0;i<=60;++i)v.push_back(i*.05);return v;}
// RecoilJets.cc runLeadIsoTightPhotonJetLoopAllRadii and fillRecoTruthJES3MatchingQA:
// highest pT after finite/pT/overlap, then eta and recoil veto, never fallback.
// Global refers to the stored collection; storage acceptance remains a separate gate.
template<class Pt,class Eta,class Phi>
long associationLeader(const Profile&p,double pe,double pp,const Pt&pt,const Eta&eta,const Phi&phi){
 long best=-1;
 for(size_t j=0;j<pt.GetSize();++j){
  if(!std::isfinite(pt[j])||!std::isfinite(eta[j])||!std::isfinite(phi[j])||pt[j]<p.jetPtMin||deltaR(eta[j],phi[j],pe,pp)<p.overlapDeltaR)continue;
  if(best<0||pt[j]>pt[best])best=j;
 }
 if(best<0||!truthRecoilJet(p,pe,pp,pt[best],eta[best],phi[best]))return -1;
 return best;
}
inline bool dataTimeCut(double time){return std::isfinite(time)&&std::fabs(time+2.0)<6.0;}
}
#endif
