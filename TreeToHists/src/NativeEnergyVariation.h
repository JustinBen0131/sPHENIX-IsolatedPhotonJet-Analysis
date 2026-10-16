/* Original-tree adapter for the shared energy variation arithmetic.
   Matching remains ResponseMatching's exact source-index calculation. */
#pragma once
#include "EnergyVariation.h"
#include "CandidateCache.h"
namespace photonjet { namespace energy_selection {
struct ChangedMomenta {
  response::Varied photons,jets;
  template<class P,class J>ChangedMomenta(const P& p,const J& j,const response::Variation& v,std::uint64_t source,std::uint64_t entry)
    :photons(p,v.photonResolution=="off"?v.photonScale*v.emScale:1.,v.photonSmear,v.seed,source,entry,0),
     jets(j,v.jetResolution=="off"?v.jetScale:1.,v.jetSmear,v.seed,source,entry,1){}
};
inline bool changesPhotons(const response::Variation& v){return v.photonScale!=1.||v.emScale!=1.||v.photonSmear!=0.||v.photonResolution!="off";}
inline bool needsTruth(const response::Variation& v){return v.photonResolution!="off"||v.jetResolution!="off";}

class NativeVariation {
  using Array=candidate_cache::Array;
  std::map<std::string,std::unique_ptr<Array>> arrays_;
  int radius_;
  const Array& get(const std::string& name)const {
    const auto it=arrays_.find(name);if(it==arrays_.end()||!it->second)throw std::runtime_error("energy variation missing array: "+name);
    return *it->second;
  }
  struct Floats {const Array& a;size_t GetSize()const{return a.size();}double operator[](size_t i)const{return a.at(i);}};
  struct Ints {const Array& a;size_t GetSize()const{return a.size();}int operator[](size_t i)const {
    const double x=a.at(i);if(!std::isfinite(x)||x!=std::floor(x)||x<std::numeric_limits<int>::min()||x>std::numeric_limits<int>::max())
      throw std::runtime_error("inexact energy-variation match index");return static_cast<int>(x);}};
public:
  NativeVariation(TTree& tree,TTreeReader& reader,const Profile& profile,const response::Variation& v):radius_(jetRadiusIndex(profile.jetRadius)) {
    if(needsTruth(v)&&!profile.simulation)throw std::runtime_error("truth-based energy resolution requires SIM");
    std::vector<std::string> names;
    if(v.emScale!=1.)names.push_back(profile.jetBranch("jet","emcal_frac"));
    if(needsTruth(v)) {
      for(const char* kind:{"jet","jetunsub","truthjet"})for(const char* radius:{"02","03","04"})names.push_back(std::string(kind)+radius+"_pt");
      for(const char* field:{"rawpt"})for(const char* kind:{"jet","jetunsub"})names.push_back(profile.jetBranch(kind,field));
      for(const char* name:{"truthpho_pt","plink_reco","plink_truth","plink_class","plink_metric","jlink_reco","jlink_truth","jlink_class","jlink_metric"})names.push_back(name);
    }
    for(const auto& name:names)arrays_[name]=candidate_cache::reader(tree,reader,name);
  }
  template<class P,class J>ChangedMomenta apply(const Profile& profile,const P& photons,const J& jets,
        const response::Variation& v,std::uint64_t source,std::uint64_t entry)const {
    ChangedMomenta result(photons,jets,v,source,entry);
    if(needsTruth(v)) {
      const auto counts=[&](const char* kind){std::array<size_t,3> out{};int r=0;
        for(const char* suffix:{"02","03","04"})out[r++]=get(std::string(kind)+suffix+"_pt").size();return out;};
      const auto ranges=response::JetRanges::of(counts("jet"),counts("jetunsub"),counts("truthjet"),radius_);
      const auto raw=Floats{get(profile.jetBranch("jet","rawpt"))},unsub=Floats{get(profile.jetBranch("jetunsub","rawpt"))};
      if(raw.GetSize()!=ranges.nJet||unsub.GetSize()!=ranges.nUnsub||jets.GetSize()!=ranges.nJet)throw std::runtime_error("energy variation jet counts differ");
      for(const char* name:{"plink_truth","plink_class","plink_metric"})if(get(name).size()!=get("plink_reco").size())throw std::runtime_error("photon match arrays differ");
      for(const char* name:{"jlink_truth","jlink_class","jlink_metric"})if(get(name).size()!=get("jlink_reco").size())throw std::runtime_error("jet match arrays differ");
      const auto match=response::eventMatches(photons.GetSize(),get("truthpho_pt").size(),Ints{get("plink_reco")},Ints{get("plink_truth")},Ints{get("plink_class")},Floats{get("plink_metric")},get("plink_reco").size(),
          ranges,raw,unsub,Ints{get("jlink_reco")},Ints{get("jlink_truth")},Ints{get("jlink_class")},Floats{get("jlink_metric")},get("jlink_reco").size());
      if(v.photonResolution!="off") {result.photons.smearPhotons(Floats{get("truthpho_pt")},match.photonTruth,v,source,entry);result.photons.scaleAll(v.photonScale*v.emScale);}
      if(v.jetResolution!="off") {result.jets.smearJets(Floats{get(profile.jetBranch("truthjet","pt"))},match.jetTruth,v,profile.jetRadius,source,entry);result.jets.scaleAll(v.jetScale);}
    }
    if(v.emScale!=1.)result.jets.scaleEmFraction(Floats{get(profile.jetBranch("jet","emcal_frac"))},v.emScale);
    return result;
  }
};
} }
