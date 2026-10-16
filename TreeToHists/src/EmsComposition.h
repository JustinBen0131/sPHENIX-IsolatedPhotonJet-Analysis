/* Coherent EM scale composition for a fully retained SIM event.
   Selection, leakage and spectra use EMS * original ET. Response/transfer
   kinematics use EMS * (original ET + the named EMR draw). Jet selection
   scales only the EM fraction; jet response first adds the named JER draw.
   This adapter owns coordinate propagation, not model or cone prescriptions.
   A caller must use the bound whole-selection evaluator at selection ET.
   Native public fills require the explicit EmsNativeSelection contract.
   Retained-record replay and scientific cone adoption remain separate gates. */
#pragma once
#include "FillPhotonJetResponse.C"
#include "MeasurementAccumulator.h"
#include "PhotonEnergySelection.h"

namespace photonjet { namespace ems {

template<class T> struct Values : std::vector<T> {
  using std::vector<T>::vector;
  size_t GetSize() const { return this->size(); }
};

inline void require(bool ok,const char* why) {
  if(!ok)throw std::runtime_error(std::string("EMS composition: ")+why);
}

// All arrays keep original event-local indices. Never construct from selected
// compact rows: incoming threshold migrants and draw identity would be lost.
struct Momenta {
  response::Varied selectionPhotons,responsePhotons,selectionJets,responseJets;
  Values<double> originalPhotons;
  response::EventMatches matches;
  std::uint64_t source,entry;
  std::string selectionJson;
  template<class P,class J,class TP,class TJ,class F>
  Momenta(const Profile& profile,const response::Variation& v,const P& photons,const J& jets,
      const TP& truthPhotons,const TJ& truthJets,const response::EventMatches& match,
      const F& jetEmFraction,std::uint64_t sourceKey,std::uint64_t entryIndex)
      :selectionPhotons(photons,1.,0.,v.seed,sourceKey,entryIndex,0),
       responsePhotons(photons,1.,0.,v.seed,sourceKey,entryIndex,0),
       selectionJets(jets,1.,0.,v.seed,sourceKey,entryIndex,1),
       responseJets(jets,1.,0.,v.seed,sourceKey,entryIndex,1),
       matches(match),source(sourceKey),entry(entryIndex),selectionJson(profile.selectionJson) {
    require(profile.simulation,"DATA must use its unchanged measurement path");
    require(profile.configured && !profile.selectionJson.empty(),"resolved selection binding is required");
    require(response::emsSelectionResponse(v),"explicit ems_selection_response mode is required");
    for(size_t i=0;i<photons.GetSize();++i)originalPhotons.push_back(photons[i]);
    // Revalidate programmatically constructed settings with the same parser.
    const auto checked=response::parseVariation(("em_scale="+number(v.emScale)+
      " photon_resolution="+v.photonResolution+" jet_resolution="+v.jetResolution+
      " photon_scale="+number(v.photonScale)+" jet_scale="+number(v.jetScale)+
      " photon_smear="+number(v.photonSmear)+" jet_smear="+number(v.jetSmear)+
      " seed="+std::to_string(v.seed)+" photon_selection="+v.photonSelection).c_str());
    selectionPhotons.scaleAll(checked.emScale);
    responsePhotons.smearPhotons(truthPhotons,match.photonTruth,checked,source,entry);
    responsePhotons.scaleAll(checked.emScale);
    selectionJets.scaleEmFraction(jetEmFraction,checked.emScale);
    responseJets.smearJets(truthJets,match.jetTruth,checked,profile.jetRadius,source,entry);
    responseJets.scaleEmFraction(jetEmFraction,checked.emScale);
  }
private:
  static std::string number(double value) {
    std::ostringstream s;s.precision(17);s<<value;return s.str();
  }
};

// Result of the existing source-bound PhotonEnergySelection evaluator and its
// explicit cone prescription. Neither a stored score nor an EMS-scaled ET is
// sufficient to create this state. Isolation arrays are R03 then R04.
struct PhotonState {
  double score=std::numeric_limits<double>::quiet_NaN();
  int preselection=0,state=0;
  std::array<double,2> cone{},threshold{},sideband{};
  std::array<int,2> method{},role{};
};

struct Selection {
  Values<double> et,score,cone,threshold,sideband;
  Values<int> state,method,role,idClass,isoClass;
  std::vector<int> preselection;
  std::uint64_t source=0,entry=0;
  std::string selectionJson,evaluatorIdentity,conePrescriptionIdentity;

  void check(const Profile& profile,const Momenta& momenta) const {
    require(profile.simulation && profile.configured,"resolved SIM selection required");
    require(selectionJson==profile.selectionJson && selectionJson==momenta.selectionJson,
            "selection configuration differs");
    require(source==momenta.source && entry==momenta.entry,"original event identity differs");
    require(!evaluatorIdentity.empty() && !conePrescriptionIdentity.empty(),"evaluator/cone identity is missing");
    const auto n=et.size();
    require(n==momenta.selectionPhotons.GetSize() && score.size()==n && state.size()==n &&
      preselection.size()==n && idClass.size()==n && isoClass.size()==n && cone.size()==2*n &&
      threshold.size()==2*n && sideband.size()==2*n && method.size()==2*n && role.size()==2*n,
      "evaluated photon arrays are not aligned");
    for(size_t i=0;i<n;++i)require(et[i]==momenta.selectionPhotons[i],"selection ET differs from EMS coordinate");
  }
};

// evaluate(index, originalET, selectionET) must adapt the existing bound
// model/domain/NPB evaluator and explicit cone prescription. Resolution ET is
// intentionally absent from the callback, preventing EMR-driven classification.
// Identity strings record provenance; they do not verify model bytes. The
// source-bound evaluator constructor retains responsibility for SHA256 checks.
template<class P,class Evaluate>
inline Selection evaluateSelection(const Profile& profile,const Momenta& momenta,const P& original,
    const std::string& evaluatorIdentity,const std::string& coneIdentity,const Evaluate& evaluate) {
  require(original.GetSize()==momenta.selectionPhotons.GetSize(),"original photon count differs");
  require(!evaluatorIdentity.empty() && !coneIdentity.empty(),"evaluator/cone identity is required");
  Selection out;out.source=momenta.source;out.entry=momenta.entry;
  out.selectionJson=profile.selectionJson;out.evaluatorIdentity=evaluatorIdentity;out.conePrescriptionIdentity=coneIdentity;
  for(size_t i=0;i<original.GetSize();++i) {
    const double et=momenta.selectionPhotons[i];
    require(std::isfinite(original[i]) && std::isfinite(et),"photon ET is nonfinite");
    require(original[i]==momenta.originalPhotons[i],"original photon ET differs from realization");
    const PhotonState p=evaluate(i,original[i],et);
    require((p.preselection==0||p.preselection==1) && (p.state==0||p.state==1),"invalid evaluator state");
    require(p.state==0 || (std::isfinite(p.score)&&p.score>=0.&&p.score<=1.),"invalid evaluated score");
    out.et.push_back(et);out.score.push_back(p.score);out.preselection.push_back(p.preselection);out.state.push_back(p.state);
    out.idClass.push_back(0);out.isoClass.push_back(0); // Configured PhotonJet recomputes classes.
    for(size_t r=0;r<2;++r) {
      require(std::isfinite(p.cone[r]) && std::isfinite(p.threshold[r]) && std::isfinite(p.sideband[r]) &&
              p.sideband[r]>=p.threshold[r],"invalid prescribed isolation witness");
      out.cone.push_back(p.cone[r]);out.threshold.push_back(p.threshold[r]);out.sideband.push_back(p.sideband[r]);
      out.method.push_back(p.method[r]);out.role.push_back(p.role[r]);
    }
  }
  out.check(profile,momenta);return out;
}

// Reuse the established NPB/domain/selected-model implementation. Its
// constructor verifies every model SHA256 and the selected feature contract.
// The original-tree caller still owns exact captured-input extraction and
// the source-specific cone prescription; no fixed cone is a default here.
template<class P,class Inputs,class ConePrescription>
inline Selection evaluateBoundSelection(const Profile& profile,const Momenta& momenta,const P& original,
    const Inputs& captured,const energy_selection::Evaluator& evaluator,
    const std::string& evaluatorIdentity,const std::string& coneIdentity,const ConePrescription& prescribedCone) {
  require(captured.size()==original.GetSize(),"captured evaluator inputs are not aligned");
  return evaluateSelection(profile,momenta,original,evaluatorIdentity,coneIdentity,
    [&](size_t i,double originalEt,double selectionEt) {
      auto input=captured[i];
      require(input.et==originalEt,"captured evaluator ET differs from original photon");
      input.prescribedCone=prescribedCone(input,selectionEt);
      const auto result=evaluator.evaluate(input,selectionEt);
      PhotonState state;state.score=result.score;state.state=result.scoreState;state.preselection=result.preselection;
      state.cone=input.prescribedCone;state.threshold=result.threshold;state.sideband=result.sideband;
      state.method=input.method;state.role=input.role;return state;
    });
}

// Inclusive photon/ABCD and recoil spectra use only the selection coordinates.
// Named EMR/JER are reserved for response and both transfer components below.
template<class F>
inline void accumulateSpectra(const Profile& profile,const CentralityClass& c,
    int run,Long64_t seq,Long64_t entry,int fileIndex,double centrality,double vz,double weight,
    Histograms& histograms,EventTable* table,const Momenta& momenta,const Selection& selected,
    const F& eta,const F& phi,const F& jetEta,const F& jetPhi) {
  selected.check(profile,momenta);
  require(entry>=0 && static_cast<std::uint64_t>(entry)==momenta.entry,"spectra entry differs");
  std::vector<int> visible(selected.et.size(),1);
  accumulateMeasurement(profile,c,run,seq,entry,fileIndex,centrality,vz,weight,histograms,table,
    selected.et,eta,phi,selected.score,selected.preselection,visible,selected.state,selected.cone,
    selected.threshold,selected.sideband,selected.method,selected.role,selected.idClass,selected.isoClass,
    momenta.selectionJets,jetEta,jetPhi,1.);
}

// Response and numerator/denominator of the simulation transfer use one shared
// varied photon/jet realization; leakage uses the evaluated selection ET.
template<class F,class I>
inline void accumulateResponse(const Profile& profile,const std::vector<CentralityClass>& classes,
    const response::ResponseGrids& grids,std::vector<std::vector<response::Accumulators>>& histograms,
    double vz,double centrality,double weight,const Momenta& momenta,const Selection& selected,
    const F& eta,const F& phi,const F& jetEta,const F& jetPhi,
    const F& truthPt,const F& truthEta,const F& truthPhi,const I& truthRole,
    const F& truthJetPt,const F& truthJetEta,const F& truthJetPhi,const response::EventMatches& match) {
  selected.check(profile,momenta);
  require(match.photonTruth.size()==selected.et.size() && match.jetTruth.size()==momenta.responseJets.GetSize(),
          "response matches are not aligned");
  require(match.photonTruth==momenta.matches.photonTruth && match.photonMetric==momenta.matches.photonMetric &&
          match.jetTruth==momenta.matches.jetTruth,"original response matches differ");
  response::accumulateResponse(profile,classes,grids,histograms,vz,centrality,weight,momenta.responsePhotons,
    eta,phi,selected.score,selected.preselection,selected.state,selected.cone,selected.threshold,selected.sideband,
    selected.method,selected.role,selected.idClass,selected.isoClass,momenta.responseJets,jetEta,jetPhi,
    truthPt,truthEta,truthPhi,truthRole,truthJetPt,truthJetEta,truthJetPhi,match,&selected.et);
}
} }
