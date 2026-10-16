#include "ScientificProducerEvidence.h"
#include "SelectionViewContract.h"
#include "PurityInputs.h"
#include "SelectionEfficiency.h"
#pragma once
/* Simulation response reduction from eventTree v1.

   Fills photon and joint photon-pT/xJ responses, misses, categorized fakes,
   leakage and combinatoric templates on the plan's native axes.
   PhotonJet.h owns selection, truth matching, source admission and weights.
   Sumw2 and source identity support independently resampled responses.
   Public execution is tree_to_hists.sh; the detailed response contract
   below defines categories and saved grids. */
#include "PhotonJet.h"
#include "ResponseMatching.h"
#include "LeadingAssociation.h"
#include "MeasurementAccumulator.h"
#include <THnSparse.h>

#ifdef __CLING__
R__LOAD_LIBRARY(libTreePlayer)
#endif

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TObjString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/* Fill the photon + recoil-jet response inputs from eventTree v1 simulation files. */

/* root -l -q 'TreeToHists/src/FillPhotonJetResponse.C("a.eventtree.root b.eventtree.root", "out.root", "nominal")' */

/* Analysis events (evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && an accepted producer class,
   acceptedEventClass in PhotonJet.h) whose vertex
   is finite and inside the profile window. The reco photon is the region-A event-leading photon (PhotonJet.h),
   the truth photon and the truth recoil jets come from PhotonJet.h, and the two are paired by the reco-truth
   matches (plink_, jlink_; see eventMatches: jets use the match inside the analysis jet view, rebuilt from the
   producer's candidate edges). Per centrality class <c> (Au+Au: cent0_20, cent20_50, cent50_80; pp: all) and jet-pT
   threshold <t> (jet5, jet7, jet10, jet12 GeV), every name below carries the suffix _<c>_<t>:
   photon_reco                             TH1D  reco pT                  region-A photon
   photon_truth                            TH1D  truth pT                 truth photon
   photon_response_truth_x_reco            TH2D  truth pT x reco pT       matched photons
   photon_fakes, photon_misses             TH1D  reco pT, truth pT        unmatched reco, truth photons
   xj_reco, xj_truth                       TH2D  pT x xJ                  reco, truth recoil jets
   xj_response_truth_global_x_reco_global  TH2D  truth cell x reco cell   matched recoil jets
   xj_fakes                                TH2D  reco pT x xJ             reco recoil jets without a match, split into
   xj_fake_photon_reco, xj_combinatoric_reco, xj_detector_fakes_reco      (RecoilClass in PhotonJet.h)
   xj_misses                               TH2D  truth pT x xJ            truth recoil jets without a match
   xj_leading_reco                         TH2D  reco pT x xJ             the leading reco recoil jet of the region-A photon
   xj_leading_matched_reco                 TH2D  reco pT x xJ             the same, when the photon is matched and the
   jet's match is the leading truth recoil jet of the truth photon (the leading-jet association of IAN Fig. 25)
   signal_leakage_truth_matched_abcd       TH2D  region x reco pT         one region leader per signal truth photon
   dphi_truthjet_<side>                    TH3D  reco pT x xJ x |dphi|    every jet of a matched region-A photon passing
   dphi_combinatoric_<side>                                               the recoil cuts except dphi, with and without
   a truth jet nearby (nearbyTruthJet); side near
   (|deta| <= 0.8) or far
   The dphi fills are the simulation closure of the data-driven combinatoric estimate (FinalAnalysis combinatoric.py):
   their sum is what FillPhotonJet.C fills on data, and the recoil-window integral of dphi_combinatoric is
   xj_combinatoric_reco.
   Axes (arguments; the defaults are the reference analysis's): reco pT 10 15 20 25 35 40 GeV; truth pT 5 10 15 20 25
   35 40 GeV; xJ in 0.05 steps to 2, then 3; dphi deltaPhiEdges() (PhotonJet.h). A response cell is pT bin *
   (xJ bins) + xJ bin (0-based). A value outside its axis is never filled (no under- or overflow) and decides the
   category as described in PhotonJet.h, so only the axis ranges change the categories. The axes used (dphi
   excepted) are written as the JSON object "grids". Weights are w_event
   times, in Au+Au, the centrality weight of the event's 5% slice (an input, CentralityWeights in PhotonJet.h; 1
   when the product is not finite), written as "centrality_weights". All histograms keep Sumw2. "sources" lists every
   input with its production input identity (meta/source input_file_sha256), which the half closure splits on.
   "variation" records the reco energy scale and resolution variation of a systematic run (Variation; "" is nominal). */

namespace photonjet {
namespace response {

/* The response axes. The response categories depend only on the outer range of each axis (a value outside it is
   never filled and decides the category as described in PhotonJet.h), not on the inner edges. So a response
   filled on fine axes sums exactly onto any coarser axes with the same ranges whose edges lie on the fine ones.
   The defaults are the reference analysis's axes. */
struct ResponseGrids {
  std::vector<double> reco{10.0, 15.0, 20.0, 25.0, 35.0, 40.0};
  std::vector<double> truth{5.0, 10.0, 15.0, 20.0, 25.0, 35.0, 40.0};
  std::vector<double> xj = xjEdges();
};

/* Whitespace-separated, strictly increasing edges; an empty string keeps the default. */
inline std::vector<double> parseEdges(const char* text, const std::vector<double>& fallback) {
  std::vector<double> edges;
  std::istringstream stream(text ? text : "");
  for (double value; stream >> value;) edges.push_back(value);
  if (edges.empty()) return fallback;
  if (edges.size() < 2) throw std::runtime_error("an axis needs at least two edges");
  for (size_t i = 1; i < edges.size(); ++i)
    if (!(edges[i] > edges[i - 1])) throw std::runtime_error("axis edges must be strictly increasing");
  return edges;
}

inline std::string jsonArray(const std::vector<double>& edges) {
  std::ostringstream s;
  s.precision(17);
  s << '[';
  for (size_t i = 0; i < edges.size(); ++i) s << (i ? "," : "") << edges[i];
  return s.str() + ']';
}

/* A variation of the reco energy scale and resolution, for the systematic uncertainties: every reco jet of the analysis collection and
   photon of the event is taken at value x scale x (1 + smear g) before any selection, g a standard Gaussian drawn per
   object from a counter-based generator keyed by (seed, the file's production input, entry, object), so the draws of
   one file do not depend on the other inputs or their order. The reco-truth jet match stays the producer's, and the
   photon classes and isolation thresholds are the stored ones. The text is "" (nominal) or whitespace-separated
   key=value pairs of jet_scale, jet_smear, photon_scale, photon_smear, seed. */
} } // photonjet::response
#include "EnergyVariation.h"
#include "EmsNativeSelection.h"
namespace photonjet { namespace response {

/* Bin of value on edges, -1 when it is not finite or outside [front, back). */
inline int binIndex(double value, const std::vector<double>& edges) {
  if (!std::isfinite(value) || value < edges.front() || value >= edges.back()) return -1;
  for (size_t i = 0; i + 1 < edges.size(); ++i)
    if (edges[i] <= value && value < edges[i + 1]) return static_cast<int>(i);
  return -1;
}

/* The JSON object of meta/collections whose name is `name` (compact JSON, Format tables.dumps). */
inline std::string collectionEntry(const std::string& collections, const std::string& name) {
  const size_t at = collections.find("\"name\":\"" + name + "\"");
  if (at == std::string::npos) return "";
  const size_t open = collections.rfind('{', at), close = collections.find('}', at);
  return open == std::string::npos || close == std::string::npos ? "" : collections.substr(open, close - open + 1);
}

/* System of a simulation file. Fails unless it is simulation with anti-kT reco and truth jet collections, the reco
   collection of the configured radius (jet04 at R = 0.4) being the analysis view of the system. */
inline System simulationSystem(TFile& f) {
  const System system = fileSystem(f);
  if (!isSimulation(f)) throw std::runtime_error(std::string("not a simulation file: ") + f.GetName());
  const std::string name = std::string("jet") + jetRadiusTag(configuredJetRadius);
  const std::string jets = collectionEntry(metaText(f, "meta/collections"), name);
  const std::string view = system == System::pp ? "towerinfo_calib" : "towerinfo_sub1";
  if (jets.find("\"algorithm\":\"antikt\"") == std::string::npos ||
      jets.find("\"input_identity\":\"" + view + "\"") == std::string::npos)
    throw std::runtime_error(name + " is not the anti-kT analysis view in " + f.GetName());
  if (metaText(f, "meta/constants").find("\"truthJets.algorithm\":\"antikt\"") == std::string::npos)
    throw std::runtime_error(std::string("truth jets are not anti-kT in ") + f.GetName());
  return system;
}

/* The accumulators of one centrality class and jet threshold. */
struct Accumulators {
  TH1D *photonReco, *photonTruth, *photonFakes, *photonMisses;
  TH2D *photonResponse, *leakage;
  TH2D *xjReco, *xjTruth, *xjResponse, *xjFakes, *xjFakePhoton, *xjDetectorFakes, *xjCombinatoric, *xjMisses;
  TH2D *xjLeading, *xjLeadingMatched, *xjIAnLeading, *xjIAnTagged;
  TH2D *xjLeadingUniform=nullptr,*xjLeadingMatchedUniform=nullptr,*xjLeadingAnyUniform=nullptr;
  std::unique_ptr<THnSparseD> pairReason;
  std::array<TH3D*, 2> dphiTruthJet{}, dphiCombinatoric{};  /* side 0 near, 1 far */
  std::vector<std::unique_ptr<TH1>> owned;  /* booking order, which is the write order */

  Accumulators(const std::string& suffix, const ResponseGrids& g, bool ppQa=false) {
    const std::vector<double>& xj = g.xj;
    const int nXj = xj.size() - 1, nReco = g.reco.size() - 1, nTruth = g.truth.size() - 1;
    const double* reco = g.reco.data();
    const double* truth = g.truth.data();
    const std::string recoAxes = ";reco p_{T}^{#gamma} [GeV];reco x_{J#gamma}";
    const std::string truthAxes = ";truth p_{T}^{#gamma} [GeV];truth x_{J#gamma}";
    auto h1 = [&](const std::string& key, const std::string& title, const std::vector<double>& edges) {
      return keep(new TH1D((key + suffix).c_str(), title.c_str(), edges.size() - 1, edges.data()));
    };
    auto h2 = [&](const std::string& key, const std::string& title, int nx, const double* x, int ny, const double* y) {
      return keep(new TH2D((key + suffix).c_str(), title.c_str(), nx, x, ny, y));
    };
    photonReco = h1("photon_reco", "region-A photons;reco p_{T}^{#gamma} [GeV]", g.reco);
    photonTruth = h1("photon_truth", "truth photons;truth p_{T}^{#gamma} [GeV]", g.truth);
    photonResponse = h2("photon_response_truth_x_reco",
                        "matched photons;truth p_{T}^{#gamma} [GeV];reco p_{T}^{#gamma} [GeV]", nTruth, truth, nReco, reco);
    photonFakes = h1("photon_fakes", "unmatched region-A photons;reco p_{T}^{#gamma} [GeV]", g.reco);
    photonMisses = h1("photon_misses", "unmatched truth photons;truth p_{T}^{#gamma} [GeV]", g.truth);
    xjReco = h2("xj_reco", "reco recoil jets" + recoAxes, nReco, reco, nXj, xj.data());
    xjTruth = h2("xj_truth", "truth recoil jets" + truthAxes, nTruth, truth, nXj, xj.data());
    std::vector<double> truthCells(nTruth * nXj + 1), recoCells(nReco * nXj + 1);
    for (size_t i = 0; i < truthCells.size(); ++i) truthCells[i] = i;
    for (size_t i = 0; i < recoCells.size(); ++i) recoCells[i] = i;
    xjResponse = h2("xj_response_truth_global_x_reco_global",
                    "matched recoil jets;truth cell (p_{T} bin #times 41 + x_{J} bin);"
                    "reco cell (p_{T} bin #times 41 + x_{J} bin)",
                    nTruth * nXj, truthCells.data(), nReco * nXj, recoCells.data());
    xjFakes = h2("xj_fakes", "unmatched reco recoil jets" + recoAxes, nReco, reco, nXj, xj.data());
    xjFakePhoton = h2("xj_fake_photon_reco", "reco recoil jets of an unmatched photon" + recoAxes, nReco, reco, nXj,
                      xj.data());
    xjDetectorFakes = h2("xj_detector_fakes_reco", "unmatched reco recoil jets with a truth jet nearby" + recoAxes,
                         nReco, reco, nXj, xj.data());
    xjMisses = h2("xj_misses", "unmatched truth recoil jets" + truthAxes, nTruth, truth, nXj, xj.data());
    xjCombinatoric = h2("xj_combinatoric_reco", "reco recoil jets without a truth jet nearby" + recoAxes, nReco, reco,
                        nXj, xj.data());
    xjLeading = h2("xj_leading_reco", "leading reco recoil jet" + recoAxes, nReco, reco, nXj, xj.data());
    xjLeadingMatched = h2("xj_leading_matched_reco", "leading reco recoil jet matched to the leading truth recoil jet" + recoAxes,
                          nReco, reco, nXj, xj.data());
    const auto ae=associationEtEdges(),ax=associationXjEdges();
    if(ppQa){
      xjLeadingUniform=h2("xj_leading_reco_u","QA only: leading reco recoil jet"+recoAxes,ae.size()-1,ae.data(),ax.size()-1,ax.data());
      xjLeadingMatchedUniform=h2("xj_leading_matched_reco_u","QA only: leading truth recoil match"+recoAxes,ae.size()-1,ae.data(),ax.size()-1,ax.data());
      xjLeadingAnyUniform=h2("xj_leading_matched_any_reco_u","QA only: any selected truth recoil match"+recoAxes,ae.size()-1,ae.data(),ax.size()-1,ax.data());
    }
    xjIAnLeading=h2("xj_leading_global_reco","global leading jet, then fiducial/recoil veto"+recoAxes,ae.size()-1,ae.data(),ax.size()-1,ax.data());
    xjIAnTagged=h2("xj_leading_global_matched_reco","matched photon and global leading jet dR <= 0.3"+recoAxes,ae.size()-1,ae.data(),ax.size()-1,ax.data());
    const int bins[3]={nReco,nXj,65536};const double lo[3]={g.reco.front(),xj.front(),-.5},hi[3]={g.reco.back(),xj.back(),65535.5};
    pairReason=std::make_unique<THnSparseD>(("reco_pair_truth_reason"+suffix).c_str(),"all selected reco pairs;reco photon E_{T};reco xJ;reason bitmask",3,bins,lo,hi);
    pairReason->GetAxis(0)->Set(nReco,reco);pairReason->GetAxis(1)->Set(nXj,xj.data());pairReason->Sumw2();
    const double regions[5] = {0, 1, 2, 3, 4};
    leakage = h2("signal_leakage_truth_matched_abcd",
                 "one region leader per signal truth photon;region (A,B,C,D);reco p_{T}^{#gamma} [GeV]", 4, regions,
                 nReco, reco);
    const std::vector<double> dphi = deltaPhiEdges();
    for (int side = 0; side < 2; ++side) {
      auto h3 = [&](const std::string& key, const std::string& title) {
        return keep(new TH3D((key + (side ? "_far" : "_near") + suffix).c_str(),
                             (title + ";reco p_{T}^{#gamma} [GeV];reco x_{J#gamma};|#Delta#phi|").c_str(), nReco, reco,
                             nXj, xj.data(), dphi.size() - 1, dphi.data()));
      };
      dphiTruthJet[side] = h3("dphi_truthjet", "jets of a matched photon with a truth jet nearby");
      dphiCombinatoric[side] = h3("dphi_combinatoric", "jets of a matched photon without a truth jet nearby");
    }
  }
  void write() {
    for (auto& h : owned) h->Write();
    pairReason->Write();
    if(!pairDefinition.is_null())TObjString(pairDefinition.dump().c_str()).Write((std::string(pairReason->GetName())+"__physics_definition").c_str());
  }
  void scientificContext(const photonjet_science::Json& context) {
    for(auto& h:owned)photonjet_science::observeMeasurementDefinition(*h,context,true);
    // Sparse reason-bit diagnostics are not yet consumed by the maintained
    // renderer. Keep their definition beside the sparse object for the writer.
    pairDefinition=photonjet_science::measurementDefinition(pairReason->GetName(),context,true);
  }
  photonjet_science::Json pairDefinition;

 private:
  template <class H>
  H* keep(H* h) {
    h->Sumw2();
    h->SetDirectory(nullptr);
    owned.emplace_back(h);
    return h;
  }
};

struct ResolutionMigration {
  std::uint64_t events=0,candidates=0,nominalAccepted=0,variedAccepted=0,incoming=0,outgoing=0,
                binChanged=0,leadersChanged=0,unavailableNominalId=0;
};
inline nlohmann::json resolutionMigrationJson(const ResolutionMigration& c) {
  return {{"events",c.events},{"finite_candidates",c.candidates},{"nominal_accepted",c.nominalAccepted},
    {"varied_accepted",c.variedAccepted},{"incoming",c.incoming},{"outgoing",c.outgoing},
    {"response_bin_changed",c.binChanged},{"region_leaders_changed",c.leadersChanged},
    {"preselected_without_nominal_id",c.unavailableNominalId}};
}

/* The event-level response accumulator is shared by native fills and cache
   replay. Inputs retain original photon indices so deterministic smearing
   has the same object identity after cache compaction. */
template<class Floats,class Ints,class Scores,class Preselection,class Energies,class JetMomenta,class States,class WitnessValues,class SelectionInts>
inline void accumulateResponse(const Profile& profile,const std::vector<CentralityClass>& classes,
    const ResponseGrids& grids,std::vector<std::vector<Accumulators>>& acc,double vz,double centrality,double w,
    const Energies& et,const Floats& eta,const Floats& phi,const Scores& score,const Preselection& presel,
    const States& state,const WitnessValues& cone,const WitnessValues& threshold,const WitnessValues& sideband,const SelectionInts& method,
    const SelectionInts& role,const SelectionInts& idClass,const SelectionInts& isoClass,const JetMomenta& jetPt,const Floats& jetEta,
    const Floats& jetPhi,const Floats& truthPt,const Floats& truthEta,const Floats& truthPhi,const Ints& truthRole,
    const Floats& truthJetPt,const Floats& truthJetEta,const Floats& truthJetPhi,const EventMatches& match,
    const std::vector<double>* classificationEt=nullptr,ResolutionMigration* migration=nullptr) {
  const auto& xj=grids.xj;const long nXj=xj.size()-1;
  const size_t nPhoton=et.GetSize(),nTruthJets=truthJetPt.GetSize();
  if(classificationEt && classificationEt->size()!=nPhoton)throw std::runtime_error("response classification energies are misaligned");
  const auto selectedLeaders=classificationEt
    ? regionLeadersWithClassification(profile,vz,centrality,nPhoton,et,*classificationEt,eta,phi,score,presel,state,cone,threshold,sideband,method,role,idClass,isoClass)
    : regionLeaders(profile,vz,centrality,nPhoton,et,eta,phi,score,presel,state,cone,threshold,sideband,method,role,idClass,isoClass);
  const auto purityLeaders=classificationEt
    ? regionLeaders(profile,vz,centrality,nPhoton,*classificationEt,eta,phi,score,presel,state,cone,threshold,sideband,method,role,idClass,isoClass)
    : selectedLeaders;
  if(classificationEt && migration) {
    ++migration->events;
    const auto nominalLeaders=regionLeaders(profile,vz,centrality,nPhoton,*classificationEt,eta,phi,score,presel,state,cone,threshold,sideband,method,role,idClass,isoClass);
    for(size_t r=0;r<4;++r)migration->leadersChanged+=nominalLeaders[r]!=selectedLeaders[r];
    for(size_t i=0;i<nPhoton;++i) {
      const double original=classificationEt->at(i);if(!std::isfinite(original))continue;
      ++migration->candidates;if(presel[i]!=1)continue;
      const auto c=photonClass(profile,i,original,centrality,score,state,cone,threshold,sideband,method,role,idClass,isoClass);
      migration->unavailableNominalId+=!c.scored;
      const bool nominal=c.region()!=kNone && photonKinematics(profile,vz,original,eta[i],phi[i]);
      const bool varied=c.region()!=kNone && photonKinematics(profile,vz,et[i],eta[i],phi[i]);
      migration->nominalAccepted+=nominal;migration->variedAccepted+=varied;
      migration->incoming+=varied&&!nominal;migration->outgoing+=nominal&&!varied;
      migration->binChanged+=nominal&&varied&&binIndex(original,grids.reco)!=binIndex(et[i],grids.reco);
    }
  }
    const long tp = truthPhoton(truthPt, truthEta, truthPhi, truthRole, truthPt.GetSize());
    const int truthBin = tp >= 0 ? binIndex(truthPt[tp], grids.truth) : -1;
    for (size_t k = 0; k < classes.size(); ++k) {
      if (profile.system == System::auau && !(classes[k].low <= centrality && centrality < classes[k].high)) continue;
      const std::array<long, 4> leader = selectedLeaders;
      const long a = leader[kA];
      const int recoBin = a >= 0 ? binIndex(et[a], grids.reco) : -1;
      const bool photonMatched = a >= 0 && tp >= 0 && match.photonTruth[a] == tp;
      const PhotonResponse photon = photonResponse(photonMatched, recoBin >= 0, truthBin >= 0);
      // This legacy response-path histogram is consumed by the ABCD correction.
      // Response-only resolution keeps its nominal leaders and unsmeared ET.
      const std::array<bool, 4> witness = classificationEt
          ? leakageWitnesses(purityLeaders,match.photonTruth,match.photonMetric,*classificationEt,truthRole)
          : leakageWitnesses(purityLeaders,match.photonTruth,match.photonMetric,et,truthRole);

      for (size_t t = 0; t < kJetThresholds.size(); ++t) {
        Accumulators& h = acc[k][t];
        Profile p = profile;
        p.jetPtMin = kJetThresholds[t];
        if (truthBin >= 0) h.photonTruth->Fill(truthPt[tp], w);
        if (recoBin >= 0) h.photonReco->Fill(et[a], w);
        if (photon.matched) h.photonResponse->Fill(truthPt[tp], et[a], w);
        if (photon.fake) h.photonFakes->Fill(et[a], w);
        if (photon.miss) h.photonMisses->Fill(truthPt[tp], w);
        for (int r = 0; r < 4; ++r) if(witness[r]) {
          const double purityEt=classificationEt?classificationEt->at(purityLeaders[r]):et[purityLeaders[r]];
          if(binIndex(purityEt,grids.reco)>=0)h.leakage->Fill(r+.5,purityEt,w);
        }

        // Figure 23: no native response-grid filter and no greedy-link substitute for the direct leading-pair test.
        if(a>=0){
          const long leadReco=associationLeader(p,eta[a],phi[a],jetPt,jetEta,jetPhi);
          if(leadReco>=0){
            h.xjIAnLeading->Fill(et[a],jetPt[leadReco]/et[a],w);
            const long linkedPhoton=match.photonTruth[a];
            if(linkedPhoton>=0 && truthRole[linkedPhoton]==kSignalSourceRole && std::isfinite(truthPt[linkedPhoton]) && truthPt[linkedPhoton]>0){
              const long leadTruth=associationLeader(p,truthEta[linkedPhoton],truthPhi[linkedPhoton],truthJetPt,truthJetEta,truthJetPhi);
              if(leadTruth>=0 && deltaR(jetEta[leadReco],jetPhi[leadReco],truthJetEta[leadTruth],truthJetPhi[leadTruth])<=.3)
                h.xjIAnTagged->Fill(et[a],jetPt[leadReco]/et[a],w);
            }
          }
        }
        /* Truth recoils of the truth photon inside the xJ axis, in orderedByPt order. */
        std::vector<int> truthXjBin(nTruthJets, -1);  /* xJ bin of each selected truth recoil */
        std::vector<long> selected;
        if (truthBin >= 0) {
          auto truthRecoil = [&](size_t j) {
            return truthRecoilJet(p, truthEta[tp], truthPhi[tp], truthJetPt[j], truthJetEta[j], truthJetPhi[j]);
          };
          for (long j : orderedByPt(nTruthJets, truthRecoil, truthJetPt)) {
            const double x = truthJetPt[j] / truthPt[tp];
            const int b = binIndex(x, xj);
            if (b < 0) continue;
            truthXjBin[j] = b;
            selected.push_back(j);
            h.xjTruth->Fill(truthPt[tp], x, w);
          }
        }

        /* Reco recoils of the region-A photon, in orderedByPt order. The order decides which of several reco jets
           linked to one truth recoil is matched: the first; the others are fakes (RecoilClass). */
        std::vector<char> taken(nTruthJets, 0);
        if (a >= 0) {
          auto recoRecoil = [&](size_t j) { return recoilJet(p, eta[a], phi[a], jetPt[j], jetEta[j], jetPhi[j]); };
          bool leading = true;  /* the first reco recoil jet inside the axes is the leading one */
          for (long j : orderedByPt(jetPt.GetSize(), recoRecoil, jetPt)) {
            const double x = jetPt[j] / et[a];
            const int b = binIndex(x, xj);
            if (recoBin < 0 || b < 0) continue;
            h.xjReco->Fill(et[a], x, w);
            const long tj = match.jetTruth[j];
            if (leading) {
              leading = false;
              h.xjLeading->Fill(et[a], x, w);
              if (photonMatched && tj >= 0 && !selected.empty() && tj == selected.front()) h.xjLeadingMatched->Fill(et[a], x, w);
              if(h.xjLeadingUniform){
                h.xjLeadingUniform->Fill(et[a],x,w);
                if(photonMatched && tj>=0 && !selected.empty() && tj==selected.front())h.xjLeadingMatchedUniform->Fill(et[a],x,w);
                if(photonMatched && tj>=0 && truthXjBin[tj]>=0)h.xjLeadingAnyUniform->Fill(et[a],x,w);
              }
            }
            const bool available = tj >= 0 && truthXjBin[tj] >= 0 && !taken[tj];
            // Positive, independent facts; absent bits expose overlapping reasons without repartitioning the response.
            const bool validLinked = tj>=0 && std::isfinite(truthJetPt[tj]) && std::isfinite(truthJetEta[tj]) && std::isfinite(truthJetPhi[tj]);
            const bool hasPhoton=tp>=0 && std::isfinite(truthPt[tp]) && truthPt[tp]>0;
            const std::array<bool,16> facts={{photonMatched,tj>=0,truthBin>=0,hasPhoton,
              validLinked&&truthJetPt[tj]>=p.jetPtMin,validLinked&&std::fabs(truthJetEta[tj])<p.jetAbsEtaMax,
              validLinked&&hasPhoton&&deltaR(truthJetEta[tj],truthJetPhi[tj],truthEta[tp],truthPhi[tp])>=p.overlapDeltaR,
              validLinked&&hasPhoton&&wrappedAbsDeltaPhi(truthJetPhi[tj],truthPhi[tp])>=p.deltaPhiMin,
              validLinked&&hasPhoton&&binIndex(truthJetPt[tj]/truthPt[tp],xj)>=0,available,tj>=0&&bool(taken[tj]),
              nearbyTruthJet(-std::numeric_limits<double>::infinity(),jetEta[j],jetPhi[j],truthJetPt,truthJetEta,truthJetPhi,nTruthJets,p.jetRadius),
              nearbyTruthJet(p.jetPtMin,jetEta[j],jetPhi[j],truthJetPt,truthJetEta,truthJetPhi,nTruthJets,p.jetRadius),
              nTruthJets==0,tj>=0&&!validLinked,hasPhoton&&truthPt[tp]>=17&&truthPt[tp]<33}};
            unsigned mask=0;for(size_t bit=0;bit<facts.size();++bit)if(facts[bit])mask|=1u<<bit;
            const double point[3]={et[a],x,double(mask)};h.pairReason->Fill(point,w);
            const RecoilClass category = recoilClass(photonMatched, available, [&] {
              return nearbyTruthJet(p.jetPtMin, jetEta[j], jetPhi[j], truthJetPt, truthJetEta, truthJetPhi,
                                    nTruthJets, p.jetRadius);
            });
            if (category == RecoilClass::kMatched) {
              h.xjResponse->Fill(truthBin * nXj + truthXjBin[tj] + 0.5, recoBin * nXj + b + 0.5, w);
              taken[tj] = 1;
              continue;
            }
            h.xjFakes->Fill(et[a], x, w);
            if (category == RecoilClass::kFakePhoton) h.xjFakePhoton->Fill(et[a], x, w);
            else if (category == RecoilClass::kCombinatoric) h.xjCombinatoric->Fill(et[a], x, w);
            else h.xjDetectorFakes->Fill(et[a], x, w);
          }
        }
        for (long j : selected)
          if (!taken[j]) h.xjMisses->Fill(truthPt[tp], truthJetPt[j] / truthPt[tp], w);

        /* Every jet of a matched photon at any Delta phi, by whether a truth jet is nearby. */
        if (photonMatched && recoBin >= 0) {
          Profile anyAngle = p;
          anyAngle.deltaPhiMin = 0.0;
          for (size_t j = 0; j < jetPt.GetSize(); ++j) {
            if (!recoilJet(anyAngle, eta[a], phi[a], jetPt[j], jetEta[j], jetPhi[j])) continue;
            const double x = jetPt[j] / et[a];
            if (binIndex(x, xj) < 0) continue;
            const int side = std::fabs(jetEta[j] - eta[a]) > kFarDeltaEta ? 1 : 0;
            const bool nearby =
                nearbyTruthJet(p.jetPtMin, jetEta[j], jetPhi[j], truthJetPt, truthJetEta, truthJetPhi, nTruthJets, p.jetRadius);
            (nearby ? h.dphiTruthJet : h.dphiCombinatoric)[side]->Fill(et[a], x, wrappedAbsDeltaPhi(jetPhi[j], phi[a]), w);
          }
        }
      }
    }
}

using ResponseAccumulators=std::vector<std::vector<Accumulators>>;
struct ResponseView {
  Variation variation;ResponseAccumulators *acc=nullptr;
  const Profile* selection=nullptr;
  Histograms* measurement=nullptr;EventTable* measurementRecords=nullptr;
  int sourceIndex=0;
  purity_inputs::Book* purity=nullptr;selection_efficiency::Book* efficiency=nullptr;
  double expectedStoredGap=std::numeric_limits<double>::quiet_NaN();
};
struct ResponseReadObserver {
  std::function<void(TFile&,TTreeReader&)> opened;
  std::function<void(long long)> completed;
};
/* One native source reader, with the existing event callback applied to each view.
   This is an entry-stream fanout; no retained record or source cache is involved. */
inline void fillFileViews(const std::string& path, const Profile& requested, const std::vector<CentralityClass>& classes,
                     const ResponseGrids& grids, const CentralityWeights& centralityWeight, const std::vector<ResponseView>& views,
                     const char* emsContractText="",
                     const std::function<std::string(const std::string&)>& digest={},
                     const ResponseReadObserver& observer={}) {
  if(views.empty())throw std::runtime_error("response fanout has no views");
  for(const auto &view:views){
    if(!view.acc)throw std::runtime_error("response fanout has no accumulator");
    requireRunnableVariation(view.variation,emsContractText && *emsContractText);
  }
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  const Profile profile = validateSelectionInput(*file, requested);
  if (simulationSystem(*file) != profile.system) throw std::runtime_error("profile system differs from " + path);
  const auto dataset=SelectionJson::parse(metaText(*file,"meta/source")).at("dataset").get<std::string>();
  const std::string purityFamily=(dataset=="isSimInclusive" || dataset=="isSimEmbeddedInclusive")?"inclusivejet":"photonjet";
  std::vector<Profile> profiles;
  for(const auto &view:views){
    validateVariationSystem(view.variation, profile.system);
    if(view.purity && view.purity->family!=purityFamily)throw std::runtime_error("SIM purity book family differs from source");
    if(view.efficiency && purityFamily!="photonjet")throw std::runtime_error("stage efficiencies are photon-family SIM only");
    if(view.selection && !selection_views::cutComposition(view.variation))
      throw std::runtime_error("cut views cannot combine energy/weight effects");
    if((view.measurement || view.purity || view.efficiency) && !selection_views::cutComposition(view.variation))throw std::runtime_error("measurement fanout accepts cut-only nominal momenta/weights");
    if(bool(view.measurement)!=bool(view.measurementRecords))throw std::runtime_error("measurement accumulator/records differ");
    if(view.selection)selection_views::compatible(profile,*view.selection,view.expectedStoredGap);
    profiles.push_back(view.selection?validateSelectionInput(*file,*view.selection):profile);
  }
  const std::uint64_t source = textKey(sourceSha(*file));
  TTreeReader reader("eventTree", file.get());
  if(observer.opened)observer.opened(*file,reader);
  applyEntryRange(reader);
  const VertexEventWeight vertexWeight(*file,reader);
  const PpSimulationGate ppSelection(*file,reader);
  const AuAuTreeWeightCopies treeCopies(*file,reader,centralityWeight);
  TTreeReaderValue<int> record(reader, "evt_record"), ready(reader, "evt_nodes_ready"), vtxValid(reader, "vtx_valid"),
      terminal(reader, "evt_terminal_status");
  std::unique_ptr<TTreeReaderValue<int>> run;std::unique_ptr<TTreeReaderValue<Long64_t>> seq;
  if(std::any_of(views.begin(),views.end(),[](const ResponseView& v){return bool(v.measurement);})){
    run=std::make_unique<TTreeReaderValue<int>>(reader,"evt_run");seq=std::make_unique<TTreeReaderValue<Long64_t>>(reader,"evt_seq");
  }
  TTreeReaderValue<double> vz(reader, "vtx_z"), centrality(reader, "cent_value"), weightValue(reader, "w_event");
  TTreeReaderArray<double> storedEt(reader, "pho_et"), eta(reader, "pho_eta"), phi(reader, "pho_phi");
  PhotonScores score(*file,reader,profile);
  PhotonPreselection presel(reader);
  std::unique_ptr<ems_native::Runtime> emsRuntime;
  std::unique_ptr<energy_selection::NativeInputs> emsInputs;
  if(std::any_of(views.begin(),views.end(),[&](const ResponseView &view){return ems_native::active(profile,view.variation);})){
    emsRuntime=std::make_unique<ems_native::Runtime>(profile,ems_native::parseContract(emsContractText),digest);
    emsInputs=std::make_unique<energy_selection::NativeInputs>(*file->Get<TTree>("eventTree"),reader,profile,score);
  }
  TTreeReaderArray<int> state(reader, profile.stateBranch.c_str());
  TTreeReaderArray<double> cone(reader, "isow_cone_sum"), threshold(reader, "isow_threshold"),
      sideband(reader, "isow_sideband_threshold");
  TTreeReaderArray<int> method(reader, "isow_subtraction_method"), role(reader, "isow_reconstructed_or_truth");
  TTreeReaderArray<int> idClass(reader, "pho_id_class"), isoClass(reader, "pho_iso_class");
  auto branch = [&profile](const char* collection, const char* leaf) { return profile.jetBranch(collection, leaf); };
  TTreeReaderArray<double> storedJetPt(reader, branch("jet", "pt").c_str()), jetEta(reader, branch("jet", "eta").c_str()),
      jetPhi(reader, branch("jet", "phi").c_str());
  std::unique_ptr<TTreeReaderArray<double>> jetEmFraction;
  if (std::any_of(views.begin(),views.end(),[](const ResponseView &view){return view.variation.emScale!=1.;})) {
    const std::string name=branch("jet","emcal_frac");
    if (!file->Get<TTree>("eventTree")->GetBranch(name.c_str())) throw std::runtime_error("EM scale needs " + name);
    jetEmFraction=std::make_unique<TTreeReaderArray<double>>(reader,name.c_str());
  }
  // Equal-distance matches use the producer pre-JES ordering across views.
  // Calibrated pT still controls the analysis cuts and response coordinates.
  TTreeReaderArray<double> rawJetPt(reader, branch("jet", "rawpt").c_str()), unsubRawPt(reader, branch("jetunsub", "rawpt").c_str());
  TTreeReaderValue<int> nJet02(reader, "n_jet02"), nJet03(reader, "n_jet03"), nJet04(reader, "n_jet04"),
      nUnsub02(reader, "n_jetunsub02"), nUnsub03(reader, "n_jetunsub03"), nUnsub04(reader, "n_jetunsub04");
  TTreeReaderArray<double> truthPt(reader, "truthpho_pt"), truthEta(reader, "truthpho_eta"),
      truthPhi(reader, "truthpho_phi");
  TTreeReaderArray<int> truthRole(reader, "truthpho_source_role"), truthG4Valid(reader, "truthpho_g4_photon_valid"),
      truthPromptClass(reader, "truthpho_prompt_class");
  TTreeReaderArray<double> truthIsolation(reader, "truthpho_truth_isolation_r03");
  TTreeReaderArray<double> truthJetPt(reader, branch("truthjet", "pt").c_str()), truthJetEta(reader, branch("truthjet", "eta").c_str()),
      truthJetPhi(reader, branch("truthjet", "phi").c_str());
  TTreeReaderValue<int> nTruthJet02(reader, "n_truthjet02"), nTruthJet03(reader, "n_truthjet03"),
      nTruthJet04(reader, "n_truthjet04");
  TTreeReaderArray<int> plinkReco(reader, "plink_reco"), plinkTruth(reader, "plink_truth"),
      plinkClass(reader, "plink_class");
  TTreeReaderArray<double> plinkMetric(reader, "plink_metric");
  TTreeReaderArray<int> jlinkReco(reader, "jlink_reco"), jlinkTruth(reader, "jlink_truth"),
      jlinkClass(reader, "jlink_class");
  TTreeReaderArray<double> jlinkMetric(reader, "jlink_metric"), unsubPt(reader, branch("jetunsub", "pt").c_str());
  const int radiusIndex = jetRadiusIndex(profile.jetRadius);
  const OwnershipGate owner(reader);

  for(std::size_t i=0;i<views.size();++i) {
    const auto& v=views[i];
    auto context=photonjet_science::measurementContext(*file,profiles[i],centralityWeight,
        photonjet_science::Json::parse(variationJson(v.variation)),
        {{"reco",grids.reco},{"truth",grids.truth},{"xj",grids.xj}},1.,false,"",score.scientificFunctions());
    if(v.variation.centralityOff)for(auto& weight:context["centrality_weights"])weight=1.;
    for(auto& perClass:*v.acc)for(auto& accumulator:perClass)accumulator.scientificContext(context);
    if(v.measurement)v.measurement->scientificContext(context);
    if(v.measurementRecords)v.measurementRecords->scientificContext(context);
  }

  const std::vector<double>& xj = grids.xj;
  const long nXj = xj.size() - 1;
  /* Sums: entries are in source order, the reference's event order. Every fill of one event carries the same weight,
     so the order of the fills inside an event (region, truth photon, recoil order) cannot change any sum. */
  while (reader.Next()) {
    treeCopies.check();
    reportProgress(reader.GetCurrentEntry()+1);
    validateSelectionEvent(profile, *centrality, storedEt, score, state, cone, threshold, sideband, method, role, idClass, isoClass);
    if (*record != 0) { /* upstream-rejected entries complete the shared skip for every view */
      if(observer.completed)observer.completed(reader.GetCurrentEntry());
      continue;
    }
    bool haveFactors=false;std::array<double,3> factors{};
    for(size_t viewIndex=0;viewIndex<views.size();++viewIndex){
    const auto& view=views[viewIndex];
    const auto& selectedProfile=profiles[viewIndex];
    const auto &variation=view.variation;auto &acc=*view.acc;
    /* the reco energies every selection below uses (Variation); the reco-truth match uses the stored ones */
    const std::uint64_t entry = static_cast<std::uint64_t>(reader.GetCurrentEntry());
    Varied et(storedEt, variation.photonResolution=="off" ? variation.photonScale*variation.emScale : 1., variation.photonSmear, variation.seed, source, entry, 0);
    Varied jetPt(storedJetPt, variation.jetResolution=="off" ? variation.jetScale : 1., variation.jetSmear, variation.seed, source, entry, 1);
    const size_t nPhoton = et.GetSize(), nTruthJets = truthJetPt.GetSize();
    if (score.GetSize() != nPhoton || (!profile.optionalID && state.GetSize() != nPhoton))
      throw std::runtime_error("photon-ID branches pho_bdt_score, pho_bdt_state are not filled in " + path);
    if (profile.augmentedClasses && (idClass.GetSize() != nPhoton || isoClass.GetSize() != nPhoton))
      throw std::runtime_error("the augmented classes pho_id_class, pho_iso_class are not filled in " + path);
    /* Matches are read for every event: the reference refuses a duplicate final match anywhere in its input. */
    auto sizes = [](int a, int b, int c) { return std::array<size_t, 3>{static_cast<size_t>(a), static_cast<size_t>(b), static_cast<size_t>(c)}; };
    const JetRanges ranges = JetRanges::of(sizes(*nJet02, *nJet03, *nJet04), sizes(*nUnsub02, *nUnsub03, *nUnsub04),
                                           sizes(*nTruthJet02, *nTruthJet03, *nTruthJet04), radiusIndex);
    if (ranges.nJet != storedJetPt.GetSize() || ranges.nUnsub != unsubPt.GetSize() || ranges.nTruthJet != nTruthJets)
      throw std::runtime_error("jet collection counters differ from their arrays in " + path);
    if (rawJetPt.GetSize() != ranges.nJet || unsubRawPt.GetSize() != ranges.nUnsub)
      throw std::runtime_error("raw jet collection sizes differ from their analysis views in " + path);
    const EventMatches match =
        eventMatches(nPhoton, truthPt.GetSize(), plinkReco, plinkTruth, plinkClass, plinkMetric, plinkReco.GetSize(),
                     ranges, rawJetPt, unsubRawPt, jlinkReco, jlinkTruth, jlinkClass, jlinkMetric, jlinkReco.GetSize());
    // Named resolution modes first form the calibrated nominal realization;
    // scale variations then move that same realization coherently.
    if(variation.photonResolution!="off") {
      et.smearPhotons(truthPt,match.photonTruth,variation,source,entry);
      et.scaleAll(variation.photonScale*variation.emScale);
    }
    if(variation.jetResolution!="off") {
      jetPt.smearJets(truthJetPt,match.jetTruth,variation,profile.jetRadius,source,entry);
      jetPt.scaleAll(variation.jetScale);
    }
    if(jetEmFraction)jetPt.scaleEmFraction(*jetEmFraction,variation.emScale);
    /* The producer's signal flag must be the analysis truth signal of the stored primitives. */
    for (size_t t = 0; t < truthPt.GetSize(); ++t)
      if (analysisTruthSignal(truthG4Valid[t], truthPromptClass[t], truthEta[t], truthIsolation[t]) !=
          (truthRole[t] == kSignalSourceRole))
        throw std::runtime_error("truth signal flag differs from its definition in " + path);
    if (*ready != 1 || *vtxValid != 1 || !acceptedEventClass(*terminal, *centrality)) continue;
    if (!ppSelection.accepts()) continue;
    if (!owner.owned()) continue;  /* a stitched sample's events outside its window belong to another sample */
    if (!responseEvent(profile, *vz)) continue;
    /* One weight reaches every response, leakage and simulated-transfer fill.
       Keep the original multiplication order for the nominal/on path. */
    if(!haveFactors){
      factors={treeCopies.baseWeight(ppSelection.weight(*weightValue)),vertexWeight(),centralityWeight(*centrality)};
      haveFactors=true;
    }
    const double w = responseEventWeight(variation,factors[0],factors[1],factors[2]);
    if (cone.GetSize() != 2 * nPhoton) throw std::runtime_error("isolation witnesses are not 2 per photon in " + path);

    selection_views::validateStoredGap(profile,view.expectedStoredGap,nPhoton,cone,threshold,sideband,method,role);
    if(view.measurement || view.purity || view.efficiency){
      size_t ci=0;
      if(profile.system==System::auau){while(ci<classes.size() && !(classes[ci].low<=*centrality && *centrality<classes[ci].high))++ci;}
      if(ci<classes.size()){
        if(view.measurement)view.measurement->h1("events_"+classes[ci].name).Fill(.5,w);
        std::vector<int> accepted(nPhoton),visible(nPhoton,1);
        for(size_t i=0;i<nPhoton;++i)accepted[i]=presel[i];
        if(view.measurement)accumulateMeasurement(selectedProfile,classes[ci],**run,**seq,reader.GetCurrentEntry(),view.sourceIndex,*centrality,*vz,w,
          *view.measurement,view.measurementRecords,storedEt,eta,phi,score,accepted,visible,state,cone,threshold,sideband,
          method,role,idClass,isoClass,storedJetPt,jetEta,jetPhi,1.);
        if(view.purity){
          std::vector<int> signal(nPhoton,0);
          for(size_t i=0;i<nPhoton;++i)signal[i]=match.photonTruth[i]>=0 && truthRole[match.photonTruth[i]]==kSignalSourceRole;
          view.purity->accumulate(selectedProfile,classes[ci],*vz,*centrality,w,storedEt,eta,phi,score,presel,state,
            cone,threshold,sideband,method,role,idClass,isoClass,storedJetPt,jetEta,jetPhi,signal,1.);
        }
        if(view.efficiency)view.efficiency->accumulate(selectedProfile,classes[ci],*vz,*centrality,w,storedEt,eta,phi,score,presel,state,
            cone,threshold,sideband,method,role,idClass,isoClass,truthPt,truthEta,truthRole,match.photonTruth);
      }
    }
    std::vector<double> classificationEt;
    if(responseOnly(variation)){classificationEt.reserve(storedEt.GetSize());for(size_t i=0;i<storedEt.GetSize();++i)classificationEt.push_back(storedEt[i]);}
    if(emsRuntime && ems_native::active(profile,variation)){
      const auto selected=emsRuntime->evaluate(*emsInputs,storedEt,*vz,*centrality,variation.emScale);
      accumulateResponse(profile,classes,grids,acc,*vz,*centrality,w,et,eta,phi,selected.score,selected.preselection,selected.state,
        selected.cone,selected.threshold,selected.sideband,method,role,idClass,isoClass,jetPt,jetEta,jetPhi,truthPt,truthEta,
        truthPhi,truthRole,truthJetPt,truthJetEta,truthJetPhi,match,&selected.et);
    }else{
    accumulateResponse(selectedProfile,classes,grids,acc,*vz,*centrality,w,et,eta,phi,score,presel,state,
        cone,threshold,sideband,method,role,idClass,isoClass,jetPt,jetEta,jetPhi,truthPt,truthEta,
        truthPhi,truthRole,truthJetPt,truthJetEta,truthJetPhi,match,responseOnly(variation)?&classificationEt:nullptr);
    }
    } // views
    if(observer.completed)observer.completed(reader.GetCurrentEntry());
  }
  /* a missing branch or a read error ends the loop like the end of the tree; only the end is accepted */
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
}

inline void fillFile(const std::string& path,const Profile& profile,const std::vector<CentralityClass>& classes,
                     const ResponseGrids& grids,const CentralityWeights& weights,const Variation& variation,
                     ResponseAccumulators& acc,const char* emsContractText="",
                     const std::function<std::string(const std::string&)>& digest={}) {
  fillFileViews(path,profile,classes,grids,weights,{{variation,&acc}},emsContractText,digest);
}

inline void writeProduct(const char* output,const Profile& profile,const ResponseGrids& grids,
                         const CentralityWeights& centralityWeight,const Variation& variation,
                         ResponseAccumulators& acc,const std::string& sources,const char* emsContractText="",
                         const std::string& gridProvenance="") {
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  for (auto& perClass : acc)
    for (auto& h : perClass) h.write();
  TObjString(profile.name.c_str()).Write("profile");
  writeSelection(profile);
  TObjString(("{\"reco_pt\":" + jsonArray(grids.reco) + ",\"truth_pt\":" + jsonArray(grids.truth) + ",\"xj\":" +
              jsonArray(grids.xj) + "}").c_str()).Write("grids");
  TObjString(("{\"per_slice\":" + jsonArray(centralityWeight.perSlice) + "}").c_str()).Write("centrality_weights");
  TObjString(variationJson(variation).c_str()).Write("variation");
  if(emsSelectionResponse(variation))TObjString(emsContractText).Write("ems_selection_contract");
  if(emsSelectionResponse(variation))TObjString("true").Write("ems_applied_to_simulation");
  TObjString(sources.c_str()).Write("sources");
  std::vector<std::string> nativeInputs;for(const auto& source:photonjet_science::Json::parse(sources))nativeInputs.push_back(source.at("path").get<std::string>());
  photonjet_science::writeProducerEvidence(*out,photonjet_science::producerInputs(nativeInputs));
  if(!gridProvenance.empty())TObjString(gridProvenance.c_str()).Write("native_response_grid");
  out->Close();
}

inline int run(const char* inputs, const char* output, const char* profileName, const char* recoEdges,
               const char* truthEdges, const char* xjAxisEdges, const char* centralityWeightText, const char* variationText,
               const char* selection,const char* emsContractText="",
               const std::function<std::string(const std::string&)>& digest={}) {
  const Variation variation = parseVariation(variationText);
  requireRunnableVariation(variation,emsContractText && *emsContractText);
  if(emsSelectionResponse(variation))ems_native::parseContract(emsContractText);
  const ResponseGrids defaults;
  ResponseGrids grids;
  grids.reco = parseEdges(recoEdges, defaults.reco);
  grids.truth = parseEdges(truthEdges, defaults.truth);
  grids.xj = parseEdges(xjAxisEdges, defaults.xj);
  const std::vector<std::string> paths = inputPaths(inputs);
  System system;
  {
    std::unique_ptr<TFile> first(TFile::Open(paths.front().c_str(), "READ"));
    if (!first || first->IsZombie()) throw std::runtime_error("cannot open " + paths.front());
    system = simulationSystem(*first);
  }
  validateVariationSystem(variation, system);
  const CentralityWeights centralityWeight = centralityWeights(centralityWeightText, system);
  const Profile profile = configuredProfile(profileName, system, selection);
  const std::vector<CentralityClass> classes = centralityClasses(system);
  std::vector<std::vector<Accumulators>> acc(classes.size());
  for (size_t k = 0; k < classes.size(); ++k)
    for (double t : kJetThresholds) acc[k].emplace_back("_" + classes[k].name + "_" + thresholdName(t), grids,profile.system==System::pp);
  for (const auto& path : paths) fillFile(path, profile, classes, grids, centralityWeight, variation, acc,emsContractText,digest);
  const std::string sources = sourceRecords(paths);  /* before the output is opened (sourceRecord) */
  writeProduct(output,profile,grids,centralityWeight,variation,acc,sources,emsContractText);
  return 0;
}

}  /* namespace response */
}  /* namespace photonjet */

/* inputs: whitespace-separated eventTree simulation paths. profile: "nominal" (the analysis), a variant of it
   (PhotonJet.h namedProfile), or "historical" (the reference reproduction only).
   Optional axes (whitespace-separated edges; empty keeps the reference axes): reco photon pT, truth photon pT, xJ.
   centralityWeights: Au+Au only and required there, 16 weights for the 5% slices of 0-80% (CentralityWeights).
   variation: "" (nominal) or the reco energy scale and resolution variation of the systematics (Variation). */
/* Simulation responses, misses and fake categories use the plan's native
   grids. FinalAnalysis restricts the measurement window and applies the
   sample's absolute normalization when it reads the saved package. */
int FillPhotonJetResponse(const char* inputs, const char* output, const char* profileName = "nominal",
                          const char* recoPtEdges = "", const char* truthPtEdges = "", const char* xjEdges = "",
                          const char* centralityWeights = "", const char* variation = "", const char* selection = "",
                          const char* emsContractText="",const std::function<std::string(const std::string&)>& digest={}) {
  return photonjet::response::run(inputs, output, profileName, recoPtEdges, truthPtEdges, xjEdges, centralityWeights,
                                  variation, selection,emsContractText,digest);
}
